# Variable-block measurement-form storage ----------------------------------

.measurement_block_index <- function(edges, left_frame, right_frame = left_frame) {
  edges <- .validate_measurement_edges(edges, left_frame, right_frame)
  table <- edges$edges
  lengths <- vapply(seq_len(nrow(table)), function(index) {
    nrow(left_frame$legs[[table$left[[index]]]]$operator) *
      nrow(right_frame$legs[[table$right[[index]]]]$operator)
  }, numeric(1))
  offsets <- c(0, utils::head(cumsum(lengths), -1L))
  value <- data.frame(
    edge_id = paste0("edge_", sprintf("%06d", seq_len(nrow(table)))),
    left = table$left,
    right = table$right,
    weight = table$weight,
    left_axis = vapply(table$left, function(node) {
      left_frame$legs[[node]]$output_space$signature
    }, character(1)),
    right_axis = vapply(table$right, function(node) {
      right_frame$legs[[node]]$output_space$signature
    }, character(1)),
    d_left = as.integer(vapply(table$left, function(node) {
      nrow(left_frame$legs[[node]]$operator)
    }, integer(1))),
    d_right = as.integer(vapply(table$right, function(node) {
      nrow(right_frame$legs[[node]]$operator)
    }, integer(1))),
    offset_elements = as.double(offsets),
    length_elements = as.double(lengths),
    stringsAsFactors = FALSE
  )
  semantic <- list(
    schema_version = 1L,
    edge_set = edges$signature,
    left_frame = left_frame$signature,
    right_frame = right_frame$signature,
    index = lapply(value, unname)
  )
  structure(value,
    edge_set = edges$signature,
    left_frame = left_frame$signature,
    right_frame = right_frame$signature,
    signature = .sha256_signature(semantic),
    class = c("effect_measurement_block_index", "data.frame")
  )
}

.validate_measurement_block_index <- function(index) {
  expected <- c("edge_id", "left", "right", "weight", "left_axis",
    "right_axis", "d_left", "d_right", "offset_elements", "length_elements")
  if (!inherits(index, "effect_measurement_block_index") ||
      !is.data.frame(index) || !identical(names(index), expected) ||
      nrow(index) < 1L || anyNA(index) ||
      anyDuplicated(index$edge_id) || any(!nzchar(index$edge_id)) ||
      !identical(index$edge_id,
        paste0("edge_", sprintf("%06d", seq_len(nrow(index))))) ||
      any(index$d_left < 1L) || any(index$d_right < 1L) ||
      any(index$length_elements !=
        as.double(index$d_left) * as.double(index$d_right)) ||
      !identical(index$offset_elements,
        c(0, utils::head(cumsum(index$length_elements), -1L))) ||
      !is.character(attr(index, "edge_set", exact = TRUE)) ||
      !is.character(attr(index, "left_frame", exact = TRUE)) ||
      !is.character(attr(index, "right_frame", exact = TRUE)) ||
      !is.character(attr(index, "signature", exact = TRUE))) {
    .input_error("Measurement block index is missing or noncanonical.")
  }
  semantic <- list(
    schema_version = 1L,
    edge_set = attr(index, "edge_set", exact = TRUE),
    left_frame = attr(index, "left_frame", exact = TRUE),
    right_frame = attr(index, "right_frame", exact = TRUE),
    index = lapply(index, unname)
  )
  expected_signature <- .sha256_signature(semantic)
  if (!identical(attr(index, "signature", exact = TRUE),
      expected_signature)) {
    .contract_error("Measurement block-index identity is inconsistent.")
  }
  index
}

.measurement_store_manifest <- function(index, complete, written,
                                        representation, expected_bytes) {
  index <- .validate_measurement_block_index(index)
  if (!.is_flag(complete) || !is.logical(written) ||
      length(written) != nrow(index) || anyNA(written)) {
    .input_error("Measurement-store completion metadata is invalid.")
  }
  list(
    schema_version = 1L,
    codec = "measurement-block-double-v1",
    representation = representation,
    complete = complete,
    written = unname(written),
    index_signature = attr(index, "signature", exact = TRUE),
    block_count = nrow(index),
    expected_bytes = as.double(expected_bytes)
  )
}

.measurement_store_signature <- function(index, manifest) {
  .sha256_signature(list(
    schema_version = 1L,
    index_signature = attr(index, "signature", exact = TRUE),
    codec = manifest$codec,
    block_count = manifest$block_count
  ))
}

.measurement_edge_position <- function(index, edge) {
  if (is.numeric(edge) && length(edge) == 1L && !is.na(edge) &&
      is.finite(edge) && edge %% 1 == 0) {
    position <- as.integer(edge)
  } else if (is.character(edge) && length(edge) == 1L && !is.na(edge)) {
    position <- match(edge, index$edge_id)
  } else {
    position <- NA_integer_
  }
  if (is.na(position) || position < 1L || position > nrow(index)) {
    .input_error(
      "Measurement edge identifier is not present in the block index."
    )
  }
  position
}

.memory_measurement_store <- function(blocks, index) {
  index <- .validate_measurement_block_index(index)
  if (!is.list(blocks) || length(blocks) != nrow(index) ||
      !identical(names(blocks), index$edge_id)) {
    .input_error("Memory measurement store requires one named block per edge.")
  }
  blocks <- lapply(seq_along(blocks), function(position) {
    value <- blocks[[position]]
    expected <- c(index$d_left[[position]], index$d_right[[position]])
    if (!.is_finite_matrix(value) || !identical(dim(value), expected)) {
      .input_error("Measurement block has invalid dimensions or values.")
    }
    unname(value)
  })
  names(blocks) <- index$edge_id
  expected_bytes <- sum(index$length_elements) * 8
  manifest <- .measurement_store_manifest(
    index, TRUE, rep(TRUE, nrow(index)), "memory", expected_bytes
  )
  structure(list(
    representation = "memory",
    codec = manifest$codec,
    index = index,
    manifest = manifest,
    path = NULL,
    manifest_path = NULL,
    read = function(edge) {
      blocks[[.measurement_edge_position(index, edge)]]
    },
    signature = .measurement_store_signature(index, manifest)
  ), class = "effect_measurement_store")
}

.measurement_manifest_path <- function(path) paste0(path, ".manifest.rds")

.file_measurement_store <- function(path, index, create = FALSE) {
  .check_string(path, "path", what = "one nonempty measurement-store path")
  index <- .validate_measurement_block_index(index)
  manifest_path <- .measurement_manifest_path(path)
  expected_bytes <- sum(index$length_elements) * 8
  if (isTRUE(create)) {
    if (file.exists(path) || file.exists(manifest_path)) {
      .input_error("Refusing to overwrite an existing measurement store.")
    }
    connection <- file(path, open = "w+b")
    on.exit(if (!is.null(connection)) try(close(connection), silent = TRUE),
      add = TRUE)
    remaining <- sum(index$length_elements)
    zero <- numeric(min(8192, remaining))
    while (remaining > 0) {
      count <- min(remaining, length(zero))
      writeBin(zero[seq_len(count)], connection, size = 8, endian = "little")
      remaining <- remaining - count
    }
    close(connection)
    connection <- NULL
    manifest <- .measurement_store_manifest(
      index, FALSE, rep(FALSE, nrow(index)), "block_backed", expected_bytes
    )
    saveRDS(list(manifest = manifest, index = index), manifest_path,
      version = 3)
  }
  if (!file.exists(path) || !file.exists(manifest_path) ||
      is.na(file.info(path)$size) || file.info(path)$size != expected_bytes) {
    .input_error(
      "Measurement store payload or manifest is missing or has the wrong size."
    )
  }
  persisted <- readRDS(manifest_path)
  if (!is.list(persisted) ||
      !identical(names(persisted), c("manifest", "index")) ||
      !identical(persisted$index, index)) {
    .contract_error(
      "Measurement store manifest does not match its block index."
    )
  }
  manifest <- persisted$manifest
  structure(list(
    representation = "block_backed",
    codec = "measurement-block-double-v1",
    index = index,
    manifest = manifest,
    path = normalizePath(path, mustWork = TRUE),
    manifest_path = normalizePath(manifest_path, mustWork = TRUE),
    read = function(edge) {
      position <- .measurement_edge_position(index, edge)
      row <- index[position, , drop = FALSE]
      connection <- file(path, open = "rb")
      on.exit(close(connection), add = TRUE)
      seek(connection, where = row$offset_elements[[1L]] * 8,
        origin = "start", rw = "read")
      values <- readBin(connection, "double",
        n = row$length_elements[[1L]], size = 8, endian = "little")
      matrix(values, row$d_left[[1L]], row$d_right[[1L]])
    },
    signature = .measurement_store_signature(index, manifest)
  ), class = "effect_measurement_store")
}

.write_measurement_block <- function(store, edge, value) {
  if (!inherits(store, "effect_measurement_store") ||
      store$representation != "block_backed") {
    .input_error("Only a block-backed measurement store can be written.")
  }
  position <- .measurement_edge_position(store$index, edge)
  row <- store$index[position, , drop = FALSE]
  if (!.is_finite_matrix(value) ||
      !identical(dim(value), c(row$d_left[[1L]], row$d_right[[1L]]))) {
    .input_error("Measurement block write has invalid dimensions or values.")
  }
  persisted <- readRDS(store$manifest_path)
  if (persisted$manifest$written[[position]]) {
    .input_error("Refusing to overwrite an already written measurement block.")
  }
  connection <- file(store$path, open = "r+b")
  on.exit(if (!is.null(connection)) try(close(connection), silent = TRUE),
    add = TRUE)
  seek(connection, where = row$offset_elements[[1L]] * 8,
    origin = "start", rw = "write")
  writeBin(as.double(value), connection, size = 8, endian = "little")
  close(connection)
  connection <- NULL
  persisted$manifest$written[[position]] <- TRUE
  persisted$manifest$complete <- all(persisted$manifest$written)
  saveRDS(persisted, store$manifest_path, version = 3)
  .file_measurement_store(store$path, store$index, create = FALSE)
}

.validate_measurement_store <- function(store, require_complete = TRUE,
                                        probe = TRUE) {
  expected <- c("representation", "codec", "index", "manifest", "path",
    "manifest_path", "read", "signature")
  if (!.sealed_fields(store, "effect_measurement_store", expected) ||
      !store$representation %in% c("memory", "block_backed") ||
      !identical(store$codec, "measurement-block-double-v1") ||
      !is.function(store$read)) {
    .input_error("Measurement-store fields are missing or noncanonical.")
  }
  index <- .validate_measurement_block_index(store$index)
  manifest <- store$manifest
  if (store$representation == "block_backed") {
    persisted <- readRDS(store$manifest_path)
    if (!identical(persisted$index, index) ||
        !identical(persisted$manifest, manifest)) {
      .contract_error(
        "Measurement-store object is stale relative to its manifest."
      )
    }
    if (!file.exists(store$path) ||
        file.info(store$path)$size != manifest$expected_bytes) {
      .input_error("Measurement-store payload is incomplete.")
    }
  }
  expected_manifest <- .measurement_store_manifest(
    index, manifest$complete, manifest$written,
    store$representation, sum(index$length_elements) * 8
  )
  if (!identical(manifest, expected_manifest) ||
      !identical(store$signature,
        .measurement_store_signature(index, manifest))) {
    .contract_error("Measurement-store manifest or identity is inconsistent.")
  }
  if (isTRUE(require_complete) && !isTRUE(manifest$complete)) {
    .input_error("Measurement store is incomplete for its declared edge set.")
  }
  if (isTRUE(probe)) {
    positions <- if (manifest$complete) seq_len(nrow(index)) else
      which(manifest$written)
    for (position in positions) {
      value <- store$read(position)
      if (!.is_finite_matrix(value) ||
          !identical(dim(value), c(index$d_left[[position]], index$d_right[[position]]))) {
        .input_error(
          "Measurement-store reader returned an invalid logical block."
        )
      }
    }
  }
  invisible(store)
}

.read_measurement_store <- function(store, edge) {
  .validate_measurement_store(store, require_complete = TRUE, probe = FALSE)
  value <- store$read(edge)
  position <- .measurement_edge_position(store$index, edge)
  if (!.is_finite_matrix(value) ||
      !identical(dim(value), c(store$index$d_left[[position]], store$index$d_right[[position]]))) {
    .input_error("Measurement-store reader returned an invalid logical block.")
  }
  value
}
