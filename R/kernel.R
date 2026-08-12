# Bounded additive contraction ----------------------------------------------

.tile_starts <- function(n, size) seq.int(1L, n, by = size)

.validate_tile_size <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x %% 1 != 0) {
    stop(sprintf("`%s` must be one positive integer.", name), call. = FALSE)
  }
  as.integer(x)
}

.tiled_contraction <- function(weights, atoms, row_tile, coordinate_tile,
                               feature_tile, write_tile = NULL) {
  if (!is.matrix(weights) || !is.numeric(weights) || any(!is.finite(weights))) {
    stop("`weights` must be a finite numeric matrix.", call. = FALSE)
  }
  if (!is.matrix(atoms) || !is.numeric(atoms) || any(!is.finite(atoms))) {
    stop("`atoms` must be a finite numeric matrix.", call. = FALSE)
  }
  if (ncol(weights) != nrow(atoms)) {
    stop("The feature dimension of `weights` and `atoms` must agree.",
      call. = FALSE)
  }
  if (nrow(weights) < 1L || ncol(weights) < 1L || ncol(atoms) < 1L) {
    stop("Contraction inputs must have positive dimensions.", call. = FALSE)
  }
  row_tile <- .validate_tile_size(row_tile, "row_tile")
  coordinate_tile <- .validate_tile_size(coordinate_tile, "coordinate_tile")
  feature_tile <- .validate_tile_size(feature_tile, "feature_tile")
  if (!is.null(write_tile) && !is.function(write_tile)) {
    stop("`write_tile` must be NULL or a function.", call. = FALSE)
  }

  measurements <- nrow(weights)
  features <- ncol(weights)
  coordinates <- ncol(atoms)
  output <- if (is.null(write_tile)) matrix(0, measurements, coordinates) else NULL
  tile_count <- 0L
  max_temporary_elements <- 0L

  for (row_start in .tile_starts(measurements, row_tile)) {
    rows <- row_start:min(row_start + row_tile - 1L, measurements)
    for (coordinate_start in .tile_starts(coordinates, coordinate_tile)) {
      coordinates_in_tile <- coordinate_start:min(
        coordinate_start + coordinate_tile - 1L, coordinates
      )
      tile <- matrix(0, length(rows), length(coordinates_in_tile))
      max_temporary_elements <- max(max_temporary_elements, length(tile))

      for (feature_start in .tile_starts(features, feature_tile)) {
        features_in_tile <- feature_start:min(
          feature_start + feature_tile - 1L, features
        )
        tile <- tile +
          weights[rows, features_in_tile, drop = FALSE] %*%
          atoms[features_in_tile, coordinates_in_tile, drop = FALSE]
      }

      if (is.null(write_tile)) {
        output[rows, coordinates_in_tile] <- tile
      } else {
        write_tile(rows, coordinates_in_tile, tile)
      }
      tile_count <- tile_count + 1L
    }
  }

  list(
    value = output,
    diagnostics = list(
      row_tile = row_tile,
      coordinate_tile = coordinate_tile,
      feature_tile = feature_tile,
      tile_count = tile_count,
      max_temporary_elements = max_temporary_elements
    )
  )
}

# Primary additive-frame lowering. Relation blocks are read exactly once per
# canonical feature block and converted to packed cross-Gram atoms before the
# sparse frame distributes them across measurement rows.
.streamed_crossgram_contraction <- function(frame, read_relation, partitions,
                                            effects, over,
                                            feature_block = 1024L,
                                            row_tile = 1024L,
                                            coordinate_tile = 256L,
                                            accumulate_tile = NULL,
                                            retain_local_relations = FALSE,
                                            query = NULL,
                                            form_total = TRUE,
                                            task_observer = NULL) {
  .validate_frame_for_compile(frame)
  if (!identical(frame$representation, "additive_diagonal")) {
    stop("The streamed cross-Gram lowering requires an additive diagonal frame.",
      call. = FALSE)
  }
  if (!is.function(read_relation)) {
    stop("`read_relation` must be a feature-block reader function.", call. = FALSE)
  }
  if (!is.character(partitions) || length(partitions) < 1L || anyNA(partitions) ||
      any(!nzchar(partitions)) || anyDuplicated(partitions)) {
    stop("`partitions` must contain unique nonempty identifiers.", call. = FALSE)
  }
  if (!is.character(effects) || length(effects) < 1L || anyNA(effects) ||
      any(!nzchar(effects)) || anyDuplicated(effects)) {
    stop("`effects` must contain unique nonempty identifiers.", call. = FALSE)
  }
  .validate_pairing(over)
  if (anyNA(match(over$left, partitions)) ||
      anyNA(match(over$right, partitions))) {
    stop("Every pairing endpoint must identify a declared partition.", call. = FALSE)
  }
  feature_block <- .validate_tile_size(feature_block, "feature_block")
  row_tile <- .validate_tile_size(row_tile, "row_tile")
  coordinate_tile <- .validate_tile_size(coordinate_tile, "coordinate_tile")
  if (!is.null(accumulate_tile) && !is.function(accumulate_tile)) {
    stop("`accumulate_tile` must be NULL or a function.", call. = FALSE)
  }
  if (!is.logical(retain_local_relations) ||
      length(retain_local_relations) != 1L || is.na(retain_local_relations)) {
    stop("`retain_local_relations` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(form_total) || length(form_total) != 1L || is.na(form_total)) {
    stop("`form_total` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!form_total && !retain_local_relations) {
    stop("Streaming must form total output, local relations, or both.",
      call. = FALSE)
  }
  if (!is.null(task_observer) && !is.function(task_observer)) {
    stop("`task_observer` must be NULL or a function.", call. = FALSE)
  }

  features <- ncol(frame$weights)
  q <- length(effects)
  packed_width <- q * (q + 1L) / 2L
  if (!is.null(query) &&
      (!is.matrix(query) || !is.numeric(query) || nrow(query) != packed_width ||
       ncol(query) < 1L || any(!is.finite(query)))) {
    stop("`query` must be a finite packed-coordinate-by-view matrix.",
      call. = FALSE)
  }
  output_width <- if (is.null(query)) packed_width else ncol(query)
  reducer <- .new_crossgram_reducer(
    frame = frame,
    partitions = partitions,
    effects = effects,
    output_width = output_width,
    row_tile = row_tile,
    coordinate_tile = coordinate_tile,
    accumulate_tile = accumulate_tile,
    retain_local_relations = retain_local_relations,
    form_total = form_total
  )

  for (feature_start in .tile_starts(features, feature_block)) {
    feature_ids <- feature_start:min(feature_start + feature_block - 1L, features)
    if (!is.null(task_observer)) task_observer("started", feature_ids)
    task <- tryCatch(
      {
        relations <- stats::setNames(lapply(partitions, function(partition) {
          value <- read_relation(partition, feature_ids)
          if (!is.matrix(value) || !is.numeric(value) ||
              !identical(dim(value), c(q, length(feature_ids))) ||
              any(!is.finite(value))) {
            stop("Relation reader returned an invalid effect-by-feature block.",
              call. = FALSE)
          }
          value
        }), partitions)
        .crossgram_feature_task(
          relations = relations,
          feature_ids = feature_ids,
          effects = effects,
          partitions = partitions,
          over = over,
          query = query,
          form_atoms = form_total
        )
      },
      error = function(error) {
        if (!is.null(task_observer)) task_observer("failed", feature_ids)
        stop(error)
      }
    )
    tryCatch(
      .reduce_crossgram_task(reducer, task),
      error = function(error) {
        if (!is.null(task_observer)) task_observer("failed", feature_ids)
        stop(error)
      }
    )
    if (!is.null(task_observer)) task_observer("completed", feature_ids)
  }

  list(
    value = reducer$output,
    local_relations = reducer$local_relations,
    diagnostics = reducer$diagnostics
  )
}

.coherent_geometry_from_local <- function(local_relations, over, mass,
                                          row_tile = 1024L,
                                          write_tile = NULL, query = NULL) {
  if (!is.array(local_relations) || length(dim(local_relations)) != 3L ||
      !is.numeric(local_relations) || any(!is.finite(local_relations))) {
    stop("`local_relations` must be a finite measurement-by-effect-by-partition array.",
      call. = FALSE)
  }
  .validate_pairing(over)
  partition_ids <- dimnames(local_relations)[[3L]]
  if (is.null(partition_ids) || anyNA(partition_ids) ||
      any(!nzchar(partition_ids)) || anyDuplicated(partition_ids)) {
    stop("Local relation partitions must have unique nonempty names.", call. = FALSE)
  }
  left_index <- match(over$left, partition_ids)
  right_index <- match(over$right, partition_ids)
  if (anyNA(left_index) || anyNA(right_index)) {
    stop("Every pairing endpoint must identify a local relation partition.",
      call. = FALSE)
  }
  measurements <- dim(local_relations)[[1L]]
  q <- dim(local_relations)[[2L]]
  if (!is.numeric(mass) || !(length(mass) %in% c(1L, measurements)) ||
      any(!is.finite(mass)) || any(mass <= 0)) {
    stop("`mass` must be one positive finite value or one per measurement.",
      call. = FALSE)
  }
  mass <- rep_len(mass, measurements)
  row_tile <- .validate_tile_size(row_tile, "row_tile")
  if (!is.null(write_tile) && !is.function(write_tile)) {
    stop("`write_tile` must be NULL or a function.", call. = FALSE)
  }
  packed_width <- q * (q + 1L) / 2L
  if (!is.null(query) &&
      (!is.matrix(query) || !is.numeric(query) || nrow(query) != packed_width ||
       ncol(query) < 1L || any(!is.finite(query)))) {
    stop("`query` must be a finite packed-coordinate-by-view matrix.",
      call. = FALSE)
  }
  output_width <- if (is.null(query)) packed_width else ncol(query)
  output <- if (is.null(write_tile)) matrix(0, measurements, output_width) else NULL
  max_work_bytes <- 0
  tile_count <- 0L

  for (row_start in .tile_starts(measurements, row_tile)) {
    rows <- row_start:min(row_start + row_tile - 1L, measurements)
    tile <- matrix(0, length(rows), packed_width)
    coordinate <- 0L
    for (column in seq_len(q)) {
      for (row in column:q) {
        coordinate <- coordinate + 1L
        work <- numeric(length(rows))
        for (edge in seq_len(nrow(over))) {
          work <- work + over$weight[[edge]] * 0.5 * (
            local_relations[rows, row, left_index[[edge]]] *
              local_relations[rows, column, right_index[[edge]]] +
            local_relations[rows, row, right_index[[edge]]] *
              local_relations[rows, column, left_index[[edge]]]
          )
        }
        work <- work / mass[rows]
        if (row != column) work <- sqrt(2) * work
        tile[, coordinate] <- work
        max_work_bytes <- max(max_work_bytes,
          as.double(utils::object.size(work)))
      }
    }
    if (!is.null(query)) tile <- tile %*% query
    if (is.null(write_tile)) {
      output[rows, ] <- tile
    } else {
      write_tile(rows, seq_len(output_width), tile)
    }
    tile_count <- tile_count + 1L
  }

  list(
    value = output,
    diagnostics = list(
      row_tile = row_tile,
      tile_count = tile_count,
      max_tile_bytes = as.double(utils::object.size(
        matrix(0, min(row_tile, measurements), output_width)
      )),
      max_work_bytes = max_work_bytes,
      durable_output_bytes = if (is.null(output)) 0 else
        as.double(utils::object.size(output))
    )
  )
}
