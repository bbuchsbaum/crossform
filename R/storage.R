# File-backed packed geometry storage ---------------------------------------

.file_geometry_store <- function(path, dim, create = FALSE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be one nonempty file path.", call. = FALSE)
  }
  if (!is.numeric(dim) || length(dim) != 2L || any(!is.finite(dim)) ||
      any(dim < 1) || any(dim %% 1 != 0)) {
    stop("`dim` must contain two positive integers.", call. = FALSE)
  }
  dim <- as.integer(dim)
  expected_bytes <- prod(as.double(dim)) * 8
  if (create) {
    if (file.exists(path)) stop("Refusing to overwrite an existing geometry store.", call. = FALSE)
    connection <- file(path, open = "w+b")
    values_remaining <- prod(as.double(dim))
    zero_block <- rep(0, min(values_remaining, 8192))
    while (values_remaining > 0) {
      count <- min(values_remaining, length(zero_block))
      writeBin(zero_block[seq_len(count)], connection, size = 8, endian = "little")
      values_remaining <- values_remaining - count
    }
    close(connection)
  }
  if (!file.exists(path) || file.info(path)$size != expected_bytes) {
    stop("Geometry store file has the wrong size.", call. = FALSE)
  }

  read_rows <- function(rows = NULL) {
    if (is.null(rows)) rows <- seq_len(dim[[1L]])
    rows <- as.integer(rows)
    out <- matrix(0, length(rows), dim[[2L]])
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    consecutive <- length(rows) < 2L || all(diff(rows) == 1L)
    for (column in seq_len(dim[[2L]])) {
      if (consecutive) {
        offset <- ((column - 1) * dim[[1L]] + rows[[1L]] - 1) * 8
        seek(connection, where = offset, origin = "start", rw = "read")
        out[, column] <- readBin(connection, "double", n = length(rows),
          size = 8, endian = "little")
      } else {
        for (i in seq_along(rows)) {
          offset <- ((column - 1) * dim[[1L]] + rows[[i]] - 1) * 8
          seek(connection, where = offset, origin = "start", rw = "read")
          out[i, column] <- readBin(connection, "double", n = 1L,
            size = 8, endian = "little")
        }
      }
    }
    out
  }

  structure(
    list(
      dim = dim,
      representation = "block_backed",
      manifest = list(schema_version = 1L, complete = TRUE, dim = dim,
        format = "packed-double-v1", expected_bytes = expected_bytes),
      path = normalizePath(path, mustWork = TRUE),
      read = read_rows
    ),
    class = "effect_geometry_store"
  )
}

.write_geometry_tile <- function(store, rows, coordinates, value) {
  if (!inherits(store, "effect_geometry_store") ||
      !identical(store$representation, "block_backed")) {
    stop("`store` must be a block-backed geometry store.", call. = FALSE)
  }
  rows <- as.integer(rows)
  coordinates <- as.integer(coordinates)
  if (length(rows) < 1L || length(coordinates) < 1L ||
      any(diff(rows) != 1L) || any(rows < 1L) || any(rows > store$dim[[1L]]) ||
      any(coordinates < 1L) || any(coordinates > store$dim[[2L]]) ||
      !is.matrix(value) || !identical(dim(value), c(length(rows), length(coordinates))) ||
      any(!is.finite(value))) {
    stop("Invalid geometry tile coordinates, shape, or values.", call. = FALSE)
  }
  connection <- file(store$path, open = "r+b")
  on.exit(close(connection), add = TRUE)
  for (j in seq_along(coordinates)) {
    offset <- ((coordinates[[j]] - 1) * store$dim[[1L]] + rows[[1L]] - 1) * 8
    seek(connection, where = offset, origin = "start", rw = "write")
    writeBin(as.double(value[, j]), connection, size = 8, endian = "little")
  }
  invisible(NULL)
}

.read_geometry_tile <- function(store, rows, coordinates) {
  if (!inherits(store, "effect_geometry_store") ||
      !identical(store$representation, "block_backed")) {
    stop("`store` must be a block-backed geometry store.", call. = FALSE)
  }
  rows <- as.integer(rows)
  coordinates <- as.integer(coordinates)
  if (length(rows) < 1L || length(coordinates) < 1L ||
      any(diff(rows) != 1L) || any(rows < 1L) || any(rows > store$dim[[1L]]) ||
      any(coordinates < 1L) || any(coordinates > store$dim[[2L]])) {
    stop("Invalid geometry tile coordinates.", call. = FALSE)
  }
  out <- matrix(0, length(rows), length(coordinates))
  connection <- file(store$path, open = "rb")
  on.exit(close(connection), add = TRUE)
  for (j in seq_along(coordinates)) {
    offset <- ((coordinates[[j]] - 1) * store$dim[[1L]] + rows[[1L]] - 1) * 8
    seek(connection, where = offset, origin = "start", rw = "read")
    out[, j] <- readBin(connection, "double", n = length(rows),
      size = 8, endian = "little")
  }
  if (any(!is.finite(out))) {
    stop("Geometry tile contains non-finite values.", call. = FALSE)
  }
  out
}
