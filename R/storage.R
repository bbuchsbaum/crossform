# File-backed packed geometry storage ---------------------------------------
#
# A geometry store is one dense column-major `dim[1] x dim[2]` block of
# little-endian float64 with no header and no sidecar. Every access goes
# through `.tile_io()` in R/primitives.R, which is also what the
# variable-block measurement store in R/measurement-storage.R transfers
# through; the two stores share that primitive and nothing else, because
# their offset arithmetic and their completion bookkeeping are different
# formats on disk (see the note on `.tile_io()`).

# Element offset of one entry of the dense column-major payload. Vectorized
# over both arguments, so one call serves a contiguous run per column and a
# scattered element-at-a-time read alike.
.geometry_tile_offsets <- function(dim, rows, columns) {
  (columns - 1) * dim[[1L]] + rows - 1
}

# Shared prologue for both tile operations: the store must be block backed and
# the tile must name a contiguous row run inside the store's shape. The write
# path additionally requires a finite matrix of exactly that shape, and
# reports the combined failure in one message.
.geometry_tile_extent <- function(store, rows, coordinates, value = NULL) {
  if (!inherits(store, "effect_geometry_store") ||
      !identical(store$representation, "block_backed")) {
    .input_error("`store` must be a block-backed geometry store.")
  }
  rows <- as.integer(rows)
  coordinates <- as.integer(coordinates)
  valid <- length(rows) >= 1L && length(coordinates) >= 1L &&
    all(diff(rows) == 1L) && all(rows >= 1L) &&
    all(rows <= store$dim[[1L]]) && all(coordinates >= 1L) &&
    all(coordinates <= store$dim[[2L]])
  if (is.null(value)) {
    if (!valid) .input_error("Invalid geometry tile coordinates.")
  } else {
    valid <- valid && is.matrix(value) &&
      identical(dim(value), c(length(rows), length(coordinates))) &&
      all(is.finite(value))
    if (!valid) {
      .input_error("Invalid geometry tile coordinates, shape, or values.")
    }
  }
  list(rows = rows, coordinates = coordinates)
}

.file_geometry_store <- function(path, dim, create = FALSE,
                                 codec = "symmetric_packed") {
  .check_string(path, "path", what = "one nonempty file path")
  if (!.is_finite_numeric(dim) || length(dim) != 2L || any(dim < 1) ||
      any(dim %% 1 != 0)) {
    .input_error("`dim` must contain two positive integers.")
  }
  dim <- as.integer(dim)
  format <- .effect_form_codec_format(codec)
  expected_bytes <- prod(as.double(dim)) * 8
  if (create) {
    if (file.exists(path)) .input_error("Refusing to overwrite an existing geometry store.")
    .tile_zero_fill(path, prod(as.double(dim)))
  }
  if (!file.exists(path) || file.info(path)$size != expected_bytes) {
    .input_error("Geometry store file has the wrong size.")
  }

  read_rows <- function(rows = NULL) {
    if (is.null(rows)) rows <- seq_len(dim[[1L]])
    rows <- as.integer(rows)
    out <- matrix(0, length(rows), dim[[2L]])
    columns <- seq_len(dim[[2L]])
    consecutive <- length(rows) < 2L || all(diff(rows) == 1L)
    # One run per column when the rows are contiguous, one run per element
    # otherwise; either way the payload is opened exactly once.
    if (consecutive) {
      runs <- .tile_io(path, .geometry_tile_offsets(dim, rows[[1L]], columns),
        n = length(rows))
      for (k in seq_along(columns)) out[, columns[[k]]] <- runs[[k]]
    } else {
      grid <- expand.grid(row = rows, column = columns)
      runs <- .tile_io(
        path, .geometry_tile_offsets(dim, grid$row, grid$column), n = 1L
      )
      out[] <- unlist(runs, use.names = FALSE)
    }
    out
  }

  structure(
    list(
      dim = dim,
      representation = "block_backed",
      manifest = list(schema_version = 1L, complete = TRUE, dim = dim,
        format = format, expected_bytes = expected_bytes),
      path = normalizePath(path, mustWork = TRUE),
      read = read_rows
    ),
    class = "effect_geometry_store"
  )
}

.write_geometry_tile <- function(store, rows, coordinates, value) {
  extent <- .geometry_tile_extent(store, rows, coordinates, value)
  .tile_io(
    store$path,
    .geometry_tile_offsets(store$dim, extent$rows[[1L]], extent$coordinates),
    mode = "write",
    values = lapply(seq_along(extent$coordinates), function(j) value[, j])
  )
}

.read_geometry_tile <- function(store, rows, coordinates) {
  extent <- .geometry_tile_extent(store, rows, coordinates)
  runs <- .tile_io(
    store$path,
    .geometry_tile_offsets(store$dim, extent$rows[[1L]], extent$coordinates),
    n = length(extent$rows)
  )
  out <- matrix(0, length(extent$rows), length(extent$coordinates))
  for (j in seq_along(extent$coordinates)) out[, j] <- runs[[j]]
  if (any(!is.finite(out))) {
    .input_error("Geometry tile contains non-finite values.")
  }
  out
}
