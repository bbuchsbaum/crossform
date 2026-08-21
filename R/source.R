# Reopenable response-source descriptors -----------------------------------

.descriptor_value_is_serializable <- function(value) {
  if (is.null(value) || is.atomic(value)) return(TRUE)
  if (!is.list(value) || is.object(value)) return(FALSE)
  all(vapply(value, .descriptor_value_is_serializable, logical(1)))
}

.source_descriptor <- function(kind, dim, access, stable_revision, spec) {
  if (!.is_string(kind)) {
    .input_error("Source descriptor `kind` must be one nonempty identifier.")
  }
  if (!.is_finite_numeric(dim) || length(dim) != 2L || anyNA(dim) ||
      any(dim < 1L) || any(dim %% 1 != 0)) {
    .input_error("Source descriptor dimensions must be two positive integers.")
  }
  access <- match.arg(access, c("coordinator", "reopenable", "shared"))
  if (!.strong_sha256(stable_revision)) {
    .input_error("Source descriptors require one strong sha256 revision.")
  }
  if (!is.list(spec) || !.descriptor_value_is_serializable(spec)) {
    .input_error(
      "Source descriptor specifications must contain only serializable values."
    )
  }
  value <- structure(list(
    kind = kind,
    dim = as.integer(dim),
    access = access,
    stable_revision = stable_revision,
    spec = spec
  ), class = "effect_source_descriptor")
  .validate_source_descriptor(value)
}

.validate_source_descriptor <- function(value) {
  expected <- c("kind", "dim", "access", "stable_revision", "spec")
  if (!.sealed_fields(value, "effect_source_descriptor", expected)) {
    .input_error("Source descriptor fields are missing or noncanonical.")
  }
  if (!.is_string(value$kind) || !.is_finite_numeric(value$dim) ||
      length(value$dim) != 2L || anyNA(value$dim) || any(value$dim < 1L) ||
      any(value$dim %% 1 != 0) ||
      !identical(value$access, match.arg(value$access, c("coordinator", "reopenable", "shared"))) ||
      !.strong_sha256(value$stable_revision) || !is.list(value$spec) ||
      !.descriptor_value_is_serializable(value$spec)) {
    .input_error("Source descriptor metadata are invalid.")
  }
  structure(as.list(value), class = "effect_source_descriptor")
}

#' Describe a read-only column-major matrix file
#'
#' The file must contain only consecutive IEEE-754 doubles in R column-major
#' order. The constructor records an absolute path and verifies a strong content
#' revision. It never opens a persistent handle.
#'
#' @param path Existing binary matrix file.
#' @param dim Two positive integers: observations by neural features.
#' @param offset_bytes Nonnegative byte offset before the matrix payload.
#' @param endian Byte order used by the file.
#' @param stable_revision Optional expected `sha256:` content revision. When
#'   omitted it is computed from the file.
#' @return An `effect_source_descriptor`: a list with `$kind`
#'   (`"file_matrix"`), the integer `$dim`, `$access` (`"reopenable"`), the
#'   `$stable_revision` content hash, and a `$spec` holding the absolute
#'   `path`, `offset_bytes`, and `endian`. Treat it as immutable.
#' @family relation planning and fitting
#' @seealso [source_capabilities()] for the capability value it implies, and
#'   [observations()] or [relation()], which accept descriptors in place of
#'   in-memory matrices.
#' @examples
#' # A 4-observation by 3-feature matrix written as raw column-major doubles.
#' path <- tempfile(fileext = ".bin")
#' writeBin(as.vector(matrix(as.double(1:12), 4L, 3L)), path, size = 8L)
#'
#' descriptor <- file_matrix_source(path, dim = c(4L, 3L))
#' descriptor$dim
#' descriptor$access
#'
#' # The content hash is recorded now and rechecked whenever the file is
#' # reopened, so a silently edited source is caught rather than used.
#' substr(descriptor$stable_revision, 1, 24)
#'
#' # Declaring a revision that no longer matches the bytes is an error.
#' try(file_matrix_source(
#'   path, dim = c(4L, 3L), stable_revision = paste0("sha256:", strrep("0", 64))
#' ))
#' unlink(path)
#' @export
file_matrix_source <- function(path, dim, offset_bytes = 0,
                               endian = .Platform$endian,
                               stable_revision = NULL) {
  if (!.is_string(path)) {
    .input_error("`path` must identify one existing matrix file.")
  }
  path <- normalizePath(path, mustWork = TRUE)
  if (!.is_finite_numeric(dim) || length(dim) != 2L || anyNA(dim) ||
      any(dim < 1L) || any(dim %% 1 != 0)) {
    .input_error("`dim` must contain two positive integers.")
  }
  dim <- as.integer(dim)
  .check_count(
    offset_bytes, "offset_bytes", what = "one nonnegative integer", min = 0L
  )
  offset_bytes <- as.double(offset_bytes)
  endian <- match.arg(endian, c("little", "big"))
  expected_size <- offset_bytes + prod(as.double(dim)) * 8
  observed_size <- file.info(path)$size
  if (is.na(observed_size) || observed_size != expected_size) {
    .input_error(sprintf(
      "Matrix file size is %.0f bytes; descriptor requires %.0f bytes.",
      observed_size, expected_size
    ))
  }
  observed_revision <- .sha256_file(path)
  if (!is.null(stable_revision) &&
      (!.strong_sha256(stable_revision) ||
       !identical(tolower(stable_revision), tolower(observed_revision)))) {
    .contract_error("The matrix file does not match `stable_revision`.")
  }
  .source_descriptor(
    kind = "file_matrix",
    dim = dim,
    access = "reopenable",
    stable_revision = observed_revision,
    spec = list(path = path, offset_bytes = offset_bytes, endian = endian)
  )
}

# Adapter-neutral shared descriptors are intentionally internal until a shard
# admission benchmark establishes a public staging API.
.shared_source_descriptor <- function(backend, token, dim, stable_revision) {
  if (!.is_string(backend)) {
    .input_error("Shared source `backend` must be one nonempty identifier.")
  }
  .source_descriptor(
    kind = "shared_matrix",
    dim = dim,
    access = "shared",
    stable_revision = stable_revision,
    spec = list(backend = backend, token = token)
  )
}

.validate_source_features <- function(features, n_features) {
  if (!.is_finite_numeric(features) || length(features) < 1L ||
      anyNA(features) || any(features %% 1 != 0) || any(features < 1L) ||
      any(features > n_features) || anyDuplicated(features)) {
    .input_error("`features` must be unique valid neural feature indices.")
  }
  as.integer(features)
}

.new_source_handle <- function(descriptor, read, close, owns_handle,
                               owns_backing = FALSE) {
  if (!is.function(read) || !is.function(close) || !.is_flag(owns_handle) ||
      !.is_flag(owns_backing)) {
    .input_error("Source handle implementation is invalid.")
  }
  structure(list(
    descriptor = .validate_source_descriptor(descriptor),
    read = read,
    close = close,
    owns_handle = owns_handle,
    owns_backing = owns_backing
  ), class = "effect_source_handle")
}

.open_file_matrix_descriptor <- function(descriptor, expected_revision) {
  path <- descriptor$spec$path
  if (!file.exists(path)) {
    .input_error("Reopenable matrix source no longer exists.")
  }
  observed <- .sha256_file(path)
  if (!identical(tolower(observed), tolower(expected_revision)) ||
      !identical(tolower(observed), tolower(descriptor$stable_revision))) {
    .contract_error("Reopenable matrix source has a stale content revision.")
  }
  connection <- file(path, open = "rb")
  state <- new.env(parent = emptyenv())
  state$closed <- FALSE
  close_handle <- function() {
    if (!state$closed) {
      close(connection)
      state$closed <- TRUE
    }
    invisible(NULL)
  }
  read_block <- function(features) {
    if (state$closed) .input_error("Source handle is closed.")
    features <- .validate_source_features(features, descriptor$dim[[2L]])
    rows <- descriptor$dim[[1L]]
    value <- matrix(NA_real_, rows, length(features))
    for (column in seq_along(features)) {
      position <- descriptor$spec$offset_bytes +
        (as.double(features[[column]]) - 1) * rows * 8
      seek(connection, where = position, origin = "start", rw = "read")
      values <- readBin(
        connection, what = double(), n = rows, size = 8,
        endian = descriptor$spec$endian
      )
      if (length(values) != rows) {
        .input_error(
          "Matrix source ended before the requested feature was read."
        )
      }
      value[, column] <- values
    }
    value
  }
  .new_source_handle(
    descriptor, read_block, close_handle,
    owns_handle = TRUE, owns_backing = FALSE
  )
}

.open_source_descriptor <- function(descriptor,
                                    expected_revision = descriptor$stable_revision,
                                    shared_opener = NULL) {
  descriptor <- .validate_source_descriptor(descriptor)
  if (!.strong_sha256(expected_revision)) {
    .input_error(
      "Opening a source requires one strong expected sha256 revision."
    )
  }
  if (!identical(tolower(expected_revision),
      tolower(descriptor$stable_revision))) {
    .contract_error(
      "Source descriptor revision does not match the expected revision."
    )
  }
  if (identical(descriptor$kind, "file_matrix")) {
    return(.open_file_matrix_descriptor(descriptor, expected_revision))
  }
  if (identical(descriptor$access, "shared")) {
    if (!is.function(shared_opener)) {
      .input_error(
        "Shared source descriptors require an explicit backend opener."
      )
    }
    opened <- shared_opener(descriptor)
    if (!is.list(opened) || !is.function(opened$read) ||
        !is.function(opened$close) ||
        !.strong_sha256(opened$stable_revision) ||
        !identical(tolower(opened$stable_revision),
          tolower(expected_revision))) {
      if (is.list(opened) && is.function(opened$close)) {
        try(opened$close(), silent = TRUE)
      }
      .input_error("Shared source opener returned an invalid or stale handle.")
    }
    return(.new_source_handle(
      descriptor, opened$read, opened$close,
      owns_handle = TRUE, owns_backing = FALSE
    ))
  }
  .input_error("Coordinator-only source descriptors cannot be reopened.")
}

.close_source_handle <- function(handle) {
  if (!inherits(handle, "effect_source_handle") || !is.list(handle) ||
      !is.function(handle$close)) {
    .input_error("`handle` must be an open source handle.")
  }
  if (isTRUE(handle$owns_handle)) handle$close()
  invisible(NULL)
}

.with_source_descriptor <- function(descriptor, code, shared_opener = NULL) {
  if (!is.function(code)) .input_error("`code` must be a function.")
  handle <- .open_source_descriptor(descriptor,
    shared_opener = shared_opener)
  on.exit(.close_source_handle(handle), add = TRUE)
  code(handle)
}

.source_descriptor_key <- function(descriptor) {
  descriptor <- .validate_source_descriptor(descriptor)
  .sha256_signature(descriptor)
}
