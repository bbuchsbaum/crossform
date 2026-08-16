# Reopenable response-source descriptors -----------------------------------

.strong_sha256 <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^sha256:[[:xdigit:]]{64}$", value)
}

.descriptor_value_is_serializable <- function(value) {
  if (is.null(value) || is.atomic(value)) return(TRUE)
  if (!is.list(value) || is.object(value)) return(FALSE)
  all(vapply(value, .descriptor_value_is_serializable, logical(1)))
}

.source_descriptor <- function(kind, dim, access, stable_revision, spec) {
  if (!is.character(kind) || length(kind) != 1L || is.na(kind) || !nzchar(kind)) {
    stop("Source descriptor `kind` must be one nonempty identifier.",
      call. = FALSE)
  }
  if (!is.numeric(dim) || length(dim) != 2L || anyNA(dim) ||
      any(!is.finite(dim)) || any(dim < 1L) || any(dim %% 1 != 0)) {
    stop("Source descriptor dimensions must be two positive integers.",
      call. = FALSE)
  }
  access <- match.arg(access, c("coordinator", "reopenable", "shared"))
  if (!.strong_sha256(stable_revision)) {
    stop("Source descriptors require one strong sha256 revision.",
      call. = FALSE)
  }
  if (!is.list(spec) || !.descriptor_value_is_serializable(spec)) {
    stop("Source descriptor specifications must contain only serializable values.",
      call. = FALSE)
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
  if (!inherits(value, "effect_source_descriptor") || !is.list(value) ||
      !identical(names(value), expected)) {
    stop("Source descriptor fields are missing or noncanonical.", call. = FALSE)
  }
  if (!is.character(value$kind) || length(value$kind) != 1L ||
      is.na(value$kind) || !nzchar(value$kind) ||
      !is.numeric(value$dim) || length(value$dim) != 2L ||
      anyNA(value$dim) || any(!is.finite(value$dim)) ||
      any(value$dim < 1L) || any(value$dim %% 1 != 0) ||
      !identical(value$access, match.arg(value$access,
        c("coordinator", "reopenable", "shared"))) ||
      !.strong_sha256(value$stable_revision) || !is.list(value$spec) ||
      !.descriptor_value_is_serializable(value$spec)) {
    stop("Source descriptor metadata are invalid.", call. = FALSE)
  }
  structure(as.list(value), class = "effect_source_descriptor")
}

.file_sha256 <- function(path) {
  paste0("sha256:", digest::digest(file = path, algo = "sha256"))
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
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must identify one existing matrix file.", call. = FALSE)
  }
  path <- normalizePath(path, mustWork = TRUE)
  if (!is.numeric(dim) || length(dim) != 2L || anyNA(dim) ||
      any(!is.finite(dim)) || any(dim < 1L) || any(dim %% 1 != 0)) {
    stop("`dim` must contain two positive integers.", call. = FALSE)
  }
  dim <- as.integer(dim)
  if (!is.numeric(offset_bytes) || length(offset_bytes) != 1L ||
      is.na(offset_bytes) || !is.finite(offset_bytes) || offset_bytes < 0 ||
      offset_bytes %% 1 != 0) {
    stop("`offset_bytes` must be one nonnegative integer.", call. = FALSE)
  }
  offset_bytes <- as.double(offset_bytes)
  endian <- match.arg(endian, c("little", "big"))
  expected_size <- offset_bytes + prod(as.double(dim)) * 8
  observed_size <- file.info(path)$size
  if (is.na(observed_size) || observed_size != expected_size) {
    stop(sprintf(
      "Matrix file size is %.0f bytes; descriptor requires %.0f bytes.",
      observed_size, expected_size
    ), call. = FALSE)
  }
  observed_revision <- .file_sha256(path)
  if (!is.null(stable_revision) &&
      (!.strong_sha256(stable_revision) ||
       !identical(tolower(stable_revision), tolower(observed_revision)))) {
    stop("The matrix file does not match `stable_revision`.", call. = FALSE)
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
  if (!is.character(backend) || length(backend) != 1L || is.na(backend) ||
      !nzchar(backend)) {
    stop("Shared source `backend` must be one nonempty identifier.",
      call. = FALSE)
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
  if (!is.numeric(features) || length(features) < 1L || anyNA(features) ||
      any(!is.finite(features)) || any(features %% 1 != 0) ||
      any(features < 1L) || any(features > n_features) ||
      anyDuplicated(features)) {
    stop("`features` must be unique valid neural feature indices.",
      call. = FALSE)
  }
  as.integer(features)
}

.new_source_handle <- function(descriptor, read, close, owns_handle,
                               owns_backing = FALSE) {
  if (!is.function(read) || !is.function(close) ||
      !is.logical(owns_handle) || length(owns_handle) != 1L ||
      is.na(owns_handle) || !is.logical(owns_backing) ||
      length(owns_backing) != 1L || is.na(owns_backing)) {
    stop("Source handle implementation is invalid.", call. = FALSE)
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
    stop("Reopenable matrix source no longer exists.", call. = FALSE)
  }
  observed <- .file_sha256(path)
  if (!identical(tolower(observed), tolower(expected_revision)) ||
      !identical(tolower(observed), tolower(descriptor$stable_revision))) {
    stop("Reopenable matrix source has a stale content revision.",
      call. = FALSE)
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
    if (state$closed) stop("Source handle is closed.", call. = FALSE)
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
        stop("Matrix source ended before the requested feature was read.",
          call. = FALSE)
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
    stop("Opening a source requires one strong expected sha256 revision.",
      call. = FALSE)
  }
  if (!identical(tolower(expected_revision),
      tolower(descriptor$stable_revision))) {
    stop("Source descriptor revision does not match the expected revision.",
      call. = FALSE)
  }
  if (identical(descriptor$kind, "file_matrix")) {
    return(.open_file_matrix_descriptor(descriptor, expected_revision))
  }
  if (identical(descriptor$access, "shared")) {
    if (!is.function(shared_opener)) {
      stop("Shared source descriptors require an explicit backend opener.",
        call. = FALSE)
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
      stop("Shared source opener returned an invalid or stale handle.",
        call. = FALSE)
    }
    return(.new_source_handle(
      descriptor, opened$read, opened$close,
      owns_handle = TRUE, owns_backing = FALSE
    ))
  }
  stop("Coordinator-only source descriptors cannot be reopened.",
    call. = FALSE)
}

.close_source_handle <- function(handle) {
  if (!inherits(handle, "effect_source_handle") || !is.list(handle) ||
      !is.function(handle$close)) {
    stop("`handle` must be an open source handle.", call. = FALSE)
  }
  if (isTRUE(handle$owns_handle)) handle$close()
  invisible(NULL)
}

.with_source_descriptor <- function(descriptor, code, shared_opener = NULL) {
  if (!is.function(code)) stop("`code` must be a function.", call. = FALSE)
  handle <- .open_source_descriptor(descriptor,
    shared_opener = shared_opener)
  on.exit(.close_source_handle(handle), add = TRUE)
  code(handle)
}

.source_descriptor_key <- function(descriptor) {
  descriptor <- .validate_source_descriptor(descriptor)
  paste0("sha256:", digest::digest(descriptor, algo = "sha256",
    serialize = TRUE, serializeVersion = 2L))
}

.open_relation_source_session <- function(relation,
                                          open_descriptor = .open_source_descriptor,
                                          shared_opener = NULL,
                                          validate = TRUE) {
  if (isTRUE(validate)) .validate_relation(relation)
  if (!is.function(open_descriptor)) {
    stop("`open_descriptor` must be a function.", call. = FALSE)
  }
  handles <- list()
  partition_keys <- stats::setNames(vector("list", length(relation$partitions)),
    relation$partitions)
  facts <- new.env(parent = emptyenv())
  facts$closed <- FALSE
  facts$close_attempts <- 0L
  facts$read_count <- stats::setNames(integer(length(relation$partitions)),
    relation$partitions)
  facts$bytes_read <- stats::setNames(numeric(length(relation$partitions)),
    relation$partitions)

  close_opened <- function() {
    if (!facts$closed) {
      facts$close_attempts <- facts$close_attempts + 1L
      failures <- character()
      for (handle in rev(handles)) {
        failure <- tryCatch({
          .close_source_handle(handle)
          NULL
        }, error = function(error) conditionMessage(error))
        if (!is.null(failure)) failures <- c(failures, failure)
      }
      facts$closed <- TRUE
      if (length(failures)) {
        stop(paste("Source-session cleanup failed:",
          paste(failures, collapse = "; ")), call. = FALSE)
      }
    }
    invisible(NULL)
  }

  for (partition in relation$partitions) {
    descriptor <- relation$sources[[partition]]$descriptor
    if (is.null(descriptor) || identical(descriptor$access, "coordinator")) {
      partition_keys[partition] <- list(NULL)
      next
    }
    key <- .source_descriptor_key(descriptor)
    partition_keys[[partition]] <- key
    if (is.null(handles[[key]])) {
      handle <- tryCatch(
        open_descriptor(descriptor,
          expected_revision = relation$capabilities[[partition]]$stable_revision,
          shared_opener = shared_opener),
        error = function(error) {
          cleanup_error <- tryCatch({
            close_opened()
            NULL
          }, error = function(cleanup) conditionMessage(cleanup))
          if (!is.null(cleanup_error)) {
            stop(paste(conditionMessage(error), cleanup_error, sep = "; "),
              call. = FALSE)
          }
          stop(error)
        }
      )
      handles[[key]] <- handle
    }
  }

  read <- function(partition, features) {
    if (facts$closed) stop("Source session is closed.", call. = FALSE)
    if (!partition %in% relation$partitions) {
      stop("Source session partition is invalid.", call. = FALSE)
    }
    features <- .validate_source_features(features, relation$n_features)
    key <- partition_keys[[partition]]
    value <- if (is.null(key)) {
      relation$sources[[partition]]$read(features)
    } else {
      handles[[key]]$read(features)
    }
    facts$read_count[[partition]] <- facts$read_count[[partition]] + 1L
    facts$bytes_read[[partition]] <- facts$bytes_read[[partition]] +
      prod(as.double(dim(value))) * 8
    value
  }

  summary <- function() {
    descriptors <- lapply(relation$sources, `[[`, "descriptor")
    access <- vapply(descriptors, function(descriptor) {
      if (is.null(descriptor)) "opaque_coordinator" else descriptor$access
    }, character(1))
    list(
      access_mode = stats::setNames(access, relation$partitions),
      distinct_owned_handles = length(handles),
      read_count = facts$read_count,
      bytes_read = facts$bytes_read,
      closed = facts$closed,
      close_attempts = facts$close_attempts
    )
  }

  structure(list(read = read, close = close_opened, summary = summary),
    class = "effect_source_session")
}

.close_source_session <- function(session) {
  if (!inherits(session, "effect_source_session") || !is.function(session$close)) {
    stop("`session` must be a crossform source session.", call. = FALSE)
  }
  session$close()
}

.open_effect_task_source_session <- function(task,
                                             open_descriptor = .open_source_descriptor,
                                             shared_opener = NULL,
                                             validate = TRUE) {
  if (isTRUE(validate)) .validate_compiled_effect_task(task)
  if (isTRUE(task$same_relation)) {
    session <- .open_relation_source_session(
      task$left_relation,
      open_descriptor = open_descriptor,
      shared_opener = shared_opener,
      validate = validate
    )
    return(structure(
      list(
        read = function(side, partition, features) {
          if (!side %in% c("left", "right")) {
            stop("Task source side must be `left` or `right`.", call. = FALSE)
          }
          session$read(partition, features)
        },
        close = session$close,
        summary = session$summary
      ),
      class = "effect_task_source_session"
    ))
  }
  left <- .open_relation_source_session(task$left_relation,
    open_descriptor = open_descriptor, shared_opener = shared_opener,
    validate = validate)
  right <- tryCatch(
    .open_relation_source_session(task$right_relation,
      open_descriptor = open_descriptor, shared_opener = shared_opener,
      validate = validate),
    error = function(error) {
      left$close()
      stop(error)
    }
  )
  closed <- FALSE
  close_both <- function() {
    if (!closed) {
      on.exit(right$close(), add = TRUE)
      left$close()
      closed <<- TRUE
    }
    invisible(NULL)
  }
  structure(
    list(
      read = function(side, partition, features) {
        if (identical(side, "left")) left$read(partition, features) else if (
          identical(side, "right")) right$read(partition, features) else
          stop("Task source side must be `left` or `right`.", call. = FALSE)
      },
      close = close_both,
      summary = function() list(left = left$summary(), right = right$summary())
    ),
    class = "effect_task_source_session"
  )
}
