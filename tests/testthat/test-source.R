write_matrix_file <- function(value, path, endian = .Platform$endian) {
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(as.double(value), connection, size = 8, endian = endian)
  invisible(path)
}

test_that("ordinary matrices remain coordinator-local strong-revision sources", {
  value <- matrix(seq_len(20), 4, 5)
  rel <- relation(list(run1 = value), effects = paste0("e", 1:4))
  descriptor <- crossform:::.relation_source_descriptors(rel)$run1

  expect_s3_class(descriptor, "effect_source_descriptor")
  expect_identical(descriptor$kind, "matrix")
  expect_identical(descriptor$access, "coordinator")
  expect_false(rel$capabilities$run1$reopenable)
  expect_error(
    crossform:::.relation_source_descriptors(rel, require_reopenable = TRUE),
    "requires reopenable"
  )
  expect_error(
    crossform:::.open_source_descriptor(descriptor),
    "cannot be reopened"
  )
})

test_that("file matrix descriptors reopen to the same relation blocks", {
  value <- matrix(seq_len(30) / 7, 5, 6)
  path <- tempfile("crossform-matrix-", fileext = ".bin")
  on.exit(unlink(path), add = TRUE)
  write_matrix_file(value, path)
  descriptor <- file_matrix_source(path, dim(value))
  serialized <- unserialize(serialize(descriptor, NULL, version = 3))

  expect_identical(serialized, descriptor)
  expect_identical(descriptor$access, "reopenable")
  expect_match(descriptor$stable_revision, "^sha256:[[:xdigit:]]{64}$")

  rel <- relation(list(run1 = descriptor), effects = paste0("e", 1:5))
  expect_true(rel$capabilities$run1$reopenable)
  expect_equal(
    unname(relation_block(rel, "run1", c(6, 2, 4))),
    value[, c(6, 2, 4), drop = FALSE],
    tolerance = 0
  )

  reopened <- crossform:::.relation_source_descriptors(
    rel, require_reopenable = TRUE
  )
  expect_identical(reopened$run1, descriptor)
})

test_that("opening fails closed when file content no longer matches revision", {
  value <- matrix(seq_len(12), 3, 4)
  path <- tempfile("crossform-stale-", fileext = ".bin")
  on.exit(unlink(path), add = TRUE)
  write_matrix_file(value, path)
  descriptor <- file_matrix_source(path, dim(value))
  value[1, 1] <- -999
  write_matrix_file(value, path)

  expect_error(
    crossform:::.open_source_descriptor(descriptor),
    "stale content revision"
  )
})

test_that("file handles close idempotently without owning the backing file", {
  value <- matrix(rnorm(20), 4, 5)
  path <- tempfile("crossform-owned-handle-", fileext = ".bin")
  on.exit(unlink(path), add = TRUE)
  write_matrix_file(value, path)
  descriptor <- file_matrix_source(path, dim(value))
  handle <- crossform:::.open_source_descriptor(descriptor)

  expect_true(handle$owns_handle)
  expect_false(handle$owns_backing)
  expect_equal(handle$read(c(1, 5)), value[, c(1, 5)], tolerance = 0)
  crossform:::.close_source_handle(handle)
  crossform:::.close_source_handle(handle)
  expect_true(file.exists(path))
  expect_error(handle$read(1), "closed")
})

test_that("shared descriptors delegate attachment but never backing ownership", {
  revision <- paste0("sha256:", paste(rep("a", 64), collapse = ""))
  descriptor <- crossform:::.shared_source_descriptor(
    backend = "test", token = list(id = "segment-1"), dim = c(2, 3),
    stable_revision = revision
  )
  closed <- 0L
  opener <- function(value) {
    expect_identical(value$spec$token$id, "segment-1")
    list(
      read = function(features) matrix(features, 2, length(features), byrow = TRUE),
      close = function() closed <<- closed + 1L,
      stable_revision = revision
    )
  }
  handle <- crossform:::.open_source_descriptor(
    descriptor, shared_opener = opener
  )

  expect_false(handle$owns_backing)
  expect_equal(handle$read(c(1, 3)), matrix(c(1, 3, 1, 3), 2, byrow = TRUE))
  crossform:::.close_source_handle(handle)
  expect_identical(closed, 1L)
})

test_that("invalid shared attachments clean up before failing", {
  revision <- paste0("sha256:", paste(rep("b", 64), collapse = ""))
  stale_revision <- paste0("sha256:", paste(rep("c", 64), collapse = ""))
  descriptor <- crossform:::.shared_source_descriptor(
    "test", "segment", c(2, 2), revision
  )
  closed <- 0L
  opener <- function(value) list(
    read = function(features) matrix(0, 2, length(features)),
    close = function() closed <<- closed + 1L,
    stable_revision = stale_revision
  )

  expect_error(
    crossform:::.open_source_descriptor(descriptor, shared_opener = opener),
    "invalid or stale"
  )
  expect_identical(closed, 1L)
})

test_that("reopenable claims without descriptors are rejected at construction", {
  revision <- paste0("sha256:", paste(rep("d", 64), collapse = ""))
  source <- function(features) matrix(1, 2, length(features))
  expect_error(
    relation(
      list(run1 = source), effects = c("a", "b"),
      source_dims = list(c(2, 3)),
      capabilities = source_capabilities(
        block_read = TRUE, reopenable = TRUE,
        stable_revision = revision
      )
    ),
    "require a reopenable descriptor"
  )
})

test_that("descriptor specifications reject closures and environments", {
  revision <- paste0("sha256:", paste(rep("e", 64), collapse = ""))
  expect_error(
    crossform:::.shared_source_descriptor(
      "test", list(callback = function() NULL), c(2, 2), revision
    ),
    "serializable values"
  )
  expect_error(
    crossform:::.shared_source_descriptor(
      "test", new.env(), c(2, 2), revision
    ),
    "serializable values"
  )
})

test_that("execution sessions deduplicate descriptors and close exactly once", {
  revision <- paste0("sha256:", paste(rep("f", 64), collapse = ""))
  descriptor <- crossform:::.shared_source_descriptor(
    "test", "same-segment", c(2, 3), revision
  )
  rel <- relation(list(run1 = descriptor, run2 = descriptor),
    effects = c("a", "b"))
  opens <- 0L
  closes <- 0L
  opener <- function(value, expected_revision, shared_opener) {
    opens <<- opens + 1L
    crossform:::.new_source_handle(
      value,
      read = function(features) matrix(features, 2, length(features), byrow = TRUE),
      close = function() closes <<- closes + 1L,
      owns_handle = TRUE
    )
  }
  session <- crossform:::.open_relation_source_session(rel,
    open_descriptor = opener)

  expect_identical(opens, 1L)
  expect_equal(session$read("run1", c(1, 3)),
    matrix(c(1, 3, 1, 3), 2, byrow = TRUE))
  expect_equal(session$read("run2", 2), matrix(2, 2, 1))
  expect_identical(session$summary()$distinct_owned_handles, 1L)
  expect_identical(unname(session$summary()$read_count), c(1L, 1L))
  expect_identical(unname(session$summary()$bytes_read), c(32, 16))
  crossform:::.close_source_session(session)
  crossform:::.close_source_session(session)
  expect_identical(closes, 1L)
  expect_identical(session$summary()$close_attempts, 1L)
})

test_that("partially opened sessions clean up when later admission fails", {
  first_revision <- paste0("sha256:", paste(rep("1", 64), collapse = ""))
  second_revision <- paste0("sha256:", paste(rep("2", 64), collapse = ""))
  first <- crossform:::.shared_source_descriptor("test", "one", c(2, 3),
    first_revision)
  second <- crossform:::.shared_source_descriptor("test", "two", c(2, 3),
    second_revision)
  rel <- relation(list(one = first, two = second), effects = c("a", "b"))
  closes <- 0L
  opener <- function(value, expected_revision, shared_opener) {
    if (identical(value$spec$token, "two")) stop("admission failed")
    crossform:::.new_source_handle(value,
      read = function(features) matrix(0, 2, length(features)),
      close = function() closes <<- closes + 1L, owns_handle = TRUE)
  }

  expect_error(crossform:::.open_relation_source_session(rel,
    open_descriptor = opener), "admission failed")
  expect_identical(closes, 1L)
})

test_that("source-session close failures are reported rather than swallowed", {
  revision <- paste0("sha256:", paste(rep("3", 64), collapse = ""))
  descriptor <- crossform:::.shared_source_descriptor(
    "test", "failing-close", c(2, 2), revision
  )
  rel <- relation(list(run = descriptor), effects = c("a", "b"))
  opener <- function(value, expected_revision, shared_opener) {
    crossform:::.new_source_handle(value,
      read = function(features) matrix(0, 2, length(features)),
      close = function() stop("detach failed"), owns_handle = TRUE)
  }
  session <- crossform:::.open_relation_source_session(rel,
    open_descriptor = opener)

  expect_error(crossform:::.close_source_session(session),
    "cleanup failed.*detach failed")
  expect_true(session$summary()$closed)
})
