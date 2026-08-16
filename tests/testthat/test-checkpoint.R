test_that("checkpoint round trip preserves payload and exact resume identity", {
  receipt <- receipt_fixture()
  path <- tempfile(pattern = "checkpoint-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  payload <- list(task = 3L, tile = matrix(c(1, 2, 3, 4), 2))

  crossform:::.write_checkpoint(payload, path, receipt, available_bytes = 1e8)
  got <- crossform:::.read_checkpoint(path, receipt)

  expect_identical(got, payload)
  expect_true(file.exists(file.path(path, "manifest.rds")))
  expect_true(file.exists(file.path(path, "payload.rds")))
  expect_false(any(grepl("\\.tmp-", list.files(dirname(path), all.files = TRUE))))
  if (.Platform$OS.type == "unix") {
    expect_equal(bitwAnd(as.integer(file.info(path)$mode), strtoi("077", base = 8L)), 0L)
  }
})

test_that("checkpoint disk preflight fails before staging", {
  receipt <- receipt_fixture()
  path <- tempfile(pattern = "checkpoint-")
  condition <- tryCatch(
    crossform:::.write_checkpoint(list(x = rnorm(100)), path, receipt,
      available_bytes = 1),
    effect_execution_error = identity
  )

  expect_s3_class(condition, "effect_execution_error")
  expect_match(condition$message, "Insufficient disk space")
  expect_false(file.exists(path))
})

test_that("checkpoint detects payload corruption", {
  receipt <- receipt_fixture()
  path <- tempfile(pattern = "checkpoint-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  crossform:::.write_checkpoint(list(x = 1:4), path, receipt,
    available_bytes = 1e8)
  connection <- file(file.path(path, "payload.rds"), open = "ab")
  writeBin(as.raw(1), connection)
  close(connection)

  condition <- tryCatch(
    crossform:::.read_checkpoint(path, receipt),
    effect_execution_error = identity
  )
  expect_match(condition$message, "checksum or size mismatch")
  expect_identical(condition$receipt, receipt)
})

test_that("resume rejects changed task, reduction, precision, and source identity", {
  receipt <- receipt_fixture()
  path <- tempfile(pattern = "checkpoint-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  crossform:::.write_checkpoint(42, path, receipt, available_bytes = 1e8)

  mutate_receipt <- function(field, value) {
    changed <- receipt
    changed[[field]] <- value
    changed
  }
  changed_source <- receipt
  changed_source$sources[[1]]$stable_revision <-
    paste0("sha256:", paste(rep("b", 64), collapse = ""))
  candidates <- list(
    mutate_receipt("task_partition_id", "other-tasks"),
    mutate_receipt("reduction_plan_id", "other-reduction"),
    mutate_receipt("precision", "single"),
    changed_source
  )
  for (candidate in candidates) {
    expect_error(crossform:::.read_checkpoint(path, candidate),
      "execution identity does not match")
  }
})

test_that("checkpoint rejects schema drift and permissive modes", {
  receipt <- receipt_fixture()
  path <- tempfile(pattern = "checkpoint-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  crossform:::.write_checkpoint(42, path, receipt, available_bytes = 1e8)
  manifest_path <- file.path(path, "manifest.rds")
  manifest <- readRDS(manifest_path)
  manifest$schema_version <- 99L
  saveRDS(manifest, manifest_path)
  expect_error(crossform:::.read_checkpoint(path, receipt), "schema version")

  manifest$schema_version <- 1L
  saveRDS(manifest, manifest_path)
  Sys.chmod(manifest_path, "0644", use_umask = FALSE)
  if (.Platform$OS.type == "unix") {
    expect_error(crossform:::.read_checkpoint(path, receipt), "permissions")
  }
})
