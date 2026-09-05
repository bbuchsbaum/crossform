test_that("[T45 T52 T61] fresh-checkout calibration restores frozen counts without re-piloting", {
  root <- .certification_source_root()
  skip_if(is.na(root), "Calibration tooling is available only in a source checkout")
  x <- certified_artifact("predictive-geometry-validation.rds", "predictive-geometry/certify.R")
  output <- tempfile("frozen-counts-")
  withr::defer(unlink(output, recursive = TRUE))
  run <- function() processx::run(file.path(R.home("bin"), "Rscript"),
    c(file.path(root, "benchmarks/predictive-geometry/restore-counts.R"), root, output),
    error_on_status = FALSE)
  expect_identical(run()$status, 0L)
  files <- file.path(output, paste0("predictive-geometry-", c("sampling", "selection"), "-pilot.rds"))
  for (i in seq_along(files)) {
    kind <- c("sampling", "selection")[i]; restored <- readRDS(files[i])
    expect_identical(restored$mode, "restored_pilot_count_contract")
    expect_identical(restored$counts, x[[kind]]$counts)
    expect_identical(restored$config_hash, x[[kind]]$config_hash)
    expect_identical(restored$harness, x[[kind]]$harness)
    expect_null(restored$results)
  }
  hashes <- vapply(files, function(path) digest::digest(file = path), "")
  expect_identical(run()$status, 0L)
  expect_identical(vapply(files, function(path) digest::digest(file = path), ""), hashes)
  corrupted <- readRDS(files[1]); corrupted$counts[1] <- corrupted$counts[1] + 1000
  saveRDS(corrupted, files[1]); before <- digest::digest(file = files[1])
  expect_identical(run()$status, 1L)
  expect_identical(digest::digest(file = files[1]), before)
})
