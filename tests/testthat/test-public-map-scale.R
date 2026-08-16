test_that("recorded public map-scale evidence passes its gates", {
  path <- testthat::test_path(
    "..", "..", "benchmark-results", "public-map-scale-gate.rds"
  )
  skip_if_not(file.exists(path), "public map-scale evidence is unavailable")
  artifact <- readRDS(path)
  expect_identical(artifact$schema_version, 1L)
  expect_gte(artifact$dimensions$features, 576L)
  expect_gte(artifact$dimensions$conditions, 12L)
  expect_gte(artifact$sampling_dimensions$frame_nodes, 576L)
  expect_gte(artifact$sampling_dimensions$covariance_dimension, 66L)
  expect_lte(
    artifact$numerical$independent_oracle_max_abs_error,
    artifact$gate$numerical_tolerance
  )
  expect_true(artifact$gate$passed)
})

test_that("the public map-scale benchmark passes its declared gate", {
  if (!identical(Sys.getenv("CROSSFORM_RUN_SCALE_TESTS"), "true")) {
    skip("Set CROSSFORM_RUN_SCALE_TESTS=true to run the public map gate.")
  }
  skip_if_not_installed("processx")
  repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  output <- file.path(tempdir(), paste0(
    "crossform-public-map-gate-", Sys.getpid()
  ))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  run <- processx::run(
    file.path(R.home("bin"), "Rscript"),
    c(file.path(repo, "benchmarks", "run-public-map-scale-gate.R"),
      repo, output),
    timeout = 5 * 60 * 1000,
    error_on_status = FALSE,
    echo = FALSE
  )
  expect_identical(run$status, 0L, info = paste(run$stdout, run$stderr))
  artifact <- readRDS(file.path(output, "public-map-scale-gate.rds"))
  expect_identical(artifact$schema_version, 1L)
  expect_gte(artifact$dimensions$features, 576L)
  expect_gte(artifact$dimensions$frame_nodes, 576L)
  expect_gte(artifact$dimensions$conditions, 12L)
  expect_gte(artifact$dimensions$partitions, 8L)
  expect_gte(nrow(artifact$timings), 6L)
  expect_gte(artifact$sampling_dimensions$frame_nodes, 576L)
  expect_gte(artifact$sampling_dimensions$conditions, 12L)
  expect_gte(artifact$sampling_dimensions$covariance_dimension, 66L)
  expect_gte(nrow(artifact$sampling$probe_timings), 3L)
  expect_true(artifact$gate$numerical_pass)
  expect_true(artifact$gate$implicit_runtime_pass)
  expect_true(artifact$gate$explicit_runtime_pass)
  expect_true(artifact$gate$relative_runtime_pass)
  expect_true(artifact$gate$sampling_numerical_pass)
  expect_true(artifact$gate$sampling_runtime_pass)
  expect_true(artifact$gate$memory_pass)
  expect_true(artifact$gate$passed)
})
