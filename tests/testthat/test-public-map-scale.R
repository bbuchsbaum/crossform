test_that("recorded public map-scale evidence passes its gates", {
  artifact <- certified_artifact(
    "public-map-scale-gate.rds", "run-public-map-scale-gate.R"
  )
  expect_identical(artifact$schema_version, 1L)
  expect_identical(artifact$provenance$runner, "run-public-map-scale-gate.R")
  expect_gte(artifact$dimensions$features, 576L)
  expect_gte(artifact$dimensions$conditions, 12L)
  expect_gte(artifact$dimensions$partitions, 8L)
  expect_gte(artifact$sampling_dimensions$frame_nodes, 576L)
  expect_gte(artifact$sampling_dimensions$covariance_dimension, 66L)
  expect_gte(nrow(artifact$timings), 6L)
  expect_gte(nrow(artifact$sampling$probe_timings), 3L)

  # The recorded booleans are re-derived from the recorded measurements, so a
  # gate that was written as passing while its own evidence contradicts it
  # fails here rather than being believed.
  expect_lte(
    artifact$numerical$independent_oracle_max_abs_error,
    artifact$gate$numerical_tolerance
  )
  expect_lte(
    artifact$numerical$explicit_implicit_max_abs_error,
    artifact$gate$numerical_tolerance
  )
  expect_lte(
    artifact$memory$incremental_peak_rss_bytes,
    artifact$gate$maximum_incremental_rss_bytes
  )
  expect_lte(
    artifact$summary$explicit_to_implicit_ratio,
    artifact$gate$maximum_explicit_to_implicit_ratio
  )
  expect_lte(
    artifact$summary$implicit_median_seconds,
    artifact$gate$maximum_path_seconds
  )
  expect_lte(
    artifact$summary$explicit_median_seconds,
    artifact$gate$maximum_path_seconds
  )
  expect_lte(
    artifact$sampling$full_sweep_seconds,
    artifact$gate$maximum_sampling_full_sweep_seconds
  )
  expect_true(artifact$sampling$full_sweep_finite_nonnegative)

  expect_true(artifact$gate$numerical_pass)
  expect_true(artifact$gate$implicit_runtime_pass)
  expect_true(artifact$gate$explicit_runtime_pass)
  expect_true(artifact$gate$relative_runtime_pass)
  expect_true(artifact$gate$sampling_numerical_pass)
  expect_true(artifact$gate$sampling_runtime_pass)
  expect_true(artifact$gate$memory_pass)
  expect_identical(
    artifact$gate$passed,
    all(unlist(artifact$gate[c("numerical_pass", "implicit_runtime_pass",
      "explicit_runtime_pass", "relative_runtime_pass",
      "sampling_numerical_pass", "sampling_runtime_pass", "memory_pass")],
      use.names = FALSE))
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
  expect_identical(artifact$provenance$package_version,
    as.character(utils::packageVersion("crossform")))
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
