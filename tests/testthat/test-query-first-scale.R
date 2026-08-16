test_that("recorded query-first scale evidence passes its gates", {
  artifact <- certified_artifact(
    "query-first-scale-gate.rds", "run-query-first-scale.R"
  )
  expect_identical(artifact$schema_version, 1L)
  expect_identical(artifact$provenance$runner, "run-query-first-scale.R")
  expect_gte(artifact$dimensions$conditions, 100L)
  expect_gte(artifact$dimensions$rdm_coordinates, 4950L)
  expect_gte(artifact$dimensions$frame_nodes, 1000L)

  # Re-derive each recorded verdict from the recorded measurements.
  expect_lte(
    artifact$summary$fused_to_materialized_ratio,
    artifact$gate$maximum_fused_to_materialized_ratio
  )
  expect_lte(
    artifact$summary$selected_median_seconds /
      artifact$summary$full_fused_median_seconds,
    artifact$gate$maximum_selected_to_full_ratio
  )
  expect_lte(
    artifact$summary$full_fused_median_seconds,
    artifact$gate$maximum_path_seconds
  )
  expect_lte(
    artifact$summary$independent_oracle_max_abs_error,
    artifact$gate$numerical_tolerance
  )
  expect_lte(
    artifact$memory$incremental_peak_r_heap_bytes,
    artifact$gate$maximum_incremental_r_heap_bytes
  )

  expect_true(artifact$gate$identity_pass)
  expect_true(artifact$gate$fused_not_slower_pass)
  expect_true(artifact$gate$selected_cheaper_pass)
  expect_identical(
    artifact$memory$scope,
    "maximum_fresh_worker_r_heap_increment_across_query_first_paths"
  )
  expect_identical(artifact$memory$measurement, "r_gc_max_used")
  expect_setequal(
    artifact$memory$path_measurements$path,
    c("selected_100", "full_fused", "rsa_fused")
  )
  expect_true(artifact$gate$memory_pass)
  expect_identical(
    artifact$gate$passed,
    all(unlist(artifact$gate[c("numerical_pass", "identity_pass",
      "fused_not_slower_pass", "selected_cheaper_pass", "runtime_pass",
      "memory_pass")], use.names = FALSE))
  )
  expect_true(artifact$gate$passed)
})

test_that("the query-first scale benchmark passes its declared gate", {
  if (!identical(Sys.getenv("CROSSFORM_RUN_SCALE_TESTS"), "true")) {
    skip("Set CROSSFORM_RUN_SCALE_TESTS=true to run the query-first gate.")
  }
  skip_if_not_installed("processx")
  repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  output <- file.path(tempdir(), paste0(
    "crossform-query-first-gate-", Sys.getpid()
  ))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  run <- processx::run(
    file.path(R.home("bin"), "Rscript"),
    c(file.path(repo, "benchmarks", "run-query-first-scale.R"),
      repo, output),
    timeout = 15 * 60 * 1000,
    error_on_status = FALSE,
    echo = FALSE
  )
  expect_identical(run$status, 0L, info = paste(run$stdout, run$stderr))
  artifact <- readRDS(file.path(output, "query-first-scale-gate.rds"))
  expect_identical(artifact$schema_version, 1L)
  expect_identical(artifact$provenance$package_version,
    as.character(utils::packageVersion("crossform")))
  expect_gte(artifact$dimensions$conditions, 100L)
  expect_true(artifact$gate$numerical_pass)
  expect_true(artifact$gate$identity_pass)
  expect_true(artifact$gate$fused_not_slower_pass)
  expect_true(artifact$gate$selected_cheaper_pass)
  expect_true(artifact$gate$runtime_pass)
  expect_identical(
    artifact$memory$scope,
    "maximum_fresh_worker_r_heap_increment_across_query_first_paths"
  )
  expect_identical(artifact$memory$measurement, "r_gc_max_used")
  expect_setequal(
    artifact$memory$path_measurements$path,
    c("selected_100", "full_fused", "rsa_fused")
  )
  expect_true(artifact$gate$memory_pass)
  expect_true(artifact$gate$passed)
})
