test_that("recorded query-first scale evidence passes its gates", {
  path <- testthat::test_path(
    "..", "..", "benchmark-results", "query-first-scale-gate.rds"
  )
  skip_if_not(file.exists(path), "query-first scale evidence is unavailable")
  artifact <- readRDS(path)
  expect_identical(artifact$schema_version, 1L)
  expect_gte(artifact$dimensions$conditions, 100L)
  expect_gte(artifact$dimensions$rdm_coordinates, 4950L)
  expect_gte(artifact$dimensions$frame_nodes, 1000L)
  expect_true(artifact$gate$identity_pass)
  expect_true(artifact$gate$fused_not_slower_pass)
  expect_true(artifact$gate$selected_cheaper_pass)
  expect_identical(
    artifact$memory$scope,
    "query_first_paths_only_after_all_route_warmup"
  )
  expect_true(artifact$gate$memory_pass)
  expect_true(artifact$gate$passed)
})

test_that("the query-first scale benchmark passes its declared gate", {
  if (!identical(Sys.getenv("EFFECTAGRAM_RUN_SCALE_TESTS"), "true")) {
    skip("Set EFFECTAGRAM_RUN_SCALE_TESTS=true to run the query-first gate.")
  }
  skip_if_not_installed("processx")
  repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  output <- file.path(tempdir(), paste0(
    "effectagram-query-first-gate-", Sys.getpid()
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
  expect_gte(artifact$dimensions$conditions, 100L)
  expect_true(artifact$gate$numerical_pass)
  expect_true(artifact$gate$identity_pass)
  expect_true(artifact$gate$fused_not_slower_pass)
  expect_true(artifact$gate$selected_cheaper_pass)
  expect_true(artifact$gate$runtime_pass)
  expect_identical(
    artifact$memory$scope,
    "query_first_paths_only_after_all_route_warmup"
  )
  expect_true(artifact$gate$memory_pass)
  expect_true(artifact$gate$passed)
})
