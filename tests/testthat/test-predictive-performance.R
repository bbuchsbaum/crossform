pg_performance_harness <- function() {
  root <- testthat::test_path("..", "..", "benchmarks", "predictive-geometry")
  skip_if_not(file.exists(file.path(root, "performance.R")), "Source benchmark harness is not installed")
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "performance-config.R"), env)
  sys.source(file.path(root, "performance.R"), env)
  env
}

test_that("[T25 T43 T56 T57] all declared workloads admit before reading neural data", {
  h <- pg_performance_harness(); config <- h$pg_performance_config()
  for (i in seq_len(nrow(config$cases))) {
    f <- h$pg_perf_fixture(config$cases[i, ], config)
    prepared <- crossform:::.geometry_fit_prepare(f$train, f$basis)
    expect_identical(prepared$plan$logical_shape, rep(f$case$model_rank, 2L))
    expect_equal(prepared$plan$packed_width, f$case$model_rank * (f$case$model_rank + 1L) / 2L)
    routes <- if (f$case$id == "large") "block" else c("memory", "block", "dense")
    for (route in routes) {
      admission <- h$pg_perf_admission(f, route)
      expect_lte(admission$planned_workspace_bytes, config$workspace_bytes)
      expect_identical(sum(f$reads$calls), 0L)
      if (route != "dense") expect_identical(admission$fit$workers, 1L)
    }
  }
})

test_that("[T24 T25 T43 T56] matched tiny computations agree with an independent dense partition oracle", {
  h <- pg_performance_harness(); config <- h$pg_performance_config()
  config$row_block <- 3L; config$block_features <- 17L
  f <- h$pg_perf_fixture(config$cases[1, ], config)
  directory <- tempfile("performance-test-"); dir.create(directory)
  withr::defer(unlink(directory, recursive = TRUE))
  public <- h$pg_perf_public(f, "block", directory)
  reads <- as.list(f$reads)
  f$reads$calls[] <- 0L; f$reads$columns[] <- 0L; f$reads$maximum_width <- 0L
  dense <- h$pg_perf_full_space(f)
  expect_identical(as.list(f$reads), reads)
  oracle <- h$pg_perf_reference(f, seq_len(f$case$nodes))
  expect_lt(h$pg_perf_error(public$forms, oracle$forms), 1e-10)
  expect_lt(h$pg_perf_error(public$table$gain, oracle$gain), 1e-10)
  for (field in c("forms", "spectra", "modes", "table", "diagnostics"))
    expect_equal(public[[field]], dense[[field]], tolerance = 1e-9)
  legacy <- h$pg_perf_legacy_overlap(f)
  expect_true(all(legacy$errors < 1e-10))
})

test_that("[T18 T53 T56] measured rank-path reuse preserves every fixed candidate gain", {
  h <- pg_performance_harness(); config <- h$pg_performance_config()
  f <- h$pg_perf_fixture(config$cases[1, ], config)
  invisible(capture.output(reuse <- h$pg_perf_reuse_comparison(f)))
  expect_true(all(reuse$errors < 1e-10))
  expect_identical(reuse$cached$calls$materializations, 2L)
  expect_identical(reuse$uncached$calls$materializations, 12L)
  expect_lt(sum(reuse$cached$reads$columns), sum(reuse$uncached$reads$columns))
  expect_lte(reuse$cached$maximum_planned_bytes, config$workspace_bytes)
})
