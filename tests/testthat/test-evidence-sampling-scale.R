test_that("committed sampling-covariance scale evidence passes its gates", {
  path <- testthat::test_path(
    "..", "..", "benchmark-results", "sampling-covariance-scale.rds"
  )
  skip_if_not(file.exists(path), "scale evidence has not been generated")
  artifact <- readRDS(path)

  expect_identical(artifact$schema_version, 1L)
  expect_identical(artifact$numerical_contract,
    "exact_factorized_hadamard_gram")
  expect_gte(artifact$dimensions$domain_features, 50000L)
  expect_gte(artifact$dimensions$experimental_conditions, 120L)
  expect_gte(artifact$dimensions$evidence_dimension, 7000L)
  expect_gte(artifact$dimensions$mean_support, 100)
  expect_true(artifact$dense_materialization_refused_at_64_mib)
  expect_true(all(unlist(artifact$checks, use.names = FALSE)))
  expect_true(all(is.finite(artifact$operations$elapsed_seconds)))
  expect_true(all(artifact$operations$elapsed_seconds >= 0))
  expect_lt(max(artifact$operations$resident_delta_bytes), 1024^3)
  expect_lt(max(artifact$operations$result_bytes), 512 * 1024^2)
})
