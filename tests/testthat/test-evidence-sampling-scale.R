test_that("committed sampling-covariance scale evidence passes its gates", {
  artifact <- certified_artifact(
    "sampling-covariance-scale.rds", "run-sampling-covariance-scale.R"
  )

  expect_identical(artifact$schema_version, 1L)
  expect_identical(artifact$provenance$runner, "run-sampling-covariance-scale.R")
  expect_identical(artifact$numerical_contract,
    "exact_factorized_hadamard_gram")
  expect_gte(artifact$dimensions$domain_features, 50000L)
  expect_gte(artifact$dimensions$experimental_conditions, 120L)
  expect_gte(artifact$dimensions$evidence_dimension, 7000L)
  expect_gte(artifact$dimensions$mean_support, 100)
  expect_true(artifact$dense_materialization_refused_at_64_mib)
  expect_gt(artifact$dimensions$dense_covariance_payload_bytes, 64 * 1024^2)
  expect_setequal(names(artifact$checks), c(
    "diagonal_finite_nonnegative", "selected_finite",
    "quadratic_finite_nonnegative", "transport_symmetric",
    "no_dense_covariance_field"
  ))
  expect_true(all(unlist(artifact$checks, use.names = FALSE)))
  expect_setequal(artifact$operations$operation, c(
    "volume_frame_compile", "factor_construction", "diagonal",
    "selected_100", "quadratic_form", "transport_8"
  ))
  expect_true(all(is.finite(artifact$operations$elapsed_seconds)))
  expect_true(all(artifact$operations$elapsed_seconds >= 0))
  expect_lt(max(artifact$operations$resident_delta_bytes), 1024^3)
  expect_lt(max(artifact$operations$result_bytes), 512 * 1024^2)

  # The factorized route never forms the dense D-by-D field it replaces.
  query_operations <- artifact$operations[
    artifact$operations$operation != "volume_frame_compile",
  ]
  expect_lt(
    max(query_operations$result_bytes),
    artifact$dimensions$dense_covariance_payload_bytes
  )
})
