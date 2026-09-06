if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

pcg_fixture <- function() {
  effects <- effect_space(c("face", "house"), basis_id = "cov-gate:v1")
  sizes <- c(s01 = 5L, s02 = 6L, s03 = 7L)
  subjects <- stats::setNames(lapply(names(sizes), function(id) {
    features <- sizes[[id]]
    domain <- abstract_domain(features,
      coordinates = cbind(x = seq_len(features) - 1),
      feature_ids = paste0("f", seq_len(features)), id = id)
    values <- matrix(seq_len(2L * features), 2L, features,
      dimnames = list(c("face", "house"), NULL)) / features
    relation <- relation(list(run1 = values, run2 = values * 1.1),
      effects = effects, domain = domain)
    plan_geometry(relation, compile_frame(voxelwise(), domain),
      cross_partitions(relation))
  }), names(sizes))
  transport <- lapply(sizes, function(features) {
    anatomical_transport(
      native_coords = cbind(seq_len(features) - 1),
      group_coords = cbind(c(0, 3)), semantics = "budget"
    )
  })
  list(subjects = subjects, transport = transport)
}

test_that("unweighted node-wise OLS does not require cross-node covariance", {
  fixture <- pcg_fixture()
  plan <- plan_population(fixture$subjects, fixture$transport)
  fit <- estimate_population(plan, rbind(`face-house` = c(1, -1)))

  expect_s3_class(fit, "effect_population_result")
  expect_identical(plan$fit, list(
    kind = "ols", weights = "subject_constant", commuting = TRUE,
    evaluation_order = "transport_then_fit"
  ))
  expect_true(any(is.finite(fit$coefficients)))
  expect_null(fit$uncertainty$within)
  expect_identical(fit$uncertainty$between$design,
    plan$model$matrix)
})

test_that("covariance-dependent precision has one actionable refusal", {
  fixture <- pcg_fixture()
  refusal <- catch_refusal(plan_population(
    fixture$subjects, fixture$transport,
    normalization = "precision_weighted"
  ))

  expect_identical(refusal$capability, "precision_weighted_normalization")
  expect_true(all(c(
    "per_subject_budget_variance_unavailable",
    "cross_node_sampling_covariance_route_absent"
  ) %in% refusal$reasons))
  expect_match(refusal$remedies,
    "design/cross-node-covariance-contract.md", all = FALSE, fixed = TRUE)
  expect_match(refusal$message, "per-measurement blocks must not be assembled")
})

test_that("the future dense law is executable but not package inference", {
  environment <- new.env(parent = globalenv())
  sys.source(test_path("..", "..", "design", "oracles",
    "cross-node-covariance.R"), envir = environment)
  oracle <- environment$cross_node_covariance_oracle()

  expect_gte(oracle$psd_minimum, 0)
  expect_identical(oracle$sparse_roundtrip, oracle$covariance)
  expect_equal(oracle$transported_budget_variance,
    oracle$full_budget_variance, tolerance = 1e-12)
  expect_gt(abs(oracle$diagonal_shortcut - oracle$full_budget_variance), 0.5)
  expect_false(exists("cross_node_sampling_covariance",
    envir = asNamespace("crossform"), inherits = FALSE))
  expect_false("cross_node_sampling_covariance" %in% getNamespaceExports(
    "crossform"))
})

test_that("future validation rejects indefiniteness and dense over-allocation", {
  environment <- new.env(parent = globalenv())
  sys.source(test_path("..", "..", "design", "oracles",
    "cross-node-covariance.R"), envir = environment)

  expect_error(environment$future_cross_node_dense(
    matrix(c(1, 2, 2, 1), 2L, 2L)), "not positive semidefinite")
  expect_error(environment$future_cross_node_dense(
    diag(20), max_dense_bytes = 100), "exceeds its declared compute budget")
  expect_error(environment$future_cross_node_dense(
    matrix(c(1, 0.2, 0.1, 1), 2L, 2L)), "not symmetric")
})

test_that("transported covariance refusal names the future contract", {
  refusal <- .population_transport_covariance_refusal()
  expect_identical(refusal$capability, "transported_sampling_covariance")
  expect_identical(refusal$future,
    .population_cross_node_future_contract())
  expect_identical(refusal$future$status, "not_implemented")
  expect_true(all(c("dense_symmetric", "sparse_symmetric") %in%
    refusal$future$representations))
  expect_true(all(c("positive_semidefinite",
    "declared_error_model") %in% refusal$future$required_checks))
  expect_true("no_implicit_dense_materialization" %in%
    refusal$future$scaling_guardrails)
})
