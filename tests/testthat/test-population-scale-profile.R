if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

test_that("analytic scale bands retain estimand and coverage identities", {
  decomposition <- do.call(population_decomposition, pdec_results())
  classical <- population_scale_profile(decomposition, "age",
    interval = "classical")
  robust <- population_scale_profile(decomposition, "age", interval = "HC3")

  expect_s3_class(robust, "effect_population_scale_profile")
  expect_identical(classical$interval$method, "classical")
  expect_identical(robust$interval$method, "HC3")
  expect_identical(robust$interval$semantics, "pointwise")
  expect_identical(robust$interval$simultaneous_coverage,
                   "not_available_unimplemented_uncalibrated")
  expect_identical(robust$data$subject_set_id, classical$data$subject_set_id)
  expect_identical(robust$data$estimate, classical$data$estimate)
  expect_true(all(c("n", "fraction", "n_eff", "subject_set_id", "gap") %in%
                    names(robust$data)))
  expect_identical(robust$receipt$interpolation,
                   "none_gaps_and_subject_set_changes_retained")
  expect_match(robust$interval$calibration_scope, "no marginal claim")
})

test_that("wild-bootstrap bands are bound to all three component results", {
  fits <- pdec_results()
  decomposition <- do.call(population_decomposition, fits)
  boot <- Map(function(result, seed) population_wild_bootstrap(
    result, contrast = "age", replicates = 99L, seed = seed
  ), fits, 8101:8103)
  profile <- population_scale_profile(decomposition, "age",
    query = "face-house", interval = "wild_bootstrap", bootstrap = boot)

  expect_identical(profile$interval$method, "wild_bootstrap")
  expect_identical(profile$receipt$uncertainty,
                   "null_imposed_wild_bootstrap_HC3_t")
  expect_true(any(is.finite(profile$data$lower)))
  expect_true(all(profile$data$lower <= profile$data$estimate |
                    profile$data$gap))
  expect_true(all(profile$data$upper >= profile$data$estimate |
                    profile$data$gap))
})

test_that("calibration artifact supports only its declared pointwise regimes", {
  results <- utils::read.csv(testthat::test_path(
    "..", "..", "inst", "extdata", "certification",
    "population-calibration-results.csv"
  ), stringsAsFactors = FALSE)
  hc3 <- results$method == "HC3" & results$marginal_claim_supported
  lower <- results$coverage[hc3] - 2 * results$coverage_mcse[hc3]
  expect_gte(min(lower), 0.88)
  expect_true(all(!results$marginal_claim_supported[
    results$coverage_regime == "informative"
  ]))
})

test_that("bootstrap bands refuse mismatched component parents", {
  fits <- pdec_results()
  decomposition <- do.call(population_decomposition, fits)
  boot <- Map(function(result, seed) population_wild_bootstrap(
    result, contrast = "age", replicates = 99L, seed = seed
  ), fits, 8201:8203)
  boot$configuration <- boot$total
  expect_error(population_scale_profile(decomposition, "age",
    interval = "wild_bootstrap", bootstrap = boot),
    class = "effect_capability_refusal")
})
