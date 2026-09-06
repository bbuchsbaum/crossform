if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

matched_baseline_environment <- function() {
  environment <- new.env(parent = globalenv())
  for (file in c(
    "00-mixture-generator.R", "01-multiscale-scenarios.R",
    "02-paired-observations.R", "03-conventional-baselines.R"
  )) {
    sys.source(testthat::test_path(
      "..", "..", "benchmarks", "matched-interpretability", file
    ), envir = environment)
  }
  environment
}

test_that("aggregate magnitude is matched across fixed-total organizations", {
  generator <- matched_baseline_environment()
  bundle <- generator$matched_multiscale_scenarios()
  truth <- generator$matched_conventional_truth(bundle)
  summaries <- truth$summaries

  expect_equal(summaries$aggregate_crossvalidated_magnitude,
               rep(bundle$metadata$total_magnitude, 3L), tolerance = 1e-12)
  expect_true(length(unique(round(
    summaries$regional_activation, 12L
  ))) == 3L)
})

test_that("a pure and mixed pair meet the ambiguity criterion", {
  generator <- matched_baseline_environment()
  truth <- generator$matched_conventional_truth()
  broad_mixed <- generator$matched_ambiguity(
    truth, "broad_coherent", "mixed_broad_fine"
  )
  fine_mixed <- generator$matched_ambiguity(
    truth, "fine_configuration", "mixed_broad_fine"
  )

  expect_true(broad_mixed$ambiguous)
  expect_true(broad_mixed$spectrum_separates)
  expect_true(broad_mixed$criterion_passes)
  expect_lte(broad_mixed$baseline_gap, 1e-12)
  expect_gte(broad_mixed$spectrum_max_gap, 0.45)
  expect_true(fine_mixed$ambiguous)
  expect_true(fine_mixed$spectrum_separates)
  expect_true(fine_mixed$criterion_passes)
})

test_that("regional activation is a sufficient negative control", {
  generator <- matched_baseline_environment()
  truth <- generator$matched_conventional_truth()
  activation <- generator$matched_ambiguity(
    truth, "broad_coherent", "mixed_broad_fine",
    baseline = "regional_activation", baseline_tolerance = 1e-12
  )

  expect_false(activation$ambiguous)
  expect_true(activation$spectrum_separates)
  expect_false(activation$criterion_passes)
  expect_gt(activation$baseline_gap, 0.1)
})

test_that("direct comparators equal public whole-territory readings", {
  generator <- matched_baseline_environment()
  simulation <- generator$matched_paired_observations()
  for (noise_regime in c("gaussian", "heteroskedastic",
                         "spatial_correlated")) {
    for (scenario in names(simulation$bundle$scenarios)) {
      fitted <- generator$matched_conventional_fit(
        simulation, noise_regime, 24L, 0.8, scenario
      )
      info <- paste(noise_regime, scenario)
      expect_equal(
        fitted$conventional$regional_activation,
        fitted$crossform_global$regional_activation,
        tolerance = 2e-12, info = info
      )
      expect_equal(
        fitted$conventional$aggregate_crossvalidated_magnitude,
        fitted$crossform_global$aggregate_crossvalidated_magnitude,
        tolerance = 2e-12, info = info
      )
      expect_equal(
        sum(fitted$spectrum$total),
        fitted$conventional$aggregate_crossvalidated_magnitude,
        tolerance = 2e-12, info = info
      )
    }
  }
})

test_that("comparator formulas use every unordered partition pair", {
  generator <- matched_baseline_environment()
  blocks <- list(
    run1 = rbind(a = c(2, 0), b = c(0, 0)),
    run2 = rbind(a = c(1, 1), b = c(0, 0)),
    run3 = rbind(a = c(0, 2), b = c(0, 0))
  )
  observed <- generator$conventional_summaries_from_blocks(
    blocks, c(a = 1, b = -1)
  )

  expect_identical(length(observed$pair_magnitudes), 3L)
  expect_equal(observed$pair_magnitudes, c(2, 0, 2), tolerance = 0)
  expect_equal(observed$aggregate_crossvalidated_magnitude, 4 / 3,
               tolerance = 1e-15)
  expect_equal(observed$regional_activation, 1, tolerance = 1e-15)
})

test_that("the ambiguity function refuses malformed comparisons", {
  generator <- matched_baseline_environment()
  truth <- generator$matched_conventional_truth()
  expect_error(generator$matched_ambiguity(
    truth, "broad_coherent", "mixed_broad_fine", baseline = "accuracy"
  ), "Unknown conventional baseline")
  expect_error(generator$matched_ambiguity(
    truth, "broad_coherent", "broad_coherent"
  ), "two distinct scenarios")
})
