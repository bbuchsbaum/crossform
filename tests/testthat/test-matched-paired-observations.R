if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

matched_observation_environment <- function() {
  environment <- new.env(parent = globalenv())
  for (file in c(
    "00-mixture-generator.R", "01-multiscale-scenarios.R",
    "02-paired-observations.R"
  )) {
    sys.source(testthat::test_path(
      "..", "..", "benchmarks", "matched-interpretability", file
    ), envir = environment)
  }
  environment
}

test_that("the manifest spans paired null low-power and recoverable cells", {
  generator <- matched_observation_environment()
  simulation <- generator$matched_paired_observations()
  manifest <- simulation$manifest

  expect_identical(simulation$schema_version,
                   "matched-paired-observations-v1")
  expect_identical(nrow(manifest), 54L)
  expect_setequal(manifest$scenario,
                  c("broad_coherent", "mixed_broad_fine",
                    "fine_configuration"))
  expect_setequal(manifest$noise_regime,
                  c("gaussian", "heteroskedastic", "spatial_correlated"))
  expect_setequal(manifest$trials_per_condition, c(6L, 24L))
  expect_setequal(manifest$snr, c(0, 0.2, 0.8))
  expect_setequal(manifest$regime_label,
                  c("null", "low_power", "recoverable"))
  expect_identical(anyDuplicated(manifest$cell_id), 0L)
  expect_identical(length(unique(manifest$base_noise_id)), 1L)
  expect_true(all(manifest$seed == simulation$metadata$seed))
})

test_that("matched organizations reuse exactly the same noise realization", {
  generator <- matched_observation_environment()
  simulation <- generator$matched_paired_observations()
  scenarios <- unique(simulation$manifest$scenario)
  for (noise_regime in unique(simulation$manifest$noise_regime)) {
    for (n in unique(simulation$manifest$trials_per_condition)) {
      for (signal_to_noise in unique(simulation$manifest$snr)) {
        residuals <- lapply(scenarios, function(scenario) {
          key <- generator$.paired_cell_key(
            noise_regime, n, signal_to_noise, scenario
          )
          cell <- simulation$cells[[key]]
          Map(function(response) response - cell$design %*% cell$effect_matrix,
              cell$responses)
        })
        for (scenario in seq_along(residuals)[-1L]) {
          expect_equal(residuals[[scenario]], residuals[[1L]],
                       tolerance = 5e-15,
                       info = paste(noise_regime, n, signal_to_noise))
        }
      }
    }
  }
})

test_that("sample-size cells are nested condition-stratified prefixes", {
  generator <- matched_observation_environment()
  simulation <- generator$matched_paired_observations()
  small_key <- generator$.paired_cell_key(
    "spatial_correlated", 6L, 0.8, "mixed_broad_fine"
  )
  large_key <- generator$.paired_cell_key(
    "spatial_correlated", 24L, 0.8, "mixed_broad_fine"
  )
  small <- simulation$cells[[small_key]]
  large <- simulation$cells[[large_key]]
  rows <- c(1:6, 24L + 1:6)

  expect_identical(small$design, large$design[rows, , drop = FALSE])
  for (partition in names(small$responses)) {
    expect_identical(small$responses[[partition]],
                     large$responses[[partition]][rows, , drop = FALSE])
  }
})

test_that("noise transforms are distinct and share unit average variance", {
  generator <- matched_observation_environment()
  simulation <- generator$matched_paired_observations()
  covariance <- simulation$source$noise_covariance

  expect_equal(covariance$gaussian, diag(17L), tolerance = 1e-15)
  expect_equal(mean(diag(covariance$heteroskedastic)), 1,
               tolerance = 1e-15)
  expect_gt(stats::sd(diag(covariance$heteroskedastic)), 0.2)
  expect_equal(mean(diag(covariance$spatial_correlated)), 1,
               tolerance = 1e-15)
  expect_equal(covariance$spatial_correlated[1L, 2L], 0.6,
               tolerance = 1e-15)
  expect_equal(covariance$spatial_correlated[1L, 3L], 0.6^2,
               tolerance = 1e-15)
})

test_that("truth is separate from observations and fixed across noise regimes", {
  generator <- matched_observation_environment()
  simulation <- generator$matched_paired_observations()
  expect_false(any(c("estimate", "standard_error") %in%
                     names(simulation$truth)))
  expect_false("truth" %in% names(simulation$cells[[1L]]))
  expect_identical(nrow(simulation$truth), 9L)
  expect_equal(simulation$truth$total[simulation$truth$snr == 0],
               rep(0, 3L), tolerance = 0)
  recoverable <- simulation$truth[simulation$truth$snr == 0.8, ]
  expect_equal(recoverable$total, rep(17 * 0.8^2, 3L), tolerance = 1e-12)
})

test_that("a compact seed rebuilds observations and truth exactly", {
  generator <- matched_observation_environment()
  set.seed(771L)
  before <- .Random.seed
  first <- generator$matched_paired_observations(seed = 4201L)
  expect_identical(.Random.seed, before)
  second <- generator$matched_paired_observations(seed = 4201L)

  expect_identical(first$manifest, second$manifest)
  expect_identical(first$truth, second$truth)
  expect_identical(first$scale_truth, second$scale_truth)
  expect_identical(first$source$base_standard_noise,
                   second$source$base_standard_noise)
  expect_identical(first$cells, second$cells)
})

test_that("public end-to-end fits execute all three organizations", {
  generator <- matched_observation_environment()
  simulation <- generator$matched_paired_observations()
  fits <- lapply(names(simulation$bundle$scenarios), function(scenario) {
    generator$matched_observation_fit(
      simulation, "gaussian", 24L, 0.8, scenario
    )
  })

  for (fit in fits) {
    expect_s3_class(fit$fit, "effect_relation_fit")
    expect_identical(nrow(fit$spectrum), 4L)
    expect_true(all(is.finite(fit$spectrum$total)))
    expect_true(all(is.finite(fit$spectrum$coherent)))
    expect_true(all(is.finite(fit$spectrum$configuration)))
    expect_identical(nrow(fit$truth), 1L)
    expect_identical(nrow(fit$scale_truth), 4L)
  }

  null_broad <- generator$matched_observation_fit(
    simulation, "gaussian", 24L, 0, "broad_coherent"
  )
  null_fine <- generator$matched_observation_fit(
    simulation, "gaussian", 24L, 0, "fine_configuration"
  )
  expect_identical(null_broad$spectrum$total, null_fine$spectrum$total)
  expect_identical(null_broad$spectrum$coherent, null_fine$spectrum$coherent)
  expect_identical(null_broad$spectrum$configuration,
                   null_fine$spectrum$configuration)
})
