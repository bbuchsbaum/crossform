test_that("the generated fMRI example is deterministic and RNG-neutral", {
  set.seed(81L)
  before <- .Random.seed
  first <- example_fmri_effects(seed = 12L, dimensions = c(4L, 4L, 3L))
  expect_identical(.Random.seed, before)
  second <- example_fmri_effects(seed = 12L, dimensions = c(4L, 4L, 3L))
  expect_identical(first$fit$signature, second$fit$signature)
  expect_identical(first$truth, second$truth)

  rm(.Random.seed, envir = .GlobalEnv)
  invisible(example_fmri_effects(seed = 13L, dimensions = c(4L, 4L, 3L)))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("the generated fixture supports the complete newcomer workflow", {
  example <- example_fmri_effects()
  expect_s3_class(example$fit, "effect_relation_fit")
  expect_s3_class(example$domain, "effect_domain")
  expect_s3_class(example$frame, "effect_frame")
  expect_identical(names(example$contrast), example$fit$relation$effects)
  expect_identical(
    rownames(example$model_rdm), example$fit$relation$effects
  )

  plan <- plan_geometry(
    example$fit$relation, example$frame,
    cross_partitions(example$fit$relation, independence = "independent")
  )
  effect <- contrast_energy(plan, example$contrast)
  distances <- rdm(plan)
  model <- rsa(plan, models = list(category = example$model_rdm))
  peak <- which.max(effect$total)

  expect_identical(dim(distances$values), c(280L, 6L))
  expect_identical(dim(model$coefficients), c(280L, 2L))
  expect_true(peak %in% example$truth$signal_measurements)
  expect_gt(
    stats::median(effect$total[example$truth$signal_measurements]),
    stats::median(effect$total[-example$truth$signal_measurements])
  )

  covariance <- rdm_sampling_covariance(
    plan, example$fit, target = "plugin", at = peak
  )
  standard_errors <- sqrt(sampling_covariance(covariance))
  expect_length(standard_errors, 6L)
  expect_true(all(is.finite(standard_errors) & standard_errors > 0))
})

test_that("the generated fixture rejects ambiguous dimensions and noise", {
  expect_error(example_fmri_effects(dimensions = c(3, 3)), "three integers",
    class = "effect_input_error")
  expect_error(example_fmri_effects(partitions = 1L), "at least 2",
    class = "effect_input_error")
  expect_error(
    example_fmri_effects(trials_per_condition = 1L), "at least 2"
  , class = "effect_input_error")
  expect_error(example_fmri_effects(noise_sd = 0), "positive",
    class = "effect_input_error")
  expect_error(example_fmri_effects(searchlight_radius = NA_real_), "positive",
    class = "effect_input_error")
})

test_that("the fixture plants one pattern block and one mean-shift block", {
  example <- example_fmri_effects()
  truth <- example$truth

  expect_length(intersect(truth$planted_features, truth$mean_features), 0L)
  expect_length(
    intersect(truth$pattern_measurements, truth$mean_measurements), 0L
  )
  expect_identical(truth$signal_measurements,
    sort(union(truth$pattern_measurements, truth$mean_measurements)))

  plan <- plan_geometry(
    example$fit$relation, example$frame,
    cross_partitions(example$fit$relation, independence = "independent",
      generalizes_over = "run")
  )
  effect <- contrast_energy(plan, example$contrast)

  # Each block demonstrates the half of the decomposition it was built from.
  expect_true(all(
    effect$configuration[truth$pattern_measurements] >
      effect$coherent[truth$pattern_measurements]
  ))
  expect_gt(max(effect$coherent[truth$mean_measurements]),
    5 * max(effect$coherent[truth$pattern_measurements]))
  expect_gt(max(effect$configuration[truth$pattern_measurements]),
    3 * max(effect$configuration[truth$mean_measurements]))

  # Nothing outside the two blocks reproduces.
  outside <- effect$total[-truth$signal_measurements]
  expect_lt(max(outside), 0.1)
  ranked <- order(effect$total, decreasing = TRUE)
  expect_true(all(
    ranked[seq_along(truth$signal_measurements)] %in%
      truth$signal_measurements
  ))
})

test_that("dropping a planted block leaves its truth empty", {
  pattern_only <- example_fmri_effects(plant = "pattern")
  expect_length(pattern_only$truth$mean_features, 0L)
  expect_length(pattern_only$truth$mean_measurements, 0L)
  expect_identical(pattern_only$truth$signal_measurements,
    sort(pattern_only$truth$pattern_measurements))

  mean_only <- example_fmri_effects(plant = "mean")
  expect_length(mean_only$truth$planted_features, 0L)
  expect_gt(length(mean_only$truth$mean_features), 0L)

  # The noise draw does not depend on what was planted, so the two runs
  # differ only by the planted signal.
  both <- example_fmri_effects()
  expect_identical(both$truth$planted_features,
    pattern_only$truth$planted_features)
  expect_identical(both$truth$mean_features, mean_only$truth$mean_features)
})
