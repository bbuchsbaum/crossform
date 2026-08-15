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
  effect <- contrast(plan, example$contrast)
  distances <- rdm(plan)
  model <- rsa(plan, models = list(category = example$model_rdm))
  peak <- which.max(effect$total)

  expect_identical(dim(distances$values), c(245L, 6L))
  expect_identical(dim(model$coefficients), c(245L, 2L))
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
  expect_error(example_fmri_effects(dimensions = c(3, 3)), "three integers")
  expect_error(example_fmri_effects(partitions = 1L), "at least 2")
  expect_error(
    example_fmri_effects(trials_per_condition = 1L), "at least 2"
  )
  expect_error(example_fmri_effects(noise_sd = 0), "positive")
  expect_error(example_fmri_effects(searchlight_radius = NA_real_), "positive")
})
