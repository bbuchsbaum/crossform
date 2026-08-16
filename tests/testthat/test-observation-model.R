test_that("observation models require an explicit sampling unit", {
  expect_error(observation_model("ols"), "declared explicitly")
  ols <- observation_model("ols", sampling_unit = "scan")

  expect_s3_class(ols, "effect_observation_model")
  expect_true(ols$capabilities$fixed_observation_model)
  expect_true(ols$capabilities$analytic_effect_covariance)
  expect_false(ols$capabilities$declared_independence)
  expect_null(ols$independence)
})

test_that("fixed and learned whiteners expose different capabilities", {
  whitener <- diag(seq(0.8, 1.2, length.out = 5L))
  fixed <- observation_model(
    "fixed_gls", "scan", whitener = whitener,
    independence = "within-run-whitened-errors"
  )
  learned <- observation_model(
    "learned_frozen_gls", "scan", whitener = whitener,
    training_revision = paste0("sha256:", paste(rep("a", 64), collapse = "")),
    training_provenance = list(method = "AR1", training_rows = 1:5)
  )

  expect_true(fixed$capabilities$fixed_observation_model)
  expect_true(fixed$capabilities$analytic_effect_covariance)
  expect_true(fixed$capabilities$declared_independence)
  expect_true(learned$capabilities$learned_observation_model)
  expect_false(learned$capabilities$fixed_observation_model)
  expect_false(learned$capabilities$analytic_effect_covariance)
  expect_false(identical(fixed$observation_model_id,
    learned$observation_model_id))
})

test_that("learned provenance is mandatory and fixed declarations are honest", {
  whitener <- diag(3)
  expect_error(observation_model(
    "learned_frozen_gls", "scan", whitener = whitener
  ), "training_revision")
  expect_error(observation_model(
    "learned_frozen_gls", "scan", whitener = whitener,
    training_revision = paste0("sha256:", paste(rep("b", 64), collapse = ""))
  ), "training_provenance")
  expect_error(observation_model(
    "fixed_gls", "scan", whitener = whitener,
    training_revision = paste0("sha256:", paste(rep("c", 64), collapse = ""))
  ), "only valid")
  expect_error(observation_model("ols", "scan", whitener = whitener),
    "must be NULL")
})

test_that("observation-model identity binds sampling and whitener values", {
  first <- observation_model("fixed_gls", "scan", whitener = diag(3))
  changed_sampling <- observation_model(
    "fixed_gls", "trial", whitener = diag(3)
  )
  changed_whitener <- observation_model(
    "fixed_gls", "scan", whitener = diag(c(1, 1, 2))
  )
  expect_length(unique(c(
    first$observation_model_id,
    changed_sampling$observation_model_id,
    changed_whitener$observation_model_id
  )), 3L)

  tampered <- first
  tampered$whitener[[1L]][1, 1] <- 2
  expect_error(crossform:::.validate_observation_model(tampered),
    "inconsistent")
})
