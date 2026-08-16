fmridesign_adapter_fixture <- function() {
  skip_if_not_installed("fmridesign", "0.6.0")
  if (!identical(as.character(utils::packageVersion("fmridesign")), "0.6.0")) {
    skip("The installed fmridesign version is outside the certified court.")
  }
  bound <- bound_study_fixture()
  study_value <- study(
    bound$observations, bound$events, bound$confounds, bound$hierarchy
  )
  event_data <- bound$events$data
  block <- match(event_data$partition, bound$fixture$partitions)
  frame <- fmridesign::sampling_frame(
    blocklens = bound$fixture$counts,
    TR = 2,
    start_time = 0
  )
  external <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition),
    data = event_data,
    block = block,
    sampling_frame = frame,
    durations = event_data$duration
  )
  model <- fmridesign_design_model(
    external,
    study_value,
    basis_id = "canonical-hrf-amplitude",
    units = "percent-signal-change"
  )
  list(bound = bound, study = study_value, external = external, model = model)
}

test_that("fmridesign compiles a study-bound semantic model", {
  fixture <- fmridesign_adapter_fixture()
  expect_s3_class(fixture$model, "effect_design_model")
  expect_identical(fixture$model$partitions, fixture$study$partitions)
  expect_identical(
    unname(vapply(fixture$model$designs, nrow, integer(1))),
    unname(fixture$bound$fixture$counts)
  )
  expect_true(fixture$model$capabilities$coding_invariant)

  coordinates <- fixture$model$condition_space$coordinates
  weights <- matrix(c(1, -1, 0), nrow = 1L,
    dimnames = list("face_minus_place", coordinates))
  effects <- effect_map(weights, fixture$model$condition_space)
  plan <- plan_relation(
    fixture$study,
    fixture$model,
    effects,
    observation_model("ols", sampling_unit = "scan")
  )
  expect_true(all(as.matrix(compiler_conformance(plan)[-1L])))
  expect_s3_class(estimate_relation(plan), "effect_relation_fit")
})

test_that("fmridesign compilation refuses a model built from other events", {
  fixture <- fmridesign_adapter_fixture()
  altered <- fixture$study$events$data
  altered$condition[[1L]] <- "place"
  altered_study <- study(
    fixture$bound$observations,
    observation_events(altered),
    fixture$bound$confounds,
    fixture$bound$hierarchy
  )
  refusal <- catch_refusal(fmridesign_design_model(
    fixture$external,
    altered_study,
    basis_id = "canonical-hrf-amplitude",
    units = "percent-signal-change"
  ))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "study_bound_compilation")
  expect_match(conditionMessage(refusal), "disagrees")
})
