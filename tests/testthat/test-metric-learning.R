test_that("public metric recipes are compact, explicit specifications", {
  identity <- identity_metric()
  diagonal <- diagonal_precision()
  shrinkage <- shrinkage_precision(0.25)

  expect_s3_class(identity, "effect_metric_recipe")
  expect_null(identity$domain)
  expect_identical(identity$kind, "identity")
  expect_true(metric_capabilities(identity)$feature_additive)
  expect_identical(diagonal$kind, "diagonal_variance_precision")
  expect_true(metric_capabilities(diagonal)$feature_additive)
  expect_identical(shrinkage$kind,
    "fixed_diagonal_shrinkage_precision")
  expect_false(metric_capabilities(shrinkage)$feature_additive)
  expect_identical(shrinkage$hyperparameters$randomness, "none")
  expect_null(shrinkage$hyperparameters$seed)
  expect_error(shrinkage_precision(0), "finite number", fixed = TRUE)
  expect_error(diagonal_precision(relative_variance_floor = 0),
    "positive relative floor")
})

test_that("evaluation-residual reuse requires an explicit policy contract", {
  disjoint <- metric_training_policy("exclude_evaluation")
  expect_false(disjoint$includes_evaluation_residuals)
  refusal <- catch_refusal(
    metric_training_policy("all_partitions_residual_orthogonality")
  )
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "evaluation_residual_reuse")
  expect_identical(refusal$namespace, "metric_learning")
  expect_identical(refusal$reasons, "residual_reuse_justification_absent")
  expect_match(refusal$remedies, "justification", all = FALSE)
  all_runs <- metric_training_policy(
    "all_partitions_residual_orthogonality",
    justification = paste(
      "GLM residuals are treated as orthogonal to fitted effects;",
      "metric-estimation uncertainty remains uncalibrated."
    )
  )
  expect_true(all_runs$includes_evaluation_residuals)
  expect_false(identical(disjoint$signature, all_runs$signature))
})

test_that("metric compilation freezes provenance but retains no factor table", {
  setup <- metric_learning_setup()
  recipe <- shrinkage_precision(0.2)
  schedule <- crossform:::compile_metric_schedule(
    recipe, setup$statistics, setup$fixture$frame, setup$over
  )
  record <- schedule$records$edge_1

  expect_s3_class(schedule, "effect_frozen_metric_schedule")
  expect_identical(schedule$recipe_specification, recipe$signature)
  expect_true(crossform:::.same_domain_reference(
    schedule$recipe$domain, setup$statistics$domain
  ))
  expect_identical(record$evaluation_left, "run1")
  expect_identical(record$evaluation_right, "run2")
  expect_identical(record$training_partitions, "run3")
  expect_identical(names(record$source_revisions), "run3")
  expect_identical(schedule$execution$local_metric_storage,
    "none_derived_on_demand")
  expect_false(schedule$execution$retained_factor_table)
  expect_false(any(c("metrics", "factors", "operators") %in% names(schedule)))
  expect_true(schedule$capabilities$calibration_requires_metric_uncertainty)
  expect_identical(setup$fixture$reads(),
    c(run1 = 0L, run2 = 0L, run3 = 0L))
  expect_silent(crossform:::.validate_frozen_metric_schedule(
    schedule, deep = TRUE
  ))
})

test_that("on-demand shrinkage precision agrees with an independent oracle", {
  setup <- metric_learning_setup()
  schedule <- crossform:::compile_metric_schedule(
    shrinkage_precision(
      shrinkage = 0.2,
      relative_variance_floor = 1e-7,
      relative_spectral_floor = 1e-9
    ),
    setup$statistics, setup$fixture$frame, setup$over
  )
  provider <- crossform:::.metric_schedule_provider(schedule, 1L)
  handle <- provider$at(10L)
  support_positions <- handle$support_positions
  raw <- oracle_local_residual_covariance(
    setup$fixture, "run3", support_positions
  )
  expected_covariance <- oracle_shrinkage_covariance(raw, schedule$recipe)
  set.seed(8302)
  left <- matrix(rnorm(3 * length(support_positions)), 3)
  right <- matrix(rnorm(2 * length(support_positions)), 2)
  expected_form <- left %*% solve(expected_covariance) %*% t(right)
  metric <- handle$materialize()
  receipt <- provider$receipt()

  expect_equal(handle$covariance, expected_covariance, tolerance = 2e-13)
  expect_equal(handle$form(left, right), expected_form, tolerance = 4e-11)
  expect_equal(metric$value, solve(expected_covariance), tolerance = 4e-11)
  expect_equal(metric$inverse_representation$value,
    expected_covariance, tolerance = 2e-13)
  expect_identical(metric$estimation, "learned_frozen")
  expect_identical(metric$provenance$training_signature,
    schedule$records$edge_1$training_signature)
  expect_identical(receipt$residual_reads_during_derivation, 0L)
  expect_identical(
    receipt$atomic_accumulation_reads,
    c(run3 = schedule$statistics$execution$atomic$run3$residual_reads)
  )
  expect_identical(receipt$nodes_derived, 1L)
  expect_false(receipt$retained_factor_table)
  # Nothing in this block reads a residual source: the provider answers from
  # the frozen atomic statistics, and the oracle above regresses the raw
  # responses directly instead of calling `residual_block()`.
  expect_identical(setup$fixture$reads(),
    c(run1 = 0L, run2 = 0L, run3 = 0L))
})

test_that("providers reduce only local pairs in the canonical global order", {
  setup <- metric_learning_setup()
  schedule <- crossform:::compile_metric_schedule(
    shrinkage_precision(0.2),
    setup$statistics, setup$fixture$frame, setup$over,
    metric_training_policy(
      "all_partitions_residual_orthogonality",
      justification = "Exercise the multi-partition canonical reduction law."
    )
  )
  provider <- crossform:::.metric_schedule_provider(schedule, 1L)
  handle <- provider$at(10L)
  scope <- crossform:::.residual_pair_scope(
    schedule$statistics,
    schedule$records$edge_1$training_partitions
  )
  global_then_gather <- crossform:::.local_residual_covariance(
    schedule$support_index,
    handle$support_positions,
    scope$covariance
  )
  gather_then_reduce <- crossform:::.local_residual_scope_covariance(
    schedule$statistics,
    schedule$support_index,
    handle$support_positions,
    schedule$records$edge_1$training_partitions
  )

  expect_identical(gather_then_reduce, global_then_gather)
  expect_false(any(c("scope", "schedule") %in%
    ls(environment(provider$at), all.names = TRUE)))
  expect_identical(provider$receipt()$residual_reads_during_derivation, 0L)
})

test_that("identity and diagonal schedules use their exact fast actions", {
  setup <- metric_learning_setup()
  identity_schedule <- crossform:::compile_metric_schedule(
    identity_metric(), setup$statistics, setup$fixture$frame, setup$over
  )
  diagonal_schedule <- crossform:::compile_metric_schedule(
    diagonal_precision(relative_variance_floor = 1e-7),
    setup$statistics, setup$fixture$frame, setup$over
  )
  identity <- crossform:::.metric_schedule_provider(
    identity_schedule, 1L
  )$at(10L)
  diagonal <- crossform:::.metric_schedule_provider(
    diagonal_schedule, 1L
  )$at(10L)
  set.seed(8303)
  left <- matrix(rnorm(2 * length(identity$support)), 2)
  right <- matrix(rnorm(4 * length(identity$support)), 4)
  raw <- oracle_local_residual_covariance(
    setup$fixture, "run3", diagonal$support_positions
  )
  variance <- diag(raw)
  floor <- 1e-7 * mean(variance[variance > 0])
  expected_variance <- pmax(variance, floor)

  expect_identical(identity$form(left, right), left %*% t(right))
  expect_equal(diagonal$covariance, diag(expected_variance),
    tolerance = 2e-13)
  expect_equal(diagonal$form(left, right),
    left %*% diag(1 / expected_variance) %*% t(right),
    tolerance = 3e-11)
  expect_true(metric_capabilities(identity$materialize())$feature_additive)
  expect_true(metric_capabilities(diagonal$materialize())$feature_additive)
})

test_that("all-partitions and disjoint training are distinct estimators", {
  setup <- metric_learning_setup()
  recipe <- shrinkage_precision(0.15)
  disjoint <- crossform:::compile_metric_schedule(
    recipe, setup$statistics, setup$fixture$frame, setup$over,
    metric_training_policy("exclude_evaluation")
  )
  all_runs <- crossform:::compile_metric_schedule(
    recipe, setup$statistics, setup$fixture$frame, setup$over,
    metric_training_policy(
      "all_partitions_residual_orthogonality",
      justification = paste(
        "Reuse follows the fitted-effect/residual orthogonality contract;",
        "calibration must propagate metric uncertainty."
      )
    )
  )

  expect_identical(disjoint$records$edge_1$training_partitions, "run3")
  expect_identical(all_runs$records$edge_1$training_partitions,
    c("run1", "run2", "run3"))
  expect_false(identical(disjoint$signature, all_runs$signature))
  expect_false(identical(
    crossform:::materialize_metric(disjoint, 10L)$signature,
    crossform:::materialize_metric(all_runs, 10L)$signature
  ))
})

test_that("workspace-invariant statistics yield identical metric schedules", {
  narrow_setup <- metric_learning_setup("narrow")
  wide_statistics <- residual_pair_statistics(
    narrow_setup$fixture$fit, narrow_setup$fixture$frame,
    workspace_bytes = narrow_setup$budgets$wider
  )
  recipe <- shrinkage_precision(0.3)
  narrow <- crossform:::compile_metric_schedule(
    recipe, narrow_setup$statistics, narrow_setup$fixture$frame,
    narrow_setup$over
  )
  wide <- crossform:::compile_metric_schedule(
    recipe, wide_statistics, narrow_setup$fixture$frame,
    narrow_setup$over
  )
  narrow_metric <- crossform:::materialize_metric(narrow, 10L)
  wide_metric <- crossform:::materialize_metric(wide, 10L)

  expect_identical(narrow$signature, wide$signature)
  expect_identical(narrow_metric$value, wide_metric$value)
  expect_identical(narrow_metric$inverse_representation$value,
    wide_metric$inverse_representation$value)
})

test_that("metric schedules refuse leakage and identity mutations", {
  setup <- metric_learning_setup()
  only_two <- residual_pair_statistics(
    setup$fixture$fit, setup$fixture$frame,
    partitions = c("run1", "run2"),
    workspace_bytes = setup$budgets$wider
  )
  expect_error(
    crossform:::compile_metric_schedule(
      shrinkage_precision(), only_two, setup$fixture$frame, setup$over
    ),
    "leaves no residual partition"
  )
  schedule <- crossform:::compile_metric_schedule(
    shrinkage_precision(), setup$statistics, setup$fixture$frame, setup$over
  )
  mutated <- schedule
  mutated$records$edge_1$training_partitions <- "run1"
  expect_error(
    crossform:::.validate_frozen_metric_schedule(mutated, deep = TRUE),
    "identity is inconsistent"
  )
})
