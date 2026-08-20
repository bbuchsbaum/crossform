# The learned local metric as an ordinary compiler lowering.
#
# `test-crossnobis-learned.R` checks the science: oracle agreement, training
# provenance, read counts. This file checks that the science is delivered by
# the geometry compiler rather than by a second plan class and a second
# driver -- that the plan is an `effect_geometry_plan`, that
# `compile()` lowers it, that the refusals are capabilities rather than class
# checks, and that the numbers are the ones the retired driver produced,
# pinned as recorded golden values now that the driver itself is deleted.

learned_geometry_plan <- function(setup, ...) {
  plan_geometry(
    setup$fixture$fit, setup$fixture$frame, setup$over,
    metric = shrinkage_precision(0.2),
    residual_workspace_bytes = setup$budgets$wider,
    ...
  )
}

test_that("a metric recipe compiles a learned local geometry plan", {
  setup <- metric_learning_setup()
  plan <- learned_geometry_plan(setup)

  expect_s3_class(plan, "effect_geometry_plan")
  expect_identical(plan$metric_schedule$kind, "learned_local_before_frame")
  expect_identical(plan$metric_schedule$materialization, "on_demand_local")
  expect_identical(plan$metric_schedule$scope, "support_local")
  expect_false(plan$metric_schedule$feature_additive)
  expect_true(plan$metric_schedule$support_dense)
  expect_null(plan$metric_schedule$metric)
  expect_null(plan$metric_schedule$metric_signature)
  expect_s3_class(plan$metric_schedule$schedule,
    "effect_frozen_metric_schedule")
  expect_identical(plan$lowering,
    "derive_then_support_streamed_pair_contraction")
  # The schedule's own declarations must be the plan's: same evaluation edges,
  # same support graph, same neural domain.
  expect_identical(plan$metric_schedule$schedule$pairing, plan$pairing)
  expect_identical(
    plan$metric_schedule$schedule$support_index$signature,
    plan$frame$support_index$signature
  )
  # The plan reads at compile time, and says so where a reader can see it.
  expect_true(plan$execution_hints$plan_time_residual_accumulation)
  expect_identical(plan$execution_hints$residual_workspace_bytes,
    as.double(setup$budgets$wider))
  expect_match(
    paste(utils::capture.output(print(plan)), collapse = "\n"),
    "learned fixed_diagonal_shrinkage_precision"
  )
})

test_that("plan_crossnobis is sugar for the same compiled estimand", {
  setup <- metric_learning_setup()
  recipe <- shrinkage_precision(0.2)
  sugar <- plan_crossnobis(
    setup$fixture$fit, setup$fixture$frame, setup$over,
    metric = recipe, residual_workspace_bytes = setup$budgets$wider
  )
  direct <- plan_geometry(
    setup$fixture$fit, setup$fixture$frame, setup$over,
    metric = recipe, residual_workspace_bytes = setup$budgets$wider
  )

  expect_s3_class(sugar, "effect_geometry_plan")
  expect_identical(sugar$scientific_plan_id, direct$scientific_plan_id)
  expect_identical(sugar$signature, direct$signature)
})

test_that("a learned recipe refuses without a residual channel or support", {
  setup <- metric_learning_setup()
  bare <- catch_refusal(plan_geometry(
    setup$fixture$fit$relation, setup$fixture$frame, setup$over,
    metric = shrinkage_precision(0.2)
  ))
  expect_s3_class(bare, "effect_capability_refusal")
  expect_identical(bare$capability, "learned_metric_input")
  expect_match(conditionMessage(bare), "beta matrices alone")

  # A learned schedule accumulates support-local residual statistics, so a
  # frame with no support index has nothing to accumulate over.
  expect_error(
    plan_geometry(
      setup$fixture$fit,
      compile_frame(voxelwise(), setup$fixture$domain), setup$over,
      metric = shrinkage_precision(0.2)
    ),
    "explicit support index",
    class = "effect_input_error"
  )

  # The workspace-hint argument is meaningless without a recipe.
  expect_error(
    plan_geometry(
      setup$fixture$fit$relation, setup$fixture$frame, setup$over,
      residual_workspace_bytes = 1024
    ),
    "admitted only with a learned",
    class = "effect_input_error"
  )
})

test_that("the compiler lowers a learned plan to the scheduled kernel", {
  setup <- metric_learning_setup()
  plan <- learned_geometry_plan(setup)
  weights <- c(condition = 1, drift = 0)
  execution <- crossform:::.compile_geometry_execution_plan(
    plan, query = bilinear_query(tcrossprod(weights)), component = "total",
    signed_query = weights
  )

  expect_s3_class(execution, "effect_geometry_execution_plan")
  expect_identical(execution$lowering,
    "support_streamed_scheduled_metric_query_contraction")
  expect_identical(execution$kernel_version,
    "support-streamed-scheduled-metric-v1")
  expect_identical(execution$component, "total")
  expect_identical(execution$output_width, 1L)
  expect_identical(execution$row_tile, 1L)
  expect_identical(execution$coordinate_tile, 1L)
  expect_identical(execution$feature_block, as.integer(max(
    plan$metric_schedule$schedule$support_index$cost$support_size
  )))
  expect_identical(execution$task_count, as.double(plan$measurements))
  expect_identical(execution$signed_query, weights)
  # The validator re-derives lowering, kernel and identity through the same
  # widened helpers, so compiler and validator cannot drift apart.
  expect_silent(crossform:::.validate_geometry_execution_plan(execution))
})

test_that("a learned plan charges its residual statistics, not a metric", {
  setup <- metric_learning_setup()
  learned <- learned_geometry_plan(setup)
  weights <- c(condition = 1, drift = 0)
  learned_execution <- crossform:::.compile_geometry_execution_plan(
    learned, query = bilinear_query(tcrossprod(weights)), component = "total",
    signed_query = weights
  )
  fixed <- plan_geometry(
    setup$fixture$fit$relation, setup$fixture$frame, setup$over,
    metric = noise_precision(
      diag(setup$fixture$domain$n_features), setup$fixture$domain
    )
  )
  fixed_execution <- crossform:::.compile_geometry_execution_plan(
    fixed, query = bilinear_query(tcrossprod(weights)), component = "total"
  )

  # The resident term a learned schedule carries is the retained residual
  # pair statistics, and there is no materialized metric object to measure at
  # all. Inheriting the fixed branch's `object.size(metric$value)` accounting
  # would silently charge nothing for the dominant payload.
  index <- learned$metric_schedule$schedule$support_index
  pair_count <- length(index$pair_pattern@i)
  partitions <- length(setup$fixture$fit$relation$partitions)
  expect_gte(
    learned_execution$memory$categories[["resident_source"]],
    pair_count * 8 * partitions
  )
  expect_gt(
    learned_execution$memory$categories[["resident_source"]],
    fixed_execution$memory$categories[["resident_source"]]
  )
})

test_that("the lowering refuses every component it cannot form", {
  setup <- metric_learning_setup()
  plan <- learned_geometry_plan(setup)
  weights <- c(condition = 1, drift = 0)
  query <- bilinear_query(tcrossprod(weights))

  refusals <- list(
    coherent = catch_refusal(evaluate_geometry(plan, query = query,
      component = "coherent")),
    configuration = catch_refusal(evaluate_geometry(plan, query = query,
      component = "configuration")),
    contrast = catch_refusal(contrast_energy(plan, weights)),
    materialization = catch_refusal(materialize_geometry(plan)),
    unsigned = catch_refusal(evaluate_geometry(plan, query = query,
      component = "total"))
  )
  for (name in names(refusals)) {
    expect_s3_class(refusals[[name]], "effect_capability_refusal")
    expect_identical(refusals[[name]]$capability,
      "scheduled_metric_component")
    expect_identical(refusals[[name]]$namespace, "geometry_views")
    expect_match(conditionMessage(refusals[[name]]), "crossnobis|learned")
  }
  expect_identical(refusals$unsigned$reasons,
    "scheduled_metric_requires_signed_contrast")
  expect_identical(refusals$materialization$reasons,
    "scheduled_metric_full_materialization_not_admitted")

  # An identity recipe is a legitimate schedule but is not a noise model, so
  # the Mahalanobis reading refuses it where the fixed route refuses the
  # implicit identity metric.
  identity_plan <- plan_geometry(
    setup$fixture$fit, setup$fixture$frame, setup$over,
    metric = identity_metric(), residual_workspace_bytes = setup$budgets$wider
  )
  identity_refusal <- catch_refusal(crossnobis(identity_plan, weights))
  expect_s3_class(identity_refusal, "effect_capability_refusal")
  expect_identical(identity_refusal$capability, "declared_noise_metric")
  expect_identical(identity_refusal$reasons,
    "identity_recipe_is_not_a_noise_model")
})

test_that("the sampling layer detects the learned metric by capability", {
  setup <- metric_learning_setup()
  plan <- learned_geometry_plan(setup)
  fixed <- plan_geometry(
    setup$fixture$fit$relation, setup$fixture$frame, setup$over,
    metric = noise_precision(
      diag(setup$fixture$domain$n_features), setup$fixture$domain
    )
  )

  expect_true(crossform:::.metric_schedule_requires_metric_uncertainty(
    plan$metric_schedule
  ))
  expect_false(crossform:::.metric_schedule_requires_metric_uncertainty(
    fixed$metric_schedule
  ))
  expect_true(
    plan$metric_schedule$schedule$capabilities$
      calibration_requires_metric_uncertainty
  )
  # The refusal reaches the sampling layer through the descriptor, with no
  # test of the plan's class anywhere on the path.
  descriptor <- crossform:::.sampling_evidence_descriptor(plan)
  expect_identical(descriptor$record$metric_status, "learned")
  refusal <- catch_refusal(rdm_sampling_covariance(
    plan, setup$fixture$fit, target = "null", at = 1L
  ))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "fixed_metric_sampling_law")
  expect_identical(refusal$reasons, "learned_metric_law_not_admitted")
})

# The numerical equality with the driver retired in B3, pinned.
#
# Before `.execute_learned_crossnobis()` and `.plan_learned_crossnobis()` were
# deleted, the retired driver was run once on this fixture under exactly the
# arguments below and its outputs recorded in
# `fixtures/learned-crossnobis-golden.rds` (`$values`, `$contrast`,
# `$estimand`, `$metric`, `$pairing`, `$index`, `$kernel_version`,
# `$lowering`, `$source_read_count`). The equality guarantee B2 established is
# not weakened by the deletion: it is now asserted against those recorded
# values rather than against a live second engine.
test_that("the new lowering reproduces the retired driver's recorded values", {
  setup <- metric_learning_setup()
  recipe <- shrinkage_precision(
    shrinkage = 0.2, relative_variance_floor = 1e-7,
    relative_spectral_floor = 1e-9
  )
  training <- metric_training_policy("exclude_evaluation")
  arguments <- list(
    setup$fixture$fit, setup$fixture$frame, setup$over,
    metric = recipe, training = training,
    compute = compute_policy(workspace_bytes = 64 * 1024^2),
    residual_workspace_bytes = setup$budgets$wider
  )
  weights <- c(condition = 1, drift = 0)
  golden <- readRDS(test_path("fixtures/learned-crossnobis-golden.rds"))

  compiled <- crossnobis(do.call(plan_crossnobis, arguments), weights)

  expect_equal(compiled$values, golden$values, tolerance = 1e-12)
  expect_identical(compiled$contrast, golden$contrast)
  expect_identical(compiled$estimand, golden$estimand)
  # The metric identity is the frozen schedule signature: a content digest, so
  # any drift in the recipe, statistics, support graph or training record
  # would move it even where the values agreed to tolerance.
  expect_identical(compiled$metric, golden$metric)
  expect_identical(compiled$pairing, golden$pairing)
  expect_identical(compiled$index, golden$index)
  expect_identical(compiled$receipt$kernel_version, golden$kernel_version)
  expect_identical(
    compiled$metadata$execution_plan$lowering, golden$lowering
  )
  expect_identical(
    compiled$metadata$source_session$read_count, golden$source_read_count
  )
  # Bit-for-bit within one session, where BLAS is fixed: the recorded pin is
  # held to a tolerance because it crosses machines, not because the route is
  # allowed to wander between two runs of itself.
  expect_identical(
    compiled$values, crossnobis(do.call(plan_crossnobis, arguments),
      weights)$values
  )
  # The retired driver's private `support_tasks` stage label becomes the
  # geometry executor's `feature_tasks`.
  expect_true(
    "feature_tasks" %in% names(compiled$receipt$observed$stage_seconds)
  )
  expect_false(
    "support_tasks" %in% names(compiled$receipt$observed$stage_seconds)
  )
  expect_identical(compiled$receipt$task_partition_id,
    sprintf("ascending-features-%d", max(
      compiled$receipt$observed$tiles$feature_block
    )))
})

test_that("training policy moves the estimand and tiles move only execution", {
  setup <- metric_learning_setup()
  recipe <- shrinkage_precision(0.2)
  disjoint <- plan_geometry(
    setup$fixture$fit, setup$fixture$frame, setup$over, metric = recipe,
    training = metric_training_policy("exclude_evaluation"),
    residual_workspace_bytes = setup$budgets$wider
  )
  all_runs <- plan_geometry(
    setup$fixture$fit, setup$fixture$frame, setup$over, metric = recipe,
    training = metric_training_policy(
      "all_partitions_residual_orthogonality",
      justification = paste(
        "Evaluation-run GLM residuals are admitted under fitted-effect and",
        "residual orthogonality; metric uncertainty remains uncalibrated."
      )
    ),
    residual_workspace_bytes = setup$budgets$wider
  )

  # The training policy is part of what the number means.
  expect_false(identical(
    disjoint$scientific_plan_id, all_runs$scientific_plan_id
  ))

  # The residual cache budget is not. It is a capacity knob, so it moves the
  # execution signature and leaves the estimand alone.
  narrow <- plan_geometry(
    setup$fixture$fit, setup$fixture$frame, setup$over, metric = recipe,
    residual_workspace_bytes = setup$budgets$minimum
  )
  wide <- plan_geometry(
    setup$fixture$fit, setup$fixture$frame, setup$over, metric = recipe,
    residual_workspace_bytes = setup$budgets$wider
  )
  expect_identical(narrow$scientific_plan_id, wide$scientific_plan_id)
  expect_false(identical(narrow$signature, wide$signature))

  weights <- c(condition = 1, drift = 0)
  query <- bilinear_query(tcrossprod(weights))
  narrow_execution <- crossform:::.compile_geometry_execution_plan(
    narrow, query = query, component = "total", signed_query = weights
  )
  wide_execution <- crossform:::.compile_geometry_execution_plan(
    wide, query = query, component = "total", signed_query = weights
  )
  expect_identical(
    narrow_execution$scientific_plan_id, wide_execution$scientific_plan_id
  )
  expect_false(identical(
    narrow_execution$signature, wide_execution$signature
  ))
  expect_identical(
    crossnobis(narrow, weights)$values, crossnobis(wide, weights)$values
  )
})

test_that("the signed hint is executor-facing, not estimand-bearing", {
  setup <- metric_learning_setup()
  plan <- learned_geometry_plan(setup)
  weights <- c(condition = 1, drift = 0)
  query <- bilinear_query(tcrossprod(weights))
  positive <- crossform:::.compile_geometry_execution_plan(
    plan, query = query, component = "total", signed_query = weights
  )
  negative <- crossform:::.compile_geometry_execution_plan(
    plan, query = query, component = "total", signed_query = -weights
  )

  # `c` and `-c` name one estimand and are executed by two hints.
  expect_identical(
    positive$scientific_plan_id, negative$scientific_plan_id
  )
  expect_false(identical(positive$signature, negative$signature))

  # A hint that does not reproduce the compiled query is refused, so it can
  # never alter what is estimated.
  expect_error(
    crossform:::.compile_geometry_execution_plan(
      plan, query = query, component = "total",
      signed_query = c(condition = 1, drift = 1)
    ),
    "signed outer product",
    class = "effect_input_error"
  )
})
