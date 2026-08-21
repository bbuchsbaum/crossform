learned_crossnobis_oracle <- function(fixture, frame, over, contrast, recipe,
                                      training) {
  fit <- fixture$fit
  weights <- as.matrix(frame$weights)
  partitions <- fit$relation$partitions
  vapply(seq_len(nrow(weights)), function(node) {
    support <- which(weights[node, ] > 0)
    root <- sqrt(weights[node, support])
    effects <- lapply(partitions, function(partition) {
      value <- relation_block(fit, partition, support)
      drop(contrast %*% value) * root
    })
    names(effects) <- partitions
    sum(vapply(seq_len(nrow(over)), function(edge) {
      metric_partitions <- if (identical(training$kind,
          "exclude_evaluation")) {
        setdiff(partitions, c(over$left[[edge]], over$right[[edge]]))
      } else {
        partitions
      }
      raw <- oracle_local_residual_covariance(
        fixture, metric_partitions, support
      )
      covariance <- oracle_shrinkage_covariance(raw, recipe)
      over$weight[[edge]] * drop(
        effects[[over$left[[edge]]]] %*%
          solve(covariance, effects[[over$right[[edge]]]])
      )
    }, numeric(1)))
  }, numeric(1))
}

test_that("learned crossnobis plans execute one exact support-streamed kernel", {
  setup <- metric_learning_setup()
  recipe <- shrinkage_precision(
    shrinkage = 0.2,
    relative_variance_floor = 1e-7,
    relative_spectral_floor = 1e-9
  )
  training <- metric_training_policy("exclude_evaluation")
  plan <- plan_crossnobis(
    setup$fixture$fit, setup$fixture$frame, setup$over,
    metric = recipe, training = training,
    compute = compute_policy(workspace_bytes = 64 * 1024^2),
    residual_workspace_bytes = setup$budgets$wider
  )
  setup$fixture$reset_reads()
  observed <- crossnobis(plan, c(condition = 1, drift = 0))
  expected <- learned_crossnobis_oracle(
    setup$fixture, setup$fixture$frame, setup$over,
    c(condition = 1, drift = 0), recipe, training
  )

  expect_s3_class(plan, "effect_geometry_plan")
  expect_s3_class(observed, "effect_crossnobis_view")
  expect_equal(observed$values, expected, tolerance = 5e-10)
  expect_identical(observed$receipt$completion_status, "complete")
  expect_identical(observed$receipt$kernel_version,
    "support-streamed-scheduled-metric-v1")
  expect_identical(observed$metadata$execution_plan$lowering,
    "support_streamed_scheduled_metric_query_contraction")
  # Kernel diagnostics now sit under the component they describe, where every
  # other geometry route puts them.
  expect_false(observed$metadata$diagnostics$total$pair_atoms_materialized)
  expect_false(observed$metadata$diagnostics$total$pair_frame_materialized)
  expect_false(
    observed$metadata$diagnostics$total$metric_factor_table_retained
  )
  expect_identical(
    observed$metadata$diagnostics$total$metric_handles_derived,
    nrow(setup$fixture$frame$weights) * nrow(setup$over)
  )
  expect_true(all(vapply(
    observed$metadata$metric_receipts,
    function(receipt) receipt$residual_reads_during_derivation == 0L,
    logical(1)
  )))
  # Execution reads only evaluation endpoints. Metric derivation consumes the
  # frozen atomic residual statistics and performs no residual-source reads.
  expect_identical(
    observed$metadata$source_session$read_count,
    c(run1 = nrow(setup$fixture$frame$weights),
      run2 = nrow(setup$fixture$frame$weights), run3 = 0L)
  )
})

test_that("residual provenance policies share execution but retain identity", {
  setup <- metric_learning_setup()
  recipe <- shrinkage_precision(0.15)
  disjoint_policy <- metric_training_policy("exclude_evaluation")
  all_policy <- metric_training_policy(
    "all_partitions_residual_orthogonality",
    justification = paste(
      "Evaluation-run GLM residuals are admitted under fitted-effect and",
      "residual orthogonality; metric uncertainty remains uncalibrated."
    )
  )
  disjoint <- plan_crossnobis(
    setup$fixture$fit, setup$fixture$frame, setup$over,
    metric = recipe, training = disjoint_policy,
    residual_workspace_bytes = setup$budgets$wider
  )
  all_runs <- plan_crossnobis(
    setup$fixture$fit, setup$fixture$frame, setup$over,
    metric = recipe, training = all_policy,
    residual_workspace_bytes = setup$budgets$wider
  )
  disjoint_value <- crossnobis(disjoint, c(condition = 1, drift = 0))
  all_value <- crossnobis(all_runs, c(condition = 1, drift = 0))

  # The same assertions, read off the geometry plan: both policies compile one
  # plan class and lower to one kernel, and differ only in identity.
  expect_identical(class(disjoint), class(all_runs))
  expect_identical(disjoint$lowering, all_runs$lowering)
  expect_identical(disjoint$lowering,
    "derive_then_support_streamed_pair_contraction")
  expect_identical(disjoint_value$receipt$kernel_version,
    all_value$receipt$kernel_version)
  expect_identical(disjoint_value$metadata$execution_plan$lowering,
    all_value$metadata$execution_plan$lowering)
  expect_false(identical(disjoint$signature, all_runs$signature))
  expect_false(identical(
    disjoint$metric_schedule$signature,
    all_runs$metric_schedule$signature
  ))
  expect_identical(
    disjoint$metric_schedule$schedule$records$edge_1$training_partitions,
    "run3"
  )
  expect_identical(
    all_runs$metric_schedule$schedule$records$edge_1$training_partitions,
    c("run1", "run2", "run3")
  )
  expect_equal(
    disjoint_value$values,
    learned_crossnobis_oracle(
      setup$fixture, setup$fixture$frame, setup$over,
      c(condition = 1, drift = 0), recipe, disjoint_policy
    ), tolerance = 5e-10
  )
  expect_equal(
    all_value$values,
    learned_crossnobis_oracle(
      setup$fixture, setup$fixture$frame, setup$over,
      c(condition = 1, drift = 0), recipe, all_policy
    ), tolerance = 5e-10
  )
  expect_true(
    all_value$metadata$metric_schedule$
      calibration_requires_metric_uncertainty
  )
})

test_that("invalid training and memory plans fail before residual reads", {
  fixture <- residual_statistics_fixture()
  over <- pairing("run1", "run2", independence = "independent")
  # Construct a valid two-partition fit rather than mutating a signed object.
  fit <- lm_relation_fit(
    fixture$response[c("run1", "run2")],
    cbind(
      intercept = 1,
      condition = rep(c(-0.5, 0.5), length.out = 37L),
      drift = seq(-1, 1, length.out = 37L)
    ),
    rbind(condition = c(0, 1, 0), drift = c(0, 0, 1)),
    domain = fixture$domain
  )
  reads <- integer(2L)
  names(reads) <- c("run1", "run2")
  original <- lapply(fit$relation$sources, `[[`, "read")
  for (partition in names(original)) {
    local({
      name <- partition
      reader <- original[[partition]]
      fit$relation$sources[[name]]$read <- function(features) {
        reads[[name]] <<- reads[[name]] + 1L
        reader(features)
      }
      fit$error_models[[name]]$residual_source$read <- function(features) {
        reads[[name]] <<- reads[[name]] + 1L
        reader(features)
      }
    })
  }
  expect_error(
    plan_crossnobis(fit, fixture$frame, over),
    "leaves no residual partition"
  , class = "effect_input_error")
  expect_identical(reads, c(run1 = 0L, run2 = 0L))

  setup <- metric_learning_setup()
  setup$fixture$reset_reads()
  expect_error(
    plan_crossnobis(
      setup$fixture$fit, setup$fixture$frame, setup$over,
      compute = compute_policy(workspace_bytes = 1),
      residual_workspace_bytes = setup$budgets$wider
    ),
    "exceeding the 1-byte workspace budget"
  , class = "effect_input_error")
  expect_identical(setup$fixture$reads(),
    c(run1 = 0L, run2 = 0L, run3 = 0L))
})

test_that("shrinkage target and calibration limits are explicit", {
  recipe <- shrinkage_precision(0.2)
  expect_identical(recipe$hyperparameters$target,
    "sample_residual_diagonal")
  expect_true(recipe$hyperparameters$target_estimated)
  expect_false("ld_t" %in% getNamespaceExports("crossform"))
  expect_false("confidence_interval" %in% getNamespaceExports("crossform"))
})

test_that("workspace changes execution identity but not the scientific plan", {
  fixture <- residual_statistics_fixture()
  budgets <- residual_statistics_budgets(fixture)
  over <- pairing("run1", "run2", independence = "independent")
  recipe <- shrinkage_precision(0.2)
  narrow <- plan_crossnobis(
    fixture$fit, fixture$frame, over, metric = recipe,
    compute = compute_policy(workspace_bytes = 64 * 1024^2),
    residual_workspace_bytes = budgets$minimum
  )
  wide <- plan_crossnobis(
    fixture$fit, fixture$frame, over, metric = recipe,
    compute = compute_policy(workspace_bytes = 128 * 1024^2),
    residual_workspace_bytes = budgets$wider
  )

  expect_identical(narrow$scientific_plan_id, wide$scientific_plan_id)
  expect_identical(narrow$metric_schedule$signature,
    wide$metric_schedule$signature)
  expect_false(identical(narrow$signature, wide$signature))
  expect_identical(
    crossnobis(narrow, c(condition = 1, drift = 0))$values,
    crossnobis(wide, c(condition = 1, drift = 0))$values
  )
})

test_that("compiled support execution does not revalidate global state per node", {
  setup <- metric_learning_setup()
  plan <- plan_crossnobis(
    setup$fixture$fit, setup$fixture$frame, setup$over,
    metric = shrinkage_precision(0.2),
    training = metric_training_policy("exclude_evaluation"),
    residual_workspace_bytes = setup$budgets$wider
  )
  frame_checks <- 0L
  support_checks <- 0L
  original_frame <- crossform:::.validate_frame_for_compile
  original_support <- crossform:::.validate_support_index
  testthat::local_mocked_bindings(
    .validate_frame_for_compile = function(...) {
      frame_checks <<- frame_checks + 1L
      original_frame(...)
    },
    .validate_support_index = function(...) {
      support_checks <<- support_checks + 1L
      original_support(...)
    },
    .package = "crossform"
  )

  value <- crossnobis(plan, c(condition = 1, drift = 0))

  expect_true(all(is.finite(value$values)))
  expect_lt(frame_checks, nrow(plan$frame$weights))
  expect_lt(support_checks, nrow(plan$frame$weights))
})
