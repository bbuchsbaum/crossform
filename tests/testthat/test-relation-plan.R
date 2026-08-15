test_that("relation plan identity is stable across coding and solver routes", {
  cell <- relation_plan_fixture("qr", "cell")
  treatment <- relation_plan_fixture("svd", "treatment")

  expect_s3_class(cell$plan, "effect_relation_plan")
  expect_identical(cell$model$design_model_id,
    treatment$model$design_model_id)
  expect_false(identical(cell$model$compilation_route_id,
    treatment$model$compilation_route_id))
  expect_identical(cell$effects$effect_map_id,
    treatment$effects$effect_map_id)
  expect_identical(cell$plan$relation_plan_id,
    treatment$plan$relation_plan_id)
  expect_false(identical(
    cell$plan$design_receipts$`run-1`$design_receipt_id,
    treatment$plan$design_receipts$`run-1`$design_receipt_id
  ))
  expect_identical(cell$plan$design_receipts$`run-1`$solver,
    "pivoted_qr")
  expect_identical(treatment$plan$design_receipts$`run-1`$solver,
    "svd_requested")
})

test_that("estimate emits the existing relation-fit contract with receipt identity", {
  cell <- relation_plan_fixture("qr", "cell")
  treatment <- relation_plan_fixture("svd", "treatment")
  cell_fit <- estimate(cell$plan)
  treatment_fit <- estimate(treatment$plan)

  expect_s3_class(cell_fit, "effect_relation_fit")
  expect_s3_class(treatment_fit, "effect_relation_fit")
  expect_identical(cell_fit$provenance$relation_plan_id,
    cell$plan$relation_plan_id)
  expect_identical(cell_fit$provenance$design_receipt_ids,
    vapply(cell$plan$design_receipts, `[[`, character(1), "design_receipt_id"))
  expect_false(identical(cell_fit$signature, treatment_fit$signature))

  for (partition in cell$plan$partitions) {
    expect_equal(
      relation_block(cell_fit, partition, 1:5),
      relation_block(treatment_fit, partition, 1:5),
      tolerance = 1e-12
    )
  }
  capabilities <- relation_fit_capabilities(cell_fit)
  expect_true(all(capabilities$residual_blocks))
  expect_true(all(capabilities$effect_covariance))
  expect_true(all(capabilities$residual_df))
})

test_that("planning is metadata-only and neural reads remain behind estimate", {
  fixture <- relation_plan_fixture(use_functions = TRUE)
  expect_identical(fixture$bound$fixture$reads$count, 0L)

  fit <- estimate(fixture$plan)
  expect_identical(fixture$bound$fixture$reads$count, 0L)
  relation_block(fit, "run-1", c(1, 3))
  expect_identical(fixture$bound$fixture$reads$count, 1L)
})

test_that("planned estimates agree with a direct censored linear oracle", {
  fixture <- relation_plan_fixture("svd", "cell", "fixed_gls")
  fit <- estimate(fixture$plan)

  for (partition in fixture$plan$partitions) {
    receipt <- fixture$plan$design_receipts[[partition]]
    rows <- fixture$plan$retained_rows[[partition]]
    response <- fixture$bound$fixture$sources[[partition]][rows, , drop = FALSE]
    whitener <- fixture$plan$whiteners[[partition]]
    direct <- receipt$lowered_target %*%
      relation_plan_inverse(whitener %*% receipt$design) %*%
      whitener %*% response
    expect_equal(relation_block(fit, partition, 1:5), direct,
      tolerance = 1e-12)
  }
})

test_that("learned observation models retain point fits but withhold analytic law", {
  fixed <- relation_plan_fixture("qr", "cell", "fixed_gls")
  learned <- relation_plan_fixture("qr", "cell", "learned_frozen_gls")
  fixed_fit <- estimate(fixed$plan)
  learned_fit <- estimate(learned$plan)

  for (partition in fixed$plan$partitions) {
    expect_equal(
      relation_block(fixed_fit, partition, 1:5),
      relation_block(learned_fit, partition, 1:5),
      tolerance = 1e-12
    )
  }
  learned_capabilities <- relation_fit_capabilities(learned_fit)
  expect_false(any(learned_capabilities$residual_blocks))
  expect_false(any(learned_capabilities$effect_covariance))
  refusal <- catch_refusal(rdm_sampling_covariance(
    plan_geometry(
      learned_fit$relation,
      compile_frame(whole_brain(), learned$bound$fixture$domain),
      cross_partitions(learned_fit$relation)
    ),
    learned_fit,
    target = "null"
  ))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_match(conditionMessage(refusal), "no error channel")
})

test_that("raw X and T plans bind parameterization and lack semantic claims", {
  cell <- relation_plan_fixture("qr", "cell")
  treatment <- relation_plan_fixture("svd", "treatment")
  cell_target <- cell$plan$lowered_effects$`run-1`
  treatment_target <- treatment$plan$lowered_effects$`run-1`
  raw_cell <- plan_relation(
    cell$plan$study,
    raw_design_model(
      cell$model$designs,
      row_ids = cell$model$row_ids,
      solver = "qr"
    ),
    raw_effect_map(
      cell_target$target,
      effects = cell_target$effect_space,
      coefficients = colnames(cell_target$target)
    ),
    cell$observation
  )
  raw_treatment <- plan_relation(
    treatment$plan$study,
    raw_design_model(
      treatment$model$designs,
      row_ids = treatment$model$row_ids,
      solver = "svd"
    ),
    raw_effect_map(
      treatment_target$target,
      effects = treatment_target$effect_space,
      coefficients = colnames(treatment_target$target)
    ),
    treatment$observation
  )

  expect_false(raw_cell$capabilities$symbolic_model)
  expect_false(raw_cell$capabilities$coding_invariant)
  expect_false(identical(raw_cell$relation_plan_id,
    raw_treatment$relation_plan_id))
  expect_equal(
    relation_block(estimate(raw_cell), "run-1", 1:5),
    relation_block(estimate(raw_treatment), "run-1", 1:5),
    tolerance = 1e-12
  )
})

test_that("plan identity changes with the request but not numerical tolerance", {
  fixture <- relation_plan_fixture()
  tolerance_route <- plan_relation(
    fixture$plan$study,
    fixture$model,
    fixture$effects,
    fixture$observation,
    tolerance = 1e-10
  )
  expect_identical(fixture$plan$relation_plan_id,
    tolerance_route$relation_plan_id)
  expect_false(identical(
    fixture$plan$design_receipts$`run-1`$design_receipt_id,
    tolerance_route$design_receipts$`run-1`$design_receipt_id
  ))

  changed_weights <- fixture$effects$weights
  changed_weights[1, ] <- c(-1, 1, 0)
  changed_effect <- plan_relation(
    fixture$plan$study,
    fixture$model,
    effect_map(changed_weights, fixture$conditions),
    fixture$observation
  )
  changed_sampling <- plan_relation(
    fixture$plan$study,
    fixture$model,
    fixture$effects,
    observation_model("ols", sampling_unit = "trial")
  )
  expect_length(unique(c(
    fixture$plan$relation_plan_id,
    changed_effect$relation_plan_id,
    changed_sampling$relation_plan_id
  )), 3L)
})

test_that("row alignment and nonestimability refuse before neural reads", {
  fixture <- relation_plan_fixture()
  model <- fixture$model
  model$row_ids$`run-1` <- rev(model$row_ids$`run-1`)
  model$compilation_route_id <- "tampered"
  expect_error(effectagram:::.validate_design_model(model),
    "disagree|inconsistent")

  reversed_designs <- fixture$model$designs
  reversed_ids <- lapply(fixture$model$row_ids, rev)
  for (partition in names(reversed_designs)) {
    rownames(reversed_designs[[partition]]) <-
      as.character(reversed_ids[[partition]])
  }
  raw <- raw_design_model(reversed_designs, row_ids = reversed_ids)
  target <- fixture$plan$lowered_effects$`run-1`
  refusal <- catch_refusal(plan_relation(
    fixture$plan$study,
    raw,
    raw_effect_map(target$target, target$effect_space,
      colnames(target$target)),
    fixture$observation
  ))
  expect_identical(refusal$capability, "aligned_observations")
  expect_match(conditionMessage(refusal), "run-1")

  deficient_designs <- fixture$model$designs
  deficient_designs$`run-1` <- cbind(
    deficient_designs$`run-1`,
    duplicate_intercept = deficient_designs$`run-1`[, 1]
  )
  raw_deficient <- raw_design_model(
    deficient_designs,
    row_ids = fixture$model$row_ids,
    solver = "auto"
  )
  bad_target <- cbind(
    fixture$plan$lowered_effects$`run-1`$target,
    duplicate_intercept = c(1, 0)
  )
  refusal <- catch_refusal(plan_relation(
    fixture$plan$study,
    raw_deficient,
    raw_effect_map(bad_target,
      fixture$plan$lowered_effects$`run-1`$effect_space,
      colnames(bad_target)),
    fixture$observation
  ))
  expect_identical(refusal$capability, "estimable_effects")
  expect_match(conditionMessage(refusal), "run-1")
})

test_that("design receipts are portable and identity guarded", {
  fixture <- relation_plan_fixture()
  receipts <- relation_plan_receipts(fixture$plan)

  expect_identical(names(receipts), fixture$plan$partitions)
  expect_true(all(vapply(receipts, function(value) {
    inherits(value, "effect_design_receipt") &&
      value$capabilities$portable_design_receipt
  }, logical(1))))
  expect_silent(unserialize(serialize(receipts, NULL, version = 3)))

  tampered <- receipts$`run-1`
  tampered$design[1, 1] <- tampered$design[1, 1] + 1
  expect_error(effectagram:::.validate_design_receipt(tampered),
    "identity is inconsistent")
})
