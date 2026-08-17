test_that("semantic relation plans are coding- and route-stable", {
  cell <- first_moment_vertical_fixture("cell", "qr", "fixed_gls")
  treatment <- first_moment_vertical_fixture(
    "treatment", "svd", "fixed_gls"
  )

  expect_identical(cell$study$study_id, treatment$study$study_id)
  expect_identical(cell$model$design_model_id,
    treatment$model$design_model_id)
  expect_identical(cell$effects$effect_map_id,
    treatment$effects$effect_map_id)
  expect_identical(cell$plan$relation_plan_id,
    treatment$plan$relation_plan_id)
  expect_false(identical(cell$model$compilation_route_id,
    treatment$model$compilation_route_id))
  expect_false(identical(
    cell$plan$design_receipts$`run-1`$design_receipt_id,
    treatment$plan$design_receipts$`run-1`$design_receipt_id
  ))
  expect_true(all(as.matrix(compiler_conformance(cell$plan)[-1L])))
  expect_identical(vapply(cell$plan$retained_rows, length, integer(1)),
    stats::setNames(rep(cell$retained_count, 4L), cell$partitions))
})

test_that("planned extraction agrees with direct E and legacy lm_relation_fit", {
  fixture <- first_moment_vertical_fixture("treatment", "svd", "fixed_gls")
  fit <- estimate_relation(fixture$plan)
  direct <- first_moment_direct_blocks(fixture)
  retained_sources <- stats::setNames(lapply(fixture$partitions,
    function(partition) {
      fixture$sources[[partition]][fixture$plan$retained_rows[[partition]], ,
        drop = FALSE]
    }), fixture$partitions)
  legacy <- lm_relation_fit(
    retained_sources,
    design = lapply(fixture$plan$design_receipts, `[[`, "design"),
    effects = lapply(fixture$plan$design_receipts, `[[`, "lowered_target"),
    observation_whitener = fixture$plan$whiteners,
    effect_names = fixture$effects$effect_space,
    domain = fixture$study$observations$domain,
    sampling_unit = "scan",
    solver = "svd"
  )

  for (partition in fixture$partitions) {
    expect_equal(relation_block(fit, partition, 1:60), direct[[partition]],
      tolerance = 2e-12)
    expect_equal(relation_block(fit, partition, 1:60),
      relation_block(legacy, partition, 1:60), tolerance = 2e-12)
  }
})

test_that("geometry, RDM, and RSA agree with independent form oracles", {
  fixture <- first_moment_vertical_fixture()
  fit <- estimate_relation(fixture$plan)
  plan <- plan_geometry(fit$relation, fixture$frame, fixture$over)
  direct_blocks <- first_moment_direct_blocks(fixture)
  direct_forms <- first_moment_direct_forms(
    direct_blocks, fixture$frame, fixture$over
  )
  direct_rdm <- first_moment_direct_rdm(direct_forms)

  distances <- rdm(plan)
  expect_equal(unname(as.matrix(distances$values)), direct_rdm,
    tolerance = 3e-11)

  weights <- c(face = 0.5, body = 0.5, house = -0.5, tool = -0.5)
  effect <- contrast_energy(plan, weights)
  direct_contrast <- vapply(direct_forms, function(form) {
    drop(crossprod(weights, form %*% weights))
  }, numeric(1))
  expect_equal(unname(effect$total), direct_contrast, tolerance = 3e-11)
  expect_equal(effect$total, effect$coherent + effect$configuration,
    tolerance = 2e-12)

  labels <- fixture$conditions$coordinates
  category <- outer(c(0, 0, 1, 1), c(0, 0, 1, 1), function(x, y) x != y)
  storage.mode(category) <- "double"
  dimnames(category) <- list(labels, labels)
  model_vector <- category[upper.tri(category)]
  expected_rsa <- t(apply(direct_rdm, 1L, function(value) {
    stats::lm.fit(cbind(intercept = 1, category = model_vector), value)$coefficients
  }))
  got_rsa <- rsa(plan, models = list(category = category))
  expect_equal(unname(got_rsa$coefficients), unname(expected_rsa),
    tolerance = 3e-11)
})

test_that("analytic covariance is route-stable when its assumptions are earned", {
  fixture <- first_moment_vertical_fixture()
  fit <- estimate_relation(fixture$plan)
  retained_sources <- stats::setNames(lapply(fixture$partitions,
    function(partition) {
      fixture$sources[[partition]][fixture$plan$retained_rows[[partition]], ,
        drop = FALSE]
    }), fixture$partitions)
  legacy <- lm_relation_fit(
    retained_sources,
    design = lapply(fixture$plan$design_receipts, `[[`, "design"),
    effects = lapply(fixture$plan$design_receipts, `[[`, "lowered_target"),
    observation_whitener = fixture$plan$whiteners,
    effect_names = fixture$effects$effect_space,
    domain = fixture$study$observations$domain,
    sampling_unit = "scan"
  )
  native_plan <- plan_geometry(fit$relation, fixture$frame, fixture$over)
  legacy_plan <- plan_geometry(legacy$relation, fixture$frame, fixture$over)
  native <- rdm_sampling_covariance(native_plan, fit, target = "plugin", at = 3L)
  reference <- rdm_sampling_covariance(
    legacy_plan, legacy, target = "plugin", at = 3L
  )
  expect_equal(sampling_covariance(native), sampling_covariance(reference),
    tolerance = 3e-11)
})

test_that("learned observation models keep points and withhold analytic law", {
  fixed <- first_moment_vertical_fixture(
    "cell", "qr", "fixed_gls"
  )
  learned <- first_moment_vertical_fixture(
    "cell", "qr", "learned_frozen_gls"
  )
  fixed_fit <- estimate_relation(fixed$plan)
  learned_fit <- estimate_relation(learned$plan)
  for (partition in fixed$partitions) {
    expect_equal(
      relation_block(fixed_fit, partition, 1:60),
      relation_block(learned_fit, partition, 1:60),
      tolerance = 2e-12
    )
  }
  learned_plan <- plan_geometry(
    learned_fit$relation, learned$frame, learned$over
  )
  refusal <- catch_refusal(rdm_sampling_covariance(
    learned_plan, learned_fit, target = "plugin", at = 1L
  ))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_match(conditionMessage(refusal), "no error channel")
})

test_that("the external fmrireg route reproduces points but not uncertainty", {
  skip_if_not_installed("fmrireg", "0.1.2")
  if (!identical(as.character(utils::packageVersion("fmrireg")), "0.1.2")) {
    skip("The installed fmrireg version is outside the certified court.")
  }
  fixture <- first_moment_vertical_fixture("cell", "qr", "ols")
  native <- estimate_relation(fixture$plan)
  external <- fmrireg_relation(fixture$plan)
  for (partition in fixture$partitions) {
    expect_equal(
      relation_block(native, partition, 1:60),
      relation_block(external, partition, 1:60),
      tolerance = 2e-12
    )
  }
  external_plan <- plan_geometry(
    external$relation, fixture$frame, fixture$over
  )
  refusal <- catch_refusal(rdm_sampling_covariance(
    external_plan, external, target = "null", at = 1L
  ))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_match(conditionMessage(refusal), "no error channel")
})

test_that("vertical-slice ambiguities refuse before scientific execution", {
  fixture <- first_moment_vertical_fixture("cell", "qr", "fixed_gls")

  shifted <- fixture$model$row_ids
  shifted$`run-2`[[1L]] <- 999L
  shifted_designs <- fixture$model$designs
  rownames(shifted_designs$`run-2`) <- as.character(shifted$`run-2`)
  misaligned_model <- design_model(
    fixture$model$specification,
    fixture$conditions,
    shifted_designs,
    fixture$model$parameterizations,
    row_ids = shifted,
    solver = fixture$model$solver,
    protocol = fixture$model$compiler$protocol,
    protocol_version = fixture$model$compiler$protocol_version,
    package = fixture$model$compiler$package,
    package_version = fixture$model$compiler$package_version
  )
  alignment <- catch_refusal(plan_relation(
    fixture$study, misaligned_model, fixture$effects, fixture$observation
  ))
  expect_s3_class(alignment, "effect_capability_refusal")
  expect_identical(alignment$capability, "aligned_observations")
  expect_match(conditionMessage(alignment), "run-2")

  undeclared <- plan_geometry(
    estimate_relation(fixture$plan)$relation,
    fixture$frame,
    cross_partitions(fixture$partitions, generalizes_over = "run")
  )
  capability <- sampling_capabilities(
    undeclared, estimate_relation(fixture$plan)
  )
  expect_false(capability$available)
  expect_true("endpoint_independence_not_declared" %in%
    capability$reasons$reason)

  design <- fixture$model$designs
  design <- lapply(design, function(value) {
    cbind(value, duplicate = value[, 1L])
  })
  target <- cbind(
    fixture$plan$design_receipts$`run-1`$lowered_target,
    duplicate = c(1, 0, 0, 0)
  )
  target[1L, 1L] <- 0
  deficient <- raw_design_model(
    design,
    row_ids = fixture$model$row_ids,
    solver = "svd"
  )
  nonestimable <- catch_refusal(plan_relation(
    fixture$study,
    deficient,
    raw_effect_map(target, fixture$effects$effect_space, colnames(target)),
    fixture$observation
  ))
  expect_s3_class(nonestimable, "effect_capability_refusal")
  expect_identical(nonestimable$capability, "estimable_effects")
  expect_match(conditionMessage(nonestimable), "run-1")
})

test_that("the committed first-moment benchmark receipt passes", {
  artifact <- certified_artifact(
    "first-moment-vertical-slice.rds", "run-first-moment-vertical-slice.R"
  )
  expect_identical(artifact$schema_version, 1L)
  expect_identical(artifact$provenance$runner,
    "run-first-moment-vertical-slice.R")
  expect_identical(artifact$fixture_version,
    "first-moment-vertical-slice:v1")
  expect_true(artifact$gate$passed)
  expect_true(all(unlist(artifact$checks[c(
    "relation_plan_identity_stable", "design_receipts_route_specific",
    "conformance_pass", "covariance_finite", "timing_complete",
    "memory_complete"
  )], use.names = FALSE)))
  expect_lt(artifact$checks$native_oracle_error, 1e-10)
  expect_lt(artifact$checks$coding_error, 1e-10)
  expect_lt(artifact$checks$fmrireg_error, 1e-10)
  expect_lt(artifact$checks$rdm_oracle_error, 1e-10)
  expect_setequal(artifact$timings$path, c(
    "plan_relation_cell_fixed", "plan_relation_treatment_fixed",
    "estimate_native_fixed", "compiler_conformance",
    "fmrireg_relation_ols", "plan_geometry", "contrast_query_first",
    "rdm_query_first", "rsa_query_first", "rdm_sampling_covariance",
    "bids_events", "bids_confounds", "bids_study",
    "fmridesign_design_model"
  ))

  # The receipt is only evidence for the adapter versions it was recorded
  # against, so it names them rather than implying version independence.
  expect_named(artifact$adapter_versions, c("fmridesign", "fmrireg"))
  expect_true(all(nzchar(artifact$adapter_versions)))
})
