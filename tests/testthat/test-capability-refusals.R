refusal_fixture <- function(partitions = 3L,
                            domain_id = "capability-refusal-domain") {
  set.seed(52107)
  design <- cbind(1, condition = rep(c(-0.5, 0.5), 8))
  effects <- rbind(level = c(1, 0), condition = c(0, 1))
  sources <- stats::setNames(lapply(seq_len(partitions), function(index) {
    matrix(rnorm(16 * 7), 16, 7)
  }), paste0("run", seq_len(partitions)))
  domain <- abstract_domain(7L, id = domain_id)
  fit <- lm_relation_fit(
    sources, design, effects, sampling_unit = "trial", domain = domain
  )
  frame <- compile_frame(whole_brain(), domain)
  list(fit = fit, frame = frame, domain = domain)
}

test_that("capability refusals are classed conditions with structured fields", {
  fixture <- refusal_fixture()
  plan <- plan_geometry(
    fixture$fit$relation, fixture$frame,
    cross_partitions(fixture$fit$relation, independence = "independent")
  )
  refusal <- catch_refusal(
    rdm_sampling_covariance(plan, fixture$fit$relation, target = "null")
  )
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "sampling_covariance")
  expect_identical(refusal$namespace, "evidence_sampling")
  expect_true("missing_error_channel" %in% refusal$reasons)
  expect_true(length(refusal$remedies) >= 1L)
  expect_match(refusal$message, "lm_relation_fit")
  expect_null(catch_refusal(rdm(plan)))
})

test_that("every unmet sampling requirement is reported, not only the first", {
  fixture <- refusal_fixture()
  partial <- plan_geometry(
    fixture$fit$relation, fixture$frame,
    pairing("run1", "run2", independence = "independent")
  )
  compiled <- crossform:::.compile_evidence_sampling_plan(
    partial, fixture$fit$relation, sampling_axis = "trial"
  )
  expect_true(all(c("missing_error_channel", "not_equal_weight_all_pairs")
    %in% compiled$unavailable_reasons))
  refusal <- catch_refusal(
    crossform:::.require_sampling_covariance(compiled)
  )
  expect_identical(refusal$reasons, compiled$unavailable_reasons)
  expect_match(refusal$message, "complete, equally weighted set")
  expect_match(refusal$message, "beta matrices alone")

  # Without a declared sampling axis, that reason joins the set and gets its
  # own sentence rather than the generic fallback.
  axis_free <- crossform:::.compile_evidence_sampling_plan(
    partial, fixture$fit$relation
  )
  expect_true("sampling_axis_missing_or_inconsistent" %in%
    axis_free$unavailable_reasons)
  axis_refusal <- catch_refusal(
    crossform:::.require_sampling_covariance(axis_free)
  )
  expect_match(axis_refusal$message, "no single sampling axis is declared")
})

test_that("undeclared endpoint independence refuses with its own sentence", {
  fixture <- refusal_fixture(partitions = 2L,
    domain_id = "refusal-independence-domain")
  dependent <- plan_geometry(
    fixture$fit$relation, fixture$frame,
    pairing("run1", "run2", independence = "not_independent")
  )
  compiled <- crossform:::.compile_evidence_sampling_plan(
    dependent, fixture$fit
  )
  expect_true("endpoint_independence_not_declared" %in%
    compiled$unavailable_reasons)
  refusal <- catch_refusal(
    crossform:::.require_sampling_covariance(compiled)
  )
  expect_match(refusal$message, "statistically\\s+independent")
  expect_true(any(grepl("independence", refusal$remedies)))
})

test_that("a learned-frozen metric is refused by admission, not by accident", {
  fixture <- refusal_fixture(domain_id = "refusal-learned-domain")
  learned_metric <- neural_metric(
    diag(7), fixture$domain,
    estimation = "learned_frozen",
    provenance = list(
      frozen = TRUE,
      training_signature = paste0("sha256:", strrep("ab", 32))
    )
  )
  plan <- plan_geometry(
    fixture$fit$relation, fixture$frame,
    cross_partitions(fixture$fit$relation, independence = "independent"),
    metric = learned_metric
  )
  descriptor <- crossform:::.sampling_evidence_descriptor(plan)
  expect_identical(descriptor$record$metric_status, "learned")

  refusal <- catch_refusal(
    rdm_sampling_covariance(plan, fixture$fit, target = "null")
  )
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "fixed_metric_sampling_law")
  expect_identical(refusal$reasons, "learned_metric_law_not_admitted")
  expect_match(refusal$message, "metric was learned")
  expect_no_match(refusal$message, "frozen training provenance")
})

test_that("an omitted calibration target refuses as a classed condition", {
  fixture <- refusal_fixture(domain_id = "refusal-target-domain")
  plan <- plan_geometry(
    fixture$fit$relation, fixture$frame,
    cross_partitions(fixture$fit$relation, independence = "independent")
  )
  refusal <- catch_refusal(rdm_sampling_covariance(plan, fixture$fit))
  expect_identical(refusal$capability, "calibration_target_declared")
  expect_match(refusal$message, "not inferred")
})

test_that("plan_crossnobis diagnoses a missing residual channel first", {
  fixture <- refusal_fixture(partitions = 2L,
    domain_id = "refusal-crossnobis-domain")
  bare_refusal <- catch_refusal(plan_crossnobis(
    fixture$fit$relation, fixture$frame,
    pairing("run1", "run2", independence = "independent")
  ))
  expect_s3_class(bare_refusal, "effect_capability_refusal")
  expect_identical(bare_refusal$capability, "learned_metric_input")
  expect_match(bare_refusal$message, "beta matrices alone.*lm_relation_fit")

  channel_free <- relation_fit(fixture$fit$relation)
  fit_refusal <- catch_refusal(plan_crossnobis(
    channel_free, fixture$frame,
    pairing("run1", "run2", independence = "independent")
  ))
  expect_s3_class(fit_refusal, "effect_capability_refusal")
  expect_match(fit_refusal$message, "beta matrices alone.*lm_relation_fit")
  expect_no_match(fit_refusal$message, "training")
})

test_that("relation-fit refusals carry the capability and partitions", {
  fixture <- refusal_fixture(domain_id = "refusal-relation-fit-domain")
  refusal <- catch_refusal(
    crossform:::.require_relation_fit_capability(
      fixture$fit$relation, "residual_blocks"
    )
  )
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "residual_blocks")
  expect_identical(refusal$namespace, "relation_fit")
  expect_identical(refusal$reasons, "missing_error_channel")
})

test_that("error capabilities are earned by their own fields", {
  fixture <- refusal_fixture(domain_id = "refusal-capability-bits-domain")
  model <- fixture$fit$error_models[[1L]]
  full <- crossform:::.error_capabilities(model)
  expect_true(full$within_participant_calibration)

  covariance_free <- model
  covariance_free$effect_covariance <- NULL
  partial <- crossform:::.error_capabilities(covariance_free)
  expect_true(partial$error_model)
  expect_false(partial$effect_covariance)
  expect_false(partial$within_participant_calibration)
  expect_true(partial$residual_blocks)
})

test_that("sampling capabilities can be asked for instead of provoked", {
  fixture <- refusal_fixture(domain_id = "refusal-introspection-domain")
  plan <- plan_geometry(
    fixture$fit$relation, fixture$frame,
    cross_partitions(fixture$fit$relation, independence = "independent")
  )
  granted <- sampling_capabilities(plan, fixture$fit)
  expect_s3_class(granted, "effect_sampling_capabilities")
  expect_true(granted$available)
  expect_identical(nrow(granted$reasons), 0L)

  bare <- sampling_capabilities(plan)
  expect_false(bare$available)
  expect_true("missing_error_channel" %in% bare$reasons$reason)
  expect_match(
    bare$reasons$remedy[bare$reasons$reason == "missing_error_channel"],
    "lm_relation_fit"
  )
  expect_output(print(bare), "unavailable")
})

test_that("the plan print surfaces metric status and generalization", {
  fixture <- refusal_fixture(domain_id = "refusal-print-domain")
  plan <- plan_geometry(
    fixture$fit$relation, fixture$frame,
    cross_partitions(
      fixture$fit$relation,
      independence = "independent",
      generalizes_over = "run"
    )
  )
  printed <- paste(utils::capture.output(print(plan)), collapse = "\n")
  expect_match(printed, "metric:")
  expect_match(printed, "over run")
  expect_match(printed, "independent")
})
