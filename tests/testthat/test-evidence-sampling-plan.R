sampling_plan_fixture <- function(partitions = 3L,
                                  domain_id = "sampling-plan-domain") {
  set.seed(81301)
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
  evidence <- plan_geometry(
    fit$relation, frame,
    cross_partitions(fit$relation, independence = "independent")
  )
  list(fit = fit, evidence = evidence)
}

test_that("sampling plans bind evidence, error, metric, and target identities", {
  fixture <- sampling_plan_fixture()
  plan <- crossform:::.compile_evidence_sampling_plan(
    fixture$evidence, fixture$fit
  )

  expect_s3_class(plan, "effect_evidence_sampling_plan")
  expect_identical(plan$contract, "evidence-sampling-v1")
  expect_identical(plan$capabilities$sampling_covariance, "available")
  expect_identical(plan$capabilities$error_channel, "relation_fit")
  expect_identical(plan$capabilities$error_model, "separable_glm")
  expect_identical(plan$capabilities$metric_role, "metric")
  expect_identical(plan$capabilities$metric_status, "fixed")
  expect_identical(plan$capabilities$metric_uncertainty, "not_applicable")
  expect_identical(plan$capabilities$partition_model, "equal")
  expect_identical(plan$capabilities$sampling_axis, "trial")
  expect_identical(plan$capabilities$calibration_target, "point_estimate")
  expect_identical(plan$capabilities$spatial_covariance, "local_marginals")
  expect_identical(plan$capabilities$requested_operation, "diagonal")
  expect_identical(plan$capabilities$materialization, "query_only")
  expect_false(plan$partition$edge_products_independent)
  expect_silent(crossform:::.require_sampling_covariance(plan))
})

test_that("same-shaped sampling requests keep distinct semantic identities", {
  fixture <- sampling_plan_fixture()
  baseline <- crossform:::.compile_evidence_sampling_plan(
    fixture$evidence, fixture$fit
  )
  null_target <- crossform:::.compile_evidence_sampling_plan(
    fixture$evidence, fixture$fit,
    target = crossform:::.sampling_target("null")
  )
  selected <- crossform:::.compile_evidence_sampling_plan(
    fixture$evidence, fixture$fit,
    operation = crossform:::.sampling_operation(
      "selected_entries", matrix(c(1L, 1L), 1L, 2L)
    )
  )
  heterogeneous <- crossform:::.compile_evidence_sampling_plan(
    fixture$evidence, fixture$fit, partition_model = "heterogeneous"
  )

  identities <- vapply(
    list(baseline, null_target, selected, heterogeneous),
    `[[`, character(1), "scientific_plan_id"
  )
  expect_length(unique(identities), 4L)
  expect_identical(
    heterogeneous$capabilities$sampling_covariance, "unavailable"
  )
  expect_true("heterogeneous_partition_model" %in%
    heterogeneous$unavailable_reasons)
})

test_that("a bare relation cannot gain a sampling law from its shape", {
  fixture <- sampling_plan_fixture()
  fitted <- crossform:::.compile_evidence_sampling_plan(
    fixture$evidence, fixture$fit
  )
  bare <- crossform:::.compile_evidence_sampling_plan(
    fixture$evidence, fixture$fit$relation, sampling_axis = "trial"
  )

  expect_identical(fitted$evidence$relation_id, bare$evidence$relation_id)
  expect_identical(bare$capabilities$error_channel, "absent")
  expect_identical(bare$capabilities$sampling_covariance, "unavailable")
  expect_true("missing_error_channel" %in% bare$unavailable_reasons)
  expect_false(identical(
    fitted$scientific_plan_id, bare$scientific_plan_id
  ))
  refusal <- catch_refusal(crossform:::.require_sampling_covariance(bare))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "sampling_covariance")
  expect_identical(refusal$namespace, "evidence_sampling")
  expect_identical(refusal$reasons, "missing_error_channel")
  expect_match(conditionMessage(refusal),
    "lm_relation_fit.*beta matrices alone")
})

test_that("error channels are identity-bound to the evidence relation", {
  first <- sampling_plan_fixture()
  second <- sampling_plan_fixture(domain_id = "different-domain")

  expect_error(
    crossform:::.compile_evidence_sampling_plan(
      first$evidence, second$fit
    ),
    "identity-bound"
  )
})

test_that("learned metric uncertainty and spatial scope are explicit gates", {
  metric <- crossform:::.sampling_metric_record(
    "sha256:learned-metric", status = "learned", uncertainty = "ignored"
  )
  bridge <- crossform:::.sampling_metric_record(
    "sha256:bridge", status = "fixed", role = "bridge"
  )

  expect_identical(metric$status, "learned")
  expect_identical(metric$uncertainty, "ignored")
  expect_identical(bridge$role, "bridge")
  expect_error(
    crossform:::.sampling_metric_record(
      "sha256:learned-metric", status = "learned"
    ),
    "must declare uncertainty"
  )

  fixture <- sampling_plan_fixture()
  spatial <- crossform:::.compile_evidence_sampling_plan(
    fixture$evidence, fixture$fit, spatial_scope = "modeled"
  )
  expect_identical(spatial$capabilities$spatial_covariance, "unavailable")
  spatial_refusal <- catch_refusal(
    crossform:::.require_sampling_covariance(spatial)
  )
  expect_s3_class(spatial_refusal, "effect_capability_refusal")
  expect_identical(spatial_refusal$capability, "sampling_covariance")
  expect_identical(spatial_refusal$namespace, "evidence_sampling")
  expect_identical(spatial_refusal$reasons,
    "spatial_covariance_model_unavailable")
  expect_match(conditionMessage(spatial_refusal),
    "cross-location sampling covariance requires an explicit spatial model")
})

test_that("sampling-plan mutation is detected by canonical validation", {
  fixture <- sampling_plan_fixture()
  plan <- crossform:::.compile_evidence_sampling_plan(
    fixture$evidence, fixture$fit
  )
  forged <- plan
  forged$capabilities$calibration_target <- "null"

  expect_error(
    crossform:::.validate_evidence_sampling_plan(forged),
    "capabilities"
  )
})
