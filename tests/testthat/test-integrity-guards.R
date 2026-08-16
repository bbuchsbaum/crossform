integrity_fixture <- function(domain = NULL,
                              domain_id = "integrity-guard-domain",
                              frame_spec = voxelwise()) {
  set.seed(90217)
  if (is.null(domain)) domain <- abstract_domain(6L, id = domain_id)
  features <- domain$n_features
  design <- cbind(1, condition = rep(c(-0.5, 0.5), 8))
  effects <- rbind(level = c(1, 0), condition = c(0, 1), extra = c(1, -1))
  sources <- stats::setNames(lapply(seq_len(3L), function(index) {
    matrix(rnorm(16 * features), 16, features)
  }), paste0("run", seq_len(3L)))
  fit <- lm_relation_fit(
    sources, design, effects, sampling_unit = "trial", domain = domain
  )
  plan <- plan_geometry(
    fit$relation, compile_frame(frame_spec, domain),
    cross_partitions(fit$relation, independence = "independent")
  )
  list(fit = fit, plan = plan, domain = domain)
}

test_that("correlation-style normalization refuses with the policy boundary", {
  fixture <- integrity_fixture()
  refusal <- catch_refusal(rdm(fixture$plan, normalize = "correlation"))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "guaranteed_psd")
  expect_match(refusal$message, "zero or negative")
  expect_match(refusal$message, "correlation-distance policy")
})

test_that("univariate removal refuses and points at the decomposition", {
  fixture <- integrity_fixture()
  weights <- c(level = 0, condition = 1, extra = 0)
  refusal <- catch_refusal(
    contrast_energy(fixture$plan, weights, remove_univariate = TRUE)
  )
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "nondestructive_decomposition")
  expect_match(refusal$message, "coherent")
  expect_match(refusal$message, "configuration")
})

test_that("rdm and rsa recompose exactly at the public surface", {
  fixture <- integrity_fixture()
  total <- rdm(fixture$plan, component = "total")$values
  coherent <- rdm(fixture$plan, component = "coherent")$values
  configuration <- rdm(fixture$plan, component = "configuration")$values
  expect_equal(total, coherent + configuration, tolerance = 1e-12)

  model <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), 3, 3,
    dimnames = list(
      c("level", "condition", "extra"), c("level", "condition", "extra")
    )
  )
  rsa_total <- rsa(fixture$plan, models = list(m = model))$coefficients
  rsa_coherent <- rsa(
    fixture$plan, models = list(m = model), component = "coherent"
  )$coefficients
  rsa_configuration <- rsa(
    fixture$plan, models = list(m = model), component = "configuration"
  )$coefficients
  expect_equal(rsa_total, rsa_coherent + rsa_configuration,
    tolerance = 1e-12)
})

test_that("conservative frames conserve total evidence; local frames do not", {
  mask <- array(TRUE, dim = c(3L, 2L, 1L))
  domain <- volume_domain(mask, id = "integrity-conservation-volume")
  conservative <- integrity_fixture(
    domain = domain, frame_spec = searchlights(1.01, "conservative")
  )
  weights <- c(level = 0, condition = 1, extra = 0)
  local_view <- contrast_energy(conservative$plan, weights)
  # The conservation comparator is the unnormalized whole-brain operator:
  # `whole_brain()`'s default local normalization divides by the feature
  # count, which is exactly the factor-of-P trap the docs warn about.
  global_plan <- plan_geometry(
    conservative$fit$relation,
    compile_frame(whole_brain("none"), domain),
    cross_partitions(
      conservative$fit$relation, independence = "independent"
    )
  )
  global_view <- contrast_energy(global_plan, weights)
  expect_equal(sum(local_view$total), global_view$total, tolerance = 1e-10)

  # The conservation law covers `total` only: coherent components are
  # defined by each node's local common mode and do not sum globally.
  expect_gt(
    abs(sum(local_view$coherent) - global_view$coherent), 1e-8
  )

  conservative_frame <- compile_frame(
    searchlights(1.01, "conservative"), domain
  )
  local_frame <- compile_frame(searchlights(1.01), domain)
  expect_true(frame_conservation(conservative_frame)$conserved)
  expect_false(frame_conservation(local_frame)$conserved)
})

test_that("unequal-length runs get a legible per-partition diagnosis", {
  domain <- abstract_domain(4L, id = "integrity-unequal-runs")
  design <- cbind(1, condition = rep(c(-0.5, 0.5), 8))
  effects <- rbind(condition = c(0, 1))
  sources <- list(
    run1 = matrix(rnorm(16 * 4), 16, 4),
    run2 = matrix(rnorm(12 * 4), 12, 4)
  )
  expect_error(
    lm_relation_fit(sources, design, effects, domain = domain),
    "Partition `run2` supplies 12 observations.*named list"
  )
  per_partition <- list(
    run1 = design,
    run2 = design[seq_len(12L), , drop = FALSE]
  )
  fit <- lm_relation_fit(sources, per_partition, effects, domain = domain)
  expect_s3_class(fit, "effect_relation_fit")
})

test_that("a covariance artifact cannot claim a plan of another dimension", {
  fixture <- integrity_fixture()
  plan <- crossform:::.compile_evidence_sampling_plan(
    fixture$plan, fixture$fit,
    target = crossform:::.sampling_target("null")
  )
  contrasts <- matrix(rnorm(5 * 120), 5, 120)
  expect_error(
    crossform:::.sampling_covariance_from_components(
      plan, contrasts,
      signal_patterns = matrix(0, 120, 4),
      effect_covariance = diag(120),
      residual_covariance = diag(4)
    ),
    "declare 120 experimental effects but the bound evidence plan has 3"
  )
})

test_that("inconsistent canonical blocks refuse instead of clipping", {
  cross <- diag(2) * 5
  left_self <- diag(2) * 0.1
  right_self <- diag(2) * 0.1
  expect_error(
    crossform:::.measurement_invariant_summary(
      cross, left_self, right_self
    ),
    "will not be silently clipped"
  )
})
