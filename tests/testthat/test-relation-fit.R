direct_gls_effects <- function(response, design, target, whitener = NULL) {
  if (is.null(whitener)) whitener <- diag(nrow(design))
  transformed_design <- whitener %*% design
  transformed_response <- whitener %*% response
  target %*% solve(
    crossprod(transformed_design),
    crossprod(transformed_design, transformed_response)
  )
}

direct_gls_residuals <- function(response, design, whitener = NULL) {
  if (is.null(whitener)) whitener <- diag(nrow(design))
  transformed_design <- whitener %*% design
  transformed_response <- whitener %*% response
  transformed_response - transformed_design %*% solve(
    crossprod(transformed_design),
    crossprod(transformed_design, transformed_response)
  )
}

test_that("raw and precomputed routes yield identical pure relation blocks", {
  set.seed(8101)
  design <- cbind(1, condition = rep(c(-0.5, 0.5), 8))
  target <- rbind(level = c(1, 0), condition = c(0, 1))
  whitener <- diag(seq(0.75, 1.25, length.out = nrow(design)))
  raw <- list(
    run1 = matrix(rnorm(16 * 9), 16, 9),
    run2 = matrix(rnorm(16 * 9), 16, 9)
  )
  fit <- lm_relation_fit(
    raw, design, target, observation_whitener = whitener
  )
  effects <- lapply(raw, direct_gls_effects,
    design = design, target = target, whitener = whitener)
  precomputed <- relation(effects, effects = c("level", "condition"))

  for (partition in names(raw)) {
    for (features in list(1:9, c(2, 5, 9))) {
      expect_equal(
        relation_block(fit, partition, features),
        relation_block(precomputed, partition, features),
        tolerance = 2e-13
      )
    }
  }
  expect_s3_class(fit$relation, "effect_relation")
  expect_false(inherits(fit$relation, "effect_relation_fit"))
})

test_that("residual blocks agree with a direct GLS residual oracle", {
  set.seed(8102)
  design <- cbind(1, linear = seq(-1, 1, length.out = 14))
  target <- rbind(intercept = c(1, 0), linear = c(0, 1))
  whitener <- diag(seq(0.6, 1.4, length.out = 14))
  response <- matrix(rnorm(14 * 11), 14, 11)
  fit <- lm_relation_fit(
    list(run = response), design, target,
    observation_whitener = whitener, sampling_unit = "time_point"
  )
  expected <- direct_gls_residuals(response, design, whitener)

  expect_equal(residual_block(fit, "run", c(1, 4, 11)),
    expected[, c(1, 4, 11)], tolerance = 2e-13)
  expect_equal(residual_block(fit, 1, 1:11), expected, tolerance = 2e-13)
  expect_identical(residual_df(fit, "run"), 12L)
  expect_identical(fit$error_models$run$sampling_unit, "time_point")
  expect_identical(fit$error_models$run$residual_source$dim, c(14L, 11L))
})

test_that("effect covariance has the declared separable scale convention", {
  design <- cbind(1, x = seq(-1, 1, length.out = 15),
    z = rep(c(-1, 0, 1), 5))
  target <- rbind(mean = c(1, 0, 0), x = c(0, 1, 0))
  whitener <- diag(seq(0.7, 1.3, length.out = 15))
  response <- matrix(seq_len(15 * 6), 15, 6)
  fit <- lm_relation_fit(
    list(run = response), design, target,
    observation_whitener = whitener
  )
  transformed <- whitener %*% design
  expected <- target %*% solve(crossprod(transformed)) %*% t(target)
  model <- fit$error_models$run

  expect_equal(effect_covariance(fit, "run"), expected, tolerance = 2e-13)
  expect_identical(model$scale_convention,
    "neural_covariance_factor_excluded")
  expect_identical(model$kind, "separable_glm")
  expect_setequal(model$assumptions, c(
    "matrix_normal_separable",
    "shared_observation_covariance_across_neural_features",
    "observation_whitener_treated_as_fixed"
  ))
  expect_match(model$source_revision, "^sha256:[[:xdigit:]]{64}$")
  expect_identical(
    model$source_revision,
    fit$relation$capabilities$run$stable_revision
  )
})

test_that("residual response sources remain lazy", {
  set.seed(8103)
  response <- matrix(rnorm(12 * 7), 12, 7)
  reads <- 0L
  source <- function(features) {
    reads <<- reads + 1L
    response[, features, drop = FALSE]
  }
  revision <- paste0("sha256:", paste(rep("a", 64), collapse = ""))
  capabilities <- source_capabilities(
    block_read = TRUE, stable_revision = revision
  )
  design <- cbind(1, x = seq_len(12))
  target <- matrix(c(0, 1), 1, dimnames = list("slope", NULL))
  fit <- lm_relation_fit(
    list(run = source), design, target,
    source_dims = list(c(12, 7)), capabilities = capabilities
  )

  expect_identical(reads, 0L)
  block <- residual_block(fit, "run", c(2, 7))
  expect_identical(reads, 1L)
  expect_identical(dim(block), c(12L, 2L))
})

test_that("pure relations advertise no recoverable error capability", {
  relation_value <- relation(
    list(run = matrix(1:12, 3, 4,
      dimnames = list(c("a", "b", "c"), NULL))),
    effects = c("a", "b", "c")
  )
  absent <- relation_fit_capabilities(relation_value)
  wrapped <- relation_fit(relation_value)

  expect_false(any(unlist(absent[-1L], use.names = FALSE)))
  expect_false(any(unlist(
    relation_fit_capabilities(wrapped)[-1L], use.names = FALSE
  )))
  expect_equal(relation_block(relation_value, "run", 1:2),
    relation_block(wrapped, "run", 1:2), tolerance = 0)
  expect_error(
    crossform:::.require_relation_fit_capability(
      relation_value, "learned_metric_input"
    ),
    "precomputed effects.*lm_relation_fit"
  )
  expect_error(residual_block(wrapped, "run", 1),
    "precomputed effects.*lm_relation_fit")
  expect_error(residual_block(relation_value, "run", 1),
    "beta matrices alone.*lm_relation_fit")
  expect_error(effect_covariance(relation_value, "run"),
    "beta matrices alone.*lm_relation_fit")
  expect_error(residual_df(relation_value, "run"),
    "beta matrices alone.*lm_relation_fit")
  expect_error(
    residual_pair_statistics(
      relation_value,
      compile_frame(searchlights(1), abstract_domain(4))
    ),
    "beta matrices alone.*lm_relation_fit"
  )
})

test_that("deep validation is a boundary operation, not a block-read rebuild", {
  design <- cbind(1, x = seq_len(10))
  target <- matrix(c(0, 1), 1, dimnames = list("slope", NULL))
  fit <- lm_relation_fit(
    list(run = matrix(seq_len(50), 10, 5)), design, target
  )
  mutated <- fit
  mutated$signature <- paste0("sha256:", paste(rep("f", 64), collapse = ""))

  expect_silent(residual_block(mutated, "run", 1:2))
  expect_error(
    crossform:::.validate_relation_fit(mutated, deep = TRUE),
    "identity"
  )
})

test_that("a saturated model cannot claim a residual error channel", {
  design <- diag(5)
  target <- matrix(1, 1, 5, dimnames = list("sum", NULL))
  response <- matrix(seq_len(20), 5, 4)

  expect_s3_class(lm_extractor(design, target), "effect_extractor")
  expect_error(
    lm_relation_fit(list(run = response), design, target),
    "residual degree of freedom"
  )
})

test_that("lm relation fits tolerate a forwarded NULL legacy whitener", {
  design <- cbind(1, x = seq_len(8))
  target <- matrix(c(0, 1), 1, dimnames = list("slope", NULL))
  response <- matrix(seq_len(32), 8, 4)

  expect_silent(lm_relation_fit(
    list(run = response), design, target, whiten = NULL
  ))
  expect_silent(lm_relation_fit(
    list(run = response), design, target,
    observation_whitener = diag(8), whiten = NULL
  ))
})

test_that("partition-specific models retain one explicit common effect space", {
  set.seed(8104)
  first_design <- cbind(1, x = seq(-1, 1, length.out = 10))
  second_design <- cbind(1, x = seq(-2, 2, length.out = 12),
    nuisance = rep(c(-1, 1), 6))
  first_target <- rbind(level = c(1, 0), slope = c(0, 1))
  second_target <- rbind(level = c(1, 0, 0), slope = c(0, 1, 0))
  response <- list(
    first = matrix(rnorm(10 * 5), 10, 5),
    second = matrix(rnorm(12 * 5), 12, 5)
  )
  fit <- lm_relation_fit(
    response,
    design = list(first = first_design, second = second_design),
    effects = list(first = first_target, second = second_target),
    sampling_unit = c("time_point", "time_point")
  )

  expect_identical(fit$relation$effects, c("level", "slope"))
  expect_identical(residual_df(fit, "first"), 8L)
  expect_identical(residual_df(fit, "second"), 9L)
  expect_equal(
    relation_block(fit, "second", 1:5),
    direct_gls_effects(
      response$second, second_design, second_target
    ),
    tolerance = 2e-13
  )
})

test_that("rank-deficient estimable fits residualize against the estimable span", {
  set.seed(8105)
  x <- seq_len(9)
  design <- cbind(x, duplicate = x)
  target <- matrix(c(1, 1), 1, dimnames = list("sum", NULL))
  response <- matrix(rnorm(9 * 4), 9, 4)
  fit <- lm_relation_fit(list(run = response), design, target)
  expected <- response - x %*% solve(crossprod(x), crossprod(x, response))

  expect_identical(residual_df(fit, "run"), 8L)
  expect_equal(residual_block(fit, "run", 1:4), expected,
    tolerance = 2e-13)
  expect_identical(
    fit$error_models$run$estimator_provenance$solver,
    "svd_estimable_fallback"
  )
})

test_that("non-estimable partition effects report rank, aliases, and remedies", {
  conditions <- rep(c("face", "house", "object", "body"), each = 3L)
  indicators <- stats::model.matrix(~ 0 + factor(conditions,
    levels = c("face", "house", "object", "body")))
  colnames(indicators) <- c("face", "house", "object", "body")
  design <- cbind(`(Intercept)` = 1, indicators)
  effects <- cbind(`(Intercept)` = 0, diag(4))
  rownames(effects) <- colnames(indicators)
  response <- matrix(seq_len(nrow(design) * 3L), nrow(design), 3L)

  refusal <- catch_refusal(lm_relation_fit(
    list(session1_run2 = response), design, effects
  ))

  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "estimable_effects")
  expect_match(conditionMessage(refusal), "Partition `session1_run2`")
  expect_match(conditionMessage(refusal), "rank 4 for 5 regressors")
  expect_match(conditionMessage(refusal), "\\(Intercept\\).+face.+house.+object.+body")
  expect_match(conditionMessage(refusal), "drop the intercept")
  expect_match(conditionMessage(refusal), "differences between conditions")
  expect_length(refusal$reasons, 3L)
  expect_length(refusal$remedies, 2L)
})
