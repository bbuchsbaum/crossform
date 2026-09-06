# Standalone oracle court for PE-D3. The oracle is source-visible and contains
# no call to a crossform internal helper; these tests bind it to the public
# population result only at the comparison boundary.

ph_oracle_path <- testthat::test_path(
  "..", "..", "design", "oracles", "population-hc3.R"
)
if (!file.exists(ph_oracle_path)) {
  testthat::skip("The source-tree HC3 oracle is not installed with the package.")
}
source(ph_oracle_path, local = TRUE)

expect_conditioned_equal <- function(actual, expected, condition_number) {
  scale <- max(1, abs(expected))
  relative_error <- max(abs(actual - expected)) / scale
  tolerance <- min(1e-3, max(
    5e-11, 100 * .Machine$double.eps * condition_number^2
  ))
  expect_lte(relative_error, tolerance)
}

test_that("the explicit oracle matches production over randomized designs", {
  for (seed in 401:405) {
    set.seed(seed)
    labels <- sprintf("s%02d", 1:11)
    x1 <- stats::rnorm(11)
    covariates <- data.frame(
      x1 = x1,
      x2 = if (seed == 405) x1 + stats::rnorm(11, sd = 1e-4) else
        stats::runif(11, -1, 1),
      row.names = labels
    )
    gains <- stats::setNames(exp(stats::rnorm(11, sd = 0.4)), labels)
    fit <- ph_fit(covariates, gains)
    production <- population_uncertainty(fit, estimator = "HC3")$between
    design <- cbind("(Intercept)" = 1,
      x1 = covariates$x1, x2 = covariates$x2)

    for (node in c("group1", "group2", "group3")) {
      response <- fit$values[node, "face-house", ]
      oracle <- population_hc3_oracle(design, response)
      expect_conditioned_equal(
        unname(fit$coefficients[node, "face-house", ]),
        unname(oracle$coefficient), oracle$condition_number
      )
      expect_conditioned_equal(
        unname(production$leverage[node, "face-house", ]),
        unname(oracle$leverage), oracle$condition_number
      )
      expect_conditioned_equal(
        unname(production$adjusted_residual[node, "face-house", ]),
        unname(oracle$adjusted_residual), oracle$condition_number
      )
      expect_conditioned_equal(
        unname(production$covariance[node, "face-house", , ]),
        unname(oracle$covariance), oracle$condition_number
      )
    }
  }
})

test_that("hand-computable intercept and coefficient subsets agree", {
  labels <- sprintf("s%02d", 1:7)
  covariates <- data.frame(
    x1 = c(-2, -1, 0, 0.5, 1, 2, 4),
    x2 = c(1, -1, 0.5, -0.5, 2, -2, 0), row.names = labels
  )
  gains <- stats::setNames(c(0.8, 1.3, 0.9, 1.6, 0.7, 1.2, 2.1), labels)

  mean_fit <- ph_fit(covariates, gains, model = ~ 1)
  mean_production <- population_uncertainty(
    mean_fit, estimator = "HC3"
  )$between
  response <- mean_fit$values["group1", "face-house", ]
  residual <- response - mean(response)
  hand_adjusted <- residual / (1 - 1 / length(response))
  hand_variance <- sum(hand_adjusted^2) / length(response)^2
  expect_equal(mean_production$covariance[
    "group1", "face-house", "(Intercept)", "(Intercept)"
  ], hand_variance, tolerance = 1e-12)

  fit <- ph_fit(covariates, gains)
  subset <- population_uncertainty(
    fit, term = "x1", estimator = "HC3"
  )$between
  design <- cbind("(Intercept)" = 1,
    x1 = covariates$x1, x2 = covariates$x2)
  oracle <- population_hc3_oracle(
    design, fit$values["group2", "face-house", ], terms = "x1"
  )
  expect_identical(dim(subset$covariance), c(4L, 1L, 1L, 1L))
  expect_equal(subset$covariance[
    "group2", "face-house", "x1", "x1"
  ], oracle$covariance["x1", "x1"], tolerance = 1e-12)
})

test_that("small n and zero residuals have documented finite behavior", {
  small_design <- cbind(
    "(Intercept)" = 1,
    x1 = c(-1, 0, 1, 2),
    x2 = c(0, 1, -1, 2)
  )
  small <- population_hc3_oracle(small_design, c(1, 2, 2.5, 5))
  expect_identical(small$residual_df, 1L)
  expect_true(all(is.finite(small$covariance)))
  expect_true(max(small$leverage) < 1)

  zero_response <- as.numeric(small_design %*% c(1, 0.5, -0.25))
  zero <- population_hc3_oracle(small_design, zero_response)
  expect_equal(zero$residual, rep(0, 4), tolerance = 1e-12)
  expect_equal(zero$adjusted_residual, rep(0, 4), tolerance = 1e-11)
  expect_equal(zero$covariance, matrix(0, 3, 3,
    dimnames = list(colnames(small_design), colnames(small_design))),
    tolerance = 1e-20)

  labels <- sprintf("s%02d", 1:6)
  subjects <- stats::setNames(lapply(labels, function(id) {
    ph_subject(id, 10L, 1)
  }), labels)
  carrier <- anatomical_transport(
    native_coords = cbind(0:9), group_coords = cbind(c(0, 4, 9)),
    semantics = "budget"
  )
  fit <- estimate_population(plan_population(
    subjects, stats::setNames(rep(list(carrier), length(labels)), labels)
  ), rbind(`face-house` = c(1, -1)))
  production <- population_uncertainty(fit, estimator = "HC3")$between
  expect_lt(max(abs(production$adjusted_residual), na.rm = TRUE), 1e-12)
  expect_lt(max(abs(production$covariance), na.rm = TRUE), 1e-22)
  expect_true(all(production$status == "degenerate_residual_scale"))
})

test_that("collinearity and leverage failures are deterministic", {
  x <- c(-2, -1, 0, 1, 2, 3)
  collinear <- cbind("(Intercept)" = 1, x1 = x, x2 = 2 * x)
  expect_error(population_hc3_oracle(collinear, seq_along(x)),
    "rank-deficient or numerically singular")
  expect_error(population_hc3_oracle(diag(3), 1:3),
    "positive residual degrees of freedom")

  leverage_one <- cbind(
    "(Intercept)" = 1,
    x1 = c(rep(0, 6), 1),
    x2 = c(-3, -2, -1, 0, 1, 2, 0)
  )
  expect_error(population_hc3_oracle(leverage_one, 1:7),
    "leverage whose complement is numerically zero")

  labels <- sprintf("s%02d", 1:7)
  covariates <- data.frame(x1 = x[c(1:6, 6)],
    x2 = 2 * x[c(1:6, 6)], row.names = labels)
  gains <- stats::setNames(seq(0.8, 1.4, length.out = 7), labels)
  refusal <- catch_refusal(ph_fit(covariates, gains))
  expect_identical(refusal$capability, "identified_group_model")
  expect_true("model_rank_deficient" %in% refusal$reasons)
})

test_that("the court detects a seeded HC2-for-HC3 defect", {
  design <- cbind(
    "(Intercept)" = 1,
    x1 = c(-1, -0.7, -0.2, 0, 0.3, 0.8, 1.1, 3.5),
    x2 = c(0.4, -0.9, 0.6, -0.3, 0.1, 1, -0.2, 0.7)
  )
  response <- c(0.5, 1.1, 0.8, 1.4, 0.7, 1.6, 1, 5.8)
  oracle <- population_hc3_oracle(design, response)

  # Seed the plausible defect: HC2 divides squared residuals by (1-h), while
  # HC3 divides by (1-h)^2. Everything else is intentionally identical.
  hc2_meat <- crossprod(design,
    design * (oracle$residual^2 / (1 - oracle$leverage)))
  mutant <- oracle$bread %*% hc2_meat %*% oracle$bread
  gap <- max(abs(mutant - oracle$covariance))

  expect_gt(gap, 1e-3)
  expect_false(isTRUE(all.equal(mutant, oracle$covariance,
    tolerance = 1e-8)))
})
