# HC3 is tested through the public population verbs and against arithmetic
# written independently here. The court never calls the package's bread,
# leverage, sandwich, or uncertainty helpers.

test_that("intercept-only HC3 has its exact closed form", {
  labels <- sprintf("s%02d", 1:7)
  covariates <- data.frame(
    x1 = seq_len(7), x2 = rev(seq_len(7)), row.names = labels
  )
  gains <- stats::setNames(c(0.7, 1.2, 0.9, 1.5, 0.8, 1.1, 1.8), labels)
  fit <- ph_fit(covariates, gains, model = ~ 1)
  robust <- population_uncertainty(fit, estimator = "HC3")
  response <- fit$values["group1", "face-house", ]
  residual <- response - mean(response)
  hand_se <- sqrt(sum((residual / (1 - 1 / length(response)))^2) /
    length(response)^2)

  expect_equal(robust$between$se[
    "group1", "face-house", "(Intercept)"
  ], hand_se, tolerance = 1e-12)
  expect_equal(unname(robust$between$leverage[
    "group1", "face-house", ]), rep(1 / length(response), length(response)),
    tolerance = 1e-14)
})

ph_hc3_reference <- function(design, response) {
  bread <- solve(crossprod(design))
  coefficient <- as.numeric(bread %*% crossprod(design, response))
  residual <- response - as.numeric(design %*% coefficient)
  leverage <- diag(design %*% bread %*% t(design))
  meat <- crossprod(design,
    design * (residual / (1 - leverage))^2)
  covariance <- bread %*% meat %*% bread
  list(
    covariance = (covariance + t(covariance)) / 2,
    leverage = leverage
  )
}

test_that("HC3 agrees with an independent randomized full-rank court", {
  for (seed in 101:105) {
    set.seed(seed)
    labels <- sprintf("s%02d", 1:10)
    covariates <- data.frame(
      x1 = stats::rnorm(10), x2 = stats::runif(10, -1, 1),
      row.names = labels
    )
    gains <- stats::setNames(exp(stats::rnorm(10, sd = 0.35)), labels)
    fit <- ph_fit(covariates, gains)
    robust <- population_uncertainty(fit, estimator = "HC3")
    design <- cbind("(Intercept)" = 1,
      x1 = covariates$x1, x2 = covariates$x2)

    expect_identical(robust$estimator, "HC3")
    expect_true("heteroskedasticity_robust_sandwich" %in%
      robust$assumptions)
    for (node in c("group1", "group2", "group3")) {
      response <- fit$values[node, "face-house", ]
      reference <- ph_hc3_reference(design, response)
      expect_equal(
        unname(robust$between$covariance[node, "face-house", , ]),
        unname(reference$covariance), tolerance = 2e-11
      )
      expect_equal(
        unname(robust$between$leverage[node, "face-house", ]),
        reference$leverage, tolerance = 2e-12
      )
      expect_equal(
        unname(robust$between$se[node, "face-house", ]),
        unname(sqrt(diag(reference$covariance))), tolerance = 2e-11
      )
    }
  }
})

test_that("a high-leverage small sample separates HC3 from classical OLS", {
  labels <- sprintf("s%02d", 1:8)
  covariates <- data.frame(
    x1 = c(-1, -0.7, -0.3, 0, 0.2, 0.5, 0.9, 3.5),
    x2 = c(0.3, -0.8, 0.7, -0.4, 0.1, 0.9, -0.2, 0.5),
    row.names = labels
  )
  gains <- stats::setNames(c(0.7, 1.1, 0.9, 1.2, 0.8, 1.3, 0.75, 2.8),
    labels)
  fit <- ph_fit(covariates, gains)
  classical <- population_uncertainty(fit, estimator = "classical")
  robust <- population_uncertainty(fit, estimator = "HC3")
  cell <- c("group1", "face-house", "x1")

  expect_gt(robust$between$max_leverage["group1", "face-house"], 0.75)
  expect_gt(robust$between$se[cell[[1L]], cell[[2L]], cell[[3L]]],
    classical$between$se[cell[[1L]], cell[[2L]], cell[[3L]]])
  expect_false(identical(robust$scientific_plan_id,
    classical$scientific_plan_id))
})

test_that("HC3 refuses leverage at one instead of dividing by zero", {
  labels <- sprintf("s%02d", 1:8)
  covariates <- data.frame(
    x1 = c(rep(0, 7), 1e9),
    x2 = c(-3, -2, -1, 0, 1, 2, 3, 0),
    row.names = labels
  )
  gains <- stats::setNames(seq(0.7, 1.4, length.out = 8), labels)
  fit <- ph_fit(covariates, gains)
  robust <- population_uncertainty(fit, estimator = "HC3")

  expect_equal(robust$between$max_leverage["group1", "face-house"], 1,
    tolerance = 1e-12)
  expect_identical(robust$between$status["group1", "face-house"],
    "refused")
  expect_identical(robust$between$reason["group1", "face-house"],
    "leverage_near_one")
  expect_true(all(is.na(robust$between$se["group1", "face-house", ])))
})

test_that("uncertainty tables retain estimator, assumptions, and diagnostics", {
  set.seed(211)
  labels <- sprintf("s%02d", 1:9)
  covariates <- data.frame(
    x1 = stats::rnorm(9), x2 = stats::runif(9), row.names = labels
  )
  gains <- stats::setNames(exp(stats::rnorm(9, sd = 0.2)), labels)
  robust <- population_uncertainty(ph_fit(covariates, gains), estimator = "HC3")
  table <- as.data.frame(robust)

  expect_true(all(c(
    "estimator", "assumptions", "n", "design_rank", "residual_df",
    "max_leverage", "max_abs_adjusted_residual", "uncertainty_status",
    "uncertainty_reason",
    "transport_conditioning", "coverage_conditioning"
  ) %in% names(table)))
  expect_identical(unique(table$estimator), "HC3")
  expect_true(all(table$transport_conditioning ==
    "conditional_on_realized_transport"))
  expect_output(print(robust), "HC3 SE")
  expect_output(print(robust), "max h")
})

test_that("the sealed design recipe detects tampering", {
  set.seed(307)
  labels <- sprintf("s%02d", 1:8)
  covariates <- data.frame(
    x1 = stats::rnorm(8), x2 = stats::runif(8), row.names = labels
  )
  gains <- stats::setNames(exp(stats::rnorm(8, sd = 0.2)), labels)
  fit <- ph_fit(covariates, gains)
  broken <- fit
  broken$uncertainty$between$design[1L, "x1"] <-
    broken$uncertainty$between$design[1L, "x1"] + 1

  expect_error(crossform:::.validate_population_result(broken),
    "recipe identity", class = "effect_contract_error")
})
