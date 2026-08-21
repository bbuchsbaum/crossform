# Null-imposed subject-level wild bootstrap (PE-D4).

pb_reference_hc3 <- function(design, response, contrast) {
  bread <- solve(crossprod(design))
  coefficient <- as.numeric(bread %*% crossprod(design, response))
  fitted <- as.numeric(design %*% coefficient)
  residual <- response - fitted
  leverage <- diag(design %*% bread %*% t(design))
  adjusted <- residual / (1 - leverage)
  covariance <- bread %*%
    crossprod(design, design * adjusted^2) %*% bread
  covariance <- (covariance + t(covariance)) / 2
  list(
    coefficient = coefficient,
    estimate = as.numeric(crossprod(contrast, coefficient)),
    fitted = fitted,
    residual = residual,
    leverage = leverage,
    adjusted = adjusted,
    bread = bread,
    covariance = covariance,
    se = sqrt(as.numeric(crossprod(contrast, covariance %*% contrast)))
  )
}

pb_slow_reference <- function(fit, node, query, contrast, null, weights) {
  available <- fit$coverage$availability[node, query, ]
  design <- fit$uncertainty$between$design[available, , drop = FALSE]
  response <- fit$values[node, query, available]
  observed <- pb_reference_hc3(design, response, contrast)
  denominator <- as.numeric(crossprod(
    contrast, observed$bread %*% contrast
  ))
  restricted <- observed$coefficient -
    as.numeric(observed$bread %*% contrast) *
    (observed$estimate - null) / denominator
  restricted_fitted <- as.numeric(design %*% restricted)
  restricted_adjusted <- (response - restricted_fitted) /
    (1 - observed$leverage)
  statistic <- vapply(seq_len(ncol(weights)), function(replicate) {
    response_star <- restricted_fitted + restricted_adjusted *
      weights[available, replicate]
    fit_star <- pb_reference_hc3(design, response_star, contrast)
    (fit_star$estimate - null) / fit_star$se
  }, numeric(1))
  observed_t <- (observed$estimate - null) / observed$se
  list(
    observed = observed,
    observed_t = observed_t,
    statistic = statistic,
    p_value = (1 + sum(abs(statistic) >= abs(observed_t))) /
      (1 + length(statistic))
  )
}

pb_fixture <- function(subjects = 8L, seed = 901L) {
  set.seed(seed)
  labels <- sprintf("s%02d", seq_len(subjects))
  covariates <- data.frame(
    x1 = stats::rnorm(subjects),
    x2 = stats::runif(subjects, -1, 1),
    row.names = labels
  )
  gains <- stats::setNames(exp(stats::rnorm(subjects, sd = 0.35)), labels)
  ph_fit(covariates, gains)
}

test_that("seeded bootstrap agrees with a slow explicit reference", {
  fit <- pb_fixture()
  contrast <- c("(Intercept)" = 0, x1 = 1, x2 = -0.5)
  null <- 0.15
  bootstrap <- population_wild_bootstrap(
    fit, contrast, null = null, replicates = 99L, seed = 1771L
  )
  reference1 <- pb_slow_reference(
    fit, "group1", "face-house", contrast, null, bootstrap$weights
  )
  reference2 <- pb_slow_reference(
    fit, "group2", "face-house", contrast, null, bootstrap$weights
  )

  expect_identical(bootstrap$contrast, contrast)
  expect_identical(bootstrap$null, null)
  expect_equal(bootstrap$observed_t["group1", "face-house"],
    reference1$observed_t, tolerance = 2e-12)
  expect_equal(unname(bootstrap$replicate_t[
    "group1", "face-house", ]), reference1$statistic, tolerance = 2e-11)
  expect_equal(unname(bootstrap$replicate_t[
    "group2", "face-house", ]), reference2$statistic, tolerance = 2e-11)
  expect_identical(bootstrap$p_value["group1", "face-house"],
    reference1$p_value)
  expect_identical(bootstrap$p_value["group2", "face-house"],
    reference2$p_value)
  ordinary <- c("group1", "group2", "group3")
  expect_true(all(bootstrap$successful_replicates[ordinary, ] == 99L))
  expect_true(all(bootstrap$failed_replicates[ordinary, ] == 0L))
  expect_identical(bootstrap$reason["<sink>", "face-house"],
    "nonpositive_standard_error")
})

test_that("one subject weight matrix is reused and fully identified", {
  fit <- pb_fixture(seed = 902L)
  before <- .Random.seed
  first <- population_wild_bootstrap(
    fit, "x1", replicates = 99L, seed = 1881L
  )
  expect_identical(.Random.seed, before)
  second <- population_wild_bootstrap(
    fit, "x1", replicates = 99L, seed = 1881L
  )

  expect_identical(first$weights, second$weights)
  expect_identical(first$replicate_t, second$replicate_t)
  expect_identical(first$scientific_plan_id, second$scientific_plan_id)
  expect_identical(rownames(first$weights), fit$coverage$planned_subjects)
  expect_identical(dim(first$weights), c(8L, 99L))
  expect_setequal(unique(as.numeric(first$weights)), c(-1, 1))
  expect_identical(first$weight_signature,
    crossform:::.sha256_signature(first$weights))

  # There is one weight for subject s03 in replicate 17. Both cells in that
  # replicate read this exact stored number; no node- or query-specific weight
  # array exists anywhere on the record.
  expect_length(first$weights["s03", "r17"], 1L)
  expect_false(any(grepl("node", names(first), fixed = TRUE) &
    grepl("weight", names(first), fixed = TRUE)))
})

test_that("local coverage subsets the shared weights without redrawing", {
  set.seed(907L)
  labels <- sprintf("s%02d", 1:8)
  covariates <- data.frame(
    x1 = stats::rnorm(8), x2 = stats::runif(8, -1, 1), row.names = labels
  )
  gains <- stats::setNames(exp(stats::rnorm(8, sd = 0.25)), labels)
  subjects <- stats::setNames(lapply(labels, function(id) {
    ph_subject(id, 10L, gains[[id]])
  }), labels)
  group_index <- c("group1", "group2")
  carriers <- stats::setNames(lapply(seq_along(labels), function(position) {
    column <- if (position <= 2L) rep(1L, 10L) else c(2L, rep(1L, 9L))
    operator <- Matrix::sparseMatrix(
      i = seq_len(10L), j = column, x = 1, dims = c(10L, 2L)
    )
    location_transport(operator,
      native_index = paste0("f", seq_len(10L)),
      group_index = group_index, semantics = "budget",
      provenance = list(method = "external", details = "coverage court"))
  }), labels)
  fit <- estimate_population(plan_population(
    subjects, carriers, model = ~ x1 + x2, data = covariates,
    coverage_policy = "available_at_node"
  ), rbind(`face-house` = c(1, -1)))
  value <- population_wild_bootstrap(
    fit, "x1", replicates = 99L, seed = 2221L
  )
  contrast <- c("(Intercept)" = 0, x1 = 1, x2 = 0)
  reference <- pb_slow_reference(
    fit, "group2", "face-house", contrast, 0, value$weights
  )

  expect_identical(fit$coverage$planned_subjects,
    rownames(value$weights))
  expect_identical(fit$coverage$n["group2", "face-house"], 6L)
  expect_identical(fit$coverage$subject_sets[[
    fit$coverage$subject_set_id["group2", "face-house"]
  ]], labels[3:8])
  expect_equal(unname(value$replicate_t[
    "group2", "face-house", ]), reference$statistic, tolerance = 2e-11)
  expect_identical(value$p_value["group2", "face-house"],
    reference$p_value)
})

test_that("Mammen weights have the declared support and moments", {
  fit <- pb_fixture(seed = 903L)
  value <- population_wild_bootstrap(
    fit, "x2", replicates = 999L, seed = 1991L, weights = "mammen"
  )
  root <- sqrt(5)
  support <- c(-(root - 1) / 2, (root + 1) / 2)

  expect_identical(value$weight_distribution, "mammen")
  expect_equal(sort(unique(as.numeric(value$weights))), sort(support),
    tolerance = 0)
  expect_lt(abs(mean(value$weights)), 0.04)
  expect_lt(abs(mean(value$weights^2) - 1), 0.04)
  expect_lt(abs(mean(value$weights^3) - 1), 0.08)
})

test_that("insufficient cells and zero residuals retain failure provenance", {
  small <- pb_fixture(subjects = 5L, seed = 904L)
  refused <- population_wild_bootstrap(
    small, "x1", replicates = 99L, seed = 2001L
  )
  expect_true(all(refused$status == "refused"))
  expect_true(all(refused$reason ==
    "insufficient_subjects_for_wild_bootstrap"))
  expect_true(all(refused$failed_replicates == 99L))
  expect_true(all(refused$replicate_failure_reason ==
    "insufficient_subjects_for_wild_bootstrap"))
  expect_true(all(is.na(refused$p_value)))
  expect_true(all(is.na(refused$reject)))

  labels <- sprintf("s%02d", 1:6)
  subjects <- stats::setNames(lapply(labels, function(id) {
    ph_subject(id, 10L, 1)
  }), labels)
  carrier <- anatomical_transport(
    native_coords = cbind(0:9), group_coords = cbind(c(0, 4, 9)),
    semantics = "budget"
  )
  zero_fit <- estimate_population(plan_population(
    subjects, stats::setNames(rep(list(carrier), length(labels)), labels)
  ), rbind(`face-house` = c(1, -1)))
  zero <- population_wild_bootstrap(
    zero_fit, "(Intercept)", replicates = 99L, seed = 2002L
  )
  expect_true(all(zero$status == "refused"))
  expect_true(all(zero$reason == "nonpositive_standard_error"))
  expect_true(all(zero$replicate_failure_reason ==
    "nonpositive_standard_error"))
})

test_that("Monte Carlo error, tables, printing, and tamper checks are explicit", {
  fit <- pb_fixture(seed = 905L)
  value <- population_wild_bootstrap(
    fit, "x1", null = 0.2, replicates = 199L, seed = 2111L
  )
  expected_mcse <- sqrt(value$p_value * (1 - value$p_value) /
    (1 + value$successful_replicates))
  expect_equal(value$monte_carlo_se, expected_mcse, tolerance = 0)
  expect_true(all(value$p_value[is.finite(value$p_value)] >= 1 / 200))

  table <- as.data.frame(value)
  expect_true(all(c(
    "contrast", "null", "observed_t", "successful_replicates",
    "failed_replicates", "p_value", "monte_carlo_se", "critical_value",
    "status", "reason", "available_subjects", "seed", "weight_signature"
  ) %in% names(table)))
  expect_identical(unique(table$seed), 2111L)
  expect_output(print(value), "one subject weight reused")
  expect_output(print(value), "Monte Carlo SE")

  broken <- value
  broken$weights[1L, 1L] <- -broken$weights[1L, 1L]
  expect_error(crossform:::.validate_population_bootstrap(broken),
    "noncanonical", class = "effect_input_error")

  broken <- value
  broken$contrast[["x1"]] <- broken$contrast[["x1"]] + 1
  expect_error(crossform:::.validate_population_bootstrap(broken),
    "identity", class = "effect_contract_error")
})

test_that("bootstrap arguments name invalid requests", {
  fit <- pb_fixture(seed = 906L)
  expect_error(population_wild_bootstrap(
    fit, "missing", replicates = 99L, seed = 1L
  ), "contrast")
  expect_error(population_wild_bootstrap(
    fit, c(0, 0, 0), replicates = 99L, seed = 1L
  ), "nonzero")
  expect_error(population_wild_bootstrap(
    fit, "x1", replicates = 98L, seed = 1L
  ), "at least 99")
  expect_error(population_wild_bootstrap(
    fit, "x1", replicates = 99L
  ), "seed")
})
