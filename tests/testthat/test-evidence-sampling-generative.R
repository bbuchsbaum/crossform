sampling_calibration_cache <- new.env(parent = emptyenv())

sampling_calibration_fixture <- function(signal = TRUE, repetitions = 4000L,
                                         seed = if (isTRUE(signal)) {
                                           81343
                                         } else {
                                           81342
                                         }) {
  cache_key <- paste(signal, repetitions, seed, sep = ":")
  if (exists(cache_key, envir = sampling_calibration_cache,
      inherits = FALSE)) {
    return(get(cache_key, envir = sampling_calibration_cache,
      inherits = FALSE))
  }
  conditions <- 3L
  features <- 40L
  partitions <- 6L
  sigma_k <- diag(0.8, conditions)
  sigma_r <- diag(features)
  set.seed(81340)
  patterns <- matrix(rnorm(conditions * features), conditions, features)
  if (isTRUE(signal)) {
    components <- sampling_oracle_components(patterns, sigma_k, sigma_r)
    patterns <- patterns * sqrt(0.5 / mean(diag(components$delta)))
  } else {
    patterns[] <- 0
  }
  value <- sampling_oracle_calibration_experiment(
    patterns, sigma_k, sigma_r, partitions, repetitions, seed,
    plugin = isTRUE(signal)
  )
  assign(cache_key, value, envir = sampling_calibration_cache)
  value
}

test_that("Eq. 13 recovers full sampling covariance and linear transports", {
  experiment <- sampling_calibration_fixture(signal = TRUE)
  empirical <- stats::cov(experiment$estimates)
  relative_error <- norm(empirical - experiment$covariance, "F") /
    norm(experiment$covariance, "F")
  diagonal_ratio <- diag(empirical) / diag(experiment$covariance)
  linear <- c(1, -0.5, 0.25)
  empirical_linear <- stats::var(drop(experiment$estimates %*% linear))
  analytic_linear <- drop(
    t(linear) %*% experiment$covariance %*% linear
  )

  expect_lt(relative_error, 0.06)
  expect_true(all(diagonal_ratio > 0.92 & diagonal_ratio < 1.08))
  expect_equal(empirical_linear / analytic_linear, 1, tolerance = 0.08)
  expect_lt(relative_error,
    norm(diag(diag(empirical)) - experiment$covariance, "F") /
      norm(experiment$covariance, "F"))
})

test_that("null and plug-in targets have explicitly different calibration", {
  null <- sampling_calibration_fixture(signal = FALSE)
  signal <- sampling_calibration_fixture(signal = TRUE)
  null_se <- sqrt(diag(null$covariance))
  null_coverage <- colMeans(abs(null$estimates) <=
    1.96 * matrix(null_se, nrow(null$estimates), length(null_se),
      byrow = TRUE))
  plugin_coverage <- colMeans(abs(sweep(
    signal$estimates, 2L, signal$truth
  )) <= 1.96 * signal$plugin_se)

  expect_true(all(null_coverage > 0.92 & null_coverage < 0.98))
  # The partition-mean plug-in target is intentionally reported as a distinct,
  # mildly conservative policy rather than silently called exact 95% coverage.
  expect_true(all(plugin_coverage > 0.95 & plugin_coverage < 0.995))
})

test_that("spread across dependent partition edges is not a valid signal SE", {
  signal <- sampling_calibration_fixture(signal = TRUE)
  analytic_se <- sqrt(diag(signal$covariance))
  analytic_coverage <- colMeans(abs(sweep(
    signal$estimates, 2L, signal$truth
  )) <= 1.96 * matrix(analytic_se, nrow(signal$estimates),
    length(analytic_se), byrow = TRUE))
  naive_coverage <- colMeans(abs(sweep(
    signal$estimates, 2L, signal$truth
  )) <= 1.96 * signal$naive_se)

  expect_true(all(analytic_coverage > 0.92 & analytic_coverage < 0.98))
  expect_true(all(naive_coverage < 0.82))
  expect_true(all(analytic_coverage - naive_coverage > 0.13))
})
