#!/usr/bin/env Rscript

# Reproducible statistical validation for Diedrichsen et al. Eq. 13 under the
# exact equal-partition, fixed-metric, separable matrix-normal capability.
# This script deliberately uses the independent base-R law court rather than
# crossform implementation code.

arguments <- commandArgs(trailingOnly = TRUE)
repetitions <- as.integer(Sys.getenv("CROSSFORM_CALIBRATION_REPS", "10000"))
if (length(arguments) > 0L && nzchar(arguments[[1L]])) {
  repetitions <- as.integer(arguments[[1L]])
}
if (is.na(repetitions) || repetitions < 1000L) {
  stop("Validation requires at least 1,000 Monte Carlo repetitions.")
}

source(file.path(
  "tests", "testthat", "helper-evidence-sampling-laws.R"
), local = TRUE)

conditions <- 3L
features <- 40L
partitions <- 6L
sigma_k <- diag(0.8, conditions)
sigma_r <- diag(features)
set.seed(81340)
patterns <- matrix(rnorm(conditions * features), conditions, features)
components <- sampling_oracle_components(patterns, sigma_k, sigma_r)
patterns <- patterns * sqrt(0.5 / mean(diag(components$delta)))

signal <- sampling_oracle_calibration_experiment(
  patterns, sigma_k, sigma_r, partitions, repetitions, seed = 81343,
  plugin = TRUE
)
null <- sampling_oracle_calibration_experiment(
  0 * patterns, sigma_k, sigma_r, partitions, repetitions, seed = 81342,
  plugin = FALSE
)

summarize_arm <- function(experiment, arm) {
  empirical <- stats::cov(experiment$estimates)
  analytic_se <- sqrt(diag(experiment$covariance))
  errors <- abs(sweep(experiment$estimates, 2L, experiment$truth))
  data.frame(
    arm = arm,
    coordinate = rownames(experiment$contrasts),
    truth = experiment$truth,
    empirical_variance = diag(empirical),
    analytic_variance = diag(experiment$covariance),
    variance_ratio = diag(empirical) / diag(experiment$covariance),
    analytic_coverage = colMeans(errors <= 1.96 * matrix(
      analytic_se, nrow(experiment$estimates), length(analytic_se),
      byrow = TRUE
    )),
    naive_edge_coverage = colMeans(errors <= 1.96 * experiment$naive_se),
    plugin_coverage = if (is.null(experiment$plugin_se)) {
      NA_real_
    } else {
      colMeans(errors <= 1.96 * experiment$plugin_se)
    },
    stringsAsFactors = FALSE
  )
}

summary <- rbind(
  summarize_arm(null, "null"),
  summarize_arm(signal, "nonzero_signal")
)
linear <- c(1, -0.5, 0.25)
linear_truth <- sum(linear * signal$truth)
linear_estimates <- drop(signal$estimates %*% linear)
linear_variance <- drop(t(linear) %*% signal$covariance %*% linear)
diagnostics <- list(
  schema_version = 1L,
  model = "fixed_metric_equal_partition_separable_matrix_normal",
  equation = "Diedrichsen_et_al_2016_Eq13",
  repetitions = repetitions,
  conditions = conditions,
  features = features,
  partitions = partitions,
  signal_full_covariance_relative_frobenius_error =
    norm(stats::cov(signal$estimates) - signal$covariance, "F") /
      norm(signal$covariance, "F"),
  null_full_covariance_relative_frobenius_error =
    norm(stats::cov(null$estimates) - null$covariance, "F") /
      norm(null$covariance, "F"),
  linear_variance_ratio = stats::var(linear_estimates) / linear_variance,
  linear_coverage = mean(abs(linear_estimates - linear_truth) <=
    1.96 * sqrt(linear_variance)),
  signal = signal,
  null = null,
  summary = summary
)

dir.create("benchmark-results", showWarnings = FALSE, recursive = TRUE)
saveRDS(diagnostics,
  file.path("benchmark-results", "sampling-covariance-validation.rds"),
  compress = "xz"
)
utils::write.csv(summary,
  file.path("benchmark-results", "sampling-covariance-validation-summary.csv"),
  row.names = FALSE
)

print(summary, row.names = FALSE)
cat(sprintf(
  paste0("signal covariance relative error: %.4f\n",
    "null covariance relative error: %.4f\n",
    "linear variance ratio: %.4f\n",
    "linear coverage: %.4f\n"),
  diagnostics$signal_full_covariance_relative_frobenius_error,
  diagnostics$null_full_covariance_relative_frobenius_error,
  diagnostics$linear_variance_ratio,
  diagnostics$linear_coverage
))
