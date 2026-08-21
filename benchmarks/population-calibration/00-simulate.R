# Paired population-calibration generator and inference references.

population_calibration_scenarios <- function() {
  data.frame(
    scenario = c(
      "gaussian_homoskedastic", "covariate_heteroskedastic",
      "heavy_tailed_influential", "unequal_transport_quality_fixed",
      "unequal_transport_quality_cross_fitted",
      "independent_coverage", "informative_coverage_fixed",
      "informative_coverage_cross_fitted"
    ),
    transport_regime = c(
      "fixed", "fixed", "cross_fitted", "fixed", "cross_fitted",
      "cross_fitted", "fixed", "cross_fitted"
    ),
    coverage_regime = c(
      "complete", "complete", "complete", "complete", "complete",
      "independent", "informative", "informative"
    ),
    marginal_claim_supported = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                                 FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}

.population_calibration_dataset <- function(seed, scenario, n = 24L) {
  set.seed(seed)
  z <- as.numeric(scale(stats::rnorm(n)))
  latent <- stats::rnorm(n)
  quality_draw <- stats::runif(n)
  coverage_draw <- stats::runif(n)
  heavy <- stats::rt(n, df = 3) / sqrt(3)
  influence <- sample.int(n, 1L)

  error <- latent
  quality <- rep(1, n)
  available <- rep(TRUE, n)
  if (identical(scenario, "covariate_heteroskedastic")) {
    scale_i <- exp(0.8 * z)
    scale_i <- scale_i / sqrt(mean(scale_i^2))
    error <- latent * scale_i
  } else if (identical(scenario, "heavy_tailed_influential")) {
    error <- heavy
    error[[influence]] <- error[[influence]] + 5 * sign(z[[influence]] + 0.1)
  } else if (identical(scenario, "unequal_transport_quality_fixed")) {
    quality <- stats::plogis(1.2 - 0.9 * z + 0.35 * quality_draw)
    error <- latent / quality
  } else if (identical(scenario,
                       "unequal_transport_quality_cross_fitted")) {
    quality <- stats::plogis(1.2 + 0.7 * stats::qnorm(quality_draw))
    error <- latent / quality
  } else if (identical(scenario, "independent_coverage")) {
    available <- coverage_draw > 0.2
  } else if (identical(scenario, "informative_coverage_fixed")) {
    quality <- stats::plogis(0.8 - 0.8 * z + 0.2 * quality_draw)
    error <- latent / quality
    available <- coverage_draw < stats::plogis(0.5 + 0.9 * z + 1.1 * error)
  } else if (identical(scenario,
                       "informative_coverage_cross_fitted")) {
    quality <- stats::plogis(1.0 + 0.6 * stats::qnorm(quality_draw))
    error <- latent / quality
    available <- coverage_draw < stats::plogis(0.5 + 1.1 * error)
  }
  if (sum(available) < 6L || length(unique(z[available])) < 2L) {
    available[order(abs(error), decreasing = FALSE)[seq_len(6L)]] <- TRUE
  }
  data.frame(
    subject = sprintf("s%02d", seq_len(n)), z = z, outcome = error,
    transport_quality = quality, available = available,
    stringsAsFactors = FALSE
  )
}

.population_calibration_fit <- function(data, weights, level = 0.95) {
  coverage_fraction <- mean(data$available)
  data <- data[data$available, , drop = FALSE]
  X <- cbind(`(Intercept)` = 1, z = data$z)
  y <- data$outcome
  n <- nrow(X)
  p <- ncol(X)
  if (n <= p || qr(X)$rank < p) return(NULL)
  bread <- solve(crossprod(X))
  beta <- drop(bread %*% crossprod(X, y))
  residual <- y - drop(X %*% beta)
  leverage <- rowSums((X %*% bread) * X)
  df <- n - p
  classical_cov <- sum(residual^2) / df * bread
  adjusted <- residual / (1 - leverage)
  hc3_cov <- bread %*% crossprod(X, X * adjusted^2) %*% bread
  ses <- c(classical = sqrt(classical_cov[2L, 2L]),
           HC3 = sqrt(hc3_cov[2L, 2L]))
  critical <- stats::qt(1 - (1 - level) / 2, df = df)

  # Null-imposed wild bootstrap for H0: beta_z = 0. The same weight columns
  # are used for this dataset by every cell and never regenerated per method.
  restricted <- mean(y)
  restricted_residual <- y - restricted
  restricted_leverage <- rep(1 / n, n)
  boot_y <- restricted +
    (restricted_residual / (1 - restricted_leverage)) * weights[seq_len(n), ]
  boot_beta <- bread %*% crossprod(X, boot_y)
  boot_residual <- boot_y - X %*% boot_beta
  influence <- drop(X %*% bread[, 2L]) / (1 - leverage)
  boot_se <- sqrt(drop(crossprod(influence^2, boot_residual^2)))
  boot_t <- boot_beta[2L, ] / boot_se
  observed_t <- beta[[2L]] / ses[["HC3"]]
  valid <- is.finite(boot_t)
  p_value <- if (any(valid)) {
    (1 + sum(abs(boot_t[valid]) >= abs(observed_t))) / (1 + sum(valid))
  } else NA_real_
  boot_critical <- if (any(valid)) {
    unname(stats::quantile(abs(boot_t[valid]), level, names = FALSE))
  } else NA_real_

  list(
    estimate = beta[[2L]], n = n, df = df,
    classical_se = ses[["classical"]], hc3_se = ses[["HC3"]],
    critical = critical, wild_p = p_value, wild_critical = boot_critical,
    wild_valid = sum(valid), quality_mean = mean(data$transport_quality),
    coverage_fraction = coverage_fraction
  )
}

population_calibration_replicates <- function(
    replications = 500L, bootstrap_replicates = 399L,
    seed = 73001L, n = 24L, level = 0.95) {
  scenarios <- population_calibration_scenarios()
  output <- list()
  for (scenario_index in seq_len(nrow(scenarios))) {
    scenario <- scenarios[scenario_index, ]
    for (replication in seq_len(replications)) {
      dataset_seed <- seed + scenario_index * 100000L + replication
      data <- .population_calibration_dataset(
        dataset_seed, scenario$scenario, n
      )
      set.seed(dataset_seed + 50000000L)
      weights <- matrix(ifelse(
        stats::runif(n * bootstrap_replicates) < 0.5, -1, 1
      ), n, bootstrap_replicates)
      fit <- .population_calibration_fit(data, weights, level)
      dataset_id <- sprintf("%s-r%04d-s%d", scenario$scenario,
                            replication, dataset_seed)
      if (is.null(fit)) {
        for (method in c("classical", "HC3", "wild_bootstrap")) {
          output[[length(output) + 1L]] <- data.frame(
            dataset_id = dataset_id, scenario = scenario$scenario,
            transport_regime = scenario$transport_regime,
            coverage_regime = scenario$coverage_regime,
            marginal_claim_supported = scenario$marginal_claim_supported,
            replication = replication, method = method, n = sum(data$available),
            estimate = NA_real_, reported_se = NA_real_, lower = NA_real_,
            upper = NA_real_, reject = NA, quality_mean = NA_real_,
            coverage_fraction = mean(data$available), failure = TRUE
          )
        }
        next
      }
      intervals <- list(
        classical = c(fit$estimate - fit$critical * fit$classical_se,
                      fit$estimate + fit$critical * fit$classical_se),
        HC3 = c(fit$estimate - fit$critical * fit$hc3_se,
                fit$estimate + fit$critical * fit$hc3_se),
        wild_bootstrap = c(fit$estimate - fit$wild_critical * fit$hc3_se,
                           fit$estimate + fit$wild_critical * fit$hc3_se)
      )
      for (method in names(intervals)) {
        reject <- if (identical(method, "wild_bootstrap")) {
          is.finite(fit$wild_p) && fit$wild_p <= 1 - level
        } else !(intervals[[method]][[1L]] <= 0 &&
                   intervals[[method]][[2L]] >= 0)
        output[[length(output) + 1L]] <- data.frame(
          dataset_id = dataset_id, scenario = scenario$scenario,
          transport_regime = scenario$transport_regime,
          coverage_regime = scenario$coverage_regime,
          marginal_claim_supported = scenario$marginal_claim_supported,
          replication = replication, method = method, n = fit$n,
          estimate = fit$estimate,
          reported_se = if (identical(method, "classical"))
            fit$classical_se else fit$hc3_se,
          lower = intervals[[method]][[1L]], upper = intervals[[method]][[2L]],
          reject = reject, quality_mean = fit$quality_mean,
          coverage_fraction = fit$coverage_fraction,
          failure = !is.finite(intervals[[method]][[1L]])
        )
      }
    }
  }
  result <- do.call(rbind, output)
  rownames(result) <- NULL
  result
}
