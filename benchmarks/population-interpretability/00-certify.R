#!/usr/bin/env Rscript

population_interpretability_truth <- function(repo) {
  x <- utils::read.csv(file.path(repo, "inst/extdata/certification",
    "matched-interpretability-truth.csv"), stringsAsFactors = FALSE)
  x[x$snr == 0.8 & is.finite(x$scale) &
      x$metric %in% c("total", "coherent_share"), ]
}

population_interpretability_fit <- function(X, y) {
  bread <- solve(crossprod(X))
  beta <- drop(bread %*% crossprod(X, y))
  residual <- y - drop(X %*% beta)
  leverage <- rowSums((X %*% bread) * X)
  adjusted <- residual / (1 - leverage)
  covariance <- bread %*% crossprod(X, X * adjusted^2) %*% bread
  se <- sqrt(diag(covariance))
  critical <- stats::qt(0.975, nrow(X) - ncol(X))
  list(beta = beta, se = se, lower = beta - critical * se,
       upper = beta + critical * se)
}

population_interpretability_replicates <- function(
    repo, replications = 200L, subjects = 24L, seed = 91001L) {
  truth <- population_interpretability_truth(repo)
  scenarios <- unique(truth$scenario)
  scales <- sort(unique(truth$scale))
  rows <- list()
  for (replication in seq_len(replications)) {
    set.seed(seed + replication)
    ids <- sprintf("s%02d", seq_len(subjects))
    x <- as.numeric(scale(stats::rnorm(subjects)))
    magnitude <- exp(0.15 * x + 0.12 * stats::rnorm(subjects))
    magnitude[[1L]] <- magnitude[[1L]] * 2.5
    quality <- stats::plogis(1.2 - 0.45 * x + 0.35 * stats::rnorm(subjects))
    organization_noise <- stats::rnorm(subjects, sd = 0.018)
    coverage_draw <- stats::runif(subjects)
    noise <- array(stats::rnorm(subjects * length(scales) * 2L, sd = 0.08),
                   c(subjects, length(scales), 2L))
    for (scale_index in seq_along(scales)) for (scenario in scenarios) {
      scale <- scales[[scale_index]]
      share <- truth$value[truth$scenario == scenario &
        truth$scale == scale & truth$metric == "coherent_share"]
      planted_total <- truth$value[truth$scenario == scenario &
        truth$scale == scale & truth$metric == "total"]
      subject_share <- pmin(1, pmax(0, share + organization_noise))
      latent <- list()
      latent$total <- planted_total * magnitude * quality
      latent$coherent <- latent$total * subject_share
      latent$configuration <- latent$total - latent$coherent
      observed <- list(
        coherent = latent$coherent + noise[, scale_index, 1L],
        configuration = latent$configuration + noise[, scale_index, 2L]
      )
      observed$total <- observed$coherent + observed$configuration
      for (coverage_regime in c("supported_complete", "failure_informative")) {
        available <- if (coverage_regime == "supported_complete") {
          rep(TRUE, subjects)
        } else {
          score <- as.numeric(scale(latent$coherent))
          coverage_draw < stats::plogis(0.25 + 0.55 * x + 0.75 * score)
        }
        if (sum(available) < 8L) {
          available[order(abs(x))[seq_len(8L)]] <- TRUE
        }
        X_all <- cbind(1, x)
        X <- X_all[available, , drop = FALSE]
        for (component in c("total", "coherent", "configuration")) {
          fit <- population_interpretability_fit(X, observed[[component]][available])
          conditional <- drop(solve(crossprod(X)) %*%
            crossprod(X, latent[[component]][available]))
          marginal <- drop(solve(crossprod(X_all)) %*%
            crossprod(X_all, latent[[component]]))
          for (term in seq_len(2L)) rows[[length(rows) + 1L]] <- data.frame(
            replication = replication, scenario = scenario, scale = scale,
            coverage_regime = coverage_regime, component = component,
            term = c("(Intercept)", "covariate_x")[[term]],
            estimate = fit$beta[[term]], se = fit$se[[term]],
            lower = fit$lower[[term]], upper = fit$upper[[term]],
            conditional_truth = conditional[[term]],
            marginal_truth = marginal[[term]], n = sum(available),
            coverage_fraction = mean(available),
            mean_transport_quality = mean(quality[available]),
            subject_set = paste(ids[available], collapse = ","),
            stringsAsFactors = FALSE)
        }
        if (scale == max(scales) && coverage_regime == "supported_complete") {
          X_all <- cbind(1, x)
          full <- population_interpretability_fit(X_all, observed$total)$beta[[1L]]
          delta <- vapply(seq_len(subjects), function(left) {
            keep <- seq_len(subjects) != left
            population_interpretability_fit(
              X_all[keep, , drop = FALSE], observed$total[keep]
            )$beta[[1L]] - full
          }, numeric(1))
          rows[[length(rows) + 1L]] <- data.frame(
            replication = replication, scenario = scenario, scale = scale,
            coverage_regime = "influence_supported", component = "total",
            term = "top_influence", estimate = which.max(abs(delta)),
            se = NA_real_, lower = NA_real_, upper = NA_real_,
            conditional_truth = 1, marginal_truth = 1, n = subjects,
            coverage_fraction = 1, mean_transport_quality = mean(quality),
            subject_set = ids[[which.max(abs(delta))]], stringsAsFactors = FALSE)
        }
      }
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

population_interpretability_results <- function(raw) {
  x <- raw[raw$term != "top_influence", ]
  groups <- split(seq_len(nrow(x)), interaction(x$scenario, x$scale,
    x$coverage_regime, x$component, x$term, drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(i) {
    z <- x[i, ]
    covered <- z$lower <= z$conditional_truth & z$conditional_truth <= z$upper
    data.frame(scenario = z$scenario[[1L]], scale = z$scale[[1L]],
      coverage_regime = z$coverage_regime[[1L]],
      component = z$component[[1L]], term = z$term[[1L]],
      replications = nrow(z), estimate_mean = mean(z$estimate),
      conditional_truth_mean = mean(z$conditional_truth),
      marginal_truth_mean = mean(z$marginal_truth),
      marginal_bias = mean(z$estimate - z$marginal_truth),
      target_shift = mean(z$conditional_truth - z$marginal_truth),
      empirical_se = stats::sd(z$estimate), mean_se = mean(z$se),
      conditional_coverage = mean(covered),
      coverage_mcse = sqrt(mean(covered) * (1 - mean(covered)) / nrow(z)),
      mean_n = mean(z$n), mean_coverage_fraction = mean(z$coverage_fraction),
      mean_transport_quality = mean(z$mean_transport_quality),
      marginal_claim_supported = z$coverage_regime[[1L]] == "supported_complete",
      stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

population_interpretability_verdicts <- function(raw, results, repo) {
  supported <- results$coverage_regime == "supported_complete" &
    results$term == "(Intercept)"
  part <- results[supported & results$component %in% c("coherent", "total"), ]
  wide <- reshape(part[, c("scenario", "scale", "component", "estimate_mean")],
    idvar = c("scenario", "scale"), timevar = "component", direction = "wide")
  wide$share <- wide$estimate_mean.coherent / wide$estimate_mean.total
  nonpoint <- wide$scale > min(wide$scale)
  ordering <- all(vapply(sort(unique(wide$scale[nonpoint])), function(s) {
    z <- wide[wide$scale == s, ]
    value <- z$share[match(c("broad_coherent", "mixed_broad_fine",
                             "fine_configuration"), z$scenario)]
    all(diff(value) < 0)
  }, logical(1)))
  total <- results[supported & results$component == "total", ]
  total_spread <- max(vapply(split(total$estimate_mean, total$scale),
    function(value) diff(range(value)), numeric(1)))
  coefficient <- raw[raw$term != "top_influence", ]
  key <- paste(coefficient$replication, coefficient$scenario, coefficient$scale,
               coefficient$coverage_regime, coefficient$term)
  conservation <- max(vapply(split(seq_len(nrow(coefficient)), key), function(i) {
    value <- stats::setNames(coefficient$estimate[i], coefficient$component[i])
    abs(value[["total"]] - value[["coherent"]] - value[["configuration"]])
  }, numeric(1)))
  cover <- results$coverage_regime == "supported_complete"
  coverage_lower <- min(results$conditional_coverage[cover] -
    2 * results$coverage_mcse[cover])
  planted <- population_interpretability_truth(repo)
  planted <- planted[planted$metric == "coherent_share" &
    planted$scale > min(planted$scale), ]
  recovered <- wide[wide$scale > min(wide$scale), ]
  profile <- merge(recovered, planted, by = c("scenario", "scale"))
  profile_correlation <- stats::cor(profile$share, profile$value)
  influence <- raw[raw$term == "top_influence", ]
  influence_rate <- mean(influence$subject_set == "s01")
  informative <- results$coverage_regime == "failure_informative" &
    results$term == "(Intercept)"
  target_shift <- max(abs(results$target_shift[informative]))
  observed <- c(ordering, total_spread, conservation, coverage_lower,
                profile_correlation, influence_rate, target_shift)
  threshold <- c(1, 0.08, 1e-10, 0.88, 0.98, 0.65, 0.05)
  comparison <- c("==", "<=", "<=", ">=", ">=", ">=", ">=")
  passes <- c(ordering, total_spread <= threshold[[2L]],
    conservation <= threshold[[3L]], coverage_lower >= threshold[[4L]],
    profile_correlation >= threshold[[5L]], influence_rate >= threshold[[6L]],
    target_shift >= threshold[[7L]])
  data.frame(gate = c("population_component_ordering", "matched_total_effects",
    "coefficient_conservation", "supported_hc3_coverage",
    "scale_profile_recovery", "injected_influence_recovery",
    "informative_coverage_limit_retained"),
    observed = observed, threshold = threshold, comparison = comparison,
    passes = passes, boundary = c(
      "broad greater than mixture greater than fine at nonpoint scales",
      "paired scenarios retain matched total coefficients",
      "total equals coherent plus configuration",
      "two-MCSE lower pointwise conditional coverage",
      "estimated versus planted coherent-share scale profile",
      "planted s01 magnitude outlier is top leave-one-out influence",
      "informative coverage shifts target and remains unsupported"),
    stringsAsFactors = FALSE)
}

population_interpretability_checksums <- function(repo) {
  paths <- c("benchmarks/population-interpretability/00-certify.R",
    "design/population-interpretability-contract.md",
    "inst/extdata/certification/matched-interpretability-truth.csv",
    "inst/extdata/certification/population-interpretability-replicates.csv",
    "inst/extdata/certification/population-interpretability-results.csv",
    "inst/extdata/certification/population-interpretability-verdicts.csv")
  full <- file.path(repo, paths)
  data.frame(schema_version = "population-interpretability-checksums-v1",
    path = paths, size_bytes = unname(file.info(full)$size),
    hash_algorithm = "md5", digest = unname(tools::md5sum(full)),
    stringsAsFactors = FALSE)
}

if (sys.nframe() == 0L) {
  script <- sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  repo <- normalizePath(file.path(dirname(script), "..", ".."))
  raw <- population_interpretability_replicates(repo)
  results <- population_interpretability_results(raw)
  verdicts <- population_interpretability_verdicts(raw, results, repo)
  output <- file.path(repo, "inst/extdata/certification")
  utils::write.csv(raw, file.path(output,
    "population-interpretability-replicates.csv"), row.names = FALSE)
  utils::write.csv(results, file.path(output,
    "population-interpretability-results.csv"), row.names = FALSE)
  utils::write.csv(verdicts, file.path(output,
    "population-interpretability-verdicts.csv"), row.names = FALSE)
  utils::write.csv(population_interpretability_checksums(repo), file.path(output,
    "population-interpretability-checksums.csv"), row.names = FALSE)
  print(verdicts, row.names = FALSE)
  if (!all(verdicts$passes)) stop("Population interpretability certification failed.")
  message("Population interpretability certification PASS: ", nrow(raw), " rows.")
}
