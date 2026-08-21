#!/usr/bin/env Rscript

population_calibration_summary <- function(replicates) {
  groups <- split(seq_len(nrow(replicates)), interaction(
    replicates$scenario, replicates$method, drop = TRUE, lex.order = TRUE
  ))
  rows <- lapply(groups, function(index) {
    x <- replicates[index, ]
    ok <- !x$failure & is.finite(x$estimate)
    estimates <- x$estimate[ok]
    coverage <- mean(x$lower[ok] <= 0 & x$upper[ok] >= 0)
    rejection <- mean(x$reject[ok])
    data.frame(
      scenario = x$scenario[[1L]], method = x$method[[1L]],
      transport_regime = x$transport_regime[[1L]],
      coverage_regime = x$coverage_regime[[1L]],
      marginal_claim_supported = x$marginal_claim_supported[[1L]],
      replications = nrow(x), successful = sum(ok),
      bias = mean(estimates), empirical_se = stats::sd(estimates),
      mean_reported_se = mean(x$reported_se[ok]), coverage = coverage,
      coverage_mcse = sqrt(coverage * (1 - coverage) / sum(ok)),
      rejection = rejection,
      rejection_mcse = sqrt(rejection * (1 - rejection) / sum(ok)),
      failure_rate = mean(!ok), mean_n = mean(x$n),
      mean_coverage_fraction = mean(x$coverage_fraction),
      mean_transport_quality = mean(x$quality_mean, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

population_calibration_verdicts <- function(summary) {
  supported <- summary$marginal_claim_supported
  hc3_supported <- summary$method == "HC3" & supported
  gaussian <- summary$scenario == "gaussian_homoskedastic"
  hetero <- summary$scenario == "covariate_heteroskedastic"
  informative <- summary$coverage_regime == "informative"
  rows <- list(
    data.frame(
      gate = "paired_methods", observed = min(table(
        summary$scenario, summary$method
      )), threshold = 1, comparison = ">=", passes =
        all(table(summary$scenario, summary$method) == 1),
      boundary = "each scenario has classical, HC3, and wild bootstrap"
    ),
    data.frame(
      gate = "gaussian_nominal_coverage",
      observed = max(abs(summary$coverage[gaussian] - 0.95)),
      threshold = 0.04, comparison = "<=",
      passes = max(abs(summary$coverage[gaussian] - 0.95)) <= 0.04,
      boundary = "absolute coverage error across all three methods"
    ),
    data.frame(
      gate = "hc3_supported_coverage",
      observed = min(summary$coverage[hc3_supported] -
                       2 * summary$coverage_mcse[hc3_supported]),
      threshold = 0.88, comparison = ">=",
      passes = min(summary$coverage[hc3_supported] -
                     2 * summary$coverage_mcse[hc3_supported]) >= 0.88,
      boundary = "two-MCSE lower bound in declared marginally supported regimes"
    ),
    data.frame(
      gate = "heteroskedastic_hc3_improves_classical",
      observed = summary$coverage[hetero & summary$method == "HC3"] -
        summary$coverage[hetero & summary$method == "classical"],
      threshold = 0.02, comparison = ">=",
      passes = summary$coverage[hetero & summary$method == "HC3"] -
        summary$coverage[hetero & summary$method == "classical"] >= 0.02,
      boundary = "paired covariate-linked heteroskedastic arm"
    ),
    data.frame(
      gate = "informative_coverage_refusal",
      observed = sum(informative & !summary$marginal_claim_supported),
      threshold = sum(informative), comparison = "==",
      passes = all(!summary$marginal_claim_supported[informative]),
      boundary = "no implemented interval is licensed for a marginal claim"
    ),
    data.frame(
      gate = "transport_regime_inventory",
      observed = length(unique(summary$transport_regime)), threshold = 2,
      comparison = "==",
      passes = setequal(summary$transport_regime, c("fixed", "cross_fitted")),
      boundary = "fixed and cross-fitted transport are both explicit"
    ),
    data.frame(
      gate = "failure_rate", observed = max(summary$failure_rate),
      threshold = 0.02, comparison = "<=",
      passes = max(summary$failure_rate) <= 0.02,
      boundary = "numerical or rank failure across certification cells"
    )
  )
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

population_calibration_parameters <- function() data.frame(
  key = c("schema_version", "seed", "replications", "subjects",
          "bootstrap_replicates", "level", "target_beta",
          "coverage_tolerance", "hc3_two_mcse_floor"),
  value = c("population-calibration-v1", "73001", "500", "24", "399",
            "0.95", "0", "0.04", "0.88"), stringsAsFactors = FALSE
)

population_calibration_checksums <- function(repo) {
  paths <- c(
    "benchmarks/population-calibration/00-simulate.R",
    "benchmarks/population-calibration/01-certify.R",
    "design/population-calibration-contract.md",
    "inst/extdata/certification/evidence-status-ledger.csv",
    "inst/extdata/certification/population-calibration-parameters.csv",
    "inst/extdata/certification/population-calibration-replicates.csv",
    "inst/extdata/certification/population-calibration-results.csv",
    "inst/extdata/certification/population-calibration-verdicts.csv"
  )
  full <- file.path(repo, paths)
  data.frame(
    schema_version = "population-calibration-checksums-v1",
    path = paths, size_bytes = unname(file.info(full)$size),
    hash_algorithm = "md5", digest = unname(tools::md5sum(full)),
    stringsAsFactors = FALSE
  )
}

if (sys.nframe() == 0L) {
  script <- sub("^--file=", "",
                grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  directory <- normalizePath(dirname(script))
  repo <- normalizePath(file.path(directory, "..", ".."))
  source(file.path(directory, "00-simulate.R"))
  replicates <- population_calibration_replicates()
  summary <- population_calibration_summary(replicates)
  verdicts <- population_calibration_verdicts(summary)
  parameters <- population_calibration_parameters()
  output <- file.path(repo, "inst", "extdata", "certification")
  utils::write.csv(parameters,
    file.path(output, "population-calibration-parameters.csv"), row.names = FALSE)
  utils::write.csv(replicates,
    file.path(output, "population-calibration-replicates.csv"), row.names = FALSE)
  utils::write.csv(summary,
    file.path(output, "population-calibration-results.csv"), row.names = FALSE)
  utils::write.csv(verdicts,
    file.path(output, "population-calibration-verdicts.csv"), row.names = FALSE)
  checksums <- population_calibration_checksums(repo)
  utils::write.csv(checksums,
    file.path(output, "population-calibration-checksums.csv"), row.names = FALSE)
  print(verdicts, row.names = FALSE)
  if (!all(verdicts$passes)) {
    stop("Population calibration failed ", sum(!verdicts$passes), " gate(s).",
         call. = FALSE)
  }
  message("Population calibration PASS: ", nrow(replicates),
          " paired method rows.")
}
