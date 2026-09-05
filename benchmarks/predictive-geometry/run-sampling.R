# Rscript benchmarks/predictive-geometry/run-sampling.R pilot|run
mode <- commandArgs(TRUE)
stopifnot(length(mode) == 1L, mode %in% c("pilot", "run"))
pkgload::load_all(".", quiet = TRUE)
source("benchmarks/provenance.R")
source("benchmarks/predictive-geometry/calibration.R")
source("benchmarks/predictive-geometry/sampling-config.R")
config <- predictive_sampling_config
config_hash <- digest::digest(config, algo = "sha256")
harness <- vapply(c("benchmarks/predictive-geometry/calibration.R",
  "benchmarks/predictive-geometry/sampling-config.R", "benchmarks/predictive-geometry/run-sampling.R"),
  function(path) digest::digest(file = path, algo = "sha256"), "")
source_before <- .crossform_source_tree_digest()
pilot_path <- "benchmark-results/predictive-geometry-sampling-pilot.rds"
results <- list()
if (mode == "run") {
  pilot <- readRDS(pilot_path)
  stopifnot(identical(pilot$config_hash, config_hash), identical(pilot$harness, harness))
}
for (i in seq_along(config$arms)) {
  arm <- config$arms[[i]]
  count <- if (mode == "pilot") config$pilot_replicates else pilot$counts[[arm]]
  seed <- (if (mode == "pilot") config$pilot_seed else config$production_seed) + 100L * i
  results[[arm]] <- pgc_run_arm(arm, count, seed, config)
}
stopifnot(identical(source_before, .crossform_source_tree_digest()))
summary <- do.call(rbind, lapply(results, pgc_summary, config = config))
record <- list(mode = mode, config = config, config_hash = config_hash, harness = harness,
  source_digest = source_before, runtime = sessionInfo(), results = results, summary = summary)
if (mode == "pilot") {
  record$counts <- vapply(results, function(result) {
    delta <- config$relative_equivalence_margin * result$design$reference_scale
    worst_sd <- max(vapply(seq_len(dim(result$values)[3]), function(j) sd(result$values[, "error", j]), 0))
    needed <- config$pilot_variance_headroom *
      (config$confidence_multiplier * worst_sd / (config$pilot_target_fraction * delta))^2
    max(config$minimum_replicates, ceiling(needed / config$count_rounding) * config$count_rounding)
  }, 0)
  saveRDS(record, pilot_path)
  print(record$counts)
} else {
  record$generator_moments <- do.call(rbind, lapply(results, function(result) {
    observed <- result$values[, "noise_moment", 1]
    error <- mean(observed) - result$design$noise_moment
    mcse <- sd(observed) / sqrt(length(observed))
    data.frame(arm = result$design$arm, exact = result$design$noise_moment, observed = mean(observed),
      error = error, mcse = mcse, pass = abs(error) <= config$confidence_multiplier * mcse + config$numerical_tolerance)
  }))
  null <- summary[summary$arm == "null", ]
  record$leakage_detected <- all(null$mean_leakage_error - config$confidence_multiplier * null$leakage_mcse >
    config$null_leakage_minimum * null$reference_scale)
  record$paired_comparisons <- do.call(rbind, lapply(results, function(result) {
    do.call(rbind, lapply(c("isotropic", "mismatched"), function(other) {
      differences <- result$values[, "gain", "model"] - result$values[, "gain", other]
      data.frame(arm = result$design$arm, comparison = paste("model", other, sep = " minus "),
        mean = mean(differences), mcse = sd(differences)/sqrt(length(differences)))
    }))
  }))
  record$pass <- all(summary$bias_pass & summary$precision_pass) && all(record$generator_moments$pass) && record$leakage_detected
  saveRDS(record, "benchmark-results/predictive-geometry-sampling.rds")
  write.csv(summary, "benchmark-results/predictive-geometry-sampling-summary.csv", row.names = FALSE)
  print(summary); print(record$generator_moments); print(record$paired_comparisons)
  stopifnot(record$pass)
}
