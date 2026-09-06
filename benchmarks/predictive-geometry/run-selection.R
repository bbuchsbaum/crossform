# Rscript benchmarks/predictive-geometry/run-selection.R pilot|run
mode <- commandArgs(TRUE); stopifnot(length(mode) == 1L, mode %in% c("pilot", "run"))
pkgload::load_all(".", quiet = TRUE)
source("benchmarks/provenance.R")
source("benchmarks/predictive-geometry/calibration.R")
source("benchmarks/predictive-geometry/selection-calibration.R")
source("benchmarks/predictive-geometry/selection-config.R")
config <- predictive_selection_config
config_hash <- digest::digest(config, algo = "sha256")
harness <- vapply(c("benchmarks/predictive-geometry/calibration.R",
  "benchmarks/predictive-geometry/selection-calibration.R", "benchmarks/predictive-geometry/selection-config.R",
  "benchmarks/predictive-geometry/run-selection.R"), function(path) digest::digest(file = path, algo = "sha256"), "")
source_before <- .crossform_source_tree_digest()
pilot_path <- "benchmark-results/predictive-geometry-selection-pilot.rds"
if (mode == "run") {
  pilot <- readRDS(pilot_path)
  stopifnot(identical(pilot$config_hash, config_hash), identical(pilot$harness, harness))
}
results <- list()
for (i in seq_along(config$arms)) {
  arm <- config$arms[[i]]
  count <- if (mode == "pilot") config$pilot_replicates else pilot$counts[[arm]]
  seed <- (if (mode == "pilot") config$pilot_seed else config$production_seed) + 100L * i
  results[[arm]] <- pgs_run_arm(arm, count, seed, config)
}
stopifnot(identical(source_before, .crossform_source_tree_digest()))
summary <- do.call(rbind, lapply(results, pgs_summary, config = config))
record <- list(mode = mode, config = config, config_hash = config_hash, harness = harness,
  source_digest = source_before, runtime = sessionInfo(), results = results, summary = summary)
if (mode == "pilot") {
  record$counts <- vapply(results, function(result) {
    delta <- config$relative_equivalence_margin * result$design$reference_scale
    needed <- config$pilot_variance_headroom *
      (config$confidence_multiplier * sd(result$values[, "error"]) / (config$pilot_target_fraction * delta))^2
    max(config$minimum_replicates, ceiling(needed/config$count_rounding) * config$count_rounding)
  }, 0)
  saveRDS(record, pilot_path); print(record$counts)
} else {
  null <- summary[summary$arm == "null", ]
  record$leakage_detected <- null$leakage_bias - config$confidence_multiplier * null$leakage_mcse >
    config$leakage_minimum * null$reference_scale
  record$pass <- all(summary$bias_pass & summary$precision_pass) && record$leakage_detected
  saveRDS(record, "benchmark-results/predictive-geometry-selection.rds")
  write.csv(summary, "benchmark-results/predictive-geometry-selection-summary.csv", row.names = FALSE)
  print(summary); stopifnot(record$pass)
}
