#!/usr/bin/env Rscript

# Isolated 50-60k-feature learned-crossnobis qualification.
#
# Run from the package root:
#
#   Rscript benchmarks/run-crossnobis-scale-gate.R . benchmark-results

if (!requireNamespace("parallel", quietly = TRUE) ||
    !requireNamespace("ps", quietly = TRUE)) {
  stop("Scale qualification requires the parallel and ps packages.")
}

arguments <- commandArgs(trailingOnly = TRUE)
repo <- if (length(arguments)) {
  normalizePath(arguments[[1L]], mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
output_dir <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(repo, "benchmark-results")
}
source(file.path(repo, "benchmarks", "provenance.R"), local = TRUE)
provenance <- crossform_benchmark_provenance(
  repo, "run-crossnobis-scale-gate.R"
)
full_dim <- if (length(arguments) >= 5L) {
  as.integer(arguments[3:5])
} else {
  c(52L, 42L, 25L)
}
active_slices <- if (length(arguments) >= 6L) {
  as.integer(arguments[[6L]])
} else {
  24L
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
result_path <- file.path(output_dir, "crossnobis-scale-gate.rds")
summary_path <- file.path(output_dir, "crossnobis-scale-gate-summary.csv")
ready_path <- file.path(output_dir, "crossnobis-scale-gate.ready")
unlink(c(result_path, ready_path))
source(file.path(repo, "benchmarks", "crossnobis-scale-worker.R"))

job <- parallel::mcparallel(
  run_crossnobis_scale_gate(
    repo, result_path, ready_path,
    full_dim = full_dim, active_slices = active_slices
  ),
  silent = TRUE
)
handle <- ps::ps_handle(job$pid)
peak_rss <- 0
measurement_started <- FALSE
collected <- NULL
repeat {
  alive <- tryCatch(ps::ps_is_running(handle), error = function(error) FALSE)
  if (!measurement_started && file.exists(ready_path)) {
    measurement_started <- TRUE
    peak_rss <- 0
  }
  if (measurement_started && alive) {
    rss <- tryCatch(
      unname(ps::ps_memory_info(handle)[["rss"]]),
      error = function(error) NA_real_
    )
    if (is.finite(rss)) peak_rss <- max(peak_rss, rss)
  }
  collected <- parallel::mccollect(job, wait = FALSE)
  if (!is.null(collected) || !alive) break
  Sys.sleep(0.01)
}
if (is.null(collected)) collected <- parallel::mccollect(job)
if (is.null(collected) || inherits(collected[[1L]], "try-error") ||
    !measurement_started || !file.exists(result_path)) {
  stop("The isolated crossnobis scale worker failed.")
}
result <- readRDS(result_path)
result$provenance <- provenance
peak_rss <- max(
  peak_rss,
  result$memory$baseline_rss_bytes,
  result$memory$rss_after_bytes,
  na.rm = TRUE
)
result$memory$os_peak_rss_bytes <- peak_rss
result$memory$incremental_peak_rss_bytes <- max(
  0, peak_rss - result$memory$baseline_rss_bytes
)
result$memory$rss_evidence <-
  "isolated_forked_child_polled_after_fixture_ready_signal"
result$gate <- list(
  maximum_incremental_rss_bytes = 4 * 1024^3,
  maximum_analysis_seconds = 30 * 60,
  memory_pass =
    result$memory$incremental_peak_rss_bytes <= 4 * 1024^3,
  runtime_pass = result$timing$analysis_seconds <= 30 * 60,
  mapping_pass = isTRUE(result$output$exact_full_index_mapping),
  structure_pass =
    !result$work$pair_atoms_materialized &&
    !result$work$pair_frame_materialized &&
    !result$work$factor_table_retained &&
    result$memory$persistent_factor_table_bytes == 0
)
result$gate$passed <- all(unlist(result$gate[c(
  "memory_pass", "runtime_pass", "mapping_pass", "structure_pass"
)]))
saveRDS(result, result_path, compress = "xz")

summary <- data.frame(
  fixture = result$fixture$id,
  features = result$fixture$features,
  support_min = result$support$min,
  support_median = result$support$median,
  support_mean = result$support$mean,
  support_max = result$support$max,
  union_pair_stored_nnz = result$support$union_pair_stored_nnz,
  structural_bytes = result$support$structural_bytes,
  local_metric_derivations = result$work$local_metric_derivations,
  factorization_units = result$work$factorization_units,
  residual_reads = result$reads$residual_total,
  evaluation_reads = result$reads$evaluation_total,
  planned_workspace_bytes = result$memory$planned_workspace_bytes,
  baseline_rss_bytes = result$memory$baseline_rss_bytes,
  os_peak_rss_bytes = result$memory$os_peak_rss_bytes,
  incremental_peak_rss_bytes = result$memory$incremental_peak_rss_bytes,
  analysis_seconds = result$timing$analysis_seconds,
  plan_seconds = result$timing$plan_and_residual_statistics_seconds,
  crossnobis_seconds = result$timing$crossnobis_seconds,
  output_mapping_seconds = result$timing$output_mapping_seconds,
  exact_full_index_mapping = result$output$exact_full_index_mapping,
  memory_pass = result$gate$memory_pass,
  runtime_pass = result$gate$runtime_pass,
  structure_pass = result$gate$structure_pass,
  passed = result$gate$passed,
  stringsAsFactors = FALSE
)
utils::write.csv(summary, summary_path, row.names = FALSE)
unlink(ready_path)
print(summary, row.names = FALSE, digits = 7)
if (!isTRUE(result$gate$passed)) {
  stop("The learned-crossnobis scale gate did not pass.")
}
