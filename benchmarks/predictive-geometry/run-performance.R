# Independent child processes, one measured workload at a time.
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
output_dir <- if (length(args) >= 2L) args[2] else file.path(repo, "benchmark-results")
selected <- if (length(args) >= 3L) args[3] else "all"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(repo, "benchmarks/predictive-geometry/performance-config.R"))
source(file.path(repo, "benchmarks/provenance.R"))
config <- pg_performance_config()
provenance <- crossform_benchmark_provenance(repo, "predictive-geometry/run-performance.R")
harness_files <- paste0("benchmarks/predictive-geometry/", c("performance-config.R", "performance.R",
  "performance-worker.R", "run-performance.R"))
harness <- stats::setNames(vapply(harness_files, function(file)
  digest::digest(file = file.path(repo, file), algo = "sha256"), ""), harness_files)
system <- list(info = Sys.info(), cores = ps::ps_cpu_count(), memory = ps::ps_system_memory(),
  load_before = ps::ps_loadavg(), cpu_before = ps::ps_system_cpu_times(),
  load_scope = "Observed host load includes other applications; no idle-machine claim is made.")
cases <- config$cases[if (selected == "all") rep(TRUE, nrow(config$cases)) else config$cases$id == selected, ]
if (!nrow(cases)) stop("Select all, tiny, medium or large")
records <- list(); summaries <- list()
for (i in seq_len(nrow(cases))) {
  id <- cases$id[i]
  routes <- if (id == "large") "block" else c("memory", "block", "dense")
  for (route in routes) {
    key <- paste(id, route, sep = "-")
    prefix <- file.path(output_dir, paste0("predictive-performance-", key))
    result_path <- paste0(prefix, ".rds")
    signals <- paste0(prefix, c(".ready", ".start", ".done"))
    unlink(c(result_path, signals))
    process <- processx::process$new(file.path(R.home("bin"), "Rscript"),
      c(file.path(repo, "benchmarks/predictive-geometry/performance-worker.R"), repo,
        id, route, result_path, signals), stdout = paste0(prefix, ".stdout.txt"),
      stderr = paste0(prefix, ".stderr.txt"), cleanup = TRUE,
      env = c("current", LC_ALL = "en_US.UTF-8"))
    observing <- FALSE; peak <- 0
    while (process$is_alive()) {
      if (!observing && file.exists(signals[1])) {
        initial <- readRDS(signals[1]); peak <- initial$baseline_rss
        observing <- TRUE; writeLines("start", signals[2])
      }
      if (observing && !file.exists(signals[3])) {
        rss <- tryCatch(process$get_memory_info()[["rss"]], error = function(e) NA_real_)
        if (is.finite(rss)) peak <- max(peak, rss)
      }
      Sys.sleep(.01)
    }
    process$wait()
    if (process$get_exit_status() != 0L || !file.exists(result_path) || !observing)
      stop("Performance worker failed: ", key, "; inspect ", prefix, ".stderr.txt")
    result <- readRDS(result_path)
    stopifnot(identical(result$source_digest, provenance$source_digest))
    result$peak_rss_bytes <- max(peak, result$final_rss_bytes)
    result$incremental_rss_bytes <- max(0, result$peak_rss_bytes - result$baseline_rss_bytes)
    result$setup_seconds <- initial$setup_seconds
    result$provenance <- provenance; result$harness_digests <- harness; result$system <- system
    result$system$load_after <- ps::ps_loadavg(); result$system$cpu_after <- ps::ps_system_cpu_times()
    saveRDS(result, result_path)
    records[[key]] <- result
    summaries[[key]] <- data.frame(case = id, route = route, seconds = result$elapsed_seconds,
      planned_bytes = result$admission$planned_workspace_bytes,
      incremental_rss_bytes = result$incremental_rss_bytes,
      largest_R_allocation = result$allocation$largest_R_allocation,
      retained_bytes = result$value$stored_bytes,
      read_calls = sum(result$source_reads$calls), read_columns = sum(result$source_reads$columns),
      maximum_read_width = result$source_reads$maximum_width,
      packed_width = if (route == "dense") cases$conditions[i] * (cases$conditions[i] + 1) / 2 else
        cases$model_rank[i] * (cases$model_rank[i] + 1) / 2,
      kernels = result$calls$kernels, eigen = result$calls$eigen,
      form_error = result$reference_errors[["form"]], gain_error = result$reference_errors[["gain"]])
    unlink(signals)
    print(summaries[[key]])
  }
  if (length(routes) > 1L) {
    reference <- records[[paste(id, "dense", sep = "-")]]
    for (route in c("memory", "block")) {
      actual <- records[[paste(id, route, sep = "-")]]
      for (field in c("forms", "spectra", "modes")) {
        err <- max(abs(actual$value[[field]] - reference$value[[field]])) /
          max(1, max(abs(reference$value[[field]])))
        if (err > config$relative_tolerance) stop("Matched output disagrees: ", id, " ", route, " ", field)
      }
      for (field in c("table", "diagnostics")) {
        err <- max(abs(as.matrix(actual$value[[field]][, -1]) - as.matrix(reference$value[[field]][, -1]))) /
          max(1, max(abs(as.matrix(reference$value[[field]][, -1]))))
        if (err > config$relative_tolerance) stop("Matched summary disagrees: ", id, " ", route, " ", field)
      }
      stopifnot(identical(actual$source_reads, reference$source_reads))
    }
  }
}
stopifnot(identical(provenance$source_digest, .crossform_source_tree_digest(repo)),
  identical(harness, stats::setNames(vapply(harness_files, function(file)
    digest::digest(file = file.path(repo, file), algo = "sha256"), ""), harness_files)))
summary <- do.call(rbind, summaries); rownames(summary) <- NULL
write.csv(summary, file.path(output_dir, paste0("predictive-performance-", selected, "-summary.csv")), row.names = FALSE)
saveRDS(list(summary = summary, provenance = provenance, harness_digests = harness,
  config = config, system = system, case_artifacts = paste0("predictive-performance-", names(records), ".rds"),
  parity = TRUE), file.path(output_dir, paste0("predictive-performance-", selected, ".rds")))
