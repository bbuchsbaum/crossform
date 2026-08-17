need_namespace <- function(package) {
  tryCatch(
    {
      loadNamespace(package)
      invisible(TRUE)
    },
    error = function(error) stop(package, " is required: ", conditionMessage(error))
  )
}
need_namespace("pkgload")
need_namespace("processx")
need_namespace("ps")

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path(repo, "benchmark-results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pkgload::load_all(repo, quiet = TRUE)
source(file.path(repo, "benchmarks", "provenance.R"), local = TRUE)
provenance <- crossform_benchmark_provenance(repo, "run-memory-benchmarks.R")
scenarios <- crossform:::.memory_benchmark_scenarios()
summary_rows <- vector("list", nrow(scenarios))

for (scenario_index in seq_len(nrow(scenarios))) {
  id <- scenarios$id[[scenario_index]]
  result_path <- file.path(output_dir, paste0(id, ".rds"))
  ready_path <- file.path(output_dir, paste0(id, ".ready"))
  unlink(c(result_path, ready_path))
  process <- processx::process$new(
    file.path(R.home("bin"), "Rscript"),
    c(file.path(repo, "benchmarks", "memory-worker.R"), repo, id, result_path,
      ready_path),
    stdout = file.path(output_dir, paste0(id, ".stdout.txt")),
    stderr = file.path(output_dir, paste0(id, ".stderr.txt")),
    cleanup = TRUE
  )
  peak_rss <- 0
  measurement_started <- FALSE
  while (process$is_alive()) {
    if (!measurement_started && file.exists(ready_path)) {
      measurement_started <- TRUE
      peak_rss <- 0
    }
    if (measurement_started) {
      memory <- tryCatch(process$get_memory_info(), error = function(error) NULL)
      if (!is.null(memory) && is.finite(memory[["rss"]])) {
        peak_rss <- max(peak_rss, unname(memory[["rss"]]))
      }
    }
    Sys.sleep(0.005)
  }
  process$wait()
  if (process$get_exit_status() != 0L || !file.exists(result_path) ||
      !measurement_started) {
    stop("benchmark worker failed: ", id)
  }
  result <- readRDS(result_path)
  result$provenance <- provenance
  peak_rss <- max(peak_rss, result$rss_before_bytes, result$rss_after_bytes,
    na.rm = TRUE)
  result$os_peak_rss_bytes <- peak_rss
  result$aggregate_peak_rss_bytes <- result$os_peak_rss_bytes + result$worker_peak_rss_bytes
  result$incremental_peak_rss_bytes <- max(0, peak_rss - result$rss_before_bytes)
  result$rss_evidence <- "isolated_child_process_polled_os_peak"
  saveRDS(result, result_path)
  summary_rows[[scenario_index]] <- data.frame(
    id = id,
    storage = result$scenario$storage,
    phase = result$scenario$phase,
    os_peak_rss_bytes = result$os_peak_rss_bytes,
    incremental_peak_rss_bytes = result$incremental_peak_rss_bytes,
    baseline_rss_bytes = result$rss_before_bytes,
    absolute_peak_rss_bytes = result$aggregate_peak_rss_bytes,
    planned_workspace_bytes = result$plan$planned_workspace_bytes,
    allocated_bytes = result$allocation$allocated_bytes,
    largest_allocation_bytes = result$allocation$largest_allocation_bytes,
    frame_bytes = result$plan$categories[["frame"]],
    resident_source_bytes = result$plan$categories[["resident_source"]],
    max_measured_live_temporary_bytes = result$total$max_live_temporary_bytes,
    stringsAsFactors = FALSE
  )
  unlink(ready_path)
}
write.csv(do.call(rbind, summary_rows),
  file.path(output_dir, "memory-benchmark-summary.csv"),
  row.names = FALSE)
