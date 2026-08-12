args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
output_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path(repo, "benchmark-results")
library_path <- if (length(args) >= 3L) normalizePath(args[[3L]], mustWork = TRUE) else
  .libPaths()[[1L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("processx", quietly = TRUE) ||
    !requireNamespace("ps", quietly = TRUE)) {
  stop("processx and ps are required")
}

tree_rss <- function(process) {
  handle <- tryCatch(ps::ps_handle(process$get_pid()), error = function(error) NULL)
  if (is.null(handle)) return(NA_real_)
  descendants <- tryCatch(ps::ps_children(handle, recursive = TRUE),
    error = function(error) list())
  handles <- c(list(handle), descendants)
  values <- vapply(handles, function(value) {
    tryCatch(unname(ps::ps_memory_info(value)[["rss"]]),
      error = function(error) NA_real_)
  }, numeric(1))
  sum(values[is.finite(values)])
}

modes <- c(
  "sequential-response",
  "shard-response-cold", "shard-response-warm",
  "shard-relation-cold", "shard-relation-warm"
)
results <- vector("list", length(modes))
names(results) <- modes

for (mode in modes) {
  result_path <- tempfile(paste0(mode, "-"), fileext = ".rds")
  ready_path <- tempfile(paste0(mode, "-"), fileext = ".ready")
  process <- processx::process$new(
    file.path(R.home("bin"), "Rscript"),
    c(
      file.path(repo, "benchmarks", "shard-admission-worker.R"),
      repo, mode, result_path, ready_path, library_path
    ),
    stdout = file.path(output_dir, paste0(mode, ".stdout.txt")),
    stderr = file.path(output_dir, paste0(mode, ".stderr.txt")),
    cleanup = TRUE
  )
  peak_tree_rss <- 0
  peak_after_ready <- 0
  while (process$is_alive()) {
    rss <- tree_rss(process)
    if (is.finite(rss)) peak_tree_rss <- max(peak_tree_rss, rss)
    if (file.exists(ready_path) && is.finite(rss)) {
      peak_after_ready <- max(peak_after_ready, rss)
    }
    Sys.sleep(0.01)
  }
  process$wait()
  if (process$get_exit_status() != 0L || !file.exists(result_path)) {
    stop("shard admission worker failed: ", mode)
  }
  value <- readRDS(result_path)
  value$aggregate_peak_tree_rss_bytes <- peak_tree_rss
  value$aggregate_peak_after_ready_bytes <- peak_after_ready
  results[[mode]] <- value
  unlink(c(result_path, ready_path))
}

# Cleanup stress is a separate isolated child because the scientific benchmark
# must not inherit a pool or shared segment from a failed dispatch.
stress_script <- tempfile("effectagram-shard-stress-", fileext = ".R")
stress_result <- tempfile("effectagram-shard-stress-", fileext = ".rds")
writeLines(c(
  sprintf(".libPaths(c(%s, .libPaths()))", dQuote(library_path)),
  "library(shard)",
  "x <- share(matrix(rnorm(2e6), 1000), backing = 'mmap', readonly = TRUE)",
  "path <- shared_info(x)$path",
  "pool <- pool_create(2, rss_limit = '2GB', rss_drift_threshold = Inf)",
  "failed <- try(shard_map(shards(8, block_size = 1, workers = 2),",
  "  function(shard, x) { if (shard$id == 3L) stop('injected failure'); sum(x[, 1]) },",
  "  borrow = list(x = x), workers = 2, max_retries = 0, autotune = FALSE),",
  "  silent = TRUE)",
  "caught <- inherits(failed, 'try-error') || !succeeded(failed)",
  "pool_stop(pool)",
  "pool_stopped <- is.null(pool_get())",
  "close(x)",
  "backing_removed <- !file.exists(path)",
  sprintf("saveRDS(list(caught = caught, pool_stopped = pool_stopped, backing_removed = backing_removed), %s)",
    dQuote(stress_result))
), stress_script)
stress <- processx::run(file.path(R.home("bin"), "Rscript"), stress_script,
  error_on_status = FALSE, timeout = 600000)
if (stress$status != 0L || !file.exists(stress_result)) {
  stop("shard cleanup stress worker failed")
}
cleanup_stress <- readRDS(stress_result)
unlink(c(stress_script, stress_result))

reference <- results[["sequential-response"]]$geometry
for (mode in modes) {
  results[[mode]]$parity <- list(
    total = isTRUE(all.equal(results[[mode]]$geometry$total,
      reference$total, tolerance = 1e-11, check.attributes = FALSE)),
    coherent = isTRUE(all.equal(results[[mode]]$geometry$coherent,
      reference$coherent, tolerance = 1e-11, check.attributes = FALSE))
  )
  results[[mode]]$geometry <- NULL
}

baseline <- results[["sequential-response"]]$elapsed_seconds
summary <- do.call(rbind, lapply(modes, function(mode) {
  value <- results[[mode]]
  data.frame(
    mode = mode,
    phase = value$phase,
    staging = value$staging,
    parallel = value$parallel,
    elapsed_seconds = value$elapsed_seconds,
    speedup_vs_sequential = baseline / value$elapsed_seconds,
    staging_seconds = value$staging_seconds,
    pool_seconds = value$pool_seconds,
    aggregate_peak_tree_rss_bytes = value$aggregate_peak_tree_rss_bytes,
    aggregate_peak_after_ready_bytes = value$aggregate_peak_after_ready_bytes,
    shared_bytes = value$shared_bytes,
    task_payload_bytes = value$task_payload_bytes,
    total_parity = value$parity$total,
    coherent_parity = value$parity$coherent,
    pool_stopped = if (is.na(value$cleanup$pool_stopped)) TRUE else
      value$cleanup$pool_stopped,
    backing_removed = if (is.na(value$cleanup$backing_removed)) TRUE else
      value$cleanup$backing_removed,
    stringsAsFactors = FALSE
  )
}))

admit <- all(summary$total_parity & summary$coherent_parity) &&
  all(summary$pool_stopped & summary$backing_removed) &&
  all(unlist(cleanup_stress, use.names = FALSE)) &&
  any(summary$speedup_vs_sequential[summary$parallel & summary$phase == "cold"] > 1.1)

artifact <- list(
  schema_version = 1L,
  benchmark_date = as.character(Sys.Date()),
  shard_version = as.character(utils::packageVersion("shard", lib.loc = library_path)),
  effectagram_version = as.character(utils::packageVersion(
    "effectagram", lib.loc = library_path)),
  gate = list(
    numerical_parity = all(summary$total_parity & summary$coherent_parity),
    cleanup = all(summary$pool_stopped & summary$backing_removed),
    cleanup_stress = all(unlist(cleanup_stress, use.names = FALSE)),
    cold_speedup_threshold = 1.1,
    cold_speedup = any(summary$speedup_vs_sequential[
      summary$parallel & summary$phase == "cold"] > 1.1),
    admitted = admit
  ),
  cleanup_stress = cleanup_stress,
  summary = summary,
  results = results
)
saveRDS(artifact, file.path(output_dir, "shard-admission.rds"), version = 3)
write.csv(summary, file.path(output_dir, "shard-admission-summary.csv"),
  row.names = FALSE)
print(summary)
print(artifact$gate)
