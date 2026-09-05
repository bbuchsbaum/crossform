args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop("usage: performance-worker.R <repo> <case> <route> <result> <ready> <start> <done>")
repo <- normalizePath(args[1], mustWork = TRUE); case_id <- args[2]; route <- args[3]
output <- args[4]; ready <- args[5]; start <- args[6]; done <- args[7]
pkgload::load_all(repo, quiet = TRUE)
source(file.path(repo, "benchmarks/predictive-geometry/performance-config.R"))
source(file.path(repo, "benchmarks/predictive-geometry/performance.R"))
source(file.path(repo, "benchmarks/provenance.R"))
config <- pg_performance_config()
case <- config$cases[config$cases$id == case_id, ]
if (nrow(case) != 1L || !route %in% c("memory", "block", "dense")) stop("Unknown performance case or route")
source_before <- .crossform_source_tree_digest(repo)
setup_started <- proc.time()[["elapsed"]]
f <- pg_perf_fixture(case, config)
admission <- pg_perf_admission(f, route)
stopifnot(sum(f$reads$calls) == 0L, admission$planned_workspace_bytes <= config$workspace_bytes)
instrument <- pg_perf_instrument(if (route != "dense") f$basis$dimension else NULL)
directory <- tempfile("prediction-records-"); dir.create(directory)
allocation_log <- tempfile("prediction-allocations-")
gc()
baseline <- unname(ps::ps_memory_info()[["rss"]])
saveRDS(list(baseline_rss = baseline, setup_seconds = proc.time()[["elapsed"]] - setup_started), ready)
deadline <- Sys.time() + 30
while (!file.exists(start)) {
  if (Sys.time() > deadline) stop("Performance monitor did not acknowledge readiness")
  Sys.sleep(.005)
}
utils::Rprofmem(allocation_log)
started <- proc.time()[["elapsed"]]
value <- if (route == "dense") pg_perf_full_space(f) else pg_perf_public(f, route, directory)
elapsed <- proc.time()[["elapsed"]] - started
utils::Rprofmem(NULL)
final_rss <- unname(ps::ps_memory_info()[["rss"]])
saveRDS(list(final_rss = final_rss), done)
counters <- as.list(instrument$counters)
instrument$cleanup()
reads <- as.list(f$reads)
allocation <- pg_perf_allocation(allocation_log)

# Reference checking follows the monitored region and does not count as
# source reads or numerical work in either measured route.
reference <- pg_perf_reference(f, value$rows)
errors <- c(form = pg_perf_error(value$forms, reference$forms),
  gain = pg_perf_error(value$table$gain[value$rows], reference$gain))
stopifnot(all(is.finite(errors)), all(errors <= config$relative_tolerance),
  identical(source_before, .crossform_source_tree_digest(repo)))
if (route != "dense") stopifnot(counters$materializations == 2L,
  all(counters$packed_widths == f$basis$dimension * (f$basis$dimension + 1L) / 2L),
  counters$paths == case$nodes)
legacy_overlap <- if (case_id == "tiny" && route == "memory") pg_perf_legacy_overlap(f) else NULL
reuse <- if (case_id %in% c("tiny", "medium") && route == "memory") pg_perf_reuse_comparison(f) else NULL
result <- list(case = case, route = route, config = config, source_digest = source_before,
  admission = admission, elapsed_seconds = elapsed, baseline_rss_bytes = baseline,
  final_rss_bytes = final_rss, source_reads = reads, calls = counters,
  allocation = allocation, reference_errors = errors,
  reference_rows = value$rows, value = value, legacy_overlap = legacy_overlap, reuse = reuse,
  baseline_scope = if (route == "dense")
    "Original-space complete executor then independent dense fit/readout; equivalent numerical outputs, without predictive record validation or durability"
    else "Public fit and score including record validation and requested persistence")
saveRDS(result, output)
unlink(c(directory, allocation_log), recursive = TRUE)
cat(case_id, route, "passed; seconds", elapsed, "reference errors", errors, "\n")
