args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: memory-worker.R <repo-root> <scenario-id> <output-rds> <ready-file>")
}
repo <- normalizePath(args[[1L]], mustWork = TRUE)
scenario_id <- args[[2L]]
output <- args[[3L]]
ready_path <- args[[4L]]

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The benchmark worker requires the suggested package 'pkgload'.")
}
pkgload::load_all(repo, quiet = TRUE)
scenarios <- effectagram:::.memory_benchmark_scenarios()
selected <- scenarios[scenarios$id == scenario_id, , drop = FALSE]
if (nrow(selected) != 1L) stop("unknown scenario: ", scenario_id)
saveRDS(effectagram:::.run_memory_benchmark_case(selected,
  ready_path = ready_path), output)
