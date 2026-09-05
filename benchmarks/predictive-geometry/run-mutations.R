args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
output <- if (length(args) >= 2L) args[2] else file.path(repo, "benchmark-results/predictive-geometry-mutations")
selected <- if (length(args) >= 3L) args[3] else "all"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source(file.path(repo, "benchmarks/predictive-geometry/mutations.R"))
source(file.path(repo, "benchmarks/provenance.R"))
catalog <- pg_mutation_catalog()
if (selected != "all") catalog <- Filter(function(x) x$id == selected, catalog)
if (!length(catalog)) stop("Unknown mutation selection")
provenance <- crossform_benchmark_provenance(repo, "predictive-geometry/run-mutations.R")
harness_files <- file.path("benchmarks/predictive-geometry", c("mutations.R", "mutation-worker.R", "run-mutations.R"))
harness <- stats::setNames(vapply(harness_files, function(file)
  digest::digest(file = file.path(repo, file), algo = "sha256"), ""), harness_files)
rows <- vector("list", length(catalog)); invalid <- character()
for (i in seq_along(catalog)) {
  mutation <- catalog[[i]]
  prefix <- file.path(output, mutation$id)
  result_file <- paste0(prefix, ".rds")
  unlink(result_file)
  process <- processx::run(file.path(R.home("bin"), "Rscript"),
    c(file.path(repo, "benchmarks/predictive-geometry/mutation-worker.R"), repo, mutation$id, result_file),
    stdout = paste0(prefix, ".stdout.txt"), stderr = paste0(prefix, ".stderr.txt"),
    error_on_status = FALSE, timeout = 180000, env = c("current", LC_ALL = "en_US.UTF-8"))
  if (!file.exists(result_file)) {
    invalid <- c(invalid, mutation$id)
    rows[[i]] <- data.frame(id = mutation$id, primary = mutation$primary, status = "INVALID",
      primary_failures = 0L, production_calls = 0L)
  } else {
    result <- readRDS(result_file)
    stopifnot(identical(result$source_digest, provenance$source_digest))
    result$provenance <- provenance; result$harness_digests <- harness
    saveRDS(result, result_file)
    rows[[i]] <- data.frame(id = result$id, primary = result$primary,
      status = if (result$killed) "KILLED" else "SURVIVED", primary_failures = nrow(result$primary_failures),
      production_calls = sum(result$mutated$hits))
  }
  print(rows[[i]])
}
stopifnot(identical(provenance$source_digest, .crossform_source_tree_digest(repo)),
  identical(harness, stats::setNames(vapply(harness_files, function(file)
    digest::digest(file = file.path(repo, file), algo = "sha256"), ""), harness_files)))
matrix <- do.call(rbind, rows); rownames(matrix) <- NULL
write.csv(matrix, file.path(output, paste0("mutation-matrix-", selected, ".csv")), row.names = FALSE)
saveRDS(list(matrix = matrix, provenance = provenance, harness_digests = harness,
  all_critical_killed = all(matrix$status == "KILLED"), scope = selected),
  file.path(output, paste0("mutation-matrix-", selected, ".rds")))
if (any(matrix$status != "KILLED")) stop("Mutation gate failed; invalid definitions and survivors are not evidence of a kill")
