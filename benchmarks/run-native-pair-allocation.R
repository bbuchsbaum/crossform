#!/usr/bin/env Rscript

# Cumulative-allocation court for the fused structured pair-query kernel.
#
# The query-first scale gate records peak R heap for a public end-to-end map.
# This narrower court answers a different legacy acceptance question: does the
# fused native kernel remove at least 30 percent of the cumulative R allocation
# incurred by the retained two-pass R oracle on the same work? Each timed call
# receives a fresh Rprofmem log after both routes have been warmed.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("The native pair-allocation court requires devtools.")
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
repetitions <- if (length(arguments) >= 3L) {
  as.integer(arguments[[3L]])
} else {
  5L
}
if (is.na(repetitions) || repetitions < 3L) {
  stop("The native pair-allocation court requires at least three repetitions.")
}

source(file.path(repo, "benchmarks", "provenance.R"), local = TRUE)
devtools::load_all(repo, quiet = TRUE)
cat("Loaded crossform\n"); flush.console()
provenance <- crossform_benchmark_provenance(
  repo, "run-native-pair-allocation.R"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260821L)
q <- 48L
p <- 512L
partition_count <- 4L
blocks <- lapply(seq_len(partition_count), function(partition) {
  matrix(stats::rnorm(q * p, mean = partition / 20), q, p)
})
pairs <- t(utils::combn(seq_len(q), 2L))
edges <- expand.grid(
  left = seq_len(partition_count),
  right = seq_len(partition_count),
  KEEP.OUT.ATTRS = FALSE
)
edges <- edges[edges$left != edges$right, , drop = FALSE]
edge_weight <- (-1)^seq_len(nrow(edges)) *
  seq_len(nrow(edges)) / sum(seq_len(nrow(edges)))

native_call <- function() {
  crossform:::.fused_pair_difference_atoms(
    blocks, blocks, edges$left, edges$right, edge_weight,
    pairs[, 1L], pairs[, 2L]
  )
}
oracle_call <- function() {
  crossform:::.pair_difference_accumulate_oracle(
    blocks, blocks, edges$left, edges$right, edge_weight,
    pairs[, 1L], pairs[, 2L]
  )
}

# Warm compilation, dispatch, and pages before either route enters the court.
invisible(native_call())
invisible(oracle_call())
invisible(gc())

native_reference <- native_call()
oracle_reference <- oracle_call()
oracle_error <- max(abs(native_reference - oracle_reference))
output_bytes <- as.double(utils::object.size(native_reference))
checksum <- sum(native_reference)
rm(native_reference, oracle_reference)
invisible(gc())

measure_once <- function(route, repetition) {
  allocation_log <- tempfile(fileext = ".Rprofmem")
  on.exit(unlink(allocation_log), add = TRUE)
  fun <- if (identical(route, "native_fused")) native_call else oracle_call
  invisible(gc())
  utils::Rprofmem(allocation_log)
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(fun(), finally = utils::Rprofmem(NULL))
  elapsed <- proc.time()[["elapsed"]] - started
  allocation <- crossform:::.allocation_summary(allocation_log)
  data.frame(
    route = route,
    repetition = as.integer(repetition),
    elapsed_seconds = elapsed,
    allocation_count = allocation$allocation_count,
    allocated_bytes = allocation$allocated_bytes,
    largest_allocation_bytes = allocation$largest_allocation_bytes,
    output_bytes = as.double(utils::object.size(value)),
    checksum = sum(value),
    stringsAsFactors = FALSE
  )
}

records <- vector("list", 2L * repetitions)
position <- 0L
for (repetition in seq_len(repetitions)) {
  order <- if (repetition %% 2L) {
    c("native_fused", "r_two_pass_oracle")
  } else {
    c("r_two_pass_oracle", "native_fused")
  }
  for (route in order) {
    position <- position + 1L
    records[[position]] <- measure_once(route, repetition)
  }
}
records <- do.call(rbind, records)

median_for <- function(column, route) {
  stats::median(records[[column]][records$route == route])
}
native_allocated <- median_for("allocated_bytes", "native_fused")
oracle_allocated <- median_for("allocated_bytes", "r_two_pass_oracle")
native_seconds <- median_for("elapsed_seconds", "native_fused")
oracle_seconds <- median_for("elapsed_seconds", "r_two_pass_oracle")
allocation_ratio <- native_allocated / oracle_allocated
runtime_ratio <- native_seconds / oracle_seconds

artifact <- list(
  schema_version = 1L,
  provenance = provenance,
  measurement = paste0(
    "Rprofmem cumulative bytes per warmed call; one fresh log per route and ",
    "repetition"
  ),
  dimensions = list(
    conditions = q,
    features = p,
    partitions = partition_count,
    ordered_edges = nrow(edges),
    pairs = nrow(pairs),
    repetitions = repetitions,
    output_bytes = output_bytes
  ),
  numerical = list(
    tolerance = 1e-12,
    max_abs_error = oracle_error,
    checksum = checksum
  ),
  records = records,
  summary = list(
    native_median_allocated_bytes = native_allocated,
    oracle_median_allocated_bytes = oracle_allocated,
    native_to_oracle_allocation_ratio = allocation_ratio,
    native_median_elapsed_seconds = native_seconds,
    oracle_median_elapsed_seconds = oracle_seconds,
    native_to_oracle_runtime_ratio = runtime_ratio
  ),
  gate = list(
    maximum_allocation_ratio = 0.70,
    numerical_pass = oracle_error <= 1e-12,
    allocation_pass = is.finite(allocation_ratio) && allocation_ratio <= 0.70
  )
)
artifact$gate$passed <- artifact$gate$numerical_pass &&
  artifact$gate$allocation_pass

summary <- data.frame(
  conditions = q,
  features = p,
  partitions = partition_count,
  ordered_edges = nrow(edges),
  pairs = nrow(pairs),
  repetitions = repetitions,
  output_bytes = output_bytes,
  oracle_max_abs_error = oracle_error,
  native_median_allocated_bytes = native_allocated,
  oracle_median_allocated_bytes = oracle_allocated,
  native_to_oracle_allocation_ratio = allocation_ratio,
  native_median_elapsed_seconds = native_seconds,
  oracle_median_elapsed_seconds = oracle_seconds,
  native_to_oracle_runtime_ratio = runtime_ratio,
  passed = artifact$gate$passed,
  stringsAsFactors = FALSE
)

saveRDS(
  artifact,
  file.path(output_dir, "native-pair-allocation.rds"),
  version = 2L
)
utils::write.csv(
  summary,
  file.path(output_dir, "native-pair-allocation-summary.csv"),
  row.names = FALSE
)
print(summary, row.names = FALSE, digits = 6)
if (!artifact$gate$passed) {
  stop("The fused pair-query kernel did not pass its allocation court.")
}
