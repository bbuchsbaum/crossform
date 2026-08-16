#!/usr/bin/env Rscript

# Deterministic public-API benchmark for explicit scalar and multivariate
# measurement edges. Run with an installed crossform package. The optional
# first argument is the iteration count (default: 3).

if (!requireNamespace("crossform", quietly = TRUE)) {
  stop("Install crossform before running measurement benchmarks.")
}

library(crossform)

arguments <- commandArgs(trailingOnly = TRUE)
iterations <- if (length(arguments)) suppressWarnings(as.integer(arguments[[1L]])) else 3L
if (length(iterations) != 1L || is.na(iterations) || iterations < 1L) {
  stop("The optional iteration count must be one positive integer.")
}

set.seed(2026081232)
q <- 24L
p <- 96L
sample_space <- effect_space(
  paste0("sample", seq_len(q)), basis_id = "benchmark:samples:v1"
)
domain <- abstract_domain(p, id = "benchmark:neural:v1")
base <- matrix(rnorm(q * p), q, p)
relations <- relation(
  list(
    run1 = base,
    run2 = 0.85 * base + matrix(rnorm(q * p, sd = 0.35), q, p)
  ),
  effects = sample_space,
  domain = domain
)
center <- diag(q) - matrix(1 / q, q, q)
by <- variation_query(
  center / (q - 1), sample_space, "sample", "joint_covariance",
  provenance = list(estimator = "benchmark-centered-within-run")
)
over <- pairing(
  relations$partitions,
  relations$partitions,
  directed = TRUE,
  self_pairs = "allow_biased",
  independence = "not_independent"
)

scalar_ids <- paste0("scalar", sprintf("%02d", seq_len(32L)))
scalar_operators <- lapply(seq_along(scalar_ids), function(index) {
  operator <- matrix(0, 1L, p)
  support <- ((index - 1L) * 3L + seq_len(3L) - 1L) %% p + 1L
  operator[1L, support] <- c(0.5, 0.3, 0.2)
  operator
})
names(scalar_operators) <- scalar_ids
scalar_frame <- measurement_frame(
  scalar_operators, domain, id = "benchmark:scalar-frame:v1"
)
# A ring plus self-blocks: 64 explicit edges, not the 1,024 possible pairs.
scalar_from <- c(scalar_ids, scalar_ids)
scalar_to <- c(scalar_ids, c(tail(scalar_ids, -1L), scalar_ids[[1L]]))
scalar_between <- edge_frame(
  scalar_from, scalar_to, scalar_frame
)

multivariate_ids <- paste0("population", sprintf("%02d", seq_len(12L)))
multivariate_operators <- lapply(seq_along(multivariate_ids), function(index) {
  columns <- ((index - 1L) * 7L + seq_len(8L) - 1L) %% p + 1L
  operator <- matrix(0, 4L, p)
  operator[cbind(rep(seq_len(4L), each = 2L), columns)] <-
    rep(c(0.75, 0.25), 4L)
  operator
})
names(multivariate_operators) <- multivariate_ids
multivariate_frame <- measurement_frame(
  multivariate_operators, domain,
  id = "benchmark:multivariate-frame:v1"
)
# Twelve self-blocks and twelve directed ring edges out of 144 possible pairs.
multivariate_from <- c(multivariate_ids, multivariate_ids)
multivariate_to <- c(
  multivariate_ids,
  c(tail(multivariate_ids, -1L), multivariate_ids[[1L]])
)
multivariate_between <- edge_frame(
  multivariate_from, multivariate_to, multivariate_frame
)

benchmark_case <- function(label, between, routes) {
  records <- lapply(routes, function(route) {
    elapsed <- numeric(iterations)
    final <- NULL
    for (iteration in seq_len(iterations)) {
      timing <- system.time({
        final <- measurement_form(
          relations, between, by, over,
          compute = compute_policy(block_features = 32L),
          route = route
        )
      })
      elapsed[[iteration]] <- unname(timing[["elapsed"]])
    }
    blocks <- effect_coupling(final)$values
    data.frame(
      frame = label,
      route = route,
      iterations = iterations,
      median_elapsed_seconds = stats::median(elapsed),
      planned_workspace_bytes =
        final$receipt$execution$memory$planned_workspace_bytes,
      output_bytes = sum(final$block_index$length_elements) * 8,
      edge_count = nrow(final$block_index),
      possible_edge_count = length(between$from_frame$node_ids) *
        length(between$to_frame$node_ids),
      scientific_plan_id = final$plan$scientific_plan_id,
      checksum = sum(vapply(blocks, sum, numeric(1))),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, records)
  if (length(unique(result$scientific_plan_id)) != 1L ||
      max(abs(result$checksum - result$checksum[[1L]])) > 1e-10) {
    stop("Execution routes changed scientific identity or numerical values.")
  }
  result
}

scalar <- benchmark_case(
  "scalar",
  scalar_between,
  c("auto", "pull_h", "forward_k", "multivariate_blocks",
    "scalar_stack", "factorized_h")
)
multivariate <- benchmark_case(
  "requested_multivariate",
  multivariate_between,
  c("auto", "pull_h", "forward_k", "multivariate_blocks", "factorized_h")
)
record <- rbind(scalar, multivariate)
stopifnot(
  all(record$edge_count < record$possible_edge_count),
  all(is.finite(record$median_elapsed_seconds)),
  all(record$median_elapsed_seconds >= 0),
  all(record$planned_workspace_bytes > record$output_bytes)
)
print(record, row.names = FALSE, digits = 6)
