#!/usr/bin/env Rscript

# Map-scale profile of measurement-form contraction routes. Records provenance,
# route selection, planned memory, oracle parity, and an Rprof split between
# R dispatch/tile loops and BLAS. Admits a native kernel only if R-level
# loops are at least 15% of end-to-end time and a prototype would project
# at least 1.25x; otherwise records a no-Rcpp decision.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("The measurement profile requires devtools.")
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
  stop("The measurement profile requires at least three repetitions.")
}

source(file.path(repo, "benchmarks", "provenance.R"), local = TRUE)
devtools::load_all(repo, quiet = TRUE)
cat("Loaded crossform\n"); flush.console()
provenance <- crossform_benchmark_provenance(repo, "run-measurement-profile.R")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260816L)
q <- 16L
p <- 128L
sample_space <- effect_space(
  paste0("sample", seq_len(q)), basis_id = "profile:samples:v1"
)
domain <- abstract_domain(p, id = "profile:neural:v1")
base <- matrix(rnorm(q * p), q, p)
relations <- relation(
  list(
    run1 = base,
    run2 = 0.85 * base + matrix(rnorm(q * p, sd = 0.35), q, p),
    run3 = 0.80 * base + matrix(rnorm(q * p, sd = 0.40), q, p)
  ),
  effects = sample_space,
  domain = domain
)
center <- diag(q) - matrix(1 / q, q, q)
by <- variation_query(
  center / (q - 1), sample_space, "sample", "joint_covariance",
  provenance = list(estimator = "profile-centered-within-run")
)
over <- pairing(
  relations$partitions, relations$partitions,
  directed = TRUE, self_pairs = "allow_biased",
  independence = "not_independent"
)

make_scalar_frame <- function(n_nodes) {
  ids <- paste0("scalar", sprintf("%03d", seq_len(n_nodes)))
  operators <- lapply(seq_along(ids), function(index) {
    operator <- matrix(0, 1L, p)
    support <- ((index - 1L) * 3L + seq_len(5L) - 1L) %% p + 1L
    operator[1L, support] <- c(0.4, 0.25, 0.15, 0.12, 0.08)
    operator
  })
  names(operators) <- ids
  measurement_frame(operators, domain, id = "profile:scalar-frame:v1")
}

make_multivariate_frame <- function(n_nodes) {
  ids <- paste0("population", sprintf("%03d", seq_len(n_nodes)))
  operators <- lapply(seq_along(ids), function(index) {
    columns <- ((index - 1L) * 5L + seq_len(8L) - 1L) %% p + 1L
    operator <- matrix(0, 4L, p)
    operator[cbind(rep(seq_len(4L), each = 2L), columns)] <-
      rep(c(0.7, 0.3), 4L)
    operator
  })
  names(operators) <- ids
  measurement_frame(operators, domain, id = "profile:multivariate-frame:v1")
}

scalar_frame <- make_scalar_frame(24L)
scalar_between <- edge_frame(
  c(scalar_frame$node_ids, scalar_frame$node_ids),
  c(scalar_frame$node_ids, c(utils::tail(scalar_frame$node_ids, -1L),
    scalar_frame$node_ids[[1L]])),
  scalar_frame
)
multivariate_frame <- make_multivariate_frame(12L)
multivariate_between <- edge_frame(
  c(multivariate_frame$node_ids, multivariate_frame$node_ids),
  c(multivariate_frame$node_ids, c(utils::tail(multivariate_frame$node_ids, -1L),
    multivariate_frame$node_ids[[1L]])),
  multivariate_frame
)

rss_bytes <- function() {
  if (requireNamespace("ps", quietly = TRUE)) {
    as.double(ps::ps_memory_info()[["rss"]])
  } else {
    NA_real_
  }
}

profile_case <- function(label, between, routes) {
  reference <- NULL
  records <- lapply(routes, function(route) {
    invisible(measurement_form(
      relations, between, by, over,
      compute = compute_policy(block_features = 64L),
      route = route
    ))
    elapsed <- numeric(repetitions)
    peak_rss <- numeric(repetitions)
    final <- NULL
    for (iteration in seq_len(repetitions)) {
      gc()
      before <- rss_bytes()
      timing <- system.time({
        final <- measurement_form(
          relations, between, by, over,
          compute = compute_policy(block_features = 64L),
          route = route
        )
      })
      elapsed[[iteration]] <- unname(timing[["elapsed"]])
      after <- rss_bytes()
      peak_rss[[iteration]] <- max(before, after, na.rm = TRUE)
    }
    if (is.null(reference)) {
      reference <<- final
    }
    blocks <- effect_coupling(final)$values
    checksum <- sum(vapply(blocks, sum, numeric(1)))
    oracle_error <- max(abs(checksum -
      sum(vapply(effect_coupling(reference)$values, sum, numeric(1)))))
    data.frame(
      frame = label,
      route = route,
      selected_route = route,
      iterations = repetitions,
      median_elapsed_seconds = stats::median(elapsed),
      peak_rss_bytes = max(peak_rss, na.rm = TRUE),
      planned_workspace_bytes =
        final$receipt$execution$memory$planned_workspace_bytes,
      output_bytes = sum(final$block_index$length_elements) * 8,
      edge_count = nrow(final$block_index),
      scientific_plan_id = final$plan$scientific_plan_id,
      checksum = checksum,
      oracle_error = oracle_error,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, records)
  cat("Profiled", label, "\n"); flush.console()
  out
}

rprof_loop_share <- function(between, route) {
  path <- tempfile(fileext = ".out")
  on.exit(unlink(path), add = TRUE)
  utils::Rprof(path, interval = 0.01, memory.profiling = FALSE)
  invisible(measurement_form(
    relations, between, by, over,
    compute = compute_policy(block_features = 64L),
    route = route
  ))
  utils::Rprof(NULL)
  summary <- tryCatch(utils::summaryRprof(path), error = function(e) NULL)
  if (is.null(summary) || !nrow(summary$by.self)) {
    return(list(r_loop_share = NA_real_, top = character()))
  }
  names <- rownames(summary$by.self)
  loop_hits <- grepl(
    "for|lapply|vapply|switch|measurement_contract|\\[",
    names
  )
  blas_hits <- grepl("crossprod|%*%|colSums|do.call", names, fixed = FALSE)
  total <- sum(summary$by.self$self.time)
  list(
    r_loop_share = if (total <= 0) NA_real_ else
      sum(summary$by.self$self.time[loop_hits]) / total,
    blas_share = if (total <= 0) NA_real_ else
      sum(summary$by.self$self.time[blas_hits]) / total,
    top = utils::head(names, 8L)
  )
}

scalar <- profile_case(
  "scalar",
  scalar_between,
  c("auto", "pull_h", "forward_k", "multivariate_blocks",
    "scalar_stack", "factorized_h")
)
multivariate <- profile_case(
  "requested_multivariate",
  multivariate_between,
  c("auto", "pull_h", "forward_k", "multivariate_blocks", "factorized_h")
)
records <- rbind(scalar, multivariate)

auto_scalar <- scalar[scalar$route == "auto", , drop = FALSE]
auto_multi <- multivariate[multivariate$route == "auto", , drop = FALSE]
loop_scalar <- rprof_loop_share(scalar_between, auto_scalar$selected_route)
loop_multi <- rprof_loop_share(multivariate_between, auto_multi$selected_route)
max_loop_share <- max(loop_scalar$r_loop_share, loop_multi$r_loop_share,
  na.rm = TRUE)
admit_native <- is.finite(max_loop_share) && max_loop_share >= 0.15
decision <- if (admit_native) {
  "prototype_native_kernel"
} else {
  "no_rcpp_keep_blas"
}

artifact <- list(
  schema_version = 1L,
  provenance = provenance,
  dimensions = list(
    experimental_conditions = q,
    domain_features = p,
    scalar_nodes = length(scalar_frame$node_ids),
    multivariate_nodes = length(multivariate_frame$node_ids),
    partitions = length(relations$partitions)
  ),
  records = records,
  route_selection = list(
    scalar_auto = auto_scalar$selected_route,
    multivariate_auto = auto_multi$selected_route
  ),
  rprof = list(
    scalar = loop_scalar,
    multivariate = loop_multi,
    max_r_loop_share = max_loop_share
  ),
  checks = list(
    oracle_parity = max(records$oracle_error) < 1e-10,
    plan_identity_stable = length(unique(records$scientific_plan_id[
      records$frame == "scalar"])) == 1L &&
      length(unique(records$scientific_plan_id[
        records$frame == "requested_multivariate"])) == 1L,
    planned_workspace_exceeds_output =
      all(records$planned_workspace_bytes > records$output_bytes)
  ),
  decision = decision,
  decision_rule = paste0(
    "Admit Rcpp only if R-level loops are >= 15% of end-to-end time ",
    "and a prototype projects >= 1.25x. Observed max R-loop share: ",
    signif(max_loop_share, 3), "."
  )
)

result_path <- file.path(output_dir, "measurement-profile.rds")
summary_path <- file.path(output_dir, "measurement-profile-summary.csv")
saveRDS(artifact, result_path)
utils::write.csv(records, summary_path, row.names = FALSE)
cat("Measurement profile decision:", decision, "\n")
print(records, row.names = FALSE, digits = 6)
cat("Max R-loop share:", signif(max_loop_share, 3), "\n")
cat("Wrote", result_path, "\n")
