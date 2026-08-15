#!/usr/bin/env Rscript

# Query-first scale gate (Gate 5): certifies that selected contrasts, selected
# RDM edges, and RSA coefficients execute at large condition count without
# materializing the full q(q+1)/2 geometry field, that the fused route is not
# slower than materialize-then-project, and that its memory stays far below
# the dense query matrix a packed representation would require
# (q = 100 implies a 5050-by-4950 dense query of about 200 MiB, plus an
# 87 MiB two-component geometry field; the fused route allocates neither).

if (!requireNamespace("parallel", quietly = TRUE) ||
    !requireNamespace("ps", quietly = TRUE) ||
    !requireNamespace("devtools", quietly = TRUE)) {
  stop("The query-first gate requires parallel, ps, and devtools.")
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
  3L
}
if (is.na(repetitions) || repetitions < 3L) {
  stop("The query-first gate requires at least three repetitions.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
result_path <- file.path(output_dir, "query-first-scale-gate.rds")
summary_path <- file.path(output_dir, "query-first-scale-gate-summary.csv")
ready_path <- file.path(output_dir, "query-first-scale-gate.ready")
unlink(c(result_path, ready_path))

run_query_first_worker <- function(repo, result_path, ready_path,
                                   repetitions, benchmark_paths) {
  suppressPackageStartupMessages(devtools::load_all(repo, quiet = TRUE))
  set.seed(2026081501L)
  dimensions <- c(12L, 10L, 9L)
  domain <- volume_domain(
    array(TRUE, dimensions), spacing = c(3, 3, 3),
    id = "query-first-scale-1080"
  )
  condition_names <- sprintf("condition%03d", seq_len(100L))
  partition_names <- sprintf("run%02d", seq_len(8L))
  q <- length(condition_names)
  design <- stats::model.matrix(~ 0 + factor(
    rep(condition_names, each = 2L), levels = condition_names
  ))
  colnames(design) <- condition_names
  targets <- diag(q)
  dimnames(targets) <- list(condition_names, condition_names)
  truth <- matrix(
    stats::rnorm(q * domain$n_features, sd = 0.3),
    q, domain$n_features, dimnames = list(condition_names, NULL)
  )
  raw_sources <- stats::setNames(lapply(partition_names, function(partition) {
    design %*% truth + matrix(
      stats::rnorm(nrow(design) * domain$n_features, sd = 0.5),
      nrow(design), domain$n_features
    )
  }), partition_names)
  fit <- lm_relation_fit(
    raw_sources, design, targets,
    effect_names = effect_space(
      condition_names, basis_id = "query-first-scale:v1",
      units = "arbitrary-BOLD"
    ),
    domain = domain, sampling_unit = "trial"
  )
  relation <- fit$relation
  frame <- compile_frame(
    searchlights(radius = 4, normalization = "local"), domain
  )
  plan <- plan_geometry(
    relation, frame,
    cross_partitions(relation, independence = "independent")
  )

  # Selected edges: one hundred pairs spread across pair space.
  all_pairs <- t(utils::combn(seq_len(q), 2L))
  selected_rows <- unique(as.integer(round(seq.int(
    1L, nrow(all_pairs), length.out = 100L
  ))))
  selected_pairs <- all_pairs[selected_rows, , drop = FALSE]

  # A fixed model RDM for the RSA readout.
  positions <- seq_len(q)
  model <- abs(outer(positions, positions, "-"))
  diag(model) <- 0
  dimnames(model) <- list(condition_names, condition_names)

  run_selected <- function() {
    rdm(plan, pairs = selected_pairs)
  }
  run_full_fused <- function() {
    rdm(plan)
  }
  run_rsa <- function() {
    rsa(plan, models = list(distance = model))
  }
  run_materialized <- function() {
    rdm(geometry(plan))
  }

  warm_selected <- run_selected()
  warm_full <- run_full_fused()
  warm_rsa <- run_rsa()
  warm_materialized <- run_materialized()

  # Independent oracle over probe nodes and the selected pairs: plain OLS
  # betas and a direct cross-partition loop, blind to every compiled path.
  effect_matrices <- stats::setNames(lapply(partition_names, function(name) {
    qr.coef(qr(design), raw_sources[[name]])
  }), partition_names)
  probe_nodes <- unique(as.integer(round(seq.int(
    1L, nrow(frame$weights), length.out = 12L
  ))))
  partition_pairs <- utils::combn(seq_along(partition_names), 2L)
  oracle <- matrix(NA_real_, length(probe_nodes), nrow(selected_pairs))
  for (node_index in seq_along(probe_nodes)) {
    node <- probe_nodes[[node_index]]
    weights <- as.numeric(frame$weights[node, ])
    support <- which(weights != 0)
    for (pair_index in seq_len(nrow(selected_pairs))) {
      pair <- selected_pairs[pair_index, ]
      edge_values <- vapply(seq_len(ncol(partition_pairs)), function(edge) {
        left <- effect_matrices[[partition_pairs[1L, edge]]]
        right <- effect_matrices[[partition_pairs[2L, edge]]]
        left_difference <- left[pair[[1L]], support] - left[pair[[2L]], support]
        right_difference <- right[pair[[1L]], support] -
          right[pair[[2L]], support]
        sum(weights[support] * left_difference * right_difference)
      }, numeric(1))
      oracle[node_index, pair_index] <- mean(edge_values)
    }
  }
  oracle_error <- max(abs(
    as.matrix(warm_selected$values)[probe_nodes, , drop = FALSE] - oracle
  ))
  selected_column_match <- max(abs(
    as.matrix(warm_selected$values) -
      as.matrix(warm_full$values)[, selected_rows, drop = FALSE]
  ))
  route_error <- max(abs(
    as.matrix(warm_full$values) - as.matrix(warm_materialized$values)
  ))
  route_identity_stable <- identical(
    warm_full$receipt$scientific_plan_id,
    warm_materialized$receipt$scientific_plan_id
  )
  if (!is.finite(oracle_error) || oracle_error > 1e-10 ||
      selected_column_match > 1e-12 || route_error > 1e-10 ||
      !route_identity_stable) {
    stop("Query-first warm-up failed its numerical or identity oracle.")
  }

  baseline_rss <- unname(ps::ps_memory_info(
    ps::ps_handle(Sys.getpid())
  )[["rss"]])
  file.create(ready_path)
  runners <- list(
    selected_100 = run_selected,
    full_fused = run_full_fused,
    rsa_fused = run_rsa,
    materialized = run_materialized
  )
  if (!is.character(benchmark_paths) || !length(benchmark_paths) ||
      any(!benchmark_paths %in% names(runners)) ||
      anyDuplicated(benchmark_paths)) {
    stop("The query-first worker received invalid benchmark paths.")
  }
  timings <- vector("list", length(benchmark_paths) * repetitions)
  cursor <- 0L
  for (iteration in seq_len(repetitions)) {
    order <- if (iteration %% 2L) benchmark_paths else rev(benchmark_paths)
    for (path in order) {
      cursor <- cursor + 1L
      elapsed <- system.time(value <- runners[[path]]())[["elapsed"]]
      timings[[cursor]] <- data.frame(
        iteration = iteration,
        path = path,
        elapsed_seconds = unname(elapsed),
        stringsAsFactors = FALSE
      )
      rm(value)
      gc()
    }
  }
  timings <- do.call(rbind, timings)

  fused_receipt <- warm_full$receipt
  result <- list(
    schema_version = 1L,
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    estimand = paste(
      "equal-weight mean cross-partition squared Euclidean RDM and fixed",
      "linear RSA readout at q = 100; query-first structured execution"
    ),
    dimensions = list(
      volume = dimensions,
      features = domain$n_features,
      frame_nodes = nrow(frame$weights),
      conditions = q,
      rdm_coordinates = nrow(all_pairs),
      selected_pairs = nrow(selected_pairs),
      partitions = length(partition_names),
      trial_rows_per_partition = nrow(design),
      dense_query_bytes_avoided = 8 * (q * (q + 1L) / 2L) * nrow(all_pairs),
      dense_geometry_bytes_avoided =
        2 * 8 * nrow(frame$weights) * (q * (q + 1L) / 2L)
    ),
    identity = list(
      plan_id = plan$scientific_plan_id,
      full_view_id = warm_full$receipt$scientific_plan_id,
      route_identity_stable = route_identity_stable,
      fused_kernel = fused_receipt$kernel_version,
      fused_lowering = fused_receipt$execution_plan$lowering,
      fused_materialization = fused_receipt$execution_plan$materialization,
      planned_workspace_bytes =
        fused_receipt$memory$planned_workspace_bytes
    ),
    numerical = list(
      independent_oracle_max_abs_error = oracle_error,
      selected_column_max_abs_error = selected_column_match,
      fused_materialized_max_abs_error = route_error
    ),
    timings = timings,
    memory = list(baseline_rss_bytes = baseline_rss)
  )
  saveRDS(result, result_path)
  TRUE
}

monitor_worker <- function(job, result_path, ready_path) {
  handle <- ps::ps_handle(job$pid)
  measurement_started <- FALSE
  peak_rss <- 0
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
    stop("The isolated query-first scale worker failed.")
  }
  result <- readRDS(result_path)
  peak_rss <- max(peak_rss, result$memory$baseline_rss_bytes, na.rm = TRUE)
  list(result = result, peak_rss = peak_rss)
}

all_paths <- c("selected_100", "full_fused", "rsa_fused", "materialized")
job <- parallel::mcparallel(
  run_query_first_worker(
    repo, result_path, ready_path, repetitions, all_paths
  ),
  silent = TRUE
)
comparison <- monitor_worker(job, result_path, ready_path)
result <- comparison$result

# Process-level RSS is allocator-sensitive when the materialized comparator
# runs in the same timing pool. Certify the public query-first paths in a
# second fresh child after all routes have been warmed, and retain the mixed
# process only as a diagnostic. This makes the memory claim match its subject.
memory_result_path <- file.path(output_dir, "query-first-scale-memory.rds")
memory_ready_path <- file.path(output_dir, "query-first-scale-memory.ready")
unlink(c(memory_result_path, memory_ready_path))
query_first_paths <- c("selected_100", "full_fused", "rsa_fused")
memory_job <- parallel::mcparallel(
  run_query_first_worker(
    repo, memory_result_path, memory_ready_path, repetitions,
    query_first_paths
  ),
  silent = TRUE
)
query_first_memory <- monitor_worker(
  memory_job, memory_result_path, memory_ready_path
)
comparison_baseline <- result$memory$baseline_rss_bytes
query_first_baseline <- query_first_memory$result$memory$baseline_rss_bytes
result$memory <- list(
  scope = "query_first_paths_only_after_all_route_warmup",
  baseline_rss_bytes = query_first_baseline,
  os_peak_rss_bytes = query_first_memory$peak_rss,
  incremental_peak_rss_bytes = max(
    0, query_first_memory$peak_rss - query_first_baseline
  ),
  comparison_process_baseline_rss_bytes = comparison_baseline,
  comparison_process_peak_rss_bytes = comparison$peak_rss,
  comparison_process_incremental_rss_bytes = max(
    0, comparison$peak_rss - comparison_baseline
  )
)
unlink(c(memory_result_path, memory_ready_path))
medians <- stats::aggregate(
  elapsed_seconds ~ path, result$timings, stats::median
)
median_of <- function(path) {
  medians$elapsed_seconds[medians$path == path]
}
selected_seconds <- median_of("selected_100")
full_seconds <- median_of("full_fused")
rsa_seconds <- median_of("rsa_fused")
materialized_seconds <- median_of("materialized")
fused_ratio <- full_seconds / max(materialized_seconds, 0.001)
result$gate <- list(
  numerical_tolerance = 1e-10,
  maximum_fused_to_materialized_ratio = 1.2,
  maximum_selected_to_full_ratio = 1.0,
  maximum_path_seconds = 120,
  # The dense packed query alone would be ~191 MiB and the materialized
  # two-component geometry field ~83 MiB more. The separately isolated
  # query-first pool must remain within this absolute map-scale budget; the
  # materialized timing comparator is intentionally outside this memory claim.
  maximum_incremental_rss_bytes = 512 * 1024^2,
  numerical_pass =
    result$numerical$independent_oracle_max_abs_error <= 1e-10 &&
    result$numerical$selected_column_max_abs_error <= 1e-12 &&
    result$numerical$fused_materialized_max_abs_error <= 1e-10,
  identity_pass = isTRUE(result$identity$route_identity_stable) &&
    identical(result$identity$fused_kernel, "additive-query-fused-v2"),
  fused_not_slower_pass = fused_ratio <= 1.2,
  selected_cheaper_pass = selected_seconds <= full_seconds,
  runtime_pass = max(
    selected_seconds, full_seconds, rsa_seconds, materialized_seconds
  ) <= 120,
  memory_pass = result$memory$incremental_peak_rss_bytes <= 512 * 1024^2
)
result$gate$passed <- all(unlist(result$gate[c(
  "numerical_pass", "identity_pass", "fused_not_slower_pass",
  "selected_cheaper_pass", "runtime_pass", "memory_pass"
)]))
result$summary <- data.frame(
  features = result$dimensions$features,
  frame_nodes = result$dimensions$frame_nodes,
  conditions = result$dimensions$conditions,
  rdm_coordinates = result$dimensions$rdm_coordinates,
  selected_pairs = result$dimensions$selected_pairs,
  selected_median_seconds = selected_seconds,
  full_fused_median_seconds = full_seconds,
  rsa_median_seconds = rsa_seconds,
  materialized_median_seconds = materialized_seconds,
  fused_to_materialized_ratio = fused_ratio,
  independent_oracle_max_abs_error =
    result$numerical$independent_oracle_max_abs_error,
  incremental_peak_rss_bytes = result$memory$incremental_peak_rss_bytes,
  planned_workspace_bytes = result$identity$planned_workspace_bytes,
  passed = result$gate$passed,
  stringsAsFactors = FALSE
)
saveRDS(result, result_path, compress = "xz")
utils::write.csv(result$summary, summary_path, row.names = FALSE)
unlink(ready_path)
print(result$summary, row.names = FALSE, digits = 7)
if (!isTRUE(result$gate$passed)) {
  stop("The query-first scale gate did not pass.")
}
