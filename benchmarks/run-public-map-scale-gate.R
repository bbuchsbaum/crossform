#!/usr/bin/env Rscript

# Public map-scale gate for the point-geometry paths taught in the newcomer
# workflow. The explicit domain-wide identity metric exercises the path that
# previously fell off a validation cliff; the implicit identity plan is its
# matched estimand. Both are checked against a direct cross-partition oracle.

if (!requireNamespace("parallel", quietly = TRUE) ||
    !requireNamespace("ps", quietly = TRUE) ||
    !requireNamespace("devtools", quietly = TRUE)) {
  stop("The public map gate requires parallel, ps, and devtools.")
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
  stop("The public map gate requires at least three repetitions.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
result_path <- file.path(output_dir, "public-map-scale-gate.rds")
summary_path <- file.path(output_dir, "public-map-scale-gate-summary.csv")
ready_path <- file.path(output_dir, "public-map-scale-gate.ready")
unlink(c(result_path, ready_path))

run_public_map_worker <- function(repo, result_path, ready_path, repetitions) {
  suppressPackageStartupMessages(devtools::load_all(repo, quiet = TRUE))
  set.seed(2026081401L)
  dimensions <- c(9L, 8L, 8L)
  domain <- volume_domain(
    array(TRUE, dimensions), spacing = c(3, 3, 3),
    id = "public-map-scale-576"
  )
  condition_names <- sprintf("condition%02d", seq_len(12L))
  partition_names <- sprintf("run%02d", seq_len(8L))
  truth <- matrix(
    stats::rnorm(length(condition_names) * domain$n_features, sd = 0.3),
    length(condition_names), domain$n_features,
    dimnames = list(condition_names, NULL)
  )
  sources <- stats::setNames(lapply(partition_names, function(partition) {
    value <- truth + matrix(
      stats::rnorm(length(truth), sd = 0.5),
      nrow(truth), ncol(truth)
    )
    rownames(value) <- condition_names
    value
  }), partition_names)
  relation <- relation(
    sources,
    effects = effect_space(
      condition_names, basis_id = "public-map-scale:v1",
      units = "arbitrary-BOLD"
    ),
    domain = domain
  )
  frame <- compile_frame(
    searchlights(radius = 4, normalization = "local"), domain
  )
  over <- cross_partitions(relation)
  implicit_plan <- plan_geometry(relation, frame, over)

  metric_elapsed <- system.time({
    explicit_metric <- neural_metric(
      diag(domain$n_features), domain, inverse = diag(domain$n_features),
      provenance = list(source = "public-map-fixed-identity")
    )
  })[["elapsed"]]
  explicit_plan <- plan_geometry(
    relation, frame, over, metric = explicit_metric
  )

  # A second relation retains raw residuals so the public analytic-covariance
  # path can be swept over a complete map at the same 12-condition scale.
  design <- stats::model.matrix(~ 0 + factor(
    rep(condition_names, each = 3L), levels = condition_names
  ))
  colnames(design) <- condition_names
  targets <- diag(length(condition_names))
  dimnames(targets) <- list(condition_names, condition_names)
  raw_sources <- stats::setNames(lapply(partition_names, function(partition) {
    design %*% truth + matrix(
      stats::rnorm(nrow(design) * domain$n_features, sd = 0.5),
      nrow(design), domain$n_features
    )
  }), partition_names)
  sampling_fit <- lm_relation_fit(
    raw_sources, design, targets,
    effect_names = effect_space(
      condition_names, basis_id = "public-map-sampling:v1",
      units = "arbitrary-BOLD"
    ),
    domain = domain, sampling_unit = "trial"
  )
  sampling_frame <- compile_frame(
    searchlights(radius = 7, normalization = "local"), domain
  )
  sampling_plan <- plan_geometry(
    sampling_fit$relation, sampling_frame,
    cross_partitions(sampling_fit$relation)
  )

  selected_nodes <- unique(as.integer(round(seq.int(
    1L, nrow(frame$weights), length.out = 16L
  ))))
  condition_pairs <- utils::combn(seq_along(condition_names), 2L)
  partition_pairs <- utils::combn(seq_along(partition_names), 2L)
  oracle <- matrix(
    NA_real_, length(selected_nodes), ncol(condition_pairs),
    dimnames = list(selected_nodes, apply(
      condition_pairs, 2L,
      function(pair) paste(condition_names[pair], collapse = " - ")
    ))
  )
  for (node_index in seq_along(selected_nodes)) {
    node <- selected_nodes[[node_index]]
    weights <- as.numeric(frame$weights[node, ])
    support <- which(weights != 0)
    for (pair_index in seq_len(ncol(condition_pairs))) {
      condition_pair <- condition_pairs[, pair_index]
      edge_values <- vapply(seq_len(ncol(partition_pairs)), function(edge) {
        left <- sources[[partition_pairs[1L, edge]]]
        right <- sources[[partition_pairs[2L, edge]]]
        left_difference <-
          left[condition_pair[1L], support] - left[condition_pair[2L], support]
        right_difference <-
          right[condition_pair[1L], support] - right[condition_pair[2L], support]
        sum(weights[support] * left_difference * right_difference)
      }, numeric(1))
      oracle[node_index, pair_index] <- mean(edge_values)
    }
  }

  run_map <- function(plan) as.matrix(rdm(plan)$values)
  run_sampling_nodes <- function(nodes) {
    values <- lapply(nodes, function(node) {
      covariance <- rdm_sampling_covariance(
        sampling_plan, sampling_fit, target = "plugin", at = node
      )
      sampling_covariance(covariance)
    })
    value <- do.call(rbind, values)
    rownames(value) <- nodes
    value
  }
  warm_implicit <- run_map(implicit_plan)
  warm_explicit <- run_map(explicit_plan)
  oracle_error <- max(abs(
    warm_implicit[selected_nodes, , drop = FALSE] - oracle
  ))
  matched_error <- max(abs(warm_explicit - warm_implicit))
  if (!is.finite(oracle_error) || oracle_error > 1e-12 ||
      !is.finite(matched_error) || matched_error > 1e-12) {
    stop("Public map warm-up failed its numerical oracle.")
  }
  sampling_probe_nodes <- unique(as.integer(round(seq.int(
    1L, nrow(sampling_frame$weights), length.out = 25L
  ))))
  warm_sampling <- run_sampling_nodes(sampling_probe_nodes)
  if (any(!is.finite(warm_sampling)) || any(warm_sampling < -1e-12)) {
    stop("Public sampling-covariance warm-up returned invalid variances.")
  }

  baseline_rss <- unname(ps::ps_memory_info(
    ps::ps_handle(Sys.getpid())
  )[["rss"]])
  file.create(ready_path)
  timings <- vector("list", 2L * repetitions)
  cursor <- 0L
  for (iteration in seq_len(repetitions)) {
    order <- if (iteration %% 2L) {
      c("implicit_identity", "explicit_identity_metric")
    } else {
      c("explicit_identity_metric", "implicit_identity")
    }
    for (path in order) {
      cursor <- cursor + 1L
      plan <- if (identical(path, "implicit_identity")) {
        implicit_plan
      } else {
        explicit_plan
      }
      elapsed <- system.time(value <- run_map(plan))[["elapsed"]]
      reference <- if (identical(path, "implicit_identity")) {
        warm_implicit
      } else {
        warm_explicit
      }
      timings[[cursor]] <- data.frame(
        iteration = iteration,
        order = match(path, order),
        path = path,
        elapsed_seconds = unname(elapsed),
        max_abs_repeat_error = max(abs(value - reference)),
        stringsAsFactors = FALSE
      )
      rm(value)
      gc()
    }
  }
  timings <- do.call(rbind, timings)
  sampling_timings <- vector("list", repetitions)
  for (iteration in seq_len(repetitions)) {
    elapsed <- system.time(value <- run_sampling_nodes(
      sampling_probe_nodes
    ))[["elapsed"]]
    sampling_timings[[iteration]] <- data.frame(
      iteration = iteration,
      nodes = length(sampling_probe_nodes),
      elapsed_seconds = unname(elapsed),
      seconds_per_node = unname(elapsed) / length(sampling_probe_nodes),
      max_abs_repeat_error = max(abs(value - warm_sampling)),
      stringsAsFactors = FALSE
    )
    rm(value)
    gc()
  }
  sampling_timings <- do.call(rbind, sampling_timings)
  full_sampling_elapsed <- system.time(
    full_sampling <- run_sampling_nodes(seq_len(nrow(sampling_frame$weights)))
  )[["elapsed"]]
  full_sampling_probe_error <- max(abs(
    full_sampling[sampling_probe_nodes, , drop = FALSE] - warm_sampling
  ))
  full_sampling_finite_nonnegative <- all(is.finite(full_sampling)) &&
    all(full_sampling >= -1e-12)
  result <- list(
    schema_version = 1L,
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    estimand = paste(
      "equal-weight mean cross-partition squared Euclidean RDM;",
      "local 1/P searchlight normalization"
    ),
    dimensions = list(
      volume = dimensions,
      features = domain$n_features,
      frame_nodes = nrow(frame$weights),
      conditions = length(condition_names),
      rdm_coordinates = ncol(condition_pairs),
      partitions = length(partition_names),
      partition_edges = ncol(partition_pairs),
      selected_oracle_nodes = length(selected_nodes),
      support_min = min(frame$support_index$cost$support_size),
      support_mean = unname(frame$support_index$cost$support_size[["mean"]]),
      support_max = max(frame$support_index$cost$support_size)
    ),
    sampling_dimensions = list(
      frame_nodes = nrow(sampling_frame$weights),
      conditions = length(condition_names),
      covariance_dimension = ncol(warm_sampling),
      partitions = length(partition_names),
      probe_nodes = length(sampling_probe_nodes),
      support_min = min(sampling_frame$support_index$cost$support_size),
      support_mean = unname(
        sampling_frame$support_index$cost$support_size[["mean"]]
      ),
      support_max = max(sampling_frame$support_index$cost$support_size)
    ),
    plan = list(
      implicit_id = implicit_plan$scientific_plan_id,
      explicit_id = explicit_plan$scientific_plan_id,
      explicit_metric_seconds = unname(metric_elapsed),
      implicit_dense_payload_bytes = implicit_plan$dense_payload_bytes,
      explicit_dense_payload_bytes = explicit_plan$dense_payload_bytes
    ),
    numerical = list(
      independent_oracle_max_abs_error = oracle_error,
      explicit_implicit_max_abs_error = matched_error,
      repeat_max_abs_error = max(timings$max_abs_repeat_error)
    ),
    timings = timings,
    sampling = list(
      target = "partition_mean_plugin",
      probe_timings = sampling_timings,
      full_sweep_seconds = unname(full_sampling_elapsed),
      full_sweep_seconds_per_node =
        unname(full_sampling_elapsed) / nrow(sampling_frame$weights),
      probe_repeat_max_abs_error = max(
        sampling_timings$max_abs_repeat_error
      ),
      full_sweep_probe_max_abs_error = full_sampling_probe_error,
      full_sweep_finite_nonnegative = full_sampling_finite_nonnegative
    ),
    memory = list(
      baseline_rss_bytes = baseline_rss
    )
  )
  saveRDS(result, result_path)
  TRUE
}

job <- parallel::mcparallel(
  run_public_map_worker(repo, result_path, ready_path, repetitions),
  silent = TRUE
)
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
  stop("The isolated public map-scale worker failed.")
}

result <- readRDS(result_path)
peak_rss <- max(peak_rss, result$memory$baseline_rss_bytes, na.rm = TRUE)
result$memory$os_peak_rss_bytes <- peak_rss
result$memory$incremental_peak_rss_bytes <- max(
  0, peak_rss - result$memory$baseline_rss_bytes
)
medians <- stats::aggregate(
  elapsed_seconds ~ path, result$timings, stats::median
)
implicit_seconds <- medians$elapsed_seconds[
  medians$path == "implicit_identity"
]
explicit_seconds <- medians$elapsed_seconds[
  medians$path == "explicit_identity_metric"
]
ratio <- explicit_seconds / max(implicit_seconds, 0.001)
result$gate <- list(
  numerical_tolerance = 1e-12,
  maximum_path_seconds = 60,
  maximum_explicit_to_implicit_ratio = 5,
  maximum_sampling_probe_seconds_per_node = 1,
  maximum_sampling_full_sweep_seconds = 10 * 60,
  maximum_incremental_rss_bytes = 1024^3,
  numerical_pass =
    result$numerical$independent_oracle_max_abs_error <= 1e-12 &&
    result$numerical$explicit_implicit_max_abs_error <= 1e-12 &&
    result$numerical$repeat_max_abs_error <= 1e-12,
  implicit_runtime_pass = implicit_seconds <= 60,
  explicit_runtime_pass = explicit_seconds <= 60,
  relative_runtime_pass = ratio <= 5,
  sampling_numerical_pass =
    result$sampling$probe_repeat_max_abs_error <= 1e-12 &&
    result$sampling$full_sweep_probe_max_abs_error <= 1e-12 &&
    isTRUE(result$sampling$full_sweep_finite_nonnegative),
  sampling_runtime_pass =
    stats::median(result$sampling$probe_timings$seconds_per_node) <= 1 &&
    result$sampling$full_sweep_seconds <= 10 * 60,
  memory_pass = result$memory$incremental_peak_rss_bytes <= 1024^3
)
result$gate$passed <- all(unlist(result$gate[c(
  "numerical_pass", "implicit_runtime_pass", "explicit_runtime_pass",
  "relative_runtime_pass", "sampling_numerical_pass",
  "sampling_runtime_pass", "memory_pass"
)]))
result$summary <- data.frame(
  features = result$dimensions$features,
  frame_nodes = result$dimensions$frame_nodes,
  conditions = result$dimensions$conditions,
  partitions = result$dimensions$partitions,
  implicit_median_seconds = implicit_seconds,
  explicit_median_seconds = explicit_seconds,
  explicit_to_implicit_ratio = ratio,
  sampling_probe_median_seconds_per_node = stats::median(
    result$sampling$probe_timings$seconds_per_node
  ),
  sampling_full_sweep_seconds = result$sampling$full_sweep_seconds,
  sampling_full_sweep_seconds_per_node =
    result$sampling$full_sweep_seconds_per_node,
  independent_oracle_max_abs_error =
    result$numerical$independent_oracle_max_abs_error,
  explicit_implicit_max_abs_error =
    result$numerical$explicit_implicit_max_abs_error,
  incremental_peak_rss_bytes = result$memory$incremental_peak_rss_bytes,
  passed = result$gate$passed,
  stringsAsFactors = FALSE
)
saveRDS(result, result_path, compress = "xz")
utils::write.csv(result$summary, summary_path, row.names = FALSE)
unlink(ready_path)
print(result$summary, row.names = FALSE, digits = 7)
if (!isTRUE(result$gate$passed)) {
  stop("The public map-scale gate did not pass.")
}
