run_crossnobis_scale_gate <- function(repo, result_path, ready_path,
                                      full_dim = c(52L, 42L, 25L),
                                      active_slices = 24L) {
  repo <- normalizePath(repo, mustWork = TRUE)
  if (!requireNamespace("pkgload", quietly = TRUE) ||
      !requireNamespace("ps", quietly = TRUE) ||
      !requireNamespace("neuroim2", quietly = TRUE)) {
    stop("Scale qualification requires pkgload, ps, and neuroim2.")
  }
  pkgload::load_all(repo, quiet = TRUE)

  seed <- 2026081308L
  set.seed(seed)
  if (!is.numeric(full_dim) || length(full_dim) != 3L ||
      anyNA(full_dim) || any(!is.finite(full_dim)) ||
      any(full_dim < 2L) || any(full_dim %% 1 != 0) ||
      !is.numeric(active_slices) || length(active_slices) != 1L ||
      is.na(active_slices) || !is.finite(active_slices) ||
      active_slices %% 1 != 0 || active_slices < 1L ||
      active_slices > full_dim[[3L]]) {
    stop("Scale-fixture dimensions and active slices are invalid.")
  }
  full_dim <- as.integer(full_dim)
  active_slices <- as.integer(active_slices)
  mask_array <- array(FALSE, full_dim)
  mask_array[, , seq_len(active_slices)] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(
    mask_array,
    neuroim2::NeuroSpace(full_dim, spacing = c(3, 3, 3))
  )
  domain <- neuroim2_volume_domain(mask, id = "crossnobis-scale-52k:v1")
  observations <- 32L
  partitions <- c("run1", "run2", "run3")
  condition <- rep(c(-0.5, 0.5), each = observations / 2L)
  design <- cbind(intercept = 1, condition = condition)
  effects <- rbind(
    baseline = c(1, 0),
    condition = c(0, 1)
  )
  signal <- numeric(domain$n_features)
  center <- ceiling(domain$n_features / 2)
  signal[pmax(1L, center - 24L):pmin(domain$n_features, center + 24L)] <-
    0.35
  raw_started <- proc.time()[["elapsed"]]
  raw <- lapply(partitions, function(partition) {
    noise <- matrix(
      stats::rnorm(observations * domain$n_features),
      observations, domain$n_features
    )
    design %*% rbind(rep(0, domain$n_features), signal) + noise
  })
  names(raw) <- partitions
  raw_seconds <- proc.time()[["elapsed"]] - raw_started

  fit_started <- proc.time()[["elapsed"]]
  fit <- lm_relation_fit(
    raw, design, effects, domain = domain,
    provenance = list(
      fixture = "crossnobis-scale-52k:v1",
      generator = "iid_gaussian_with_local_effect"
    )
  )
  fit_seconds <- proc.time()[["elapsed"]] - fit_started
  frame_started <- proc.time()[["elapsed"]]
  frame <- compile_frame(
    searchlights(radius = 6.1, normalization = "local"),
    domain
  )
  frame_seconds <- proc.time()[["elapsed"]] - frame_started
  over <- pairing("run1", "run2", independence = "independent")
  recipe <- shrinkage_precision(
    shrinkage = 0.2,
    relative_variance_floor = 1e-8,
    relative_spectral_floor = 1e-10
  )
  training <- metric_training_policy("exclude_evaluation")
  compute <- compute_policy(workspace_bytes = 4 * 1024^3)

  gc()
  baseline_rss <- effectagram:::.current_rss_bytes()
  writeLines("ready", ready_path, useBytes = TRUE)

  plan_started <- proc.time()[["elapsed"]]
  plan <- plan_crossnobis(
    fit, frame, over,
    metric = recipe,
    training = training,
    compute = compute,
    residual_workspace_bytes = 4 * 1024^3
  )
  plan_seconds <- proc.time()[["elapsed"]] - plan_started
  evaluation_started <- proc.time()[["elapsed"]]
  value <- crossnobis(plan, c(baseline = 0, condition = 1))
  evaluation_seconds <- proc.time()[["elapsed"]] - evaluation_started
  mapping_started <- proc.time()[["elapsed"]]
  mapped <- as_neurovol(
    value$values, mask, domain, fill = NA_real_,
    label = "learned local crossnobis"
  )
  mapping_seconds <- proc.time()[["elapsed"]] - mapping_started
  mapped_array <- as.array(mapped)
  mapping_exact <- identical(
    unname(mapped_array[domain$feature_ids]),
    unname(as.double(value$values))
  ) && all(is.na(mapped_array[-domain$feature_ids])) &&
    identical(neuroim2::space(mapped), neuroim2::space(mask))
  gc()
  rss_after <- effectagram:::.current_rss_bytes()

  index <- frame$support_index
  support_sizes <- diff(index$ptr)
  residual_reads <- vapply(
    plan$metric_schedule$statistics$execution$atomic,
    `[[`, integer(1), "residual_reads"
  )
  evaluation_reads <- value$metadata$source_session$read_count
  analysis_seconds <- fit_seconds + frame_seconds + plan_seconds +
    evaluation_seconds + mapping_seconds
  result <- list(
    schema_version = 1L,
    fixture = list(
      id = "crossnobis-scale-52k:v1",
      seed = seed,
      full_dim = full_dim,
      active_slices = active_slices,
      features = domain$n_features,
      spacing_mm = c(3, 3, 3),
      radius_mm = 6.1,
      observations_per_partition = observations,
      partitions = partitions,
      evaluation_edges = nrow(over),
      training_policy = training$kind,
      metric = recipe$kind
    ),
    support = list(
      min = min(support_sizes),
      median = stats::median(support_sizes),
      mean = mean(support_sizes),
      max = max(support_sizes),
      memberships = length(index$members),
      union_pair_stored_nnz = length(index$pair_pattern@i),
      union_degree = index$cost$union_degree,
      structural_bytes = index$cost$estimated_structural_bytes
    ),
    work = list(
      local_metric_derivations =
        as.double(length(index$node_ids) * nrow(over)),
      dense_metric_entries =
        as.double(index$cost$dense_metric_entries * nrow(over)),
      factorization_units =
        as.double(index$cost$one_dense_factorization_pass_units *
          nrow(over)),
      pair_atoms_materialized =
        value$metadata$diagnostics$pair_atoms_materialized,
      pair_frame_materialized =
        value$metadata$diagnostics$pair_frame_materialized,
      factor_table_retained =
        value$metadata$diagnostics$metric_factor_table_retained
    ),
    reads = list(
      residual_by_partition = residual_reads,
      residual_total = sum(residual_reads),
      evaluation_by_partition = evaluation_reads,
      evaluation_total = sum(evaluation_reads)
    ),
    timing = list(
      raw_fixture_seconds = raw_seconds,
      fit_seconds = fit_seconds,
      frame_seconds = frame_seconds,
      plan_and_residual_statistics_seconds = plan_seconds,
      crossnobis_seconds = evaluation_seconds,
      output_mapping_seconds = mapping_seconds,
      analysis_seconds = analysis_seconds
    ),
    memory = list(
      baseline_rss_bytes = baseline_rss,
      rss_after_bytes = rss_after,
      planned_workspace_bytes = plan$memory$planned_workspace_bytes,
      budget_bytes = plan$memory$budget_bytes,
      frame_object_bytes = as.double(utils::object.size(frame)),
      residual_statistics_object_bytes = as.double(utils::object.size(
        plan$metric_schedule$statistics
      )),
      result_object_bytes = as.double(utils::object.size(value)),
      persistent_factor_table_bytes = 0
    ),
    output = list(
      exact_full_index_mapping = mapping_exact,
      finite_values = all(is.finite(value$values)),
      negative_values_retained = any(value$values < 0),
      output_values = length(value$values),
      compact_checksum = sum(value$values)
    ),
    execution = list(
      lowering = plan$lowering,
      kernel_version = plan$kernel_version,
      scientific_plan_id = plan$scientific_plan_id,
      completion_status = value$receipt$completion_status
    ),
    session = utils::sessionInfo()
  )
  saveRDS(result, result_path, compress = FALSE)
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (!length(arguments) %in% c(3L, 7L)) {
    stop(paste(
      "usage: crossnobis-scale-worker.R",
      "<repo-root> <result-rds> <ready-file>"
    ))
  }
  dims <- if (length(arguments) == 7L) {
    as.integer(arguments[4:6])
  } else {
    c(52L, 42L, 25L)
  }
  slices <- if (length(arguments) == 7L) {
    as.integer(arguments[[7L]])
  } else {
    24L
  }
  run_crossnobis_scale_gate(
    arguments[[1L]], arguments[[2L]], arguments[[3L]],
    full_dim = dims, active_slices = slices
  )
}
