# Faithful conventional summaries and an operational ambiguity criterion.
# Source 00, 01, and 02 before this file.

conventional_summaries_from_blocks <- function(blocks, contrast) {
  if (!is.list(blocks) || length(blocks) < 2L ||
      !all(vapply(blocks, is.matrix, logical(1)))) {
    stop("blocks must contain at least two effect-by-feature matrices.",
         call. = FALSE)
  }
  dimensions <- vapply(blocks, function(value) paste(dim(value), collapse = "x"),
                       character(1))
  if (length(unique(dimensions)) != 1L ||
      any(!is.finite(unlist(blocks, use.names = FALSE)))) {
    stop("blocks must have aligned finite dimensions.", call. = FALSE)
  }
  effect_names <- rownames(blocks[[1L]])
  if (is.null(effect_names) || any(vapply(blocks, function(value) {
    !identical(rownames(value), effect_names)
  }, logical(1)))) {
    stop("blocks must share named effect rows.", call. = FALSE)
  }
  if (!is.numeric(contrast) || is.null(names(contrast)) ||
      !setequal(names(contrast), effect_names) || any(!is.finite(contrast))) {
    stop("contrast must be a finite named vector aligned to the effect rows.",
         call. = FALSE)
  }
  contrast <- contrast[effect_names]
  patterns <- lapply(blocks, function(value) drop(contrast %*% value))
  mean_pattern <- Reduce(`+`, patterns) / length(patterns)
  pairs <- utils::combn(seq_along(patterns), 2L)
  pair_magnitudes <- vapply(seq_len(ncol(pairs)), function(edge) {
    sum(patterns[[pairs[1L, edge]]] * patterns[[pairs[2L, edge]]])
  }, numeric(1))
  list(
    regional_activation = mean(mean_pattern),
    aggregate_crossvalidated_magnitude = mean(pair_magnitudes),
    partition_patterns = patterns,
    pair_magnitudes = pair_magnitudes,
    definitions = c(
      regional_activation = "mean feature value of the partition-averaged signed contrast",
      aggregate_crossvalidated_magnitude = "uniform mean of all unordered cross-partition Euclidean inner products"
    )
  )
}

matched_conventional_truth <- function(bundle = matched_multiscale_scenarios()) {
  summaries <- do.call(rbind, lapply(names(bundle$scenarios), function(scenario) {
    fixture <- bundle$scenarios[[scenario]]
    baseline <- conventional_summaries_from_blocks(
      list(run1 = fixture$effect_matrix, run2 = fixture$effect_matrix),
      fixture$contrast
    )
    data.frame(
      scenario = scenario,
      regional_activation = baseline$regional_activation,
      aggregate_crossvalidated_magnitude =
        baseline$aggregate_crossvalidated_magnitude,
      stringsAsFactors = FALSE
    )
  }))
  rownames(summaries) <- NULL
  list(summaries = summaries, spectrum = bundle$expected)
}

matched_ambiguity <- function(truth, left, right,
                              baseline = "aggregate_crossvalidated_magnitude",
                              baseline_tolerance = 1e-12,
                              spectrum_min_separation = 0.2) {
  if (!baseline %in% names(truth$summaries)) {
    stop("Unknown conventional baseline `", baseline, "`.", call. = FALSE)
  }
  if (!all(c(left, right) %in% truth$summaries$scenario) ||
      identical(left, right)) {
    stop("left and right must name two distinct scenarios.", call. = FALSE)
  }
  rows <- match(c(left, right), truth$summaries$scenario)
  baseline_values <- truth$summaries[[baseline]][rows]
  left_spectrum <- truth$spectrum[
    truth$spectrum$scenario == left, c("scale", "configuration_share")
  ]
  right_spectrum <- truth$spectrum[
    truth$spectrum$scenario == right, c("scale", "configuration_share")
  ]
  right_spectrum <- right_spectrum[
    match(left_spectrum$scale, right_spectrum$scale), , drop = FALSE
  ]
  spectrum_gap <- max(abs(
    left_spectrum$configuration_share - right_spectrum$configuration_share
  ))
  baseline_gap <- abs(diff(baseline_values))
  data.frame(
    left = left,
    right = right,
    baseline = baseline,
    left_value = baseline_values[[1L]],
    right_value = baseline_values[[2L]],
    baseline_gap = unname(baseline_gap),
    baseline_tolerance = baseline_tolerance,
    ambiguous = unname(baseline_gap) <= baseline_tolerance,
    spectrum_max_gap = spectrum_gap,
    spectrum_min_separation = spectrum_min_separation,
    spectrum_separates = spectrum_gap >= spectrum_min_separation,
    criterion_passes = unname(baseline_gap) <= baseline_tolerance &&
      spectrum_gap >= spectrum_min_separation,
    stringsAsFactors = FALSE
  )
}

matched_conventional_fit <- function(simulation, noise_regime,
                                     trials_per_condition, snr, scenario) {
  fitted <- matched_observation_fit(
    simulation, noise_regime, trials_per_condition, snr, scenario
  )
  features <- seq_len(simulation$bundle$metadata$n_features)
  blocks <- lapply(fitted$fit$relation$partitions, function(partition) {
    relation_block(fitted$fit, partition, features)
  })
  names(blocks) <- fitted$fit$relation$partitions
  conventional <- conventional_summaries_from_blocks(
    blocks, simulation$bundle$scenarios[[scenario]]$contrast
  )
  global_frame <- compile_frame(
    whole_brain(normalization = "none"), simulation$bundle$domain
  )
  global_plan <- plan_geometry(
    fitted$fit$relation, global_frame,
    cross_partitions(
      fitted$fit$relation, independence = "independent",
      generalizes_over = "run"
    )
  )
  global <- contrast_energy(
    global_plan, simulation$bundle$scenarios[[scenario]]$contrast
  )
  list(
    cell_id = fitted$cell_id,
    conventional = conventional,
    crossform_global = data.frame(
      regional_activation = unname(global$signed),
      aggregate_crossvalidated_magnitude = unname(global$total)
    ),
    spectrum = fitted$spectrum,
    truth = fitted$truth,
    scale_truth = fitted$scale_truth
  )
}
