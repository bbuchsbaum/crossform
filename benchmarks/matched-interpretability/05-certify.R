#!/usr/bin/env Rscript
# Full matched-interpretability certification.
#
# The exhaustive grid uses direct balanced-design OLS and independent frame
# formulas. The ordinary test suite separately runs public-path smoke cells.

.certification_patterns <- function(cell) {
  lapply(cell$responses, function(response) {
    first <- cell$design[, "condition_a"] == 1
    second <- cell$design[, "condition_b"] == 1
    colMeans(response[first, , drop = FALSE]) -
      colMeans(response[second, , drop = FALSE])
  })
}

.certification_cell <- function(simulation, manifest_row, seed) {
  cell <- simulation$cells[[manifest_row$cell_id]]
  patterns <- .certification_patterns(cell)
  pairs <- utils::combn(seq_along(patterns), 2L)
  mean_pattern <- Reduce(`+`, patterns) / length(patterns)
  pair_total <- vapply(seq_len(ncol(pairs)), function(edge) {
    sum(patterns[[pairs[1L, edge]]] * patterns[[pairs[2L, edge]]])
  }, numeric(1))
  activation <- mean(mean_pattern)
  aggregate <- mean(pair_total)
  scale_truth <- simulation$scale_truth[
    simulation$scale_truth$scenario == manifest_row$scenario &
      simulation$scale_truth$snr == manifest_row$snr,
    , drop = FALSE
  ]
  global_truth <- simulation$truth[
    simulation$truth$scenario == manifest_row$scenario &
      simulation$truth$snr == manifest_row$snr,
    , drop = FALSE
  ]
  truth_pattern <- drop(
    simulation$bundle$scenarios[[manifest_row$scenario]]$contrast %*%
      cell$effect_matrix
  )

  rows <- lapply(seq_along(simulation$bundle$frames), function(scale) {
    weights <- as.matrix(simulation$bundle$frames[[scale]]$weights)
    row_mass <- rowSums(weights)
    totals <- vapply(seq_len(ncol(pairs)), function(edge) {
      left <- patterns[[pairs[1L, edge]]]
      right <- patterns[[pairs[2L, edge]]]
      sum(weights %*% (left * right))
    }, numeric(1))
    coherents <- vapply(seq_len(ncol(pairs)), function(edge) {
      left <- drop(weights %*% patterns[[pairs[1L, edge]]])
      right <- drop(weights %*% patterns[[pairs[2L, edge]]])
      sum(left * right / row_mass)
    }, numeric(1))
    total <- manifest_row$partitions * 0 # typed numeric zero
    total <- simulation$bundle$metadata$alpha[[scale]] * mean(totals)
    coherent <- simulation$bundle$metadata$alpha[[scale]] * mean(coherents)
    configuration <- total - coherent
    expected <- scale_truth[scale_truth$scale ==
                              simulation$bundle$metadata$radii[[scale]], ]
    data.frame(
      seed = seed,
      cell_id = manifest_row$cell_id,
      scenario = manifest_row$scenario,
      noise_regime = manifest_row$noise_regime,
      trials_per_condition = manifest_row$trials_per_condition,
      snr = manifest_row$snr,
      regime_label = manifest_row$regime_label,
      scale = simulation$bundle$metadata$radii[[scale]],
      activation = activation,
      aggregate = aggregate,
      total = total,
      coherent = coherent,
      configuration = configuration,
      # Crossvalidated components are unbiased signed estimates. In finite
      # samples either component may cross zero, so truncating before taking
      # the ratio would selectively discard the hardest replications.
      coherent_share = if (is.finite(total) && abs(total) > .Machine$double.eps)
        coherent / total else NA_real_,
      truth_activation = mean(truth_pattern),
      truth_aggregate = global_truth$total[[1L]],
      truth_total = expected$total[[1L]],
      truth_coherent = expected$coherent[[1L]],
      truth_configuration = expected$configuration[[1L]],
      truth_coherent_share = expected$coherent_share[[1L]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

matched_certification_raw <- function(seeds = 5001:5048) {
  if (!is.numeric(seeds) || !length(seeds) || any(!is.finite(seeds)) ||
      any(seeds != as.integer(seeds)) || anyDuplicated(seeds)) {
    stop("seeds must be unique finite whole numbers.", call. = FALSE)
  }
  output <- list()
  for (seed in as.integer(seeds)) {
    simulation <- matched_paired_observations(seed = seed)
    for (row in seq_len(nrow(simulation$manifest))) {
      output[[length(output) + 1L]] <- .certification_cell(
        simulation, simulation$manifest[row, , drop = FALSE], seed
      )
    }
  }
  value <- do.call(rbind, output)
  rownames(value) <- NULL
  value
}

.certification_long_summary <- function(raw) {
  metric_map <- c(
    activation = "truth_activation",
    aggregate = "truth_aggregate",
    total = "truth_total",
    coherent = "truth_coherent",
    configuration = "truth_configuration",
    coherent_share = "truth_coherent_share"
  )
  rows <- list()
  group <- interaction(
    raw$scenario, raw$noise_regime, raw$trials_per_condition, raw$snr,
    raw$scale, drop = TRUE, lex.order = TRUE
  )
  groups <- split(seq_len(nrow(raw)), group)
  for (indices in groups) {
    first <- raw[indices[[1L]], ]
    for (metric in names(metric_map)) {
      values <- raw[[metric]][indices]
      finite <- values[is.finite(values)]
      rows[[length(rows) + 1L]] <- data.frame(
        scenario = first$scenario,
        noise_regime = first$noise_regime,
        trials_per_condition = first$trials_per_condition,
        snr = first$snr,
        regime_label = first$regime_label,
        scale = first$scale,
        metric = metric,
        truth = first[[metric_map[[metric]]]],
        estimate_mean = if (length(finite)) mean(finite) else NA_real_,
        estimate_sd = if (length(finite) > 1L) stats::sd(finite) else NA_real_,
        interval_lower = if (length(finite)) unname(stats::quantile(
          finite, 0.025, names = FALSE
        )) else NA_real_,
        interval_upper = if (length(finite)) unname(stats::quantile(
          finite, 0.975, names = FALSE
        )) else NA_real_,
        n_valid = length(finite),
        n_replications = length(unique(raw$seed[indices])),
        stringsAsFactors = FALSE
      )
    }
  }
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

.certification_metric_row <- function(metric, group, observed, threshold,
                                      comparison, passes, boundary) {
  data.frame(
    metric = metric, group = group, observed = observed,
    threshold = threshold, comparison = comparison, passes = passes,
    boundary = boundary, stringsAsFactors = FALSE
  )
}

.certification_pair_difference <- function(data, left, right, value) {
  left_rows <- data[data$scenario == left, ]
  right_rows <- data[data$scenario == right, ]
  key <- function(x) paste(x$seed, x$noise_regime,
                           x$trials_per_condition, x$snr, x$scale, sep = "::")
  right_rows <- right_rows[match(key(left_rows), key(right_rows)), ]
  left_rows[[value]] - right_rows[[value]]
}

matched_certification_metrics <- function(raw) {
  thresholds <- list(
    conservation_abs = 1e-10,
    recomposition_abs = 1e-10,
    null_false_separation_abs = 1e-12,
    recoverable_total_relative_bias = 0.12,
    recoverable_share_mae = 0.12,
    recoverable_ordering_rate = 0.80,
    ambiguity_equivalence_fraction = 0.15,
    spectrum_separation = 0.35,
    activation_negative_control = 0.15
  )
  rows <- list()
  add <- function(...) rows[[length(rows) + 1L]] <<-
    .certification_metric_row(...)

  cell_groups <- split(seq_len(nrow(raw)), paste(raw$seed, raw$cell_id))
  conservation <- vapply(cell_groups, function(indices) {
    abs(sum(raw$total[indices]) - raw$aggregate[indices[[1L]]])
  }, numeric(1))
  recomposition <- abs(raw$coherent + raw$configuration - raw$total)
  add("conservation_max_abs", "all cells", max(conservation),
      thresholds$conservation_abs, "<=",
      max(conservation) <= thresholds$conservation_abs,
      "sum of alpha-weighted scale totals equals global aggregate")
  add("recomposition_max_abs", "all rows", max(recomposition),
      thresholds$recomposition_abs, "<=",
      max(recomposition) <= thresholds$recomposition_abs,
      "coherent plus configuration equals total")

  null <- raw[raw$snr == 0, ]
  null_groups <- split(seq_len(nrow(null)), paste(
    null$seed, null$noise_regime, null$trials_per_condition, null$scale
  ))
  null_range <- max(vapply(null_groups, function(indices) {
    max(vapply(c("activation", "aggregate", "total", "coherent",
                 "configuration"), function(metric) {
      diff(range(null[[metric]][indices]))
    }, numeric(1)))
  }, numeric(1)))
  add("null_false_separation_max_abs", "all paired null scenarios", null_range,
      thresholds$null_false_separation_abs, "<=",
      null_range <= thresholds$null_false_separation_abs,
      "identical null observations cannot separate organization labels")

  recoverable <- raw[
    raw$snr == 0.8 & raw$trials_per_condition == 24L & raw$scale == 4.01,
  ]
  for (noise in unique(recoverable$noise_regime)) {
    for (scenario in unique(recoverable$scenario)) {
      part <- recoverable[
        recoverable$noise_regime == noise & recoverable$scenario == scenario,
      ]
      relative_bias <- abs(mean(part$aggregate) - part$truth_aggregate[[1L]]) /
        part$truth_aggregate[[1L]]
      add("recoverable_total_relative_bias", paste(noise, scenario),
          relative_bias, thresholds$recoverable_total_relative_bias, "<=",
          relative_bias <= thresholds$recoverable_total_relative_bias,
          "n=24 SNR=0.8 crossvalidated aggregate")

      all_scales <- raw[
        raw$snr == 0.8 & raw$trials_per_condition == 24L &
          raw$noise_regime == noise & raw$scenario == scenario &
          raw$scale > 0.01 & is.finite(raw$coherent_share),
      ]
      share_mae <- mean(abs(
        all_scales$coherent_share - all_scales$truth_coherent_share
      ))
      add("recoverable_share_mae", paste(noise, scenario), share_mae,
          thresholds$recoverable_share_mae, "<=",
          share_mae <= thresholds$recoverable_share_mae,
          "nonpoint coherence share at n=24 SNR=0.8")
    }

    part <- recoverable[recoverable$noise_regime == noise, ]
    seed_groups <- split(seq_len(nrow(part)), part$seed)
    ordered <- vapply(seed_groups, function(indices) {
      value <- part[indices, ]
      share <- value$coherent_share[match(
        c("broad_coherent", "mixed_broad_fine", "fine_configuration"),
        value$scenario
      )]
      all(is.finite(share)) && share[[1L]] > share[[2L]] &&
        share[[2L]] > share[[3L]]
    }, logical(1))
    ordering_rate <- mean(ordered)
    add("recoverable_ordering_rate", noise, ordering_rate,
        thresholds$recoverable_ordering_rate, ">=",
        ordering_rate >= thresholds$recoverable_ordering_rate,
        "broad greater than mixed greater than fine coherent share")

    difference <- .certification_pair_difference(
      part, "broad_coherent", "mixed_broad_fine", "aggregate"
    )
    equivalence_bound <- abs(mean(difference)) +
      1.96 * stats::sd(difference) / sqrt(length(difference))
    margin <- thresholds$ambiguity_equivalence_fraction *
      unique(part$truth_aggregate)[[1L]]
    add("aggregate_ambiguity_95_bound", noise, equivalence_bound, margin,
        "<=", equivalence_bound <= margin,
        "paired broad-minus-mixed aggregate mean plus 1.96 Monte Carlo SE")

    share_difference <- .certification_pair_difference(
      part, "broad_coherent", "mixed_broad_fine", "coherent_share"
    )
    separation <- mean(share_difference, na.rm = TRUE)
    add("spectrum_separation", noise, separation,
        thresholds$spectrum_separation, ">=",
        separation >= thresholds$spectrum_separation,
        "wide-scale broad-minus-mixed coherent share")

    activation_difference <- .certification_pair_difference(
      part, "broad_coherent", "mixed_broad_fine", "activation"
    )
    activation_gap <- mean(activation_difference)
    add("activation_negative_control", noise, activation_gap,
        thresholds$activation_negative_control, ">=",
        activation_gap >= thresholds$activation_negative_control,
        "regional activation correctly separates broad from mixed")
  }
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  attr(output, "thresholds") <- thresholds
  output
}

matched_certification_truth <- function() {
  simulation <- matched_paired_observations(seed = 5001L)
  global <- do.call(rbind, lapply(c("total", "coherent", "configuration"),
    function(metric) data.frame(
      scenario = simulation$truth$scenario,
      snr = simulation$truth$snr,
      scale = NA_real_, metric = metric,
      value = simulation$truth[[metric]], stringsAsFactors = FALSE
    )))
  scale <- do.call(rbind, lapply(
    c("total", "coherent", "configuration", "coherent_share",
      "configuration_share"),
    function(metric) data.frame(
      scenario = simulation$scale_truth$scenario,
      snr = simulation$scale_truth$snr,
      scale = simulation$scale_truth$scale,
      metric = metric, value = simulation$scale_truth[[metric]],
      stringsAsFactors = FALSE
    )
  ))
  output <- rbind(global, scale)
  output <- output[order(output$scenario, output$snr, output$scale,
                         output$metric), ]
  rownames(output) <- NULL
  output
}

matched_certification_parameters <- function(seeds, metrics) {
  thresholds <- attr(metrics, "thresholds")
  values <- c(
    schema_version = "matched-interpretability-certification-v1",
    evidence_claim = "CF-H10",
    tier = "full",
    first_seed = min(seeds), last_seed = max(seeds),
    replications = length(seeds),
    scenarios = "broad_coherent,mixed_broad_fine,fine_configuration",
    noise_regimes = "gaussian,heteroskedastic,spatial_correlated",
    trials_per_condition = "6,24",
    snr = "0,0.2,0.8",
    frame_radii = "0.01,1.01,2.01,4.01",
    frame_alpha = "0.25,0.25,0.25,0.25",
    vapply(thresholds, format, character(1), digits = 17,
           scientific = FALSE, trim = TRUE)
  )
  data.frame(key = names(values), value = unname(values),
             stringsAsFactors = FALSE)
}

matched_certification_smoke <- function(seeds = c(6201L, 6202L)) {
  raw <- matched_certification_raw(seeds)
  cell_groups <- split(seq_len(nrow(raw)), paste(raw$seed, raw$cell_id))
  conservation <- max(vapply(cell_groups, function(indices) {
    abs(sum(raw$total[indices]) - raw$aggregate[indices[[1L]]])
  }, numeric(1)))
  recomposition <- max(abs(raw$coherent + raw$configuration - raw$total))
  null <- raw[raw$snr == 0, ]
  null_groups <- split(seq_len(nrow(null)), paste(
    null$seed, null$noise_regime, null$trials_per_condition, null$scale
  ))
  null_separation <- max(vapply(null_groups, function(indices) {
    diff(range(null$total[indices]))
  }, numeric(1)))
  list(
    passes = conservation <= 1e-10 && recomposition <= 1e-10 &&
      null_separation <= 1e-12,
    conservation = conservation,
    recomposition = recomposition,
    null_separation = null_separation,
    cells = length(unique(raw$cell_id)),
    seeds = as.integer(seeds)
  )
}

.certification_checksums <- function(repo) {
  entries <- data.frame(
    role = c(rep("source", 6L), "contract", "evidence_ledger",
             rep("numeric_artifact", 6L), "semantic_plot"),
    path = c(
      paste0("benchmarks/matched-interpretability/0", 0:5, c(
        "-mixture-generator.R", "-multiscale-scenarios.R",
        "-paired-observations.R", "-conventional-baselines.R",
        "-canonical-figure.R", "-certify.R"
      )),
      "design/matched-interpretability-simulation-contract.md",
      "inst/extdata/certification/evidence-status-ledger.csv",
      "inst/extdata/certification/matched-interpretability-parameters.csv",
      "inst/extdata/certification/matched-interpretability-truth.csv",
      "inst/extdata/certification/matched-interpretability-summary.csv",
      "inst/extdata/certification/matched-interpretability-metrics.csv",
      "inst/extdata/certification/matched-interpretability-figure-data.csv",
      "vignettes/matched-interpretability.Rmd",
      "inst/extdata/certification/matched-interpretability-figure.png"
    ),
    stringsAsFactors = FALSE
  )
  paths <- file.path(repo, entries$path)
  if (any(!file.exists(paths))) {
    stop("Certification checksum input missing: ",
         paste(entries$path[!file.exists(paths)], collapse = ", "),
         call. = FALSE)
  }
  semantic <- entries$role == "semantic_plot"
  data.frame(
    schema_version = "matched-interpretability-checksums-v1",
    role = entries$role,
    path = entries$path,
    size_bytes = unname(file.info(paths)$size),
    hash_algorithm = ifelse(semantic, "semantic_render_check", "md5"),
    digest = ifelse(semantic, "", unname(tools::md5sum(paths))),
    stringsAsFactors = FALSE
  )
}

if (sys.nframe() == 0L) {
  script <- sub("^--file=", "",
                grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  directory <- normalizePath(dirname(script))
  repo <- normalizePath(file.path(directory, "..", ".."))
  for (file in c(
    "00-mixture-generator.R", "01-multiscale-scenarios.R",
    "02-paired-observations.R", "03-conventional-baselines.R"
  )) source(file.path(directory, file))
  suppressMessages(pkgload::load_all(repo, quiet = TRUE, export_all = FALSE))
  seeds <- 5001:5048
  raw <- matched_certification_raw(seeds)
  metrics <- matched_certification_metrics(raw)
  summary <- .certification_long_summary(raw)
  truth <- matched_certification_truth()
  parameters <- matched_certification_parameters(seeds, metrics)
  output <- file.path(repo, "inst", "extdata", "certification")
  utils::write.csv(parameters,
    file.path(output, "matched-interpretability-parameters.csv"),
    row.names = FALSE)
  utils::write.csv(truth,
    file.path(output, "matched-interpretability-truth.csv"), row.names = FALSE)
  utils::write.csv(summary,
    file.path(output, "matched-interpretability-summary.csv"), row.names = FALSE)
  utils::write.csv(metrics,
    file.path(output, "matched-interpretability-metrics.csv"), row.names = FALSE)
  checksums <- .certification_checksums(repo)
  utils::write.csv(checksums,
    file.path(output, "matched-interpretability-checksums.csv"),
    row.names = FALSE)
  print(metrics, row.names = FALSE)
  if (!all(metrics$passes)) {
    stop("Matched interpretability certification failed ",
         sum(!metrics$passes), " metric(s).", call. = FALSE)
  }
  message("Matched interpretability certification PASS: ", length(seeds),
          " paired replications, ", nrow(raw), " scale rows.")
}
