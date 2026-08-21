#!/usr/bin/env Rscript
# Canonical six-panel matched-interpretability figure.
#
# The script sources 00--03 when run directly. When sourced by tests, source
# those dependencies first and call the functions below.

.canonical_scenario_labels <- c(
  broad_coherent = "Broad coherent",
  mixed_broad_fine = "50/50 mixed",
  fine_configuration = "Fine configuration"
)

.canonical_scenario_colors <- c(
  broad_coherent = "#0072B2",
  mixed_broad_fine = "#CC79A7",
  fine_configuration = "#D55E00"
)

.canonical_summary <- function(values) {
  finite <- values[is.finite(values)]
  if (!length(finite)) {
    return(c(mean = NA_real_, sd = NA_real_, lower = NA_real_,
             upper = NA_real_, n_valid = 0))
  }
  c(
    mean = mean(finite),
    sd = if (length(finite) > 1L) stats::sd(finite) else 0,
    lower = unname(stats::quantile(finite, 0.025, names = FALSE)),
    upper = unname(stats::quantile(finite, 0.975, names = FALSE)),
    n_valid = length(finite)
  )
}

canonical_figure_data <- function(seeds = 4301:4324,
                                  noise_regime = "gaussian",
                                  trials_per_condition = 24L,
                                  snr = 0.8) {
  if (!is.numeric(seeds) || length(seeds) < 12L || any(!is.finite(seeds)) ||
      any(seeds != as.integer(seeds)) || anyDuplicated(seeds)) {
    stop("seeds must contain at least 12 unique finite whole numbers.",
         call. = FALSE)
  }
  seeds <- as.integer(seeds)
  estimates <- list()
  truth <- NULL
  patterns <- NULL

  for (replication in seq_along(seeds)) {
    simulation <- matched_paired_observations(seed = seeds[[replication]])
    if (is.null(truth)) {
      truth <- simulation$scale_truth[simulation$scale_truth$snr == snr, ]
      global_truth <- simulation$truth[simulation$truth$snr == snr, ]
      pattern_rows <- list()
      for (scenario in names(simulation$bundle$scenarios)) {
        key <- .paired_cell_key(
          noise_regime, trials_per_condition, snr, scenario
        )
        cell <- simulation$cells[[key]]
        pattern <- drop(
          simulation$bundle$scenarios[[scenario]]$contrast %*%
            cell$effect_matrix
        )
        pattern_rows[[length(pattern_rows) + 1L]] <- data.frame(
          scenario = scenario,
          metric = "ground_truth_pattern",
          scale = seq_along(pattern),
          truth = pattern,
          stringsAsFactors = FALSE
        )
      }
      patterns <- do.call(rbind, pattern_rows)
    }

    for (scenario in names(simulation$bundle$scenarios)) {
      fitted <- matched_conventional_fit(
        simulation, noise_regime, trials_per_condition, snr, scenario
      )
      estimates[[length(estimates) + 1L]] <- data.frame(
        seed = seeds[[replication]], scenario = scenario,
        metric = c("regional_activation",
                   "aggregate_crossvalidated_magnitude"),
        scale = NA_real_,
        estimate = c(
          fitted$conventional$regional_activation,
          fitted$conventional$aggregate_crossvalidated_magnitude
        ),
        stringsAsFactors = FALSE
      )
      spectrum <- fitted$spectrum
      for (metric in c("total", "coherent", "configuration",
                       "coherence_fraction")) {
        estimates[[length(estimates) + 1L]] <- data.frame(
          seed = seeds[[replication]], scenario = scenario,
          metric = metric,
          scale = spectrum$scale,
          estimate = spectrum[[metric]],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  estimates <- do.call(rbind, estimates)
  rownames(estimates) <- NULL

  truth_rows <- list()
  for (scenario in unique(estimates$scenario)) {
    global <- global_truth[global_truth$scenario == scenario, ]
    fixture <- matched_multiscale_scenarios()$scenarios[[scenario]]
    signal_scale <- global$signal_scale[[1L]]
    truth_rows[[length(truth_rows) + 1L]] <- data.frame(
      scenario = scenario,
      metric = c("regional_activation",
                 "aggregate_crossvalidated_magnitude"),
      scale = NA_real_,
      truth = c(
        mean(signal_scale * fixture$effect_pattern),
        global$total[[1L]]
      ),
      stringsAsFactors = FALSE
    )
    scale_truth <- truth[truth$scenario == scenario, ]
    for (metric in c("total", "coherent", "configuration")) {
      truth_rows[[length(truth_rows) + 1L]] <- data.frame(
        scenario = scenario, metric = metric, scale = scale_truth$scale,
        truth = scale_truth[[metric]], stringsAsFactors = FALSE
      )
    }
    truth_rows[[length(truth_rows) + 1L]] <- data.frame(
      scenario = scenario, metric = "coherence_fraction",
      scale = scale_truth$scale, truth = scale_truth$coherent_share,
      stringsAsFactors = FALSE
    )
  }
  truth_rows <- do.call(rbind, truth_rows)

  groups <- split(
    seq_len(nrow(estimates)),
    paste(estimates$scenario, estimates$metric,
          ifelse(is.na(estimates$scale), "NA", estimates$scale), sep = "::")
  )
  summary_rows <- lapply(groups, function(rows) {
    first <- estimates[rows[[1L]], ]
    stats <- .canonical_summary(estimates$estimate[rows])
    expected <- truth_rows$truth[
      truth_rows$scenario == first$scenario &
        truth_rows$metric == first$metric &
        ((is.na(truth_rows$scale) & is.na(first$scale)) |
           truth_rows$scale == first$scale)
    ][[1L]]
    data.frame(
      scenario = first$scenario,
      scenario_label = unname(.canonical_scenario_labels[[first$scenario]]),
      metric = first$metric,
      scale = first$scale,
      truth = expected,
      estimate_mean = unname(stats[["mean"]]),
      estimate_sd = unname(stats[["sd"]]),
      interval_lower = unname(stats[["lower"]]),
      interval_upper = unname(stats[["upper"]]),
      n_valid = as.integer(stats[["n_valid"]]),
      n_replications = length(seeds),
      noise_regime = noise_regime,
      trials_per_condition = as.integer(trials_per_condition),
      snr = snr,
      first_seed = min(seeds),
      last_seed = max(seeds),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL
  pattern_summary <- data.frame(
    scenario = patterns$scenario,
    scenario_label = unname(.canonical_scenario_labels[patterns$scenario]),
    metric = patterns$metric,
    scale = patterns$scale,
    truth = patterns$truth,
    estimate_mean = NA_real_, estimate_sd = NA_real_,
    interval_lower = NA_real_, interval_upper = NA_real_,
    n_valid = NA_integer_, n_replications = length(seeds),
    noise_regime = noise_regime,
    trials_per_condition = as.integer(trials_per_condition), snr = snr,
    first_seed = min(seeds), last_seed = max(seeds),
    stringsAsFactors = FALSE
  )
  output <- rbind(pattern_summary, summary)
  output$scenario <- factor(
    output$scenario, levels = names(.canonical_scenario_labels)
  )
  output <- output[order(output$metric, output$scenario, output$scale), ]
  output$scenario <- as.character(output$scenario)
  rownames(output) <- NULL
  output
}

.canonical_error_plot <- function(data, metric, title, ylab, ylim = NULL) {
  value <- data[data$metric == metric & is.na(data$scale), ]
  value <- value[match(names(.canonical_scenario_labels), value$scenario), ]
  if (is.null(ylim)) {
    limits <- range(c(value$truth, value$interval_lower,
                      value$interval_upper), finite = TRUE)
    padding <- max(diff(limits) * 0.12, 0.05)
    ylim <- limits + c(-padding, padding)
  }
  x <- seq_len(nrow(value))
  plot(x, value$estimate_mean, type = "n", xaxt = "n", xlab = "Scenario",
       ylab = ylab, main = title, ylim = ylim)
  axis(1, at = x, labels = c("Broad", "Mixed", "Fine"))
  abline(h = 0, col = "grey85", lty = 3)
  segments(x, value$interval_lower, x, value$interval_upper,
           col = .canonical_scenario_colors[value$scenario], lwd = 2)
  points(x, value$estimate_mean, pch = 16,
         col = .canonical_scenario_colors[value$scenario], cex = 1.1)
  points(x, value$truth, pch = 4, lwd = 2, cex = 1.2)
}

render_canonical_figure <- function(data, path, width = 1800L,
                                    height = 1200L, res = 150L) {
  required <- c(
    "ground_truth_pattern", "regional_activation",
    "aggregate_crossvalidated_magnitude", "coherent", "configuration",
    "coherence_fraction"
  )
  if (!all(required %in% data$metric)) {
    stop("Figure data is missing a required panel metric.", call. = FALSE)
  }
  grDevices::png(path, width = width, height = height, res = res,
                 type = "cairo")
  on.exit(grDevices::dev.off(), add = TRUE)
  par(mfrow = c(2, 3), mar = c(4.2, 4.4, 3.0, 1.0),
      oma = c(1.2, 0.5, 2.2, 0.5), las = 1)

  pattern <- data[data$metric == "ground_truth_pattern", ]
  pattern_ylim <- range(pattern$truth)
  plot(range(pattern$scale), pattern_ylim, type = "n",
       xlab = "Feature position", ylab = "Planted contrast",
       main = "A  Ground truth")
  for (scenario in names(.canonical_scenario_labels)) {
    part <- pattern[pattern$scenario == scenario, ]
    lines(part$scale, part$truth,
          col = .canonical_scenario_colors[[scenario]], lwd = 2)
  }
  legend("topright", legend = unname(.canonical_scenario_labels),
         col = .canonical_scenario_colors, lwd = 2, bty = "n", cex = 0.78)

  .canonical_error_plot(
    data, "regional_activation", "B  Regional activation",
    "Mean signed contrast"
  )
  .canonical_error_plot(
    data, "aggregate_crossvalidated_magnitude",
    "C  Aggregate multivariate magnitude", "Crossvalidated magnitude"
  )

  widest <- max(data$scale[
    data$metric == "coherent" & is.finite(data$scale)
  ])
  components <- data[data$metric %in% c("coherent", "configuration") &
                       data$scale == widest, ]
  component_ylim <- range(c(0, components$truth, components$interval_lower,
                            components$interval_upper), finite = TRUE)
  component_ylim[[2L]] <- component_ylim[[2L]] * 1.08
  coherent <- data[data$metric == "coherent" & data$scale == widest, ]
  coherent$scale <- NA_real_
  temporary <- rbind(data, transform(coherent, metric = "coherent_wide"))
  .canonical_error_plot(
    temporary, "coherent_wide", "D  Coherent magnitude",
    paste0("Magnitude at radius ", widest), ylim = component_ylim
  )
  configuration <- data[
    data$metric == "configuration" & data$scale == widest, ]
  configuration$scale <- NA_real_
  temporary <- rbind(data, transform(configuration,
                                     metric = "configuration_wide"))
  .canonical_error_plot(
    temporary, "configuration_wide", "E  Configuration magnitude",
    paste0("Magnitude at radius ", widest), ylim = component_ylim
  )

  spectrum <- data[data$metric == "coherence_fraction", ]
  plot(range(spectrum$scale), c(0, 1), type = "n",
       xlab = "Searchlight radius (features)",
       ylab = "Coherent share", main = "F  Coherence spectrum")
  for (scenario in names(.canonical_scenario_labels)) {
    part <- spectrum[spectrum$scenario == scenario, ]
    part <- part[order(part$scale), ]
    polygon(c(part$scale, rev(part$scale)),
            c(pmax(0, part$interval_lower), rev(pmin(1, part$interval_upper))),
            col = grDevices::adjustcolor(.canonical_scenario_colors[[scenario]],
                                        alpha.f = 0.16),
            border = NA)
    lines(part$scale, part$estimate_mean,
          col = .canonical_scenario_colors[[scenario]], lwd = 2)
    points(part$scale, part$truth, pch = 4,
           col = .canonical_scenario_colors[[scenario]], lwd = 1.5)
  }
  mtext(sprintf(
    "Crosses: truth; points/lines: mean; whiskers/bands: empirical 95%% interval; %d paired replications",
    unique(data$n_replications)[[1L]]
  ), outer = TRUE, side = 3, line = 0.5, cex = 0.82)

  invisible(list(
    panels = c("ground_truth", "regional_activation", "aggregate_magnitude",
               "coherent_magnitude", "configuration_magnitude",
               "coherence_spectrum"),
    dimensions = c(width = width, height = height, resolution = res),
    component_ylim = component_ylim,
    widest_scale = widest,
    scenarios = names(.canonical_scenario_labels)
  ))
}

if (sys.nframe() == 0L) {
  script <- sub("^--file=", "",
                grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  directory <- normalizePath(dirname(script))
  for (file in c(
    "00-mixture-generator.R", "01-multiscale-scenarios.R",
    "02-paired-observations.R", "03-conventional-baselines.R"
  )) source(file.path(directory, file))
  suppressMessages(pkgload::load_all(
    normalizePath(file.path(directory, "..", "..")), quiet = TRUE,
    export_all = FALSE
  ))
  output_directory <- normalizePath(file.path(
    directory, "..", "..", "inst", "extdata", "certification"
  ), mustWork = FALSE)
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  figure_data <- canonical_figure_data()
  utils::write.csv(
    figure_data,
    file.path(output_directory, "matched-interpretability-figure-data.csv"),
    row.names = FALSE
  )
  render_canonical_figure(
    figure_data,
    file.path(output_directory, "matched-interpretability-figure.png")
  )
  message("Wrote canonical figure data and PNG to ", output_directory)
}
