# Paired end-to-end observations for the matched multiscale scenarios.
# Source 00-mixture-generator.R and 01-multiscale-scenarios.R first.

.paired_standard_noise <- function(partitions, observations, features, seed) {
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  stats::setNames(lapply(seq_len(partitions), function(partition) {
    matrix(stats::rnorm(observations * features), observations, features)
  }), paste0("run", seq_len(partitions)))
}

.paired_cell_key <- function(noise_regime, trials_per_condition, snr,
                             scenario) {
  paste(noise_regime, paste0("n", trials_per_condition),
        paste0("snr", format(snr, trim = TRUE, scientific = FALSE)),
        scenario, sep = "::")
}

matched_paired_observations <- function(
    bundle = matched_multiscale_scenarios(),
    partitions = 4L,
    trials_per_condition = c(6L, 24L),
    snr = c(0, 0.2, 0.8),
    seed = 20260823L) {
  if (!inherits(bundle$family, "effect_frame") ||
      !identical(bundle$schema_version, "matched-multiscale-v1")) {
    stop("bundle must come from matched_multiscale_scenarios().",
         call. = FALSE)
  }
  if (!is.numeric(partitions) || length(partitions) != 1L ||
      partitions != as.integer(partitions) || partitions < 2L) {
    stop("partitions must be one whole number of at least 2.", call. = FALSE)
  }
  partitions <- as.integer(partitions)
  if (!is.numeric(trials_per_condition) ||
      length(trials_per_condition) < 2L ||
      any(trials_per_condition != as.integer(trials_per_condition)) ||
      any(trials_per_condition < 3L) ||
      is.unsorted(trials_per_condition, strictly = TRUE)) {
    stop("trials_per_condition must contain at least two strictly increasing whole numbers of at least 3.",
         call. = FALSE)
  }
  trials_per_condition <- as.integer(trials_per_condition)
  if (!is.numeric(snr) || length(snr) < 3L || any(!is.finite(snr)) ||
      any(snr < 0) || is.unsorted(snr, strictly = TRUE) || snr[[1L]] != 0) {
    stop("snr must contain null, low-power, and recoverable nonnegative levels in strictly increasing order, starting at zero.",
         call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("seed must be one finite whole number.", call. = FALSE)
  }
  seed <- as.integer(seed)

  features <- bundle$metadata$n_features
  maximum_trials <- max(trials_per_condition)
  maximum_observations <- 2L * maximum_trials
  base_noise <- .paired_standard_noise(
    partitions, maximum_observations, features, seed
  )

  heteroskedastic_scale <- seq(0.55, 1.45, length.out = features)
  heteroskedastic_scale <- heteroskedastic_scale /
    sqrt(mean(heteroskedastic_scale^2))
  spatial_covariance <- toeplitz(0.6^(seq_len(features) - 1L))
  transforms <- list(
    gaussian = diag(features),
    heteroskedastic = diag(heteroskedastic_scale),
    spatial_correlated = chol(spatial_covariance)
  )
  covariance <- lapply(transforms, function(factor) crossprod(factor))
  stopifnot(all(vapply(covariance, function(value) {
    abs(mean(diag(value)) - 1) < 1e-12
  }, logical(1))))

  designs <- stats::setNames(lapply(trials_per_condition, function(n) {
    value <- rbind(
      condition_a = cbind(condition_a = rep(1, n), condition_b = rep(0, n)),
      condition_b = cbind(condition_a = rep(0, n), condition_b = rep(1, n))
    )
    matrix(value, nrow = 2L * n, ncol = 2L,
           dimnames = list(NULL, c("condition_a", "condition_b")))
  }), paste0("n", trials_per_condition))

  snr_labels <- c(
    "null", rep("low_power", length(snr) - 2L), "recoverable"
  )
  cells <- list()
  manifest_rows <- list()
  truth_rows <- list()
  scale_truth_rows <- list()

  for (scenario in names(bundle$scenarios)) {
    base_fixture <- bundle$scenarios[[scenario]]
    base_rms <- sqrt(base_fixture$truth$total / features)
    for (signal_to_noise in snr) {
      signal_scale <- if (signal_to_noise == 0) 0 else signal_to_noise / base_rms
      effect_matrix <- signal_scale * base_fixture$effect_matrix
      global_truth <- signal_scale^2 * unlist(base_fixture$truth[c(
        "total", "coherent", "configuration"
      )])
      truth_rows[[length(truth_rows) + 1L]] <- data.frame(
        scenario = scenario,
        snr = signal_to_noise,
        regime_label = snr_labels[[match(signal_to_noise, snr)]],
        signal_scale = signal_scale,
        total = unname(global_truth[["total"]]),
        coherent = unname(global_truth[["coherent"]]),
        configuration = unname(global_truth[["configuration"]]),
        stringsAsFactors = FALSE
      )
      base_scale_truth <- bundle$expected[
        bundle$expected$scenario == scenario,
        c("scenario", "family", "scale", "alpha", "total", "coherent",
          "configuration", "coherent_share", "configuration_share")
      ]
      base_scale_truth$snr <- signal_to_noise
      base_scale_truth$total <- signal_scale^2 * base_scale_truth$total
      base_scale_truth$coherent <- signal_scale^2 * base_scale_truth$coherent
      base_scale_truth$configuration <-
        signal_scale^2 * base_scale_truth$configuration
      if (signal_to_noise == 0) {
        base_scale_truth$coherent_share <- NA_real_
        base_scale_truth$configuration_share <- NA_real_
      }
      scale_truth_rows[[length(scale_truth_rows) + 1L]] <- base_scale_truth

      for (noise_regime in names(transforms)) {
        transformed <- lapply(base_noise, `%*%`, transforms[[noise_regime]])
        for (n in trials_per_condition) {
          rows <- c(seq_len(n), maximum_trials + seq_len(n))
          design <- designs[[paste0("n", n)]]
          responses <- lapply(transformed, function(noise) {
            design %*% effect_matrix + noise[rows, , drop = FALSE]
          })
          key <- .paired_cell_key(noise_regime, n, signal_to_noise, scenario)
          cells[[key]] <- list(
            responses = responses,
            design = design,
            effect_matrix = effect_matrix
          )
          manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
            cell_id = key,
            scenario = scenario,
            noise_regime = noise_regime,
            trials_per_condition = n,
            snr = signal_to_noise,
            regime_label = snr_labels[[match(signal_to_noise, snr)]],
            partitions = partitions,
            n_features = features,
            seed = seed,
            design_id = paste0("balanced-two-condition:n", n),
            base_noise_id = paste0(
              "paired-standard-normal-v1:seed", seed, ":runs", partitions,
              ":nmax", maximum_trials, ":p", features
            ),
            truth_id = paste0(
              scenario, ":snr",
              format(signal_to_noise, trim = TRUE, scientific = FALSE)
            ),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  list(
    schema_version = "matched-paired-observations-v1",
    bundle = bundle,
    cells = cells,
    truth = do.call(rbind, truth_rows),
    scale_truth = do.call(rbind, scale_truth_rows),
    manifest = do.call(rbind, manifest_rows),
    source = list(
      base_standard_noise = base_noise,
      designs = designs,
      noise_transforms = transforms,
      noise_covariance = covariance
    ),
    metadata = list(
      seed = seed,
      partitions = partitions,
      trials_per_condition = trials_per_condition,
      snr = snr,
      snr_definition = "contrast_pattern_RMS / mean_feature_noise_SD",
      sample_size_relation = "smaller cells are condition-stratified prefixes of the maximum cell",
      pairing = "same design and base standard-normal draw across organizations, SNRs, sample sizes, and noise transforms",
      truth_storage = "truth and scale_truth are separate from cells and estimator output"
    )
  )
}

matched_observation_fit <- function(simulation, noise_regime,
                                    trials_per_condition, snr, scenario) {
  key <- .paired_cell_key(
    noise_regime, trials_per_condition, snr, scenario
  )
  if (!key %in% names(simulation$cells)) {
    stop("No simulation cell named `", key, "`.", call. = FALSE)
  }
  cell <- simulation$cells[[key]]
  effects <- diag(2L)
  dimnames(effects) <- list(colnames(cell$design), colnames(cell$design))
  fit <- lm_relation_fit(
    cell$responses, cell$design, effects,
    effect_names = colnames(cell$design), sampling_unit = "trial",
    domain = simulation$bundle$domain
  )
  plan <- plan_geometry(
    fit$relation, simulation$bundle$family,
    cross_partitions(
      fit$relation, independence = "independent", generalizes_over = "run"
    )
  )
  spectrum <- coherence_spectrum(
    plan, simulation$bundle$scenarios[[scenario]]$contrast
  )
  list(
    cell_id = key,
    fit = fit,
    plan = plan,
    spectrum = as.data.frame(spectrum),
    truth = simulation$truth[
      simulation$truth$scenario == scenario & simulation$truth$snr == snr,
      , drop = FALSE
    ],
    scale_truth = simulation$scale_truth[
      simulation$scale_truth$scenario == scenario &
        simulation$scale_truth$snr == snr,
      , drop = FALSE
    ]
  )
}
