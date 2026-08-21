#!/usr/bin/env Rscript
# 01-fixture.R -- build the shared fixture, fit it in crossform, and export
# everything the Python arm needs.
#
# What leaves this script (all CSV, no reticulate):
#   betas.csv              per-run condition patterns, one row per (run, cond)
#   residuals.csv          per-run OLS residuals, one row per (run, obs)
#   precision.csv          the fixed 40 x 40 noise precision both sides use
#   covariance.csv         its inverse, for rsatoolbox's `noise=` documentation
#   regions.csv            region label per voxel (frame support definition)
#   model-rdms.csv         the two model RDMs, vectorised in pair order
#   pairs.csv              the condition pair order crossform reports
#   crossform-rdm.csv      crossform's fixed-metric crossnobis RDM
#   crossform-rsa.csv      crossform's linear-RSA coefficients
#   fixture-meta.csv       scalars the other scripts must not re-derive
#
# The RDM is read with `rdm()` on a `plan_geometry(metric = noise_precision())`
# plan. `crossnobis()` on the same plan is the named Mahalanobis reading of
# the same compiled estimand; 04-extension.R shows they agree exactly.

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
source(file.path(exemplar_dir, "00-common.R"))
mode <- load_crossform(exemplar_dir)
results <- file.path(exemplar_dir, "results")
dir.create(results, showWarnings = FALSE, recursive = TRUE)
message("crossform loaded from ", mode)

## ---- Fixture ------------------------------------------------------------
fixture <- build_fixture()
noise <- pooled_precision(fixture)
message("fixture: ", fixture$q, " conditions x ", N_RUNS, " runs, ",
        N_VOXELS, " voxels, ", fixture$n_obs, " observations per run; ",
        "residual df = ", noise$residual_df)
message("residual covariance condition number: ",
        format(kappa(noise$covariance, exact = TRUE), digits = 4),
        "  (an identity metric is not equivalent here)")

## ---- Fit in crossform ---------------------------------------------------
domain <- abstract_domain(
  N_VOXELS,
  coordinates = cbind(x = seq_len(N_VOXELS), y = 0, z = 0),
  id = "rsatoolbox-parity-domain",
  coordinate_units = "mm"
)

fit <- lm_relation_fit(
  fixture$responses, fixture$design, fixture$effects,
  effect_names = CONDITIONS, sampling_unit = "trial", domain = domain
)
relation <- fit$relation
over <- cross_partitions(relation, independence = "independent",
                         generalizes_over = "run")
stopifnot(nrow(over) == choose(N_RUNS, 2L), abs(sum(over$weight) - 1) < 1e-12)

# The betas crossform will contract are the relation blocks. Assert they are
# the ordinary per-run OLS coefficients before exporting them, so the Python
# arm is demonstrably reading the same estimates and not a re-derivation.
betas <- lapply(relation$partitions, function(partition) {
  relation_block(fit, partition, seq_len(N_VOXELS))
})
names(betas) <- relation$partitions
manual <- lapply(fixture$responses, function(y) {
  solve(crossprod(fixture$design), crossprod(fixture$design, y))
})
beta_gap <- max(vapply(seq_along(betas), function(i) {
  max(abs(betas[[i]] - manual[[i]]))
}, numeric(1)))
message("relation blocks vs plain per-run OLS: max abs diff = ",
        format(beta_gap, digits = 3))
stopifnot(beta_gap < 1e-12)

## ---- Fixed metric and frames --------------------------------------------
metric <- noise_precision(
  noise$precision, domain, covariance = noise$covariance,
  provenance = list(
    source = "pooled within-run residual covariance",
    estimator = "sum_r E_r' E_r / nu, nu = sum_r (n_r - rank(X_r))"
  )
)

frames <- list(
  whole = compile_frame(whole_brain(normalization = "local"), domain),
  regions = compile_frame(regions(REGION_LABELS, normalization = "local"),
                          domain)
)

plans <- lapply(frames, function(frame) {
  plan_geometry(relation, at = frame, over = over, metric = metric)
})

## ---- crossform RDM and RSA ----------------------------------------------
models <- model_rdms()

rdm_rows <- list()
rsa_rows <- list()
pair_frame <- NULL
for (nm in names(plans)) {
  view <- rdm(plans[[nm]])
  values <- as.matrix(view$values)
  measurements <- as.character(view$index)
  if (is.null(pair_frame)) pair_frame <- view$pairs
  stopifnot(identical(view$pairs, pair_frame))
  for (i in seq_along(measurements)) {
    rdm_rows[[length(rdm_rows) + 1L]] <- data.frame(
      frame = nm, measurement = measurements[i],
      pair = seq_len(ncol(values)),
      left = pair_frame$left, right = pair_frame$right,
      crossform = as.numeric(values[i, ]),
      stringsAsFactors = FALSE
    )
  }
  fitted <- rsa(plans[[nm]], models = models)
  coefficients <- as.matrix(fitted$coefficients)
  for (i in seq_along(measurements)) {
    rsa_rows[[length(rsa_rows) + 1L]] <- data.frame(
      frame = nm, measurement = measurements[i],
      term = colnames(coefficients),
      crossform = as.numeric(coefficients[i, ]),
      stringsAsFactors = FALSE
    )
  }
  # The no-intercept fit is the like-for-like comparator for rsatoolbox's
  # `ModelWeighted` + `fit_regress`, which carries no constant column.
  bare <- rsa(plans[[nm]], models = models, intercept = FALSE)
  bare_coefficients <- as.matrix(bare$coefficients)
  for (i in seq_along(measurements)) {
    rsa_rows[[length(rsa_rows) + 1L]] <- data.frame(
      frame = nm, measurement = measurements[i],
      term = paste0("nointercept:", colnames(bare_coefficients)),
      crossform = as.numeric(bare_coefficients[i, ]),
      stringsAsFactors = FALSE
    )
  }
}
crossform_rdm <- do.call(rbind, rdm_rows)
crossform_rsa <- do.call(rbind, rsa_rows)

## ---- Export -------------------------------------------------------------
beta_long <- do.call(rbind, lapply(names(betas), function(run) {
  m <- betas[[run]]
  data.frame(run = run, condition = rownames(m), as.data.frame(unname(m)),
             stringsAsFactors = FALSE)
}))
utils::write.csv(beta_long, file.path(results, "betas.csv"), row.names = FALSE)

residual_long <- do.call(rbind, lapply(names(noise$residuals), function(run) {
  m <- noise$residuals[[run]]
  data.frame(run = run, observation = seq_len(nrow(m)),
             as.data.frame(unname(m)), stringsAsFactors = FALSE)
}))
utils::write.csv(residual_long, file.path(results, "residuals.csv"),
                 row.names = FALSE)

utils::write.csv(as.data.frame(unname(noise$precision)),
                 file.path(results, "precision.csv"), row.names = FALSE)
utils::write.csv(as.data.frame(unname(noise$covariance)),
                 file.path(results, "covariance.csv"), row.names = FALSE)
utils::write.csv(
  data.frame(voxel = seq_len(N_VOXELS), region = REGION_LABELS,
             stringsAsFactors = FALSE),
  file.path(results, "regions.csv"), row.names = FALSE)
utils::write.csv(
  data.frame(pair = seq_len(nrow(pair_frame)), left = pair_frame$left,
             right = pair_frame$right,
             category = rdm_pair_vector(models$category),
             animacy = rdm_pair_vector(models$animacy),
             stringsAsFactors = FALSE),
  file.path(results, "model-rdms.csv"), row.names = FALSE)
utils::write.csv(crossform_rdm, file.path(results, "crossform-rdm.csv"),
                 row.names = FALSE)
utils::write.csv(crossform_rsa, file.path(results, "crossform-rsa.csv"),
                 row.names = FALSE)
fixture_meta <- c(
  seed = SEED,
  n_runs = N_RUNS,
  n_voxels = N_VOXELS,
  n_conditions = fixture$q,
  n_obs_per_run = fixture$n_obs,
  residual_df = noise$residual_df,
  n_pairs = nrow(pair_frame),
  beta_vs_ols_max_abs_diff = beta_gap,
  pairing_edges = nrow(over),
  partition_weight = unique(over$weight),
  covariance_condition_number = kappa(noise$covariance, exact = TRUE),
  metric_role = "fixed_noise_precision",
  metric_estimator = "inverse_pooled_within_run_residual_covariance",
  metric_normalization = "frame_local_divide_by_support_size",
  effect_centering = "none_pair_differences_are_zero_sum",
  partition_scheme = "uniform_unordered_cross_run_pairs",
  pair_order = "row_major_upper_triangle",
  rsa_objective = "fixed_ols_on_vectorized_rdm",
  claim_scope = "crossnobis_and_fixed_linear_rsa_only"
)
utils::write.csv(
  data.frame(key = names(fixture_meta), value = unname(fixture_meta),
             stringsAsFactors = FALSE),
  file.path(results, "fixture-meta.csv"), row.names = FALSE)

saveRDS(list(fixture = fixture, noise = noise, domain = domain, fit = fit,
             over = over, metric = metric, frames = frames, plans = plans,
             models = models, pairs = pair_frame),
        file.path(results, "fixture.rds"))

message("Wrote ", length(list.files(results)), " files to ", results)
message("crossform RDM rows: ", nrow(crossform_rdm),
        "; RSA rows: ", nrow(crossform_rsa))
