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
  noncv_pairing = "biased_self_pairing_single_partition",
  correlation_route = "downstream_from_psd_within_sample_self_form",
  correlation_support_scope = "one_plan_per_support_not_frame_composable",
  comparison_route = "downstream_base_r_no_crossform_export",
  claim_scope = "standard_workflow_distances_comparisons_and_fixed_linear_rsa"
)
utils::write.csv(
  data.frame(key = names(fixture_meta), value = unname(fixture_meta),
             stringsAsFactors = FALSE),
  file.path(results, "fixture-meta.csv"), row.names = FALSE)

## ---- Non-crossvalidated distances (the within-sample arm) ---------------
# rsatoolbox's calc_rdm(method = "euclidean" / "mahalanobis" / "correlation")
# averages every observation of a condition across the whole dataset. With a
# balanced indicator design repeated identically in each run, that grand mean
# is exactly the mean of the per-run OLS blocks, so both sides contract the
# same patterns. crossform expresses "no cross-validation" as a declared
# biased self-pairing rather than by dropping the generalization axis.
grand <- Reduce(`+`, betas) / length(betas)
within <- relation(list(all = grand), domain = domain)
within_over <- pairing("all", "all", self_pairs = "allow_biased",
                       independence = "not_independent")
stopifnot(identical(attr(within_over, "estimate"), "self_product_biased"))

distance_rows <- list()
for (nm in names(frames)) {
  frame <- frames[[nm]]

  euclidean_plan <- plan_geometry(within, at = frame, over = within_over)
  mahalanobis_plan <- plan_geometry(within, at = frame, over = within_over,
                                    metric = metric)
  euclidean_view <- rdm(euclidean_plan)
  mahalanobis_view <- rdm(mahalanobis_plan)
  measurements <- as.character(euclidean_view$index)
  stopifnot(identical(euclidean_view$pairs, pair_frame))

  for (i in seq_along(measurements)) {
    distance_rows[[length(distance_rows) + 1L]] <- data.frame(
      frame = nm, measurement = measurements[i],
      pair = seq_len(nrow(pair_frame)),
      left = pair_frame$left, right = pair_frame$right,
      euclidean = as.numeric(as.matrix(euclidean_view$values)[i, ]),
      mahalanobis = as.numeric(as.matrix(mahalanobis_view$values)[i, ]),
      stringsAsFactors = FALSE
    )
  }
}
crossform_distances <- do.call(rbind, distance_rows)
message("non-crossvalidated distances: ", nrow(crossform_distances),
        " rows over ", length(unique(crossform_distances$measurement)),
        " measurements")

## ---- Correlation distance, one support at a time ------------------------
# Euclidean and Mahalanobis are frame-composable: one frame carrying several
# supports yields all of them in a single pass, because the distance is a
# fixed bilinear query and `normalization = "local"` supplies the channel
# count per support.
#
# Correlation distance is not. Both its centering and its normalizer depend on
# which channels are in the support -- rsatoolbox's `calc_rdm_correlation`
# calls `_parse_input(..., remove_mean = TRUE)` on the *restricted* dataset --
# so a multi-support frame cannot produce it in one pass. That is a concrete
# consequence of its being a nonlinear view rather than a bilinear query, and
# it is exactly the kind of structure the correlation-distance policy is
# about. Each support therefore gets its own domain, its own support-centered
# patterns, and its own plan.
supports <- c(list(whole_brain = seq_len(N_VOXELS)),
              split(seq_len(N_VOXELS), REGION_LABELS))
correlation_rows <- list()
for (nm in names(supports)) {
  support <- supports[[nm]]
  patterns <- grand[, support, drop = FALSE]
  patterns <- patterns - rowMeans(patterns)
  support_domain <- abstract_domain(
    length(support), id = paste0("rsatoolbox-parity-support-", nm)
  )
  support_relation <- relation(list(all = patterns), domain = support_domain)
  support_plan <- plan_geometry(
    support_relation, at = compile_frame(whole_brain(), support_domain),
    over = within_over
  )
  distances <- as.numeric(as.matrix(rdm(support_plan)$values)[1L, ])
  diagonals <- vapply(seq_along(CONDITIONS), function(i) {
    w <- setNames(rep(0, length(CONDITIONS)), CONDITIONS)
    w[i] <- 1
    as.numeric(contrast_energy(support_plan, w)$total)[1L]
  }, numeric(1))
  correlation_rows[[length(correlation_rows) + 1L]] <- data.frame(
    measurement = nm, pair = seq_len(nrow(pair_frame)),
    left = pair_frame$left, right = pair_frame$right,
    correlation = correlation_distance(distances, diagonals,
                                       length(CONDITIONS)),
    stringsAsFactors = FALSE
  )
}
crossform_correlation <- do.call(rbind, correlation_rows)
message("correlation distances: ", nrow(crossform_correlation), " rows over ",
        length(supports), " supports (one plan each; not frame-composable)")

## ---- Downstream comparison statistics -----------------------------------
# Computed from crossform's crossnobis RDM with base R. No crossform export is
# involved: `rdm()$values` is already in row-major upper-triangle order, which
# is what numpy's triu_indices produces, so the two vectors align elementwise.
all_models <- c(models, list(graded = graded_model()))
model_vectors <- lapply(all_models, rdm_pair_vector)

similarity_rows <- list()
for (m in unique(crossform_rdm$measurement)) {
  part <- crossform_rdm[crossform_rdm$measurement == m, ]
  part <- part[order(part$pair), ]
  for (model_name in names(model_vectors)) {
    stats <- similarity_statistics(part$crossform, model_vectors[[model_name]],
                                   length(CONDITIONS))
    similarity_rows[[length(similarity_rows) + 1L]] <- data.frame(
      measurement = m, model = model_name, method = names(stats),
      crossform = unname(stats), stringsAsFactors = FALSE
    )
  }
}
crossform_similarity <- do.call(rbind, similarity_rows)
message("downstream similarity statistics: ", nrow(crossform_similarity),
        " rows (", length(model_vectors), " models x 8 methods x ",
        length(unique(crossform_rdm$measurement)), " measurements)")

## ---- Diagonal noise precision -------------------------------------------
# Same pattern as the full precision: the estimator is the exemplar's, applied
# once in R, and 03-compare.R checks that rsatoolbox's own prec_from_residuals
# reproduces it.
diagonal_covariance <- diag(diag(noise$covariance))
diagonal_precision_matrix <- diag(1 / diag(noise$covariance))

utils::write.csv(crossform_distances,
                 file.path(results, "crossform-distances.csv"),
                 row.names = FALSE)
utils::write.csv(crossform_correlation,
                 file.path(results, "crossform-correlation.csv"),
                 row.names = FALSE)
utils::write.csv(crossform_similarity,
                 file.path(results, "crossform-similarity.csv"),
                 row.names = FALSE)
utils::write.csv(as.data.frame(unname(diagonal_precision_matrix)),
                 file.path(results, "precision-diagonal.csv"),
                 row.names = FALSE)
utils::write.csv(
  data.frame(condition = CONDITIONS, as.data.frame(unname(grand)),
             stringsAsFactors = FALSE),
  file.path(results, "grand-means.csv"), row.names = FALSE)
utils::write.csv(
  data.frame(pair = seq_len(nrow(pair_frame)), left = pair_frame$left,
             right = pair_frame$right,
             category = model_vectors$category,
             animacy = model_vectors$animacy,
             graded = model_vectors$graded,
             stringsAsFactors = FALSE),
  file.path(results, "model-rdms-extended.csv"), row.names = FALSE)

saveRDS(list(fixture = fixture, noise = noise, domain = domain, fit = fit,
             over = over, metric = metric, frames = frames, plans = plans,
             models = models, pairs = pair_frame),
        file.path(results, "fixture.rds"))

message("Wrote ", length(list.files(results)), " files to ", results)
message("crossform RDM rows: ", nrow(crossform_rdm),
        "; RSA rows: ", nrow(crossform_rsa))
