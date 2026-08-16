#!/usr/bin/env Rscript
# 02-crossform-searchlight.R -- the crossform arm.
#
# Consumes data/prepared-smoke.rds and produces, per VT searchlight centre:
#   * rdm_cv     the 28-element cross-run RDM (see below)
#   * rsa_score  the shared Spearman score of that RDM against the model RDM
#   * rsa_ols    crossform's own native rsa() coefficient (OLS in RDM space)
#
# WHAT THE RDM ACTUALLY IS, stated precisely because the comparison depends
# on it. `cross_partitions()` declares the 66 unordered run pairs; geometry
# forms cross-run feature atoms; `rdm()` reads the fixed difference query. The
# result is the CROSS-VALIDATED SQUARED EUCLIDEAN DISTANCE under the identity
# metric with local (mean-over-voxel) frame normalisation:
#
#   d2(i,j) = mean over run pairs r<s of
#               (1/P) * sum_v ( b_i^r[v] - b_j^r[v] ) * ( b_i^s[v] - b_j^s[v] )
#
# It is unbiased around zero, so individual entries may be negative; the
# package deliberately preserves that rather than truncating. This script
# asserts that identity against a direct loop at a sample of centres, so the
# claim in the report is checked, not assumed.
#
# NOT correlation distance. Correlation distance needs each pattern rescaled
# to unit norm inside each sphere, which is nonlinear in the patterns and
# therefore outside crossform's bilinear core. That divergence from the
# README estimand is real and is carried into 04-compare.R.

suppressMessages(library(neuroim2))

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
source(file.path(exemplar_dir, "models.R"))
crossform_source <- load_crossform(exemplar_dir)
paths <- exemplar_paths(exemplar_dir)
dir.create(paths$results, showWarnings = FALSE, recursive = TRUE)

prep <- readRDS(paths$prepared)
stopifnot(identical(prep$conditions, CONDITIONS))

## ---- Domain and frame ---------------------------------------------------
vt_vol <- neuroim2::read_vol(prep$mask_file)
vt_mask <- neuroim2::LogicalNeuroVol(as.array(vt_vol) > 0,
                                     neuroim2::space(vt_vol))
domain <- neuroim2_volume_domain(vt_mask, id = "haxby-subj1-vt")
stopifnot(identical(domain$feature_ids, prep$vt_index))

frame <- neuroim2_searchlights(vt_mask, radius = RADIUS_MM, domain = domain,
                               normalization = "local")
n_centers <- nrow(frame$weights)
message("crossform frame: ", n_centers, " searchlight rows over ",
        length(domain$feature_ids), " VT voxels at r = ", RADIUS_MM, " mm")

## ---- Relation, pairing, plan -------------------------------------------
effects <- effect_space(CONDITIONS, basis_id = "haxby-condition-means:v1",
                        units = "within-run-z")
rel <- relation(prep$runs, effects = effects, domain = domain)
over <- cross_partitions(rel)              # off-diagonal run pairs only
n_run <- length(prep$runs)
stopifnot(length(over$left) == choose(n_run, 2L),
          !any(over$left == over$right))   # no run paired with itself
message("cross_partitions: ", length(over$left), " unordered run pairs, ",
        "no within-run pair")

plan <- plan_geometry(rel, at = frame, over = over)

## ---- The RDM view -------------------------------------------------------
t0 <- Sys.time()
D <- rdm(plan)
elapsed_rdm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
message("rdm() over ", n_centers, " centres in ", round(elapsed_rdm, 2), " s")

# The pair order crossform reports must be the canonical combn order the
# shared model vector uses. Assert rather than assume.
pair_idx <- rdm_pair_index()
stopifnot(identical(D$pairs$left, pair_idx$left),
          identical(D$pairs$right, pair_idx$right))

rdm_cv <- as.matrix(D$values)              # n_centers x 28
centers <- as.integer(D$index)             # full-volume indices of the centres
stopifnot(identical(centers, prep$vt_index))

## ---- Verify the estimand against a direct loop --------------------------
# Independent, deliberately naive recomputation at a sample of centres. This
# checks that the compiled sparse contraction really is the cross-validated
# squared distance defined in the header.
verify_centers <- unique(round(seq(1, n_centers, length.out = 8)))
W <- as.matrix(frame$weights)
nbhd <- neuroim2::searchlight_indices(vt_mask, radius = RADIUS_MM, nonzero = TRUE)
runs <- prep$runs
R <- length(runs)

direct_rdm <- function(center) {
  members <- match(nbhd[[center]], domain$feature_ids)
  w <- W[center, members]
  vapply(seq_len(nrow(pair_idx)), function(k) {
    i <- pair_idx$i[k]; j <- pair_idx$j[k]
    acc <- 0; n <- 0L
    for (r in seq_len(R - 1L)) for (s in seq(r + 1L, R)) {
      dr <- runs[[r]][i, members] - runs[[r]][j, members]
      ds <- runs[[s]][i, members] - runs[[s]][j, members]
      acc <- acc + sum(w * dr * ds); n <- n + 1L
    }
    acc / n
  }, numeric(1))
}

verify_max_abs <- max(vapply(verify_centers, function(cc) {
  max(abs(direct_rdm(cc) - rdm_cv[cc, ]))
}, numeric(1)))
message("self-verification vs direct loop at ", length(verify_centers),
        " centres: max abs difference = ", format(verify_max_abs, digits = 3))
if (verify_max_abs > TOLERANCE) {
  stop("crossform rdm() does not match the direct cross-validated distance.")
}

## ---- Shared second-order score -----------------------------------------
model_vec <- model_vector(MODEL_RDMS$animacy)
rsa_score <- apply(rdm_cv, 1L, spearman_rsa, model_vec = model_vec)

## ---- crossform's own native RSA view (OLS in RDM space) --------------
# Reported alongside, NOT as the matched quantity: rsa() fits an ordinary
# least-squares regression of the RDM on the model RDM (plus intercept),
# which is a linear estimand, not the rank-based one above.
t0 <- Sys.time()
fit <- rsa(plan, models = list(animacy = MODEL_RDMS$animacy))
elapsed_rsa <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
rsa_ols <- as.numeric(as.matrix(fit$coefficients)[, "animacy"])
message("rsa() over ", n_centers, " centres in ", round(elapsed_rsa, 2), " s")

## ---- Save ---------------------------------------------------------------
out <- list(
  arm = "crossform",
  estimand = paste("cross-validated squared Euclidean distance over 66",
                   "off-diagonal run pairs, identity metric, local frame",
                   "normalisation"),
  centers = centers,
  n_voxels = vapply(nbhd, length, integer(1)),
  rdm = rdm_cv,
  pairs = D$pairs,
  rsa_score = rsa_score,          # shared Spearman statistic
  rsa_ols = rsa_ols,              # crossform-native OLS coefficient
  model = "animacy",
  radius_mm = RADIUS_MM,
  verify_max_abs = verify_max_abs,
  timing = list(rdm_seconds = elapsed_rdm, rsa_seconds = elapsed_rsa),
  session = list(crossform = as.character(utils::packageVersion("crossform")),
                 neuroim2 = as.character(utils::packageVersion("neuroim2")),
                 when = Sys.time())
)
saveRDS(out, file.path(paths$results, "smoke-crossform.rds"))
message("Wrote ", file.path(paths$results, "smoke-crossform.rds"))
message("  Spearman RSA: median ", round(median(rsa_score, na.rm = TRUE), 4),
        "  range [", round(min(rsa_score, na.rm = TRUE), 4), ", ",
        round(max(rsa_score, na.rm = TRUE), 4), "]",
        "  NA: ", sum(is.na(rsa_score)))
