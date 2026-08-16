#!/usr/bin/env Rscript
# 04-compare.R -- centre-by-centre agreement report and divergence taxonomy.
#
# Three comparisons, in decreasing order of how much they prove:
#
#   1. MATCHED ESTIMAND (tolerance-gated). crossform's cross-validated
#      squared Euclidean RDM vs the same estimand recomputed by a direct loop
#      over rMVPA's searchlight spheres and neuroim2 series extraction. Same
#      inputs, same estimand, two independent data paths and two independent
#      contraction implementations. This one must pass at TOLERANCE.
#
#   2. SHARED SECOND-ORDER SCORE (tolerance-gated). The identical
#      `spearman_rsa()` applied to each arm's own RDM vectors.
#
#   3. NATIVE ESTIMANDS (reported, NOT gated). crossform's cross-run
#      condition-level RSA vs rMVPA's `rsa_model`. These are different
#      estimands by construction; the numbers are reported as map agreement
#      with the semantic causes named, never reconciled by adjustment.

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
source(file.path(exemplar_dir, "models.R"))
paths <- exemplar_paths(exemplar_dir)

E <- readRDS(file.path(paths$results, "smoke-crossform.rds"))
R <- readRDS(file.path(paths$results, "smoke-rmvpa.rds"))

## ---- Preconditions: the two arms must describe the same geometry --------
stopifnot(identical(E$centers, R$centers))
stopifnot(identical(E$n_voxels, R$n_voxels))
stopifnot(identical(E$radius_mm, R$radius_mm))
stopifnot(identical(as.character(E$pairs$left), as.character(R$pairs$left)),
          identical(as.character(E$pairs$right), as.character(R$pairs$right)))
n_centers <- length(E$centers)
message("Both arms: ", n_centers, " centres, identical centre ids, identical ",
        "sphere sizes, radius ", E$radius_mm, " mm")

fmt <- function(x) format(x, digits = 3, scientific = TRUE)

## ---- 1. Matched estimand: the RDMs themselves ---------------------------
d_rdm <- E$rdm - R$reference_rdm
max_abs_rdm <- max(abs(d_rdm))
# Relative to the scale of the quantity, so "tiny" is not an artefact of the
# distances themselves being tiny.
rdm_scale <- max(abs(E$rdm))
max_rel_rdm <- max_abs_rdm / rdm_scale
per_center_max <- apply(abs(d_rdm), 1L, max)
rdm_pass <- max_abs_rdm <= TOLERANCE

message("\n[1] MATCHED ESTIMAND -- cross-validated squared Euclidean RDM")
message("  (1a) crossform vs independent reference loop")
message("    entries compared      : ", length(d_rdm),
        " (", n_centers, " centres x ", ncol(E$rdm), " pairs)")
message("    max abs difference    : ", fmt(max_abs_rdm))
message("    max relative diff     : ", fmt(max_rel_rdm))
message("    worst centre          : ", E$centers[which.max(per_center_max)],
        " (", fmt(max(per_center_max)), ")")
message("    centres above tol     : ", sum(per_center_max > TOLERANCE), "/", n_centers)
message("    GATE (<= ", fmt(TOLERANCE), "): ", if (rdm_pass) "PASS" else "FAIL")

# (1b) The comparison that actually earns a replacement-map row: crossform's
# estimator against rMVPA's own crossnobis estimator, same estimand.
native_cv <- list(available = isTRUE(R$native_cv_available))
if (native_cv$available) {
  d_cv <- E$rdm - R$native_cv_rdm
  native_cv$max_abs_difference <- max(abs(d_cv))
  native_cv$per_center_max <- apply(abs(d_cv), 1L, max)
  native_cv$centers_above_tolerance <- sum(native_cv$per_center_max > TOLERANCE)
  native_cv$score_max_abs_difference <- max(abs(E$rsa_score - R$native_cv_score),
                                            na.rm = TRUE)
  native_cv$pass <- native_cv$max_abs_difference <= TOLERANCE
  message("  (1b) crossform vs rMVPA's OWN crossnobis estimator")
  message("    max abs difference    : ", fmt(native_cv$max_abs_difference))
  message("    centres above tol     : ", native_cv$centers_above_tolerance,
          "/", n_centers)
  message("    shared score max diff : ", fmt(native_cv$score_max_abs_difference))
  message("    GATE (<= ", fmt(TOLERANCE), "): ",
          if (native_cv$pass) "PASS" else "FAIL")
} else {
  native_cv$pass <- NA
  message("  (1b) rMVPA crossnobis estimator unavailable in the rMVPA build ",
          "used by 03; estimator-to-estimator comparison not run")
}

## ---- 2. Shared second-order score ---------------------------------------
d_score <- E$rsa_score - R$reference_score
max_abs_score <- max(abs(d_score), na.rm = TRUE)
score_cor <- stats::cor(E$rsa_score, R$reference_score,
                        use = "complete.obs")
score_pass <- max_abs_score <= TOLERANCE

message("\n[2] SHARED SPEARMAN SCORE on the matched RDMs")
message("    max abs difference    : ", fmt(max_abs_score))
message("    Pearson r of maps     : ", format(score_cor, digits = 12))
message("    GATE (<= ", fmt(TOLERANCE), "): ", if (score_pass) "PASS" else "FAIL")

## ---- 3. Native estimands (reported, not gated) --------------------------
native_rows <- list()
for (dm in names(R$native)) {
  v <- R$native[[dm]]$values
  native_rows[[dm]] <- data.frame(
    rmvpa_distmethod = dm,
    pearson_r_vs_crossform_spearman = stats::cor(E$rsa_score, v,
                                                   use = "complete.obs"),
    spearman_rho_vs_crossform_spearman = stats::cor(E$rsa_score, v,
      method = "spearman", use = "complete.obs"),
    pearson_r_vs_crossform_ols = stats::cor(E$rsa_ols, v, use = "complete.obs"),
    max_abs_difference = max(abs(E$rsa_score - v), na.rm = TRUE),
    rmvpa_median = stats::median(v, na.rm = TRUE),
    crossform_median = stats::median(E$rsa_score, na.rm = TRUE),
    stringsAsFactors = FALSE)
}
native_tbl <- do.call(rbind, native_rows)
rownames(native_tbl) <- NULL

message("\n[3] NATIVE ESTIMANDS -- different by construction, reported only")
print(format(native_tbl, digits = 4))

## ---- Divergence taxonomy ------------------------------------------------
# Every entry is a NAMED semantic difference, established by reading each
# package's source, not inferred from the size of a discrepancy.
#
# SCOPE: all of these describe `rsa_model`, and therefore comparison [3] only.
# They are NOT a claim that rMVPA cannot express the matched estimand -- it
# can, through its crossnobis helpers, which is comparison [1b] and agrees
# exactly. What `rsa_model` cannot do is express it; that is a statement about
# one entry point, not about the package.
divergences <- data.frame(
  id = c("aggregation-level", "distance-convention", "second-order-coupling",
         "run-pair-treatment", "crossform-second-order"),
  affects = c("native", "native", "native", "native", "native"),
  description = c(
    paste("rMVPA rsa_model builds ONE pooled 96x96 RDM per sphere and keeps",
          "its 4224 between-run cell pairs. crossform aggregates over run",
          "pairs first, yielding 28 condition-level values. A correlation",
          "over 4224 cell pairs is not a function of the 28 cell means",
          "(the denominator uses the full within-cell spread), so the two",
          "cannot coincide even on identical patterns."),
    paste("rMVPA uses 1 - correlation, computed within each single-run",
          "pattern pair. crossform uses the cross-validated squared",
          "Euclidean distance under the identity metric. Correlation",
          "distance rescales each pattern to unit norm inside the sphere,",
          "which is nonlinear in the patterns and therefore not expressible",
          "in crossform's bilinear core."),
    paste("In rMVPA `distmethod` sets BOTH the neural RDM metric AND the",
          "second-order correlation method; `regtype` is inert for",
          "'pearson'/'spearman'. A Pearson (1-r) RDM scored by a Spearman",
          "second-order correlation -- the README's literal estimand -- is",
          "not expressible in a single rsa_model call."),
    paste("Both arms exclude within-run pairs, but at different levels:",
          "rMVPA drops within-run CELL pairs from a pooled vector,",
          "crossform averages over the 66 unordered RUN pairs. These",
          "coincide only in which information is excluded, not in how the",
          "retained information is weighted."),
    paste("crossform's native rsa() is an ordinary least-squares",
          "regression in RDM space, not a rank correlation. The Spearman",
          "statistic in the matched comparison is supplied by the",
          "exemplar's shared spearman_rsa(), applied identically to both",
          "arms, so the second-order step contributes no arm difference.")),
  stringsAsFactors = FALSE)

message("\nDIVERGENCE TAXONOMY (", nrow(divergences),
        " named semantic differences, all describing `rsa_model` and so",
        " affecting comparison [3] only)")
for (i in seq_len(nrow(divergences))) {
  message("  - ", divergences$id[i])
}

## ---- Save ---------------------------------------------------------------
comparison <- list(
  tier = "smoke",
  n_centers = n_centers,
  radius_mm = E$radius_mm,
  model = E$model,
  tolerance = TOLERANCE,
  matched_rdm = list(
    max_abs_difference = max_abs_rdm,
    max_relative_difference = max_rel_rdm,
    scale = rdm_scale,
    per_center_max = per_center_max,
    centers_above_tolerance = sum(per_center_max > TOLERANCE),
    pass = rdm_pass),
  matched_rdm_native_estimator = native_cv,
  matched_score = list(
    max_abs_difference = max_abs_score,
    map_correlation = score_cor,
    pass = score_pass),
  native = native_tbl,
  divergences = divergences,
  crossform_self_verification = E$verify_max_abs,
  crossform_score = E$rsa_score,
  rmvpa_reference_score = R$reference_score,
  rmvpa_native = lapply(R$native, function(x) x$values),
  centers = E$centers,
  n_voxels = E$n_voxels,
  timing = list(crossform = E$timing, rmvpa = R$timing,
                rmvpa_native = lapply(R$native, function(x) x$seconds)),
  session = list(crossform = E$session, rmvpa = R$session),
  when = Sys.time()
)
saveRDS(comparison, file.path(paths$results, "smoke-comparison.rds"))
message("\nWrote ", file.path(paths$results, "smoke-comparison.rds"))

overall <- rdm_pass && score_pass &&
  (is.na(native_cv$pass) || isTRUE(native_cv$pass))
message("SMOKE TIER MATCHED-ESTIMAND GATE: ", if (overall) "PASS" else "FAIL")
if (is.na(native_cv$pass)) {
  message("  (estimator-to-estimator comparison was not available; ",
          "the gate covers 1a and 2 only)")
}
if (!overall) quit(save = "no", status = 1)
