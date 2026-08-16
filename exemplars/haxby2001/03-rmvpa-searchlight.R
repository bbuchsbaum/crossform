#!/usr/bin/env Rscript
# 03-rmvpa-searchlight.R -- the rMVPA arm.
#
# Consumes the SAME data/prepared-smoke.rds as 02. The 12 x 8 per-run
# condition means are packed into one 96-row NeuroVec (row order: run-major,
# condition-minor) and handed to rMVPA's own searchlight machinery.
#
# Three things are computed. (A) is a different estimand from (B) and (C);
# (B) and (C) are the same estimand by two independent routes.
#
#  (A) rMVPA NATIVE RSA -- `rsa_model` + `rsa_design(block_var = run)`.
#      rMVPA builds ONE pooled 96 x 96 correlation-distance matrix per sphere
#      and keeps its 4224 between-run cell pairs, then correlates that vector
#      against the model RDM expanded to 96 x 96. It never aggregates to the
#      8 x 8 condition level and never averages over run pairs. Run in both
#      distmethod settings because in rMVPA `distmethod` sets BOTH the neural
#      distance and the second-order correlation method (`regtype` is inert
#      for "pearson"/"spearman"; see 04's divergence notes).
#
#  (B) MATCHED-ESTIMAND REFERENCE -- the same cross-validated squared
#      Euclidean RDM that 02 computes, recomputed here by a plain nested loop
#      written from the estimand definition, fed by rMVPA's own searchlight
#      spheres. Depends on neither package's estimator, so it arbitrates
#      between them.
#
#  (C) rMVPA's OWN CROSSNOBIS ESTIMATOR -- the matched estimand computed by
#      rMVPA's crossnobis machinery. `rsa_model` cannot express it, but
#      `compute_crossvalidated_means_sl` + `compute_crossnobis_distances_sl`
#      can, which makes the headline comparison estimator-to-estimator rather
#      than estimator-to-reference-loop. Needs a current rMVPA; see the
#      capability check below.
#
# Neither arm's numbers are tuned to the other.

suppressMessages({
  library(neuroim2)
  library(rMVPA)
})

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
source(file.path(exemplar_dir, "models.R"))
paths <- exemplar_paths(exemplar_dir)
dir.create(paths$results, showWarnings = FALSE, recursive = TRUE)

prep <- readRDS(paths$prepared)
stopifnot(identical(prep$conditions, CONDITIONS))

## ---- Pack the shared prepared matrices into a 96-row NeuroVec -----------
vt_vol <- neuroim2::read_vol(prep$mask_file)
vt_mask <- neuroim2::LogicalNeuroVol(as.array(vt_vol) > 0,
                                     neuroim2::space(vt_vol))
vt_index <- which(as.array(vt_mask))
stopifnot(identical(vt_index, prep$vt_index))

n_run <- length(prep$runs)
n_cond <- length(CONDITIONS)
n_obs <- n_run * n_cond
run_of <- rep(seq_len(n_run), each = n_cond)
cond_of <- rep(seq_len(n_cond), times = n_run)

M <- matrix(NA_real_, n_obs, length(vt_index))
for (r in seq_len(n_run)) {
  M[(r - 1L) * n_cond + seq_len(n_cond), ] <- prep$runs[[r]]
}
stopifnot(!anyNA(M))

# Extend the mask's own 3D space to 4D with add_dim so spacing, origin, and
# orientation are carried over. Building NeuroSpace(c(dim, n_obs)) directly
# drops them, which older rMVPA accepted and current rMVPA correctly rejects
# with "Spatial geometry mismatch between `train_data` and `mask`".
svec <- neuroim2::SparseNeuroVec(
  M, neuroim2::add_dim(neuroim2::space(vt_mask), n_obs), mask = vt_mask)
dset <- rMVPA::mvpa_dataset(svec, mask = vt_mask)

# The pack must round-trip through neuroim2 exactly, or every number below is
# about a different dataset than 02 used.
probe <- neuroim2::series(svec, vt_index[seq_len(min(25, length(vt_index)))])
stopifnot(max(abs(probe - M[, seq_len(ncol(probe)), drop = FALSE])) == 0)
message("96-row dataset packed and round-trip verified against prepared data")

## ---- (A) rMVPA native RSA searchlight -----------------------------------
model_rdm <- MODEL_RDMS$animacy[CONDITIONS, CONDITIONS, drop = FALSE]
R96 <- model_rdm[cond_of, cond_of]          # expand to 96 x 96, symmetric, 0 diag
stopifnot(isSymmetric(R96), all(diag(R96) == 0))

rsa_des <- rMVPA::rsa_design(~ animacy, data = list(animacy = R96),
                             block_var = run_of)
n_between <- length(rsa_des$model_mat$animacy)
message("rsa_design: ", n_between, " between-run cell pairs kept of ",
        choose(n_obs, 2), " (dropped ", n_run * choose(n_cond, 2),
        " within-run pairs)")
stopifnot(n_between == choose(n_obs, 2) - n_run * choose(n_cond, 2))

native <- list()
for (dm in c("pearson", "spearman")) {
  mspec <- rMVPA::rsa_model(dset, rsa_des, distmethod = dm, regtype = "pearson")
  t0 <- Sys.time()
  res <- rMVPA::run_searchlight(mspec, radius = RADIUS_MM, method = "standard")
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  vals <- res$results$animacy[vt_index]
  native[[dm]] <- list(values = as.numeric(vals), seconds = el)
  message("rMVPA rsa_model distmethod=", dm, ": ", round(el, 1), " s, median ",
          round(median(vals, na.rm = TRUE), 4), ", NA ", sum(is.na(vals)))
}

## ---- (B) Matched-estimand reference over rMVPA's data path --------------
# Cross-validated squared Euclidean distance, mean over the 66 unordered run
# pairs, with local (1/P) normalisation -- the README estimand, written as a
# direct loop with no reference to how crossform computes it.
sl <- rMVPA::get_searchlight(dset, "standard", radius = RADIUS_MM)
center_ids <- rMVPA::get_center_ids(dset)
stopifnot(identical(as.integer(center_ids), as.integer(vt_index)))

pair_idx <- rdm_pair_index()
n_pairs <- nrow(pair_idx)
run_pairs <- utils::combn(n_run, 2L)

cv_rdm <- matrix(NA_real_, length(center_ids), n_pairs)
sphere_n <- integer(length(center_ids))

t0 <- Sys.time()
for (ci in seq_along(center_ids)) {
  # sl[[ci]] is a sphere over the MASK geometry; values() on it returns mask
  # values, not data. Pull the patterns with the sphere's full-volume indices.
  members <- neuroim2::indices(sl[[ci]])
  X <- neuroim2::series(svec, members)      # 96 x P, rows in 4th-dim order
  P <- ncol(X)
  sphere_n[ci] <- P
  # difference patterns per run: [n_pairs x P] for each run
  acc <- numeric(n_pairs)
  for (pp in seq_len(ncol(run_pairs))) {
    r <- run_pairs[1L, pp]; s <- run_pairs[2L, pp]
    Xr <- X[(r - 1L) * n_cond + seq_len(n_cond), , drop = FALSE]
    Xs <- X[(s - 1L) * n_cond + seq_len(n_cond), , drop = FALSE]
    Dr <- Xr[pair_idx$i, , drop = FALSE] - Xr[pair_idx$j, , drop = FALSE]
    Ds <- Xs[pair_idx$i, , drop = FALSE] - Xs[pair_idx$j, , drop = FALSE]
    acc <- acc + rowSums(Dr * Ds) / P
  }
  cv_rdm[ci, ] <- acc / ncol(run_pairs)
  if (ci %% 100 == 0) message("  reference RDM centre ", ci, "/", length(center_ids))
}
elapsed_ref <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
message("matched-estimand reference over rMVPA spheres in ",
        round(elapsed_ref, 1), " s")

model_vec <- model_vector(MODEL_RDMS$animacy)
ref_score <- apply(cv_rdm, 1L, spearman_rsa, model_vec = model_vec)

## ---- (C) rMVPA's OWN crossnobis estimator -------------------------------
# `rsa_model` cannot express the matched estimand, but rMVPA's crossnobis
# machinery can: compute_crossvalidated_means_sl() builds per-run condition
# means and compute_crossnobis_distances_sl() forms the cross-validated
# squared distance, dividing by P once and averaging over ordered partition
# pairs (equivalent to unordered, since the summand is symmetric in r,s).
# This is what makes the comparison estimator-to-estimator rather than
# estimator-to-reference-loop.
#
# Requires a current rMVPA: the 2026-04 build both demands a whitening matrix
# for crossnobis and uses leave-one-run-out TRAINING means rather than
# per-run means, so it cannot produce this estimand at all. Capability is
# checked rather than assumed, and the arm is skipped (not faked) if absent.
have_crossnobis <-
  exists("compute_crossvalidated_means_sl", where = asNamespace("rMVPA")) &&
  exists("compute_crossnobis_distances_sl", where = asNamespace("rMVPA"))

native_cv_rdm <- NULL
native_cv_score <- NULL
elapsed_native_cv <- NA_real_
if (!have_crossnobis) {
  message("rMVPA crossnobis helpers unavailable in this build; arm (C) skipped")
} else {
  mv_design <- rMVPA::mvpa_design(
    data.frame(y = factor(rep(CONDITIONS, times = n_run), levels = CONDITIONS),
               block = run_of),
    y_train = ~ y, block_var = ~ block)
  cv_spec <- rMVPA::blocked_cross_validation(run_of)

  # Verify the SEMANTICS, not the version. The 2026-04 rMVPA builds crossnobis
  # "folds" from leave-one-run-out TRAINING means; those overlap across folds,
  # so cross-fold inner products do not cancel noise and the distances come out
  # badly positively biased. Current rMVPA uses each run's own mean. Check that
  # directly on a small slice rather than trusting a version string, because a
  # silently biased arm is exactly the failure this exemplar exists to catch.
  probe_cols <- seq_len(min(20L, ncol(M)))
  U_probe <- rMVPA::compute_crossvalidated_means_sl(
    M[, probe_cols, drop = FALSE], mv_design, cv_spec,
    estimation_method = "crossnobis", return_folds = TRUE)$fold_estimates
  per_run <- array(NA_real_, dim(U_probe))
  for (r in seq_len(n_run)) {
    per_run[, , r] <- M[(r - 1L) * n_cond + seq_len(n_cond), probe_cols]
  }
  fold_gap <- max(abs(U_probe - per_run))
  if (fold_gap > 1e-12) {
    stop("rMVPA's crossnobis fold estimates are not this run's per-run ",
         "condition means (max abs difference ", format(fold_gap, digits = 3),
         "). This build derives folds from leave-one-run-out training means, ",
         "which yields positively biased distances. Reinstall rMVPA from a ",
         "current source tree before trusting arm (C).")
  }
  message("crossnobis fold construction verified: folds are per-run means ",
          "(max abs difference ", format(fold_gap, digits = 3), ")")

  native_cv_rdm <- matrix(NA_real_, length(center_ids), n_pairs)
  t0 <- Sys.time()
  for (ci in seq_along(center_ids)) {
    X <- neuroim2::series(svec, neuroim2::indices(sl[[ci]]))
    U <- rMVPA::compute_crossvalidated_means_sl(
      X, mv_design, cv_spec, estimation_method = "crossnobis",
      return_folds = TRUE)$fold_estimates
    d <- rMVPA:::compute_crossnobis_distances_sl(U, P_voxels = ncol(X))
    if (ci == 1L) {
      # rMVPA labels pairs "<right>_vs_<left>" in lower-triangle order, which
      # is the same order as combn(); assert rather than trust it.
      stopifnot(identical(names(d),
                          paste0(pair_idx$right, "_vs_", pair_idx$left)))
    }
    native_cv_rdm[ci, ] <- as.numeric(d)
  }
  elapsed_native_cv <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  native_cv_score <- apply(native_cv_rdm, 1L, spearman_rsa, model_vec = model_vec)
  message("rMVPA native crossnobis RDM over ", length(center_ids),
          " centres in ", round(elapsed_native_cv, 1), " s")
}

## ---- Save ---------------------------------------------------------------
out <- list(
  arm = "rMVPA",
  centers = as.integer(center_ids),
  n_voxels = sphere_n,
  native = native,               # (A) rMVPA's own rsa_model, per distmethod
  native_estimand = paste("pooled 96x96 correlation-distance RDM per sphere,",
                          "between-run cell pairs only (4224 of 4560),",
                          "correlated with the expanded 96x96 model RDM"),
  reference_rdm = cv_rdm,        # (B) matched estimand over rMVPA's data path
  reference_score = ref_score,
  reference_estimand = paste("cross-validated squared Euclidean distance over",
                             "66 off-diagonal run pairs, local normalisation"),
  native_cv_rdm = native_cv_rdm,      # (C) rMVPA's own crossnobis estimator
  native_cv_score = native_cv_score,
  native_cv_available = have_crossnobis,
  native_cv_estimand = paste("compute_crossvalidated_means_sl(estimation_method",
                             "= 'crossnobis') + compute_crossnobis_distances_sl,",
                             "identity metric, per-run partition means"),
  pairs = pair_idx[, c("left", "right")],
  model = "animacy",
  radius_mm = RADIUS_MM,
  timing = list(reference_seconds = elapsed_ref,
                native_cv_seconds = elapsed_native_cv),
  session = list(rMVPA = as.character(utils::packageVersion("rMVPA")),
                 rMVPA_build = read.dcf(file.path(find.package("rMVPA"),
                                                  "DESCRIPTION"))[1, "Built"],
                 rMVPA_path = find.package("rMVPA"),
                 neuroim2 = as.character(utils::packageVersion("neuroim2")),
                 when = Sys.time())
)
saveRDS(out, file.path(paths$results, "smoke-rmvpa.rds"))
message("Wrote ", file.path(paths$results, "smoke-rmvpa.rds"))
