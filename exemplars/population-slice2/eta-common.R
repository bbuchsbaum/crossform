#!/usr/bin/env Rscript
# eta-common.R -- the machinery both eta scripts share.
#
# 04 measures eta across TASK: P^F fitted on `sharedreward`, evaluated on
# `trust`. 05 measures it across RUN: P^F fitted on `trust` run 1, evaluated on
# `trust` runs 2-5. The two differ only in where the fingerprint comes from and
# which runs are held out, so everything else lives here -- and lives here
# rather than being copied, because the null band is only a band around the
# honest estimator if the honest estimator and the 200 permuted ones are the
# same arithmetic.
#
# Sourced after 00-common.R and after `load_crossform()`.

## ---- Held-out native forms ------------------------------------------------
# `materialize_geometry()` gives the complete packed form at every native
# searchlight: a V x 45 matrix of svec'd 9 x 9 condition geometries. Transport
# is LINEAR in these, so once they exist any transport's Z is one sparse
# product away. That is the whole reason a 200-draw null band is affordable:
# the geometry is computed once and only the operator changes.
held_out_forms <- function(subjects, native_nodes, partitions, ref_space,
                           reporter = message) {
  dim3 <- as.integer(dim(ref_space))
  eff <- effect_space(TRUST_CONDITIONS,
                      basis_id = "ds003745-trust-condition-means:v1",
                      units = "percent-signal-change")
  out <- list()
  for (s in subjects) {
    t0 <- Sys.time()
    f <- readRDS(path_betas(s))
    mask_arr <- array(FALSE, dim3)
    mask_arr[f$voxel_linear] <- TRUE
    dom <- neuroim2_volume_domain(
      neuroim2::LogicalNeuroVol(mask_arr, ref_space),
      id = paste0("ds003745-", s))
    stopifnot(identical(as.character(dom$feature_ids), native_nodes[[s]]))
    frame <- compile_frame(
      searchlights(radius = SEARCHLIGHT_RADIUS_MM,
                   normalization = "conservative"), dom)
    part <- lapply(partitions, function(runs) {
      rs <- lapply(f$trust[runs], function(B) {
        dimnames(B) <- list(TRUST_CONDITIONS, as.character(dom$feature_ids))
        B
      })
      rel <- relation(rs, effects = eff, domain = dom)
      over <- cross_partitions(rel, independence = "independent",
                               generalizes_over = "run")
      M <- geometry_component(
        materialize_geometry(plan_geometry(rel, at = frame, over = over)),
        "total")
      stopifnot(nrow(M) == length(dom$feature_ids))
      M
    })
    out[[s]] <- part
    rm(frame, dom, f); gc(verbose = FALSE)
    reporter(sprintf("  %-9s V=%6d packed %d  %.0fs", s, nrow(part[[1L]]),
                     ncol(part[[1L]]), as.numeric(difftime(Sys.time(), t0,
                                                           units = "secs"))))
  }
  out
}

## ---- V^C, V^W, R (contract section 7.1, implemented literally) ------------
# Zbar_i is the average of the subject's two partition forms. The contract
# writes it without spelling the averaging out; the mean is the reading that
# makes V^C and V^W estimate the same ||mu||^2 up to the noise term, which is
# what makes their ratio a share.
#
# Sum over ORDERED pairs i != i', which is what the N(N-1) normalizer is for;
# it is computed as ||sum_i Zbar_i||^2 - sum_i ||Zbar_i||^2 rather than as a
# double loop, because at 12 subjects x 1241 nodes x 45 coordinates the double
# loop is the only slow part of a null draw.
consensus_share <- function(Pmats, forms, subjects, vw_floor = 0) {
  N <- length(subjects)
  ZA <- lapply(subjects, function(s)
    as.matrix(Matrix::crossprod(Pmats[[s]], forms[[s]][[1L]])))
  ZB <- lapply(subjects, function(s)
    as.matrix(Matrix::crossprod(Pmats[[s]], forms[[s]][[2L]])))
  # V^W is a MEAN over participants of a cross-partition inner product, so it
  # has a standard error and the per-participant terms are kept. "V^W > 0" is a
  # statement about a noisy average; without its spread there is no way to tell
  # a V^W that clears the floor from one that merely landed above it, and the
  # share R inherits whatever that average is worth.
  vw_i <- vapply(seq_len(N), function(i) sum(ZA[[i]] * ZB[[i]]), numeric(1))
  names(vw_i) <- subjects
  vw <- mean(vw_i)
  Zbar <- lapply(seq_len(N), function(i) (ZA[[i]] + ZB[[i]]) / 2)
  S <- Reduce(`+`, Zbar)
  self <- sum(vapply(Zbar, function(z) sum(z * z), numeric(1)))
  vc <- (sum(S * S) - self) / (N * (N - 1))
  # Contract 7.1: the share is formed only when V^W clears a declared floor.
  # V^W is itself a cross-partition product and can be zero or negative when
  # nothing reproduces; the share is then NA, not a large number.
  list(V_C = vc, V_W = vw, V_W_per_subject = vw_i,
       V_W_se = stats::sd(vw_i) / sqrt(N),
       R = if (is.finite(vw) && vw > vw_floor) vc / vw else NA_real_)
}

## ---- Leave-one-subject-out group fingerprint atlas ------------------------
# Per subject and group node, the mean fingerprint over that subject's voxels
# within the fingerprint radius; then the atlas a subject sees is the mean over
# the OTHER subjects only, centred and unit-normed so a dot product between two
# fingerprints is their correlation.
#
# Leaving the subject out matters even though the fitting data is already
# disjoint from the evaluation data: without it, a subject's own transport
# would be built partly from its own idiosyncrasy, and the group node it is
# steered toward would be the one it already resembles.
build_loo_atlas <- function(subjects, fingerprints, fp_supports, n_group) {
  q <- ncol(fingerprints[[subjects[1L]]])
  per <- array(0, c(length(subjects), n_group, q),
               dimnames = list(subjects, NULL, NULL))
  cnt <- matrix(0L, length(subjects), n_group, dimnames = list(subjects, NULL))
  for (s in subjects) {
    fs <- fp_supports[[s]]
    acc <- rowsum(fingerprints[[s]][fs[, "row"], , drop = FALSE], fs[, "col"],
                  reorder = TRUE)
    counts <- as.vector(table(factor(fs[, "col"], levels = seq_len(n_group))))
    hit <- as.integer(rownames(acc))
    per[s, hit, ] <- acc / counts[hit]
    cnt[s, ] <- counts
  }
  total <- apply(per, c(2L, 3L), sum)
  reach <- colSums(cnt > 0L)
  structure(lapply(stats::setNames(subjects, subjects), function(s) {
    n_other <- reach - as.integer(cnt[s, ] > 0L)
    A <- total - per[s, , ]
    A[n_other > 0L, ] <- A[n_other > 0L, , drop = FALSE] / n_other[n_other > 0L]
    A[n_other == 0L, ] <- 0
    A <- A - rowMeans(A)
    nrm <- sqrt(rowSums(A^2))
    ok <- nrm > .Machine$double.eps
    A[ok, ] <- A[ok, , drop = FALSE] / nrm[ok]
    A[!ok, ] <- 0
    list(atlas = A, defined = ok & n_other > 0L)
  }), reach = reach)
}

## ---- One softmax transport, optionally under a permuted atlas -------------
# `perm` maps group node j to the fingerprint it is given. The identity
# permutation must reproduce the operator `location_transport()` sealed, and
# both eta scripts assert exactly that before drawing a single null.
softmax_transport_matrix <- function(support, fingerprint, atlas, defined,
                                     degenerate, n_native, n_group,
                                     temperature, perm = NULL) {
  cols <- if (is.null(perm)) support[, "col"] else perm[support[, "col"]]
  sim <- rowSums(fingerprint[support[, "row"], , drop = FALSE] *
                 atlas[cols, , drop = FALSE])
  sim[!defined[cols]] <- 0
  sim[degenerate[support[, "row"]]] <- 0
  Matrix::sparseMatrix(
    i = support[, "row"], j = support[, "col"],
    x = softmax_rows(sim, support[, "row"], temperature),
    dims = c(n_native, n_group))
}

## ---- The null band --------------------------------------------------------
# One draw = one permutation of the group-node fingerprint atlas ACROSS group
# nodes: node j is given node pi(j)'s fingerprint. The permuted transport keeps
# the same support, the same softmax, the same temperature, the same sink and
# the same amount of spreading -- everything except the correspondence between
# a voxel's function and its destination's function.
#
# Which randomization is right is open maintainer decision 14.5; this is the
# one this slice picks, and it is picked because it is the only one that holds
# the smoothing fixed. Permuting subject labels or held-out partitions would
# leave P^F's extra spreading unpriced, and P^F's extra spreading is precisely
# the thing an anatomical baseline does not have.
#
# The draw loop itself lives in each eta script, because the closure that turns
# a permutation into a set of transports is the one thing the two scripts do
# not share -- 04's fingerprint is a `sharedreward` condition vector, 05's is a
# `trust` run-1 condition vector. This function summarizes the result.
summarize_null <- function(null, eta) {
  finite <- if (is.null(null) || !nrow(null)) numeric(0) else
    null$eta[is.finite(null$eta)]
  # An empty band is a legitimate state, not an error: it is what a V^W that
  # does not clear the floor produces, and it must report itself as absent
  # rather than as a band of NaNs with a p-value attached.
  if (!length(finite) || !is.finite(eta)) {
    return(list(mean = NA_real_, sd = NA_real_, q95 = NA_real_,
                max = NA_real_, negative_share = NA_real_,
                draws = length(finite), rank = NA_real_, p = NA_real_))
  }
  list(
    mean = mean(finite), sd = stats::sd(finite),
    q95 = unname(stats::quantile(finite, 0.95)), max = max(finite),
    negative_share = mean(finite < 0), draws = length(finite),
    rank = sum(finite < eta),
    p = (1 + sum(finite >= eta)) / (1 + length(finite))
  )
}
