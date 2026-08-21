#!/usr/bin/env Rscript
# 01-glm.R -- per-subject, per-run GLM on the shared MNI lattice.
#
# Reads the fMRIPrep MNI BOLD, the raw `events.tsv`, and the fMRIPrep confound
# table, and writes one RDS per subject holding
#
#   trust  : 5 runs x 9 conditions x V voxels  (the native evidence)
#   shared : 1 matrix, 7 conditions x V voxels (the P^F fingerprint source)
#
# on that subject's own coverage mask. Nothing here knows about crossform:
# this is the boundary where fMRI becomes a per-run condition-mean matrix, and
# everything downstream consumes only the RDS.
#
# WHAT THE MODEL IS, in full, because a GLM described loosely is a GLM that
# cannot be checked:
#
#   response   percent signal change: 100 * (y_t - mean_t y) / mean_t y, per
#              voxel per run. Puts every subject in the same units, which
#              matters because the population layer sums across subjects.
#   conditions the trial_type levels of 00-common.R, as boxcars at the
#              recorded onsets with the recorded durations, convolved with the
#              SPM canonical double-gamma HRF and sampled at the volume times
#              t = TR * (0:(T-1)).
#   nuisance   EVERY trial type in the events file that is not a declared
#              condition -- `missed_trial` (risk 4: nuisance, never a
#              condition) and, for `sharedreward`, its six block-level levels
#              and two too-sparse neutral event types; the 6 rigid-body motion
#              parameters and their temporal derivatives (12 columns); every
#              `cosine\d+` column fMRIPrep emits, which IS the high-pass
#              filter and cannot be dropped; `framewise_displacement`; and one
#              indicator column per volume whose FD exceeds 0.5 mm.
#   selection  confound columns are chosen BY NAME PATTERN, never by position
#              or count -- the tables run 154 to 364 columns wide across the 84
#              runs (risk 5). aCompCor is deliberately NOT used: its component
#              count is run-specific, so including "all of them" would give
#              different runs different model complexity for no stated reason.
#   estimator  OLS by QR on the full design; the 9 (or 7) condition columns are
#              read off and everything else is projected out.
#
# WHAT IS NOT MODELLED: slice timing (the derivatives are not slice-time
# corrected and no reference slice is declared, so the volume time is used
# as-is), HRF derivatives, and any autocorrelation model. OLS betas under
# whitening-free estimation are unbiased; their standard errors would not be,
# which is one reason this slice reports uncertainty as uncalibrated.
#
# Idempotent: a subject whose RDS already exists AND was fitted under the
# current `GLM_MODEL_ID` is skipped. Change the design, bump the id, and every
# cached fit refits itself rather than silently surviving the change.
#
# Environment: SLICE2_DIR, SUBJECTS (comma/space separated subset).
# Output: data/derived/<subject>_betas.rds  (git-ignored)

source(file.path(
  if (nzchar(Sys.getenv("SLICE2_DIR"))) Sys.getenv("SLICE2_DIR") else
    normalizePath(dirname(sub("^--file=", "",
      grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))),
  "00-common.R"))

dir.create(DERIVED_DIR, showWarnings = FALSE, recursive = TRUE)

FD_SPIKE_MM <- 0.5

## ---- HRF ------------------------------------------------------------------
# SPM's canonical: difference of two gamma densities, peak 6 s, undershoot 16 s,
# undershoot weight 1/6, normalized to unit peak. Written out rather than taken
# from a package so that the exemplar has no HRF dependency and the reader can
# see exactly what was convolved.
spm_hrf <- function(t) {
  h <- stats::dgamma(t, shape = 6, rate = 1) -
       stats::dgamma(t, shape = 16, rate = 1) / 6
  h / max(h)
}

# Convolve a set of (onset, duration) boxcars with the HRF and sample at the
# volume acquisition times. Done on a 16x-oversampled grid so that onsets,
# which are recorded to microsecond precision, are not rounded to the TR.
convolve_events <- function(onsets, durations, n_vol, tr, oversample = 16L) {
  dt <- tr / oversample
  n_fine <- as.integer(ceiling(n_vol * tr / dt)) + as.integer(ceiling(32 / dt))
  stick <- numeric(n_fine)
  for (i in seq_along(onsets)) {
    a <- as.integer(floor(onsets[i] / dt)) + 1L
    b <- as.integer(floor((onsets[i] + max(durations[i], dt)) / dt)) + 1L
    a <- max(a, 1L); b <- min(max(b, a), n_fine)
    stick[a:b] <- stick[a:b] + 1
  }
  kern <- spm_hrf(seq(0, 32, by = dt))
  conv <- stats::convolve(stick, rev(kern), type = "open")[seq_len(n_fine)]
  # Volume times: t = TR * (0:(T-1)), i.e. the first volume is at t = 0.
  idx <- as.integer(round(tr * seq.int(0L, n_vol - 1L) / dt)) + 1L
  conv[pmin(idx, n_fine)]
}

## ---- Confounds ------------------------------------------------------------
read_confounds <- function(path, n_vol) {
  tab <- utils::read.delim(path, sep = "\t", header = TRUE, na.strings = "n/a",
                           check.names = FALSE)
  if (nrow(tab) != n_vol) {
    stop(basename(path), ": ", nrow(tab), " confound rows for ", n_vol,
         " volumes")
  }
  nm <- names(tab)
  motion <- grep("^(trans|rot)_[xyz](_derivative1)?$", nm, value = TRUE)
  cosine <- grep("^cosine[0-9]+$", nm, value = TRUE)
  fd     <- intersect("framewise_displacement", nm)
  if (length(motion) != 12L) {
    stop(basename(path), ": expected 12 motion columns, found ", length(motion))
  }
  if (!length(cosine)) {
    stop(basename(path), ": no cosine drift basis; refusing to fit unfiltered")
  }
  X <- as.matrix(tab[, c(motion, cosine, fd), drop = FALSE])
  # fMRIPrep leaves the first sample of every derivative and of FD as n/a.
  # Zero is the right fill for a derivative at t = 1 and for FD at t = 1: both
  # are differences that do not exist yet, not missing measurements.
  X[!is.finite(X)] <- 0
  spikes <- NULL
  if (length(fd)) {
    hit <- which(X[, fd] > FD_SPIKE_MM)
    if (length(hit)) {
      spikes <- matrix(0, n_vol, length(hit),
                       dimnames = list(NULL, sprintf("spike%02d", seq_along(hit))))
      spikes[cbind(hit, seq_along(hit))] <- 1
    }
  }
  list(matrix = cbind(X, spikes), n_motion = length(motion),
       n_cosine = length(cosine), n_spike = if (is.null(spikes)) 0L
                                            else ncol(spikes),
       n_fd_over = if (length(fd)) sum(X[, fd] > FD_SPIKE_MM) else NA_integer_,
       mean_fd = if (length(fd)) mean(X[, fd]) else NA_real_,
       max_fd  = if (length(fd)) max(X[, fd])  else NA_real_,
       n_columns_available = length(nm))
}

## ---- One run --------------------------------------------------------------
fit_run <- function(subject, task, run, conditions, voxel_linear, n_vol,
                    volume_length) {
  ev <- utils::read.delim(path_events(subject, task, run), sep = "\t",
                          header = TRUE, na.strings = "n/a",
                          stringsAsFactors = FALSE)
  present <- intersect(conditions, unique(ev$trial_type))
  if (!identical(sort(present), sort(conditions))) {
    stop(subject, " ", task, " ", raw_run(run), ": missing condition(s) ",
         paste(setdiff(conditions, present), collapse = ", "))
  }
  # EVERY trial type in the file gets a regressor, and betas are read off only
  # for the declared conditions. An event that happens and is not modelled does
  # not vanish -- its variance is absorbed by whatever regressor it correlates
  # with, which for this dataset is a regressor of interest.
  #
  # This matters most for `sharedreward`, whose 15 levels include six
  # block-level regressors and two neutral event types too sparse to use as
  # conditions (DECISION.md risk 7). Measured on `sub-104`: fitting the seven
  # conditions alone versus fitting all fifteen levels changes the resulting
  # seven-condition fingerprint beyond recognition -- the correlation between
  # the two fingerprints has a MEDIAN of 0.13 across 63,860 voxels, with a
  # quarter of voxels below -0.20. Omitting a modelled-away nuisance is not a
  # simplification here, it is a different operator.
  #
  # For `trust` the two are the same model: its only levels are the nine
  # conditions plus `missed_trial`, so nothing is added.
  model_levels <- union(conditions, sort(unique(ev$trial_type)))
  model_levels <- model_levels[!is.na(model_levels) & nzchar(model_levels)]
  Xc <- vapply(model_levels, function(lv) {
    k <- ev$trial_type == lv
    convolve_events(ev$onset[k], ev$duration[k], n_vol, TR)
  }, numeric(n_vol))
  colnames(Xc) <- model_levels

  cf <- read_confounds(path_confounds(subject, task, run), n_vol)
  X  <- cbind(Xc, cf$matrix, `(Intercept)` = 1)

  # Read the BOLD and go straight to the mask. The full 4D array is
  # 66*78*61*T doubles (~545 MB at T = 217); it is dropped before the solve.
  # Slicing by linear index one volume at a time keeps the extraction copy at
  # one volume rather than materializing a T*V index matrix.
  arr <- RNifti::readNifti(path_bold(subject, task, run))
  if (!identical(dim(arr)[4L], n_vol)) {
    stop(subject, " ", task, " ", deriv_run(run), ": ", dim(arr)[4L],
         " volumes, expected ", n_vol)
  }
  nv <- length(voxel_linear)
  Y <- matrix(0, n_vol, nv)                    # T x V
  off <- (seq_len(n_vol) - 1L) * volume_length
  for (t in seq_len(n_vol)) Y[t, ] <- arr[voxel_linear + off[t]]
  rm(arr); gc(verbose = FALSE)

  # fMRIPrep's brain mask is generous at the edge: a handful of voxels inside
  # it are identically zero for the whole run, so percent signal change is
  # undefined there. They are not an error and they are not imputed -- they are
  # flagged, zeroed so the solve stays finite, and dropped from the subject's
  # native territory once every run has voted (see `degenerate` below). A voxel
  # is native territory only if it carries signal in ALL SEVEN runs.
  mu <- colMeans(Y)
  bad <- !is.finite(mu) | mu == 0
  mu[bad] <- 1
  Y <- 100 * sweep(sweep(Y, 2L, mu, "-"), 2L, mu, "/")   # percent signal change
  if (any(bad)) Y[, bad] <- 0

  qrX <- qr(X)
  if (qrX$rank < ncol(X)) {
    stop(subject, " ", task, " ", deriv_run(run), ": design is rank ",
         qrX$rank, " of ", ncol(X), " columns")
  }
  B <- qr.coef(qrX, Y)[conditions, , drop = FALSE]
  resid_df <- n_vol - qrX$rank
  list(beta = B, resid_df = resid_df, n_columns = ncol(X),
       bad = bad, confounds = cf,
       model_levels = model_levels,
       n_nuisance_levels = length(setdiff(model_levels, conditions)),
       n_missed = sum(ev$trial_type %in% NUISANCE_TRIAL_TYPES),
       trial_counts = table(factor(ev$trial_type, levels = conditions)))
}

## ---- Coverage mask --------------------------------------------------------
# A subject's native territory is the intersection of that subject's own seven
# run brain masks -- five `trust`, two `sharedreward`. Intersecting across
# tasks matters: the fingerprint that builds P^F and the evidence that P^F
# carries have to live on the same voxels, or the transport would be defined
# on rows the data does not reach.
coverage_mask <- function(subject, ref_dim) {
  acc <- NULL
  for (r in TRUST_RUNS) {
    a <- as.array(neuroim2::read_vol(path_boldmask(subject, "trust", r))) > 0
    acc <- if (is.null(acc)) a else acc & a
  }
  for (r in SHARED_RUNS) {
    a <- as.array(neuroim2::read_vol(path_boldmask(subject, "sharedreward", r))) > 0
    acc <- acc & a
  }
  stopifnot(identical(dim(acc), ref_dim))
  acc
}

## ---- Main -----------------------------------------------------------------
requested <- Sys.getenv("SUBJECTS")
requested <- if (nzchar(requested)) {
  trimws(strsplit(requested, "[,[:space:]]+")[[1L]])
} else {
  SUBJECTS
}
requested <- requested[nzchar(requested)]

ref_space <- reference_space()
ref_dim   <- as.integer(dim(ref_space))
script_t0 <- Sys.time()

for (subject in requested) {
  out <- path_betas(subject)
  if (file.exists(out)) {
    # Idempotent, but not blindly so. A cached fit is reused only if it was
    # produced by THIS model. A cache that does not notice the design changing
    # under it is how a pipeline silently keeps serving numbers that nobody
    # computes any more -- and this one did, once: an earlier version modelled
    # only the seven `sharedreward` conditions and left the other eight trial
    # types out, which changed the fingerprint beyond recognition. Bump
    # GLM_MODEL_ID whenever the design changes.
    cached <- tryCatch(readRDS(out)$model_id, error = function(e) NULL)
    if (identical(cached, GLM_MODEL_ID)) {
      say("%-9s betas present (%s), skipping", subject, GLM_MODEL_ID); next
    }
    say("%-9s betas were fitted by \"%s\"; refitting under \"%s\"",
        subject, if (is.null(cached)) "an unrecorded model" else cached,
        GLM_MODEL_ID)
  }
  t0 <- Sys.time()
  needed <- c(vapply(TRUST_RUNS, function(r) path_bold(subject, "trust", r), ""),
              vapply(SHARED_RUNS, function(r) path_bold(subject, "sharedreward", r), ""))
  absent <- needed[!file.exists(needed)]
  if (length(absent)) {
    say("%-9s SKIPPED: %d BOLD file(s) absent (run ./fetch.sh)",
        subject, length(absent))
    next
  }

  mask_arr <- coverage_mask(subject, ref_dim)
  # `which(arr)` is column-major, which is the order
  # `neuroim2_volume_domain()` reports features in. 02 asserts that alignment
  # against the domain's own coordinates rather than trusting this comment.
  lin <- which(mask_arr)
  vox <- arrayInd(lin, ref_dim)
  volume_length <- prod(ref_dim)

  trust <- lapply(TRUST_RUNS, function(r)
    fit_run(subject, "trust", r, TRUST_CONDITIONS, lin, 217L, volume_length))
  names(trust) <- paste0("run", TRUST_RUNS)

  # The fingerprint source. Both `sharedreward` runs are fitted separately and
  # averaged: a single joint fit would need a run factor, and averaging two
  # per-run beta maps is the same estimate with less code to get wrong.
  shared <- lapply(SHARED_RUNS, function(r)
    fit_run(subject, "sharedreward", r, SHARED_CONDITIONS, lin, 202L,
            volume_length))
  shared_beta <- Reduce(`+`, lapply(shared, `[[`, "beta")) / length(shared)

  # Voxels that do not vary in some run carry no geometry and would make a
  # searchlight singular; they are dropped from the native territory here, so
  # that the domain, the frame and the transport all see one voxel set. The
  # union is taken over all seven runs of both tasks.
  degenerate <- Reduce(`|`, c(lapply(trust, `[[`, "bad"),
                              lapply(shared, `[[`, "bad")))
  degenerate <- degenerate | Reduce(`|`, lapply(trust, function(x)
    !is.finite(colSums(x$beta)) | apply(x$beta, 2L, stats::var) == 0))
  degenerate <- degenerate | !is.finite(colSums(shared_beta))
  keep <- !degenerate

  saveRDS(list(
    model_id = GLM_MODEL_ID,
    subject = subject,
    group = if (subject %in% SUBJECTS_YOUNGER) "younger" else "older",
    dim = ref_dim,
    voxel_linear = lin[keep],
    voxel_index = vox[keep, , drop = FALSE],
    conditions = TRUST_CONDITIONS,
    shared_conditions = SHARED_CONDITIONS,
    trust = lapply(trust, function(x)
      x$beta[, keep, drop = FALSE]),
    shared = shared_beta[, keep, drop = FALSE],
    n_coverage_voxels = length(lin),
    n_degenerate = sum(degenerate),
    resid_df = vapply(trust, `[[`, numeric(1), "resid_df"),
    design_columns = vapply(trust, `[[`, numeric(1), "n_columns"),
    n_missed = vapply(trust, `[[`, numeric(1), "n_missed"),
    trust_nuisance_levels = vapply(trust, `[[`, numeric(1), "n_nuisance_levels"),
    shared_nuisance_levels = vapply(shared, `[[`, numeric(1), "n_nuisance_levels"),
    shared_model_levels = shared[[1L]]$model_levels,
    n_spike = vapply(trust, function(x) x$confounds$n_spike, numeric(1)),
    n_cosine = vapply(trust, function(x) x$confounds$n_cosine, numeric(1)),
    confound_columns_available = vapply(trust, function(x)
      x$confounds$n_columns_available, numeric(1)),
    mean_fd = vapply(trust, function(x) x$confounds$mean_fd, numeric(1)),
    max_fd = vapply(trust, function(x) x$confounds$max_fd, numeric(1)),
    trial_counts = do.call(rbind, lapply(trust, `[[`, "trial_counts")),
    shared_resid_df = vapply(shared, `[[`, numeric(1), "resid_df"),
    fit_seconds = elapsed(t0)
  ), out)

  say("%-9s V=%6d (dropped %3d degenerate)  df=%s  spikes=%s  %.0fs",
      subject, sum(keep), sum(degenerate),
      paste(vapply(trust, `[[`, numeric(1), "resid_df"), collapse = "/"),
      paste(vapply(trust, function(x) x$confounds$n_spike, numeric(1)),
            collapse = "/"),
      elapsed(t0))
}

say("\n01-glm.R done in %.1f min", elapsed(script_t0) / 60)
