#!/usr/bin/env Rscript
# 00-common.R -- shared design, simulation parameters, loaders, and the
# closed-form planted geometry for the rectangular encoding-retrieval RSA
# exemplar.
#
# Sourced by 01-simulate.R, 02-analyze.R, and 03-recover.R. Every constant
# that more than one script needs lives here so the three cannot drift apart.
#
# Scale knobs for a fast smoke run (they change the numbers, so the committed
# results are always produced with the defaults):
#   ER_SUBJECTS=3 Rscript 01-simulate.R

exemplar_dir_from_args <- function() {
  if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) return(Sys.getenv("EXEMPLAR_DIR"))
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) return(normalizePath(dirname(sub("^--file=", "", f[1]))))
  normalizePath(".")
}

load_crossform <- function(exemplar_dir) {
  repo <- normalizePath(file.path(exemplar_dir, "..", ".."), mustWork = FALSE)
  if (file.exists(file.path(repo, "DESCRIPTION"))) {
    suppressMessages(pkgload::load_all(repo, quiet = TRUE, export_all = FALSE))
    return("source")
  }
  suppressMessages(library(crossform))
  "installed"
}

## ---- Experimental design ------------------------------------------------
# Three study-test cycles over one 36-item list. Every cycle is one scanning
# run holding that cycle's encoding block and then its retrieval block, so
# encoding and retrieval effects for a cycle come out of the *same* GLM and
# their estimation errors are correlated within a run. That is what makes the
# cross-run pairing below a scientific choice rather than a formality.
#
# Repeated study-test cycles are also what crossform's relation contract
# needs: every partition of a relation estimates the whole shared effect
# space, so each run must contain every encoding item and every retrieval
# probe. A single-cycle design with items split across runs is a different
# partitioning problem and is out of scope here.

SEED <- 20260820L
N_SUBJECTS <- as.integer(Sys.getenv("ER_SUBJECTS", "12"))
N_RUNS <- 3L

CATEGORIES <- c("face", "scene", "object")
ITEMS_PER_CATEGORY <- 12L
RETRIEVED_PER_CATEGORY <- 8L   # encoded items that come back as old probes
LURES_PER_CATEGORY <- 2L       # retrieval probes that were never encoded

# Neural domain: three regions with different planted structure.
REGION_SIZES <- c(regionA = 40L, regionB = 30L, regionC = 20L)
REGION_LABELS <- rep(names(REGION_SIZES), REGION_SIZES)
N_VOXELS <- sum(REGION_SIZES)

# Acquisition.
TR <- 2
N_TR <- 300L                   # 600 s per run
ISI_JITTER <- c(3, 4, 5, 6, 7)
RETRIEVAL_DURATION <- 2.5
STUDY_DURATIONS <- seq(1.5, 4.0, by = 0.5)   # the trial covariate, in seconds

## ---- Planted effect amplitudes ------------------------------------------
# regionA carries item-specific reinstatement (plus a category component, so
# the naive control set is genuinely confounded).
# regionB carries category structure only.
# regionC carries no task structure at all.
AMP_ITEM_A <- 1.0        # encoding item pattern, regionA
AMP_CATEGORY_A <- 0.6    # encoding/retrieval category pattern, regionA
AMP_CATEGORY_B <- 1.2    # encoding/retrieval category pattern, regionB
AMP_LURE_A <- 0.8        # a lure's own (never-encoded) pattern, regionA
REINSTATEMENT_INTERCEPT <- 0.55   # rho0: reinstatement at mean study duration
REINSTATEMENT_SLOPE <- 0.18       # rho1: extra reinstatement per second
STATE_SD <- 0.5          # same-run, item-specific state pattern (all voxels)

# Noise and nuisance.
NOISE_SD <- 1.0
NOISE_AR <- 0.6          # spatial correlation within a region
DRIFT_SD <- 4.0
MOTION_SD <- 2.0
BASELINE_MEAN <- 100

TOLERANCE <- 1e-10

## ---- Item table ---------------------------------------------------------
er_items <- function() {
  item <- sprintf("item%02d", seq_len(length(CATEGORIES) * ITEMS_PER_CATEGORY))
  category <- rep(CATEGORIES, each = ITEMS_PER_CATEGORY)
  retrieved <- unlist(lapply(seq_along(CATEGORIES), function(g) {
    offset <- (g - 1L) * ITEMS_PER_CATEGORY
    seq_len(ITEMS_PER_CATEGORY) + offset <= offset + RETRIEVED_PER_CATEGORY
  }))
  set.seed(SEED + 1L)
  duration <- sample(STUDY_DURATIONS, length(item), replace = TRUE)
  data.frame(
    item = item, category = category, retrieved = retrieved,
    study_duration = duration, stringsAsFactors = FALSE
  )
}

er_probes <- function(items) {
  old <- items[items$retrieved, c("item", "category")]
  old$probe <- paste0("ret_", old$item)
  old$old <- TRUE
  old$source_item <- old$item
  lures <- do.call(rbind, lapply(seq_along(CATEGORIES), function(g) {
    index <- (g - 1L) * LURES_PER_CATEGORY + seq_len(LURES_PER_CATEGORY)
    data.frame(
      item = NA_character_, category = CATEGORIES[g],
      probe = sprintf("ret_lure%02d", index), old = FALSE,
      source_item = NA_character_, stringsAsFactors = FALSE
    )
  }))
  probes <- rbind(old, lures)
  probes <- probes[order(match(probes$category, CATEGORIES), !probes$old,
    probes$probe), ]
  rownames(probes) <- NULL
  probes
}

## ---- Haemodynamic design ------------------------------------------------
hrf_canonical <- function(t) {
  positive <- ifelse(t > 0, t^5 * exp(-t) / gamma(6), 0)
  negative <- ifelse(t > 0, t^15 * exp(-t) / gamma(16), 0)
  positive - negative / 6
}

# One unit-peak regressor for a boxcar trial. Unit-peak scaling keeps the
# beta interpretable as a pattern amplitude, so the study-duration covariate
# modulates later reinstatement rather than the encoding beta's scale.
trial_regressor <- function(onset, duration, times, micro = 0.5) {
  offsets <- seq(0, duration, by = micro)
  value <- rowSums(vapply(offsets, function(o) {
    hrf_canonical(times - (onset + o))
  }, numeric(length(times))))
  peak <- max(value)
  if (peak > 0) value <- value / peak
  value
}

legendre_drift <- function(n) {
  x <- seq(-1, 1, length.out = n)
  cbind(drift1 = x, drift2 = (3 * x^2 - 1) / 2)
}

#' Per-run design matrices, trial orders, and effect target matrices.
#'
#' The design is fixed across simulated subjects: subjects differ only in
#' their noise realization, so the across-subject spread reported later is
#' pure estimation variance rather than design variability.
er_design <- function(items, probes) {
  set.seed(SEED + 2L)
  encoding_names <- paste0("enc_", items$item)
  retrieval_names <- probes$probe
  times <- (seq_len(N_TR) - 1L) * TR
  drift <- legendre_drift(N_TR)

  runs <- lapply(seq_len(N_RUNS), function(run) {
    encoding_order <- sample(seq_len(nrow(items)))
    retrieval_order <- sample(seq_len(nrow(probes)))

    onset <- 12
    encoding_onsets <- numeric(nrow(items))
    for (k in seq_along(encoding_order)) {
      trial <- encoding_order[[k]]
      encoding_onsets[[trial]] <- onset
      onset <- onset + items$study_duration[[trial]] + sample(ISI_JITTER, 1L)
    }
    onset <- onset + 25
    retrieval_onsets <- numeric(nrow(probes))
    for (k in seq_along(retrieval_order)) {
      trial <- retrieval_order[[k]]
      retrieval_onsets[[trial]] <- onset
      onset <- onset + RETRIEVAL_DURATION + sample(ISI_JITTER, 1L)
    }
    stopifnot(onset + 20 < N_TR * TR)

    encoding_columns <- vapply(seq_len(nrow(items)), function(trial) {
      trial_regressor(encoding_onsets[[trial]], items$study_duration[[trial]],
        times)
    }, numeric(N_TR))
    retrieval_columns <- vapply(seq_len(nrow(probes)), function(trial) {
      trial_regressor(retrieval_onsets[[trial]], RETRIEVAL_DURATION, times)
    }, numeric(N_TR))
    motion <- cumsum(rnorm(N_TR, sd = 0.35))
    motion <- (motion - mean(motion)) / stats::sd(motion)

    X <- cbind(encoding_columns, retrieval_columns,
      `(baseline)` = 1, drift, motion = motion)
    colnames(X) <- c(encoding_names, retrieval_names,
      "(baseline)", "drift1", "drift2", "motion")
    list(
      design = X,
      last_onset = onset,
      encoding_onsets = encoding_onsets,
      retrieval_onsets = retrieval_onsets
    )
  })
  names(runs) <- paste0("run", seq_len(N_RUNS))

  n_nuisance <- 4L
  n_encoding <- length(encoding_names)
  n_retrieval <- length(retrieval_names)
  n_columns <- n_encoding + n_retrieval + n_nuisance
  column_names <- colnames(runs[[1L]]$design)

  encoding_target <- matrix(0, n_encoding, n_columns,
    dimnames = list(encoding_names, column_names))
  encoding_target[cbind(seq_len(n_encoding), seq_len(n_encoding))] <- 1
  retrieval_target <- matrix(0, n_retrieval, n_columns,
    dimnames = list(retrieval_names, column_names))
  retrieval_target[cbind(seq_len(n_retrieval),
    n_encoding + seq_len(n_retrieval))] <- 1

  list(
    runs = runs,
    encoding_names = encoding_names,
    retrieval_names = retrieval_names,
    encoding_target = encoding_target,
    retrieval_target = retrieval_target,
    residual_df = N_TR - n_columns
  )
}

## ---- Planted truth ------------------------------------------------------
region_index <- function(region) which(REGION_LABELS == region)

#' Noiseless per-run encoding and retrieval patterns.
#'
#' The same-run item state `m[k, i, ]` is added to item i's encoding pattern
#' and to that item's retrieval pattern *within the same run*. It stands in
#' for anything that couples an item's two trials inside one scanning run.
#' Cross-run pairing removes it in expectation; same-run pairing does not.
er_truth <- function(items, probes) {
  set.seed(SEED + 3L)
  A <- region_index("regionA")
  B <- region_index("regionB")
  n_items <- nrow(items)
  n_probes <- nrow(probes)

  item_pattern <- matrix(0, n_items, N_VOXELS,
    dimnames = list(items$item, NULL))
  item_pattern[, A] <- matrix(rnorm(n_items * length(A)), n_items, length(A))

  category_A <- matrix(rnorm(length(CATEGORIES) * length(A)),
    length(CATEGORIES), length(A), dimnames = list(CATEGORIES, NULL))
  category_B <- matrix(rnorm(length(CATEGORIES) * length(B)),
    length(CATEGORIES), length(B), dimnames = list(CATEGORIES, NULL))

  lure_pattern <- matrix(0, n_probes, N_VOXELS,
    dimnames = list(probes$probe, NULL))
  lure_rows <- which(!probes$old)
  lure_pattern[lure_rows, A] <- matrix(
    rnorm(length(lure_rows) * length(A)), length(lure_rows), length(A))

  duration <- setNames(items$study_duration, items$item)
  matched_items <- probes$source_item[probes$old]
  duration_center <- mean(duration[matched_items])
  gain <- REINSTATEMENT_INTERCEPT +
    REINSTATEMENT_SLOPE * (duration - duration_center)
  stopifnot(all(gain > 0))

  base_encoding <- matrix(0, n_items, N_VOXELS,
    dimnames = list(paste0("enc_", items$item), NULL))
  base_encoding[, A] <- AMP_ITEM_A * item_pattern[, A] +
    AMP_CATEGORY_A * category_A[items$category, , drop = FALSE]
  base_encoding[, B] <- AMP_CATEGORY_B * category_B[items$category, ,
    drop = FALSE]

  base_retrieval <- matrix(0, n_probes, N_VOXELS,
    dimnames = list(probes$probe, NULL))
  base_retrieval[, A] <- AMP_CATEGORY_A * category_A[probes$category, ,
    drop = FALSE] + AMP_LURE_A * lure_pattern[, A]
  old_rows <- which(probes$old)
  base_retrieval[old_rows, A] <- base_retrieval[old_rows, A] +
    AMP_ITEM_A * gain[probes$source_item[old_rows]] *
      item_pattern[probes$source_item[old_rows], A]
  base_retrieval[, B] <- AMP_CATEGORY_B * category_B[probes$category, ,
    drop = FALSE]

  state <- array(rnorm(N_RUNS * n_items * N_VOXELS, sd = STATE_SD),
    dim = c(N_RUNS, n_items, N_VOXELS))
  lure_state <- array(rnorm(N_RUNS * length(lure_rows) * N_VOXELS,
    sd = STATE_SD), dim = c(N_RUNS, length(lure_rows), N_VOXELS))

  encoding <- lapply(seq_len(N_RUNS), function(run) {
    base_encoding + state[run, , ]
  })
  retrieval <- lapply(seq_len(N_RUNS), function(run) {
    value <- base_retrieval
    value[old_rows, ] <- value[old_rows, ] +
      state[run, match(probes$source_item[old_rows], items$item), ]
    value[lure_rows, ] <- value[lure_rows, ] + lure_state[run, , ]
    value
  })
  names(encoding) <- names(retrieval) <- paste0("run", seq_len(N_RUNS))

  list(
    encoding = encoding, retrieval = retrieval,
    base_encoding = base_encoding, base_retrieval = base_retrieval,
    item_pattern = item_pattern, gain = gain,
    duration_center = duration_center
  )
}

#' Spatial noise factor: AR-like correlation inside a region, independent
#' across regions.
er_noise_factor <- function() {
  covariance <- matrix(0, N_VOXELS, N_VOXELS)
  for (region in names(REGION_SIZES)) {
    index <- region_index(region)
    covariance[index, index] <- toeplitz(NOISE_AR^(seq_along(index) - 1L))
  }
  chol(NOISE_SD^2 * covariance + diag(1e-8, N_VOXELS))
}

#' One subject: per-run response matrices, from the fixed design and truth.
er_simulate_subject <- function(design, truth, subject, factor) {
  responses <- lapply(seq_len(N_RUNS), function(run) {
    set.seed(SEED + 1000L * subject + run)
    X <- design$runs[[run]]$design
    effects <- rbind(truth$encoding[[run]], truth$retrieval[[run]])
    n_effects <- nrow(effects)
    nuisance <- rbind(
      rnorm(N_VOXELS, mean = BASELINE_MEAN, sd = 5),
      rnorm(N_VOXELS, sd = DRIFT_SD),
      rnorm(N_VOXELS, sd = DRIFT_SD),
      rnorm(N_VOXELS, sd = MOTION_SD)
    )
    signal <- X[, seq_len(n_effects), drop = FALSE] %*% effects +
      X[, n_effects + seq_len(4L), drop = FALSE] %*% nuisance
    noise <- matrix(rnorm(N_TR * N_VOXELS), N_TR, N_VOXELS) %*% factor
    signal + noise
  })
  names(responses) <- names(design$runs)
  responses
}

## ---- Closed-form planted geometry ---------------------------------------
#' The pair value a plan estimates, computed from the planted truth.
#'
#' A rectangular plan contracts, per measurement m with frame weights w,
#'   G[i, j] = sum_edges weight_e * sum_v w_v * B_enc[k(e)][i, v] * B_ret[l(e)][j, v]
#' Cross-run edges are the ordered pairs k != l with uniform weight; the
#' same-run comparator is the diagonal k == l. Both are computed exactly here
#' from the noiseless per-run patterns, with no Monte Carlo.
er_planted_geometry <- function(truth, weights, edges = c("cross", "same")) {
  edges <- match.arg(edges)
  same <- Reduce(`+`, lapply(seq_len(N_RUNS), function(run) {
    truth$encoding[[run]] %*% (weights * t(truth$retrieval[[run]]))
  }))
  if (edges == "same") return(same / N_RUNS)
  total_encoding <- Reduce(`+`, truth$encoding)
  total_retrieval <- Reduce(`+`, truth$retrieval)
  outer_all <- total_encoding %*% (weights * t(total_retrieval))
  (outer_all - same) / (N_RUNS * (N_RUNS - 1L))
}

#' The planted value of a compiled pair query at one measurement.
er_planted_value <- function(query, planted) {
  H <- as.matrix(query$operator)
  stopifnot(identical(dim(H), dim(planted)))
  sum(H * planted)
}
