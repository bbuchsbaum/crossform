#!/usr/bin/env Rscript
# 01-prepare-data.R -- the single shared preparation for both arms.
#
# Produces one 8 x P condition-mean matrix per run (P = VT mask voxels),
# computed once. Both the crossform arm and the rMVPA arm consume exactly
# this object, so no preprocessing, HRF, or averaging choice can differ
# between them. The comparison isolates pairing/distance/RSA semantics only.
#
# Steps, in order:
#   1. read bold.nii.gz, restrict to the VT mask voxels
#   2. within each run, z-score every voxel's time series (mean 0, sd 1)
#   3. shift the block labels forward by HRF_SHIFT_TR to approximate lag
#   4. average the shifted volumes of each condition block within each run
#
# Output: data/prepared-smoke.rds  (gitignored; regenerate with this script)

suppressMessages(library(neuroim2))

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
source(file.path(exemplar_dir, "models.R"))
paths <- exemplar_paths(exemplar_dir)
dir.create(paths$results, showWarnings = FALSE, recursive = TRUE)

if (file.exists(paths$prepared)) {
  message("Prepared data already present: ", paths$prepared,
          "\nDelete it to rebuild.")
  prep <- readRDS(paths$prepared)
  message("  runs: ", length(prep$runs), "  conditions: ",
          nrow(prep$runs[[1]]), "  voxels: ", ncol(prep$runs[[1]]))
  quit(save = "no", status = 0)
}

## ---- Labels -------------------------------------------------------------
labels <- utils::read.table(file.path(paths$subj, "labels.txt"), header = TRUE,
                            stringsAsFactors = FALSE)
stopifnot(identical(names(labels), c("labels", "chunks")))
runs <- sort(unique(labels$chunks))
stopifnot(length(runs) == N_RUNS)
stopifnot(setequal(setdiff(unique(labels$labels), "rest"), CONDITIONS))

## ---- BOLD, restricted to VT ---------------------------------------------
vt_vol <- neuroim2::read_vol(file.path(paths$subj, "mask4_vt.nii.gz"))
vt_mask <- neuroim2::LogicalNeuroVol(as.array(vt_vol) > 0,
                                     neuroim2::space(vt_vol))
vt_index <- which(as.array(vt_mask))          # stable full-volume indices
message("VT mask voxels: ", length(vt_index))

message("Reading bold.nii.gz (this takes a minute) ...")
bold <- neuroim2::read_vec(file.path(paths$subj, "bold.nii.gz"))
stopifnot(dim(bold)[4] == nrow(labels))

# series() returns time x voxel for the requested full-volume indices.
ts_all <- neuroim2::series(bold, vt_index)
storage.mode(ts_all) <- "double"
stopifnot(nrow(ts_all) == nrow(labels), ncol(ts_all) == length(vt_index))
rm(bold); invisible(gc())

## ---- Per-run z-scoring, HRF shift, condition means ----------------------
# Shifting labels forward by k TR is equivalent to averaging the volumes at
# t+k for each labelled block volume t; done inside the run so no shift ever
# crosses a run boundary.
run_means <- vector("list", length(runs))
names(run_means) <- paste0("run", sprintf("%02d", runs + 1L))
block_counts <- matrix(NA_integer_, length(runs), length(CONDITIONS),
                       dimnames = list(names(run_means), CONDITIONS))

# The error-channel arm (05) needs the responses the means were computed
# from, not just the means. Keep the labelled, shifted, z-scored volumes and
# the indicator design whose OLS solution IS the condition mean, so both arms
# and the uncertainty arm share one preparation rather than two.
response_runs <- vector("list", length(runs))
names(response_runs) <- names(run_means)
response_design <- NULL

for (ri in seq_along(runs)) {
  rows <- which(labels$chunks == runs[ri])
  X <- ts_all[rows, , drop = FALSE]
  lab <- labels$labels[rows]

  # 2. within-run z-score per voxel
  mu <- colMeans(X)
  sdv <- apply(X, 2L, stats::sd)
  if (any(sdv == 0)) {
    stop("Run ", runs[ri], " has ", sum(sdv == 0),
         " constant voxel time series inside the VT mask; ",
         "z-scoring is undefined for them.")
  }
  Z <- sweep(sweep(X, 2L, mu, `-`), 2L, sdv, `/`)

  # 3. HRF shift: volume t is attributed to the label at t - HRF_SHIFT_TR
  shifted <- c(rep(NA_character_, HRF_SHIFT_TR),
               lab[seq_len(length(lab) - HRF_SHIFT_TR)])

  # 4. condition means over shifted, non-rest volumes
  M <- matrix(NA_real_, length(CONDITIONS), ncol(Z),
              dimnames = list(CONDITIONS, NULL))
  for (ci in seq_along(CONDITIONS)) {
    keep <- which(!is.na(shifted) & shifted == CONDITIONS[ci])
    if (!length(keep)) stop("No volumes for ", CONDITIONS[ci], " in run ", runs[ri])
    block_counts[ri, ci] <- length(keep)
    M[ci, ] <- colMeans(Z[keep, , drop = FALSE])
  }
  if (anyNA(M)) stop("NA in condition means for run ", runs[ri])
  run_means[[ri]] <- M

  # Labelled volumes only, ordered by condition, plus the matching 0/1
  # indicator design. With disjoint indicators and no intercept the OLS
  # coefficient for condition c is exactly the mean of its volumes, so this
  # design reproduces M above and additionally leaves residuals behind.
  keep_all <- which(!is.na(shifted) & shifted %in% CONDITIONS)
  keep_all <- keep_all[order(match(shifted[keep_all], CONDITIONS))]
  Y <- Z[keep_all, , drop = FALSE]
  Xd <- stats::model.matrix(~ 0 + factor(shifted[keep_all], levels = CONDITIONS))
  colnames(Xd) <- CONDITIONS
  if (is.null(response_design)) {
    response_design <- Xd
  } else if (!identical(unname(Xd), unname(response_design))) {
    stop("Runs do not share one common indicator design; run ", runs[ri])
  }
  # The two routes must agree exactly, or 05 is describing different effects.
  beta <- qr.coef(qr(Xd), Y)
  if (max(abs(beta - M)) > 1e-10) {
    stop("Indicator-design OLS does not reproduce the condition means in run ",
         runs[ri], " (max abs difference ", max(abs(beta - M)), ")")
  }
  response_runs[[ri]] <- Y
}

message("Volumes averaged per condition per run: ",
        paste(sort(unique(as.vector(block_counts))), collapse = "/"))

## ---- Save ---------------------------------------------------------------
prepared <- list(
  runs = run_means,                 # named list of 8 x P matrices, run order
  response_runs = response_runs,    # labelled volumes x P, for the error channel
  response_design = response_design,# shared indicator design, volumes x 8
  conditions = CONDITIONS,
  run_names = names(run_means),
  vt_index = vt_index,              # full-volume indices of the P columns
  mask_dim = dim(vt_mask),
  mask_spacing = neuroim2::spacing(vt_mask),
  mask_file = file.path(paths$subj, "mask4_vt.nii.gz"),
  block_counts = block_counts,
  prep = list(tr_seconds = TR_SECONDS, hrf_shift_tr = HRF_SHIFT_TR,
              zscore = "within-run, per-voxel, mean 0 sd 1",
              built = Sys.time())
)
saveRDS(prepared, paths$prepared)
message("Wrote ", paths$prepared)
message("  ", length(run_means), " runs x ", nrow(run_means[[1]]),
        " conditions x ", ncol(run_means[[1]]), " voxels")
