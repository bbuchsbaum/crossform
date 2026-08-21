## subjects.R -- the per-subject counterpart of models.R.
##
## `01-prepare-data.R` prepares subject 1 and writes `data/prepared-smoke.rds`,
## which scripts 02-07 consume. Ticket E10 needs the same preparation run over
## all six Haxby subjects, so the preparation is factored here as one function
## that scripts 08 and 09 call in a loop. It is deliberately the SAME recipe:
## VT mask, within-run per-voxel z-score, 2 TR haemodynamic shift, condition
## means per run. Any divergence between this and 01 would make the six-subject
## panel a different analysis from the one the rest of the exemplar reports, so
## `08-prepare-subjects.R` re-derives subject 1 here and checks it against the
## committed `prepared-smoke.rds` rather than assuming agreement.
##
## Sourced after `models.R`, whose CONDITIONS / HRF_SHIFT_TR / TR_SECONDS this
## file uses. Requires neuroim2.

SUBJECT_IDS <- paste0("subj", 1:6)

#' Paths for one subject.
#'
#' Mirrors `exemplar_paths()` but parameterized. Subject 1's prepared object
#' here is `prepared-subj1.rds`, NOT the pre-existing `prepared-smoke.rds`:
#' the older file is left exactly as it is so scripts 02-07 keep reading the
#' artifact they were run against.
subject_paths <- function(exemplar_dir, subject) {
  stopifnot(length(subject) == 1L, subject %in% SUBJECT_IDS)
  list(
    subject  = subject,
    data     = file.path(exemplar_dir, "data"),
    subj     = file.path(exemplar_dir, "data", subject),
    results  = file.path(exemplar_dir, "results"),
    prepared = file.path(exemplar_dir, "data",
                         sprintf("prepared-%s.rds", subject)),
    mask     = file.path(exemplar_dir, "data", subject, "mask4_vt.nii.gz")
  )
}

#' Which subjects have a usable raw directory on disk.
#'
#' The six-subject scripts run over whatever `00-download.R` actually managed
#' to fetch and never over a subject that is absent. Missing subjects are
#' reported by the caller, not silently skipped.
available_subjects <- function(exemplar_dir, subjects = SUBJECT_IDS) {
  required <- c("bold.nii.gz", "mask4_vt.nii.gz", "labels.txt")
  keep <- vapply(subjects, function(s) {
    d <- file.path(exemplar_dir, "data", s)
    dir.exists(d) && all(file.exists(file.path(d, required)))
  }, logical(1))
  subjects[keep]
}

#' The VT mask as a logical volume, minus any voxels the preparation dropped.
#'
#' Both the preparation and the analysis must see exactly the same feature set,
#' or the domain the frame is compiled on will not match the columns of the
#' pattern matrices. The mask file is the single source of truth; `drop` is the
#' subject-specific exclusion the preparation recorded.
vt_mask_volume <- function(mask_file, drop = integer(0)) {
  vol <- neuroim2::read_vol(mask_file)
  arr <- as.array(vol) > 0
  if (length(drop)) arr[drop] <- FALSE
  neuroim2::LogicalNeuroVol(arr, neuroim2::space(vol))
}

#' Prepare one subject: per-run condition means over the VT mask.
#'
#' Identical in recipe to `01-prepare-data.R`, with three differences forced by
#' running over six subjects rather than one:
#'
#'   * the run count is READ, not asserted to be 12. Distributions of this
#'     dataset differ in whether every subject has all twelve runs; the number
#'     found is recorded per subject instead of being a precondition.
#'   * a condition absent from a run is a recorded quirk that drops the run,
#'     not an error. A run cannot contribute a condition-mean matrix it has no
#'     volumes for, and dropping the run is the honest handling; the remaining
#'     runs still pair.
#'   * a voxel whose time series is constant in ANY run is DROPPED from the
#'     mask for that subject, and the count is recorded. `01` stops on this
#'     because subject 1 has none. Z-scoring is undefined for such a voxel, and
#'     the alternative -- keeping it with an arbitrary scale -- would put a
#'     fabricated pattern into the ledger.
#'
#' `response_runs` / `response_design` (which only `05` consumes) are not built
#' here; the six-subject arm reads condition means only.
#'
#' Returns the prepared list invisibly and writes it to `paths$prepared`.
prepare_subject <- function(exemplar_dir, subject, overwrite = FALSE,
                            verbose = TRUE) {
  paths <- subject_paths(exemplar_dir, subject)
  say <- function(...) if (verbose) message(...)

  if (file.exists(paths$prepared) && !overwrite) {
    prep <- readRDS(paths$prepared)
    say("  ", subject, ": prepared object already present (",
        length(prep$runs), " runs x ", nrow(prep$runs[[1L]]),
        " conditions x ", ncol(prep$runs[[1L]]), " voxels)")
    return(invisible(prep))
  }

  t0 <- Sys.time()
  quirks <- character(0)

  ## ---- Labels ------------------------------------------------------------
  labels <- utils::read.table(file.path(paths$subj, "labels.txt"),
                              header = TRUE, stringsAsFactors = FALSE)
  stopifnot(identical(names(labels), c("labels", "chunks")))
  runs <- sort(unique(labels$chunks))
  n_runs_found <- length(runs)
  if (n_runs_found != N_RUNS) {
    quirks <- c(quirks, sprintf("%d runs in labels.txt, not the %d subject 1 has",
                                n_runs_found, N_RUNS))
  }
  missing_conditions <- setdiff(CONDITIONS, unique(labels$labels))
  if (length(missing_conditions)) {
    stop(subject, ": labels.txt has no volumes at all for ",
         paste(missing_conditions, collapse = ", "),
         "; the eight-condition basis is not available for this subject.")
  }

  ## ---- BOLD, restricted to VT --------------------------------------------
  vt_vol <- neuroim2::read_vol(paths$mask)
  vt_mask <- neuroim2::LogicalNeuroVol(as.array(vt_vol) > 0,
                                       neuroim2::space(vt_vol))
  vt_index_full <- which(as.array(vt_mask))
  say("  ", subject, ": VT mask voxels ", length(vt_index_full),
      ", runs ", n_runs_found, ", volumes ", nrow(labels))

  bold <- neuroim2::read_vec(file.path(paths$subj, "bold.nii.gz"))
  if (dim(bold)[4] != nrow(labels)) {
    stop(subject, ": bold.nii.gz has ", dim(bold)[4],
         " volumes but labels.txt has ", nrow(labels), " rows.")
  }
  ts_all <- neuroim2::series(bold, vt_index_full)
  storage.mode(ts_all) <- "double"
  stopifnot(nrow(ts_all) == nrow(labels),
            ncol(ts_all) == length(vt_index_full))
  rm(bold); invisible(gc())

  ## ---- Degenerate voxels, found before any run is z-scored ---------------
  # A voxel is usable only if it varies in EVERY run, because it must carry a
  # comparable pattern in every run the pairing uses.
  sd_by_run <- vapply(runs, function(r) {
    apply(ts_all[labels$chunks == r, , drop = FALSE], 2L, stats::sd)
  }, numeric(ncol(ts_all)))
  usable <- apply(sd_by_run, 1L, function(s) all(is.finite(s) & s > 0))
  n_dropped <- sum(!usable)
  dropped_index <- vt_index_full[!usable]
  if (n_dropped) {
    quirks <- c(quirks, sprintf(
      "%d of %d VT voxels constant in at least one run and dropped (z-scoring undefined)",
      n_dropped, length(vt_index_full)))
    say("    dropped ", n_dropped, " constant VT voxel(s)")
  }
  vt_index <- vt_index_full[usable]
  ts_all <- ts_all[, usable, drop = FALSE]
  if (!ncol(ts_all)) stop(subject, ": no usable VT voxels remain.")

  ## ---- Per-run z-scoring, HRF shift, condition means ---------------------
  run_means <- vector("list", n_runs_found)
  run_names <- paste0("run", sprintf("%02d", runs + 1L))
  names(run_means) <- run_names
  block_counts <- matrix(NA_integer_, n_runs_found, length(CONDITIONS),
                         dimnames = list(run_names, CONDITIONS))
  dropped_runs <- character(0)
  # A logical keep vector, NOT `run_means[[ri]] <- NULL`: assigning NULL to a
  # list element DELETES it, which would shift every later run into the wrong
  # name. The matrices stay where they are and the vector says which survive.
  keep_runs <- rep(TRUE, n_runs_found)

  for (ri in seq_along(runs)) {
    rows <- which(labels$chunks == runs[ri])
    X <- ts_all[rows, , drop = FALSE]
    lab <- labels$labels[rows]

    mu <- colMeans(X)
    sdv <- apply(X, 2L, stats::sd)
    Z <- sweep(sweep(X, 2L, mu, `-`), 2L, sdv, `/`)

    shifted <- c(rep(NA_character_, HRF_SHIFT_TR),
                 lab[seq_len(length(lab) - HRF_SHIFT_TR)])

    M <- matrix(NA_real_, length(CONDITIONS), ncol(Z),
                dimnames = list(CONDITIONS, NULL))
    empty <- character(0)
    for (ci in seq_along(CONDITIONS)) {
      keep <- which(!is.na(shifted) & shifted == CONDITIONS[ci])
      if (!length(keep)) { empty <- c(empty, CONDITIONS[ci]); next }
      block_counts[ri, ci] <- length(keep)
      M[ci, ] <- colMeans(Z[keep, , drop = FALSE])
    }
    if (length(empty)) {
      dropped_runs <- c(dropped_runs, run_names[ri])
      quirks <- c(quirks, sprintf(
        "run %s carries no labelled volumes (%s absent) and is dropped",
        run_names[ri],
        if (length(empty) == length(CONDITIONS)) "all eight conditions"
        else paste(empty, collapse = "/")))
      keep_runs[ri] <- FALSE
      next
    }
    if (anyNA(M)) stop(subject, ": NA in condition means for run ", runs[ri])
    run_means[[ri]] <- M
  }

  run_means <- run_means[keep_runs]
  block_counts <- block_counts[keep_runs, , drop = FALSE]
  if (length(run_means) < 2L) {
    stop(subject, ": fewer than two usable runs; nothing to cross-pair.")
  }

  prepared <- list(
    subject = subject,
    runs = run_means,
    conditions = CONDITIONS,
    run_names = names(run_means),
    n_runs = length(run_means),
    dropped_runs = dropped_runs,
    vt_index = vt_index,
    vt_index_full = vt_index_full,
    dropped_index = dropped_index,
    n_vt_full = length(vt_index_full),
    n_vt_used = length(vt_index),
    mask_dim = dim(vt_mask),
    mask_spacing = neuroim2::spacing(vt_mask),
    mask_file = paths$mask,
    block_counts = block_counts,
    quirks = quirks,
    prep = list(tr_seconds = TR_SECONDS, hrf_shift_tr = HRF_SHIFT_TR,
                zscore = "within-run, per-voxel, mean 0 sd 1",
                recipe = "identical to 01-prepare-data.R",
                seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
                built = Sys.time())
  )
  saveRDS(prepared, paths$prepared)
  say("    wrote ", basename(paths$prepared), " -- ", length(run_means),
      " runs x ", nrow(run_means[[1L]]), " conditions x ",
      ncol(run_means[[1L]]), " voxels (",
      round(prepared$prep$seconds), " s)")
  invisible(prepared)
}
