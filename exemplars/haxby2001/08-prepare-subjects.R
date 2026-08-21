#!/usr/bin/env Rscript
# 08-prepare-subjects.R -- run 01's preparation over every downloaded subject.
#
# `01-prepare-data.R` prepares subject 1 alone. This script calls the factored
# recipe in `subjects.R` for each subject `00-download.R` fetched, writing one
# `data/prepared-<subj>.rds` per subject. It does not touch
# `data/prepared-smoke.rds`, which scripts 02-07 read.
#
# The first thing it does is re-derive subject 1 through the factored path and
# compare it, element by element, against the committed `prepared-smoke.rds`.
# If the two disagree, the six-subject panel would be a different analysis from
# the one the rest of the exemplar reports, and the script stops. This check is
# the reason to run 01 first even though 08 does not otherwise need it.
#
# Idempotent: a subject whose prepared object exists is skipped (delete the
# file, or set OVERWRITE=1, to rebuild). Runtime is dominated by reading each
# 300 MB bold.nii.gz -- roughly half a minute per subject.
#
# Environment: SUBJECTS, OVERWRITE, EXEMPLAR_DIR.
#
# Output: data/prepared-subj1.rds ... data/prepared-subj6.rds (gitignored),
#         data/prepared-subjects-summary.rds.

suppressMessages(library(neuroim2))

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
source(file.path(exemplar_dir, "models.R"))
source(file.path(exemplar_dir, "subjects.R"))
paths <- exemplar_paths(exemplar_dir)

requested <- Sys.getenv("SUBJECTS")
requested <- if (nzchar(requested)) {
  trimws(strsplit(requested, "[,[:space:]]+")[[1L]])
} else {
  SUBJECT_IDS
}
requested <- requested[nzchar(requested)]
overwrite <- identical(Sys.getenv("OVERWRITE"), "1")

subjects <- available_subjects(exemplar_dir, requested)
absent <- setdiff(requested, subjects)
if (length(absent)) {
  message("Not downloaded, skipping: ", paste(absent, collapse = ", "),
          "\n  (run 00-download.R; it reports which subjects failed and why)")
}
if (!length(subjects)) {
  stop("No subject directories found under ", paths$data,
       "\nRun 00-download.R first.")
}
message("Preparing: ", paste(subjects, collapse = ", "))

script_t0 <- Sys.time()

## ---- Subject 1 against the committed preparation -------------------------
# Only meaningful if subject 1 is in scope and 01 has already been run.
if ("subj1" %in% subjects && file.exists(paths$prepared)) {
  message("\n== cross-check: factored recipe vs 01-prepare-data.R ==")
  re <- prepare_subject(exemplar_dir, "subj1", overwrite = overwrite)
  old <- readRDS(paths$prepared)
  if (!identical(as.integer(re$vt_index), as.integer(old$vt_index))) {
    stop("Factored preparation selects a different VT voxel set for subject 1 ",
         "than data/prepared-smoke.rds does.")
  }
  if (!identical(names(re$runs), names(old$runs))) {
    stop("Factored preparation produces different run names for subject 1.")
  }
  deltas <- vapply(names(re$runs), function(nm)
    max(abs(re$runs[[nm]] - old$runs[[nm]])), numeric(1))
  worst <- max(deltas)
  if (worst > 0) {
    stop("Factored preparation does not reproduce data/prepared-smoke.rds ",
         "for subject 1 (max abs difference ", format(worst, digits = 3), ").")
  }
  message("  subject 1 reproduced exactly (max abs difference ", worst,
          " over ", length(deltas), " runs x ", nrow(re$runs[[1L]]),
          " conditions x ", ncol(re$runs[[1L]]), " voxels)")
  cross_check <- list(checked = TRUE, max_abs_difference = worst,
                      runs = length(deltas))
} else {
  message("\n(no cross-check: ", if (!("subj1" %in% subjects))
    "subject 1 not in scope" else
    "data/prepared-smoke.rds absent -- run 01-prepare-data.R", ")")
  cross_check <- list(checked = FALSE, max_abs_difference = NA_real_,
                      runs = NA_integer_)
}

## ---- Every subject -------------------------------------------------------
message("\n== preparing ==")
prepared <- list()
failures <- list()
for (subj in subjects) {
  out <- tryCatch(prepare_subject(exemplar_dir, subj, overwrite = overwrite),
                  error = function(e) e)
  if (inherits(out, "error")) {
    message("  ", subj, ": FAILED -- ", conditionMessage(out))
    failures[[subj]] <- conditionMessage(out)
  } else {
    prepared[[subj]] <- out
  }
}

## ---- Summary -------------------------------------------------------------
summary_rows <- lapply(names(prepared), function(subj) {
  p <- prepared[[subj]]
  data.frame(
    subject = subj, n_runs = p$n_runs, n_conditions = nrow(p$runs[[1L]]),
    n_vt_full = p$n_vt_full, n_vt_used = p$n_vt_used,
    n_vt_dropped = p$n_vt_full - p$n_vt_used,
    volumes_per_condition_min = min(p$block_counts, na.rm = TRUE),
    volumes_per_condition_max = max(p$block_counts, na.rm = TRUE),
    quirks = if (length(p$quirks)) paste(p$quirks, collapse = "; ") else "",
    stringsAsFactors = FALSE
  )
})
summary_df <- if (length(summary_rows)) do.call(rbind, summary_rows) else NULL

script_seconds <- as.numeric(difftime(Sys.time(), script_t0, units = "secs"))
saveRDS(list(summary = summary_df, cross_check = cross_check,
             failures = failures, requested = requested,
             absent = absent, seconds = script_seconds, when = Sys.time()),
        file.path(paths$data, "prepared-subjects-summary.rds"))

message("\n== summary ==")
if (!is.null(summary_df)) print(summary_df, row.names = FALSE)
if (length(failures)) {
  message("\nFailed: ", paste(names(failures), collapse = ", "))
}
message("\nTotal: ", round(script_seconds, 1), " s for ", length(prepared),
        " subject(s)")
