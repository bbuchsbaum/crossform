#!/usr/bin/env Rscript
# 00-common.R -- shared constants, paths and helpers for WS-E slice 2.
#
# Sourced by every numbered script. Defines nothing that touches the network
# and nothing that writes: it is safe to source in an interactive session to
# inspect the settings a run used.
#
# The dataset is OpenNeuro ds003745 snapshot 2.1.1 (CC0), fMRIPrep 21.0.2,
# space MNI152NLin2009cAsym. See DECISION.md for why, and for the risk list
# that the constants below are answers to.

suppressMessages({
  library(neuroim2)
  library(Matrix)
})

SLICE_DIR <- if (nzchar(Sys.getenv("SLICE2_DIR"))) {
  Sys.getenv("SLICE2_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])))
}
REPO_DIR    <- normalizePath(file.path(SLICE_DIR, "..", ".."))
DATA_DIR    <- file.path(SLICE_DIR, "data")
FMRIPREP    <- file.path(DATA_DIR, "derivatives", "fmriprep")
DERIVED_DIR <- file.path(DATA_DIR, "derived")     # git-ignored intermediates
RESULTS_DIR <- file.path(SLICE_DIR, "results")    # committed receipts

load_crossform <- function() {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required to load crossform from source")
  }
  suppressMessages(pkgload::load_all(REPO_DIR, quiet = TRUE, export_all = FALSE))
  invisible(as.character(utils::packageVersion("crossform")))
}

## ---- The cohort -----------------------------------------------------------
# E12's deliberate 6 + 6 age split (DECISION.md section 2). n = 6 per group
# supports description, never inference: no younger-vs-older contrast is
# reported anywhere in this slice, and the age column exists only so that the
# coverage and sink numbers can be read against it.
SUBJECTS_YOUNGER <- c("sub-104", "sub-105", "sub-107", "sub-108",
                      "sub-112", "sub-113")
SUBJECTS_OLDER   <- c("sub-111", "sub-127", "sub-128", "sub-129",
                      "sub-130", "sub-131")
SUBJECTS <- c(SUBJECTS_YOUNGER, SUBJECTS_OLDER)

TRUST_RUNS  <- 1:5
SHARED_RUNS <- 1:2
TR <- 2.02                 # from task-trust_bold.json / task-sharedreward_bold.json

## ---- The conditions -------------------------------------------------------
# task `trust`: 9 levels, all present in all 60 chosen runs (DECISION.md 3b).
# `missed_trial` is present in 25 of the 60 runs and is a nuisance regressor of
# no interest -- never a tenth condition (risk 4).
TRUST_CONDITIONS <- c(
  "choice_computer", "choice_friend", "choice_stranger",
  "outcome_computer_defect", "outcome_computer_recip",
  "outcome_friend_defect",   "outcome_friend_recip",
  "outcome_stranger_defect", "outcome_stranger_recip"
)
# Trial types that are modelled as regressors of no interest wherever they
# appear, in either task. Every OTHER level present in an events file is also
# modelled (01-glm.R) -- this constant only names the one that would otherwise
# be mistaken for a condition.
NUISANCE_TRIAL_TYPES <- "missed_trial"

# Stamped into every fitted-beta file and checked before a cached fit is
# reused. Bump it whenever 01-glm.R's design changes, so that a stale cache
# refits instead of silently surviving the change.
#   psc-allevents-mot12-cosine-fd-spike:v1
#     percent signal change; every trial type in the events file modelled;
#     6 motion + 6 derivatives; the full fMRIPrep cosine basis;
#     framewise_displacement; one indicator per volume with FD > 0.5 mm.
GLM_MODEL_ID <- "psc-allevents-mot12-cosine-fd-spike:v1"

# task `sharedreward`: only the 7 event-level conditions that appear in ALL 24
# runs (DECISION.md risk 7). `event_computer_neutral` is missing from 4 runs
# and `event_stranger_neutral` from 1, so neither is used; the six block-level
# levels are not used either. This task is the P^F fitting source and nothing
# else -- no `trust` number is ever computed from it.
SHARED_CONDITIONS <- c(
  "event_computer_punish", "event_computer_reward",
  "event_friend_punish",   "event_friend_reward",
  "event_stranger_punish", "event_stranger_reward",
  "event_friend_neutral"
)

## ---- The query bank -------------------------------------------------------
# Five contrasts read at every group node in one pass. Every row is zero-sum,
# which is what the distance basis requires of a contrast and what keeps the
# whole bank inside the span a pairwise-difference readout can reach.
#
# `friend-vs-computer-choice` is a bare pairwise difference, i.e. exactly the
# RDM edge between two conditions. It is in the bank so that the commutation
# acceptance (03, identity 6) has an edge to compare against; the other four
# are the scientific readouts.
.named_contrast <- function(positive, negative) {
  w <- stats::setNames(numeric(length(TRUST_CONDITIONS)), TRUST_CONDITIONS)
  w[positive] <-  1 / length(positive)
  w[negative] <- -1 / length(negative)
  w
}
.CHOICE  <- grep("^choice_",  TRUST_CONDITIONS, value = TRUE)
.OUTCOME <- grep("^outcome_", TRUST_CONDITIONS, value = TRUE)
.RECIP   <- grep("_recip$",   TRUST_CONDITIONS, value = TRUE)
.DEFECT  <- grep("_defect$",  TRUST_CONDITIONS, value = TRUE)
.FRIEND  <- grep("friend",    TRUST_CONDITIONS, value = TRUE)
.STRANGER<- grep("stranger",  TRUST_CONDITIONS, value = TRUE)

QUERY_BANK <- rbind(
  `choice-vs-outcome`         = .named_contrast(.CHOICE, .OUTCOME),
  `social-vs-computer-choice` = .named_contrast(
                                  c("choice_friend", "choice_stranger"),
                                  "choice_computer"),
  `recip-vs-defect`           = .named_contrast(.RECIP, .DEFECT),
  `friend-vs-stranger`        = .named_contrast(.FRIEND, .STRANGER),
  `friend-vs-computer-choice` = .named_contrast("choice_friend",
                                                "choice_computer")
)
stopifnot(max(abs(rowSums(QUERY_BANK))) < 1e-12)
HEADLINE_QUERY  <- "choice-vs-outcome"
# The commutation acceptance needs a query that is literally an RDM edge.
COMMUTATION_QUERY <- "friend-vs-computer-choice"
COMMUTATION_EDGE  <- c("choice_friend", "choice_computer")

## ---- Geometry settings ----------------------------------------------------
# Searchlight radius in MILLIMETRES. The voxels are 2.973 x 2.973 x 3.22 mm --
# anisotropic (risk 3) -- so a radius in voxels would build spheres that are
# squashed along z. `compile_frame()` measures in the domain's own coordinates,
# which `neuroim2_volume_domain()` reports in mm, so a millimetre radius here
# is a millimetre radius in the brain.
SEARCHLIGHT_RADIUS_MM <- 8

# The group grid: every GRID_STEP-th voxel of the shared lattice in each
# direction, i.e. 8.92 x 8.92 x 9.66 mm spacing.
GRID_STEP <- 3L
# Grey-matter gate on the group nodes, applied to the across-subject mean of
# the fMRIPrep `label-GM_probseg` resampled to the functional grid.
GM_THRESHOLD <- 0.25

# P^A's assignment radius in mm. Half the group-grid cell diagonal is
# sqrt(4.46^2 + 4.46^2 + 4.83^2) = 7.95 mm, so every voxel inside the group
# grid's own territory has a node within 8 mm; 9 mm leaves one voxel of slack
# and still sinks everything that is genuinely outside the consensus grey
# matter. THIS RADIUS DEFINES THE SINK for both transports: a native row is
# placed if and only if P^A places it, and P^F copies P^A's sink mass exactly.
TRANSPORT_RADIUS_MM <- 9

# P^F's redistribution radius in mm, deliberately wider. At 9 mm the mean
# support is only 1.8 group nodes -- the grid is grey-matter gated and
# consensus restricted, so it is sparser than its nominal spacing suggests --
# and a transport that can choose between fewer than two destinations cannot
# express a functional preference. 12 mm gives P^F somewhere to move WITHOUT
# giving it anything to discard, because the sink is still P^A's.
#
# P^F therefore has strictly more freedom than P^A, and that asymmetry is
# exactly what the section 7.3 null band controls for: a permuted-fingerprint
# P^F has this same wider support and this same spreading, so anything eta
# shows above the null band is attributable to the fingerprint correspondence
# and not to the extra room.
PF_SUPPORT_RADIUS_MM <- 12

# Softmax temperature for P^F. Fixed rather than tuned: tuning it on the
# held-out `trust` data would be exactly the circularity section 7.2 measures
# at 3.15x the honest gain.
PF_TEMPERATURE <- 6

# The eta null band (contract section 7.3 / 14.5).
ETA_NULL_DRAWS <- 200L
ETA_NULL_SEED  <- 20260820L

# Group-node subject-coverage floor (contract section 7.5 / 14.3, an open
# maintainer decision). Declared here, reported, never used to silently drop a
# node.
COVERAGE_FLOOR <- length(SUBJECTS)

## ---- Tolerances -----------------------------------------------------------
POPULATION_TOLERANCE <- 1e-12   # contract section 11
FRAME_TOLERANCE      <- 1e-10   # conservative-frame column mass

## ---- Filename helpers -----------------------------------------------------
# Raw events are zero-padded (`run-01`), fMRIPrep derivatives are not
# (`run-1`). This is DECISION.md risk 1, "the single most likely source of a
# silent mispairing", so the two spellings get two functions and neither is
# ever built by pasting a run number into a format string at the call site.
raw_run  <- function(r) sprintf("run-%02d", as.integer(r))
deriv_run<- function(r) sprintf("run-%d",   as.integer(r))

path_events <- function(subject, task, run) {
  file.path(DATA_DIR, subject, "func",
            sprintf("%s_task-%s_%s_events.tsv", subject, task, raw_run(run)))
}
path_bold <- function(subject, task, run) {
  file.path(FMRIPREP, subject, "func",
    sprintf("%s_task-%s_%s_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz",
            subject, task, deriv_run(run)))
}
path_boldmask <- function(subject, task, run) {
  file.path(FMRIPREP, subject, "func",
    sprintf("%s_task-%s_%s_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz",
            subject, task, deriv_run(run)))
}
path_confounds <- function(subject, task, run) {
  file.path(FMRIPREP, subject, "func",
            sprintf("%s_task-%s_%s_desc-confounds_timeseries.tsv",
                    subject, task, deriv_run(run)))
}
path_gm_probseg <- function(subject) {
  file.path(FMRIPREP, subject, "anat",
    sprintf("%s_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz", subject))
}
path_anat_mask <- function(subject) {
  file.path(FMRIPREP, subject, "anat",
    sprintf("%s_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz", subject))
}
path_betas <- function(subject) {
  file.path(DERIVED_DIR, sprintf("%s_betas.rds", subject))
}

## ---- Receipts -------------------------------------------------------------
# Slice 1's shape, verbatim: one row is one number, `subject` is the leading
# key and "group" marks rows belonging to the population fit rather than to a
# participant. `tolerance` is finite only on rows asserting an identity, and on
# those rows `passes` is exactly `abs(value) <= tolerance`.
new_receipts <- function() {
  e <- new.env(parent = emptyenv())
  e$rows <- list()
  e
}
record <- function(store, subject, quantity, value,
                   tolerance = NA_real_, note = "") {
  if (length(value) != 1L || length(tolerance) != 1L) {
    stop("receipt \"", quantity, "\" for ", subject, ": value has length ",
         length(value), " and tolerance length ", length(tolerance),
         "; both must be scalars.")
  }
  passes <- if (is.na(tolerance)) NA else abs(value) <= tolerance
  store$rows[[length(store$rows) + 1L]] <- data.frame(
    subject = subject, quantity = quantity, value = as.numeric(value),
    tolerance = as.numeric(tolerance), passes = passes, note = note,
    stringsAsFactors = FALSE
  )
  invisible(value)
}
receipts_frame <- function(store) do.call(rbind, store$rows)

## ---- Small utilities ------------------------------------------------------
say <- function(...) message(sprintf(...))

elapsed <- function(t0) {
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

# The shared functional lattice, read from a file rather than hard-coded, and
# asserted against E12's measurement. Every one of the 84 run masks agrees
# (DECISION.md check 13), so any file will do; using sub-104 trust run 1 means
# the assertion below fails loudly if the fetch ever repoints.
reference_space <- function() {
  v <- neuroim2::read_vol(path_boldmask(SUBJECTS[1L], "trust", 1L))
  sp <- neuroim2::space(v)
  stopifnot(identical(as.integer(dim(sp)), c(66L, 78L, 61L)))
  stopifnot(max(abs(neuroim2::spacing(sp) -
                    c(2.973, 2.973, 3.22))) < 1e-4)
  sp
}

## ---- Shared P^F arithmetic ------------------------------------------------
# Row-normalized softmax over a (row, col) support. Lives here because 02
# builds P^F once and 04 rebuilds it several hundred times under permuted
# fingerprints: if the two used different arithmetic the null band would be a
# band around a different estimator than the one it is meant to price.
#
# No max-subtraction is needed for stability. The similarity is a correlation,
# so `temperature * sim` lies in [-6, +6] and its exponential in
# [0.0025, 403]. Saying so is cheaper than a group-max pass that would only
# guard against a range this construction cannot produce.
softmax_rows <- function(sim, rows, temperature) {
  w <- exp(temperature * sim)
  denom <- as.vector(rowsum(w, rows, reorder = TRUE))
  w / denom[match(rows, sort(unique(rows)))]
}
