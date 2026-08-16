## models.R -- model RDMs and shared constants for the Haxby 2001 exemplar.
##
## Sourced by every numbered script. Defines the experimental basis, the
## searchlight geometry, the model RDM, and the one shared second-order
## statistic that both arms are scored with.

## ---- Experimental basis -------------------------------------------------
## Ordered once here and used everywhere. Row order of every per-run
## condition-mean matrix is exactly this order.
CONDITIONS <- c("bottle", "cat", "chair", "face",
                "house", "scissors", "scrambledpix", "shoe")

## Haxby's animate categories among the eight stimulus classes.
ANIMATE <- c("cat", "face")

## ---- Acquisition / preparation constants --------------------------------
TR_SECONDS   <- 2.5   # Haxby 2001 acquisition TR
HRF_SHIFT_TR <- 2L    # 2 TR = 5 s haemodynamic lag applied to block onsets
N_RUNS       <- 12L

## ---- Searchlight geometry -----------------------------------------------
## neuroim2 searchlight radii are millimetres. The Haxby subject-1 volume has
## 3.5 x 3.75 x 3.75 mm voxels, so 11.25 mm is the README's "radius 3 voxels"
## measured along the 3.75 mm in-plane axes. Both arms are handed this same
## number and, verified in 04-compare.R, receive bit-identical sphere
## membership from neuroim2 at every centre.
RADIUS_MM <- 11.25

## ---- Agreement tolerance ------------------------------------------------
## Gate for "these two systems computed the same number". Both arms run in
## double precision over the same inputs; anything above this is a semantic
## difference, not floating-point drift.
TOLERANCE <- 1e-8

SEED <- 20260814L

## ---- Model RDMs ---------------------------------------------------------
#' Binary animate/inanimate model RDM over CONDITIONS.
#'
#' Zero for pairs sharing an animacy class, one for pairs that cross it, zero
#' diagonal. This is the canonical Haxby category-structure model: the eight
#' stimulus categories are the rows, and the modelled dissimilarity is the
#' animate/inanimate partition over them.
animacy_rdm <- function(conditions = CONDITIONS, animate = ANIMATE) {
  is_animate <- conditions %in% animate
  m <- outer(is_animate, is_animate, function(a, b) as.numeric(a != b))
  dimnames(m) <- list(conditions, conditions)
  diag(m) <- 0
  m
}

MODEL_RDMS <- list(animacy = animacy_rdm())

## ---- The one shared second-order statistic ------------------------------
#' Spearman RSA score of an empirical RDM against a model RDM.
#'
#' Both arms are scored with this identical function applied to their own
#' 28-element lower-triangle RDM vector, so the second-order step contributes
#' no difference between arms. `rdm_vec` and the model must use the same pair
#' ordering, which `rdm_pair_index()` fixes.
#'
#' Returns NA when the empirical RDM is constant (Spearman undefined).
spearman_rsa <- function(rdm_vec, model_vec) {
  if (anyNA(rdm_vec) || stats::sd(rdm_vec) == 0) return(NA_real_)
  suppressWarnings(stats::cor(rdm_vec, model_vec, method = "spearman"))
}

#' Canonical ordering of the 28 unique condition pairs.
#'
#' `utils::combn(q, 2)` column order -- the same order crossform's `rdm()`
#' view reports its pairs in, which 02 asserts against rather than assumes.
rdm_pair_index <- function(conditions = CONDITIONS) {
  p <- utils::combn(seq_along(conditions), 2L)
  data.frame(i = p[1L, ], j = p[2L, ],
             left = conditions[p[1L, ]], right = conditions[p[2L, ]],
             stringsAsFactors = FALSE)
}

#' Lower-triangle vector of a q x q model RDM in `rdm_pair_index()` order.
model_vector <- function(model, conditions = CONDITIONS) {
  idx <- rdm_pair_index(conditions)
  model <- model[conditions, conditions, drop = FALSE]
  model[cbind(idx$i, idx$j)]
}

## ---- Paths --------------------------------------------------------------
exemplar_paths <- function(exemplar_dir) {
  list(
    data    = file.path(exemplar_dir, "data"),
    subj    = file.path(exemplar_dir, "data", "subj1"),
    results = file.path(exemplar_dir, "results"),
    prepared = file.path(exemplar_dir, "data", "prepared-smoke.rds")
  )
}

#' Attach crossform, whether installed or only present as a source tree.
#'
#' The exemplar lives inside the package repository, so a developer checkout
#' normally has no installed copy. Fall back to devtools::load_all on the
#' repository root two levels up.
load_crossform <- function(exemplar_dir) {
  if (requireNamespace("crossform", quietly = TRUE)) {
    suppressMessages(library(crossform))
    return("installed")
  }
  repo <- normalizePath(file.path(exemplar_dir, "..", ".."), mustWork = FALSE)
  if (!file.exists(file.path(repo, "DESCRIPTION"))) {
    stop("crossform is neither installed nor found as a source tree at ", repo)
  }
  suppressMessages(devtools::load_all(repo, quiet = TRUE))
  "source"
}

#' Resolve the exemplar directory whether run by Rscript or sourced.
exemplar_dir_from_args <- function() {
  if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) return(Sys.getenv("EXEMPLAR_DIR"))
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) return(normalizePath(dirname(sub("^--file=", "", f[1]))))
  normalizePath(".")
}
