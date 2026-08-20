#!/usr/bin/env Rscript
# 00-common.R -- shared fixture definition, paths, and crossform loading.
#
# Sourced by 01-fixture.R, 03-compare.R, and 04-extension.R. It holds nothing
# that is specific to one arm; the fixture constants live here so the three
# R scripts cannot drift apart.

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

## ---- Fixture constants --------------------------------------------------
# Condition names are chosen so that alphabetical order (what rsatoolbox's
# `sort_by()` and `np.unique()` impose) is the same as the declared effect
# order in crossform. That removes one whole class of silent mismatch.
SEED <- 20260817L
CONDITIONS <- c("cond1_faceA", "cond2_faceB", "cond3_houseA",
                "cond4_houseB", "cond5_toolA", "cond6_toolB")
N_RUNS <- 4L
N_VOXELS <- 40L
TRIALS_PER_CONDITION <- 8L
REGION_LABELS <- rep(c("roiA", "roiB", "roiC"), times = c(16L, 14L, 10L))

# The contrast used by the strict-extension arm: mean(face) - mean(house).
CONTRAST <- setNames(c(0.5, 0.5, -0.5, -0.5, 0, 0), CONDITIONS)

#' Vectorise a square RDM in crossform's declared pair order.
#'
#' crossform's `.rdm_vector()` and `.pair_difference_query()` both use
#' `utils::combn(seq_len(q), 2)`, i.e. the row-major upper triangle
#' (1,2), (1,3), ..., (1,q), (2,3), ... This is the same order numpy's
#' `np.triu_indices(q, 1)` produces, which is what makes the two vectorised
#' RDMs comparable element by element without a permutation.
rdm_pair_vector <- function(value) {
  pairs <- utils::combn(seq_len(nrow(value)), 2L)
  value[cbind(pairs[1L, ], pairs[2L, ])]
}

#' Two model RDMs over the six conditions.
model_rdms <- function() {
  category <- c(face = 1, face = 1, house = 2, house = 2, tool = 3, tool = 3)
  animacy <- c(1, 1, 0, 0, 0, 0)
  build <- function(code) {
    m <- outer(code, code, function(a, b) as.numeric(a != b))
    dimnames(m) <- list(CONDITIONS, CONDITIONS)
    m
  }
  list(category = build(category), animacy = build(animacy))
}

#' The deterministic fixture: raw per-run responses, design, truth.
#'
#' Residual noise is drawn with a known non-spherical covariance so the noise
#' metric is not a decoration: an identity metric would give different
#' numbers, and both packages must agree under the same non-identity one.
build_fixture <- function() {
  set.seed(SEED)
  q <- length(CONDITIONS)
  n_obs <- q * TRIALS_PER_CONDITION

  design <- model.matrix(~ 0 + factor(rep(seq_len(q), TRIALS_PER_CONDITION)))
  colnames(design) <- CONDITIONS
  attr(design, "assign") <- NULL
  attr(design, "contrasts") <- NULL
  design <- matrix(as.numeric(design), n_obs, q,
                   dimnames = list(NULL, CONDITIONS))

  effects <- diag(q)
  dimnames(effects) <- list(CONDITIONS, CONDITIONS)

  # Non-spherical residual covariance: AR-like correlation across voxels
  # times heterogeneous per-voxel scale, plus a low-rank bump so it is not
  # exactly Toeplitz.
  correlation <- toeplitz(0.65^(seq_len(N_VOXELS) - 1L))
  bump <- outer(cos(seq_len(N_VOXELS) / 3), cos(seq_len(N_VOXELS) / 3))
  scale <- seq(0.6, 1.6, length.out = N_VOXELS)
  covariance <- diag(scale) %*% (correlation + 0.35 * bump) %*% diag(scale)
  covariance <- (covariance + t(covariance)) / 2
  covariance <- covariance + diag(0.05, N_VOXELS)
  factor <- chol(covariance)

  # Planted condition patterns: a category axis, an animacy axis, and
  # exemplar-specific detail, each on a different voxel block.
  truth <- matrix(0, q, N_VOXELS, dimnames = list(CONDITIONS, NULL))
  truth[, 1:12] <- rep(c(0.8, 0.8, -0.4, -0.4, -0.4, -0.4), 12) *
    rep(seq(1, 0.4, length.out = 12), each = q)
  truth[, 13:26] <- outer(c(0.5, 0.5, 0.5, 0.5, -0.9, -0.9),
                          seq(0.9, 0.3, length.out = 14))
  truth[, 27:34] <- outer(c(0.6, -0.6, 0.5, -0.5, 0.4, -0.4),
                          seq(0.8, 0.2, length.out = 8))

  responses <- setNames(lapply(seq_len(N_RUNS), function(run) {
    noise <- matrix(rnorm(n_obs * N_VOXELS), n_obs, N_VOXELS) %*% factor
    design %*% truth + noise
  }), paste0("run", seq_len(N_RUNS)))

  list(
    responses = responses, design = design, effects = effects,
    truth = truth, covariance = covariance,
    conditions = CONDITIONS, n_obs = n_obs, q = q
  )
}

#' Pooled residual covariance / precision, computed once in R.
#'
#' `noise_precision()` in crossform is a *constructor for a fixed metric*, not
#' an estimator: it takes the precision the analyst supplies and records that
#' it was fixed before effect evaluation. So the estimator below is the
#' exemplar's, not crossform's, and the identical matrix is handed to
#' rsatoolbox. `03-compare.R` additionally checks that rsatoolbox's own
#' `prec_from_residuals(..., method = "full", dof = nu)` reproduces it.
pooled_precision <- function(fixture) {
  x <- fixture$design
  hat <- diag(nrow(x)) - x %*% solve(crossprod(x), t(x))
  residuals <- lapply(fixture$responses, function(y) hat %*% y)
  crossproducts <- Reduce(`+`, lapply(residuals, crossprod))
  df <- length(residuals) * (nrow(x) - qr(x)$rank)
  covariance <- crossproducts / df
  list(
    residuals = residuals,
    residual_df = df,
    covariance = covariance,
    precision = solve(covariance)
  )
}

TOLERANCE <- 1e-10
