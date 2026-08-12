# Explicit observation-to-effect extractors --------------------------------

#' Construct an explicit linear effect extractor
#'
#' An extractor is the declared map `E` in `B = E Y`. It contains no neural
#' data and may therefore be reused across feature blocks.
#'
#' @param map A finite effect-by-observation numeric matrix.
#' @param effects An `effect_space()` or unique names used as shorthand for an
#'   unspecified-basis effect space.
#' @param estimator Short estimator identity.
#' @param diagnostics Optional estimator diagnostics.
#' @return An `effect_extractor`.
#' @export
effect_extractor <- function(map, effects = rownames(map),
                             estimator = "explicit", diagnostics = list()) {
  if (!is.matrix(map) || !is.numeric(map) || any(dim(map) < 1L) ||
      any(!is.finite(map))) {
    stop("`map` must be a finite nonempty effect-by-observation matrix.",
      call. = FALSE)
  }
  effects <- .as_effect_space(effects, nrow(map))
  if (!is.character(estimator) || length(estimator) != 1L ||
      is.na(estimator) || !nzchar(estimator)) {
    stop("`estimator` must be one nonempty identifier.", call. = FALSE)
  }
  if (!is.list(diagnostics)) {
    stop("`diagnostics` must be a list.", call. = FALSE)
  }
  rownames(map) <- effects$coordinates
  structure(
    list(
      map = map,
      effect_space = effects,
      effects = effects$coordinates,
      n_observations = ncol(map),
      estimator = estimator,
      diagnostics = diagnostics
    ),
    class = "effect_extractor"
  )
}

#' Compile a supplied linear model into an effect extractor
#'
#' The function accepts an already constructed design matrix; it does not own
#' event timing, HRF construction, nuisance selection, or a formula language.
#' Full-rank designs use pivoted QR. Rank-deficient designs use an SVD
#' pseudoinverse only after every requested effect is proven estimable.
#'
#' @param design Finite observation-by-coefficient design matrix.
#' @param effects Finite effect-by-coefficient target matrix.
#' @param whiten Optional finite square observation whitener `L`.
#' @param effect_names Optional names or an `effect_space()` for target effects.
#' @param tolerance Positive rank and estimability tolerance.
#' @return An `effect_extractor`.
#' @export
lm_extractor <- function(design, effects, whiten = NULL,
                         effect_names = rownames(effects),
                         tolerance = sqrt(.Machine$double.eps)) {
  if (!is.matrix(design) || !is.numeric(design) || any(dim(design) < 1L) ||
      any(!is.finite(design))) {
    stop("`design` must be a finite nonempty observation-by-coefficient matrix.",
      call. = FALSE)
  }
  if (!is.matrix(effects) || !is.numeric(effects) || nrow(effects) < 1L ||
      ncol(effects) != ncol(design) || any(!is.finite(effects))) {
    stop("`effects` must be a finite effect-by-coefficient matrix matching the design.",
      call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || is.na(tolerance) ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop("`tolerance` must be one positive finite number.", call. = FALSE)
  }
  n <- nrow(design)
  if (is.null(whiten)) whiten <- diag(n)
  if (!is.matrix(whiten) || !is.numeric(whiten) ||
      !identical(dim(whiten), c(n, n)) || any(!is.finite(whiten))) {
    stop("`whiten` must be NULL or a finite square observation matrix.",
      call. = FALSE)
  }
  effect_names <- .as_effect_space(effect_names, nrow(effects))
  coordinate_names <- effect_names$coordinates
  whitened_design <- whiten %*% design
  decomposition <- qr(whitened_design, tol = tolerance, LAPACK = FALSE)
  rank <- decomposition$rank
  coefficients <- ncol(design)

  if (rank == coefficients) {
    coefficient_map <- qr.solve(whitened_design, whiten, tol = tolerance)
    solver <- "pivoted_qr"
    estimability_error <- rep(0, nrow(effects))
  } else {
    singular <- svd(whitened_design, nu = min(dim(whitened_design)),
      nv = coefficients)
    cutoff <- tolerance * max(singular$d)
    keep <- singular$d > cutoff
    rank <- sum(keep)
    if (rank < 1L) {
      stop("The whitened design has zero estimable rank.", call. = FALSE)
    }
    basis <- singular$v[, keep, drop = FALSE]
    projected <- effects %*% basis %*% t(basis)
    scale <- pmax(1, sqrt(rowSums(effects^2)))
    estimability_error <- sqrt(rowSums((effects - projected)^2)) / scale
    if (any(estimability_error > tolerance * 10)) {
      bad <- coordinate_names[estimability_error > tolerance * 10]
      stop(sprintf("Requested effects are not estimable: %s.",
        paste(bad, collapse = ", ")), call. = FALSE)
    }
    inverse <- singular$v[, keep, drop = FALSE] %*%
      (t(singular$u[, keep, drop = FALSE]) / singular$d[keep])
    coefficient_map <- inverse %*% whiten
    solver <- "svd_estimable_fallback"
  }

  effect_extractor(
    effects %*% coefficient_map,
    effects = effect_names,
    estimator = "linear_model",
    diagnostics = list(
      solver = solver,
      observations = n,
      coefficients = coefficients,
      rank = as.integer(rank),
      rank_deficient = rank < coefficients,
      estimability_error = stats::setNames(estimability_error, coordinate_names),
      tolerance = tolerance
    )
  )
}

.validate_effect_names <- function(effects, expected) {
  if (is.null(effects)) effects <- paste0("effect", seq_len(expected))
  if (!is.character(effects) || length(effects) != expected || anyNA(effects) ||
      any(!nzchar(effects)) || anyDuplicated(effects)) {
    stop("Effect coordinates must have unique nonempty names.", call. = FALSE)
  }
  effects
}

.validate_effect_extractor <- function(x) {
  if (!inherits(x, "effect_extractor") || !is.list(x) ||
      !identical(names(x), c("map", "effect_space", "effects",
        "n_observations", "estimator", "diagnostics"))) {
    stop("Extractor fields are missing or noncanonical.", call. = FALSE)
  }
  rebuilt <- effect_extractor(x$map, x$effect_space, x$estimator, x$diagnostics)
  if (!identical(x$effects, rebuilt$effects)) {
    stop("Extractor coordinate labels are inconsistent with its effect space.",
      call. = FALSE)
  }
  if (!identical(x$n_observations, rebuilt$n_observations)) {
    stop("Extractor observation metadata is inconsistent with its map.",
      call. = FALSE)
  }
  rebuilt
}
