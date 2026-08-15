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
#' @param observation_whitener Optional finite square observation whitener `L`.
#' @param effect_names Optional names or an `effect_space()` for target effects.
#' @param tolerance Positive rank and estimability tolerance.
#' @param whiten Deprecated alias for `observation_whitener`. Supplying both is
#'   an error.
#' @return An `effect_extractor`.
#' @export
lm_extractor <- function(design, effects, observation_whitener = NULL,
                         effect_names = rownames(effects),
                         tolerance = sqrt(.Machine$double.eps),
                         whiten = NULL) {
  resolved <- .resolve_observation_whitener(
    observation_whitener, whiten, nrow(design),
    legacy_supplied = !missing(whiten) && !is.null(whiten)
  )
  .compile_lm_estimator(
    design, effects, resolved, effect_names, tolerance
  )$extractor
}

.resolve_observation_whitener <- function(observation_whitener, whiten,
                                          observations, legacy_supplied) {
  if (isTRUE(legacy_supplied)) {
    if (!is.null(observation_whitener)) {
      stop("Supply only `observation_whitener`; `whiten` is its deprecated alias.",
        call. = FALSE)
    }
    warning("`whiten` is deprecated; use `observation_whitener`.",
      call. = FALSE)
    observation_whitener <- whiten
  }
  if (!is.numeric(observations) || length(observations) != 1L ||
      is.na(observations) || observations < 1L || observations %% 1 != 0) {
    stop("The observation count must be one positive integer.", call. = FALSE)
  }
  observations <- as.integer(observations)
  if (is.null(observation_whitener)) {
    semantic <- list(kind = "identity", dim = c(observations, observations))
    return(structure(list(
      matrix = NULL,
      identity = TRUE,
      descriptor = c(semantic, list(signature = paste0(
        "sha256:", digest::digest(semantic, algo = "sha256", serialize = TRUE)
      )))
    ), class = "effect_observation_whitener"))
  }
  if (!is.matrix(observation_whitener) ||
      !is.numeric(observation_whitener) ||
      !identical(dim(observation_whitener), c(observations, observations)) ||
      any(!is.finite(observation_whitener))) {
    stop(paste0(
      "`observation_whitener` must be NULL or a finite square observation ",
      "matrix."
    ), call. = FALSE)
  }
  semantic <- list(
    kind = "explicit",
    dim = as.integer(dim(observation_whitener)),
    matrix_revision = paste0("sha256:", digest::digest(
      observation_whitener, algo = "sha256", serialize = TRUE
    ))
  )
  structure(list(
    matrix = observation_whitener,
    identity = FALSE,
    descriptor = c(semantic, list(signature = paste0(
      "sha256:", digest::digest(semantic, algo = "sha256", serialize = TRUE)
    )))
  ), class = "effect_observation_whitener")
}

.compile_lm_estimator <- function(design, effects, observation_whitener,
                                  effect_names, tolerance, partition = NULL) {
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
  if (!inherits(observation_whitener, "effect_observation_whitener") ||
      !is.list(observation_whitener) ||
      !identical(names(observation_whitener),
        c("matrix", "identity", "descriptor"))) {
    stop("Observation-whitener compilation metadata are invalid.",
      call. = FALSE)
  }
  effect_names <- .as_effect_space(effect_names, nrow(effects))
  coordinate_names <- effect_names$coordinates
  whitener <- observation_whitener$matrix
  whitened_design <- if (isTRUE(observation_whitener$identity)) {
    design
  } else {
    whitener %*% design
  }
  decomposition <- qr(whitened_design, tol = tolerance, LAPACK = FALSE)
  rank <- decomposition$rank
  coefficients <- ncol(design)

  if (rank == coefficients) {
    residual_basis <- qr.Q(decomposition, complete = FALSE)[,
      seq_len(rank), drop = FALSE]
    triangular <- qr.R(decomposition)[seq_len(rank), seq_len(rank), drop = FALSE]
    pivoted_map <- backsolve(triangular, t(residual_basis))
    coefficient_whitened_map <- matrix(0, coefficients, n)
    coefficient_whitened_map[decomposition$pivot, ] <- pivoted_map
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
      null_basis <- singular$v[, !keep, drop = FALSE]
      coefficient_names <- colnames(design)
      if (is.null(coefficient_names)) {
        coefficient_names <- paste0("coefficient", seq_len(coefficients))
      }
      alias_groups <- unique(vapply(seq_len(ncol(null_basis)), function(index) {
        values <- null_basis[, index]
        active <- abs(values) > tolerance * 10 * max(1, max(abs(values)))
        paste(coefficient_names[active], collapse = ", ")
      }, character(1)))
      alias_groups <- alias_groups[nzchar(alias_groups)]
      location <- if (is.null(partition)) {
        "The design"
      } else {
        sprintf("Partition `%s`", partition)
      }
      reasons <- c(
        sprintf("%s has rank %d for %d regressors.",
          location, rank, coefficients),
        if (length(alias_groups)) {
          sprintf("Aliased regressor set%s: %s.",
            if (length(alias_groups) == 1L) "" else "s",
            paste(alias_groups, collapse = "; "))
        },
        sprintf("Requested effect%s outside the estimable row space: %s.",
          if (length(bad) == 1L) "" else "s",
          paste(bad, collapse = ", "))
      )
      remedies <- c(
        paste0(
          "Remove a redundant design column. For a complete set of condition ",
          "indicators, drop the intercept or one condition indicator."
        ),
        paste0(
          "Request estimable contrasts, such as differences between ",
          "conditions, instead of unidentified individual coefficients."
        )
      )
      message <- paste0(
        "Requested effects are not estimable: ", paste(bad, collapse = ", "),
        ". ", paste(reasons, collapse = " "), " ",
        paste(remedies, collapse = " ")
      )
      .capability_refusal(
        message,
        capability = "estimable_effects",
        namespace = "relation_fit",
        reasons = reasons,
        remedies = remedies
      )
    }
    inverse <- singular$v[, keep, drop = FALSE] %*%
      (t(singular$u[, keep, drop = FALSE]) / singular$d[keep])
    coefficient_whitened_map <- inverse
    residual_basis <- singular$u[, keep, drop = FALSE]
    solver <- "svd_estimable_fallback"
  }

  coefficient_map <- if (isTRUE(observation_whitener$identity)) {
    coefficient_whitened_map
  } else {
    coefficient_whitened_map %*% whitener
  }
  effect_whitened_map <- effects %*% coefficient_whitened_map
  effect_covariance <- tcrossprod(effect_whitened_map)
  dimnames(effect_covariance) <- list(coordinate_names, coordinate_names)
  residual_df <- as.integer(n - rank)
  extractor <- effect_extractor(
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
      tolerance = tolerance,
      observation_whitener = observation_whitener$descriptor
    )
  )
  residualize <- function(response) {
    if (!is.matrix(response) || !is.numeric(response) ||
        nrow(response) != n || any(!is.finite(response))) {
      stop("Residualization requires a finite observation-by-feature block.",
        call. = FALSE)
    }
    transformed <- if (isTRUE(observation_whitener$identity)) {
      response
    } else {
      whitener %*% response
    }
    transformed - residual_basis %*% crossprod(residual_basis, transformed)
  }
  semantic <- list(
    schema_version = 1L,
    design_revision = paste0("sha256:", digest::digest(
      design, algo = "sha256", serialize = TRUE
    )),
    target_revision = paste0("sha256:", digest::digest(
      effects, algo = "sha256", serialize = TRUE
    )),
    observation_whitener = observation_whitener$descriptor,
    effect_space = effect_names,
    solver = solver,
    rank = as.integer(rank),
    tolerance = tolerance
  )
  list(
    extractor = extractor,
    effect_covariance = effect_covariance,
    effect_covariance_convention =
      "neural_covariance_factor_excluded",
    residual_df = residual_df,
    residualize = residualize,
    observation_whitener = observation_whitener$descriptor,
    estimator_provenance = c(semantic, list(signature = paste0(
      "sha256:", digest::digest(semantic, algo = "sha256", serialize = TRUE)
    )))
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
