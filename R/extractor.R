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
#' @return An `effect_extractor`: a list with the effect-by-observation `$map`,
#'   its `$effect_space` and `$effects` labels, `$n_observations`, the
#'   `$estimator` identity, and `$diagnostics`.
#' @family relation planning and fitting
#' @seealso [lm_extractor()] to compile one from a design and target, and
#'   [relation()], which applies one extractor per partition.
#' @examples
#' # Average four observations, and contrast the first two with the last two.
#' map <- rbind(
#'   mean = c(0.25, 0.25, 0.25, 0.25),
#'   difference = c(0.5, 0.5, -0.5, -0.5)
#' )
#' extractor <- effect_extractor(map, estimator = "hand-specified")
#' extractor$effects
#' extractor$n_observations
#'
#' # The extractor holds no neural data, so the same map is reused across every
#' # feature block a relation reads.
#' set.seed(1)
#' extractor$map %*% matrix(rnorm(8), 4L, 2L)
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
#' @param solver Numerical factorization route: automatic, pivoted QR, or SVD.
#'   This changes execution provenance, not the requested effect.
#' @return An `effect_extractor` whose `$map` is `T (L X)^+ L`, with
#'   `$diagnostics` recording `solver`, `solver_policy`, `observations`,
#'   `coefficients`, `rank`, `rank_deficient`, per-effect
#'   `estimability_error`, `tolerance`, and the `observation_whitener`
#'   descriptor.
#' @family relation planning and fitting
#' @seealso [effect_extractor()] for a hand-specified map, and
#'   [lm_relation_fit()], which additionally keeps the residual error channel.
#' @examples
#' condition <- factor(rep(c("face", "body"), each = 3L))
#' design <- cbind(
#'   stats::model.matrix(~ 0 + condition), drift = seq(-1, 1, length.out = 6L)
#' )
#' colnames(design)[1:2] <- c("face", "body")
#'
#' # Ask for one contrast on the design's coefficient axis.
#' target <- rbind(`face-body` = c(1, -1, 0))
#' colnames(target) <- colnames(design)
#' extractor <- lm_extractor(design, target)
#' extractor$diagnostics$solver
#' round(extractor$map, 3)
#'
#' # Adding an intercept aliases the two condition columns, so a request for
#' # `face` alone is refused and the aliased regressors are named.
#' aliased <- cbind(intercept = 1, design)
#' request <- matrix(0, 1L, ncol(aliased),
#'   dimnames = list("face", colnames(aliased)))
#' request[, "face"] <- 1
#' refusal <- catch_refusal(lm_extractor(aliased, request))
#' refusal$capability
#' refusal$reasons
#' @export
lm_extractor <- function(design, effects, observation_whitener = NULL,
                         effect_names = rownames(effects),
                         tolerance = sqrt(.Machine$double.eps),
                         whiten = NULL,
                         solver = c("auto", "qr", "svd")) {
  solver <- match.arg(solver)
  resolved <- .resolve_observation_whitener(
    observation_whitener, whiten, nrow(design),
    legacy_supplied = !missing(whiten) && !is.null(whiten)
  )
  .compile_lm_estimator(
    design, effects, resolved, effect_names, tolerance, solver = solver
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
        "sha256:", digest::digest(semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L)
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
      observation_whitener, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    ))
  )
  structure(list(
    matrix = observation_whitener,
    identity = FALSE,
    descriptor = c(semantic, list(signature = paste0(
      "sha256:", digest::digest(semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L)
    )))
  ), class = "effect_observation_whitener")
}

.compile_lm_estimator <- function(design, effects, observation_whitener,
                                  effect_names, tolerance, partition = NULL,
                                  solver = c("auto", "qr", "svd")) {
  solver <- match.arg(solver)
  where <- if (is.null(partition)) "" else sprintf("Partition `%s`: ", partition)
  if (!is.matrix(design) || !is.numeric(design) || any(dim(design) < 1L) ||
      any(!is.finite(design))) {
    stop(sprintf(paste0(
      "%s`design` must be a finite nonempty observation-by-coefficient ",
      "matrix, one row per observation; received %s."
    ), where, .msg_value(design)), call. = FALSE)
  }
  if (!is.matrix(effects) || !is.numeric(effects) || nrow(effects) < 1L ||
      any(!is.finite(effects))) {
    stop(sprintf(paste0(
      "%s`effects` must be a finite effect-by-coefficient matrix saying which ",
      "linear combination of design coefficients each effect is; received %s."
    ), where, .msg_value(effects)), call. = FALSE)
  }
  if (ncol(effects) != ncol(design)) {
    stop(sprintf(paste0(
      "%s`effects` has %s but the design has %s, and every effect must be a ",
      "combination of the design's coefficients. Supply one `effects` column ",
      "per design column, in design column order."
    ), where, .msg_count(ncol(effects), "column"),
      .msg_count(ncol(design), "coefficient column")), call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || is.na(tolerance) ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop(sprintf("`tolerance` must be one positive finite number; received %s.",
      .msg_value(tolerance)), call. = FALSE)
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

  if (identical(solver, "qr") && rank < coefficients) {
    location <- if (is.null(partition)) "The design" else {
      sprintf("Partition `%s`", partition)
    }
    .capability_refusal(
      sprintf(
        "%s has rank %d for %d regressors, so the requested QR route is unavailable.",
        location, rank, coefficients
      ),
      capability = "full_rank_design",
      namespace = "relation_fit",
      reasons = sprintf("%s is rank deficient.", location),
      remedies = c(
        "Use `solver = \"auto\"` or `solver = \"svd\"` for estimable effects.",
        "Remove redundant design columns."
      )
    )
  }

  use_qr <- rank == coefficients && solver %in% c("auto", "qr")
  if (use_qr) {
    residual_basis <- qr.Q(decomposition, complete = FALSE)[,
      seq_len(rank), drop = FALSE]
    triangular <- qr.R(decomposition)[seq_len(rank), seq_len(rank), drop = FALSE]
    pivoted_map <- backsolve(triangular, t(residual_basis))
    coefficient_whitened_map <- matrix(0, coefficients, n)
    coefficient_whitened_map[decomposition$pivot, ] <- pivoted_map
    solver_used <- "pivoted_qr"
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
    solver_used <- if (identical(solver, "svd")) {
      "svd_requested"
    } else {
      "svd_estimable_fallback"
    }
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
      solver = solver_used,
      solver_policy = solver,
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
      design, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    )),
    target_revision = paste0("sha256:", digest::digest(
      effects, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    )),
    observation_whitener = observation_whitener$descriptor,
    effect_space = effect_names,
    solver = solver_used,
    solver_policy = solver,
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
      "sha256:", digest::digest(semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L)
    )))
  )
}

.validate_effect_names <- function(effects, expected) {
  if (is.null(effects)) effects <- paste0("effect", seq_len(expected))
  if (!is.character(effects)) {
    stop(sprintf(
      "Effect coordinates must be a character vector; received %s.",
      .msg_value(effects)), call. = FALSE)
  }
  if (length(effects) != expected) {
    stop(sprintf(
      "Effect coordinates must name %s; received %s (%s).",
      .msg_count(expected, "coordinate"),
      .msg_count(length(effects), "name"), .msg_names(effects)),
      call. = FALSE)
  }
  if (anyNA(effects) || any(!nzchar(effects)) || anyDuplicated(effects)) {
    stop(sprintf(
      "Effect coordinates must have unique nonempty names%s.",
      if (anyDuplicated(effects)) {
        sprintf("; %s appears more than once",
          .msg_names(unique(effects[duplicated(effects)])))
      } else {
        "; some are missing or empty"
      }), call. = FALSE)
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
