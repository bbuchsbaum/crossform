# Statistical envelopes for effect relations -------------------------------

.residual_source <- function(dim, read, read_padded, stable_revision) {
  if (!is.numeric(dim) || length(dim) != 2L || anyNA(dim) ||
      any(!is.finite(dim)) || any(dim < 1L) || any(dim %% 1 != 0) ||
      !is.function(read) || !is.function(read_padded) ||
      !.strong_sha256(stable_revision)) {
    stop("Residual-source metadata are invalid.", call. = FALSE)
  }
  structure(list(
    dim = as.integer(dim),
    read = read,
    read_padded = read_padded,
    stable_revision = stable_revision
  ), class = "effect_residual_source")
}

.validate_whitener_descriptor <- function(x, deep = TRUE) {
  if (!is.list(x) || !identical(names(x), if (identical(x$kind, "identity")) {
    c("kind", "dim", "signature")
  } else {
    c("kind", "dim", "matrix_revision", "signature")
  }) || !x$kind %in% c("identity", "explicit") ||
      !is.integer(x$dim) || length(x$dim) != 2L || x$dim[[1L]] < 1L ||
      x$dim[[1L]] != x$dim[[2L]] || !.strong_sha256(x$signature) ||
      (identical(x$kind, "explicit") &&
       !.strong_sha256(x$matrix_revision))) {
    stop("Observation-whitener identity is invalid.", call. = FALSE)
  }
  if (isTRUE(deep)) {
    semantic <- x[names(x) != "signature"]
    expected <- .sha256_signature(semantic)
    if (!identical(x$signature, expected)) {
      stop("Observation-whitener identity is inconsistent.", call. = FALSE)
    }
  }
  x
}

.error_model_signature <- function(model) {
  semantic <- list(
    schema_version = 1L,
    kind = model$kind,
    assumptions = model$assumptions,
    scale_convention = model$scale_convention,
    effect_covariance = model$effect_covariance,
    residual_df = model$residual_df,
    sampling_unit = model$sampling_unit,
    observation_whitener = model$observation_whitener,
    source_revision = model$source_revision,
    residual_revision = model$residual_source$stable_revision,
    residual_dim = model$residual_source$dim,
    estimator_provenance = model$estimator_provenance,
    provenance = model$provenance
  )
  .sha256_signature(semantic)
}

.separable_glm_error_model <- function(effect_covariance, residual_df,
                                       sampling_unit,
                                       observation_whitener,
                                       source_revision, residual_source,
                                       estimator_provenance,
                                       provenance = list()) {
  model <- structure(list(
    kind = "separable_glm",
    assumptions = c(
      "matrix_normal_separable",
      "shared_observation_covariance_across_neural_features",
      "observation_whitener_treated_as_fixed"
    ),
    scale_convention = "neural_covariance_factor_excluded",
    effect_covariance = effect_covariance,
    residual_df = as.integer(residual_df),
    sampling_unit = sampling_unit,
    observation_whitener = observation_whitener,
    source_revision = source_revision,
    residual_source = residual_source,
    estimator_provenance = estimator_provenance,
    provenance = provenance,
    signature = NA_character_
  ), class = "effect_error_model")
  model$signature <- .error_model_signature(model)
  .validate_error_model(model, deep = TRUE)
}

.validate_error_model <- function(model, deep = TRUE) {
  expected <- c("kind", "assumptions", "scale_convention",
    "effect_covariance", "residual_df", "sampling_unit",
    "observation_whitener", "source_revision", "residual_source",
    "estimator_provenance", "provenance", "signature")
  if (!inherits(model, "effect_error_model") || !is.list(model) ||
      !identical(names(model), expected) ||
      !identical(model$kind, "separable_glm") ||
      !identical(model$assumptions, c(
        "matrix_normal_separable",
        "shared_observation_covariance_across_neural_features",
        "observation_whitener_treated_as_fixed"
      )) ||
      !identical(model$scale_convention,
        "neural_covariance_factor_excluded") ||
      !is.matrix(model$effect_covariance) ||
      !is.numeric(model$effect_covariance) ||
      nrow(model$effect_covariance) < 1L ||
      nrow(model$effect_covariance) != ncol(model$effect_covariance) ||
      any(!is.finite(model$effect_covariance)) ||
      max(abs(model$effect_covariance - t(model$effect_covariance))) > 1e-10 ||
      !is.integer(model$residual_df) || length(model$residual_df) != 1L ||
      is.na(model$residual_df) || model$residual_df < 1L ||
      !is.character(model$sampling_unit) ||
      length(model$sampling_unit) != 1L || is.na(model$sampling_unit) ||
      !nzchar(model$sampling_unit) || !.strong_sha256(model$source_revision) ||
      !is.list(model$estimator_provenance) || !is.list(model$provenance) ||
      !.strong_sha256(model$signature)) {
    stop("Error-model fields are missing or noncanonical.", call. = FALSE)
  }
  if (isTRUE(deep)) {
    eigenvalues <- eigen(model$effect_covariance,
      symmetric = TRUE, only.values = TRUE)$values
    covariance_scale <- max(1, max(abs(diag(model$effect_covariance))))
    if (min(eigenvalues) < -1e-10 * covariance_scale) {
      stop("Effect-coordinate covariance must be positive semidefinite.",
        call. = FALSE)
    }
  }
  whitener <- .validate_whitener_descriptor(
    model$observation_whitener, deep = deep
  )
  residual <- model$residual_source
  if (!inherits(residual, "effect_residual_source") || !is.list(residual) ||
      !identical(names(residual),
        c("dim", "read", "read_padded", "stable_revision")) ||
      !is.integer(residual$dim) || length(residual$dim) != 2L ||
      residual$dim[[1L]] != whitener$dim[[1L]] ||
      residual$dim[[2L]] < 1L || !is.function(residual$read) ||
      !is.function(residual$read_padded) ||
      !.strong_sha256(residual$stable_revision)) {
    stop("Error-model residual source is invalid.", call. = FALSE)
  }
  if (isTRUE(deep) && !identical(model$signature,
      .error_model_signature(model))) {
    stop("Error-model identity is inconsistent.", call. = FALSE)
  }
  model
}

.error_capabilities <- function(model) {
  # Each capability is earned by the field that actually supplies it, so a
  # future error-model variant cannot inherit a capability it lacks. Under
  # the current separable-GLM model every field is mandatory and the bits
  # coincide; the derivations keep that an invariant, not an assumption.
  present <- !is.null(model)
  separable <- present && identical(model$kind, "separable_glm")
  residuals <- present &&
    inherits(model$residual_source, "effect_residual_source")
  covariance <- present && is.matrix(model$effect_covariance)
  degrees <- present && is.integer(model$residual_df) &&
    length(model$residual_df) == 1L && !is.na(model$residual_df) &&
    model$residual_df >= 1L
  structure(list(
    error_model = present,
    residual_blocks = residuals,
    effect_covariance = covariance,
    residual_df = degrees,
    separable_error = separable,
    learned_metric_input = residuals,
    within_participant_calibration = separable && residuals &&
      covariance && degrees
  ), class = "effect_error_capabilities")
}

.relation_fit_signature <- function(relation, error_models, provenance) {
  semantic <- list(
    schema_version = 1L,
    relation = .relation_family_identity(relation),
    errors = vapply(error_models, function(model) {
      if (is.null(model)) NA_character_ else model$signature
    }, character(1)),
    provenance = provenance
  )
  .sha256_signature(semantic)
}

#' Attach statistical error information to a relation
#'
#' `relation_fit()` keeps the mathematical `effect_relation` intact and adds a
#' separate, capability-bearing error channel. A relation can remain valid
#' without this envelope; operations that learn a neural metric or calibrate
#' within-participant uncertainty must require the corresponding capability.
#'
#' @param relation A validated `effect_relation`.
#' @param error_models Optional named list of internal error-model values, one
#'   per relation partition. Omitting it records an explicitly absent channel.
#' @param provenance Compact fit-level provenance.
#' @return An `effect_relation_fit`: a list with the untouched `$relation`, one
#'   `$error_models` entry per partition (`NULL` where absent), the derived
#'   per-partition `$capabilities`, `$provenance`, and a `$signature` binding
#'   the relation and every error model.
#' @family relation planning and fitting
#' @seealso [lm_relation_fit()], the route that installs a real residual
#'   channel, and [relation_fit_capabilities()] to inspect the result.
#' @examples
#' set.seed(1)
#' domain <- abstract_domain(4L, id = "relation-fit-example")
#' betas <- matrix(rnorm(8), 2L, 4L,
#'   dimnames = list(c("face", "body"), NULL))
#' point <- relation(list(`run-1` = betas), domain = domain)
#'
#' # Wrapping a point relation records an explicitly absent error channel,
#' # which is different from an unstated one.
#' fit <- relation_fit(point)
#' relation_fit_capabilities(fit)
#'
#' # The relation itself is unchanged, so point geometry still works.
#' identical(fit$relation, point)
#' @export
relation_fit <- function(relation, error_models = NULL, provenance = list()) {
  .validate_relation(relation)
  .relation_source_capabilities(relation)
  if (is.null(error_models)) {
    error_models <- stats::setNames(
      rep(list(NULL), length(relation$partitions)), relation$partitions
    )
  }
  if (!is.list(error_models) ||
      !identical(names(error_models), relation$partitions)) {
    stop("`error_models` must provide one named value per relation partition.",
      call. = FALSE)
  }
  error_models <- lapply(error_models, function(model) {
    if (is.null(model)) NULL else .validate_error_model(model, deep = TRUE)
  })
  names(error_models) <- relation$partitions
  .validate_effect_provenance(provenance, "relation-fit provenance")
  capabilities <- lapply(error_models, .error_capabilities)
  names(capabilities) <- relation$partitions
  value <- structure(list(
    relation = relation,
    error_models = error_models,
    capabilities = capabilities,
    provenance = provenance,
    signature = .relation_fit_signature(relation, error_models, provenance)
  ), class = "effect_relation_fit")
  .validate_relation_fit(value, deep = TRUE)
}

.validate_relation_fit <- function(x, deep = TRUE) {
  if (.validated_before(x, "relation_fit", deep)) return(invisible(x))
  expected <- c("relation", "error_models", "capabilities", "provenance",
    "signature")
  if (!inherits(x, "effect_relation_fit") || !is.list(x) ||
      !identical(names(x), expected)) {
    stop("Relation-fit fields are missing or noncanonical.", call. = FALSE)
  }
  .validate_relation(x$relation, deep = deep)
  partitions <- x$relation$partitions
  if (!is.list(x$error_models) ||
      !identical(names(x$error_models), partitions) ||
      !is.list(x$capabilities) ||
      !identical(names(x$capabilities), partitions) ||
      !is.list(x$provenance) || !.strong_sha256(x$signature)) {
    stop("Relation-fit metadata are inconsistent.", call. = FALSE)
  }
  source_capabilities <- x$relation$capabilities
  q <- length(x$relation$effects)
  for (partition in partitions) {
    model <- x$error_models[[partition]]
    expected_capability <- .error_capabilities(model)
    if (!identical(x$capabilities[[partition]], expected_capability)) {
      stop("Relation-fit error capabilities are inconsistent.",
        call. = FALSE)
    }
    if (is.null(model)) next
    model <- .validate_error_model(model, deep = deep)
    if (!identical(dim(model$effect_covariance), c(q, q)) ||
        !identical(rownames(model$effect_covariance), x$relation$effects) ||
        !identical(colnames(model$effect_covariance), x$relation$effects) ||
        model$residual_source$dim[[2L]] != x$relation$n_features ||
        is.null(source_capabilities) ||
        !identical(tolower(model$source_revision),
          tolower(source_capabilities[[partition]]$stable_revision))) {
      stop("Relation-fit error model is incompatible with its relation.",
        call. = FALSE)
    }
  }
  if (isTRUE(deep)) {
    .validate_effect_provenance(x$provenance, "relation-fit provenance")
    expected_signature <- .relation_fit_signature(
      x$relation, x$error_models, x$provenance
    )
    if (!identical(x$signature, expected_signature)) {
      stop("Relation-fit identity is inconsistent.", call. = FALSE)
    }
  }
  .record_validated(x, "relation_fit", deep)
  invisible(x)
}

.partition_values <- function(value, partitions, what, allow_null = FALSE) {
  if (is.null(value)) {
    if (!allow_null) {
      stop(sprintf(paste0(
        "`%s` is required: supply one value shared by every partition, or a ",
        "named list with one entry per partition (%s)."
      ), what, .msg_names(partitions)), call. = FALSE)
    }
    return(stats::setNames(rep(list(NULL), length(partitions)), partitions))
  }
  if (!is.list(value) || is.data.frame(value) ||
      inherits(value, "effect_space")) {
    return(stats::setNames(rep(list(value), length(partitions)), partitions))
  }
  if (length(value) != length(partitions)) {
    stop(sprintf(paste0(
      "`%s` is a list of %d, but there %s %s (%s). Supply one entry per ",
      "partition, or a single value shared by all of them."
    ), what, length(value), if (length(partitions) == 1L) "is" else "are",
      .msg_count(length(partitions), "partition"), .msg_names(partitions)),
      call. = FALSE)
  }
  if (!is.null(names(value))) {
    if (!setequal(names(value), partitions) || anyDuplicated(names(value))) {
      stop(sprintf(paste0(
        "Named `%s` entries must match the partitions exactly: %s were ",
        "supplied for partitions %s."
      ), what, .msg_names(names(value)), .msg_names(partitions)),
        call. = FALSE)
    }
    value <- value[partitions]
  }
  names(value) <- partitions
  value
}

#' Fit a partitioned linear-model relation with a residual error channel
#'
#' The adapter compiles each supplied design once. Its pure effect map forms the
#' underlying relation; its orthogonal residual projection, unscaled
#' effect-coordinate covariance, residual degrees of freedom, and source
#' revision form a separate separable-GLM error capability. Residual blocks are
#' produced lazily without constructing a dense observation residualizer.
#' This error channel is required by the admitted fixed-metric analytic RDM
#' covariance in [rdm_sampling_covariance()]. The separable GLM is a declared
#' capability: it assumes the supplied observation model yields a common
#' effect-coordinate covariance structure across neural features and, for the
#' equal-partition specialization, across partitions.
#'
#' @section Independent observations, or a whitener:
#' The error channel this function records --- the effect covariance
#' \eqn{\Xi=A A^\top} and the residual degrees of freedom \eqn{\nu=n-\mathrm{rank}(X)}
#' --- describes observations that are independent given the design. fMRI
#' residuals are not: they are temporally autocorrelated. Fitting without
#' `observation_whitener` therefore leaves \eqn{\Xi} as the plain OLS factor,
#' which is not the covariance of the estimates under correlated errors, and
#' leaves \eqn{\nu} counting observations rather than independent ones, so
#' \eqn{\nu} overstates the residual information available.
#'
#' The consequence for [rdm_sampling_covariance()] runs in *either*
#' direction and is not small. Under AR(1) errors with \eqn{\rho=0.75}, 32
#' trials, four conditions, six runs and 50 features, the ratio of the true
#' spread of the distance estimator to the standard error crossform reports
#' is 0.50 for a randomly interleaved trial order (the standard error is twice
#' too large) and 5.10 for a blocked order (five times too small). Supplying
#' `observation_whitener = L` with \eqn{L^\top L=\Sigma_t^{-1}} brings the same
#' blocked design to 1.03. crossform applies the \eqn{L} you give it and
#' records its identity; it cannot check that \eqn{L} matches the
#' autocorrelation actually present in your data.
#'
#' @param sources Raw observation-by-feature sources accepted by `relation()`.
#' @param design One design matrix, or one per partition.
#' @param effects One effect target matrix, or one per partition.
#' @param observation_whitener NULL, one observation whitener \eqn{L}, or one
#'   per partition. A finite square observation-by-observation matrix; the fit
#'   is carried out on \eqn{LX} and \eqn{Ly}, so \eqn{L} should satisfy
#'   \eqn{L^\top L=\Sigma_t^{-1}} for the within-partition error covariance
#'   \eqn{\Sigma_t}. Leaving it `NULL` asserts that the observations are
#'   already independent given the design; see the section above.
#' @param effect_names One common effect space/name vector, or one per partition.
#' @param tolerance Positive rank and estimability tolerance.
#' @param source_dims,partitions,domain,domain_id,capabilities Passed to
#'   `relation()`.
#' @param sampling_unit One nonempty sampling-unit label, or one per partition.
#' @param provenance Compact estimator provenance.
#' @param whiten Deprecated alias for `observation_whitener`.
#' @param solver One numerical route, or one per partition: automatic, QR, or
#'   SVD. The route is recorded in estimator provenance.
#' @return An `effect_relation_fit` whose `$relation` carries the fitted effect
#'   map and whose `$error_models` carry, per partition, the unscaled
#'   `effect_covariance`, `residual_df`, the lazy `residual_source`, the
#'   `observation_whitener` descriptor, and `estimator_provenance`. Every
#'   partition reports `within_participant_calibration`.
#' @family relation planning and fitting
#' @seealso [relation()] for precomputed effects with no error channel,
#'   [estimate_relation()] to reach this object from a validated
#'   [plan_relation()], and [residual_block()], [effect_covariance()],
#'   [residual_df()] to read the error channel.
#' @examples
#' set.seed(20260815)
#' domain <- abstract_domain(6L, id = "lm-relation-fit-example")
#' condition <- factor(rep(c("face", "body", "tool"), each = 4L))
#' design <- stats::model.matrix(~ 0 + condition)
#' colnames(design) <- c("face", "body", "tool")
#'
#' # Ask for the three condition means themselves.
#' targets <- diag(3)
#' dimnames(targets) <- list(colnames(design), colnames(design))
#' runs <- lapply(c("run-1", "run-2"), function(partition) {
#'   design %*% matrix(rnorm(18), 3L, 6L) + matrix(rnorm(72), 12L, 6L)
#' })
#' names(runs) <- c("run-1", "run-2")
#'
#' fit <- lm_relation_fit(
#'   runs, design, targets, domain = domain, sampling_unit = "trial"
#' )
#' relation_fit_capabilities(fit)$within_participant_calibration
#'
#' # Twelve trials minus three estimated means leaves nine residual df, which
#' # is what later analytic uncertainty is calibrated against.
#' residual_df(fit, "run-1")
#' @export
lm_relation_fit <- function(sources, design, effects,
                            observation_whitener = NULL,
                            effect_names = if (is.matrix(effects)) {
                              rownames(effects)
                            } else {
                              NULL
                            },
                            tolerance = sqrt(.Machine$double.eps),
                            source_dims = NULL, partitions = NULL,
                            domain = NULL, domain_id = "abstract",
                            capabilities = NULL,
                            sampling_unit = "observation",
                            provenance = list(), whiten = NULL,
                            solver = "auto") {
  if (is.matrix(sources)) sources <- list(sources)
  if (!is.list(sources) || length(sources) < 1L) {
    stop("`sources` must be a matrix or nonempty list.", call. = FALSE)
  }
  if (is.null(partitions)) partitions <- names(sources)
  if (is.null(partitions) || any(!nzchar(partitions))) {
    partitions <- paste0("partition", seq_along(sources))
  }
  if (!is.character(partitions) || length(partitions) != length(sources) ||
      anyNA(partitions) || any(!nzchar(partitions)) || anyDuplicated(partitions)) {
    stop("Partitions must have unique nonempty names.", call. = FALSE)
  }
  names(sources) <- partitions
  if (!missing(whiten) && !is.null(whiten)) {
    if (!is.null(observation_whitener)) {
      stop("Supply only `observation_whitener`; `whiten` is its deprecated alias.",
        call. = FALSE)
    }
    warning("`whiten` is deprecated; use `observation_whitener`.",
      call. = FALSE)
    observation_whitener <- whiten
  }
  if (missing(design)) {
    stop(paste0(
      "`design` is required: pass the observation-by-coefficient design ",
      "matrix, or a named list with one design per partition when runs ",
      "differ in length."
    ), call. = FALSE)
  }
  if (missing(effects)) {
    stop(paste0(
      "`effects` is required: pass the effect-by-coefficient matrix saying ",
      "which linear combination of design coefficients each experimental ",
      "effect is (`diag(ncol(design))` when the coefficients are the ",
      "effects)."
    ), call. = FALSE)
  }
  designs <- .partition_values(design, partitions, "design")
  targets <- .partition_values(effects, partitions, "effects")
  whiteners <- .partition_values(
    observation_whitener, partitions, "observation_whitener",
    allow_null = TRUE
  )
  names_by_partition <- .partition_values(
    effect_names, partitions, "effect_names", allow_null = TRUE
  )
  if (!is.character(sampling_unit) || anyNA(sampling_unit) ||
      any(!nzchar(sampling_unit)) ||
      !length(sampling_unit) %in% c(1L, length(partitions))) {
    stop("`sampling_unit` must provide one nonempty label per partition.",
      call. = FALSE)
  }
  if (length(sampling_unit) == 1L) {
    sampling_unit <- rep(sampling_unit, length(partitions))
  }
  names(sampling_unit) <- partitions
  if (!is.character(solver) || anyNA(solver) ||
      !length(solver) %in% c(1L, length(partitions)) ||
      any(!solver %in% c("auto", "qr", "svd"))) {
    stop("`solver` must provide `auto`, `qr`, or `svd` per partition.",
      call. = FALSE)
  }
  if (length(solver) == 1L) solver <- rep(solver, length(partitions))
  names(solver) <- partitions
  for (partition in partitions) {
    source_value <- sources[[partition]]
    design_value <- designs[[partition]]
    if (is.matrix(source_value) && is.matrix(design_value) &&
        nrow(source_value) != nrow(design_value)) {
      stop(sprintf(paste0(
        "Partition `%s` supplies %d observations but its design has %d ",
        "rows. Runs of unequal length are supported: pass `design` (and ",
        "any `observation_whitener`) as a named list with one entry per ",
        "partition."
      ), partition, nrow(source_value), nrow(design_value)), call. = FALSE)
    }
  }
  compiled <- lapply(partitions, function(partition) {
    design_value <- designs[[partition]]
    name_value <- names_by_partition[[partition]]
    if (is.null(name_value)) {
      name_value <- rownames(targets[[partition]])
    }
    resolved <- .resolve_observation_whitener(
      whiteners[[partition]], NULL, nrow(design_value),
      legacy_supplied = FALSE
    )
    .compile_lm_estimator(
      design_value, targets[[partition]], resolved,
      name_value, tolerance, partition = partition,
      solver = solver[[partition]]
    )
  })
  names(compiled) <- partitions
  if (any(vapply(compiled, function(value) value$residual_df < 1L,
    logical(1)))) {
    stop("Every fitted partition requires at least one residual degree of freedom.",
      call. = FALSE)
  }
  relation_value <- relation(
    sources, extract = lapply(compiled, `[[`, "extractor"),
    source_dims = source_dims, partitions = partitions,
    domain = domain, domain_id = domain_id, capabilities = capabilities,
    provenance = c(provenance, list(adapter = "lm_relation_fit"))
  )
  source_capabilities <- .relation_source_capabilities(relation_value)
  error_models <- lapply(partitions, function(partition) {
    source <- relation_value$sources[[partition]]
    estimate <- compiled[[partition]]
    source_revision <- source_capabilities[[partition]]$stable_revision
    residual_semantic <- list(
      source_revision = source_revision,
      estimator = estimate$estimator_provenance$signature,
      kind = "whitened_orthogonal_residuals"
    )
    residual_revision <- .sha256_signature(residual_semantic)
    residual_source <- .residual_source(
      source$dim,
      read = function(features) {
        estimate$residualize(source$read(features))
      },
      read_padded = function(features, width) {
        features <- .validate_source_features(features, source$dim[[2L]])
        if (!is.numeric(width) || length(width) != 1L || is.na(width) ||
            !is.finite(width) || width %% 1 != 0 ||
            width < length(features)) {
          stop("A padded residual width must be an integer at least as large as the feature block.",
            call. = FALSE)
        }
        raw <- source$read(features)
        if (!is.matrix(raw) || !is.numeric(raw) ||
            !identical(dim(raw), c(source$dim[[1L]], length(features))) ||
            any(!is.finite(raw))) {
          stop("Response source returned an invalid observation-by-feature block.",
            call. = FALSE)
        }
        padded <- matrix(0, source$dim[[1L]], as.integer(width))
        padded[, seq_along(features)] <- raw
        estimate$residualize(padded)
      },
      stable_revision = residual_revision
    )
    .separable_glm_error_model(
      effect_covariance = estimate$effect_covariance,
      residual_df = estimate$residual_df,
      sampling_unit = sampling_unit[[partition]],
      observation_whitener = estimate$observation_whitener,
      source_revision = source_revision,
      residual_source = residual_source,
      estimator_provenance = estimate$estimator_provenance,
      provenance = list(partition = partition)
    )
  })
  names(error_models) <- partitions
  relation_fit(
    relation_value, error_models,
    provenance = c(provenance, list(adapter = "lm_relation_fit"))
  )
}

#' Inspect statistical capabilities of a relation or relation fit
#'
#' Call this before requesting anything that needs residuals, so a missing
#' error channel surfaces as an explicit `FALSE` rather than a later refusal.
#'
#' @param x An `effect_relation` or `effect_relation_fit`.
#' @return A data frame with one row per partition: `partition` plus the
#'   logical columns `error_model`, `residual_blocks`, `effect_covariance`,
#'   `residual_df`, `separable_error`, `learned_metric_input`, and
#'   `within_participant_calibration`.
#' @family relation planning and fitting
#' @seealso [lm_relation_fit()] to obtain the capabilities,
#'   [sampling_capabilities()] for the uncertainty-side report, and
#'   [catch_refusal()] to inspect the refusal raised when one is missing.
#' @examples
#' example <- example_fmri_effects()
#'
#' # A fit built from raw responses carries the full error channel.
#' relation_fit_capabilities(example$fit)[
#'   , c("partition", "residual_blocks", "within_participant_calibration")
#' ]
#'
#' # The bare relation underneath reports every statistical capability FALSE:
#' # the point geometry is intact, the uncertainty channel is not there.
#' relation_fit_capabilities(example$fit$relation)[1L, ]
#' @export
relation_fit_capabilities <- function(x) {
  if (inherits(x, "effect_relation")) {
    .validate_relation(x, deep = FALSE)
    capabilities <- lapply(x$partitions, function(partition) {
      .error_capabilities(NULL)
    })
    partitions <- x$partitions
  } else {
    .validate_relation_fit(x, deep = FALSE)
    capabilities <- x$capabilities
    partitions <- x$relation$partitions
  }
  values <- do.call(rbind, lapply(capabilities, function(value) {
    as.data.frame(unclass(value), stringsAsFactors = FALSE)
  }))
  rownames(values) <- NULL
  data.frame(partition = partitions, values, check.names = FALSE,
    stringsAsFactors = FALSE)
}

.require_relation_fit_capability <- function(x, capability,
                                             partitions = NULL) {
  if (!is.character(capability) || length(capability) != 1L ||
      is.na(capability) || !nzchar(capability)) {
    stop("A capability requirement must be one nonempty identifier.",
      call. = FALSE)
  }
  table <- relation_fit_capabilities(x)
  if (!capability %in% names(table)) {
    stop(sprintf("Unknown relation-fit capability `%s`.", capability),
      call. = FALSE)
  }
  if (is.null(partitions)) partitions <- table$partition
  if (!is.character(partitions) || any(!partitions %in% table$partition)) {
    stop("Capability partitions must identify fitted relation partitions.",
      call. = FALSE)
  }
  selected <- match(partitions, table$partition)
  if (!all(table[[capability]][selected])) {
    if (inherits(x, "effect_relation") ||
        (inherits(x, "effect_relation_fit") &&
         all(vapply(x$error_models, is.null, logical(1))))) {
      .capability_refusal(paste0(
        "This relation was built from precomputed effects and has no ",
        "residual error channel. `", capability, "` cannot be recovered ",
        "from beta matrices alone. Start from raw responses with ",
        "`lm_relation_fit()` when this capability is required."
      ),
        capability = capability,
        namespace = "relation_fit",
        reasons = "missing_error_channel",
        remedies = "Refit raw observations with `lm_relation_fit()`."
      )
    }
    missing_partitions <- partitions[!table[[capability]][selected]]
    .capability_refusal(sprintf(
      paste0(
        "This operation requires relation-fit capability `%s` for every ",
        "partition; it is absent for: %s."
      ),
      capability, paste(missing_partitions, collapse = ", ")
    ),
      capability = capability,
      namespace = "relation_fit",
      reasons = paste0("missing_in_partition:", missing_partitions),
      remedies = paste0(
        "Refit the listed partitions with `lm_relation_fit()` so every ",
        "partition carries the capability."
      )
    )
  }
  invisible(TRUE)
}

.fit_partition <- function(x, partition) {
  partitions <- x$relation$partitions
  if (is.numeric(partition) && length(partition) == 1L && !is.na(partition) &&
      partition %% 1 == 0 && partition >= 1L && partition <= length(partitions)) {
    partition <- partitions[[as.integer(partition)]]
  }
  if (!is.character(partition) || length(partition) != 1L ||
      is.na(partition) || !partition %in% partitions) {
    stop("`partition` must identify one fitted relation partition.",
      call. = FALSE)
  }
  partition
}

#' Read fitted residuals for a neural feature block
#'
#' Residuals are produced lazily, one feature block at a time, without ever
#' materializing a dense observation residualizer. This is the input a learned
#' neural metric or an analytic RDM covariance draws on.
#'
#' @param x An `effect_relation_fit` with residual-block capability.
#' @param partition One partition name or index.
#' @param features Unique neural feature indices.
#' @return A whitened residual matrix with one row per observation in the
#'   partition and one column per requested feature.
#' @family relation planning and fitting
#' @seealso [lm_relation_fit()] which installs the residual channel,
#'   [residual_df()] for the matching degrees of freedom, and
#'   [noise_precision()], which learns a metric from these blocks.
#' @examples
#' example <- example_fmri_effects()
#'
#' # Whitened residuals for the first two neural features of one run.
#' residuals <- residual_block(example$fit, "run1", 1:2)
#' dim(residuals)
#'
#' # They are orthogonal to the fitted design, so the condition means have
#' # already been projected out.
#' round(colMeans(residuals), 8)
#'
#' # Divide by residual_df(), not nrow(), to estimate the noise variance.
#' round(colSums(residuals^2) / residual_df(example$fit, "run1"), 3)
#' @export
residual_block <- function(x, partition, features) {
  if (inherits(x, "effect_relation")) {
    .require_relation_fit_capability(x, "residual_blocks")
  }
  .validate_relation_fit(x, deep = FALSE)
  partition <- .fit_partition(x, partition)
  .require_relation_fit_capability(x, "residual_blocks", partition)
  features <- .validate_source_features(features, x$relation$n_features)
  source <- x$error_models[[partition]]$residual_source
  value <- source$read(features)
  if (!is.matrix(value) || !is.numeric(value) ||
      !identical(dim(value), c(source$dim[[1L]], length(features))) ||
      any(!is.finite(value))) {
    stop("Residual source returned an invalid observation-by-feature block.",
      call. = FALSE)
  }
  value
}

.residual_padded_block <- function(x, partition, features, width) {
  .validate_relation_fit(x, deep = FALSE)
  partition <- .fit_partition(x, partition)
  .require_relation_fit_capability(x, "residual_blocks", partition)
  features <- .validate_source_features(features, x$relation$n_features)
  if (!is.numeric(width) || length(width) != 1L || is.na(width) ||
      !is.finite(width) || width %% 1 != 0 || width < length(features)) {
    stop("`width` must be an integer at least as large as `features`.",
      call. = FALSE)
  }
  source <- x$error_models[[partition]]$residual_source
  value <- source$read_padded(features, as.integer(width))
  if (!is.matrix(value) || !is.numeric(value) ||
      !identical(dim(value), c(source$dim[[1L]], as.integer(width))) ||
      any(!is.finite(value))) {
    stop("Residual source returned an invalid fixed-width residual block.",
      call. = FALSE)
  }
  value
}

#' Read the unscaled effect-coordinate covariance of a fitted partition
#'
#' This is the design-side factor of the separable error model. By convention
#' the neural residual covariance is excluded, so the value depends on the
#' design and whitener but not on the data.
#'
#' @inheritParams residual_block
#' @return A symmetric effect-by-effect matrix, with rows and columns named by
#'   the relation's effect coordinates, excluding the neural residual
#'   covariance factor.
#' @family neural metrics
#' @seealso [residual_df()] and [residual_block()] for the other two pieces of
#'   the error channel, and [rdm_sampling_covariance()], which combines them.
#' @examples
#' example <- example_fmri_effects()
#'
#' # Four condition means estimated from eight trials each: the design factor
#' # is diagonal with entries 1/8, because the conditions are orthogonal.
#' covariance <- effect_covariance(example$fit, "run1")
#' round(covariance, 4)
#'
#' # Scaling it by a residual variance gives an effect standard error.
#' residuals <- residual_block(example$fit, "run1", 1L)
#' variance <- sum(residuals^2) / residual_df(example$fit, "run1")
#' round(sqrt(diag(covariance) * variance), 3)
#' @export
effect_covariance <- function(x, partition) {
  if (inherits(x, "effect_relation")) {
    .require_relation_fit_capability(x, "effect_covariance")
  }
  .validate_relation_fit(x, deep = FALSE)
  partition <- .fit_partition(x, partition)
  .require_relation_fit_capability(x, "effect_covariance", partition)
  x$error_models[[partition]]$effect_covariance
}

#' Read residual degrees of freedom from a fitted partition
#'
#' The divisor for any noise-variance estimate built from
#' [residual_block()]: observations minus the numerical rank of the whitened
#' design, not minus the number of design columns.
#'
#' @inheritParams residual_block
#' @return One positive integer.
#' @family relation planning and fitting
#' @seealso [residual_block()] and [effect_covariance()] for the rest of the
#'   error channel, and [relation_fit_capabilities()] to check availability.
#' @examples
#' example <- example_fmri_effects()
#'
#' # 32 trials per run minus the four estimated condition means.
#' residual_df(example$fit, "run1")
#'
#' # Each independent run contributes its own degrees of freedom.
#' vapply(example$fit$relation$partitions,
#'   function(partition) residual_df(example$fit, partition), integer(1))
#'
#' # A relation built from precomputed betas has none to report.
#' catch_refusal(residual_df(example$fit$relation, "run1"))$capability
#' @export
residual_df <- function(x, partition) {
  if (inherits(x, "effect_relation")) {
    .require_relation_fit_capability(x, "residual_df")
  }
  .validate_relation_fit(x, deep = FALSE)
  partition <- .fit_partition(x, partition)
  .require_relation_fit_capability(x, "residual_df", partition)
  x$error_models[[partition]]$residual_df
}
