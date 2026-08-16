# Observation-model assumptions -------------------------------------------

.normalize_whitener_declaration <- function(value, kind) {
  if (identical(kind, "ols")) {
    if (!is.null(value)) {
      stop("`whitener` must be NULL for `kind = \"ols\"`.", call. = FALSE)
    }
    return(NULL)
  }
  if (is.matrix(value)) value <- list(value)
  if (!is.list(value) || length(value) < 1L) {
    stop("GLS observation models require at least one whitener matrix.",
      call. = FALSE)
  }
  for (index in seq_along(value)) {
    matrix <- value[[index]]
    if (!is.matrix(matrix) || !is.numeric(matrix) || nrow(matrix) < 1L ||
        nrow(matrix) != ncol(matrix) || any(!is.finite(matrix))) {
      stop("Every declared whitener must be a finite nonempty square matrix.",
        call. = FALSE)
    }
  }
  value
}

.whitener_revisions <- function(value) {
  if (is.null(value)) return(NULL)
  revisions <- vapply(value, function(matrix) {
    .semantic_digest("sha256:", unname(matrix))
  }, character(1))
  names(revisions) <- names(value)
  revisions
}

#' Declare the first-moment observation model
#'
#' States the estimation assumptions the data tables cannot prove: the error
#' law, the sampling unit, any whitener, and whether independence was asserted.
#' Nothing here is inferred from events or confounds, and the declaration is
#' what determines whether analytic uncertainty is later available.
#'
#' @param kind `"ols"`, `"fixed_gls"`, or `"learned_frozen_gls"`.
#' @param sampling_unit One sampling-unit axis name, or a named value per
#'   partition. It is required and is never inferred from events.
#' @param whitener For GLS, one finite square matrix or a named list per
#'   partition.
#' @param independence Optional explicit independence declaration. `NULL`
#'   records that no independence guarantee was asserted.
#' @param training_revision Strong training revision required for learned GLS.
#' @param training_provenance Portable learned-model provenance.
#' @param assumptions Additional portable declared assumptions.
#' @param provenance Portable model provenance.
#' @return An `effect_observation_model`: a list with `$kind`,
#'   `$sampling_unit`, the `$whitener` and its `$whitener_revisions`,
#'   `$independence`, learned-model `$training_revision` and
#'   `$training_provenance`, `$assumptions`, `$provenance`, an
#'   `$observation_model_id`, and `$capabilities` with
#'   `fixed_observation_model`, `learned_observation_model`,
#'   `separable_glm_law`, `analytic_effect_covariance`, and
#'   `declared_independence`.
#' @family studies and effect maps
#' @seealso [plan_relation()], which consumes this declaration, and
#'   [rdm_sampling_covariance()], whose analytic route requires
#'   `analytic_effect_covariance`.
#' @examples
#' # OLS with scans as the sampling unit. Independence must be declared; it is
#' # never read off the event table.
#' ols <- observation_model(
#'   "ols", sampling_unit = "scan",
#'   independence = "runs independently acquired conditional on the model"
#' )
#' ols$capabilities$analytic_effect_covariance
#' ols$capabilities$declared_independence
#'
#' # A frozen learned whitener still yields point estimates, but the analytic
#' # error channel is withheld because its training is not accounted for.
#' learned <- observation_model(
#'   "learned_frozen_gls", sampling_unit = "scan", whitener = diag(6),
#'   training_revision = paste0("sha256:", strrep("a", 64)),
#'   training_provenance = list(method = "AR1", fitted_on = "held-out runs")
#' )
#' learned$capabilities$analytic_effect_covariance
#' @export
observation_model <- function(
    kind = c("ols", "fixed_gls", "learned_frozen_gls"), sampling_unit,
    whitener = NULL, independence = NULL, training_revision = NULL,
    training_provenance = list(), assumptions = list(), provenance = list()) {
  kind <- match.arg(kind)
  if (missing(sampling_unit) || !is.character(sampling_unit) ||
      length(sampling_unit) < 1L || anyNA(sampling_unit) ||
      any(!nzchar(sampling_unit))) {
    stop("`sampling_unit` must be declared explicitly.", call. = FALSE)
  }
  if (!is.null(independence)) {
    independence <- .validate_nonempty_id(independence, "independence")
  }
  whitener <- .normalize_whitener_declaration(whitener, kind)
  learned <- identical(kind, "learned_frozen_gls")
  if (learned) {
    if (!.strong_sha256(training_revision)) {
      stop("Learned GLS requires one strong `training_revision`.",
        call. = FALSE)
    }
    if (!is.list(training_provenance) || !length(training_provenance)) {
      stop("Learned GLS requires nonempty `training_provenance`.",
        call. = FALSE)
    }
  } else if (!is.null(training_revision) || length(training_provenance)) {
    stop("Training provenance is only valid for learned GLS.", call. = FALSE)
  }
  .validate_effect_provenance(training_provenance, "training provenance")
  .validate_effect_provenance(assumptions, "observation-model assumptions")
  .validate_effect_provenance(provenance, "observation-model provenance")
  semantic <- list(
    schema_version = 1L,
    kind = kind,
    sampling_unit = unname(sampling_unit),
    sampling_unit_names = names(sampling_unit),
    whitener_revisions = .whitener_revisions(whitener),
    independence = independence,
    training_revision = training_revision,
    training_provenance = training_provenance,
    assumptions = assumptions,
    provenance = provenance
  )
  structure(list(
    kind = kind,
    sampling_unit = sampling_unit,
    whitener = whitener,
    whitener_revisions = semantic$whitener_revisions,
    independence = independence,
    training_revision = training_revision,
    training_provenance = training_provenance,
    assumptions = assumptions,
    provenance = provenance,
    observation_model_id = .semantic_digest(
      "observation-model-sha256:", semantic
    ),
    capabilities = list(
      fixed_observation_model = !learned,
      learned_observation_model = learned,
      separable_glm_law = !learned,
      analytic_effect_covariance = !learned,
      declared_independence = !is.null(independence)
    )
  ), class = "effect_observation_model")
}

.validate_observation_model <- function(value) {
  expected <- c(
    "kind", "sampling_unit", "whitener", "whitener_revisions",
    "independence", "training_revision", "training_provenance",
    "assumptions", "provenance", "observation_model_id", "capabilities"
  )
  if (!inherits(value, "effect_observation_model") || !is.list(value) ||
      !identical(names(value), expected)) {
    stop("Observation-model fields are missing or noncanonical.",
      call. = FALSE)
  }
  rebuilt <- observation_model(
    kind = value$kind,
    sampling_unit = value$sampling_unit,
    whitener = value$whitener,
    independence = value$independence,
    training_revision = value$training_revision,
    training_provenance = value$training_provenance,
    assumptions = value$assumptions,
    provenance = value$provenance
  )
  if (!identical(value, rebuilt)) {
    stop("Observation-model metadata or identity is inconsistent.",
      call. = FALSE)
  }
  rebuilt
}

.observation_model_for_partitions <- function(value, partitions, full_rows,
                                              retained_rows) {
  value <- .validate_observation_model(value)
  sampling <- value$sampling_unit
  if (length(sampling) == 1L) {
    sampling <- rep(sampling, length(partitions))
  } else if (!is.null(names(sampling))) {
    if (!setequal(names(sampling), partitions)) {
      stop("Named sampling units must cover every study partition.",
        call. = FALSE)
    }
    sampling <- sampling[partitions]
  }
  if (length(sampling) != length(partitions)) {
    stop("`sampling_unit` must provide one value per study partition.",
      call. = FALSE)
  }
  names(sampling) <- partitions

  whiteners <- value$whitener
  if (is.null(whiteners)) {
    whiteners <- stats::setNames(rep(list(NULL), length(partitions)), partitions)
  } else {
    if (length(whiteners) == 1L && length(partitions) > 1L &&
        is.null(names(whiteners))) {
      whiteners <- rep(whiteners, length(partitions))
    } else if (!is.null(names(whiteners))) {
      if (!setequal(names(whiteners), partitions)) {
        stop("Named whiteners must cover every study partition.",
          call. = FALSE)
      }
      whiteners <- whiteners[partitions]
    }
    if (length(whiteners) != length(partitions)) {
      stop("`whitener` must provide one matrix per study partition.",
        call. = FALSE)
    }
    names(whiteners) <- partitions
    for (partition in partitions) {
      matrix <- whiteners[[partition]]
      retained <- retained_rows[[partition]]
      if (nrow(matrix) == full_rows[[partition]]) {
        matrix <- matrix[retained, retained, drop = FALSE]
      } else if (nrow(matrix) != length(retained)) {
        .study_refusal(
          sprintf("Partition `%s` whitener does not match its observation rows.",
            partition),
          capability = "aligned_observations",
          reasons = sprintf(
            "Whitener dimension is %d; expected %d full or %d retained rows.",
            nrow(matrix), full_rows[[partition]], length(retained)
          ),
          remedies = "Supply a whitener on the complete or retained observation axis."
        )
      }
      whiteners[[partition]] <- matrix
    }
  }
  list(sampling_unit = sampling, whiteners = whiteners)
}
