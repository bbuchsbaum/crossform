# Same-space neural metrics -------------------------------------------------

.metric_capability_record <- function(identity, native_diagonal,
                                      positive_semidefinite,
                                      positive_definite, fixed,
                                      learned_frozen, learned_recipe,
                                      materialized, inverse_quadratic,
                                      inverse_quadratic_shortcut,
                                      inverse_quadratic_recipe = FALSE) {
  structure(list(
    role = "same_space_metric",
    identity = isTRUE(identity),
    native_diagonal = isTRUE(native_diagonal),
    feature_additive = isTRUE(native_diagonal),
    support_dense = !isTRUE(native_diagonal),
    positive_semidefinite = isTRUE(positive_semidefinite),
    positive_definite = isTRUE(positive_definite),
    fixed = isTRUE(fixed),
    learned_frozen = isTRUE(learned_frozen),
    learned_recipe = isTRUE(learned_recipe),
    materialized = isTRUE(materialized),
    inverse_quadratic = isTRUE(inverse_quadratic),
    inverse_quadratic_shortcut = isTRUE(inverse_quadratic_shortcut),
    inverse_quadratic_recipe = isTRUE(inverse_quadratic_recipe)
  ), class = "effect_metric_capabilities")
}

.metric_support <- function(support, domain, dimension) {
  domain <- .domain_reference(domain)
  if (is.null(support)) {
    if (dimension != domain$n_features) {
      .input_error(
        "A support-local metric requires explicit ordered feature identities."
      )
    }
    support <- domain$feature_ids
  }
  if (length(support) != dimension || anyNA(support) ||
      anyDuplicated(support)) {
    .input_error(
      "Metric support must uniquely identify every local matrix coordinate."
    )
  }
  positions <- match(support, domain$feature_ids)
  if (anyNA(positions) ||
      !identical(support, domain$feature_ids[positions])) {
    .input_error("Metric support must belong to its exact neural domain.")
  }
  list(domain = domain, support = support, positions = as.integer(positions))
}

.canonical_symmetric_metric <- function(value, tolerance, label) {
  matrix_like <- (is.matrix(value) && is.numeric(value)) ||
    inherits(value, "Matrix")
  if (!matrix_like || length(dim(value)) != 2L || any(dim(value) < 1L) ||
      nrow(value) != ncol(value) || any(!is.finite(value))) {
    .input_error(
      sprintf("%s must be one finite nonempty square matrix.", label)
    )
  }
  value <- as.matrix(value)
  storage.mode(value) <- "double"
  scale <- max(1, max(abs(value)))
  if (max(abs(value - t(value))) > tolerance * scale) {
    .input_error(
      sprintf("%s must be symmetric within the declared tolerance.", label)
    )
  }
  unname((value + t(value)) / 2)
}

.metric_inverse_representation <- function(inverse, value, tolerance) {
  if (is.null(inverse)) {
    semantic <- list(kind = "none")
    return(structure(c(semantic, list(
      signature = .sha256_signature(semantic)
    )), class = "effect_metric_inverse_representation"))
  }
  inverse <- .canonical_symmetric_metric(
    inverse, tolerance, "`inverse`"
  )
  if (!identical(dim(inverse), dim(value))) {
    .contract_error(
      "A retained inverse metric must match the metric dimensions."
    )
  }
  identity <- diag(nrow(value))
  residual <- max(
    abs(value %*% inverse - identity),
    abs(inverse %*% value - identity)
  )
  scale <- max(1, max(abs(value)) * max(abs(inverse)))
  if (!is.finite(residual) || residual > tolerance * scale * nrow(value)) {
    .contract_error(
      "The retained inverse is inconsistent with the neural metric."
    )
  }
  semantic <- list(kind = "retained_inverse_metric", value = inverse)
  structure(c(semantic, list(
    signature = .sha256_signature(semantic)
  )), class = "effect_metric_inverse_representation")
}

.metric_spectrum <- function(value, tolerance) {
  eigenvalues <- eigen(value, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(1, max(abs(eigenvalues)), max(abs(diag(value))))
  list(
    values = eigenvalues,
    scale = scale,
    positive_semidefinite = min(eigenvalues) >= -tolerance * scale,
    positive_definite = min(eigenvalues) > tolerance * scale
  )
}

#' Construct a support-local same-space neural metric
#'
#' A neural metric is the PSD same-space role of the evidence-pairing operator
#' `K`. It is distinct from a cross-space measurement bridge. The stored
#' matrix has width equal to its local support, never the full neural domain
#' unless the metric is genuinely global. Supplying `inverse` retains a
#' mathematically equivalent inverse action, such as the covariance from which
#' a precision metric was derived, so inverse quadratic forms need no second
#' factorization.
#'
#' @param value A finite symmetric positive-semidefinite local metric matrix.
#' @param domain Exact neural feature domain.
#' @param support Ordered feature identities for the local matrix coordinates.
#' @param inverse Optional finite symmetric inverse metric.
#' @param estimation Whether the materialized metric was fixed a priori or
#'   learned and frozen before evaluation.
#' @param tolerance Positive numerical validation tolerance.
#' @param provenance Compact metric provenance. Learned-frozen metrics require
#'   `frozen = TRUE` and a strong `training_signature`.
#' @return An `effect_neural_metric` carrying the canonicalized `$value`, its
#'   `$domain` and `$support`, the `$estimation` status, any retained
#'   `$inverse_representation`, the derived `$capabilities` record, and a
#'   content-addressed `$signature` bound into every plan that uses it.
#' @seealso [noise_precision()] for the semantically specific
#'   inverse-noise-covariance constructor, [metric_capabilities()] to inspect
#'   what a metric admits, and [plan_geometry()] to compile it into a plan.
#' @family neural metrics
#' @examples
#' # A fixed diagonal metric: each feature is reweighted independently, so
#' # the metric stays feature-additive and survives searchlight restriction.
#' domain <- abstract_domain(3, id = "metric-example")
#' metric <- neural_metric(diag(c(1, 2, 3)), domain)
#' metric_capabilities(metric)$feature_additive
#'
#' # Retaining the inverse lets inverse quadratic forms be read without a
#' # second factorization.
#' full <- neural_metric(
#'   matrix(c(2, 0.5, 0, 0.5, 2, 0.5, 0, 0.5, 2), 3, 3), domain,
#'   inverse = solve(matrix(c(2, 0.5, 0, 0.5, 2, 0.5, 0, 0.5, 2), 3, 3))
#' )
#' full$inverse_representation$kind
#'
#' # A same-space metric must be positive semidefinite.
#' indefinite <- try(
#'   neural_metric(diag(c(1, -1, 1)), domain), silent = TRUE
#' )
#' conditionMessage(attr(indefinite, "condition"))
#' @export
neural_metric <- function(value, domain, support = NULL, inverse = NULL,
                          estimation = c("fixed", "learned_frozen"),
                          tolerance = 1e-10, provenance = list()) {
  .check_number(tolerance, "tolerance", positive = TRUE)
  estimation <- match.arg(estimation)
  .validate_effect_provenance(provenance, "neural-metric provenance")
  value <- .canonical_symmetric_metric(value, tolerance, "`value`")
  identified <- .metric_support(support, domain, nrow(value))
  spectrum <- .metric_spectrum(value, tolerance)
  if (!spectrum$positive_semidefinite) {
    .input_error("A same-space neural metric must be positive semidefinite.")
  }
  if (estimation == "learned_frozen" &&
      (!identical(provenance$frozen, TRUE) ||
       !.strong_sha256(provenance$training_signature))) {
    .input_error(paste0(
      "Learned metrics require frozen training provenance and a strong ",
      "training signature."
    ))
  }
  inverse_representation <- .metric_inverse_representation(
    inverse, value, tolerance
  )
  if (!identical(inverse_representation$kind, "none") &&
      !spectrum$positive_definite) {
    .input_error(
      "A retained inverse requires a positive-definite neural metric."
    )
  }
  off_diagonal <- row(value) != col(value)
  native_diagonal <- all(value[off_diagonal] == 0)
  identity <- native_diagonal && all(diag(value) == 1)
  positive_diagonal <- native_diagonal &&
    all(diag(value) > tolerance * spectrum$scale)
  inverse_shortcut <- positive_diagonal ||
    identical(inverse_representation$kind, "retained_inverse_metric")
  capabilities <- .metric_capability_record(
    identity = identity,
    native_diagonal = native_diagonal,
    positive_semidefinite = spectrum$positive_semidefinite,
    positive_definite = spectrum$positive_definite,
    fixed = TRUE,
    learned_frozen = estimation == "learned_frozen",
    learned_recipe = FALSE,
    materialized = TRUE,
    inverse_quadratic = spectrum$positive_definite,
    inverse_quadratic_shortcut = inverse_shortcut
  )
  semantic <- list(
    schema_version = 1L,
    role = "same_space_metric",
    domain = identified$domain,
    support = identified$support,
    value = value,
    estimation = estimation,
    tolerance = tolerance,
    provenance = provenance
  )
  structure(list(
    role = semantic$role,
    domain = semantic$domain,
    support = semantic$support,
    positions = identified$positions,
    value = semantic$value,
    estimation = semantic$estimation,
    tolerance = semantic$tolerance,
    provenance = semantic$provenance,
    inverse_representation = inverse_representation,
    capabilities = capabilities,
    signature = .sha256_signature(semantic)
  ), class = "effect_neural_metric")
}

.validate_metric_capabilities <- function(x) {
  expected <- c("role", "identity", "native_diagonal", "feature_additive",
    "support_dense", "positive_semidefinite", "positive_definite", "fixed",
    "learned_frozen", "learned_recipe", "materialized",
    "inverse_quadratic", "inverse_quadratic_shortcut",
    "inverse_quadratic_recipe")
  if (!.sealed_fields(x, "effect_metric_capabilities", expected) ||
      !identical(x$role, "same_space_metric") ||
      !all(vapply(x[-1L], function(value) {     is.logical(value) && length(value) == 1L && !is.na(value) }, logical(1))) ||
      !identical(x$feature_additive, x$native_diagonal) ||
      identical(x$support_dense, x$native_diagonal)) {
    .input_error("Metric capabilities are missing or inconsistent.")
  }
  x
}

.validate_metric_inverse_representation <- function(x, value, tolerance,
                                                    deep = TRUE) {
  if (!inherits(x, "effect_metric_inverse_representation") || !is.list(x)) {
    .input_error("Metric inverse representation is invalid.")
  }
  expected <- if (identical(x$kind, "none")) {
    c("kind", "signature")
  } else {
    c("kind", "value", "signature")
  }
  if (!identical(names(x), expected) ||
      !x$kind %in% c("none", "retained_inverse_metric") ||
      !.strong_sha256(x$signature)) {
    .input_error("Metric inverse representation is invalid.")
  }
  if (isTRUE(deep)) {
    rebuilt <- .metric_inverse_representation(
      if (identical(x$kind, "none")) NULL else x$value,
      value, tolerance
    )
    if (!identical(x, rebuilt)) {
      .contract_error("Metric inverse representation is inconsistent.")
    }
  }
  x
}

.validate_neural_metric <- function(x, deep = TRUE) {
  if (.validated_before(x, "neural_metric", deep)) return(x)
  expected <- c("role", "domain", "support", "positions", "value",
    "estimation", "tolerance", "provenance", "inverse_representation",
    "capabilities", "signature")
  if (!.sealed_fields(x, "effect_neural_metric", expected) ||
      !identical(x$role, "same_space_metric") || !.is_finite_matrix(x$value) ||
      any(dim(x$value) < 1L) || nrow(x$value) != ncol(x$value) ||
      !.is_string(x$estimation, allow_empty = TRUE) ||
      !x$estimation %in% c("fixed", "learned_frozen") ||
      !.is_number(x$tolerance) || x$tolerance <= 0 || !is.list(x$provenance) ||
      !.strong_sha256(x$signature)) {
    .input_error("Neural-metric fields are missing or noncanonical.")
  }
  domain <- .validate_domain_reference(x$domain)
  identified <- .metric_support(x$support, domain, nrow(x$value))
  if (!identical(x$positions, identified$positions)) {
    .contract_error("Neural-metric support positions are inconsistent.")
  }
  inverse <- .validate_metric_inverse_representation(
    x$inverse_representation, x$value, x$tolerance, deep = deep
  )
  capabilities <- .validate_metric_capabilities(x$capabilities)
  if (isTRUE(deep)) {
    rebuilt <- neural_metric(
      x$value, domain, x$support,
      inverse = if (identical(inverse$kind, "none")) NULL else inverse$value,
      estimation = x$estimation,
      tolerance = x$tolerance,
      provenance = x$provenance
    )
    if (!identical(x, rebuilt) || !identical(capabilities, rebuilt$capabilities)) {
      .contract_error("Neural-metric identity is inconsistent.")
    }
  }
  .record_validated(x, "neural_metric", deep)
  x
}

.metric_recipe <- function(kind, domain,
                           native_diagonal = FALSE,
                           positive_definite = TRUE,
                           inverse_quadratic_recipe = FALSE,
                           hyperparameters = list(), provenance = list()) {
  if (!.is_string(kind) || !.is_flag(native_diagonal) ||
      !.is_flag(positive_definite) || !.is_flag(inverse_quadratic_recipe) ||
      !is.list(hyperparameters)) {
    .input_error("Metric-recipe metadata are invalid.")
  }
  domain <- if (is.null(domain)) NULL else .domain_reference(domain)
  .validate_effect_provenance(hyperparameters, "metric-recipe hyperparameters")
  .validate_effect_provenance(provenance, "metric-recipe provenance")
  capabilities <- .metric_capability_record(
    identity = FALSE,
    native_diagonal = native_diagonal,
    positive_semidefinite = TRUE,
    positive_definite = positive_definite,
    fixed = FALSE,
    learned_frozen = FALSE,
    learned_recipe = TRUE,
    materialized = FALSE,
    inverse_quadratic = FALSE,
    inverse_quadratic_shortcut = FALSE,
    inverse_quadratic_recipe = inverse_quadratic_recipe
  )
  semantic <- list(
    schema_version = 1L,
    role = "same_space_metric",
    kind = kind,
    domain = domain,
    support_local = TRUE,
    capabilities = unclass(capabilities),
    hyperparameters = hyperparameters,
    provenance = provenance
  )
  structure(list(
    role = semantic$role,
    kind = kind,
    domain = domain,
    support_local = TRUE,
    capabilities = capabilities,
    hyperparameters = hyperparameters,
    provenance = provenance,
    signature = .sha256_signature(semantic)
  ), class = "effect_metric_recipe")
}

.validate_metric_recipe <- function(x) {
  expected <- c("role", "kind", "domain", "support_local", "capabilities",
    "hyperparameters", "provenance", "signature")
  if (!.sealed_fields(x, "effect_metric_recipe", expected) ||
      !identical(x$role, "same_space_metric") ||
      !identical(x$support_local, TRUE) || !.strong_sha256(x$signature)) {
    .input_error("Metric-recipe fields are missing or noncanonical.")
  }
  if (!is.null(x$domain)) .validate_domain_reference(x$domain)
  rebuilt <- .metric_recipe(
    x$kind, x$domain,
    native_diagonal = x$capabilities$native_diagonal,
    positive_definite = x$capabilities$positive_definite,
    inverse_quadratic_recipe = x$capabilities$inverse_quadratic_recipe,
    hyperparameters = x$hyperparameters,
    provenance = x$provenance
  )
  if (!identical(x, rebuilt)) {
    .contract_error("Metric-recipe identity is inconsistent.")
  }
  x
}

#' Inspect exact neural-metric capabilities
#'
#' `metric_capabilities()` answers what a metric or recipe admits before a
#' plan provokes a refusal: whether it is diagonal, feature-additive,
#' positive definite, already materialized, and whether inverse quadratic
#' forms are available.
#'
#' @param x A `neural_metric()` or an on-demand metric recipe. A cross-space
#'   measurement bridge is refused because cross-space bridges are not
#'   same-space metrics.
#' @return An `effect_metric_capabilities` record of logical flags, including
#'   `$identity`, `$native_diagonal`, `$feature_additive`, `$support_dense`,
#'   `$positive_definite`, `$fixed`, `$learned_recipe`, `$materialized`, and
#'   `$inverse_quadratic`.
#' @seealso [neural_metric()] and [noise_precision()] for fixed metrics;
#'   [identity_metric()], [diagonal_precision()], and [shrinkage_precision()]
#'   for on-demand recipes.
#' @family neural metrics
#' @examples
#' # A materialized diagonal metric is feature-additive, so it restricts
#' # cleanly to any searchlight support.
#' domain <- abstract_domain(3, id = "capability-example")
#' fixed <- metric_capabilities(neural_metric(diag(c(1, 2, 3)), domain))
#' fixed[c("native_diagonal", "feature_additive", "materialized")]
#'
#' # A recipe is not yet a matrix: it promises a local metric per support.
#' recipe <- metric_capabilities(shrinkage_precision(0.2))
#' recipe[c("learned_recipe", "materialized", "support_dense")]
#'
#' # Anything that is not a metric or a recipe is refused rather than
#' # coerced, so a bare matrix or a cross-space object cannot pass as one.
#' refused <- try(metric_capabilities(diag(c(1, 2, 3))), silent = TRUE)
#' conditionMessage(attr(refused, "condition"))
#' @export
metric_capabilities <- function(x) {
  if (inherits(x, "effect_measurement_bridge")) {
    .input_error(paste0(
      "A measurement bridge relates two neural spaces; it is not a ",
      "same-space PSD metric."
    ))
  }
  if (inherits(x, "effect_neural_metric")) {
    return(.validate_neural_metric(x, deep = FALSE)$capabilities)
  }
  if (inherits(x, "effect_metric_recipe")) {
    return(.validate_metric_recipe(x)$capabilities)
  }
  .input_error("`x` must be a neural metric or metric recipe.")
}

.metric_handle <- function(metric) {
  metric <- .validate_neural_metric(metric, deep = FALSE)
  capabilities <- metric$capabilities
  state <- new.env(parent = emptyenv())
  state$factorizations <- 0L
  apply_K <- function(value) {
    if (is.atomic(value) && is.null(dim(value))) {
      if (!.is_finite_numeric(value) || length(value) != nrow(metric$value)) {
        .input_error("Metric application requires one finite local vector.")
      }
      return(drop(metric$value %*% value))
    }
    if (!.is_finite_matrix(value) || nrow(value) != nrow(metric$value)) {
      .input_error(
        "Metric application requires finite local vectors in columns."
      )
    }
    metric$value %*% value
  }
  inverse <- metric$inverse_representation
  if (identical(inverse$kind, "retained_inverse_metric")) {
    inverse_mode <- "retained_inverse_metric"
    quadform_Kinv <- function(value) {
      if (!.is_finite_numeric(value) || is.matrix(value) ||
          length(value) != nrow(metric$value)) {
        .input_error("Inverse quadratic forms require one finite local vector.")
      }
      drop(crossprod(value, inverse$value %*% value))
    }
  } else if (capabilities$native_diagonal &&
      capabilities$positive_definite) {
    inverse_mode <- "diagonal_reciprocal"
    diagonal_inverse <- 1 / diag(metric$value)
    quadform_Kinv <- function(value) {
      if (!.is_finite_numeric(value) || is.matrix(value) ||
          length(value) != nrow(metric$value)) {
        .input_error("Inverse quadratic forms require one finite local vector.")
      }
      sum(value^2 * diagonal_inverse)
    }
  } else if (capabilities$positive_definite) {
    inverse_mode <- "cholesky_solve"
    factor <- chol(metric$value)
    state$factorizations <- state$factorizations + 1L
    quadform_Kinv <- function(value) {
      if (!.is_finite_numeric(value) || is.matrix(value) ||
          length(value) != nrow(metric$value)) {
        .input_error("Inverse quadratic forms require one finite local vector.")
      }
      solved <- backsolve(factor, forwardsolve(t(factor), value))
      sum(value * solved)
    }
  } else {
    inverse_mode <- "unavailable_singular_metric"
    quadform_Kinv <- function(value) {
      .input_error(paste0(
        "This metric has no positive-definite inverse-quadratic capability."
      ))
    }
  }
  structure(list(
    metric = metric$signature,
    apply_K = apply_K,
    quadform_Kinv = quadform_Kinv,
    diagnostics = function() list(
      inverse_mode = inverse_mode,
      factorizations = state$factorizations,
      shortcut = inverse_mode %in% c(
        "retained_inverse_metric", "diagonal_reciprocal"
      )
    )
  ), class = "effect_metric_handle")
}

.metric_lowering <- function(metric) {
  capabilities <- metric_capabilities(metric)
  if (capabilities$learned_recipe) {
    return(if (capabilities$feature_additive) {
      "derive_then_additive_contraction"
    } else {
      "derive_then_support_streamed_pair_contraction"
    })
  }
  if (!capabilities$fixed || !capabilities$materialized) {
    .input_error(
      "Only fixed materialized metrics can enter an evidence lowering."
    )
  }
  if (capabilities$feature_additive) {
    "additive_contraction"
  } else {
    "support_streamed_pair_contraction"
  }
}

# The field layout of a geometry metric schedule depends on its kind: the
# learned kind carries a frozen `$schedule` where the other two carry a
# materialized `$metric`. Stated once so the constructors in the plan layer
# and the validator below cannot disagree about the seal.
.geometry_metric_schedule_fields <- function(kind) {
  base <- c("role", "kind", "frame_composition", "feature_additive",
    "support_dense", "materialization", "scope", "lowering",
    "metric_signature", "metric")
  if (identical(kind, "learned_local_before_frame")) {
    return(c(base, "schedule", "signature"))
  }
  c(base, "signature")
}

# The semantic digest of a schedule. The frozen schedule enters by its own
# signature rather than by value: the payload is the whole residual
# sufficient-statistics block, and its signature already digests the recipe,
# the statistics, the support index, the pairing, the training policy, and
# every per-edge training record.
.geometry_metric_schedule_semantic <- function(x) {
  semantic <- c(list(schema_version = 1L), unclass(x[
    !names(x) %in% c("metric", "signature")
  ]))
  if (identical(x$kind, "learned_local_before_frame")) {
    semantic$schedule <- x$schedule$signature
  }
  semantic
}

# The geometry metric schedule is built in the plan layer, but it is a
# statement about a metric, so its validator lives beside the metric it
# constrains and the plan calls down into it.
.validate_geometry_metric_schedule <- function(x, deep = TRUE) {
  expected <- .geometry_metric_schedule_fields(x$kind)
  if (!.sealed_fields(x, "effect_metric_schedule", expected) ||
      !identical(x$role, "same_space_metric_schedule") ||
      !x$kind %in% c("implicit_identity_before_frame",
        "fixed_metric_before_frame", "learned_local_before_frame") ||
      !identical(x$frame_composition, "sqrt_weight_congruence") ||
      !.is_flag(x$feature_additive) || !.is_flag(x$support_dense) ||
      identical(x$feature_additive, x$support_dense) ||
      !x$materialization %in% c("implicit", "fixed_metric",
        "on_demand_local") ||
      !x$scope %in% c("domain_operator", "single_node", "support_local") ||
      !x$lowering %in% c("additive_contraction",
        "support_streamed_pair_contraction",
        "derive_then_support_streamed_pair_contraction") ||
      !.strong_sha256(x$signature)) {
    .input_error("Geometry metric-schedule fields are missing or noncanonical.")
  }
  if (identical(x$kind, "implicit_identity_before_frame")) {
    if (!identical(x$feature_additive, TRUE) ||
        !identical(x$support_dense, FALSE) ||
        !identical(x$materialization, "implicit") ||
        !identical(x$lowering, "additive_contraction") ||
        !is.null(x$metric_signature) || !is.null(x$metric)) {
      .input_error("The implicit identity metric schedule is inconsistent.")
    }
  } else if (identical(x$kind, "learned_local_before_frame")) {
    # A learned local metric is derived per support and per evaluation edge,
    # so the schedule is not feature additive and needs the dense local
    # support even when the recipe itself is diagonal: the declaration is
    # about the schedule's admitted lowering, not about the recipe's shape.
    #
    # The frozen schedule is checked structurally here, not by calling
    # `.validate_frozen_metric_schedule()`: that validator lives in
    # metric-learning.R, which sits above this file in the value order
    # (metric-learning consumes recipes this file defines), and calling it
    # from here closed a four-file cycle. Deep validation of the frozen
    # record happens where the record is made and where it is spent -- at
    # construction (`.geometry_learned_metric_schedule()`,
    # `compile_metric_schedule()`), at plan validation
    # (`.validate_geometry_plan()`), and again in the kernel before any
    # contraction -- and the signature check below binds
    # `x$schedule$signature` into this schedule's own sealed identity.
    schedule <- x$schedule
    if (!inherits(schedule, "effect_frozen_metric_schedule") ||
        !.strong_sha256(schedule$signature)) {
      .contract_error(
        "A learned metric schedule must carry a frozen metric schedule."
      )
    }
    if (!identical(x$materialization, "on_demand_local") ||
        !identical(x$scope, "support_local") ||
        !identical(x$feature_additive, FALSE) ||
        !identical(x$support_dense, TRUE) ||
        !identical(x$lowering,
          "derive_then_support_streamed_pair_contraction") ||
        !is.null(x$metric_signature) || !is.null(x$metric) ||
        !isTRUE(schedule$capabilities$provenance_frozen) ||
        !identical(schedule$capabilities$materialized, FALSE)) {
      .contract_error("The learned local metric schedule is inconsistent.")
    }
    # Risk 1 of design section 16.7: a training assignment that varied with
    # the spatial support would be a different estimator per node, because
    # the compiler makes one tile one support. The records are frozen per
    # evaluation edge before any tiling exists; assert the shape that makes
    # that true rather than trusting the construction order.
    if (length(schedule$records) != nrow(schedule$pairing)) {
      .contract_error(
        "Learned metric training records must be one per evaluation edge."
      )
    }
  } else {
    # The `metric_signature` comparison below binds the metric to the schedule
    # identity at either depth; only the metric-internal rebuild is gated. That
    # rebuild is O(p^3) in a domain-wide metric, and every plan has already
    # passed it at its compile boundary.
    metric <- .validate_neural_metric(x$metric, deep = deep)
    if (!identical(x$materialization, "fixed_metric") ||
        !identical(x$metric_signature, metric$signature) ||
        !identical(x$feature_additive,
          metric$capabilities$feature_additive) ||
        !identical(x$support_dense, metric$capabilities$support_dense) ||
        !identical(x$lowering, .metric_lowering(metric))) {
      .contract_error("The fixed neural metric schedule is inconsistent.")
    }
  }
  expected_signature <- .sha256_signature(
    .geometry_metric_schedule_semantic(x)
  )
  .check_signature(
    x$signature, expected_signature,
    "Geometry metric-schedule identity is inconsistent."
  )
  x
}

.metric_additive_frame <- function(frame, metric_schedule) {
  .validate_frame_for_compile(frame)
  schedule <- .validate_geometry_metric_schedule(metric_schedule)
  if (!isTRUE(schedule$feature_additive)) {
    .input_error(
      "Only a compiler-proven native-diagonal metric is feature additive."
    )
  }
  if (identical(schedule$kind, "implicit_identity_before_frame")) {
    return(frame)
  }
  metric <- .validate_neural_metric(schedule$metric, deep = FALSE)
  diagonal <- numeric(frame$domain$n_features)
  diagonal[metric$positions] <- diag(metric$value)
  effective <- frame$weights %*% Matrix::Diagonal(x = diagonal)
  # The fold scales column `v` by `diagonal[v]`, so the composed weights no
  # longer satisfy the declared column-mass ("conservative") or row-mass
  # ("local") law, and declaring one would make the frame validator refuse a
  # numerically correct operator. The composed frame therefore declares
  # `normalization = "none"` -- the truth about the weights it carries -- and
  # records the declared normalization as provenance under `$metric_folded`,
  # so a reader can still tell what the weights were before composition and
  # that they are post-composition now.
  #
  # Conservation survives the fold. For a conservative frame
  # `sum_x w_xv d_v = d_v`, which is exactly the whole-support mass of feature
  # `v` read under the same diagonal metric, so `sum_x G_x = G_Omega` still
  # holds against the metric-folded global comparator. `reference_mass` is
  # that per-feature target, which is what `frame_conservation()` and
  # `.metric_frame_conservation()` certify against.
  folded <- additive_frame(effective, normalization = "none",
    domain = frame$domain)
  folded$metric_folded <- list(
    folded = TRUE,
    declared_normalization = frame$normalization,
    metric_kind = if (isTRUE(metric$capabilities$native_diagonal)) {
      "native_diagonal"
    } else {
      "feature_additive"
    },
    metric_signature = metric$signature,
    schedule_kind = schedule$kind,
    composition = "diagonal_metric_fold",
    reference_mass = diagonal
  )
  .validate_frame_for_compile(folded)
  folded
}

#' Identify an oriented coherent neural functional
#'
#' The vector `a` defines the scientifically oriented amplitude `B a`. It is
#' separate from the metric that normalizes the energy of that amplitude.
#'
#' @param value One finite nonzero local covector.
#' @param domain Exact neural feature domain.
#' @param support Ordered feature identities for `value`.
#' @param label Nonempty scientific identity for the functional.
#' @param provenance Compact provenance.
#' @return An immutable `effect_coherent_functional`.
#' @keywords internal
coherent_functional <- function(value, domain, support = NULL,
                                label = "raw_weighted_mean",
                                provenance = list()) {
  if (!.is_finite_numeric(value) || is.matrix(value) || length(value) < 1L ||
      all(value == 0)) {
    .input_error("A coherent functional must be one finite nonzero vector.")
  }
  if (!.is_string(label)) {
    .input_error("A coherent functional requires one nonempty label.")
  }
  .validate_effect_provenance(provenance,
    "coherent-functional provenance")
  identified <- .metric_support(support, domain, length(value))
  semantic <- list(
    schema_version = 1L,
    domain = identified$domain,
    support = identified$support,
    value = unname(as.double(value)),
    label = label,
    orientation = "fixed",
    provenance = provenance
  )
  structure(c(semantic[c("domain", "support", "value", "label",
    "orientation", "provenance")], list(
    signature = .sha256_signature(semantic)
  )), class = "effect_coherent_functional")
}

.validate_coherent_functional <- function(x) {
  expected <- c("domain", "support", "value", "label", "orientation",
    "provenance", "signature")
  if (!.sealed_fields(x, "effect_coherent_functional", expected) ||
      !identical(x$orientation, "fixed") || !.strong_sha256(x$signature)) {
    .input_error("Coherent-functional fields are missing or noncanonical.")
  }
  rebuilt <- coherent_functional(
    x$value, x$domain, x$support, x$label, x$provenance
  )
  if (!identical(x, rebuilt)) {
    .contract_error("Coherent-functional identity is inconsistent.")
  }
  x
}

# The two frame-metric node readers below -- the validating one and the
# trusted accessor the kernels use -- open the same way and close the same
# way. Both admit the same frame and derive the same measurement labels, and
# both report the same five fields for a node. Only the middle differs: one
# indexes a single row of whatever the frame happens to be, the other pays
# once for an Rsparse view and then does CSR slot arithmetic per node. These
# two helpers are that shared opening and closing.
.additive_frame_node_ids <- function(frame) {
  .validate_frame_for_compile(frame)
  if (!identical(frame$representation, "additive_diagonal")) {
    .input_error(
      "Frame-metric composition requires an additive localization frame."
    )
  }
  if (!is.null(frame$index) && "measurement" %in% names(frame$index)) {
    frame$index$measurement
  } else {
    seq_len(nrow(frame$weights))
  }
}

.frame_metric_node_value <- function(frame, node_ids, position,
                                     support_positions, weight) {
  list(
    position = position,
    id = node_ids[[position]],
    weight = weight,
    support_positions = as.integer(support_positions),
    support = frame$domain$feature_ids[support_positions]
  )
}

# The dense fallback for a frame whose weights are an ordinary matrix: the
# support is wherever the row is strictly positive.
.dense_frame_row_support <- function(row) {
  dense <- as.numeric(row)
  support_positions <- which(dense > 0)
  list(support_positions = support_positions, weight = dense[support_positions])
}

.frame_metric_node <- function(frame, node) {
  node_ids <- .additive_frame_node_ids(frame)
  if (is.numeric(node) && length(node) == 1L && !is.na(node) &&
      node %% 1 == 0 && node >= 1L && node <= nrow(frame$weights)) {
    position <- as.integer(node)
  } else {
    position <- match(node, node_ids)
    if (length(position) != 1L || is.na(position)) {
      .input_error("`node` must identify one additive-frame measurement.")
    }
  }
  row <- frame$weights[position, , drop = FALSE]
  if (inherits(row, "sparseMatrix")) {
    entries <- Matrix::summary(row)
    keep <- entries$x > 0
    support_positions <- as.integer(entries$j[keep])
    weight <- as.numeric(entries$x[keep])
    ordering <- order(support_positions)
    support_positions <- support_positions[ordering]
    weight <- weight[ordering]
  } else {
    dense <- .dense_frame_row_support(row)
    support_positions <- dense$support_positions
    weight <- dense$weight
  }
  .frame_metric_node_value(
    frame, node_ids, position, support_positions, weight
  )
}

# Trusted compiled-frame accessor for numerical kernels. Public and compiler
# boundaries validate the frame once; this function performs only local CSR
# indexing. Searchlight frames use the support index as the authoritative
# support order, while row weights are read from a one-time Rsparse view.
.frame_metric_node_accessor <- function(frame, argument = "at") {
  node_ids <- .additive_frame_node_ids(frame)
  row_weights <- if (inherits(frame$weights, "sparseMatrix")) {
    methods::as(frame$weights, "RsparseMatrix")
  } else {
    NULL
  }
  support_index <- frame$support_index
  function(position) {
    position <- .msg_measurement_index(
      position, nrow(frame$weights), argument = argument, subject = "frame"
    )
    if (!is.null(support_index)) {
      support_positions <- .support_index_support_trusted(
        support_index, position
      )
      first <- row_weights@p[[position]] + 1L
      last <- row_weights@p[[position + 1L]]
      slots <- seq.int(first, last)
      columns <- row_weights@j[slots] + 1L
      weights <- row_weights@x[slots]
      matched <- match(support_positions, columns)
      if (anyNA(matched)) {
        .contract_error(
          "Compiled frame weights do not match their support index."
        )
      }
      weight <- as.numeric(weights[matched])
    } else if (!is.null(row_weights)) {
      first <- row_weights@p[[position]] + 1L
      last <- row_weights@p[[position + 1L]]
      slots <- seq.int(first, last)
      support_positions <- as.integer(row_weights@j[slots] + 1L)
      weight <- as.numeric(row_weights@x[slots])
    } else {
      dense <- .dense_frame_row_support(frame$weights[position, ])
      support_positions <- dense$support_positions
      weight <- dense$weight
    }
    .frame_metric_node_value(
      frame, node_ids, position, support_positions, weight
    )
  }
}

.compose_frame_metric <- function(frame, metric = NULL, node) {
  node_value <- .frame_metric_node(frame, node)
  frame_signature <- .additive_frame_signature(frame)
  if (is.null(metric)) {
    metric <- neural_metric(
      diag(length(node_value$support)), frame$domain,
      support = node_value$support,
      provenance = list(construction = "identity_local_metric")
    )
  }
  metric <- .validate_neural_metric(metric)
  if (!.same_domain_reference(metric$domain, frame$domain) ||
      !identical(metric$support, node_value$support)) {
    .contract_error(
      "The local metric support must exactly match its frame node."
    )
  }
  root_weight <- sqrt(node_value$weight)
  effective_value <- if (metric$capabilities$native_diagonal) {
    # This is the exact diagonal specialization of
    # D(sqrt(w)) K D(sqrt(w)).  Besides avoiding a dense outer product, it
    # preserves the existing additive weights without a sqrt/square round trip.
    diag(node_value$weight * diag(metric$value))
  } else {
    metric$value * tcrossprod(root_weight)
  }
  inverse <- metric$inverse_representation
  effective_inverse <- if (identical(inverse$kind,
      "retained_inverse_metric")) {
    inverse$value * tcrossprod(1 / root_weight)
  } else {
    NULL
  }
  provenance <- list(
    construction = "sqrt_weight_congruence",
    law = "D(sqrt(w)) K D(sqrt(w))",
    frame = frame_signature,
    node = as.character(node_value$id),
    base_metric = metric$signature,
    base_provenance = metric$provenance
  )
  if (metric$estimation == "learned_frozen") {
    provenance$frozen <- TRUE
    provenance$training_signature <- metric$provenance$training_signature
  }
  effective <- neural_metric(
    effective_value, frame$domain, node_value$support,
    inverse = effective_inverse,
    estimation = metric$estimation,
    tolerance = metric$tolerance,
    provenance = provenance
  )
  coherent <- coherent_functional(
    node_value$weight / sum(node_value$weight),
    frame$domain, node_value$support,
    label = "normalized_raw_frame_weight",
    provenance = list(
      frame = frame_signature,
      node = as.character(node_value$id),
      metric_independent_orientation = TRUE
    )
  )
  semantic <- list(
    schema_version = 1L,
    composition = "sqrt_weight_congruence",
    frame = frame_signature,
    node = as.character(node_value$id),
    base_metric = metric$signature,
    effective_metric = effective$signature,
    coherent_functional = coherent$signature
  )
  structure(list(
    composition = "sqrt_weight_congruence",
    frame = frame_signature,
    node = node_value$id,
    weight = node_value$weight,
    base_metric = metric,
    metric = effective,
    coherent = coherent,
    signature = .sha256_signature(semantic)
  ), class = "effect_composed_metric")
}

#' Split a neural metric into coherent and configuration components
#'
#' For SPD `K` and an oriented covector `a`, the exact metric-aware split is
#' `K_coh = a a^T / (a^T K^-1 a)` and `K_cfg = K - K_coh`. The coherent
#' amplitude remains `B a`; the metric controls only its normalization.
#'
#' @param metric A materialized `neural_metric()` or an internal composed
#'   frame-metric value.
#' @param coherent A matching `coherent_functional()` or numeric vector. It may
#'   be omitted for a composed frame metric.
#' @return An `effect_metric_components` value.
#' @examples
#' domain <- abstract_domain(3, id = "component-example")
#' metric <- neural_metric(diag(c(1, 2, 3)), domain)
#' mean_functional <- crossform:::coherent_functional(rep(1 / 3, 3), domain)
#' components <- crossform:::metric_components(metric, mean_functional)
#' all.equal(components$coherent + components$configuration, metric$value)
#' @keywords internal
metric_components <- function(metric, coherent = NULL) {
  if (inherits(metric, "effect_measurement_bridge")) {
    .input_error(
      "Cross-space bridges do not admit same-space metric components."
    )
  }
  if (inherits(metric, "effect_metric_recipe")) {
    .input_error(paste0(
      "An on-demand metric recipe must be derived and frozen before ",
      "decomposition."
    ))
  }
  if (inherits(metric, "effect_composed_metric")) {
    if (is.null(coherent)) coherent <- metric$coherent
    metric <- metric$metric
  }
  metric <- .validate_neural_metric(metric)
  if (is.null(coherent)) {
    .input_error(
      "Metric decomposition requires an oriented coherent functional."
    )
  }
  if (is.numeric(coherent) && !is.matrix(coherent)) {
    coherent <- coherent_functional(
      coherent, metric$domain, metric$support
    )
  }
  coherent <- .validate_coherent_functional(coherent)
  if (!.same_domain_reference(metric$domain, coherent$domain) ||
      !identical(metric$support, coherent$support)) {
    .contract_error(
      "The metric and coherent functional must share one ordered support."
    )
  }
  if (!metric$capabilities$positive_definite ||
      !metric$capabilities$inverse_quadratic) {
    .input_error(paste0(
      "Coherent/configuration decomposition requires an SPD metric with ",
      "inverse-quadratic capability."
    ))
  }
  handle <- .metric_handle(metric)
  denominator <- handle$quadform_Kinv(coherent$value)
  if (!is.finite(denominator) || denominator <= 0) {
    .invariant_error(
      "The coherent inverse-metric norm must be positive and finite."
    )
  }
  coherent_metric <- tcrossprod(coherent$value) / denominator
  configuration_metric <- metric$value - coherent_metric
  tolerance <- metric$tolerance
  spectrum <- .metric_spectrum(configuration_metric, tolerance)
  if (!spectrum$positive_semidefinite) {
    .invariant_error("The derived configuration metric failed its PSD law.")
  }
  coherent_spectrum <- .metric_spectrum(coherent_metric, tolerance)
  coherent_rank <- sum(coherent_spectrum$values >
    tolerance * coherent_spectrum$scale)
  if (coherent_rank != 1L) {
    .invariant_error("The derived coherent metric failed its rank-one law.")
  }
  diagnostic <- handle$diagnostics()
  semantic <- list(
    schema_version = 1L,
    metric = metric$signature,
    coherent_functional = coherent$signature,
    coherent_metric = coherent_metric,
    configuration_metric = configuration_metric
  )
  structure(list(
    metric = metric,
    coherent_functional = coherent,
    coherent = coherent_metric,
    configuration = configuration_metric,
    denominator = denominator,
    coherent_rank = as.integer(coherent_rank),
    configuration_psd = TRUE,
    inverse_quadratic_mode = diagnostic$inverse_mode,
    factorization_count = diagnostic$factorizations,
    signature = .sha256_signature(semantic)
  ), class = "effect_metric_components")
}

.metric_frame_conservation <- function(frame, metrics = NULL,
                                       tolerance = 1e-10) {
  .validate_frame_for_compile(frame)
  nodes <- nrow(frame$weights)
  if (is.null(metrics)) metrics <- rep(list(NULL), nodes)
  if (!is.list(metrics) || length(metrics) != nodes) {
    .input_error(
      "Conservation analysis requires one base metric per frame node."
    )
  }
  composed <- lapply(seq_len(nodes), function(node) {
    .compose_frame_metric(frame, metrics[[node]], node)
  })
  additive <- all(vapply(composed, function(value) {
    value$metric$capabilities$feature_additive
  }, logical(1)))
  # The composed diagonals must sum to the mass the frame's columns are meant
  # to carry: one for a declared frame, and the folded metric diagonal for a
  # frame that already has a diagonal metric folded in, whose global
  # comparator is read under that same metric. This keeps the certificate
  # agreeing with `frame_conservation()` on the same frame.
  fold <- frame$metric_folded
  reference <- .frame_conservation_reference(frame)
  global_diagonal <- NULL
  identity_conservation <- FALSE
  if (additive) {
    global_diagonal <- numeric(frame$domain$n_features)
    for (value in composed) {
      global_diagonal[value$metric$positions] <-
        global_diagonal[value$metric$positions] + diag(value$metric$value)
    }
    identity_conservation <- max(abs(global_diagonal - reference)) <= tolerance
  }
  semantic <- list(
    schema_version = 1L,
    frame = .additive_frame_signature(frame),
    metrics = vapply(composed, function(value) {
      value$base_metric$signature
    }, character(1)),
    feature_additive = additive,
    global_diagonal = global_diagonal,
    reference_mass = reference,
    metric_folded = !is.null(fold),
    declared_normalization = if (is.null(fold)) {
      frame$normalization
    } else {
      fold$declared_normalization
    },
    identity_conservation = identity_conservation,
    tolerance = tolerance
  )
  structure(list(
    feature_additive = additive,
    global_metric_kind = if (additive) "native_diagonal" else
      "support_pair_operator",
    global_diagonal = global_diagonal,
    reference_mass = reference,
    metric_folded = !is.null(fold),
    declared_normalization = if (is.null(fold)) {
      frame$normalization
    } else {
      fold$declared_normalization
    },
    identity_conservation = identity_conservation,
    reason = if (!additive) {
      paste0(
        "Conservative frame weights alone do not conserve a non-diagonal ",
        "metric schedule."
      )
    } else if (!identity_conservation) {
      if (is.null(fold)) {
        "The summed diagonal metric is not the native identity metric."
      } else {
        "The summed diagonal metric is not the folded reference mass."
      }
    } else if (is.null(fold)) {
      "The composed feature-additive metrics resolve the native identity."
    } else {
      paste0(
        "The composed feature-additive metrics resolve the folded metric ",
        "diagonal."
      )
    },
    signature = .sha256_signature(semantic)
  ), class = "effect_metric_conservation")
}

.require_metric_conservation <- function(x,
                                         target = c("feature_additive",
                                                    "identity")) {
  target <- match.arg(target)
  if (!inherits(x, "effect_metric_conservation") || !is.list(x)) {
    .input_error("`x` must be a metric-conservation certificate.")
  }
  admitted <- if (target == "feature_additive") {
    x$feature_additive
  } else {
    x$identity_conservation
  }
  if (!isTRUE(admitted)) {
    .input_error(sprintf("Metric conservation is not certified: %s", x$reason))
  }
  invisible(TRUE)
}
