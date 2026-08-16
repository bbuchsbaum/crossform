# Frozen-provenance, on-demand metric learning -----------------------------

.validate_metric_floor <- function(relative_floor, absolute_floor) {
  if (!is.numeric(relative_floor) || length(relative_floor) != 1L ||
      is.na(relative_floor) || !is.finite(relative_floor) ||
      relative_floor <= 0 || !is.numeric(absolute_floor) ||
      length(absolute_floor) != 1L || is.na(absolute_floor) ||
      !is.finite(absolute_floor) || absolute_floor < 0) {
    stop("Metric variance floors must be finite, with a positive relative floor and a nonnegative absolute floor.",
      call. = FALSE)
  }
  list(
    relative_variance_floor = as.double(relative_floor),
    absolute_variance_floor = as.double(absolute_floor)
  )
}

#' Specify an on-demand identity metric
#'
#' These constructors describe how a local neural metric is to be produced.
#' They do not allocate a domain-wide matrix. Use `compile_metric_schedule()`
#' to bind a recipe to residual sufficient statistics, spatial supports, and
#' evaluation edges.
#'
#' @param domain Optional exact neural domain. Omitting it defers domain binding
#'   until schedule compilation.
#' @return An `effect_metric_recipe` with `$kind`, the optional bound
#'   `$domain`, a `$capabilities` record, the fixed `$hyperparameters`
#'   (estimator, randomness, seed), and a `$signature`. It allocates no
#'   matrix.
#' @seealso [diagonal_precision()] and [shrinkage_precision()] for
#'   residual-derived recipes, [metric_capabilities()] to inspect one, and
#'   [plan_crossnobis()], which binds a recipe to residual statistics.
#' @family neural metrics
#' @examples
#' # A recipe is a declaration, not a matrix: no domain-wide operator exists
#' # until a schedule binds it to supports.
#' recipe <- identity_metric()
#' recipe$kind
#' recipe$hyperparameters$estimator
#' metric_capabilities(recipe)$materialized
#'
#' # Bind the domain up front when you want a domain mismatch caught before
#' # schedule compilation rather than during it.
#' domain <- abstract_domain(3, id = "identity-metric-example")
#' identity_metric(domain)$domain$id
#'
#' # Unlike the residual-derived recipes, this one estimates nothing, so it
#' # stays diagonal and needs no residual channel.
#' c(identity = metric_capabilities(recipe)$native_diagonal,
#'   shrinkage = metric_capabilities(shrinkage_precision())$native_diagonal)
#' @export
identity_metric <- function(domain = NULL) {
  .metric_recipe(
    "identity", domain,
    native_diagonal = TRUE,
    positive_definite = TRUE,
    inverse_quadratic_recipe = TRUE,
    hyperparameters = list(
      estimator = "fixed_identity",
      randomness = "none",
      seed = NULL
    )
  )
}

#' Specify on-demand diagonal residual-variance precision
#'
#' The local metric is the inverse of the residual variance of each feature,
#' floored so that a near-silent feature cannot dominate. Use it for
#' univariate noise normalization when a full local covariance is not wanted
#' or not estimable.
#'
#' @inheritParams identity_metric
#' @param relative_variance_floor Positive floor relative to the mean positive
#'   local residual variance.
#' @param absolute_variance_floor Nonnegative floor in squared response units.
#' @return An `effect_metric_recipe` whose `$hyperparameters` record the
#'   `residual_diagonal_inverse` estimator and both variance floors, with
#'   `$capabilities$native_diagonal` true and `$capabilities$materialized`
#'   false until a schedule binds it.
#' @seealso [shrinkage_precision()] for a full local covariance shrunk to its
#'   diagonal, and [plan_crossnobis()], which consumes the recipe.
#' @family neural metrics
#' @examples
#' # The estimator and both floors are fixed by the recipe, so they are part
#' # of the plan identity rather than a runtime tuning choice.
#' recipe <- diagonal_precision(relative_variance_floor = 1e-6)
#' recipe$kind
#' recipe$hyperparameters[c("estimator", "relative_variance_floor")]
#'
#' # The floor is a guard against dividing by a near-zero residual variance.
#' diagonal_precision(absolute_variance_floor = 0.01)$
#'   hyperparameters$absolute_variance_floor
#'
#' # A nonpositive relative floor would remove that guard, so it is refused.
#' refused <- try(diagonal_precision(relative_variance_floor = 0), silent = TRUE)
#' conditionMessage(attr(refused, "condition"))
#' @export
diagonal_precision <- function(relative_variance_floor = 1e-8,
                               absolute_variance_floor = 0,
                               domain = NULL) {
  floors <- .validate_metric_floor(
    relative_variance_floor, absolute_variance_floor
  )
  .metric_recipe(
    "diagonal_variance_precision", domain,
    native_diagonal = TRUE,
    positive_definite = TRUE,
    inverse_quadratic_recipe = TRUE,
    hyperparameters = c(list(
      estimator = "residual_diagonal_inverse",
      randomness = "none",
      seed = NULL
    ), floors)
  )
}

#' Specify on-demand shrinkage-to-diagonal precision
#'
#' The local covariance estimator is
#' `(1 - shrinkage) * S + shrinkage * diag(diag(S))`, followed by the declared
#' variance and spectral floors. `shrinkage` is fixed by the recipe; it is not
#' tuned on evaluation effects.
#'
#' @inheritParams diagonal_precision
#' @param shrinkage Fixed number in `(0, 1]`.
#' @param relative_spectral_floor Positive minimum eigenvalue relative to the
#'   local covariance scale.
#' @return An `effect_metric_recipe` whose `$hyperparameters` record the
#'   `fixed_shrinkage_to_residual_diagonal` estimator, the fixed
#'   `shrinkage`, and the variance and spectral floors. It is not diagonal,
#'   so it needs a dense local support.
#' @seealso [diagonal_precision()] for the diagonal-only recipe,
#'   [metric_training_policy()] for which partitions may train it, and
#'   [plan_crossnobis()] to compile it.
#' @family neural metrics
#' @examples
#' # The default shrinks the local residual covariance 10% toward its own
#' # diagonal, which keeps small searchlights invertible.
#' recipe <- shrinkage_precision()
#' recipe$hyperparameters$shrinkage
#'
#' # Shrinkage is fixed by the recipe, never tuned on evaluation effects, so
#' # changing it names a different estimand.
#' shrinkage_precision(0.3)$hyperparameters$shrinkage
#'
#' # Unlike a diagonal recipe, this one needs the full local support.
#' metric_capabilities(recipe)[c("native_diagonal", "support_dense")]
#'
#' # Shrinkage must lie in (0, 1]; zero would be an unregularized covariance.
#' refused <- try(shrinkage_precision(0), silent = TRUE)
#' conditionMessage(attr(refused, "condition"))
#' @export
shrinkage_precision <- function(shrinkage = 0.1,
                                relative_variance_floor = 1e-8,
                                absolute_variance_floor = 0,
                                relative_spectral_floor = 1e-10,
                                domain = NULL) {
  if (!is.numeric(shrinkage) || length(shrinkage) != 1L ||
      is.na(shrinkage) || !is.finite(shrinkage) ||
      shrinkage <= 0 || shrinkage > 1) {
    stop("`shrinkage` must be one finite number in (0, 1].",
      call. = FALSE)
  }
  if (!is.numeric(relative_spectral_floor) ||
      length(relative_spectral_floor) != 1L ||
      is.na(relative_spectral_floor) ||
      !is.finite(relative_spectral_floor) ||
      relative_spectral_floor <= 0) {
    stop("`relative_spectral_floor` must be one positive finite number.",
      call. = FALSE)
  }
  floors <- .validate_metric_floor(
    relative_variance_floor, absolute_variance_floor
  )
  .metric_recipe(
    "fixed_diagonal_shrinkage_precision", domain,
    native_diagonal = FALSE,
    positive_definite = TRUE,
    inverse_quadratic_recipe = TRUE,
    hyperparameters = c(list(
      estimator = "fixed_shrinkage_to_residual_diagonal",
      target = "sample_residual_diagonal",
      target_estimated = TRUE,
      shrinkage = as.double(shrinkage),
      relative_spectral_floor = as.double(relative_spectral_floor),
      randomness = "none",
      seed = NULL
    ), floors)
  )
}

.bind_metric_recipe_domain <- function(recipe, domain) {
  recipe <- .validate_metric_recipe(recipe)
  domain <- .domain_reference(domain)
  if (!is.null(recipe$domain)) {
    if (!.same_domain_reference(recipe$domain, domain)) {
      stop("The metric recipe and residual statistics use different neural domains.",
        call. = FALSE)
    }
    return(recipe)
  }
  .metric_recipe(
    recipe$kind, domain,
    native_diagonal = recipe$capabilities$native_diagonal,
    positive_definite = recipe$capabilities$positive_definite,
    inverse_quadratic_recipe =
      recipe$capabilities$inverse_quadratic_recipe,
    hyperparameters = recipe$hyperparameters,
    provenance = c(recipe$provenance, list(
      unbound_recipe_signature = recipe$signature
    ))
  )
}

#' Declare which residual partitions may train a metric
#'
#' `metric_training_policy()` fixes the partition discipline used when a
#' metric recipe is estimated from residuals, so metric training cannot
#' silently borrow the partitions whose products are being evaluated. Pass it
#' to [plan_crossnobis()].
#'
#' @param kind `"exclude_evaluation"` estimates each edge metric without its
#'   two evaluation partitions. `"all_partitions_residual_orthogonality"`
#'   admits evaluation-partition GLM residuals under an explicit orthogonality
#'   justification.
#' @param justification Required for the all-partitions policy. It records why
#'   residual reuse is scientifically admitted; it is not treated as proof.
#' @return An `effect_metric_training_policy` recording `$kind`,
#'   `$includes_evaluation_residuals`, the named `$assumption` it rests on,
#'   the stored `$justification`, and a `$signature` bound into plan identity.
#' @seealso [plan_crossnobis()], which enforces this policy, and
#'   [shrinkage_precision()] for the recipe it governs.
#' @family neural metrics
#'
#' @section Refusal:
#' Requesting `"all_partitions_residual_orthogonality"` without a
#' `justification` signals an `effect_capability_refusal` carrying capability
#' `"evaluation_residual_reuse"` in namespace `"metric_learning"`, with reason
#' `"residual_reuse_justification_absent"`. Inspect it with [catch_refusal()].
#' @examples
#' # The default trains each edge's metric without its two evaluation
#' # partitions, so the metric and the products stay disjoint.
#' policy <- metric_training_policy()
#' policy$kind
#' c(reuses_evaluation = policy$includes_evaluation_residuals,
#'   assumption = policy$assumption)
#'
#' # Reusing evaluation-partition residuals is admitted only with an explicit
#' # written justification, which is recorded, not verified.
#' permissive <- metric_training_policy(
#'   "all_partitions_residual_orthogonality",
#'   justification = "GLM residuals are orthogonal to the fitted effects."
#' )
#' permissive$includes_evaluation_residuals
#'
#' # Omitting that justification is refused, and the refusal is classed, so a
#' # caller can branch on the capability rather than on the prose.
#' refusal <- catch_refusal(
#'   metric_training_policy("all_partitions_residual_orthogonality")
#' )
#' refusal$capability
#' refusal$reasons
#' @export
metric_training_policy <- function(
    kind = c("exclude_evaluation",
             "all_partitions_residual_orthogonality"),
    justification = NULL) {
  kind <- match.arg(kind)
  if (identical(kind, "all_partitions_residual_orthogonality")) {
    if (!is.character(justification) || length(justification) != 1L ||
        is.na(justification) || !nzchar(justification)) {
      .capability_refusal(paste0(
        "Reusing evaluation-partition residuals to train a metric requires ",
        "one explicit written `justification`, which is recorded in the plan ",
        "identity and never treated as proof. Without it the metric would be ",
        "estimated from the same residuals whose products are being ",
        "evaluated, and the resulting distances would be optimistically ",
        "biased by an amount this package cannot quantify."
      ),
        capability = "evaluation_residual_reuse",
        namespace = "metric_learning",
        reasons = "residual_reuse_justification_absent",
        remedies = paste0(
          "Keep the default `metric_training_policy(\"exclude_evaluation\")`, ",
          "or pass `justification = ` stating why residual reuse is admitted ",
          "for this design."
        )
      )
    }
  } else if (!is.null(justification) &&
      (!is.character(justification) || length(justification) != 1L ||
       is.na(justification) || !nzchar(justification))) {
    stop("`justification` must be NULL or one nonempty string.",
      call. = FALSE)
  }
  semantic <- list(
    schema_version = 1L,
    kind = kind,
    includes_evaluation_residuals = identical(
      kind, "all_partitions_residual_orthogonality"
    ),
    assumption = if (identical(
      kind, "all_partitions_residual_orthogonality"
    )) {
      "fitted_effect_and_glm_residual_orthogonality"
    } else {
      "partition_disjoint_metric_training"
    },
    justification = justification
  )
  structure(c(semantic[-1L], list(
    signature = .sha256_signature(semantic)
  )), class = "effect_metric_training_policy")
}

.validate_metric_training_policy <- function(x) {
  expected <- c("kind", "includes_evaluation_residuals", "assumption",
    "justification", "signature")
  if (!inherits(x, "effect_metric_training_policy") || !is.list(x) ||
      !identical(names(x), expected) || !.strong_sha256(x$signature)) {
    stop("Metric-training-policy fields are missing or noncanonical.",
      call. = FALSE)
  }
  rebuilt <- metric_training_policy(x$kind, x$justification)
  if (!identical(x, rebuilt)) {
    stop("Metric-training-policy identity is inconsistent.", call. = FALSE)
  }
  x
}

.preflight_metric_training <- function(recipe, partitions, over, policy) {
  recipe <- .validate_metric_recipe(recipe)
  policy <- .validate_metric_training_policy(policy)
  .validate_pairing(over)
  if (!is.character(partitions) || length(partitions) < 1L ||
      anyNA(partitions) || any(!nzchar(partitions)) ||
      anyDuplicated(partitions) || any(!over$left %in% partitions) ||
      any(!over$right %in% partitions)) {
    stop("Metric-training partitions and evaluation edges are inconsistent.",
      call. = FALSE)
  }
  if (identical(recipe$kind, "identity")) {
    return(invisible(TRUE))
  }
  if (identical(policy$kind, "exclude_evaluation")) {
    available <- vapply(seq_len(nrow(over)), function(edge) {
      length(setdiff(partitions, c(over$left[[edge]], over$right[[edge]])))
    }, integer(1))
    if (any(available < 1L)) {
      stop(sprintf(
        "Evaluation edge %d leaves no residual partition for metric training.",
        which(available < 1L)[[1L]]
      ), call. = FALSE)
    }
  }
  invisible(TRUE)
}

.metric_pairing_identity <- function(over) {
  .validate_pairing(over)
  identity <- list(
    edges = unclass(as.data.frame(over)),
    directed = attr(over, "directed", exact = TRUE),
    self_pairs = attr(over, "self_pairs", exact = TRUE),
    independence = attr(over, "independence", exact = TRUE),
    estimate = attr(over, "estimate", exact = TRUE)
  )
  # Included only when declared so undeclared-axis pairings keep their
  # pre-existing identities; a declared axis is estimand-bearing and must
  # move the digest.
  axis <- attr(over, "generalizes_over", exact = TRUE)
  if (!is.null(axis)) identity$generalizes_over <- axis
  identity
}

.metric_training_record <- function(edge, over, recipe, statistics, policy) {
  left <- over$left[[edge]]
  right <- over$right[[edge]]
  training <- if (identical(recipe$kind, "identity")) {
    character()
  } else if (identical(policy$kind, "exclude_evaluation")) {
    setdiff(statistics$partitions, c(left, right))
  } else {
    statistics$partitions
  }
  if (!identical(recipe$kind, "identity") && length(training) < 1L) {
    stop(sprintf(
      "Evaluation edge %d leaves no residual partition for metric training.",
      edge
    ), call. = FALSE)
  }
  atomic <- statistics$atomic[training]
  atomic_signatures <- vapply(atomic, `[[`, character(1), "signature")
  source_revisions <- vapply(atomic, `[[`, character(1), "source_revision")
  residual_revisions <- vapply(
    atomic, `[[`, character(1), "residual_revision"
  )
  semantic <- list(
    schema_version = 1L,
    edge = as.integer(edge),
    evaluation_left = left,
    evaluation_right = right,
    recipe = recipe$signature,
    statistics = statistics$signature,
    policy = policy$signature,
    training_partitions = training,
    atomic_signatures = atomic_signatures,
    source_revisions = source_revisions,
    residual_revisions = residual_revisions
  )
  structure(c(semantic[c("edge", "evaluation_left", "evaluation_right",
    "training_partitions", "atomic_signatures", "source_revisions",
    "residual_revisions")], list(
    training_signature = .sha256_signature(semantic)
  )), class = "effect_metric_training_record")
}

.metric_schedule_signature <- function(x) {
  .sha256_signature(list(
    schema_version = 1L,
    role = x$role,
    recipe_specification = x$recipe_specification,
    recipe = x$recipe$signature,
    statistics = x$statistics$signature,
    support_index = x$support_index$signature,
    pairing = .metric_pairing_identity(x$pairing),
    training_policy = x$training_policy$signature,
    records = vapply(x$records, `[[`, character(1), "training_signature"),
    capabilities = unclass(x$capabilities)
  ))
}

.validate_frozen_metric_schedule <- function(x, deep = FALSE) {
  expected <- c("role", "recipe_specification", "recipe", "statistics",
    "support_index", "pairing", "training_policy", "records",
    "capabilities", "execution", "signature")
  if (!inherits(x, "effect_frozen_metric_schedule") || !is.list(x) ||
      !identical(names(x), expected) ||
      !identical(x$role, "same_space_metric_schedule") ||
      !.strong_sha256(x$recipe_specification) ||
      !is.list(x$records) || length(x$records) < 1L ||
      !is.list(x$capabilities) || !is.list(x$execution) ||
      !.strong_sha256(x$signature)) {
    stop("Frozen metric-schedule fields are missing or noncanonical.",
      call. = FALSE)
  }
  recipe <- .validate_metric_recipe(x$recipe)
  statistics <- x$statistics
  .validate_residual_pair_statistics(statistics)
  index <- .validate_support_index(x$support_index)
  .validate_pairing(x$pairing)
  policy <- .validate_metric_training_policy(x$training_policy)
  if (is.null(recipe$domain) ||
      !.same_domain_reference(recipe$domain, statistics$domain) ||
      !.same_domain_reference(index$domain, statistics$domain) ||
      !identical(index$signature, statistics$support_index) ||
      length(x$records) != nrow(x$pairing)) {
    stop("Frozen metric-schedule domains, supports, or edges are inconsistent.",
      call. = FALSE)
  }
  expected_capabilities <- list(
    feature_additive = recipe$capabilities$feature_additive,
    support_dense = recipe$capabilities$support_dense,
    learned = !identical(recipe$kind, "identity"),
    provenance_frozen = TRUE,
    materialized = FALSE,
    on_demand = TRUE,
    calibration_requires_metric_uncertainty =
      !identical(recipe$kind, "identity")
  )
  if (!identical(x$capabilities, expected_capabilities)) {
    stop("Frozen metric-schedule capabilities are inconsistent.",
      call. = FALSE)
  }
  if (isTRUE(deep)) {
    expected_records <- lapply(seq_len(nrow(x$pairing)), function(edge) {
      .metric_training_record(
        edge, x$pairing, recipe, statistics, policy
      )
    })
    names(expected_records) <- paste0("edge_", seq_along(expected_records))
    if (!identical(x$records, expected_records) ||
        !identical(x$signature, .metric_schedule_signature(x))) {
      stop("Frozen metric-schedule identity is inconsistent.", call. = FALSE)
    }
  }
  invisible(x)
}

#' Compile a frozen-provenance, on-demand metric schedule
#'
#' Compilation binds a metric recipe to canonical residual sufficient
#' statistics, spatial supports, and evaluation edges. What is frozen is the
#' estimator identity, training assignment, source revisions, and
#' regularization policy. Local covariance and precision matrices are derived
#' when requested and are never retained as a node-by-fold table.
#'
#' @param recipe A recipe from `identity_metric()`, `diagonal_precision()`, or
#'   `shrinkage_precision()`.
#' @param statistics Canonical `residual_pair_statistics()`.
#' @param at The same support-index-backed frame used to accumulate
#'   `statistics`.
#' @param over An evaluation `pairing()`.
#' @param training A `metric_training_policy()`.
#' @return An `effect_frozen_metric_schedule`.
#' @keywords internal
compile_metric_schedule <- function(
    recipe, statistics, at, over,
    training = metric_training_policy("exclude_evaluation")) {
  recipe_specification <- .validate_metric_recipe(recipe)$signature
  .validate_residual_pair_statistics(statistics, deep = TRUE)
  index <- .residual_statistics_support_index(at, statistics$domain)
  if (!identical(index$signature, statistics$support_index)) {
    stop("Metric compilation requires the exact support index used by the residual statistics.",
      call. = FALSE)
  }
  .validate_pairing(over)
  if (any(!over$left %in% statistics$partitions) ||
      any(!over$right %in% statistics$partitions)) {
    stop("Every metric evaluation endpoint must identify a residual-statistics partition.",
      call. = FALSE)
  }
  policy <- .validate_metric_training_policy(training)
  .preflight_metric_training(
    recipe, statistics$partitions, over, policy
  )
  recipe <- .bind_metric_recipe_domain(recipe, statistics$domain)
  supported <- c(
    "identity", "diagonal_variance_precision",
    "fixed_diagonal_shrinkage_precision"
  )
  if (!recipe$kind %in% supported) {
    stop("This metric recipe has no admitted on-demand learner.",
      call. = FALSE)
  }
  coordinates <- .residual_pair_coordinates(index)
  if (!identical(coordinates$i, statistics$pair_i) ||
      !identical(coordinates$j, statistics$pair_j)) {
    stop("Residual statistics do not match the support pair graph order.",
      call. = FALSE)
  }
  records <- lapply(seq_len(nrow(over)), function(edge) {
    .metric_training_record(edge, over, recipe, statistics, policy)
  })
  names(records) <- paste0("edge_", seq_along(records))
  capabilities <- list(
    feature_additive = recipe$capabilities$feature_additive,
    support_dense = recipe$capabilities$support_dense,
    learned = !identical(recipe$kind, "identity"),
    provenance_frozen = TRUE,
    materialized = FALSE,
    on_demand = TRUE,
    calibration_requires_metric_uncertainty =
      !identical(recipe$kind, "identity")
  )
  value <- structure(list(
    role = "same_space_metric_schedule",
    recipe_specification = recipe_specification,
    recipe = recipe,
    statistics = statistics,
    support_index = index,
    pairing = over,
    training_policy = policy,
    records = records,
    capabilities = capabilities,
    execution = list(
      residual_accumulation = statistics$execution,
      local_metric_storage = "none_derived_on_demand",
      retained_factor_table = FALSE
    ),
    signature = NA_character_
  ), class = "effect_frozen_metric_schedule")
  value$signature <- .metric_schedule_signature(value)
  .validate_frozen_metric_schedule(value, deep = TRUE)
}

.metric_schedule_edge <- function(schedule, edge) {
  if (is.character(edge) && length(edge) == 1L &&
      !is.na(edge) && edge %in% names(schedule$records)) {
    return(match(edge, names(schedule$records)))
  }
  if (!is.numeric(edge) || length(edge) != 1L || is.na(edge) ||
      !is.finite(edge) || edge %% 1 != 0 || edge < 1L ||
      edge > length(schedule$records)) {
    stop("`edge` must identify one compiled metric evaluation edge.",
      call. = FALSE)
  }
  as.integer(edge)
}

.metric_schedule_node_trusted <- function(index, node) {
  if (is.numeric(node) && length(node) == 1L && !is.na(node) &&
      is.finite(node) && node %% 1 == 0 && node >= 1L &&
      node <= length(index$node_ids)) {
    return(as.integer(node))
  }
  position <- match(node, index$node_ids)
  if (length(position) != 1L || is.na(position)) {
    stop("`node` must identify one compiled metric support.", call. = FALSE)
  }
  as.integer(position)
}

.local_residual_covariance <- function(index, support, covariance) {
  index <- .validate_support_index(index)
  .local_residual_covariance_trusted(index, support, covariance)
}

.local_residual_covariance_trusted <- function(index, support, covariance) {
  pattern <- index$pair_pattern
  if (!identical(pattern@uplo, "U") ||
      length(covariance) != length(pattern@i)) {
    stop("Local covariance extraction requires the canonical upper pair graph.",
      call. = FALSE)
  }
  dimension <- length(support)
  value <- matrix(0, dimension, dimension)
  for (column in seq_len(dimension)) {
    global_column <- support[[column]]
    first <- pattern@p[[global_column]] + 1L
    last <- pattern@p[[global_column + 1L]]
    slots <- seq.int(first, last)
    rows <- pattern@i[slots] + 1L
    desired <- support[seq_len(column)]
    matched <- match(desired, rows)
    if (anyNA(matched)) {
      stop("A support pair is absent from its declared union pair graph.",
        call. = FALSE)
    }
    pair_values <- covariance[slots[matched]]
    value[seq_len(column), column] <- pair_values
    value[column, seq_len(column)] <- pair_values
  }
  value
}

.local_residual_scope_covariance <- function(statistics, index, support,
                                             partitions) {
  if (!is.character(partitions) || length(partitions) < 1L ||
      anyNA(partitions) || any(!nzchar(partitions)) ||
      anyDuplicated(partitions) ||
      any(!partitions %in% statistics$partitions)) {
    stop("A local residual scope must select unique atomic partitions.",
      call. = FALSE)
  }
  # Match .canonical_reduce() exactly, but gather only the requested support.
  # This avoids an edge-specific vector over the complete union pair graph
  # while retaining the same elementwise floating-point reduction order.
  partitions <- partitions[order(partitions, method = "radix")]
  combined <- matrix(0, length(support), length(support))
  residual_df <- 0
  for (partition in partitions) {
    atomic <- statistics$atomic[[partition]]
    local <- .local_residual_covariance_trusted(
      index, support, atomic$cross_products
    )
    combined[] <- combined + local
    rm(local)
    residual_df <- residual_df + atomic$residual_df
  }
  combined / residual_df
}

.metric_covariance_diagnostics <- function(covariance, tolerance = 1e-10) {
  spectrum <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(abs(spectrum), max(abs(diag(covariance))),
    .Machine$double.xmin)
  positive <- spectrum > tolerance * scale
  list(
    eigenvalues = spectrum,
    scale = scale,
    estimated_rank = as.integer(sum(positive)),
    condition = if (all(positive)) max(spectrum) / min(spectrum) else Inf
  )
}

.metric_variance_floor <- function(covariance, recipe) {
  variance <- diag(covariance)
  positive <- variance[variance > 0]
  if (!length(positive)) {
    stop("Residual covariance has no positive local variance.", call. = FALSE)
  }
  hyper <- recipe$hyperparameters
  floor <- max(
    hyper$absolute_variance_floor,
    hyper$relative_variance_floor * mean(positive)
  )
  list(value = pmax(variance, floor), floor = floor,
    floored = as.integer(sum(variance < floor)))
}

.scheduled_metric_handle <- function(recipe, domain, support_index,
                                     schedule_signature,
                                     calibration_requires_uncertainty,
                                     record, statistics, node) {
  node <- .metric_schedule_node_trusted(support_index, node)
  support_positions <- .support_index_support_trusted(
    support_index, node
  )
  support <- domain$feature_ids[support_positions]
  dimension <- length(support)
  raw_covariance <- if (identical(recipe$kind, "identity")) {
    diag(dimension)
  } else {
    .local_residual_scope_covariance(
      statistics, support_index, support_positions,
      record$training_partitions
    )
  }
  raw_diagnostics <- .metric_covariance_diagnostics(raw_covariance)
  variance <- NULL
  ridge <- 0
  if (identical(recipe$kind, "identity")) {
    covariance <- raw_covariance
    precision_diagonal <- rep(1, dimension)
    factor <- NULL
    estimator <- "fixed_identity"
  } else if (identical(recipe$kind, "diagonal_variance_precision")) {
    variance <- .metric_variance_floor(raw_covariance, recipe)
    covariance <- diag(variance$value, dimension)
    precision_diagonal <- 1 / variance$value
    factor <- NULL
    estimator <- "residual_diagonal_inverse"
  } else {
    variance <- .metric_variance_floor(raw_covariance, recipe)
    alpha <- recipe$hyperparameters$shrinkage
    covariance <- (1 - alpha) * raw_covariance +
      alpha * diag(variance$value, dimension)
    before_ridge <- .metric_covariance_diagnostics(covariance)
    target_minimum <- recipe$hyperparameters$relative_spectral_floor *
      before_ridge$scale
    ridge <- max(0, target_minimum - min(before_ridge$eigenvalues))
    if (ridge > 0) covariance <- covariance + diag(ridge, dimension)
    precision_diagonal <- NULL
    factor <- chol(covariance)
    estimator <- "fixed_shrinkage_to_residual_diagonal"
  }
  final_diagnostics <- .metric_covariance_diagnostics(covariance)
  if (!is.finite(final_diagnostics$condition)) {
    stop("Regularized local covariance is not positive definite.",
      call. = FALSE)
  }
  apply_K <- if (!is.null(precision_diagonal)) {
    function(value) {
      if (is.atomic(value) && is.null(dim(value))) {
        if (!is.numeric(value) || length(value) != dimension ||
            any(!is.finite(value))) {
          stop("Metric application requires one finite local vector.",
            call. = FALSE)
        }
        return(precision_diagonal * value)
      }
      if (!is.matrix(value) || !is.numeric(value) ||
          nrow(value) != dimension || any(!is.finite(value))) {
        stop("Metric application requires finite local vectors in columns.",
          call. = FALSE)
      }
      precision_diagonal * value
    }
  } else {
    function(value) {
      if (is.atomic(value) && is.null(dim(value))) {
        if (!is.numeric(value) || length(value) != dimension ||
            any(!is.finite(value))) {
          stop("Metric application requires one finite local vector.",
            call. = FALSE)
        }
      } else if (!is.matrix(value) || !is.numeric(value) ||
          nrow(value) != dimension || any(!is.finite(value))) {
        stop("Metric application requires finite local vectors in columns.",
          call. = FALSE)
      }
      backsolve(factor, forwardsolve(t(factor), value))
    }
  }
  quadform_Kinv <- function(value) {
    if (!is.numeric(value) || is.matrix(value) ||
        length(value) != dimension || any(!is.finite(value))) {
      stop("Inverse quadratic forms require one finite local vector.",
        call. = FALSE)
    }
    drop(crossprod(value, covariance %*% value))
  }
  provenance <- list(
    frozen = !identical(recipe$kind, "identity"),
    training_signature = record$training_signature,
    schedule = schedule_signature,
    recipe = recipe$signature,
    edge = record$edge,
    node = as.character(support_index$node_ids[[node]]),
    training_partitions = record$training_partitions,
    source_revisions = record$source_revisions,
    residual_revisions = record$residual_revisions,
    estimator = estimator
  )
  diagnostics <- list(
    estimator = estimator,
    support_size = as.integer(dimension),
    estimated_rank = raw_diagnostics$estimated_rank,
    condition_before = raw_diagnostics$condition,
    condition_after = final_diagnostics$condition,
    variance_floor = if (is.null(variance)) 0 else variance$floor,
    variance_entries_floored = if (is.null(variance)) 0L else variance$floored,
    spectral_ridge = ridge,
    shrinkage = if (identical(recipe$kind,
      "fixed_diagonal_shrinkage_precision")) {
      recipe$hyperparameters$shrinkage
    } else {
      0
    },
    structural_rank_guarantee = FALSE,
    calibration_requires_metric_uncertainty =
      calibration_requires_uncertainty
  )
  materialize <- function() {
    precision <- if (!is.null(precision_diagonal)) {
      diag(precision_diagonal, dimension)
    } else {
      chol2inv(factor)
    }
    neural_metric(
      precision, domain, support,
      inverse = covariance,
      estimation = if (identical(recipe$kind, "identity")) {
        "fixed"
      } else {
        "learned_frozen"
      },
      provenance = provenance
    )
  }
  form <- function(left, right) {
    if (!is.matrix(left) || !is.numeric(left) || ncol(left) != dimension ||
        any(!is.finite(left)) || !is.matrix(right) || !is.numeric(right) ||
        ncol(right) != dimension || any(!is.finite(right))) {
      stop("Metric forms require finite matrices sharing the local feature axis.",
        call. = FALSE)
    }
    left %*% apply_K(t(right))
  }
  semantic <- list(
    schema_version = 1L,
    schedule = schedule_signature,
    training = record$training_signature,
    node = as.character(support_index$node_ids[[node]]),
    support = support,
    covariance = covariance,
    diagnostics = diagnostics
  )
  structure(list(
    support = support,
    support_positions = support_positions,
    covariance = covariance,
    apply_K = apply_K,
    quadform_Kinv = quadform_Kinv,
    form = form,
    materialize = materialize,
    diagnostics = diagnostics,
    provenance = provenance,
    signature = .sha256_signature(semantic)
  ), class = "effect_scheduled_metric_handle")
}

.new_metric_schedule_provider <- function(
    recipe, domain, support_index, schedule_signature,
    calibration_requires_uncertainty, record, statistics,
    selected_execution, edge) {
  state <- new.env(parent = emptyenv())
  state$nodes_derived <- 0L
  at <- function(node) {
    state$nodes_derived <- state$nodes_derived + 1L
    .scheduled_metric_handle(
      recipe, domain, support_index, schedule_signature,
      calibration_requires_uncertainty, record, statistics, node
    )
  }
  structure(list(
    edge = edge,
    record = record,
    at = at,
    receipt = function() list(
      schedule = schedule_signature,
      edge = edge,
      training_signature = record$training_signature,
      training_partitions = record$training_partitions,
      atomic_statistics = record$atomic_signatures,
      atomic_accumulation_reads = vapply(
        selected_execution, `[[`, integer(1), "residual_reads"
      ),
      residual_reads_during_derivation = 0L,
      nodes_derived = state$nodes_derived,
      retained_factor_table = FALSE
    )
  ), class = "effect_metric_schedule_provider")
}

.metric_schedule_provider <- function(schedule, edge = 1L) {
  .validate_frozen_metric_schedule(schedule, deep = FALSE)
  edge <- .metric_schedule_edge(schedule, edge)
  record <- schedule$records[[edge]]
  selected_execution <- schedule$statistics$execution$atomic[
    record$training_partitions
  ]
  .new_metric_schedule_provider(
    schedule$recipe,
    schedule$statistics$domain,
    schedule$support_index,
    schedule$signature,
    schedule$capabilities$calibration_requires_metric_uncertainty,
    record,
    schedule$statistics,
    selected_execution,
    edge
  )
}

#' Materialize one metric from an on-demand schedule
#'
#' This is an inspection operation. Execution kernels should request one
#' edge-scoped provider and reuse its combined sufficient statistics while
#' deriving local solve handles support by support.
#'
#' @param schedule A compiled metric schedule.
#' @param node One support position or node identifier.
#' @param edge One evaluation-edge position or name.
#' @return A support-local `neural_metric()`.
#' @keywords internal
materialize_metric <- function(schedule, node, edge = 1L) {
  provider <- .metric_schedule_provider(schedule, edge)
  provider$at(node)$materialize()
}
