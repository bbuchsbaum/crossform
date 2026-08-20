# Crossvalidated Mahalanobis contrast views -------------------------------

#' Construct a fixed neural noise-precision metric
#'
#' `noise_precision()` is the semantically specific constructor for a known or
#' externally estimated precision that is fixed before effect evaluation. It
#' uses the same `neural_metric()` representation while recording that the
#' operator has inverse-noise-covariance meaning.
#'
#' @param value A finite symmetric positive-semidefinite precision matrix.
#' @param domain Exact neural feature domain.
#' @param support Ordered support identities for a local matrix.
#' @param covariance Optional retained inverse covariance.
#' @param tolerance Positive numerical validation tolerance.
#' @param provenance Additional compact provenance. Reserved semantic fields
#'   cannot be replaced.
#' @return An `effect_neural_metric` whose `$provenance$metric_role` is
#'   `"noise_precision"`, carrying the canonicalized `$value`, its `$domain`
#'   and `$support`, `$capabilities`, and a `$signature`. The recorded role is
#'   what lets [crossnobis()] name the result a Mahalanobis distance.
#' @seealso [crossnobis()] and [plan_geometry()] (its `metric` argument);
#'   [shrinkage_precision()] when the precision must instead be learned from
#'   residuals.
#' @family neural metrics
#' @examples
#' # A precision estimated outside crossform and fixed before evaluation.
#' domain <- abstract_domain(3, id = "noise-precision-example")
#' covariance <- matrix(c(1, 0.4, 0.1, 0.4, 1, 0.3, 0.1, 0.3, 1), 3, 3)
#' metric <- noise_precision(
#'   solve(covariance), domain, covariance = covariance,
#'   provenance = list(source = "resting-run residual covariance")
#' )
#' metric$provenance$metric_role
#' metric$capabilities$positive_definite
#'
#' # The role is what a crossnobis plan checks for, so a plain
#' # `neural_metric()` with the same numbers is deliberately not equivalent.
#' identical(
#'   neural_metric(solve(covariance), domain)$provenance$metric_role,
#'   metric$provenance$metric_role
#' )
#'
#' # The semantic provenance fields are reserved and cannot be overwritten.
#' refused <- try(
#'   noise_precision(diag(3), domain, provenance = list(metric_role = "other")),
#'   silent = TRUE
#' )
#' conditionMessage(attr(refused, "condition"))
#' @export
noise_precision <- function(value, domain, support = NULL,
                            covariance = NULL, tolerance = 1e-10,
                            provenance = list()) {
  .validate_effect_provenance(provenance, "noise-precision provenance")
  reserved <- c("metric_role", "estimator_status")
  if (!is.null(names(provenance)) && any(names(provenance) %in% reserved)) {
    .input_error("Noise-precision semantic provenance fields are reserved.")
  }
  neural_metric(
    value, domain, support = support, inverse = covariance,
    estimation = "fixed", tolerance = tolerance,
    provenance = c(list(
      metric_role = "noise_precision",
      estimator_status = "fixed_before_effect_evaluation"
    ), provenance)
  )
}

.crossnobis_plan_metric <- function(x) {
  .validate_geometry_plan(x)
  schedule <- x$metric_schedule
  if (!identical(schedule$kind, "fixed_metric_before_frame")) {
    .capability_refusal(paste0(
      "Crossnobis is a Mahalanobis reading, so it requires an explicit ",
      "noise-precision metric; this plan carries the implicit identity ",
      "metric, which is Euclidean and declares nothing about the noise."
    ),
      capability = "declared_noise_metric",
      namespace = "geometry_views",
      reasons = "implicit_identity_metric_is_not_a_noise_model",
      remedies = paste0(
        "Compile the plan with `plan_geometry(..., metric = ",
        "noise_precision(...))`, use `plan_crossnobis()` to learn the ",
        "precision from residuals, or read the same estimand without a ",
        "noise claim through `contrast_energy()`."
      )
    )
  }
  metric <- .validate_neural_metric(schedule$metric, deep = FALSE)
  if (!identical(metric$provenance$metric_role, "noise_precision")) {
    .capability_refusal(paste0(
      "Crossnobis interpretation requires a metric declared as a noise ",
      "precision; this plan's fixed metric was not constructed by ",
      "`noise_precision()` and carries no such role, so calling its output a ",
      "Mahalanobis distance would overstate what the metric means."
    ),
      capability = "declared_noise_metric",
      namespace = "geometry_views",
      reasons = "metric_role_is_not_noise_precision",
      remedies = paste0(
        "Build the metric with `noise_precision()`, or read the plan with ",
        "`contrast_energy()`, which makes no noise-model claim."
      )
    )
  }
  metric
}

# The learned counterpart of `.crossnobis_plan_metric()`. A recipe earns the
# Mahalanobis reading by being a residual-derived precision; the identity
# recipe is a legitimate schedule but declares nothing about the noise, so it
# gets the same refusal a plan with the implicit identity metric gets.
.crossnobis_learned_metric <- function(x) {
  schedule <- x$metric_schedule$schedule
  if (identical(schedule$recipe$kind, "identity")) {
    .capability_refusal(paste0(
      "Crossnobis is a Mahalanobis reading, so it requires a residual-",
      "derived precision; this plan's recipe is `identity_metric()`, which ",
      "is Euclidean and declares nothing about the noise."
    ),
      capability = "declared_noise_metric",
      namespace = "geometry_views",
      reasons = "identity_recipe_is_not_a_noise_model",
      remedies = paste0(
        "Compile the plan with `shrinkage_precision()` or ",
        "`diagonal_precision()`, or read the same estimand without a noise ",
        "claim through `contrast_energy()`."
      )
    )
  }
  schedule
}

.require_crossnobis_pairing <- function(over) {
  .validate_pairing(over)
  if (!identical(attr(over, "independence", exact = TRUE), "independent") ||
      !identical(attr(over, "estimate", exact = TRUE),
        "cross_generalized") || any(over$left == over$right)) {
    .input_error(sprintf(paste0(
      "Crossnobis requires cross-partition edges declared independent and ",
      "containing no self-products; this pairing declares independence `%s` ",
      "with estimate `%s`%s. The crossvalidated distance is unbiased only ",
      "because the two endpoints carry independent noise."
    ), attr(over, "independence", exact = TRUE),
      attr(over, "estimate", exact = TRUE),
      if (any(over$left == over$right)) {
        sprintf(" and contains %s",
          .msg_count(sum(over$left == over$right), "self-product"))
      } else {
        ""
      }))
  }
  invisible(TRUE)
}

.crossnobis_scientific_plan_id <- function(relation_fit_signature, task,
                                            frame, pairing,
                                            metric_schedule) {
  .sha256_signature(list(
    schema_version = 1L,
    relation_fit = relation_fit_signature,
    evidence_task = task$task_id,
    frame = .additive_frame_signature(frame),
    pairing = .metric_pairing_identity(pairing),
    metric_schedule = metric_schedule$signature
  ), "crossnobis-sha256:")
}

.crossnobis_plan_signature <- function(fields) {
  .sha256_signature(list(
    schema_version = 1L,
    scientific_plan_id = fields$scientific_plan_id,
    compute = unclass(fields$compute),
    lowering = fields$lowering,
    kernel_version = fields$kernel_version,
    memory = unclass(fields$memory)
  ))
}

.learned_crossnobis_memory_plan <- function(fit, frame, support_index, over,
                                            compute) {
  relation <- fit$relation
  support_index <- .support_index_materialize_pair_pattern(
    .validate_support_index(support_index)
  )
  .validate_pairing(over)
  support_sizes <- support_index$cost$support_size
  k <- as.double(support_sizes[["max"]])
  q <- length(relation$effects)
  endpoints <- unique(c(over$left, over$right))
  maximum_observations <- max(vapply(
    relation$sources, function(source) source$dim[[1L]], integer(1)
  ))
  pair_count <- length(support_index$pair_pattern@i)
  partitions <- length(relation$partitions)
  support_entries <- length(support_index$members)
  measurements <- length(support_index$node_ids)
  # Payload formulas only: preflight never allocates the buffers it sizes.
  # Statistics retain pair coordinates and one cross-product vector per
  # partition. The support index retains CSR membership plus the canonical
  # symmetric pair-pattern CSC slots; these are charged separately instead of
  # disappearing into an unmeasured object header.
  resident_statistics <- pair_count * (8 * partitions + 8) +
    support_entries * 4 + (measurements + 1) * 8 +
    length(support_index$pair_pattern@i) * 4 +
    length(support_index$pair_pattern@p) * 4
  canonical_weights <- .canonical_additive_weights(frame$weights)
  frame_payload <- if (inherits(canonical_weights, "sparseMatrix")) {
    length(canonical_weights@x) * 8 + length(canonical_weights@i) * 4 +
      length(canonical_weights@p) * 4
  } else {
    prod(as.double(dim(canonical_weights))) * 8
  }
  source_shared <- sum(vapply(relation$sources, function(source) {
    if (identical(source$kind, "matrix")) {
      prod(as.double(source$dim)) * 8
    } else {
      0
    }
  }, numeric(1)))
  distinct_handles <- unique(vapply(relation$sources, function(source) {
    descriptor <- source$descriptor
    if (is.null(descriptor) || identical(descriptor$access, "coordinator")) {
      return(NA_character_)
    }
    .source_descriptor_key(descriptor)
  }, character(1)))
  distinct_handles <- sum(!is.na(distinct_handles))
  plan <- memory_plan(
    frame_bytes = frame_payload,
    resident_source_bytes = source_shared + resident_statistics,
    source_handle_bytes = distinct_handles * 4096,
    source_block_bytes = 8 * maximum_observations * k,
    relation_block_bytes = 8 * length(endpoints) * q * k,
    atom_block_bytes = 0,
    local_state_bytes = 0,
    output_bytes = 8 * measurements,
    # At most two local covariance matrices are live while atomics are reduced;
    # regularization then reuses that storage before factorization. Weighted
    # endpoint patterns and one relation block complete the active task.
    contraction_bytes = 8 * (2 * k * k + 2 * length(endpoints) * k + q * k),
    replacement_copy_bytes = 8 * max(k * k, measurements),
    workers = 1L,
    n_active = 1L,
    budget_bytes = compute$workspace_bytes
  )
  if (identical(plan$fits_budget, FALSE)) {
    .input_error(sprintf(
      paste0("Learned crossnobis requires %.0f bytes, exceeding the ",
        "%.0f-byte workspace budget."),
      plan$planned_workspace_bytes, compute$workspace_bytes
    ))
  }
  plan
}

#' Compile an on-demand learned-metric crossnobis plan
#'
#' The plan freezes the residual-statistics identity, training policy,
#' regularization recipe, evaluation pairing, and spatial support graph. It
#' does not retain one covariance or precision matrix per node and fold; each
#' local solve handle is derived on demand from canonical pair statistics.
#' The learned metric is itself random. Version 0.1 therefore returns signed
#' crossnobis point estimates but does not apply the fixed-metric analytic
#' covariance law in [rdm_sampling_covariance()] or claim calibrated intervals
#' for this plan; doing so requires propagation of metric-estimation
#' uncertainty (for example, an admitted LD-t specialization).
#'
#' `plan_crossnobis()` is sugar: it names the intent, and compiles exactly
#' `plan_geometry(x, at, over, metric = <recipe>, training = )`. The plan it
#' returns is an ordinary `effect_geometry_plan` whose metric schedule has
#' kind `learned_local_before_frame`, lowered by the geometry compiler like
#' any other. The pairing contract crossnobis requires -- independent,
#' cross-partition, self-product-free -- is enforced by [crossnobis()] when
#' the plan is read, which is where the fixed-metric route enforces it.
#'
#' @param x An `effect_relation_fit` with residual-block capability.
#' @param at A support-index-backed compiled spatial frame.
#' @param over Independent cross-partition evaluation edges.
#' @param metric An on-demand metric recipe such as `shrinkage_precision()`.
#' @param training A `metric_training_policy()`.
#' @param compute A sequential `compute_policy()`.
#' @param residual_workspace_bytes Positive budget used while accumulating
#'   canonical residual pair sufficient statistics. It changes cache capacity,
#'   never the canonical numerical tile shape. Defaults to
#'   `compute$workspace_bytes`, or 512 MiB when the policy declares none.
#' @return An `effect_geometry_plan` reusable across fixed contrasts. It
#'   carries the learned `$metric_schedule` wrapping the frozen schedule, the
#'   `$frame`, `$pairing`, and `$execution_hints`, and a
#'   `$scientific_plan_id` that identifies the estimand independently of
#'   execution choices.
#' @seealso [plan_geometry()], which this wraps; [crossnobis()] to read a
#'   contrast from this plan,
#'   [shrinkage_precision()] and [metric_training_policy()] for the metric
#'   declarations it freezes, and [residual_pair_statistics()] for the
#'   sufficient statistics it compiles.
#' @family geometry plans and views
#' @examples
#' set.seed(7)
#' domain <- abstract_domain(
#'   3, coordinates = cbind(x = 0:2, y = 0),
#'   id = "learned-crossnobis-example"
#' )
#' design <- cbind(intercept = 1, condition = rep(c(-0.5, 0.5), 4))
#' coefficients <- rbind(
#'   intercept = c(0, 0, 0),
#'   condition = c(0.6, -0.4, 0.2)
#' )
#' responses <- setNames(lapply(seq_len(3), function(run) {
#'   design %*% coefficients + matrix(rnorm(8 * 3, sd = 0.25), 8, 3)
#' }), paste0("run", seq_len(3)))
#' fit <- lm_relation_fit(
#'   responses, design, rbind(condition = c(0, 1)), domain = domain
#' )
#' plan <- plan_crossnobis(
#'   fit,
#'   compile_frame(searchlights(1.01), domain),
#'   pairing("run1", "run2", independence = "independent"),
#'   metric = shrinkage_precision(0.2)
#' )
#' result <- crossnobis(plan, c(condition = 1))
#' result
#' as.data.frame(result)
#' @export
plan_crossnobis <- function(
    x, at, over, metric = shrinkage_precision(),
    training = metric_training_policy("exclude_evaluation"),
    compute = compute_policy(),
    residual_workspace_bytes = NULL) {
  if (missing(at)) {
    .input_error(paste0(
      "`at` is required: pass a compiled frame from `compile_frame()`, for ",
      "example `compile_frame(searchlights(8), domain)`."
    ))
  }
  if (missing(over)) {
    .input_error(paste0(
      "`over` is required: pass a pairing from `cross_partitions()` or ",
      "`pairing()` declaring which partition products may be formed."
    ))
  }
  if (inherits(metric, "effect_metric_recipe") &&
      identical(.validate_metric_recipe(metric)$kind, "identity")) {
    .input_error(paste0(
      "A learned crossnobis plan requires a residual-derived precision ",
      "recipe; use `noise_precision(diag(...))` for a fixed identity metric."
    ))
  }
  plan_geometry(
    x, at, over, compute = compute, metric = metric, training = training,
    residual_workspace_bytes = residual_workspace_bytes
  )
}

# The retired learned-crossnobis plan constructor. Nothing public reaches it:
# `plan_crossnobis()` above compiles an ordinary geometry plan, and the
# geometry compiler lowers it. It is kept for one ticket so that the
# old-versus-new numerical equality can be asserted against the route it
# replaces, and is deleted in B3 together with `effect_crossnobis_plan`,
# `.execute_learned_crossnobis()`, and their identity and receipt helpers.
.plan_learned_crossnobis <- function(
    x, at, over, metric = shrinkage_precision(),
    training = metric_training_policy("exclude_evaluation"),
    compute = compute_policy(),
    residual_workspace_bytes = if (is.null(compute$workspace_bytes)) {
      512 * 1024^2
    } else {
      compute$workspace_bytes
    }) {
  # Diagnose a missing residual channel before shape validation or any
  # metric-training preflight: a bare relation or a channel-free fit must get
  # the capability refusal, not a field-shape error or a training-partition
  # shortage that misstates the real problem.
  .require_relation_fit_capability(x, "learned_metric_input")
  .validate_relation_fit(x, deep = FALSE)
  compute <- .validate_compute_policy(compute)
  .require_crossnobis_pairing(over)
  .validate_geometry_plan_inputs(x$relation, at, over)
  metric <- .validate_metric_recipe(metric)
  if (identical(metric$kind, "identity")) {
    .input_error(paste0(
      "A learned crossnobis plan requires a residual-derived precision ",
      "recipe; use `noise_precision(diag(...))` for a fixed identity metric."
    ))
  }
  training <- .validate_metric_training_policy(training)
  .preflight_metric_training(
    metric, x$relation$partitions, over, training
  )
  capabilities <- .execution_preflight(compute, function() {
    .relation_source_capabilities(x$relation)
  })$source_capabilities
  relation <- x$relation
  relation$capabilities <- capabilities
  support_index <- .residual_statistics_support_index(
    at, x$relation$domain
  )
  memory <- .learned_crossnobis_memory_plan(
    x, at, support_index, over, compute
  )
  statistics <- residual_pair_statistics(
    x, at, workspace_bytes = residual_workspace_bytes
  )
  schedule <- compile_metric_schedule(metric, statistics, at, over, training)
  task <- .compile_effect_evidence_task(relation, over)
  lowering <- "support_streamed_scheduled_metric_query_contraction"
  kernel_version <- "support-streamed-scheduled-metric-v1"
  scientific_plan_id <- .crossnobis_scientific_plan_id(
    x$signature, task, at, over, schedule
  )
  fields <- list(
    relation_fit_signature = x$signature,
    task = task,
    frame = at,
    pairing = over,
    metric_schedule = schedule,
    compute = compute,
    memory = memory,
    lowering = lowering,
    kernel_version = kernel_version,
    scientific_plan_id = scientific_plan_id
  )
  value <- structure(c(fields, list(
    signature = .crossnobis_plan_signature(fields)
  )), class = "effect_crossnobis_plan")
  .validate_crossnobis_plan(value, deep = FALSE)
  value
}

.validate_crossnobis_plan <- function(x, deep = TRUE) {
  expected <- c("relation_fit_signature", "task", "frame", "pairing",
    "metric_schedule", "compute", "memory", "lowering", "kernel_version",
    "scientific_plan_id", "signature")
  if (!.sealed_fields(x, "effect_crossnobis_plan", expected) ||
      !.strong_sha256(x$relation_fit_signature) ||
      !identical(x$lowering, "support_streamed_scheduled_metric_query_contraction") ||
      !identical(x$kernel_version, "support-streamed-scheduled-metric-v1") ||
      !grepl("^crossnobis-sha256:[[:xdigit:]]{64}$", x$scientific_plan_id) ||
      !.strong_sha256(x$signature)) {
    .input_error("Crossnobis-plan fields are missing or noncanonical.")
  }
  .validate_evidence_task(x$task)
  .validate_frame_for_compile(x$frame)
  .require_crossnobis_pairing(x$pairing)
  .validate_frozen_metric_schedule(x$metric_schedule, deep = deep)
  .validate_compute_policy(x$compute)
  .validate_memory_plan_for_receipt(x$memory)
  if (!isTRUE(x$task$same_relation) ||
      !identical(x$metric_schedule$pairing, x$pairing) ||
      !identical(x$metric_schedule$support_index$signature,
        x$frame$support_index$signature) ||
      !.same_domain_reference(x$frame$domain,
        x$task$left_relation$domain)) {
    .contract_error(
      "Crossnobis-plan relation, support, pairing, or metric is inconsistent."
    )
  }
  if (isTRUE(deep)) {
    fields <- x[names(x) != "signature"]
    expected_scientific <- .crossnobis_scientific_plan_id(
      x$relation_fit_signature, x$task, x$frame, x$pairing,
      x$metric_schedule
    )
    if (!identical(x$scientific_plan_id, expected_scientific) ||
        !identical(x$signature, .crossnobis_plan_signature(fields))) {
      .contract_error("Crossnobis-plan identity is inconsistent.")
    }
  }
  invisible(x)
}

#' @export
print.effect_crossnobis_plan <- function(x, ...) {
  .validate_crossnobis_plan(x, deep = FALSE)
  payload <- structure(x$memory$planned_workspace_bytes, class = "object_size")
  cat("<effect_crossnobis_plan>\n", sep = "")
  cat("  measurements: ", nrow(x$frame$weights), "\n", sep = "")
  cat("  metric:       ", x$metric_schedule$recipe$kind, "\n", sep = "")
  cat("  training:     ", x$metric_schedule$training_policy$kind, "\n", sep = "")
  cat("  lowering:     ", x$lowering, "\n", sep = "")
  cat("  workspace:    ", format(payload, units = "auto"), "\n", sep = "")
  invisible(x)
}
