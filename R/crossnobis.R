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
#' @return A fixed `effect_neural_metric` with noise-precision identity.
#' @export
noise_precision <- function(value, domain, support = NULL,
                            covariance = NULL, tolerance = 1e-10,
                            provenance = list()) {
  .validate_effect_provenance(provenance, "noise-precision provenance")
  reserved <- c("metric_role", "estimator_status")
  if (!is.null(names(provenance)) && any(names(provenance) %in% reserved)) {
    stop("Noise-precision semantic provenance fields are reserved.",
      call. = FALSE)
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
    stop(paste0(
      "Crossnobis requires an explicit noise-precision metric schedule; ",
      "an implicit Euclidean geometry has no declared noise metric."
    ), call. = FALSE)
  }
  metric <- .validate_neural_metric(schedule$metric, deep = FALSE)
  if (!identical(metric$provenance$metric_role, "noise_precision")) {
    stop(paste0(
      "Crossnobis interpretation requires a metric constructed as ",
      "`noise_precision()` or an admitted learned precision schedule."
    ), call. = FALSE)
  }
  metric
}

.require_crossnobis_pairing <- function(over) {
  .validate_pairing(over)
  if (!identical(attr(over, "independence", exact = TRUE), "independent") ||
      !identical(attr(over, "estimate", exact = TRUE),
        "cross_generalized") || any(over$left == over$right)) {
    stop(paste0(
      "Crossnobis requires cross-partition evaluation edges declared ",
      "independent and containing no self-products."
    ), call. = FALSE)
  }
  invisible(TRUE)
}

.crossnobis_scientific_plan_id <- function(relation_fit_signature, task,
                                            frame, pairing,
                                            metric_schedule) {
  paste0("crossnobis-sha256:", digest::digest(list(
    schema_version = 1L,
    relation_fit = relation_fit_signature,
    evidence_task = task$task_id,
    frame = .additive_frame_signature(frame),
    pairing = .metric_pairing_identity(pairing),
    metric_schedule = metric_schedule$signature
  ), algo = "sha256", serialize = TRUE))
}

.crossnobis_plan_signature <- function(fields) {
  paste0("sha256:", digest::digest(list(
    schema_version = 1L,
    scientific_plan_id = fields$scientific_plan_id,
    compute = unclass(fields$compute),
    lowering = fields$lowering,
    kernel_version = fields$kernel_version,
    memory = unclass(fields$memory)
  ), algo = "sha256", serialize = TRUE))
}

.learned_crossnobis_memory_plan <- function(fit, frame, support_index, over,
                                            compute) {
  relation <- fit$relation
  support_index <- .validate_support_index(support_index)
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
    stop(sprintf(
      paste0("Learned crossnobis requires %.0f bytes, exceeding the ",
        "%.0f-byte workspace budget."),
      plan$planned_workspace_bytes, compute$workspace_bytes
    ), call. = FALSE)
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
#' @param x An `effect_relation_fit` with residual-block capability.
#' @param at A support-index-backed compiled spatial frame.
#' @param over Independent cross-partition evaluation edges.
#' @param metric An on-demand metric recipe such as `shrinkage_precision()`.
#' @param training A `metric_training_policy()`.
#' @param compute A sequential `compute_policy()`.
#' @param residual_workspace_bytes Positive budget used while accumulating
#'   canonical residual pair sufficient statistics. It changes cache capacity,
#'   never the canonical numerical tile shape.
#' @return An `effect_crossnobis_plan` reusable across fixed contrasts.
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
  .validate_compiler_inputs(x$relation, at, over)
  metric <- .validate_metric_recipe(metric)
  if (identical(metric$kind, "identity")) {
    stop(paste0(
      "A learned crossnobis plan requires a residual-derived precision ",
      "recipe; use `noise_precision(diag(...))` for a fixed identity metric."
    ), call. = FALSE)
  }
  training <- .validate_metric_training_policy(training)
  .preflight_metric_training(
    metric, x$relation$partitions, over, training
  )
  capabilities <- .execution_preflight(compute, function() {
    .compiler_capabilities(x$relation)
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
  if (!inherits(x, "effect_crossnobis_plan") || !is.list(x) ||
      !identical(names(x), expected) ||
      !.strong_sha256(x$relation_fit_signature) ||
      !identical(x$lowering,
        "support_streamed_scheduled_metric_query_contraction") ||
      !identical(x$kernel_version,
        "support-streamed-scheduled-metric-v1") ||
      !grepl("^crossnobis-sha256:[[:xdigit:]]{64}$",
        x$scientific_plan_id) ||
      !.strong_sha256(x$signature)) {
    stop("Crossnobis-plan fields are missing or noncanonical.", call. = FALSE)
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
    stop("Crossnobis-plan relation, support, pairing, or metric is inconsistent.",
      call. = FALSE)
  }
  if (isTRUE(deep)) {
    fields <- x[names(x) != "signature"]
    expected_scientific <- .crossnobis_scientific_plan_id(
      x$relation_fit_signature, x$task, x$frame, x$pairing,
      x$metric_schedule
    )
    if (!identical(x$scientific_plan_id, expected_scientific) ||
        !identical(x$signature, .crossnobis_plan_signature(fields))) {
      stop("Crossnobis-plan identity is inconsistent.", call. = FALSE)
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

.planned_crossnobis_receipt <- function(plan, contrast) {
  execution_receipt(
    scientific_plan_id = paste0("crossnobis-sha256:", digest::digest(list(
      schema_version = 1L,
      parent = plan$scientific_plan_id,
      contrast = contrast,
      query_role = "effect"
    ), algo = "sha256", serialize = TRUE)),
    compute = plan$compute,
    sources = plan$task$left_relation$capabilities,
    memory = plan$memory,
    kernel_version = plan$kernel_version,
    task_partition_id = "ascending-supports-one-live-node",
    reduction_plan_id = plan$task$stages$signature,
    numeric_contract = numerical_contract(),
    completion_status = "planned",
    task_count = as.double(nrow(plan$frame$weights)),
    completed_task_count = 0L,
    blas = .compiler_blas(),
    domain_signature = plan$task$left_relation$domain$signature
  )
}

.execute_learned_crossnobis <- function(plan, weights, reporter = NULL) {
  .validate_crossnobis_plan(plan)
  weights <- .align_contrast(
    weights, plan$task$left_relation$effect_space$coordinates
  )
  task <- .as_compiled_effect_task(plan$task, validate = FALSE)
  planned_receipt <- .planned_crossnobis_receipt(plan, weights)
  final_receipt <- new.env(parent = emptyenv())
  final_receipt$value <- planned_receipt
  observed <- .empty_execution_observations()
  observed$task_counts[["planned"]] <- nrow(plan$frame$weights)
  observed$tiles <- list(
    feature_block = as.integer(max(plan$metric_schedule$support_index$cost$support_size)),
    row_tile = 1L,
    coordinate_tile = 1L
  )
  observed_state <- new.env(parent = emptyenv())
  observed_state$value <- observed
  task_observer <- function(event, features) {
    observed_state$value$task_counts[[event]] <-
      observed_state$value$task_counts[[event]] + 1
    if (identical(event, "completed")) {
      observed_state$value$features_completed <-
        observed_state$value$features_completed + length(features)
    }
    invisible(NULL)
  }

  result <- .execute_guarded(
    compute = function() {
      admission_started <- proc.time()[["elapsed"]]
      source_session <- .open_effect_task_source_session(task, validate = FALSE)
      on.exit(source_session$close(), add = TRUE)
      observed_state$value$stage_seconds[["source_admission"]] <-
        max(0, proc.time()[["elapsed"]] - admission_started)
      observed_state$value$source_access <-
        source_session$summary()$access_mode
      read_relation <- function(partition, features) {
        value <- .relation_block_with_reader(
          task$left_relation, partition, features,
          function(partition, features) {
            source_session$read("left", partition, features)
          }, validate = FALSE
        )
        observed_state$value$bytes_read <-
          sum(source_session$summary()$bytes_read)
        value
      }
      kernel_started <- proc.time()[["elapsed"]]
      evaluated <- .support_streamed_scheduled_crossnobis(
        frame = plan$frame,
        metric_schedule = plan$metric_schedule,
        read_relation = read_relation,
        partitions = task$left_relation$partitions,
        effects = task$left_space$coordinates,
        ordered_edges = task$ordered_edges,
        contrast = weights,
        task_observer = task_observer
      )
      observed_state$value$stage_seconds[["support_tasks"]] <-
        max(0, proc.time()[["elapsed"]] - kernel_started)
      source_session$close()
      source_summary <- source_session$summary()
      observed_state$value$bytes_read <- sum(source_summary$bytes_read)
      list(evaluated = evaluated, source_summary = source_summary)
    },
    receipt = planned_receipt,
    reporter = .compiler_reporter(reporter, final_receipt),
    observations = function() observed_state$value,
    receipt_sink = function(receipt) final_receipt$value <- receipt
  )
  metadata <- list(
    frame = list(
      representation = plan$frame$representation,
      normalization = plan$frame$normalization,
      domain = plan$frame$domain
    ),
    metric_schedule = list(
      signature = plan$metric_schedule$signature,
      recipe = plan$metric_schedule$recipe$signature,
      recipe_kind = plan$metric_schedule$recipe$kind,
      training_policy = plan$metric_schedule$training_policy,
      records = lapply(plan$metric_schedule$records, function(record) {
        list(
          evaluation_left = record$evaluation_left,
          evaluation_right = record$evaluation_right,
          training_partitions = record$training_partitions,
          training_signature = record$training_signature
        )
      }),
      local_metric_storage = "none_derived_on_demand",
      retained_factor_table = FALSE,
      calibration_requires_metric_uncertainty =
        plan$metric_schedule$capabilities$calibration_requires_metric_uncertainty
    ),
    diagnostics = result$evaluated$diagnostics,
    metric_receipts = result$evaluated$metric_receipts,
    source_session = result$source_summary,
    execution_plan = list(
      signature = plan$signature,
      lowering = plan$lowering,
      kernel_version = plan$kernel_version,
      query_fused = TRUE,
      materialization = "direct_crossnobis_contrast"
    ),
    scientific_plan_id = final_receipt$value$scientific_plan_id
  )
  structure(list(
    values = result$evaluated$values,
    contrast = weights,
    estimand = "crossvalidated_squared_mahalanobis_contrast",
    metric = plan$metric_schedule$signature,
    pairing = .metric_pairing_identity(plan$pairing),
    index = .compiler_index(plan$frame),
    receipt = final_receipt$value,
    metadata = metadata
  ), class = "effect_crossnobis_view")
}

#' Evaluate a signed local crossnobis contrast
#'
#' The result is the query-first evidence
#' `tr(c c^T B_a K B_b^T)` aggregated over the plan's independent partition
#' edges. This function is a validating view over the ordinary geometry
#' compiler; it does not introduce a second numerical engine. Negative finite
#' estimates are retained.
#'
#' @param x Either an `effect_geometry_plan` carrying an explicit fixed
#'   `noise_precision()` metric, or an `effect_crossnobis_plan` carrying a
#'   provenance-frozen on-demand metric schedule.
#' @param weights One finite contrast weight per experimental effect.
#' @return An `effect_crossnobis_view` with one signed value per spatial
#'   measurement and the executed plan receipt.
#'
#' @section One estimand, two views:
#' On a fixed-metric geometry plan, `crossnobis(x, weights)` is the named
#' Mahalanobis reading of exactly `contrast(x, weights)$total`: the same
#' compiled estimand, exposed as a single signed value. When the analysis
#' also needs the signed endpoint marginals or the exact
#' coherent/configuration decomposition of the same quantity, call
#' [contrast()] on the same plan; every component it returns inherits the
#' plan's fixed metric.
#' @examples
#' domain <- abstract_domain(3, id = "crossnobis-example")
#' run1 <- rbind(a = c(1, 0, 0), b = c(0, 1, 0))
#' run2 <- rbind(a = c(1.1, 0, 0), b = c(0, 0.9, 0))
#' relation <- relation(list(run1 = run1, run2 = run2), domain = domain)
#' plan <- plan_geometry(
#'   relation,
#'   compile_frame(whole_brain(), domain),
#'   cross_partitions(relation, independence = "independent"),
#'   metric = noise_precision(diag(3), domain)
#' )
#' result <- crossnobis(plan, c(a = 1, b = -1))
#' result
#' as.data.frame(result)
#' @export
crossnobis <- function(x, weights) {
  if (inherits(x, "effect_crossnobis_plan")) {
    return(.execute_learned_crossnobis(x, weights))
  }
  if (!inherits(x, "effect_geometry_plan")) {
    stop("`x` must be a compiled geometry or crossnobis plan.", call. = FALSE)
  }
  metric <- .crossnobis_plan_metric(x)
  .require_crossnobis_pairing(x$pairing)
  weights <- .align_contrast(
    weights, x$task$left_relation$effect_space$coordinates
  )
  query <- bilinear_query(tcrossprod(weights))
  evaluated <- evaluate_geometry(x, query = query, component = "total")
  values <- drop(evaluated$values)
  structure(list(
    values = values,
    contrast = weights,
    estimand = "crossvalidated_squared_mahalanobis_contrast",
    metric = metric$signature,
    pairing = .metric_pairing_identity(x$pairing),
    index = evaluated$index,
    receipt = evaluated$receipt,
    metadata = evaluated$metadata
  ), class = "effect_crossnobis_view")
}
