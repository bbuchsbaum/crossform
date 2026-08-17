# The crossnobis runtime and its public entry --------------------------------
#
# `R/crossnobis.R` is the plan: `noise_precision()`, `plan_crossnobis()`, the
# plan validator, its print method. It used to be the executor as well, which
# made a layer-3 plan file call `.run_geometry_compiler()` in
# R/execution-driver.R and `.support_streamed_scheduled_crossnobis()` in
# R/kernel.R -- the last two entries in the architecture register, and the same
# category error the compiler had before the executor was split out of it in
# wave 4B.
#
# The split follows that precedent exactly. What lives here is everything that
# runs: the planned receipt, the learned-metric runtime, and the exported
# `crossnobis()` entry that dispatches between the learned route and the
# ordinary geometry compiler. The plan constructs nothing here and this file
# constructs no plan; the direction is plan -> runtime, one way.

.planned_crossnobis_receipt <- function(plan, contrast) {
  execution_receipt(
    scientific_plan_id = .sha256_signature(list(
      schema_version = 1L,
      parent = plan$scientific_plan_id,
      contrast = contrast,
      query_role = "effect"
    ), "crossnobis-sha256:"),
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
    blas = .execution_blas_record(),
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
  # The same observer and the same stage clock the geometry executor uses; now
  # that both runtimes are in layer 4 this is a sideways call rather than a
  # second copy of the accounting.
  task_observer <- .report_execution_event(observed_state)

  result <- .execute_guarded(
    compute = function() {
      admission_started <- proc.time()[["elapsed"]]
      source_session <- .open_effect_task_source_session(task, validate = FALSE)
      on.exit(source_session$close(), add = TRUE)
      .record_execution_stage(
        observed_state, "source_admission", admission_started
      )
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
      .record_execution_stage(observed_state, "support_tasks", kernel_started)
      source_session$close()
      source_summary <- source_session$summary()
      observed_state$value$bytes_read <- sum(source_summary$bytes_read)
      list(evaluated = evaluated, source_summary = source_summary)
    },
    receipt = planned_receipt,
    reporter = .execution_reporter(reporter, final_receipt),
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
    index = .execution_measurement_index(plan$frame),
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
#' @return An `effect_crossnobis_view` whose `$values` holds one signed
#'   crossvalidated squared Mahalanobis value per spatial measurement, with
#'   the aligned `$contrast`, the named `$estimand`, the `$metric` and
#'   `$pairing` identities, `$index`, and the executed `$receipt`.
#' @section Structure:
#' One signed value per spatial measurement, alongside the declarations that
#' make it a Mahalanobis reading.
#'
#' - `$values`: the signed crossvalidated squared Mahalanobis contrast, one
#'   per measurement, in `$index` order. Negative estimates are retained.
#' - `$contrast`: the weights, reordered to the relation's effect order and
#'   named.
#' - `$estimand`: the named estimand,
#'   `"crossvalidated_squared_mahalanobis_contrast"`.
#' - `$metric`: the signature of the noise-precision metric the values were
#'   read under.
#' - `$pairing`: the identity of the independent partition edges they were
#'   generalized over.
#' - `$index`: the measurement identifiers, one per value, carried from the
#'   frame's `$index$measurement`.
#' - `$receipt`: the execution receipt for the run that produced the values.
#'
#' The `$metadata` block and any other element not listed here are internal
#' and may change.
#' @seealso [noise_precision()] for the fixed metric a geometry plan needs,
#'   [plan_crossnobis()] for the learned-metric route, and
#'   [contrast_energy()] for the decomposed reading of the same estimand.
#' @family geometry plans and views
#'
#' @section Refusal:
#' A plan carrying the implicit identity metric, or a fixed metric that was not
#' built by [noise_precision()], signals an `effect_capability_refusal` with
#' capability `"declared_noise_metric"` in namespace `"geometry_views"`.
#' Inspect it with [catch_refusal()].
#'
#' @section One estimand, two views:
#' On a fixed-metric geometry plan, `crossnobis(x, weights)` is the named
#' Mahalanobis reading of exactly `contrast_energy(x, weights)$total`: the same
#' compiled estimand, exposed as a single signed value. When the analysis
#' also needs the signed endpoint marginals or the exact
#' coherent/configuration decomposition of the same quantity, call
#' [contrast_energy()] on the same plan; every component it returns inherits the
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
    .input_error(sprintf(paste0(
      "`x` must be an `effect_geometry_plan` from `plan_geometry()` or an ",
      "`effect_crossnobis_plan` from `plan_crossnobis()`; received %s."
    ), .msg_value(x)))
  }
  if (missing(weights)) {
    .input_error(paste0(
      "`weights` is required: pass one finite weight per experimental ",
      "effect, for example `crossnobis(plan, c(face = 1, house = -1))`."
    ))
  }
  metric <- .crossnobis_plan_metric(x)
  .require_crossnobis_pairing(x$pairing)
  weights <- .align_contrast(
    weights, x$task$left_relation$effect_space$coordinates
  )
  query <- bilinear_query(tcrossprod(weights))
  # The internal runner rather than `evaluate_geometry()`: this is the plan
  # layer, and the public entry point sits above the views it feeds.
  evaluated <- .run_geometry_compiler(x, query = query, component = "total")
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
