# Admission for model fitting/readout workspace (layer 4) ------------------
# The ordinary geometry compiler still plans and executes neural work.
# Reserve model/readout buffers first, then give that compiler the remainder.

.geometry_prediction_extra_memory <- function(prepared, budget, storage, rows,
                                              operation, fit = NULL, modes = FALSE,
                                              components = FALSE, frozen_prediction_bytes = NULL) {
  r <- as.double(prepared$basis$dimension)
  n <- as.double(prepared$plan$measurements)
  k <- as.double(budget)
  h <- r * (r + 1) / 2
  width <- if (operation == "fit") r * k + k + 2 * r else if (modes) 4 * k else 0
  # Factors are never stored as redundant model-coordinate square matrices.
  retained <- 8 * n * width * (storage == "memory") + 512 * n
  model <- 2 * as.double(utils::object.size(prepared$basis))
  frozen <- if (is.null(fit)) 0 else 512 * n +
    if (identical(fit$values$representation, "memory")) 8 * prod(as.double(fit$values$dim)) else 0
  if (!is.null(frozen_prediction_bytes)) frozen <- max(frozen, frozen_prediction_bytes)
  scratch <- 8 * (32 * r^2 + 20 * r + 12 * r * k) + 65536
  blocks <- rows * (8 * ((if (components) 4 else 2) * h +
    4 * r * k + 12 * k + 2 * r) + 4096)
  # Serialization and copy-on-modify overlap include metadata and output.
  c(model = model, frozen_prediction = frozen, retained_output = retained,
    spectral_scratch = scratch, row_buffers = blocks,
    serialization_overlap = model + retained)
}

.geometry_prediction_workspace <- function(prepared, budget, storage, row_block,
                                           operation, fit = NULL, modes = FALSE,
                                           components = FALSE, cache_bytes = 0,
                                           frozen_prediction_bytes = NULL) {
  .check_number(cache_bytes, "cache_bytes", nonnegative = TRUE)
  if (!is.null(frozen_prediction_bytes)) .check_number(frozen_prediction_bytes, "frozen_prediction_bytes", nonnegative = TRUE)
  plan <- prepared$plan
  policy <- .validate_compute_policy(plan$compute)
  if (storage == "block" && !identical(plan$metric_schedule$kind, "implicit_identity_before_frame")) {
    .capability_refusal("The complete fixed-metric executor currently supports memory output only.",
      capability = "predictive_storage_route", namespace = "predictive_geometry",
      reasons = "fixed_metric_block_storage_not_admitted",
      remedies = "Use storage = 'memory' with an adequate workspace budget for a fixed neural metric.")
  }
  maximum <- min(row_block, plan$measurements)
  choices <- if (is.null(policy$workspace_bytes)) maximum else .descending_tiles(maximum)
  for (rows in choices) {
    categories <- c(.geometry_prediction_extra_memory(prepared, budget, storage, rows,
      operation, fit, modes, components, frozen_prediction_bytes), training_cache = cache_bytes)
    extra <- 1.25 * sum(categories)
    if (!is.finite(extra) || extra > 2^53) .input_error("Predictive workspace dimensions overflow exact byte accounting.")
    remaining <- if (is.null(policy$workspace_bytes)) NULL else policy$workspace_bytes - extra
    if (!is.null(remaining) && remaining <= 0) next
    candidate <- if (is.null(remaining)) plan else plan_geometry(plan$task$left_relation,
      plan$frame, plan$pairing, metric = plan$metric_schedule$metric,
      compute = compute_policy(block_features = policy$block_features, workspace_bytes = remaining))
    execution <- tryCatch(.compile_geometry_execution_plan(candidate, storage = storage),
      effect_input_error = function(error) {
        if (startsWith(conditionMessage(error), "The conservative memory plan requires")) NULL else stop(error)
      })
    if (is.null(execution)) next
    record <- list(operation = operation, storage = storage, row_block = as.integer(rows),
      categories = categories, safety_factor = 1.25, extra_bytes = extra,
      geometry = execution$memory, planned_workspace_bytes = extra + execution$memory$planned_workspace_bytes,
      budget_bytes = policy$workspace_bytes, workers = policy$workers)
    record$signature <- .sha256_signature(record, "prediction-workspace-sha256:")
    prepared$plan <- candidate
    prepared$workspace <- record
    return(prepared)
  }
  .capability_refusal("The declared workspace cannot hold the reduced geometry and predictive buffers, even one measurement row at a time.",
    capability = "predictive_workspace", namespace = "predictive_geometry",
    reasons = "predictive_workspace_budget_exceeded",
    remedies = "Increase workspace_bytes, reduce model/response ranks, or request admitted block storage.")
}

.validate_prediction_workspace <- function(x, receipt) {
  fields <- c("operation", "storage", "row_block", "categories", "safety_factor", "extra_bytes",
    "geometry", "planned_workspace_bytes", "budget_bytes", "workers", "signature")
  if (!is.list(x) || !identical(names(x), fields) ||
      !x$operation %in% c("fit", "score") || !x$storage %in% c("memory", "block") ||
      !.is_count(x$row_block) || !.is_finite_numeric(x$categories) || any(x$categories < 0) ||
      !identical(x$workers, 1L) || !identical(x$safety_factor, 1.25) ||
      !identical(x$extra_bytes, x$safety_factor * sum(x$categories)) ||
      !identical(x$geometry, receipt$memory) ||
      !identical(x$planned_workspace_bytes, x$extra_bytes + x$geometry$planned_workspace_bytes) ||
      (!is.null(x$budget_bytes) && (!.is_number(x$budget_bytes) ||
        x$planned_workspace_bytes > x$budget_bytes))) {
    .contract_error("Predictive workspace accounting is incomplete or inconsistent with execution.")
  }
  .validate_memory_plan_for_receipt(x$geometry)
  .check_signature(x$signature, .sha256_signature(x[setdiff(fields, "signature")], "prediction-workspace-sha256:"),
    "Predictive workspace identity is inconsistent.")
  invisible(x)
}
