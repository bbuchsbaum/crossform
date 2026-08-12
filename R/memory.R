# Conservative owned-workspace planning -------------------------------------

#' Construct an effectagram-owned workspace plan
#'
#' The hard budget covers package-owned live objects and conservative temporary
#' overlap. Process baseline RSS and absolute RSS are observations, not hidden
#' additions to the workspace model.
#'
#' @param frame_bytes Resident dense or sparse frame storage.
#' @param resident_source_bytes Resident source objects owned by effectagram.
#' @param source_handle_bytes Execution-owned source handles and bookkeeping.
#' @param source_block_bytes,relation_block_bytes,atom_block_bytes Per-active-task
#'   data blocks.
#' @param local_state_bytes Durable component-dependent local relation state.
#' @param output_bytes Durable in-memory result components.
#' @param contraction_bytes,replacement_copy_bytes Per-task contraction and
#'   copy-on-modify temporaries.
#' @param serialization_overlap_bytes,reorder_buffer_bytes,checkpoint_buffer_bytes
#'   Optional per-task execution buffers.
#' @param workers,n_active Positive worker and active-task counts.
#' @param safety_factor Multiplicative headroom, at least one.
#' @param budget_bytes Optional hard owned-workspace budget.
#' @param measured_workspace_bytes Optional observed package-owned live bytes.
#' @param baseline_rss_bytes,peak_rss_bytes Optional process RSS observations.
#' @return An immutable-by-convention `effect_memory_plan`.
#' @export
memory_plan <- function(frame_bytes = 0,
                        resident_source_bytes = 0,
                        source_handle_bytes = 0,
                        source_block_bytes = 0,
                        relation_block_bytes = 0,
                        atom_block_bytes = 0,
                        local_state_bytes = 0,
                        output_bytes = 0,
                        contraction_bytes = 0,
                        replacement_copy_bytes = 0,
                        serialization_overlap_bytes = 0,
                        reorder_buffer_bytes = 0,
                        checkpoint_buffer_bytes = 0,
                        workers = 1L,
                        n_active = workers,
                        safety_factor = 1.25,
                        budget_bytes = NULL,
                        measured_workspace_bytes = NULL,
                        baseline_rss_bytes = NULL,
                        peak_rss_bytes = NULL) {
  categories <- c(
    frame = frame_bytes,
    resident_source = resident_source_bytes,
    source_handles = source_handle_bytes,
    source_block = source_block_bytes,
    relation_block = relation_block_bytes,
    atom_block = atom_block_bytes,
    local_state = local_state_bytes,
    output = output_bytes,
    contraction = contraction_bytes,
    replacement_copy = replacement_copy_bytes,
    serialization_overlap = serialization_overlap_bytes,
    reorder_buffer = reorder_buffer_bytes,
    checkpoint_buffer = checkpoint_buffer_bytes
  )
  max_exact <- 2^53
  whole <- function(value, name, positive = FALSE) {
    lower <- if (positive) 1 else 0
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < lower || value %% 1 != 0 ||
        value > max_exact) {
      stop(sprintf("`%s` must be one %sfinite whole scalar no greater than 2^53.",
        name, if (positive) "positive " else "nonnegative "), call. = FALSE)
    }
    value
  }
  for (name in names(categories)) whole(categories[[name]], name)
  workers <- whole(workers, "workers", positive = TRUE)
  n_active <- whole(n_active, "n_active", positive = TRUE)
  if (n_active > workers) stop("`n_active` cannot exceed `workers`.", call. = FALSE)
  if (!is.numeric(safety_factor) || length(safety_factor) != 1L ||
      is.na(safety_factor) || !is.finite(safety_factor) || safety_factor < 1) {
    stop("`safety_factor` must be one finite number greater than or equal to one.",
      call. = FALSE)
  }
  optional <- function(value, name) {
    if (!is.null(value) && (!is.numeric(value) || length(value) != 1L ||
        is.na(value) || !is.finite(value) || value < 0)) {
      stop(sprintf("`%s` must be NULL or one nonnegative finite byte count.",
        name), call. = FALSE)
    }
    value
  }
  optional(budget_bytes, "budget_bytes")
  if (!is.null(budget_bytes) && budget_bytes == 0) {
    stop("`budget_bytes` must be NULL or positive.", call. = FALSE)
  }
  optional(measured_workspace_bytes, "measured_workspace_bytes")
  optional(baseline_rss_bytes, "baseline_rss_bytes")
  optional(peak_rss_bytes, "peak_rss_bytes")
  if (!is.null(baseline_rss_bytes) && !is.null(peak_rss_bytes) &&
      peak_rss_bytes < baseline_rss_bytes) {
    stop("Peak RSS cannot be smaller than baseline RSS.", call. = FALSE)
  }

  persistent_names <- c("frame", "resident_source", "source_handles",
    "local_state", "output")
  task_names <- c("source_block", "relation_block", "atom_block",
    "contraction", "replacement_copy", "serialization_overlap",
    "reorder_buffer", "checkpoint_buffer")
  persistent <- sum(categories[persistent_names])
  per_active <- sum(categories[task_names])
  active <- n_active * per_active
  if (!is.finite(active) || active > max_exact ||
      !is.finite(persistent + active) || persistent + active > max_exact) {
    stop("Workspace byte accounting overflows exact representation.",
      call. = FALSE)
  }
  modeled <- persistent + active
  conservative <- ceiling(modeled * safety_factor)
  incremental_rss <- if (is.null(baseline_rss_bytes) || is.null(peak_rss_bytes)) {
    NULL
  } else {
    peak_rss_bytes - baseline_rss_bytes
  }

  structure(list(
    categories = categories,
    workers = as.integer(workers),
    n_active = as.integer(n_active),
    persistent_workspace_bytes = persistent,
    task_workspace_per_active_bytes = per_active,
    active_task_workspace_bytes = active,
    modeled_workspace_bytes = modeled,
    safety_factor = safety_factor,
    planned_workspace_bytes = conservative,
    budget_bytes = budget_bytes,
    fits_budget = if (is.null(budget_bytes)) NA else conservative <= budget_bytes,
    measured_workspace_bytes = measured_workspace_bytes,
    measured_workspace_within_plan = if (is.null(measured_workspace_bytes)) NA else
      measured_workspace_bytes <= conservative,
    baseline_rss_bytes = baseline_rss_bytes,
    incremental_peak_rss_bytes = incremental_rss,
    absolute_peak_rss_bytes = peak_rss_bytes,
    prediction_kind = "effectagram_owned_workspace_upper_bound"
  ), class = "effect_memory_plan")
}
