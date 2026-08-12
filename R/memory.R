# Conservative execution-memory planning ------------------------------------

#' Construct a conservative peak-memory plan
#'
#' The plan names every allocation category required by the execution design.
#' It is an upper-bound model, not a promise of exact portable prediction: R,
#' BLAS, allocators, and operating systems may add unobserved overhead. Plans
#' should therefore be checked against measured aggregate peak memory.
#'
#' @param shared_source_bytes Source memory shared across workers.
#' @param runtime_reserve_bytes Fixed process-level reserve for R, Matrix/BLAS
#'   lazy initialization, allocator arenas, and other measured runtime pages
#'   that are not attributable to one scientific array. The default 64 MiB is
#'   intentionally conservative and must remain benchmark-validated.
#' @param private_source_bytes Source memory privately held by each worker.
#' @param source_block_bytes,relation_block_bytes,atom_block_bytes Working
#'   blocks retained during one contraction task.
#' @param output_bytes Durable in-memory output or active output tiles.
#' @param contraction_bytes Largest dense sparse-contraction result temporary.
#' @param replacement_copy_bytes Allowance for R copy-on-modify/replacement.
#' @param serialization_overlap_bytes Data simultaneously present while being
#'   serialized or transferred.
#' @param reorder_buffer_bytes Maximum bounded out-of-order result buffer.
#' @param checkpoint_buffer_bytes Maximum checkpoint staging buffer.
#' @param workers Positive worker count represented by the plan.
#' @param n_active Maximum concurrently active tasks. It cannot exceed workers.
#' @param safety_factor Multiplicative allocator/implementation headroom, at
#'   least one.
#' @param budget_bytes Optional positive memory budget.
#' @param measured_peak_bytes Optional observed aggregate peak RSS for
#'   validation.
#' @return A categorized conservative memory plan.
#' @export
memory_plan <- function(shared_source_bytes = 0,
                        runtime_reserve_bytes = 64 * 1024^2,
                        private_source_bytes = 0,
                        source_block_bytes = 0,
                        relation_block_bytes = 0,
                        atom_block_bytes = 0,
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
                        measured_peak_bytes = NULL) {
  categories <- c(
    shared_source = shared_source_bytes,
    runtime_reserve = runtime_reserve_bytes,
    private_source_per_worker = private_source_bytes,
    source_block = source_block_bytes,
    relation_block = relation_block_bytes,
    atom_block = atom_block_bytes,
    output = output_bytes,
    contraction = contraction_bytes,
    replacement_copy = replacement_copy_bytes,
    serialization_overlap = serialization_overlap_bytes,
    reorder_buffer = reorder_buffer_bytes,
    checkpoint_buffer = checkpoint_buffer_bytes
  )
  max_exact <- 2^53
  validate_whole <- function(value, name, positive = FALSE) {
    lower <- if (positive) 1 else 0
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < lower || value %% 1 != 0 || value > max_exact) {
      stop(sprintf("`%s` must be one %sfinite whole scalar no greater than 2^53.",
        name, if (positive) "positive " else "nonnegative "), call. = FALSE)
    }
    value
  }
  for (name in names(categories)) validate_whole(categories[[name]], name)
  workers <- validate_whole(workers, "workers", positive = TRUE)
  n_active <- validate_whole(n_active, "n_active", positive = TRUE)
  if (n_active > workers) {
    stop("`n_active` cannot exceed `workers`.", call. = FALSE)
  }
  if (!is.numeric(safety_factor) || length(safety_factor) != 1L ||
      is.na(safety_factor) || !is.finite(safety_factor) || safety_factor < 1) {
    stop("`safety_factor` must be one finite number greater than or equal to one.",
      call. = FALSE)
  }
  validate_optional_bytes <- function(value, name) {
    if (!is.null(value) &&
        (!is.numeric(value) || length(value) != 1L || is.na(value) ||
         !is.finite(value) || value <= 0)) {
      stop(sprintf("`%s` must be NULL or one positive finite byte count.", name),
        call. = FALSE)
    }
  }
  validate_optional_bytes(budget_bytes, "budget_bytes")
  validate_optional_bytes(measured_peak_bytes, "measured_peak_bytes")

  safe_product <- function(left, right, name) {
    value <- left * right
    if (!is.finite(value) || value > max_exact) {
      stop(sprintf("Memory product `%s` overflows exact byte accounting.", name),
        call. = FALSE)
    }
    value
  }
  task_names <- c(
    "source_block", "relation_block", "atom_block", "contraction",
    "replacement_copy", "serialization_overlap", "reorder_buffer",
    "checkpoint_buffer"
  )
  task_private_per_active <- sum(categories[task_names])
  if (!is.finite(task_private_per_active) || task_private_per_active > max_exact) {
    stop("Task-private byte categories overflow exact accounting.", call. = FALSE)
  }
  worker_private <- safe_product(workers,
    categories[["private_source_per_worker"]], "workers * private_source")
  active_task_total <- safe_product(n_active, task_private_per_active,
    "n_active * task_private")
  modeled <- categories[["shared_source"]] + categories[["runtime_reserve"]] +
    categories[["output"]] +
    worker_private + active_task_total
  if (!is.finite(modeled) || modeled > max_exact) {
    stop("Modeled peak overflows exact byte accounting.", call. = FALSE)
  }
  conservative <- ceiling(modeled * safety_factor)

  structure(
    list(
      categories = categories,
      workers = as.integer(workers),
      n_active = as.integer(n_active),
      worker_private_total_bytes = worker_private,
      task_private_per_active_bytes = task_private_per_active,
      active_task_total_bytes = active_task_total,
      modeled_peak_bytes = modeled,
      safety_factor = safety_factor,
      conservative_peak_bytes = conservative,
      budget_bytes = budget_bytes,
      fits_budget = if (is.null(budget_bytes)) NA else conservative <= budget_bytes,
      measured_peak_bytes = measured_peak_bytes,
      measurement_within_plan = if (is.null(measured_peak_bytes)) NA else
        measured_peak_bytes <= conservative,
      prediction_kind = "conservative_upper_bound_to_validate"
    ),
    class = "effect_memory_plan"
  )
}
