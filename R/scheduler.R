# In-process bounded scheduling model ---------------------------------------

.simulate_bounded_schedule <- function(values, task_id, result_bytes,
                                       completion_order = task_id,
                                       max_inflight = 1L,
                                       max_reorder_bytes = Inf,
                                       fail_task_id = NULL) {
  n <- length(values)
  if (!is.list(values) || n < 1L || length(task_id) != n ||
      length(result_bytes) != n || anyNA(task_id) || anyDuplicated(task_id)) {
    stop("Tasks require unique ids, values, and result-byte estimates.",
      call. = FALSE)
  }
  max_inflight <- .validate_tile_size(max_inflight, "max_inflight")
  if (!is.numeric(result_bytes) || anyNA(result_bytes) ||
      any(!is.finite(result_bytes)) || any(result_bytes < 0)) {
    stop("`result_bytes` must be finite nonnegative estimates.", call. = FALSE)
  }
  if (!is.numeric(max_reorder_bytes) || length(max_reorder_bytes) != 1L ||
      is.na(max_reorder_bytes) || max_reorder_bytes <= 0) {
    stop("`max_reorder_bytes` must be one positive number.", call. = FALSE)
  }
  if (any(result_bytes > max_reorder_bytes)) {
    stop("A task result exceeds `max_reorder_bytes` and cannot be dispatched.",
      call. = FALSE)
  }
  if (!setequal(completion_order, task_id) || length(completion_order) != n) {
    stop("`completion_order` must be a permutation of task ids.", call. = FALSE)
  }
  valid_values <- vapply(values, function(x) is.numeric(x) && all(is.finite(x)), logical(1))
  if (!all(valid_values)) stop("Task values must be finite numeric objects.", call. = FALSE)

  canonical <- order(task_id, method = "radix")
  task_id <- task_id[canonical]
  values <- values[canonical]
  result_bytes <- result_bytes[canonical]
  priority <- match(task_id, completion_order)
  next_dispatch <- 1L
  next_reduce <- 1L
  inflight <- integer()
  buffer <- list()
  reduced <- NULL
  release_order <- task_id[0]
  max_inflight_seen <- 0L
  max_reserved_bytes <- 0
  max_buffer_bytes <- 0
  backpressure_count <- 0L

  repeat {
    dispatched <- FALSE
    while (next_dispatch <= n && length(inflight) < max_inflight) {
      occupied <- sum(result_bytes[inflight]) +
        sum(result_bytes[as.integer(names(buffer))])
      if (occupied + result_bytes[[next_dispatch]] > max_reorder_bytes) {
        backpressure_count <- backpressure_count + 1L
        break
      }
      inflight <- c(inflight, next_dispatch)
      next_dispatch <- next_dispatch + 1L
      dispatched <- TRUE
      max_inflight_seen <- max(max_inflight_seen, length(inflight))
      max_reserved_bytes <- max(max_reserved_bytes,
        sum(result_bytes[inflight]) + sum(result_bytes[as.integer(names(buffer))]))
    }

    while (as.character(next_reduce) %in% names(buffer)) {
      value <- buffer[[as.character(next_reduce)]]
      reduced <- if (is.null(reduced)) value else reduced + value
      buffer[[as.character(next_reduce)]] <- NULL
      release_order <- c(release_order, task_id[[next_reduce]])
      next_reduce <- next_reduce + 1L
    }
    if (next_reduce > n) break
    if (length(inflight) == 0L && next_dispatch <= n) next
    if (length(inflight) == 0L) {
      stop("Scheduler deadlock: budget prevents progress.", call. = FALSE)
    }

    finished <- inflight[[which.min(priority[inflight])]]
    inflight <- inflight[inflight != finished]
    if (!is.null(fail_task_id) && identical(task_id[[finished]], fail_task_id)) {
      condition <- structure(
        list(
          message = paste0("simulated task failure: ", fail_task_id),
          call = NULL,
          cleanup_status = list(
            success = TRUE,
            inflight_released = length(inflight),
            buffered_released = length(buffer)
          ),
          diagnostics = list(max_inflight_seen = max_inflight_seen,
            max_reserved_bytes = max_reserved_bytes)
        ),
        class = c("effect_scheduler_error", "error", "condition")
      )
      stop(condition)
    }
    buffer[[as.character(finished)]] <- values[[finished]]
    max_buffer_bytes <- max(max_buffer_bytes,
      sum(result_bytes[as.integer(names(buffer))]))
  }

  list(
    value = reduced,
    diagnostics = list(
      max_inflight_seen = max_inflight_seen,
      max_reserved_bytes = max_reserved_bytes,
      max_buffer_bytes = max_buffer_bytes,
      backpressure_count = backpressure_count,
      release_order = release_order,
      dispatched_all = next_dispatch > n
    )
  )
}
