# Conservative owned-workspace planning -------------------------------------

# Exact size of a base-R double matrix on the current R runtime.
# Keeping this arithmetic here is deliberate: a memory preflight must not
# allocate the buffer whose affordability it is deciding.  The fixed overhead
# is measured once from a zero-payload skeleton rather than assumed.
.dense_double_matrix_bytes <- function(rows, columns) {
  dimensions <- c(rows = rows, columns = columns)
  if (!is.numeric(dimensions) || length(dimensions) != 2L ||
      anyNA(dimensions) || any(!is.finite(dimensions)) ||
      any(dimensions < 0) || any(dimensions %% 1 != 0)) {
    stop("Matrix dimensions must be finite nonnegative whole scalars.",
      call. = FALSE)
  }
  max_exact <- 2^53
  overhead <- as.double(utils::object.size(
    structure(numeric(), dim = c(0L, 0L))
  ))
  if (rows != 0 && columns > floor((max_exact - overhead) / (8 * rows))) {
    stop("Matrix byte accounting overflows exact representation.",
      call. = FALSE)
  }
  8 * rows * columns + overhead
}

.dense_double_vector_bytes <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      !is.finite(n) || n < 0 || n %% 1 != 0) {
    stop("Vector length must be one finite nonnegative whole scalar.",
      call. = FALSE)
  }
  max_exact <- 2^53
  overhead <- as.double(utils::object.size(numeric()))
  if (n > floor((max_exact - overhead) / 8)) {
    stop("Vector byte accounting overflows exact representation.",
      call. = FALSE)
  }
  8 * n + overhead
}

# An array's values may be enormous while its axis metadata is small.  Build
# only a zero-payload skeleton so named-array accounting remains exact without
# materializing the values being planned.
.named_double_array_bytes <- function(dimensions, dimnames = NULL) {
  if (!is.numeric(dimensions) || length(dimensions) < 1L ||
      anyNA(dimensions) || any(!is.finite(dimensions)) ||
      any(dimensions < 0) || any(dimensions %% 1 != 0)) {
    stop("Array dimensions must be finite nonnegative whole scalars.",
      call. = FALSE)
  }
  max_exact <- 2^53
  elements <- prod(dimensions)
  if (!is.finite(elements) || elements > floor(max_exact / 8)) {
    stop("Array byte accounting overflows exact representation.",
      call. = FALSE)
  }
  skeleton_dimensions <- dimensions
  skeleton_dimensions[[1L]] <- 0
  skeleton <- array(numeric(), skeleton_dimensions, dimnames = dimnames)
  bytes <- as.double(utils::object.size(skeleton)) + 8 * elements
  if (!is.finite(bytes) || bytes > max_exact) {
    stop("Array byte accounting overflows exact representation.",
      call. = FALSE)
  }
  bytes
}

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
#' @keywords internal
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

# Conservative plan for the universal rectangular/packed streaming primitive.
# Matrix object sizes, rather than payload-only byte counts, keep the plan
# comparable to the kernel's named-live-object diagnostics.
.effect_form_kernel_memory_plan <- function(
    frame, left_effects, right_effects, left_partitions, right_partitions,
    codec = c("rectangular", "symmetric_packed"), query = NULL,
    same_relation = FALSE, feature_block = 1024L, row_tile = 1024L,
    coordinate_tile = 256L, storage = c("memory", "block"),
    retain_first_moments = FALSE, form_total = TRUE) {
  .validate_frame_for_compile(frame)
  codec <- match.arg(codec)
  storage <- match.arg(storage)
  feature_block <- .validate_tile_size(feature_block, "feature_block")
  row_tile <- .validate_tile_size(row_tile, "row_tile")
  coordinate_tile <- .validate_tile_size(
    coordinate_tile, "coordinate_tile"
  )
  if (!is.logical(same_relation) || length(same_relation) != 1L ||
      is.na(same_relation)) {
    stop("`same_relation` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(retain_first_moments) ||
      length(retain_first_moments) != 1L || is.na(retain_first_moments) ||
      !is.logical(form_total) || length(form_total) != 1L ||
      is.na(form_total)) {
    stop("First-moment and total planning flags must be TRUE or FALSE.",
      call. = FALSE)
  }
  if (!retain_first_moments && !form_total) {
    stop("A plan must retain first moments, form total output, or both.",
      call. = FALSE)
  }
  .validate_effect_names(left_effects, length(left_effects))
  .validate_effect_names(right_effects, length(right_effects))
  if (!is.character(left_partitions) || length(left_partitions) < 1L ||
      anyNA(left_partitions) || any(!nzchar(left_partitions)) ||
      anyDuplicated(left_partitions) ||
      !is.character(right_partitions) || length(right_partitions) < 1L ||
      anyNA(right_partitions) || any(!nzchar(right_partitions)) ||
      anyDuplicated(right_partitions)) {
    stop("Effect-form partition axes must be unique nonempty identifiers.",
      call. = FALSE)
  }
  if (same_relation && (!identical(left_partitions, right_partitions) ||
      !identical(left_effects, right_effects))) {
    stop("A shared relation plan requires identical partition and effect axes.",
      call. = FALSE)
  }

  q_left <- length(left_effects)
  q_right <- length(right_effects)
  physical_width <- if (codec == "rectangular") {
    q_left * q_right
  } else {
    if (!same_relation || !identical(left_effects, right_effects)) {
      stop("Symmetric-packed planning requires a self-form.", call. = FALSE)
    }
    q_left * (q_left + 1L) / 2L
  }
  if (!is.null(query) && (!is.matrix(query) || !is.numeric(query) ||
      nrow(query) != physical_width || ncol(query) < 1L ||
      any(!is.finite(query)))) {
    stop("`query` must match the finite physical form coordinates.",
      call. = FALSE)
  }

  features <- ncol(frame$weights)
  measurements <- nrow(frame$weights)
  output_width <- if (is.null(query)) physical_width else ncol(query)
  f <- min(feature_block, features)
  rows <- min(row_tile, measurements)
  coordinates <- min(coordinate_tile, output_width)
  relation_block <- length(left_partitions) *
    .dense_double_matrix_bytes(q_left, f) +
    if (same_relation) 0 else
      length(right_partitions) *
        .dense_double_matrix_bytes(q_right, f)
  atom_matrix <- if (form_total) {
    .dense_double_matrix_bytes(f, output_width)
  } else 0
  atom_work <- max(
    if (form_total) .dense_double_vector_bytes(f) else 0,
    if (!form_total || is.null(query)) 0 else
      .dense_double_matrix_bytes(q_left, q_right)
  )
  weight_slice <- .dense_double_matrix_bytes(rows, f)
  atom_slice <- if (form_total) {
    .dense_double_matrix_bytes(f, coordinates)
  } else 0
  product <- if (form_total) {
    .dense_double_matrix_bytes(rows, coordinates)
  } else 0
  first_product <- if (retain_first_moments) {
    max(.dense_double_matrix_bytes(rows, q_left),
      .dense_double_matrix_bytes(rows, q_right))
  } else {
    0
  }
  first_state <- if (retain_first_moments) {
    .named_double_array_bytes(
      c(measurements, q_left, length(left_partitions)),
      list(NULL, left_effects, left_partitions)
    ) + if (same_relation) 0 else .named_double_array_bytes(
      c(measurements, q_right, length(right_partitions)),
      list(NULL, right_effects, right_partitions)
    )
  } else {
    0
  }
  output <- if (form_total && storage == "memory") {
    .dense_double_matrix_bytes(measurements, output_width)
  } else {
    0
  }

  memory_plan(
    frame_bytes = as.double(utils::object.size(frame$weights)),
    relation_block_bytes = relation_block,
    atom_block_bytes = atom_matrix + atom_work,
    local_state_bytes = first_state,
    output_bytes = output,
    contraction_bytes = weight_slice + atom_slice + product + first_product,
    replacement_copy_bytes =
      (if (form_total && storage == "memory") 2 * product else 0) +
      (if (retain_first_moments) 2 * first_product else 0),
    safety_factor = 1.25
  )
}
