# Memory-plan declarations ---------------------------------------------------
#
# Layer 2 (values). A memory plan is a statement about bytes, computed from
# shapes and policies; nothing here allocates the buffers it is deciding on.
# That is why values, plans, and the executor may all build one.

# Exact size of a base-R double matrix on the current R runtime.
# Keeping this arithmetic here is deliberate: a memory preflight must not
# allocate the buffer whose affordability it is deciding.  The fixed overhead
# is measured once from a zero-payload skeleton rather than assumed.
.dense_double_matrix_bytes <- function(rows, columns) {
  dimensions <- c(rows = rows, columns = columns)
  if (!.is_finite_numeric(dimensions) || length(dimensions) != 2L ||
      anyNA(dimensions) || any(dimensions < 0) || any(dimensions %% 1 != 0)) {
    .input_error("Matrix dimensions must be finite nonnegative whole scalars.")
  }
  max_exact <- 2^53
  overhead <- as.double(utils::object.size(
    structure(numeric(), dim = c(0L, 0L))
  ))
  if (rows != 0 && columns > floor((max_exact - overhead) / (8 * rows))) {
    .invariant_error("Matrix byte accounting overflows exact representation.")
  }
  8 * rows * columns + overhead
}

.dense_double_vector_bytes <- function(n) {
  if (!.is_number(n) || n < 0 || n %% 1 != 0) {
    .input_error("Vector length must be one finite nonnegative whole scalar.")
  }
  max_exact <- 2^53
  overhead <- as.double(utils::object.size(numeric()))
  if (n > floor((max_exact - overhead) / 8)) {
    .invariant_error("Vector byte accounting overflows exact representation.")
  }
  8 * n + overhead
}

# An array's values may be enormous while its axis metadata is small.  Build
# only a zero-payload skeleton so named-array accounting remains exact without
# materializing the values being planned.
.named_double_array_bytes <- function(dimensions, dimnames = NULL) {
  if (!.is_finite_numeric(dimensions) || length(dimensions) < 1L ||
      anyNA(dimensions) || any(dimensions < 0) || any(dimensions %% 1 != 0)) {
    .input_error("Array dimensions must be finite nonnegative whole scalars.")
  }
  max_exact <- 2^53
  elements <- prod(dimensions)
  if (!is.finite(elements) || elements > floor(max_exact / 8)) {
    .invariant_error("Array byte accounting overflows exact representation.")
  }
  skeleton_dimensions <- dimensions
  skeleton_dimensions[[1L]] <- 0
  skeleton <- array(numeric(), skeleton_dimensions, dimnames = dimnames)
  bytes <- as.double(utils::object.size(skeleton)) + 8 * elements
  if (!is.finite(bytes) || bytes > max_exact) {
    .invariant_error("Array byte accounting overflows exact representation.")
  }
  bytes
}

#' Construct a crossform-owned workspace plan
#'
#' The hard budget covers package-owned live objects and conservative temporary
#' overlap. Process baseline RSS and absolute RSS are observations, not hidden
#' additions to the workspace model.
#'
#' @param frame_bytes Resident dense or sparse frame storage.
#' @param resident_source_bytes Resident source objects owned by crossform.
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
    if (!.is_number(value) || value < lower || value %% 1 != 0 ||
        value > max_exact) {
      .input_error(sprintf(
        "`%s` must be one %sfinite whole scalar no greater than 2^53.",
        name, if (positive) "positive " else "nonnegative "))
    }
    value
  }
  for (name in names(categories)) whole(categories[[name]], name)
  workers <- whole(workers, "workers", positive = TRUE)
  n_active <- whole(n_active, "n_active", positive = TRUE)
  if (n_active > workers) .input_error("`n_active` cannot exceed `workers`.")
  if (!.is_number(safety_factor) || safety_factor < 1) {
    .input_error(
      "`safety_factor` must be one finite number greater than or equal to one."
    )
  }
  optional <- function(value, name) {
    if (!is.null(value) && (!is.numeric(value) || length(value) != 1L ||
        is.na(value) || !is.finite(value) || value < 0)) {
      .input_error(sprintf(
        "`%s` must be NULL or one nonnegative finite byte count.",
        name))
    }
    value
  }
  optional(budget_bytes, "budget_bytes")
  if (!is.null(budget_bytes) && budget_bytes == 0) {
    .input_error("`budget_bytes` must be NULL or positive.")
  }
  optional(measured_workspace_bytes, "measured_workspace_bytes")
  optional(baseline_rss_bytes, "baseline_rss_bytes")
  optional(peak_rss_bytes, "peak_rss_bytes")
  if (!is.null(baseline_rss_bytes) && !is.null(peak_rss_bytes) &&
      peak_rss_bytes < baseline_rss_bytes) {
    .input_error("Peak RSS cannot be smaller than baseline RSS.")
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
    .invariant_error(
      "Workspace byte accounting overflows exact representation."
    )
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
    prediction_kind = "crossform_owned_workspace_upper_bound"
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
  .check_flag(same_relation, "same_relation")
  if (!.is_flag(retain_first_moments) || !.is_flag(form_total)) {
    .input_error("First-moment and total planning flags must be TRUE or FALSE.")
  }
  if (!retain_first_moments && !form_total) {
    .input_error(
      "A plan must retain first moments, form total output, or both."
    )
  }
  .validate_effect_names(left_effects, length(left_effects))
  .validate_effect_names(right_effects, length(right_effects))
  if (!.is_strings(left_partitions, unique = TRUE) ||
      length(left_partitions) < 1L ||
      !.is_strings(right_partitions, unique = TRUE) ||
      length(right_partitions) < 1L) {
    .input_error(
      "Effect-form partition axes must be unique nonempty identifiers."
    )
  }
  if (same_relation && (!identical(left_partitions, right_partitions) ||
      !identical(left_effects, right_effects))) {
    .input_error(
      "A shared relation plan requires identical partition and effect axes."
    )
  }

  q_left <- length(left_effects)
  q_right <- length(right_effects)
  physical_width <- if (codec == "rectangular") {
    q_left * q_right
  } else {
    if (!same_relation || !identical(left_effects, right_effects)) {
      .input_error("Symmetric-packed planning requires a self-form.")
    }
    q_left * (q_left + 1L) / 2L
  }
  structured_query <- !is.null(query) && .is_pair_difference_query(query)
  if (!is.null(query) && !structured_query &&
      (!is.matrix(query) || !is.numeric(query) ||
      nrow(query) != physical_width || ncol(query) < 1L ||
      any(!is.finite(query)))) {
    .contract_error("`query` must match the finite physical form coordinates.")
  }

  features <- ncol(frame$weights)
  measurements <- nrow(frame$weights)
  output_width <- if (is.null(query)) {
    physical_width
  } else {
    .query_output_width(query)
  }
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
    if (!form_total || is.null(query)) 0 else if (structured_query) {
      # Structured evaluation tiles the pair axis. Account conservatively for
      # the accumulator plus indexing, difference, product, weighting, and
      # addition temporaries for one live tile, plus the query payload.
      6 * .dense_double_matrix_bytes(
        min(64L, length(query$pair_left)), f
      ) +
        .query_payload_bytes(query)
    } else {
      .dense_double_matrix_bytes(q_left, q_right) +
        .query_payload_bytes(query)
    }
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

# Conservative plans for the streaming geometry compiler ---------------------
#
# These are declarations too: they compute byte counts from shapes and
# policies and allocate nothing. The compiler selects among them before it
# opens a source.

.compiler_memory_plan <- function(x, at, compute, feature_block, row_tile,
                                  coordinate_tile, output_width, storage,
                                  requirements, query = NULL,
                                  right_relation = NULL) {
  q <- length(x$effects)
  h <- q * (q + 1L) / 2L
  p <- x$n_features
  m <- nrow(at$weights)
  r <- length(x$partitions)
  f <- min(feature_block, p)
  rows <- min(row_tile, m)
  coordinates <- min(coordinate_tile, output_width)
  max_observations <- max(vapply(x$sources, function(source) source$dim[[1L]],
    integer(1)))
  source_shared <- sum(vapply(x$sources, function(source) {
    if (identical(source$kind, "matrix")) {
      prod(as.double(source$dim)) * 8 +
        as.double(utils::object.size(source))
    } else {
      0
    }
  }, numeric(1)))
  distinct_handles <- unique(vapply(x$sources, function(source) {
    descriptor <- source$descriptor
    if (is.null(descriptor) || identical(descriptor$access, "coordinator")) {
      return(NA_character_)
    }
    .source_descriptor_key(descriptor)
  }, character(1)))
  distinct_handles <- sum(!is.na(distinct_handles))
  source_block <- max_observations * f * 8
  relation_block <- r * q * f * 8 +
    if (is.null(right_relation)) 0 else
      length(right_relation$partitions) *
        length(right_relation$effects) * f * 8
  atom_coordinates <- if (identical(
      requirements$materialization, "full_geometry")) {
    h
  } else {
    output_width
  }
  # Query-fused execution forms only the requested per-feature coordinates.
  # The live query representation is accounted explicitly: a structured
  # pair-difference query stores index vectors plus its coefficient map and
  # a pair-by-feature-block workspace; a dense physical query stores its
  # full matrix plus one reconstructed q-by-q operator at a time.
  query_operator <- if (requirements$total &&
      !identical(requirements$materialization, "full_geometry")) {
    if (!is.null(query) && .is_pair_difference_query(query)) {
      .query_payload_bytes(query) + 8 * length(query$pair_left) * f
    } else {
      .query_payload_bytes(query) + q * q * 8
    }
  } else {
    0
  }
  atom_block <- if (requirements$total) {
    f * atom_coordinates * 8 + query_operator
  } else {
    0
  }
  local_bytes <- if (requirements$coherent) m * q * r * 8 else 0
  marginal_bytes <- if (requirements$marginals) 2 * m * q * 8 else 0
  component_count <- as.integer(requirements$total) +
    as.integer(requirements$coherent)
  durable_geometry <- if (storage == "memory") {
    component_count * m * output_width * 8
  } else {
    0
  }
  output_bytes <- marginal_bytes + durable_geometry
  total_contraction <- if (requirements$total) (
    rows * f + f * coordinates + rows * coordinates + rows * q +
      rows * h
  ) * 8 else 0
  coherent_contraction <- if (requirements$coherent) (
    rows * f + rows * q * 3 + rows * h + rows * output_width
  ) * 8 else 0
  contraction <- max(total_contraction, coherent_contraction)
  replacement <- if (requirements$total) {
    (2 * rows * coordinates + 2 * rows * q) * 8
  } else {
    2 * rows * q * 8
  }

  memory_plan(
    frame_bytes = as.double(utils::object.size(at$weights)),
    resident_source_bytes = source_shared,
    source_handle_bytes = distinct_handles * 4096,
    source_block_bytes = source_block,
    relation_block_bytes = relation_block,
    atom_block_bytes = atom_block,
    local_state_bytes = local_bytes,
    output_bytes = output_bytes,
    contraction_bytes = contraction,
    replacement_copy_bytes = replacement,
    workers = 1L,
    n_active = 1L,
    budget_bytes = compute$workspace_bytes
  )
}

.descending_tiles <- function(maximum, initial = maximum) {
  value <- min(as.integer(maximum), as.integer(initial))
  out <- integer()
  while (value > 1L) {
    out <- c(out, value)
    value <- max(1L, as.integer(floor(value / 2)))
  }
  unique(c(out, 1L))
}

.select_compiler_memory_plan <- function(x, at, compute, output_width, storage,
                                         requirements, query = NULL,
                                         right_relation = NULL) {
  if (is.null(compute$workspace_bytes)) {
    feature <- min(if (is.null(compute$block_features)) 1024L else
      compute$block_features, x$n_features)
    rows <- min(256L, nrow(at$weights))
    coordinates <- min(64L, output_width)
    return(list(feature_block = feature, row_tile = rows,
      coordinate_tile = coordinates,
      memory = .compiler_memory_plan(x, at, compute, feature, rows,
        coordinates, output_width, storage, requirements, query = query,
        right_relation = right_relation)))
  }
  feature_candidates <- if (is.null(compute$block_features)) {
    .descending_tiles(x$n_features)
  } else {
    min(as.integer(compute$block_features), x$n_features)
  }
  row_candidates <- .descending_tiles(nrow(at$weights))
  coordinate_candidates <- .descending_tiles(output_width)
  smallest <- NULL
  for (feature in feature_candidates) {
    for (rows in row_candidates) {
      for (coordinates in coordinate_candidates) {
        memory <- .compiler_memory_plan(x, at, compute, feature, rows,
          coordinates, output_width, storage, requirements, query = query,
          right_relation = right_relation)
        smallest <- list(feature_block = feature, row_tile = rows,
          coordinate_tile = coordinates, memory = memory)
        if (isTRUE(memory$fits_budget)) return(smallest)
      }
    }
  }
  smallest
}

.support_metric_memory_plan <- function(x, at, metric_schedule, compute,
                                        output_width, storage, requirements) {
  schedule <- .validate_geometry_metric_schedule(metric_schedule)
  if (!identical(schedule$materialization, "fixed_metric")) {
    .input_error("Support-metric planning requires a fixed metric schedule.")
  }
  support_sizes <- if (!is.null(at$support_index)) {
    diff(at$support_index$ptr)
  } else if (inherits(at$weights, "Matrix")) {
    as.numeric(Matrix::rowSums(at$weights != 0))
  } else {
    rowSums(at$weights != 0)
  }
  k <- max(support_sizes)
  q <- length(x$effects)
  r <- length(x$partitions)
  m <- nrow(at$weights)
  components <- as.integer(requirements$total) +
    as.integer(requirements$coherent)
  output <- if (storage == "memory") 8 * components * m * output_width else 0
  first <- if (requirements$marginals) 8 * m * q * r else 0
  metric_bytes <- as.double(utils::object.size(schedule$metric$value))
  plan <- memory_plan(
    frame_bytes = as.double(utils::object.size(at$weights)),
    resident_source_bytes = metric_bytes,
    source_block_bytes = 8 * max(vapply(
      x$sources, function(source) source$dim[[1L]], integer(1)
    )) * k,
    relation_block_bytes = 8 * r * q * k,
    atom_block_bytes = 0,
    local_state_bytes = first,
    output_bytes = output,
    contraction_bytes = 8 * (2 * k * k + 3 * q * k + q * q),
    replacement_copy_bytes = 8 * max(output_width, q * q),
    workers = 1L,
    n_active = 1L,
    budget_bytes = compute$workspace_bytes
  )
  list(
    feature_block = as.integer(k),
    row_tile = 1L,
    coordinate_tile = as.integer(min(output_width, 64L)),
    memory = plan
  )
}
