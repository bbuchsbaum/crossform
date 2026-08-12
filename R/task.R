# Canonical feature tasks ---------------------------------------------------

.measured_object_bytes <- function(x) as.double(utils::object.size(x))

.validate_relation_task_inputs <- function(relations, feature_ids, effects,
                                           partitions, over, query) {
  if (!is.list(relations) || length(relations) < 1L || is.null(names(relations)) ||
      anyNA(names(relations)) || any(!nzchar(names(relations))) ||
      anyDuplicated(names(relations))) {
    stop("`relations` must be a uniquely named list of relation matrices.",
      call. = FALSE)
  }
  if (!is.character(partitions) || length(partitions) < 1L || anyNA(partitions) ||
      any(!nzchar(partitions)) || anyDuplicated(partitions)) {
    stop("`partitions` must contain unique nonempty identifiers.", call. = FALSE)
  }
  if (!setequal(names(relations), partitions)) {
    stop("Relation names must equal the declared partition identifiers.",
      call. = FALSE)
  }
  if (!is.numeric(feature_ids) || length(feature_ids) < 1L ||
      anyNA(feature_ids) || any(!is.finite(feature_ids)) ||
      any(feature_ids < 1) || any(feature_ids %% 1 != 0) ||
      is.unsorted(feature_ids, strictly = TRUE)) {
    stop("`feature_ids` must be strictly increasing positive integers.",
      call. = FALSE)
  }
  feature_ids <- as.integer(feature_ids)
  if (!is.character(effects) || length(effects) < 1L || anyNA(effects) ||
      any(!nzchar(effects)) || anyDuplicated(effects)) {
    stop("`effects` must contain unique nonempty identifiers.", call. = FALSE)
  }
  q <- length(effects)
  expected_dim <- c(q, length(feature_ids))
  for (partition in partitions) {
    value <- relations[[partition]]
    if (!is.matrix(value) || !is.numeric(value) ||
        !identical(dim(value), expected_dim) || any(!is.finite(value))) {
      stop("Every relation task input must be a finite effect-by-feature matrix.",
        call. = FALSE)
    }
  }
  .validate_pairing(over)
  if (anyNA(match(over$left, partitions)) ||
      anyNA(match(over$right, partitions))) {
    stop("Every pairing endpoint must identify a declared partition.",
      call. = FALSE)
  }
  packed_width <- q * (q + 1L) / 2L
  if (!is.null(query) &&
      (!is.matrix(query) || !is.numeric(query) ||
       nrow(query) != packed_width || ncol(query) < 1L ||
       any(!is.finite(query)))) {
    stop("`query` must be a finite packed-coordinate-by-view matrix.",
      call. = FALSE)
  }
  list(feature_ids = feature_ids, packed_width = packed_width)
}

# A feature task is deliberately unaware of sources, frames, accumulators, and
# executors. The same pure function is used in-process and by any later worker.
.crossgram_feature_task <- function(relations, feature_ids, effects,
                                    partitions = names(relations), over,
                                    query = NULL, form_atoms = TRUE) {
  validated <- .validate_relation_task_inputs(
    relations, feature_ids, effects, partitions, over, query
  )
  feature_ids <- validated$feature_ids
  relations <- relations[partitions]
  left_index <- match(over$left, partitions)
  right_index <- match(over$right, partitions)
  q <- length(effects)

  if (!is.logical(form_atoms) || length(form_atoms) != 1L || is.na(form_atoms)) {
    stop("`form_atoms` must be TRUE or FALSE.", call. = FALSE)
  }

  atoms <- if (form_atoms) {
    matrix(0, length(feature_ids), validated$packed_width)
  } else {
    NULL
  }
  packed_coordinate <- 0L
  max_atom_work_bytes <- 0
  if (form_atoms) {
    for (column in seq_len(q)) {
      for (row in column:q) {
        packed_coordinate <- packed_coordinate + 1L
        atom_work <- numeric(length(feature_ids))
        for (edge in seq_len(nrow(over))) {
          left <- relations[[left_index[[edge]]]]
          right <- relations[[right_index[[edge]]]]
          atom_work <- atom_work + over$weight[[edge]] * 0.5 *
            (left[row, ] * right[column, ] +
              right[row, ] * left[column, ])
        }
        if (row != column) atom_work <- sqrt(2) * atom_work
        if (any(!is.finite(atom_work))) {
          stop("Cross-Gram atom formation produced non-finite values.",
            call. = FALSE)
        }
        atoms[, packed_coordinate] <- atom_work
        max_atom_work_bytes <- max(
          max_atom_work_bytes, .measured_object_bytes(atom_work)
        )
      }
    }
  }

  relation_bytes <- sum(vapply(
    relations, .measured_object_bytes, numeric(1)
  ))
  atom_bytes <- if (is.null(atoms)) 0 else .measured_object_bytes(atoms)
  query_atom_bytes <- 0
  if (!is.null(query) && form_atoms) {
    atoms <- atoms %*% query
    query_atom_bytes <- .measured_object_bytes(atoms)
  }

  structure(list(
    feature_ids = feature_ids,
    partitions = partitions,
    effects = effects,
    relations = relations,
    atoms = atoms,
    projected = !is.null(query) && form_atoms,
    atoms_formed = form_atoms,
    packed_width = validated$packed_width,
    diagnostics = list(
      relation_bytes = relation_bytes,
      atom_bytes = atom_bytes,
      query_atom_bytes = query_atom_bytes,
      max_atom_work_bytes = max_atom_work_bytes
    )
  ), class = "effect_feature_task_result")
}

.new_crossgram_reducer <- function(frame, partitions, effects, output_width,
                                   row_tile, coordinate_tile,
                                   accumulate_tile = NULL,
                                   retain_local_relations = FALSE,
                                   form_total = TRUE) {
  .validate_frame_for_compile(frame)
  row_tile <- .validate_tile_size(row_tile, "row_tile")
  coordinate_tile <- .validate_tile_size(coordinate_tile, "coordinate_tile")
  if (!is.logical(form_total) || length(form_total) != 1L || is.na(form_total)) {
    stop("`form_total` must be TRUE or FALSE.", call. = FALSE)
  }
  if (form_total) {
    output_width <- .validate_tile_size(output_width, "output_width")
  } else {
    output_width <- 0L
    if (!retain_local_relations) {
      stop("A reducer must form total output, local relations, or both.",
        call. = FALSE)
    }
  }
  if (!is.null(accumulate_tile) && !is.function(accumulate_tile)) {
    stop("`accumulate_tile` must be NULL or a function.", call. = FALSE)
  }
  if (!is.logical(retain_local_relations) ||
      length(retain_local_relations) != 1L || is.na(retain_local_relations)) {
    stop("`retain_local_relations` must be TRUE or FALSE.", call. = FALSE)
  }

  measurements <- nrow(frame$weights)
  q <- length(effects)
  state <- new.env(parent = emptyenv())
  state$weights <- frame$weights
  state$partitions <- partitions
  state$effects <- effects
  state$output_width <- output_width
  state$row_tile <- row_tile
  state$coordinate_tile <- coordinate_tile
  state$accumulate_tile <- accumulate_tile
  state$retain_local_relations <- retain_local_relations
  state$form_total <- form_total
  state$next_feature <- 1L
  state$output <- if (form_total && is.null(accumulate_tile)) {
    matrix(0, measurements, output_width)
  } else {
    NULL
  }
  state$local_relations <- if (retain_local_relations) {
    array(0, dim = c(measurements, q, length(partitions)),
      dimnames = list(NULL, effects, partitions))
  } else {
    NULL
  }
  state$diagnostics <- list(
    feature_blocks = 0L,
    relation_reads = 0L,
    atom_count = 0L,
    contraction_tiles = 0L,
    max_relation_bytes = 0,
    max_atom_bytes = 0,
    max_query_atom_bytes = 0,
    max_atom_work_bytes = 0,
    max_weight_slice_bytes = 0,
    max_atom_slice_bytes = 0,
    max_product_bytes = 0,
    max_existing_slice_bytes = 0,
    max_replacement_bytes = 0,
    max_local_product_bytes = 0,
    max_local_existing_bytes = 0,
    max_local_replacement_bytes = 0,
    max_live_temporary_bytes = 0,
    durable_output_bytes = if (is.null(state$output)) 0 else
      .measured_object_bytes(state$output),
    durable_local_relation_bytes = if (is.null(state$local_relations)) 0 else
      .measured_object_bytes(state$local_relations),
    measurement_kind = "named-live-R-object-size",
    external_accumulator_memory_measured = is.null(accumulate_tile)
  )
  class(state) <- "effect_crossgram_reducer"
  state
}

.reduce_crossgram_task <- function(reducer, task) {
  if (!inherits(reducer, "effect_crossgram_reducer")) {
    stop("`reducer` must be a cross-Gram reducer.", call. = FALSE)
  }
  if (!inherits(task, "effect_feature_task_result")) {
    stop("`task` must be a canonical feature-task result.", call. = FALSE)
  }
  if (!identical(task$partitions, reducer$partitions) ||
      !identical(task$effects, reducer$effects)) {
    stop("Task partition and effect identities must match the reducer.",
      call. = FALSE)
  }
  if (task$feature_ids[[1L]] != reducer$next_feature ||
      !identical(task$feature_ids,
        seq.int(reducer$next_feature,
          length.out = length(task$feature_ids)))) {
    stop("Feature tasks must be reduced once in canonical contiguous order.",
      call. = FALSE)
  }
  if (max(task$feature_ids) > ncol(reducer$weights) ||
      !identical(isTRUE(task$atoms_formed), reducer$form_total) ||
      (reducer$form_total && ncol(task$atoms) != reducer$output_width)) {
    stop("Task dimensions do not match the reducer.", call. = FALSE)
  }

  diagnostics <- reducer$diagnostics
  task_diag <- task$diagnostics
  diagnostics$feature_blocks <- diagnostics$feature_blocks + 1L
  diagnostics$relation_reads <- diagnostics$relation_reads +
    length(task$partitions)
  if (reducer$form_total) {
    diagnostics$atom_count <- diagnostics$atom_count + length(task$feature_ids)
  }
  diagnostics$max_relation_bytes <- max(
    diagnostics$max_relation_bytes, task_diag$relation_bytes
  )
  diagnostics$max_atom_bytes <- max(
    diagnostics$max_atom_bytes, task_diag$atom_bytes
  )
  diagnostics$max_query_atom_bytes <- max(
    diagnostics$max_query_atom_bytes, task_diag$query_atom_bytes
  )
  diagnostics$max_atom_work_bytes <- max(
    diagnostics$max_atom_work_bytes, task_diag$max_atom_work_bytes
  )
  base_live_bytes <- task_diag$relation_bytes + task_diag$atom_bytes +
    task_diag$query_atom_bytes + task_diag$max_atom_work_bytes

  measurements <- nrow(reducer$weights)
  for (row_start in .tile_starts(measurements, reducer$row_tile)) {
    rows <- row_start:min(row_start + reducer$row_tile - 1L, measurements)
    weight_slice <- reducer$weights[rows, task$feature_ids, drop = FALSE]
    weight_slice_bytes <- .measured_object_bytes(weight_slice)
    if (reducer$retain_local_relations) {
      for (partition_index in seq_along(task$partitions)) {
        local_product <- as.matrix(
          weight_slice %*% t(task$relations[[partition_index]])
        )
        local_existing <- matrix(
          reducer$local_relations[rows, , partition_index, drop = FALSE],
          nrow = length(rows), ncol = length(task$effects)
        )
        local_replacement <- local_existing + local_product
        reducer$local_relations[rows, , partition_index] <- local_replacement
        local_product_bytes <- .measured_object_bytes(local_product)
        local_existing_bytes <- .measured_object_bytes(local_existing)
        local_replacement_bytes <- .measured_object_bytes(local_replacement)
        diagnostics$max_local_product_bytes <- max(
          diagnostics$max_local_product_bytes, local_product_bytes
        )
        diagnostics$max_local_existing_bytes <- max(
          diagnostics$max_local_existing_bytes, local_existing_bytes
        )
        diagnostics$max_local_replacement_bytes <- max(
          diagnostics$max_local_replacement_bytes, local_replacement_bytes
        )
        diagnostics$max_live_temporary_bytes <- max(
          diagnostics$max_live_temporary_bytes,
          base_live_bytes + weight_slice_bytes + local_product_bytes +
            local_existing_bytes + local_replacement_bytes
        )
      }
    }

    if (reducer$form_total) for (coordinate_start in .tile_starts(
      reducer$output_width, reducer$coordinate_tile
    )) {
      coordinates <- coordinate_start:min(
        coordinate_start + reducer$coordinate_tile - 1L,
        reducer$output_width
      )
      atom_slice <- task$atoms[, coordinates, drop = FALSE]
      product <- as.matrix(weight_slice %*% atom_slice)
      if (any(!is.finite(product))) {
        stop("Sparse frame contraction produced non-finite values.",
          call. = FALSE)
      }
      existing <- replacement <- NULL
      if (is.null(reducer$accumulate_tile)) {
        existing <- reducer$output[rows, coordinates, drop = FALSE]
        replacement <- existing + product
        reducer$output[rows, coordinates] <- replacement
      } else {
        reducer$accumulate_tile(rows, coordinates, product)
      }

      atom_slice_bytes <- .measured_object_bytes(atom_slice)
      product_bytes <- .measured_object_bytes(product)
      existing_bytes <- if (is.null(existing)) 0 else
        .measured_object_bytes(existing)
      replacement_bytes <- if (is.null(replacement)) 0 else
        .measured_object_bytes(replacement)
      diagnostics$max_weight_slice_bytes <- max(
        diagnostics$max_weight_slice_bytes, weight_slice_bytes
      )
      diagnostics$max_atom_slice_bytes <- max(
        diagnostics$max_atom_slice_bytes, atom_slice_bytes
      )
      diagnostics$max_product_bytes <- max(
        diagnostics$max_product_bytes, product_bytes
      )
      diagnostics$max_existing_slice_bytes <- max(
        diagnostics$max_existing_slice_bytes, existing_bytes
      )
      diagnostics$max_replacement_bytes <- max(
        diagnostics$max_replacement_bytes, replacement_bytes
      )
      diagnostics$max_live_temporary_bytes <- max(
        diagnostics$max_live_temporary_bytes,
        base_live_bytes + weight_slice_bytes + atom_slice_bytes +
          product_bytes + existing_bytes + replacement_bytes
      )
      diagnostics$contraction_tiles <- diagnostics$contraction_tiles + 1L
    }
  }

  reducer$next_feature <- max(task$feature_ids) + 1L
  reducer$diagnostics <- diagnostics
  invisible(reducer)
}

.reduce_crossgram_tasks <- function(tasks, frame, partitions, effects,
                                    output_width, row_tile = 1024L,
                                    coordinate_tile = 256L,
                                    accumulate_tile = NULL,
                                    retain_local_relations = FALSE,
                                    form_total = TRUE) {
  if (!is.list(tasks) || length(tasks) < 1L ||
      any(!vapply(tasks, inherits, logical(1), "effect_feature_task_result"))) {
    stop("`tasks` must be a nonempty list of canonical feature-task results.",
      call. = FALSE)
  }
  starts <- vapply(tasks, function(task) task$feature_ids[[1L]], integer(1))
  if (anyDuplicated(starts)) {
    stop("Feature tasks must have unique canonical starting features.",
      call. = FALSE)
  }
  tasks <- tasks[order(starts)]
  reducer <- .new_crossgram_reducer(
    frame, partitions, effects, output_width, row_tile, coordinate_tile,
    accumulate_tile, retain_local_relations, form_total
  )
  for (task in tasks) .reduce_crossgram_task(reducer, task)
  if (reducer$next_feature != ncol(frame$weights) + 1L) {
    stop("Feature tasks must cover the frame exactly once without gaps.",
      call. = FALSE)
  }
  list(
    value = reducer$output,
    local_relations = reducer$local_relations,
    diagnostics = reducer$diagnostics
  )
}
