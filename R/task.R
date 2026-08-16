# Canonical feature tasks ---------------------------------------------------

.validate_effect_form_task_inputs <- function(left_relations, right_relations,
                                              feature_ids, left_effects,
                                              right_effects, left_partitions,
                                              right_partitions, ordered_edges,
                                              codec, query, same_relation) {
  validate_family <- function(relations, partitions, effects, side) {
    if (!is.list(relations) || length(relations) < 1L ||
        is.null(names(relations)) || anyNA(names(relations)) ||
        any(!nzchar(names(relations))) || anyDuplicated(names(relations)) ||
        !is.character(partitions) || length(partitions) < 1L ||
        anyNA(partitions) || any(!nzchar(partitions)) ||
        anyDuplicated(partitions) || !setequal(names(relations), partitions)) {
      stop(sprintf("`%s_relations` must match one named partition family.", side),
        call. = FALSE)
    }
    .validate_effect_names(effects, length(effects))
    expected <- c(length(effects), length(feature_ids))
    for (partition in partitions) {
      value <- relations[[partition]]
      if (!is.matrix(value) || !is.numeric(value) ||
          !identical(dim(value), expected) || any(!is.finite(value))) {
        stop(sprintf(
          "Every %s relation input must be a finite effect-by-feature matrix.",
          side
        ), call. = FALSE)
      }
    }
    invisible(TRUE)
  }
  if (!is.numeric(feature_ids) || length(feature_ids) < 1L ||
      anyNA(feature_ids) || any(!is.finite(feature_ids)) ||
      any(feature_ids < 1) || any(feature_ids %% 1 != 0) ||
      is.unsorted(feature_ids, strictly = TRUE)) {
    stop("`feature_ids` must be strictly increasing positive integers.",
      call. = FALSE)
  }
  feature_ids <- as.integer(feature_ids)
  validate_family(left_relations, left_partitions, left_effects, "left")
  validate_family(right_relations, right_partitions, right_effects, "right")
  .validate_ordered_partition_edges(
    ordered_edges, left_partitions, right_partitions, same_relation
  )
  codec <- match.arg(codec, c("rectangular", "symmetric_packed"))
  if (codec == "symmetric_packed" &&
      (!same_relation || !identical(left_effects, right_effects) ||
       !identical(attr(ordered_edges, "expansion"),
         "self_adjoint_half_edges"))) {
    stop("Symmetric-packed atoms require a self-adjoint self-form task.",
      call. = FALSE)
  }
  logical_width <- length(left_effects) * length(right_effects)
  # Both codecs report `$physical_width` as an integer. `q * (q + 1L) / 2L` is
  # exact but double-typed in R, so it is coerced here rather than leaving one
  # codec reporting `6` and the other `6L` for the same quantity. This field
  # reaches no signature: it is compared against `nrow(query)`, used as an
  # output width, and returned. The same expression is still computed as a
  # double for local widths in R/memory.R and R/kernel.R, which are not
  # reported anywhere; R/geometry-plan.R already coerces its `packed_width`.
  physical_width <- if (codec == "rectangular") {
    logical_width
  } else {
    as.integer(length(left_effects) * (length(left_effects) + 1L) / 2L)
  }
  if (!is.null(query)) {
    if (.is_pair_difference_query(query)) {
      .validate_pair_difference_for_task(
        query, left_effects, right_effects, same_relation
      )
    } else if (!is.matrix(query) || !is.numeric(query) ||
        nrow(query) != physical_width || ncol(query) < 1L ||
        any(!is.finite(query))) {
      stop("`query` must match the finite physical form coordinates.",
        call. = FALSE)
    }
  }
  list(
    feature_ids = feature_ids,
    codec = codec,
    logical_width = logical_width,
    physical_width = physical_width
  )
}

# The universal numerical primitive. Every edge forms one ordered outer
# product; symmetry is supplied only by an explicit reverse half-edge.
.effect_form_feature_task <- function(left_relations, right_relations,
                                      feature_ids, left_effects, right_effects,
                                      left_partitions = names(left_relations),
                                      right_partitions = names(right_relations),
                                      ordered_edges,
                                      codec = c("rectangular", "symmetric_packed"),
                                      query = NULL, form_atoms = TRUE,
                                      same_relation = FALSE) {
  validated <- .validate_effect_form_task_inputs(
    left_relations, right_relations, feature_ids, left_effects, right_effects,
    left_partitions, right_partitions, ordered_edges, codec, query,
    same_relation
  )
  if (!is.logical(form_atoms) || length(form_atoms) != 1L || is.na(form_atoms)) {
    stop("`form_atoms` must be TRUE or FALSE.", call. = FALSE)
  }
  feature_ids <- validated$feature_ids
  left_relations <- left_relations[left_partitions]
  right_relations <- right_relations[right_partitions]
  left_index <- match(ordered_edges$left, left_partitions)
  right_index <- match(ordered_edges$right, right_partitions)
  q_left <- length(left_effects)
  q_right <- length(right_effects)
  structured_query <- !is.null(query) && .is_pair_difference_query(query)
  output_width <- if (is.null(query)) {
    validated$physical_width
  } else {
    .query_output_width(query)
  }
  tiled_pair_output <- form_atoms && structured_query &&
    is.null(query$coefficients)
  atoms <- if (!form_atoms || tiled_pair_output) {
    NULL
  } else {
    matrix(0, length(feature_ids), output_width)
  }
  feature_count <- length(feature_ids)
  pair_tile_size <- if (structured_query) {
    min(64L, length(query$pair_left))
  } else {
    0L
  }
  max_atom_work_bytes <- if (!form_atoms) {
    0
  } else if (is.null(query)) {
    8 * feature_count
  } else if (structured_query) {
    # One accumulator and the difference/product temporaries for one bounded
    # pair tile are live; the full pair-by-feature matrix is never formed.
    8 * (6 * pair_tile_size * feature_count)
  } else {
    8 * (feature_count + q_left * q_right)
  }

  if (form_atoms && structured_query) {
    atom_tiles <- if (tiled_pair_output) {
      vector("list", ceiling(length(query$pair_left) / pair_tile_size))
    } else {
      NULL
    }
    tile_index <- 0L
    for (pair_start in seq.int(1L, length(query$pair_left),
        by = pair_tile_size)) {
      tile_index <- tile_index + 1L
      pair_indices <- pair_start:min(
        pair_start + pair_tile_size - 1L, length(query$pair_left)
      )
      pair_work <- matrix(0, length(pair_indices), feature_count)
      for (edge in seq_len(nrow(ordered_edges))) {
        edge_product <- .pair_difference_edge_products(
          query,
          left_relations[[left_index[[edge]]]],
          right_relations[[right_index[[edge]]]],
          pair_indices
        )
        pair_work <- pair_work + ordered_edges$weight[[edge]] * edge_product
      }
      if (any(!is.finite(pair_work))) {
        stop("Direct effect-form querying produced non-finite values.",
          call. = FALSE)
      }
      if (is.null(query$coefficients)) {
        atom_tiles[[tile_index]] <- t(pair_work)
      } else {
        atoms <- atoms + t(
          query$coefficients[, pair_indices, drop = FALSE] %*% pair_work
        )
      }
      rm(pair_work, edge_product)
      invisible(gc(full = FALSE))
    }
    if (tiled_pair_output) {
      atoms <- do.call(cbind, atom_tiles)
    }
  } else if (form_atoms && is.null(query)) {
    coordinate <- 0L
    for (column in seq_len(q_right)) {
      rows <- if (identical(validated$codec, "symmetric_packed")) {
        column:q_left
      } else {
        seq_len(q_left)
      }
      for (row in rows) {
        coordinate <- coordinate + 1L
        work <- numeric(length(feature_ids))
        for (edge in seq_len(nrow(ordered_edges))) {
          left <- left_relations[[left_index[[edge]]]]
          right <- right_relations[[right_index[[edge]]]]
          work <- work + ordered_edges$weight[[edge]] *
            left[row, ] * right[column, ]
        }
        if (identical(validated$codec, "symmetric_packed") && row != column) {
          work <- sqrt(2) * work
        }
        if (any(!is.finite(work))) {
          stop("Effect-form atom formation produced non-finite values.",
            call. = FALSE)
        }
        atoms[, coordinate] <- work
      }
    }
  } else if (form_atoms) {
    for (view in seq_len(ncol(query))) {
      work <- numeric(length(feature_ids))
      operator <- .physical_query_operator(
        query[, view], q_left, q_right, validated$codec
      )
      for (edge in seq_len(nrow(ordered_edges))) {
        left <- left_relations[[left_index[[edge]]]]
        right <- right_relations[[right_index[[edge]]]]
        work <- work + ordered_edges$weight[[edge]] *
          colSums(left * (operator %*% right))
      }
      if (any(!is.finite(work))) {
        stop("Direct effect-form querying produced non-finite values.",
          call. = FALSE)
      }
      atoms[, view] <- work
    }
  }

  relation_bytes <- 8 * feature_count * (
    length(left_partitions) * q_left +
      if (same_relation) 0 else length(right_partitions) * q_right
  )
  atom_bytes <- if (is.null(atoms)) 0 else 8 * length(atoms)
  structure(list(
    feature_ids = feature_ids,
    left_partitions = left_partitions,
    right_partitions = right_partitions,
    left_effects = left_effects,
    right_effects = right_effects,
    left_relations = left_relations,
    right_relations = right_relations,
    ordered_edges = ordered_edges,
    atoms = atoms,
    projected = !is.null(query) && form_atoms,
    atoms_formed = form_atoms,
    codec = validated$codec,
    logical_shape = as.integer(c(q_left, q_right)),
    physical_width = validated$physical_width,
    diagnostics = list(
      relation_bytes = relation_bytes,
      atom_bytes = if (is.null(query)) atom_bytes else 0,
      query_atom_bytes = if (is.null(query)) 0 else atom_bytes,
      max_atom_work_bytes = max_atom_work_bytes
    )
  ), class = "effect_form_feature_task_result")
}
