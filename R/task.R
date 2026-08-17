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
        !.is_strings(partitions, unique = TRUE) || length(partitions) < 1L ||
        !setequal(names(relations), partitions)) {
      .input_error(
        sprintf("`%s_relations` must match one named partition family.", side)
      )
    }
    .validate_effect_names(effects, length(effects))
    expected <- c(length(effects), length(feature_ids))
    for (partition in partitions) {
      value <- relations[[partition]]
      if (!.is_finite_matrix(value) || !identical(dim(value), expected)) {
        .input_error(sprintf(
          "Every %s relation input must be a finite effect-by-feature matrix.",
          side
        ))
      }
    }
    invisible(TRUE)
  }
  if (!.is_finite_numeric(feature_ids) || length(feature_ids) < 1L ||
      anyNA(feature_ids) || any(feature_ids < 1) || any(feature_ids %% 1 != 0) ||
      is.unsorted(feature_ids, strictly = TRUE)) {
    .input_error("`feature_ids` must be strictly increasing positive integers.")
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
    .input_error(
      "Symmetric-packed atoms require a self-adjoint self-form task."
    )
  }
  logical_width <- length(left_effects) * length(right_effects)
  # Both codecs report `$physical_width` as an integer. `q * (q + 1L) / 2L` is
  # exact but double-typed in R, so it is coerced here rather than leaving one
  # codec reporting `6` and the other `6L` for the same quantity. This field
  # reaches no signature: it is compared against `nrow(query)`, used as an
  # output width, and returned. The same expression is still computed as a
  # double for local widths in R/memory-plan.R and R/kernel.R, which are not
  # reported anywhere; R/geometry-plan.R already coerces its `packed_width`.
  physical_width <- if (codec == "rectangular") {
    logical_width
  } else {
    as.integer(length(left_effects) * (length(left_effects) + 1L) / 2L)
  }
  .validate_task_query(
    query, physical_width, left_effects, right_effects, same_relation
  )
  list(
    feature_ids = feature_ids,
    codec = codec,
    logical_width = logical_width,
    physical_width = physical_width
  )
}

.packed_effect_form_atoms_oracle <- function(left_relations, right_relations,
                                             left_index, right_index,
                                             edge_weight, codec) {
  q_left <- nrow(left_relations[[1L]])
  q_right <- nrow(right_relations[[1L]])
  n_features <- ncol(left_relations[[1L]])
  packed <- identical(codec, "symmetric_packed")
  n_coords <- if (packed) {
    as.integer(q_left * (q_left + 1L) / 2L)
  } else {
    q_left * q_right
  }
  atoms <- matrix(0, n_features, n_coords)
  coordinate <- 0L
  for (column in seq_len(q_right)) {
    rows <- if (packed) column:q_left else seq_len(q_left)
    for (row in rows) {
      coordinate <- coordinate + 1L
      work <- numeric(n_features)
      for (edge in seq_along(edge_weight)) {
        work <- work + edge_weight[[edge]] *
          left_relations[[left_index[[edge]]]][row, ] *
          right_relations[[right_index[[edge]]]][column, ]
      }
      if (packed && row != column) {
        work <- sqrt(2) * work
      }
      atoms[, coordinate] <- work
    }
  }
  atoms
}

.packed_effect_form_atoms <- function(left_relations, right_relations,
                                      left_index, right_index,
                                      edge_weight, codec) {
  .packed_effect_form_atoms_cpp(
    unname(left_relations),
    unname(right_relations),
    as.integer(left_index),
    as.integer(right_index),
    as.numeric(edge_weight),
    identical(codec, "symmetric_packed")
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
  .check_flag(form_atoms, "form_atoms")
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
  atoms <- if (!form_atoms || structured_query) {
    NULL
  } else {
    matrix(0, length(feature_ids), output_width)
  }
  feature_count <- length(feature_ids)
  max_atom_work_bytes <- if (!form_atoms) {
    0
  } else if (is.null(query)) {
    # Native packed/rectangular accumulation writes directly to atoms.
    0
  } else if (structured_query) {
    # Native fused evaluation writes directly to atoms. Difference and
    # product tiles are never allocated, and there is no per-tile gc().
    0
  } else {
    8 * (feature_count + q_left * q_right)
  }

  if (form_atoms && structured_query) {
    atoms <- .fused_pair_difference_atoms(
      left_relations, right_relations, left_index, right_index,
      ordered_edges$weight, query$pair_left, query$pair_right,
      query$coefficients
    )
    if (any(!is.finite(atoms))) {
      .input_error(
        paste0("Direct effect-form querying produced non-finite values.",
        " Finite inputs overflowed double precision during the computation; rescale the responses (for example to unit variance) before building the relation.")
      )
    }
  } else if (form_atoms && is.null(query)) {
    atoms <- .packed_effect_form_atoms(
      left_relations, right_relations, left_index, right_index,
      ordered_edges$weight, validated$codec
    )
    if (any(!is.finite(atoms))) {
      .input_error(
        paste0("Effect-form atom formation produced non-finite values.",
      " Finite inputs overflowed double precision during the computation; rescale the responses (for example to unit variance) before building the relation.")
      )
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
        .input_error(
          paste0("Direct effect-form querying produced non-finite values.",
          " Finite inputs overflowed double precision during the computation; rescale the responses (for example to unit variance) before building the relation.")
        )
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
