# Bounded additive contraction ----------------------------------------------

.tiled_contraction <- function(weights, atoms, row_tile, coordinate_tile,
                               feature_tile, write_tile = NULL) {
  if (!is.matrix(weights) || !is.numeric(weights) || any(!is.finite(weights))) {
    stop("`weights` must be a finite numeric matrix.", call. = FALSE)
  }
  if (!is.matrix(atoms) || !is.numeric(atoms) || any(!is.finite(atoms))) {
    stop("`atoms` must be a finite numeric matrix.", call. = FALSE)
  }
  if (ncol(weights) != nrow(atoms)) {
    stop("The feature dimension of `weights` and `atoms` must agree.",
      call. = FALSE)
  }
  if (nrow(weights) < 1L || ncol(weights) < 1L || ncol(atoms) < 1L) {
    stop("Contraction inputs must have positive dimensions.", call. = FALSE)
  }
  row_tile <- .validate_tile_size(row_tile, "row_tile")
  coordinate_tile <- .validate_tile_size(coordinate_tile, "coordinate_tile")
  feature_tile <- .validate_tile_size(feature_tile, "feature_tile")
  if (!is.null(write_tile) && !is.function(write_tile)) {
    stop("`write_tile` must be NULL or a function.", call. = FALSE)
  }

  measurements <- nrow(weights)
  features <- ncol(weights)
  coordinates <- ncol(atoms)
  output <- if (is.null(write_tile)) matrix(0, measurements, coordinates) else NULL
  tile_count <- 0L
  max_temporary_elements <- 0L

  for (row_start in .tile_starts(measurements, row_tile)) {
    rows <- row_start:min(row_start + row_tile - 1L, measurements)
    for (coordinate_start in .tile_starts(coordinates, coordinate_tile)) {
      coordinates_in_tile <- coordinate_start:min(
        coordinate_start + coordinate_tile - 1L, coordinates
      )
      tile <- matrix(0, length(rows), length(coordinates_in_tile))
      max_temporary_elements <- max(max_temporary_elements, length(tile))

      for (feature_start in .tile_starts(features, feature_tile)) {
        features_in_tile <- feature_start:min(
          feature_start + feature_tile - 1L, features
        )
        tile <- tile +
          weights[rows, features_in_tile, drop = FALSE] %*%
          atoms[features_in_tile, coordinates_in_tile, drop = FALSE]
      }

      if (is.null(write_tile)) {
        output[rows, coordinates_in_tile] <- tile
      } else {
        write_tile(rows, coordinates_in_tile, tile)
      }
      tile_count <- tile_count + 1L
    }
  }

  list(
    value = output,
    diagnostics = list(
      row_tile = row_tile,
      coordinate_tile = coordinate_tile,
      feature_tile = feature_tile,
      tile_count = tile_count,
      max_temporary_elements = max_temporary_elements
    )
  )
}

# First moments, flat and named ----------------------------------------------
#
# A partition family of effect-by-feature blocks contracts against the frame
# one partition at a time only if the code asks it to. Stacking the
# transposed blocks side by side gives a feature-by-(effect within partition)
# panel, so the whole family contracts in one product against a measurement
# tile of the frame. The stacking order matches the column-major storage of
# the measurement-by-effect-by-partition array the kernel returns, so the
# accumulator can be a flat matrix and become that array by relabelling.

.first_moment_feature_panel <- function(relations) {
  if (length(relations) == 1L) return(t(relations[[1L]]))
  do.call(cbind, lapply(relations, t))
}

.named_first_moment_array <- function(flat, measurements, effects,
                                      partitions) {
  array(flat, dim = c(measurements, length(effects), length(partitions)),
    dimnames = list(NULL, effects, partitions))
}

# Stream the universal total effect form. The feature task owns ordered
# products and optional direct querying; this coordinator owns only bounded
# reads and spatial contraction.
.streamed_effect_form_contraction <- function(
    frame, read_left, read_right = read_left,
    left_partitions, right_partitions, left_effects, right_effects,
    ordered_edges, codec = c("rectangular", "symmetric_packed"),
    same_relation = FALSE, query = NULL,
    feature_block = 1024L, row_tile = 1024L, coordinate_tile = 256L,
    accumulate_tile = NULL, retain_first_moments = FALSE,
    form_total = TRUE, task_observer = NULL) {
  .validate_frame_for_compile(frame)
  if (!identical(frame$representation, "additive_diagonal")) {
    stop("The streamed effect-form lowering requires an additive diagonal frame.",
      call. = FALSE)
  }
  if (!is.function(read_left) || !is.function(read_right)) {
    stop("Effect-form relation readers must be functions.", call. = FALSE)
  }
  if (!is.logical(same_relation) || length(same_relation) != 1L ||
      is.na(same_relation)) {
    stop("`same_relation` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(retain_first_moments) ||
      length(retain_first_moments) != 1L || is.na(retain_first_moments)) {
    stop("`retain_first_moments` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(form_total) || length(form_total) != 1L ||
      is.na(form_total)) {
    stop("`form_total` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!form_total && !retain_first_moments) {
    stop("Streaming must form total output, first moments, or both.",
      call. = FALSE)
  }
  if (same_relation && (!identical(left_partitions, right_partitions) ||
      !identical(left_effects, right_effects))) {
    stop("A shared relation reader requires identical partition and effect axes.",
      call. = FALSE)
  }
  .validate_effect_names(left_effects, length(left_effects))
  .validate_effect_names(right_effects, length(right_effects))
  .validate_ordered_partition_edges(
    ordered_edges, left_partitions, right_partitions, same_relation
  )
  codec <- match.arg(codec)
  q_left <- length(left_effects)
  q_right <- length(right_effects)
  physical_width <- if (codec == "rectangular") {
    q_left * q_right
  } else {
    if (!same_relation || !identical(left_effects, right_effects) ||
        !identical(attr(ordered_edges, "expansion"),
          "self_adjoint_half_edges")) {
      stop("Symmetric-packed streaming requires a self-adjoint self form.",
        call. = FALSE)
    }
    q_left * (q_left + 1L) / 2L
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
  feature_block <- .validate_tile_size(feature_block, "feature_block")
  row_tile <- .validate_tile_size(row_tile, "row_tile")
  coordinate_tile <- .validate_tile_size(coordinate_tile, "coordinate_tile")
  if (!is.null(accumulate_tile) && !is.function(accumulate_tile)) {
    stop("`accumulate_tile` must be NULL or a function.", call. = FALSE)
  }
  if (!form_total && !is.null(accumulate_tile)) {
    stop("`accumulate_tile` requires total-form contraction.", call. = FALSE)
  }
  if (!is.null(task_observer) && !is.function(task_observer)) {
    stop("`task_observer` must be NULL or a function.", call. = FALSE)
  }

  features <- ncol(frame$weights)
  measurements <- nrow(frame$weights)
  output_width <- if (is.null(query)) {
    physical_width
  } else {
    .query_output_width(query)
  }
  output <- if (form_total && is.null(accumulate_tile)) {
    matrix(0, measurements, output_width)
  } else {
    NULL
  }
  # First moments accumulate in a flat measurement-by-(effect within
  # partition) matrix and are reshaped to the named three-way array once, at
  # the end. The flat layout is the same storage the array uses (measurement
  # fastest, then effect, then partition), so the reshape is a relabelling;
  # what it buys is that the per-tile update is a matrix subassignment
  # rather than a three-way array subassignment, which duplicated the whole
  # durable array on every tile.
  left_first <- if (retain_first_moments) {
    matrix(0, measurements, q_left * length(left_partitions))
  } else {
    NULL
  }
  right_first <- if (retain_first_moments && !same_relation) {
    matrix(0, measurements, q_right * length(right_partitions))
  } else {
    NULL
  }
  max_features <- min(feature_block, features)
  max_rows <- min(row_tile, measurements)
  max_coordinates <- min(coordinate_tile, output_width)
  # Two-index subsetting of a column-compressed frame scans the whole row
  # axis on every call, so `weights[rows, features]` costs the same whether
  # one measurement tile is asked for or all of them. Taking the feature
  # block's columns once and slicing measurement tiles out of that block is
  # the identical submatrix at a fraction of the cost; it is worth one block
  # copy only when more than one measurement tile reads it.
  hoist_feature_columns <- length(.tile_starts(measurements, row_tile)) > 1L
  frame_block_bytes <- if (hoist_feature_columns) {
    as.double(utils::object.size(frame$weights))
  } else {
    0
  }
  relation_bytes <- 8 * max_features * (
    length(left_partitions) * q_left +
      if (same_relation) 0 else length(right_partitions) * q_right
  )
  atom_bytes <- if (form_total) 8 * max_features * output_width else 0
  atom_work_bytes <- if (!form_total) {
    0
  } else if (is.null(query)) {
    8 * max_features
  } else {
    8 * (max_features + q_left * q_right)
  }
  weight_slice_bytes <- 8 * max_rows * max_features
  atom_slice_bytes <- if (form_total) {
    8 * max_features * max_coordinates
  } else {
    0
  }
  product_bytes <- if (form_total) 8 * max_rows * max_coordinates else 0
  replacement_bytes <- if (form_total && is.null(accumulate_tile)) {
    product_bytes
  } else {
    0
  }
  # One fused product per side covers every partition at once, so the live
  # product is measurement-tile by (effect within partition) and the feature
  # panel it multiplies is one transposed copy of the relation blocks.
  left_first_bytes <- if (retain_first_moments) {
    8 * max_rows * q_left * length(left_partitions)
  } else {
    0
  }
  right_first_bytes <- if (retain_first_moments && !same_relation) {
    8 * max_rows * q_right * length(right_partitions)
  } else {
    0
  }
  first_panel_bytes <- if (retain_first_moments) relation_bytes else 0
  base_live_bytes <- frame_block_bytes + relation_bytes + atom_bytes +
    atom_work_bytes
  total_live_bytes <- if (form_total) {
    base_live_bytes + weight_slice_bytes + atom_slice_bytes +
      product_bytes + 2 * replacement_bytes
  } else {
    0
  }
  first_live_bytes <- if (retain_first_moments) {
    base_live_bytes + first_panel_bytes + weight_slice_bytes + 3 * max(
      left_first_bytes, right_first_bytes
    )
  } else {
    0
  }
  diagnostics <- list(
    feature_blocks = 0L,
    relation_reads = 0L,
    atom_count = 0L,
    contraction_tiles = 0L,
    max_relation_bytes = relation_bytes,
    max_atom_bytes = if (is.null(query)) atom_bytes else 0,
    max_query_atom_bytes = if (is.null(query)) 0 else atom_bytes,
    max_atom_work_bytes = atom_work_bytes,
    max_frame_block_bytes = frame_block_bytes,
    max_weight_slice_bytes = weight_slice_bytes,
    max_atom_slice_bytes = atom_slice_bytes,
    max_product_bytes = product_bytes,
    max_existing_slice_bytes = replacement_bytes,
    max_replacement_bytes = replacement_bytes,
    max_left_first_product_bytes = left_first_bytes,
    max_left_first_existing_bytes = left_first_bytes,
    max_left_first_replacement_bytes = left_first_bytes,
    max_right_first_product_bytes = right_first_bytes,
    max_right_first_existing_bytes = right_first_bytes,
    max_right_first_replacement_bytes = right_first_bytes,
    max_live_temporary_bytes = max(total_live_bytes, first_live_bytes),
    durable_output_bytes = if (is.null(output)) 0 else 8 * length(output),
    durable_left_first_moment_bytes = if (is.null(left_first)) 0 else
      8 * length(left_first),
    durable_right_first_moment_bytes = if (is.null(right_first)) 0 else
      8 * length(right_first),
    first_moment_sides_share_storage = retain_first_moments && same_relation,
    measurement_kind = "static-owned-buffer-accounting",
    external_accumulator_memory_measured = is.null(accumulate_tile)
  )

  for (feature_start in .tile_starts(features, feature_block)) {
    feature_ids <- feature_start:min(feature_start + feature_block - 1L, features)
    if (!is.null(task_observer)) task_observer("started", feature_ids)
    task <- tryCatch({
      read_family <- function(partitions, effects, reader, side) {
        stats::setNames(lapply(partitions, function(partition) {
          value <- reader(partition, feature_ids)
          if (!is.matrix(value) || !is.numeric(value) ||
              !identical(dim(value), c(length(effects), length(feature_ids))) ||
              any(!is.finite(value))) {
            stop(sprintf(
              "%s relation reader returned an invalid effect-by-feature block.",
              side
            ), call. = FALSE)
          }
          value
        }), partitions)
      }
      left <- read_family(left_partitions, left_effects, read_left, "Left")
      right <- if (same_relation) left else
        read_family(right_partitions, right_effects, read_right, "Right")
      .effect_form_feature_task(
        left, right, feature_ids, left_effects, right_effects,
        left_partitions, right_partitions, ordered_edges,
        codec = codec, query = query, form_atoms = form_total,
        same_relation = same_relation
      )
    }, error = function(error) {
      if (!is.null(task_observer)) task_observer("failed", feature_ids)
      stop(error)
    })

    tryCatch({
      # The feature panel is the transposed relation family for this block,
      # built once per block rather than once per measurement tile. Its
      # column order is effect within partition, which is exactly the
      # storage order of the three-way first-moment array, so one fused
      # product per side replaces one product per partition and the fill is
      # a plain matrix subassignment.
      left_panel <- if (retain_first_moments) {
        .first_moment_feature_panel(task$left_relations)
      } else {
        NULL
      }
      right_panel <- if (retain_first_moments && !same_relation) {
        .first_moment_feature_panel(task$right_relations)
      } else {
        NULL
      }
      feature_weights <- if (hoist_feature_columns) {
        frame$weights[, feature_ids, drop = FALSE]
      } else {
        NULL
      }
      for (row_start in .tile_starts(measurements, row_tile)) {
        rows <- row_start:min(row_start + row_tile - 1L, measurements)
        weight_slice <- if (hoist_feature_columns) {
          feature_weights[rows, , drop = FALSE]
        } else {
          frame$weights[rows, feature_ids, drop = FALSE]
        }

        if (retain_first_moments) {
          left_first[rows, ] <- left_first[rows, , drop = FALSE] +
            as.matrix(weight_slice %*% left_panel)
          if (!same_relation) {
            right_first[rows, ] <- right_first[rows, , drop = FALSE] +
              as.matrix(weight_slice %*% right_panel)
          }
        }

        if (form_total) for (coordinate_start in .tile_starts(
          output_width, coordinate_tile
        )) {
          coordinates <- coordinate_start:min(
            coordinate_start + coordinate_tile - 1L, output_width
          )
          atom_slice <- task$atoms[, coordinates, drop = FALSE]
          product <- as.matrix(weight_slice %*% atom_slice)
          if (any(!is.finite(product))) {
            stop("Effect-form spatial contraction produced non-finite values.",
              call. = FALSE)
          }
          existing <- replacement <- NULL
          if (is.null(accumulate_tile)) {
            existing <- output[rows, coordinates, drop = FALSE]
            replacement <- existing + product
            output[rows, coordinates] <- replacement
          } else {
            accumulate_tile(rows, coordinates, product)
          }
          diagnostics$contraction_tiles <- diagnostics$contraction_tiles + 1L
        }
      }
    }, error = function(error) {
      if (!is.null(task_observer)) task_observer("failed", feature_ids)
      stop(error)
    })
    diagnostics$feature_blocks <- diagnostics$feature_blocks + 1L
    diagnostics$relation_reads <- diagnostics$relation_reads +
      length(left_partitions) + if (same_relation) 0L else length(right_partitions)
    if (form_total) {
      diagnostics$atom_count <- diagnostics$atom_count + length(feature_ids)
    }
    if (!is.null(task_observer)) task_observer("completed", feature_ids)
  }

  if (retain_first_moments) {
    left_first <- .named_first_moment_array(
      left_first, measurements, left_effects, left_partitions
    )
    right_first <- if (same_relation) {
      left_first
    } else {
      .named_first_moment_array(
        right_first, measurements, right_effects, right_partitions
      )
    }
  }
  list(
    value = output,
    first_moments = if (retain_first_moments) {
      list(left = left_first, right = right_first)
    } else {
      NULL
    },
    mass = Matrix::rowSums(frame$weights),
    codec = codec,
    logical_shape = as.integer(c(q_left, q_right)),
    diagnostics = diagnostics
  )
}

# Support-streamed lowering for a fixed non-diagonal metric. Unlike the
# feature-additive route, this kernel never materializes feature-pair atoms or
# an m-by-p-squared pair frame. One support, its local metric, and the relation
# blocks needed for that support are live at a time.
.support_streamed_metric_contraction <- function(
    frame, metric_schedule, read_relation, partitions, effects,
    ordered_edges, query = NULL, form_total = TRUE, form_coherent = FALSE,
    retain_first_moments = FALSE, task_observer = NULL) {
  .validate_frame_for_compile(frame)
  schedule <- .validate_geometry_metric_schedule(metric_schedule)
  if (!identical(schedule$kind, "fixed_metric_before_frame") ||
      !is.function(read_relation)) {
    stop("Support-streamed execution requires a fixed metric and relation reader.",
      call. = FALSE)
  }
  metric <- .validate_neural_metric(schedule$metric, deep = FALSE)
  .validate_effect_names(effects, length(effects))
  .validate_ordered_partition_edges(
    ordered_edges, partitions, partitions, TRUE
  )
  if (!is.logical(form_total) || length(form_total) != 1L ||
      is.na(form_total) || !is.logical(form_coherent) ||
      length(form_coherent) != 1L || is.na(form_coherent) ||
      !is.logical(retain_first_moments) ||
      length(retain_first_moments) != 1L || is.na(retain_first_moments) ||
      (!form_total && !form_coherent && !retain_first_moments)) {
    stop("Support-streamed output requirements are invalid.", call. = FALSE)
  }
  if (!is.null(task_observer) && !is.function(task_observer)) {
    stop("`task_observer` must be NULL or a function.", call. = FALSE)
  }
  q <- length(effects)
  packed_width <- q * (q + 1L) / 2L
  structured_query <- !is.null(query) && .is_pair_difference_query(query)
  if (structured_query) {
    .validate_pair_difference_for_task(query, effects, effects, TRUE)
  } else if (!is.null(query) && (!is.matrix(query) || !is.numeric(query) ||
      nrow(query) != packed_width || ncol(query) < 1L ||
      any(!is.finite(query)))) {
    stop("`query` must match the finite packed effect coordinates.",
      call. = FALSE)
  }
  operators <- if (is.null(query) || structured_query) {
    NULL
  } else {
    .physical_query_operators(query, q, q, "symmetric_packed")
  }
  measurements <- nrow(frame$weights)
  output_width <- if (is.null(query)) {
    packed_width
  } else {
    .query_output_width(query)
  }
  total <- if (form_total) matrix(0, measurements, output_width) else NULL
  coherent <- if (form_coherent) {
    matrix(0, measurements, output_width)
  } else {
    NULL
  }
  first <- if (retain_first_moments) {
    array(0, c(measurements, q, length(partitions)),
      dimnames = list(NULL, effects, partitions))
  } else {
    NULL
  }
  mass <- numeric(measurements)
  left_index <- match(ordered_edges$left, partitions)
  right_index <- match(ordered_edges$right, partitions)
  metric_positions <- stats::setNames(
    seq_along(metric$support), as.character(metric$support)
  )
  diagnostics <- list(
    support_tasks = 0L,
    relation_reads = 0L,
    pair_atoms_materialized = FALSE,
    pair_frame_materialized = FALSE,
    max_support_size = 0L,
    max_relation_block_bytes = 0,
    max_metric_bytes = 0,
    max_query_work_bytes = 0,
    metric_factorizations = 0L,
    durable_output_bytes =
      (if (is.null(total)) 0 else 8 * length(total)) +
      (if (is.null(coherent)) 0 else 8 * length(coherent)),
    durable_first_moment_bytes = if (is.null(first)) 0 else 8 * length(first),
    measurement_kind = "static-owned-buffer-accounting"
  )

  project_pairs <- function(values) {
    if (is.null(query$coefficients)) values else
      drop(query$coefficients %*% values)
  }
  contract_metric <- function(left, right, K) {
    if (structured_query) {
      dl <- left[query$pair_left, , drop = FALSE] -
        left[query$pair_right, , drop = FALSE]
      dr <- right[query$pair_left, , drop = FALSE] -
        right[query$pair_right, , drop = FALSE]
      return(project_pairs(rowSums((dl %*% K) * dr)))
    }
    if (is.null(operators)) {
      return(left %*% K %*% t(right))
    }
    vapply(operators, function(H) {
      sum((t(H) %*% left %*% K) * right)
    }, numeric(1))
  }
  encode <- function(value) {
    if (is.null(query)) .svec_symmetric(value) else value
  }

  for (node_index in seq_len(measurements)) {
    node <- .frame_metric_node(frame, node_index)
    if (!is.null(task_observer)) {
      task_observer("started", node$support_positions)
    }
    tryCatch({
      local_metric_index <- unname(metric_positions[as.character(node$support)])
      if (anyNA(local_metric_index)) {
        stop("A frame support falls outside the fixed metric operator.",
          call. = FALSE)
      }
      base <- metric$value[
        local_metric_index, local_metric_index, drop = FALSE
      ]
      root_weight <- sqrt(node$weight)
      K <- if (metric$capabilities$native_diagonal) {
        diag(node$weight * diag(base), nrow = length(node$weight))
      } else {
        base * tcrossprod(root_weight)
      }
      relations <- stats::setNames(lapply(partitions, function(partition) {
        value <- read_relation(partition, node$support_positions)
        if (!is.matrix(value) || !is.numeric(value) ||
            !identical(dim(value), c(q, length(node$support_positions))) ||
            any(!is.finite(value))) {
          stop("Relation reader returned an invalid local effect block.",
            call. = FALSE)
        }
        value
      }), partitions)
      mass[[node_index]] <- sum(node$weight)
      if (retain_first_moments || form_coherent) {
        amplitudes <- lapply(relations, function(value) {
          drop(value %*% (node$weight / mass[[node_index]]))
        })
        if (retain_first_moments) {
          for (partition_index in seq_along(partitions)) {
            first[node_index, , partition_index] <-
              amplitudes[[partition_index]] * mass[[node_index]]
          }
        }
      }
      if (form_total) {
        value <- if (is.null(query)) matrix(0, q, q) else
          numeric(output_width)
        for (edge in seq_len(nrow(ordered_edges))) {
          value <- value + ordered_edges$weight[[edge]] * contract_metric(
            relations[[left_index[[edge]]]],
            relations[[right_index[[edge]]]], K
          )
        }
        total[node_index, ] <- encode(value)
      }
      if (form_coherent) {
        if (!metric$capabilities$positive_definite) {
          stop(paste0(
            "Metric-aware coherent/configuration output requires an SPD ",
            "metric on every support."
          ), call. = FALSE)
        }
        a <- node$weight / mass[[node_index]]
        factor <- tryCatch(chol(K), error = function(error) NULL)
        if (is.null(factor)) {
          stop("The composed local metric is not positive definite.",
            call. = FALSE)
        }
        solved <- backsolve(factor, forwardsolve(t(factor), a))
        denominator <- sum(a * solved)
        if (!is.finite(denominator) || denominator <= 0) {
          stop("The coherent inverse-metric norm is not positive and finite.",
            call. = FALSE)
        }
        value <- if (is.null(query)) matrix(0, q, q) else
          numeric(output_width)
        for (edge in seq_len(nrow(ordered_edges))) {
          left_amplitude <- amplitudes[[left_index[[edge]]]]
          right_amplitude <- amplitudes[[right_index[[edge]]]]
          edge_value <- if (structured_query) {
            dl <- left_amplitude[query$pair_left] -
              left_amplitude[query$pair_right]
            dr <- right_amplitude[query$pair_left] -
              right_amplitude[query$pair_right]
            project_pairs(dl * dr / denominator)
          } else if (is.null(operators)) {
            tcrossprod(left_amplitude, right_amplitude) / denominator
          } else {
            vapply(operators, function(H) {
              sum((t(H) %*% left_amplitude) * right_amplitude) /
                denominator
            }, numeric(1))
          }
          value <- value + ordered_edges$weight[[edge]] * edge_value
        }
        coherent[node_index, ] <- encode(value)
        diagnostics$metric_factorizations <-
          diagnostics$metric_factorizations + 1L
      }
      k <- length(node$support_positions)
      diagnostics$support_tasks <- diagnostics$support_tasks + 1L
      diagnostics$relation_reads <- diagnostics$relation_reads +
        length(partitions)
      diagnostics$max_support_size <- max(diagnostics$max_support_size, k)
      diagnostics$max_relation_block_bytes <- max(
        diagnostics$max_relation_block_bytes,
        8 * length(partitions) * q * k
      )
      diagnostics$max_metric_bytes <- max(
        diagnostics$max_metric_bytes, 8 * k * k
      )
      diagnostics$max_query_work_bytes <- max(
        diagnostics$max_query_work_bytes,
        8 * (q * k + k * k + if (is.null(operators)) q * q else q * k)
      )
    }, error = function(error) {
      if (!is.null(task_observer)) {
        task_observer("failed", node$support_positions)
      }
      stop(error)
    })
    if (!is.null(task_observer)) {
      task_observer("completed", node$support_positions)
    }
  }
  list(
    value = total,
    coherent = coherent,
    first_moments = if (retain_first_moments) {
      list(left = first, right = first)
    } else {
      NULL
    },
    mass = mass,
    codec = "symmetric_packed",
    logical_shape = as.integer(c(q, q)),
    diagnostics = diagnostics
  )
}

# Query-fused crossnobis for a provenance-frozen, on-demand metric schedule.
# The metric can vary by support and evaluation edge, so feature additivity is
# unavailable. The kernel streams one support and derives one local solve
# handle at a time; it never materializes pair atoms, a p-squared frame, or a
# node-by-edge factor table.
.support_streamed_scheduled_crossnobis <- function(
    frame, metric_schedule, read_relation, partitions, effects,
    ordered_edges, contrast, task_observer = NULL) {
  .validate_frame_for_compile(frame)
  schedule <- metric_schedule
  .validate_frozen_metric_schedule(schedule, deep = FALSE)
  if (!is.function(read_relation)) {
    stop("Scheduled crossnobis requires one relation reader.", call. = FALSE)
  }
  .validate_effect_names(effects, length(effects))
  .validate_ordered_partition_edges(
    ordered_edges, partitions, partitions, TRUE
  )
  contrast <- .align_contrast(contrast, effects)
  if (!is.null(task_observer) && !is.function(task_observer)) {
    stop("`task_observer` must be NULL or a function.", call. = FALSE)
  }
  for (edge in seq_len(nrow(ordered_edges))) {
    input <- ordered_edges$input_edge[[edge]]
    declared <- c(schedule$pairing$left[[input]],
      schedule$pairing$right[[input]])
    observed <- c(ordered_edges$left[[edge]], ordered_edges$right[[edge]])
    if (!identical(observed, declared) &&
        !identical(observed, rev(declared))) {
      stop("Scheduled metric records do not match ordered evaluation edges.",
        call. = FALSE)
    }
  }
  endpoints <- partitions[partitions %in%
    unique(c(ordered_edges$left, ordered_edges$right))]
  providers <- lapply(seq_len(nrow(schedule$pairing)), function(edge) {
    .metric_schedule_provider(schedule, edge)
  })
  read_node <- .frame_metric_node_accessor(frame)
  measurements <- nrow(frame$weights)
  values <- numeric(measurements)
  diagnostics <- list(
    support_tasks = 0L,
    relation_reads = 0L,
    metric_handles_derived = 0L,
    pair_atoms_materialized = FALSE,
    pair_frame_materialized = FALSE,
    metric_factor_table_retained = FALSE,
    max_support_size = 0L,
    max_relation_block_bytes = 0,
    max_local_covariance_bytes = 0,
    durable_output_bytes = 8 * measurements,
    measurement_kind = "static-owned-buffer-accounting"
  )

  for (node_index in seq_len(measurements)) {
    node <- read_node(node_index)
    if (!is.null(task_observer)) {
      task_observer("started", node$support_positions)
    }
    tryCatch({
      root_weight <- sqrt(node$weight)
      patterns <- stats::setNames(lapply(endpoints, function(partition) {
        relation <- read_relation(partition, node$support_positions)
        if (!is.matrix(relation) || !is.numeric(relation) ||
            !identical(dim(relation),
              c(length(effects), length(node$support_positions))) ||
            any(!is.finite(relation))) {
          stop("Relation reader returned an invalid local effect block.",
            call. = FALSE)
        }
        drop(contrast %*% relation) * root_weight
      }), endpoints)
      value <- 0
      # A scalar form with symmetric K is self-adjoint: the reverse half-edge
      # equals the forward half-edge exactly in the mathematical estimand.
      # Contract each declared edge once, avoiding a duplicate local solve.
      for (edge in seq_len(nrow(schedule$pairing))) {
        provider <- providers[[edge]]
        handle <- provider$at(node_index)
        edge_value <- drop(handle$form(
          matrix(patterns[[schedule$pairing$left[[edge]]]], nrow = 1L),
          matrix(patterns[[schedule$pairing$right[[edge]]]], nrow = 1L)
        ))
        value <- value + schedule$pairing$weight[[edge]] * edge_value
      }
      if (!is.finite(value)) {
        stop("Scheduled crossnobis produced a non-finite value.",
          call. = FALSE)
      }
      values[[node_index]] <- value
      k <- length(node$support_positions)
      diagnostics$support_tasks <- diagnostics$support_tasks + 1L
      diagnostics$relation_reads <- diagnostics$relation_reads +
        length(endpoints)
      diagnostics$metric_handles_derived <-
        diagnostics$metric_handles_derived + nrow(schedule$pairing)
      diagnostics$max_support_size <- max(
        diagnostics$max_support_size, k
      )
      diagnostics$max_relation_block_bytes <- max(
        diagnostics$max_relation_block_bytes,
        8 * length(endpoints) * length(effects) * k
      )
      diagnostics$max_local_covariance_bytes <- max(
        diagnostics$max_local_covariance_bytes, 8 * k * k
      )
      if (!is.null(task_observer)) {
        task_observer("completed", node$support_positions)
      }
    }, error = function(error) {
      if (!is.null(task_observer)) {
        task_observer("failed", node$support_positions)
      }
      stop(error)
    })
  }
  list(
    values = values,
    diagnostics = diagnostics,
    metric_receipts = lapply(providers, function(provider) provider$receipt()),
    endpoints_read = endpoints
  )
}

# Component-aware lowering for two-sided forms. Configuration is never a
# separate feature statistic: it is the total contraction minus the coherent
# bilinear correction built from the two retained first-moment families.
.streamed_effect_form_components <- function(
    frame, read_left, read_right = read_left,
    left_partitions, right_partitions, left_effects, right_effects,
    ordered_edges, codec = c("rectangular", "symmetric_packed"),
    same_relation = FALSE, query = NULL,
    component = c("complete", "total", "coherent", "configuration"),
    feature_block = 1024L, row_tile = 1024L, coordinate_tile = 256L,
    task_observer = NULL) {
  codec <- match.arg(codec)
  component <- match.arg(component)
  if (component == "complete" && !is.null(query)) {
    stop("A queried effect-form lowering must name one result component.",
      call. = FALSE)
  }
  needs_total <- component %in% c("complete", "total", "configuration")
  needs_coherent <- component %in% c(
    "complete", "coherent", "configuration"
  )
  memory <- .effect_form_kernel_memory_plan(
    frame,
    left_effects,
    right_effects,
    left_partitions,
    right_partitions,
    codec = codec,
    query = query,
    same_relation = same_relation,
    feature_block = feature_block,
    row_tile = row_tile,
    coordinate_tile = coordinate_tile,
    storage = "memory",
    retain_first_moments = needs_coherent,
    form_total = needs_total
  )
  streamed <- .streamed_effect_form_contraction(
    frame = frame,
    read_left = read_left,
    read_right = read_right,
    left_partitions = left_partitions,
    right_partitions = right_partitions,
    left_effects = left_effects,
    right_effects = right_effects,
    ordered_edges = ordered_edges,
    codec = codec,
    same_relation = same_relation,
    query = query,
    feature_block = feature_block,
    row_tile = row_tile,
    coordinate_tile = coordinate_tile,
    retain_first_moments = needs_coherent,
    form_total = needs_total,
    task_observer = task_observer
  )
  coherent <- if (needs_coherent) {
    .effect_form_coherent_from_first_moments(
      streamed$first_moments$left,
      streamed$first_moments$right,
      ordered_edges,
      streamed$mass,
      codec = codec,
      same_relation = same_relation,
      row_tile = row_tile,
      query = query
    )
  } else {
    NULL
  }
  value <- switch(component,
    complete = NULL,
    total = streamed$value,
    coherent = coherent$value,
    configuration = streamed$value - coherent$value
  )
  structure(list(
    value = value,
    total = if (needs_total) streamed$value else NULL,
    coherent = if (needs_coherent) coherent$value else NULL,
    first_moments = streamed$first_moments,
    mass = streamed$mass,
    codec = codec,
    logical_shape = streamed$logical_shape,
    component = component,
    memory = memory,
    diagnostics = list(
      stream = streamed$diagnostics,
      coherent = if (needs_coherent) coherent$diagnostics else NULL
    )
  ), class = "effect_form_component_contraction")
}

# Primary additive-frame lowering. Relation blocks are read exactly once per
# canonical feature block and converted to packed cross-Gram atoms before the
# sparse frame distributes them across measurement rows.
.streamed_crossgram_contraction <- function(frame, read_relation, partitions,
                                            effects, over,
                                            feature_block = 1024L,
                                            row_tile = 1024L,
                                            coordinate_tile = 256L,
                                            accumulate_tile = NULL,
                                            retain_local_relations = FALSE,
                                            query = NULL,
                                            form_total = TRUE,
                                            task_observer = NULL) {
  ordered_edges <- .ordered_partition_edges(
    over, partitions, partitions, same_relation = TRUE
  )
  streamed <- .streamed_effect_form_contraction(
    frame = frame,
    read_left = read_relation,
    left_partitions = partitions,
    right_partitions = partitions,
    left_effects = effects,
    right_effects = effects,
    ordered_edges = ordered_edges,
    codec = "symmetric_packed",
    same_relation = TRUE,
    query = query,
    feature_block = feature_block,
    row_tile = row_tile,
    coordinate_tile = coordinate_tile,
    accumulate_tile = accumulate_tile,
    retain_first_moments = retain_local_relations,
    form_total = form_total,
    task_observer = task_observer
  )
  diagnostics <- streamed$diagnostics
  diagnostics$max_local_product_bytes <-
    diagnostics$max_left_first_product_bytes
  diagnostics$max_local_existing_bytes <-
    diagnostics$max_left_first_existing_bytes
  diagnostics$max_local_replacement_bytes <-
    diagnostics$max_left_first_replacement_bytes
  diagnostics$durable_local_relation_bytes <-
    diagnostics$durable_left_first_moment_bytes
  list(
    value = streamed$value,
    local_relations = if (retain_local_relations) {
      streamed$first_moments$left
    } else {
      NULL
    },
    diagnostics = diagnostics
  )
}

.effect_form_coherent_from_first_moments <- function(
    left_first, right_first, ordered_edges, mass,
    codec = c("rectangular", "symmetric_packed"), same_relation = FALSE,
    row_tile = 1024L, write_tile = NULL, query = NULL) {
  validate_first <- function(value, side) {
    if (!is.array(value) || length(dim(value)) != 3L ||
        !is.numeric(value) || any(!is.finite(value))) {
      stop(sprintf(
        "`%s_first` must be a finite measurement-by-effect-by-partition array.",
        side
      ), call. = FALSE)
    }
    effects <- dimnames(value)[[2L]]
    partitions <- dimnames(value)[[3L]]
    if (is.null(effects) || anyNA(effects) || any(!nzchar(effects)) ||
        anyDuplicated(effects) || is.null(partitions) || anyNA(partitions) ||
        any(!nzchar(partitions)) || anyDuplicated(partitions)) {
      stop("First-moment effect and partition axes must be uniquely named.",
        call. = FALSE)
    }
    list(effects = effects, partitions = partitions)
  }
  left_axis <- validate_first(left_first, "left")
  right_axis <- validate_first(right_first, "right")
  measurements <- dim(left_first)[[1L]]
  if (dim(right_first)[[1L]] != measurements) {
    stop("Left and right first moments must share a measurement axis.",
      call. = FALSE)
  }
  .validate_ordered_partition_edges(
    ordered_edges, left_axis$partitions, right_axis$partitions, same_relation
  )
  codec <- match.arg(codec)
  if (codec == "symmetric_packed" &&
      (!same_relation || !identical(left_axis, right_axis) ||
       !identical(attr(ordered_edges, "expansion"),
         "self_adjoint_half_edges"))) {
    stop("Symmetric-packed coherent forms require a self-adjoint self form.",
      call. = FALSE)
  }
  if (!is.numeric(mass) || !(length(mass) %in% c(1L, measurements)) ||
      any(!is.finite(mass)) || any(mass <= 0)) {
    stop("`mass` must be one positive finite value or one per measurement.",
      call. = FALSE)
  }
  mass <- rep_len(mass, measurements)
  row_tile <- .validate_tile_size(row_tile, "row_tile")
  if (!is.null(write_tile) && !is.function(write_tile)) {
    stop("`write_tile` must be NULL or a function.", call. = FALSE)
  }
  q_left <- dim(left_first)[[2L]]
  q_right <- dim(right_first)[[2L]]
  physical_width <- if (codec == "rectangular") {
    q_left * q_right
  } else {
    q_left * (q_left + 1L) / 2L
  }
  structured_query <- !is.null(query) && .is_pair_difference_query(query)
  if (structured_query) {
    .validate_pair_difference_for_task(
      query, left_axis$effects, right_axis$effects, same_relation
    )
  } else if (!is.null(query) && (!is.matrix(query) || !is.numeric(query) ||
      nrow(query) != physical_width || ncol(query) < 1L ||
      any(!is.finite(query)))) {
    stop("`query` must match the finite physical form coordinates.",
      call. = FALSE)
  }
  output_width <- if (is.null(query)) {
    physical_width
  } else {
    .query_output_width(query)
  }
  output <- if (is.null(write_tile)) {
    matrix(0, measurements, output_width)
  } else {
    NULL
  }
  left_index <- match(ordered_edges$left, left_axis$partitions)
  right_index <- match(ordered_edges$right, right_axis$partitions)
  max_rows <- min(row_tile, measurements)
  max_work_bytes <- 8 * max_rows
  max_operand_bytes <- if (is.null(query)) {
    0
  } else if (structured_query) {
    8 * (2 * max_rows * length(query$pair_left))
  } else {
    8 * (max_rows * q_left + max_rows * q_right + q_left * q_right)
  }
  tile_count <- 0L

  for (row_start in .tile_starts(measurements, row_tile)) {
    rows <- row_start:min(row_start + row_tile - 1L, measurements)
    tile <- matrix(0, length(rows), output_width)
    if (is.null(query)) {
      coordinate <- 0L
      for (column in seq_len(q_right)) {
        form_rows <- if (codec == "symmetric_packed") column:q_left else
          seq_len(q_left)
        for (row in form_rows) {
          coordinate <- coordinate + 1L
          work <- numeric(length(rows))
          for (edge in seq_len(nrow(ordered_edges))) {
            work <- work + ordered_edges$weight[[edge]] *
              left_first[rows, row, left_index[[edge]]] *
              right_first[rows, column, right_index[[edge]]]
          }
          work <- work / mass[rows]
          if (codec == "symmetric_packed" && row != column) {
            work <- sqrt(2) * work
          }
          tile[, coordinate] <- work
        }
      }
    } else if (structured_query) {
      pair_tile <- matrix(0, length(rows), length(query$pair_left))
      for (edge in seq_len(nrow(ordered_edges))) {
        left <- matrix(
          left_first[rows, , left_index[[edge]], drop = FALSE],
          nrow = length(rows), ncol = q_left
        )
        right <- matrix(
          right_first[rows, , right_index[[edge]], drop = FALSE],
          nrow = length(rows), ncol = q_right
        )
        dl <- left[, query$pair_left, drop = FALSE] -
          left[, query$pair_right, drop = FALSE]
        dr <- right[, query$pair_left, drop = FALSE] -
          right[, query$pair_right, drop = FALSE]
        pair_tile <- pair_tile + ordered_edges$weight[[edge]] * (dl * dr)
      }
      tile[, ] <- if (is.null(query$coefficients)) {
        pair_tile / mass[rows]
      } else {
        (pair_tile %*% t(query$coefficients)) / mass[rows]
      }
    } else {
      for (view in seq_len(ncol(query))) {
        work <- numeric(length(rows))
        operator <- .physical_query_operator(
          query[, view], q_left, q_right, codec
        )
        for (edge in seq_len(nrow(ordered_edges))) {
          left <- matrix(
            left_first[rows, , left_index[[edge]], drop = FALSE],
            nrow = length(rows), ncol = q_left
          )
          right <- matrix(
            right_first[rows, , right_index[[edge]], drop = FALSE],
            nrow = length(rows), ncol = q_right
          )
          work <- work + ordered_edges$weight[[edge]] *
            rowSums(left * (right %*% t(operator)))
        }
        tile[, view] <- work / mass[rows]
      }
    }
    if (any(!is.finite(tile))) {
      stop("Coherent effect-form contraction produced non-finite values.",
        call. = FALSE)
    }
    if (is.null(write_tile)) {
      output[rows, ] <- tile
    } else {
      write_tile(rows, seq_len(output_width), tile)
    }
    tile_count <- tile_count + 1L
  }

  list(
    value = output,
    diagnostics = list(
      row_tile = row_tile,
      tile_count = tile_count,
      max_tile_bytes = 8 * max_rows * output_width,
      max_work_bytes = max_work_bytes,
      max_operand_bytes = max_operand_bytes,
      durable_output_bytes = if (is.null(output)) 0 else 8 * length(output),
      measurement_kind = "static-owned-buffer-accounting"
    )
  )
}

.coherent_geometry_from_local <- function(local_relations, over, mass,
                                          row_tile = 1024L,
                                          write_tile = NULL, query = NULL) {
  if (!is.array(local_relations) || length(dim(local_relations)) != 3L ||
      !is.numeric(local_relations) || any(!is.finite(local_relations))) {
    stop("`local_relations` must be a finite measurement-by-effect-by-partition array.",
      call. = FALSE)
  }
  partitions <- dimnames(local_relations)[[3L]]
  effects <- dimnames(local_relations)[[2L]]
  if (is.null(partitions) || is.null(effects)) {
    stop("Local relation effect and partition axes must be named.",
      call. = FALSE)
  }
  edges <- .ordered_partition_edges(over, partitions, partitions, TRUE)
  .effect_form_coherent_from_first_moments(
    local_relations, local_relations, edges, mass,
    codec = "symmetric_packed", same_relation = TRUE,
    row_tile = row_tile, write_tile = write_tile, query = query
  )
}
