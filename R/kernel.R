# Bounded additive contraction ----------------------------------------------

.tiled_contraction <- function(weights, atoms, row_tile, coordinate_tile,
                               feature_tile, write_tile = NULL) {
  .check_matrix(weights, "weights", what = "a finite numeric matrix")
  .check_matrix(atoms, "atoms", what = "a finite numeric matrix")
  if (ncol(weights) != nrow(atoms)) {
    .input_error("The feature dimension of `weights` and `atoms` must agree.")
  }
  if (nrow(weights) < 1L || ncol(weights) < 1L || ncol(atoms) < 1L) {
    .input_error("Contraction inputs must have positive dimensions.")
  }
  row_tile <- .validate_tile_size(row_tile, "row_tile")
  coordinate_tile <- .validate_tile_size(coordinate_tile, "coordinate_tile")
  feature_tile <- .validate_tile_size(feature_tile, "feature_tile")
  if (!is.null(write_tile) && !is.function(write_tile)) {
    .input_error("`write_tile` must be NULL or a function.")
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
    .input_error(
      "The streamed effect-form lowering requires an additive diagonal frame."
    )
  }
  if (!is.function(read_left) || !is.function(read_right)) {
    .input_error("Effect-form relation readers must be functions.")
  }
  .check_flag(same_relation, "same_relation")
  .check_flag(retain_first_moments, "retain_first_moments")
  .check_flag(form_total, "form_total")
  if (!form_total && !retain_first_moments) {
    .input_error("Streaming must form total output, first moments, or both.")
  }
  if (same_relation && (!identical(left_partitions, right_partitions) ||
      !identical(left_effects, right_effects))) {
    .input_error(
      "A shared relation reader requires identical partition and effect axes."
    )
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
      .input_error(
        "Symmetric-packed streaming requires a self-adjoint self form."
      )
    }
    q_left * (q_left + 1L) / 2L
  }
  .validate_task_query(
    query, physical_width, left_effects, right_effects, same_relation
  )
  feature_block <- .validate_tile_size(feature_block, "feature_block")
  row_tile <- .validate_tile_size(row_tile, "row_tile")
  coordinate_tile <- .validate_tile_size(coordinate_tile, "coordinate_tile")
  if (!is.null(accumulate_tile) && !is.function(accumulate_tile)) {
    .input_error("`accumulate_tile` must be NULL or a function.")
  }
  if (!form_total && !is.null(accumulate_tile)) {
    .input_error("`accumulate_tile` requires total-form contraction.")
  }
  if (!is.null(task_observer) && !is.function(task_observer)) {
    .input_error("`task_observer` must be NULL or a function.")
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
          if (!.is_finite_matrix(value) ||
              !identical(dim(value), c(length(effects), length(feature_ids)))) {
            .input_error(sprintf(
              "%s relation reader returned an invalid effect-by-feature block.",
              side
            ))
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
            .input_error(
              paste0("Effect-form spatial contraction produced non-finite values.",
          " Finite inputs overflowed double precision during the computation; rescale the responses (for example to unit variance) before building the relation.")
            )
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

# Grouped-partner Gram route for support-streamed pair differences ----------
#
# `.support_streamed_metric_contraction()` evaluates one row-difference product
# per requested effect pair per ordered edge per support. That is the right
# association when a query names a few pairs, and the wrong one when it names
# most of them: the same numbers are the entries of one effect-by-effect matrix
# per support, formed by grouping the edge sum on its left partition. The
# grouping is a re-association of a sum of like-signed terms, so it introduces
# no cancellation the estimand does not already carry, and it factors no
# metric, so an indefinite local metric is contracted rather than refused --
# which is what the reference loop does too. `src/support-metric-gram.cpp`
# carries the derivation.
#
# Everything below is admission. The route declines by returning NULL, which is
# not an error and not a refusal: the reference loop then runs and produces the
# behaviour, the diagnostics, and any refusal the caller would have seen
# anyway. Declining is therefore always safe, and every condition it cannot
# cheaply prove is a decline.

# Reading whole partition blocks once, rather than one overlapping local block
# per support, is strictly less I/O but it makes the blocks resident, alongside
# the partner-sum staging they are folded into. This is the ceiling on that
# residency; above it the route declines and the reference loop keeps its
# per-support bound.
.support_gram_block_budget_bytes <- function() 256 * 1024^2

# The one resource in which the Gram route is genuinely worse than the
# difference route: it forms a q-by-q accumulator per measurement, where the
# difference route's live working set is only the pair-by-support block. At
# large q that is real memory, and it is what bounds the route in practice --
# on the frames measured below, this ceiling is reached before the cost
# comparison ever prefers the difference route. Above it the route declines.
.support_gram_accumulator_budget_bytes <- function() 64 * 1024^2

# Cost of the interpreted difference route for one measurement, in operation
# equivalents. An operation count alone gets this wrong, because two of the
# three terms are not flops at all. Measured on a 1000-node searchlight frame,
# four partitions, mean support 26.75, R 4.5.1 with Accelerate and threads
# pinned: a fixed per-edge dispatch and allocation cost of about 1.7e5
# equivalents; a term of about 64 equivalents per element of the q-by-support
# block, which is the two row differences the route takes before it contracts
# anything and the reason a one-pair query does not get cheaper as q grows; and
# the 2 k m^2 the flop count already had. Fitting only the third term is what
# sent a 28-of-4950 query down a route measured 15x slower than the one it
# declined.
.support_gram_difference_operations <- function(edges, q, pairs, support) {
  edges * (1.7e5 + 64 * q * support + 2 * pairs * support^2)
}

# The difference route is interpreted R and the Gram route is compiled, so a
# unit of work does not cost the same in each and their operation counts are
# not comparable as they stand. This is the measured ratio, from the marginal
# seconds per operation equivalent of each route on the frame above: 2.8e-10 s
# for the difference route against 4.5e-12 s for the Gram route, a ratio of 62,
# taken down to 32 so the constant is not the most favourable reading of the
# measurement. It reproduces every routing decision measured: across q in
# 8, 32, 100, 200, 400, 800, 1200 and 1600 against one, 28 and all pairs, the
# Gram route was faster by between 4.1x (q = 1600, one pair) and 230x
# (q = 100, all pairs), and this comparison admits all of them, with the margin
# narrowing to 1.4x exactly where the accumulator budget above cuts in.
.support_gram_interpreted_penalty <- function() 32

.support_gram_admits_metric <- function(metric) {
  !isTRUE(metric$capabilities$native_diagonal) &&
    is.matrix(metric$value) && is.double(metric$value) &&
    nrow(metric$value) == ncol(metric$value) && nrow(metric$value) >= 1L
}

# Two layouts of the same arithmetic. Below the threshold the gather is one
# contiguous copy per support feature and the collapse is one small product per
# partition; above it the wider per-partition gather pays for itself, because
# the collapse becomes a single symmetric rank-2k update whose inner dimension
# is the whole support-by-partition extent. The threshold is measured on one
# BLAS rather than derived, and the two layouts agree to rounding, so a wrong
# call here costs a few per cent and nothing else.
.support_gram_layout <- function(q) if (q >= 64L) 2L else 1L

# The half-edge expansion carries each unordered partition pair twice, forward
# then reversed at one weight, so the odd rows name every unordered pair once.
# The grouping needs each pair oriented on its smaller endpoint; the ordered
# sum is symmetric in the pair, so orienting it is free. A self pair has no
# smaller endpoint and is declined.
.support_gram_unordered_edges <- function(ordered_edges, partitions) {
  forward <- seq.int(1L, nrow(ordered_edges), by = 2L)
  left <- match(ordered_edges$left[forward], partitions)
  right <- match(ordered_edges$right[forward], partitions)
  if (anyNA(left) || anyNA(right) || any(left == right)) return(NULL)
  list(left = pmin(left, right), right = pmax(left, right),
    weight = ordered_edges$weight[forward])
}

.support_metric_gram_route <- function(frame, metric, read_relation,
                                       partitions, effects, ordered_edges,
                                       query, structured_query, form_total,
                                       form_coherent, retain_first_moments,
                                       task_observer, measurements, q,
                                       output_width) {
  if (!form_total || form_coherent || retain_first_moments) return(NULL)
  if (is.null(query) || !isTRUE(structured_query)) return(NULL)
  if (!is.null(query$coefficients)) return(NULL)
  if (!identical(attr(ordered_edges, "expansion", exact = TRUE),
      "self_adjoint_half_edges")) {
    return(NULL)
  }
  if (!.support_gram_admits_metric(metric)) return(NULL)
  weights <- frame$weights
  if (!inherits(weights, "sparseMatrix")) return(NULL)
  n_features <- ncol(weights)
  n_partitions <- length(partitions)
  n_pairs <- length(query$pair_left)
  if (n_pairs != output_width || n_partitions < 2L || q < 2L) return(NULL)

  # The resident blocks plus the partner-sum staging they are folded into.
  block_bytes <- 8 * (n_partitions + 2 * (n_partitions - 1L)) * q * n_features
  if (block_bytes > .support_gram_block_budget_bytes()) return(NULL)
  accumulator_bytes <- 8 * as.numeric(q) * q
  if (accumulator_bytes > .support_gram_accumulator_budget_bytes()) {
    return(NULL)
  }

  edges <- .support_gram_unordered_edges(ordered_edges, partitions)
  if (is.null(edges)) return(NULL)

  rows <- methods::as(weights, "RsparseMatrix")
  if (!inherits(rows, "dgRMatrix") || length(rows@x) < 1L) return(NULL)
  if (any(rows@x <= 0)) return(NULL)
  if (length(rows@p) != measurements + 1L) return(NULL)
  support_sizes <- diff(rows@p)
  if (any(support_sizes < 1L)) return(NULL)

  # Which route is cheaper, in operation equivalents corrected for the fact
  # that one of them is interpreted. The Gram route is 2 (P-1) q m^2 for the
  # metric products plus 2 (P-1) q^2 m for the collapse; the difference route
  # carries two terms beyond its flops, which is why it is modelled rather than
  # counted. The q^2 term is what eventually makes the Gram route the worse
  # choice, and the accumulator budget above reaches that point first.
  mean_support <- length(rows@x) / measurements
  partners <- n_partitions - 1L
  difference_cost <- .support_gram_difference_operations(
    nrow(ordered_edges), q, n_pairs, mean_support
  )
  gram_cost <- 2 * partners * q * mean_support^2 +
    2 * partners * q^2 * mean_support
  if (gram_cost >= difference_cost * .support_gram_interpreted_penalty()) {
    return(NULL)
  }

  metric_row <- if (identical(metric$support, frame$domain$feature_ids)) {
    seq_len(n_features)
  } else {
    match(as.character(frame$domain$feature_ids),
      as.character(metric$support))
  }
  if (length(metric_row) != n_features || anyNA(metric_row)) return(NULL)

  feature_ids <- seq_len(n_features)
  blocks <- vector("list", n_partitions)
  for (index in seq_len(n_partitions)) {
    value <- read_relation(partitions[[index]], feature_ids)
    if (!.is_finite_matrix(value) ||
        !identical(dim(value), c(q, n_features))) {
      return(NULL)
    }
    blocks[[index]] <- value
  }

  native <- .support_metric_gram_pairs_cpp(
    rows@p, rows@j, rows@x, metric$value, as.integer(metric_row - 1L), blocks,
    as.integer(edges$left - 1L), as.integer(edges$right - 1L),
    as.numeric(edges$weight), as.integer(query$pair_left - 1L),
    as.integer(query$pair_right - 1L), .support_gram_layout(q)
  )

  if (!is.null(task_observer)) {
    supports <- split(as.integer(rows@j + 1L),
      rep.int(seq_len(measurements), support_sizes))
    for (node_index in seq_len(measurements)) {
      support_positions <- supports[[node_index]]
      task_observer("started", support_positions)
      task_observer("completed", support_positions)
    }
  }

  max_support <- as.integer(native$max_support)
  list(
    value = native$value,
    coherent = NULL,
    first_moments = NULL,
    mass = native$mass,
    codec = "symmetric_packed",
    logical_shape = as.integer(c(q, q)),
    diagnostics = list(
      support_tasks = measurements,
      relation_reads = n_partitions,
      pair_atoms_materialized = FALSE,
      pair_frame_materialized = FALSE,
      max_support_size = max_support,
      max_relation_block_bytes = block_bytes,
      max_metric_bytes = 8 * max_support^2,
      # The live working set of one support: the gathered partition and
      # partner slabs, the composed local metric, the metric products, and the
      # staged matrices the kernel fills before writing.
      max_query_work_bytes = 8 * (
        (native$stage_rows + 2 * partners * q) *
          as.numeric(max_support) +
          max_support^2 +
          native$node_block * as.numeric(q) * q
      ),
      metric_factorizations = 0L,
      durable_output_bytes = 8 * length(native$value),
      durable_first_moment_bytes = 0,
      measurement_kind = "static-owned-buffer-accounting"
    )
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
    .input_error(
      "Support-streamed execution requires a fixed metric and relation reader."
    )
  }
  metric <- .validate_neural_metric(schedule$metric, deep = FALSE)
  .validate_effect_names(effects, length(effects))
  .validate_ordered_partition_edges(
    ordered_edges, partitions, partitions, TRUE
  )
  if (!.is_flag(form_total) || !.is_flag(form_coherent) ||
      !.is_flag(retain_first_moments) ||
      (!form_total && !form_coherent && !retain_first_moments)) {
    .input_error("Support-streamed output requirements are invalid.")
  }
  if (!is.null(task_observer) && !is.function(task_observer)) {
    .input_error("`task_observer` must be NULL or a function.")
  }
  q <- length(effects)
  packed_width <- q * (q + 1L) / 2L
  structured_query <- .validate_task_query(
    query, packed_width, effects, effects, TRUE,
    what = "packed effect coordinates"
  )
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
  # One execution route, not a second estimand: when the query asks for most
  # of the effect pairs the same numbers come out of one small Gram matrix per
  # support instead of one row-difference product per pair per edge. The
  # attempt declines, silently and without side effects, on anything it does
  # not recognize, and the reference loop below then runs unchanged.
  native <- .support_metric_gram_route(
    frame, metric, read_relation, partitions, effects, ordered_edges,
    query = query, structured_query = structured_query,
    form_total = form_total, form_coherent = form_coherent,
    retain_first_moments = retain_first_moments,
    task_observer = task_observer, measurements = measurements, q = q,
    output_width = output_width
  )
  if (!is.null(native)) return(native)

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
        .contract_error(
          "A frame support falls outside the fixed metric operator."
        )
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
        if (!.is_finite_matrix(value) ||
            !identical(dim(value), c(q, length(node$support_positions)))) {
          .input_error(
            "Relation reader returned an invalid local effect block."
          )
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
          .input_error(paste0(
            "Metric-aware coherent/configuration output requires an SPD ",
            "metric on every support."
          ))
        }
        a <- node$weight / mass[[node_index]]
        factor <- tryCatch(chol(K), error = function(error) NULL)
        if (is.null(factor)) {
          .input_error("The composed local metric is not positive definite.")
        }
        solved <- backsolve(factor, forwardsolve(t(factor), a))
        denominator <- sum(a * solved)
        if (!is.finite(denominator) || denominator <= 0) {
          .invariant_error(
            "The coherent inverse-metric norm is not positive and finite."
          )
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
    .input_error("Scheduled crossnobis requires one relation reader.")
  }
  .validate_effect_names(effects, length(effects))
  .validate_ordered_partition_edges(
    ordered_edges, partitions, partitions, TRUE
  )
  contrast <- .align_contrast(contrast, effects)
  if (!is.null(task_observer) && !is.function(task_observer)) {
    .input_error("`task_observer` must be NULL or a function.")
  }
  for (edge in seq_len(nrow(ordered_edges))) {
    input <- ordered_edges$input_edge[[edge]]
    declared <- c(schedule$pairing$left[[input]],
      schedule$pairing$right[[input]])
    observed <- c(ordered_edges$left[[edge]], ordered_edges$right[[edge]])
    if (!identical(observed, declared) &&
        !identical(observed, rev(declared))) {
      .contract_error(
        "Scheduled metric records do not match ordered evaluation edges."
      )
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
        if (!.is_finite_matrix(relation) ||
            !identical(dim(relation), c(length(effects), length(node$support_positions)))) {
          .input_error(
            "Relation reader returned an invalid local effect block."
          )
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
        .invariant_error("Scheduled crossnobis produced a non-finite value.")
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
    .input_error(
      "A queried effect-form lowering must name one result component."
    )
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
        !.is_finite_numeric(value)) {
      .input_error(sprintf(
        "`%s_first` must be a finite measurement-by-effect-by-partition array.",
        side
      ))
    }
    effects <- dimnames(value)[[2L]]
    partitions <- dimnames(value)[[3L]]
    if (is.null(effects) || anyNA(effects) || any(!nzchar(effects)) ||
        anyDuplicated(effects) || is.null(partitions) || anyNA(partitions) ||
        any(!nzchar(partitions)) || anyDuplicated(partitions)) {
      .input_error(
        "First-moment effect and partition axes must be uniquely named."
      )
    }
    list(effects = effects, partitions = partitions)
  }
  left_axis <- validate_first(left_first, "left")
  right_axis <- validate_first(right_first, "right")
  measurements <- dim(left_first)[[1L]]
  if (dim(right_first)[[1L]] != measurements) {
    .input_error("Left and right first moments must share a measurement axis.")
  }
  .validate_ordered_partition_edges(
    ordered_edges, left_axis$partitions, right_axis$partitions, same_relation
  )
  codec <- match.arg(codec)
  if (codec == "symmetric_packed" &&
      (!same_relation || !identical(left_axis, right_axis) ||
       !identical(attr(ordered_edges, "expansion"),
         "self_adjoint_half_edges"))) {
    .input_error(
      "Symmetric-packed coherent forms require a self-adjoint self form."
    )
  }
  if (!.is_finite_numeric(mass) || !(length(mass) %in% c(1L, measurements)) ||
      any(mass <= 0)) {
    .input_error(
      "`mass` must be one positive finite value or one per measurement."
    )
  }
  mass <- rep_len(mass, measurements)
  row_tile <- .validate_tile_size(row_tile, "row_tile")
  if (!is.null(write_tile) && !is.function(write_tile)) {
    .input_error("`write_tile` must be NULL or a function.")
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
    .contract_error("`query` must match the finite physical form coordinates.")
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
      tile <- .coherent_effect_form_atoms_cpp(
        left_first, right_first,
        as.integer(left_index), as.integer(right_index),
        as.numeric(ordered_edges$weight), as.numeric(mass),
        as.integer(rows[[1L]]), as.integer(length(rows)),
        identical(codec, "symmetric_packed")
      )
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
      .input_error(
        paste0("Coherent effect-form contraction produced non-finite values.",
          " Finite inputs overflowed double precision during the computation; rescale the responses (for example to unit variance) before building the relation.")
      )
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
      !.is_finite_numeric(local_relations)) {
    .input_error(paste0(
      "`local_relations` must be a finite measurement-by-effect-by-partition ",
      "array."
    ))
  }
  partitions <- dimnames(local_relations)[[3L]]
  effects <- dimnames(local_relations)[[2L]]
  if (is.null(partitions) || is.null(effects)) {
    .input_error("Local relation effect and partition axes must be named.")
  }
  edges <- .ordered_partition_edges(over, partitions, partitions, TRUE)
  .effect_form_coherent_from_first_moments(
    local_relations, local_relations, edges, mass,
    codec = "symmetric_packed", same_relation = TRUE,
    row_tile = row_tile, write_tile = write_tile, query = query
  )
}
