# Requested-edge measurement-form contraction ------------------------------

.validate_raw_measurement_task <- function(task) {
  .validate_evidence_task(task)
  if (task$materialization$kind != "measurement_form" ||
      task$experimental_boundary$state != "closed" ||
      task$neural_boundary$state != "open") {
    .input_error(paste0(
      "Raw measurement contraction requires closed experimental and open ",
      "neural boundaries."
    ))
  }
  normalizer <- task$stages$normalization$operation
  transform <- task$stages$transform$operation
  if (!identical(normalizer$kind, "inner_product") ||
      !identical(transform$kind, "identity")) {
    .input_error(paste0(
      "Raw measurement blocks require bilinear inner-product stages; ",
      "normalization or nonlinear transformation must remain explicit."
    ))
  }
  if (!is.null(task$materialization$projection)) {
    .input_error(
      "Complete raw measurement blocks cannot carry an effect-form projection."
    )
  }
  task
}

# The raw-measurement half of the two-sided source session. The task is
# admitted here, in the layer that owns it, and the session wiring itself is
# `.open_two_sided_source_session()` in R/source.R.
.open_evidence_task_source_session <- function(
    task, open_descriptor = .open_source_descriptor, shared_opener = NULL) {
  .open_two_sided_source_session(
    .validate_raw_measurement_task(task),
    "effect_evidence_source_session", "Evidence",
    open_descriptor = open_descriptor, shared_opener = shared_opener
  )
}

.measurement_requested_nodes <- function(task, side = c("left", "right")) {
  side <- match.arg(side)
  edges <- task$neural_boundary$edges$edges
  unique(edges[[side]])
}

.measurement_frame_width <- function(frame, nodes) {
  sum(vapply(frame$legs[nodes], function(leg) nrow(leg$operator), integer(1)))
}

.measurement_h_factorization <- function(h, tolerance = 1e-12) {
  if (!.is_number(tolerance) || tolerance <= 0) {
    .input_error("Factorization tolerance must be one positive finite value.")
  }
  decomposition <- svd(h, nu = min(dim(h)), nv = min(dim(h)))
  threshold <- if (!length(decomposition$d)) 0 else
    max(decomposition$d) * tolerance
  retained <- which(decomposition$d > threshold)
  if (!length(retained)) {
    return(list(
      left = matrix(0, nrow(h), 0L),
      right = matrix(0, ncol(h), 0L),
      rank = 0L,
      singular_values = decomposition$d,
      tolerance = tolerance
    ))
  }
  root <- sqrt(decomposition$d[retained])
  list(
    left = sweep(decomposition$u[, retained, drop = FALSE], 2L, root, `*`),
    right = sweep(decomposition$v[, retained, drop = FALSE], 2L, root, `*`),
    rank = as.integer(length(retained)),
    singular_values = decomposition$d,
    tolerance = tolerance
  )
}

.measurement_route_plan <- function(
    task,
    route = c("auto", "forward_k", "pull_h", "factorized_h",
              "scalar_stack", "multivariate_blocks"),
    feature_block = 1024L, edge_tile = 128L, rank_tolerance = 1e-12) {
  task <- .validate_raw_measurement_task(task)
  route <- match.arg(route)
  feature_block <- .validate_tile_size(feature_block, "feature_block")
  edge_tile <- .validate_tile_size(edge_tile, "edge_tile")
  h <- as.matrix(task$experimental_boundary$query$operator)
  factorization <- .measurement_h_factorization(h, rank_tolerance)
  left_frame <- task$neural_boundary$left_frame
  right_frame <- task$neural_boundary$right_frame
  measurement_edges <- task$neural_boundary$edges$edges
  left_nodes <- unique(measurement_edges$left)
  right_nodes <- unique(measurement_edges$right)
  left_ranks <- vapply(left_frame$legs[left_nodes], function(x) {
    nrow(x$operator)
  }, integer(1))
  right_ranks <- vapply(right_frame$legs[right_nodes], function(x) {
    nrow(x$operator)
  }, integer(1))
  scalar_eligible <- all(left_ranks == 1L) && all(right_ranks == 1L)
  blocks <- nrow(measurement_edges)
  partition_edges <- nrow(task$ordered_partition_products)
  q_left <- nrow(h)
  q_right <- ncol(h)
  block_coordinates <- sum(vapply(seq_len(blocks), function(index) {
    left_ranks[[match(measurement_edges$left[[index]], left_nodes)]] *
      right_ranks[[match(measurement_edges$right[[index]], right_nodes)]]
  }, numeric(1)))
  costs <- c(
    forward_k = partition_edges * block_coordinates * q_left * q_right,
    pull_h = partition_edges * sum(vapply(seq_len(blocks), function(index) {
      right_rank <- right_ranks[[match(
        measurement_edges$right[[index]], right_nodes
      )]]
      left_rank <- left_ranks[[match(
        measurement_edges$left[[index]], left_nodes
      )]]
      q_left * q_right * right_rank + q_left * left_rank * right_rank
    }, numeric(1))),
    factorized_h = partition_edges * factorization$rank *
      (q_left * sum(left_ranks) + q_right * sum(right_ranks) +
       block_coordinates),
    scalar_stack = if (scalar_eligible) {
      partition_edges * (q_left * q_right * length(right_nodes) +
        q_left * blocks)
    } else {
      Inf
    },
    multivariate_blocks = partition_edges * block_coordinates *
      (q_left + q_right)
  )
  eligible <- c(
    forward_k = TRUE,
    pull_h = TRUE,
    factorized_h = factorization$rank < min(dim(h)),
    scalar_stack = scalar_eligible,
    multivariate_blocks = TRUE
  )
  chosen <- if (route == "auto") {
    candidates <- names(costs)[eligible]
    candidates[[which.min(costs[candidates])]]
  } else {
    if (!isTRUE(eligible[[route]])) {
      .input_error(sprintf(
        "Measurement route `%s` is not eligible for this task.",
        route))
    }
    route
  }
  semantic <- list(
    schema_version = 1L,
    task_id = task$task_id,
    route = chosen,
    eligible = eligible,
    costs = costs,
    factor_rank = factorization$rank,
    rank_tolerance = rank_tolerance,
    feature_block = feature_block,
    edge_tile = edge_tile
  )
  structure(c(semantic, list(
    signature = .sha256_signature(semantic)
  )), class = "effect_measurement_route_plan")
}

.measurement_kernel_memory_plan <- function(task, route_plan,
                                            storage = c("memory", "block"),
                                            budget_bytes = NULL) {
  task <- .validate_raw_measurement_task(task)
  if (!inherits(route_plan, "effect_measurement_route_plan")) {
    .input_error("Measurement memory planning requires a route plan.")
  }
  storage <- match.arg(storage)
  left_frame <- task$neural_boundary$left_frame
  right_frame <- task$neural_boundary$right_frame
  left_nodes <- .measurement_requested_nodes(task, "left")
  right_nodes <- .measurement_requested_nodes(task, "right")
  left_width <- .measurement_frame_width(left_frame, left_nodes)
  right_width <- .measurement_frame_width(right_frame, right_nodes)
  q_left <- length(task$spaces$experimental_left$coordinates)
  q_right <- length(task$spaces$experimental_right$coordinates)
  f_left <- min(route_plan$feature_block, task$left_relation$n_features)
  f_right <- min(route_plan$feature_block, task$right_relation$n_features)
  frame_bytes <- sum(vapply(c(
    left_frame$legs[left_nodes], right_frame$legs[right_nodes]
  ), function(x) {
    .dense_double_matrix_bytes(nrow(x$operator), ncol(x$operator))
  }, numeric(1)))
  relation_buffer <- .dense_double_matrix_bytes(q_left, f_left) +
    if (task$same_relation) 0 else
      .dense_double_matrix_bytes(q_right, f_right)
  factor_slice <- .dense_double_matrix_bytes(left_width, f_left) +
    .dense_double_matrix_bytes(right_width, f_right)
  projected_state <- length(task$left_relation$partitions) *
    .dense_double_matrix_bytes(q_left, left_width) +
    if (task$same_relation && identical(left_nodes, right_nodes) &&
        identical(left_frame$signature, right_frame$signature)) {
      0
    } else {
      length(task$right_relation$partitions) *
        .dense_double_matrix_bytes(q_right, right_width)
    }
  h_factor <- if (route_plan$route == "factorized_h") {
    .dense_double_matrix_bytes(q_left + q_right, route_plan$factor_rank)
  } else {
    .dense_double_matrix_bytes(q_left, q_right)
  }
  edge_table <- task$neural_boundary$edges$edges
  block_lengths <- vapply(seq_len(nrow(edge_table)), function(index) {
    nrow(left_frame$legs[[edge_table$left[[index]]]]$operator) *
      nrow(right_frame$legs[[edge_table$right[[index]]]]$operator)
  }, numeric(1))
  output_buffer <- sum(block_lengths) * 8
  largest_block <- max(block_lengths) * 8
  edge_buffer <- min(route_plan$edge_tile, nrow(edge_table)) *
    max(q_left, q_right, 1L) * 8 + largest_block
  plan <- memory_plan(
    frame_bytes = frame_bytes,
    resident_source_bytes = as.double(utils::object.size(
      list(task$left_relation$sources,
        if (task$same_relation) NULL else task$right_relation$sources)
    )),
    relation_block_bytes = relation_buffer,
    atom_block_bytes = factor_slice + projected_state + h_factor,
    output_bytes = if (storage == "memory") output_buffer else 0,
    contraction_bytes = edge_buffer,
    replacement_copy_bytes = if (storage == "memory") 2 * largest_block else
      largest_block,
    workers = 1L,
    n_active = 1L,
    safety_factor = 1.25,
    budget_bytes = budget_bytes
  )
  structure(list(
    plan = plan,
    buffers = c(
      relation = relation_buffer,
      factor_slice = factor_slice,
      projected_state = projected_state,
      query_or_factor = h_factor,
      edge = edge_buffer,
      output = if (storage == "memory") output_buffer else 0,
      replacement = if (storage == "memory") 2 * largest_block else largest_block
    ),
    storage = storage,
    route = route_plan$route
  ), class = "effect_measurement_memory_plan")
}

.project_measurement_family <- function(task, session, side = c("left", "right"),
                                        nodes, feature_block) {
  side <- match.arg(side)
  relation <- if (side == "left") task$left_relation else task$right_relation
  frame <- if (side == "left") task$neural_boundary$left_frame else
    task$neural_boundary$right_frame
  nodes <- unique(nodes)
  legs <- frame$legs[nodes]
  widths <- vapply(legs, function(x) nrow(x$operator), integer(1))
  ends <- cumsum(widths)
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  ranges <- Map(seq.int, starts, ends)
  names(ranges) <- nodes
  operator <- do.call(rbind, lapply(legs, `[[`, "operator"))
  q <- length(relation$effect_space$coordinates)
  largest_feature_block <- min(feature_block, relation$n_features)
  projected <- stats::setNames(vector("list", length(relation$partitions)),
    relation$partitions)
  diagnostics <- list(
    max_relation_block_bytes = .dense_double_matrix_bytes(
      q, largest_feature_block
    ),
    max_factor_slice_bytes = .dense_double_matrix_bytes(
      nrow(operator), largest_feature_block
    ),
    projected_state_bytes = 0,
    blocks_read = 0L
  )
  for (partition in relation$partitions) {
    value <- matrix(0, q, nrow(operator))
    for (start in seq.int(1L, relation$n_features, by = feature_block)) {
      features <- start:min(start + feature_block - 1L, relation$n_features)
      relation_block <- session$read(side, partition, features)
      factor_slice <- operator[, features, drop = FALSE]
      value <- value + relation_block %*% t(factor_slice)
      diagnostics$blocks_read <- diagnostics$blocks_read + 1L
    }
    projected[[partition]] <- lapply(ranges, function(columns) {
      value[, columns, drop = FALSE]
    })
    diagnostics$projected_state_bytes <-
      diagnostics$projected_state_bytes +
      as.double(utils::object.size(projected[[partition]]))
  }
  list(values = projected, diagnostics = diagnostics)
}

.measurement_empty_blocks <- function(task) {
  edges <- task$neural_boundary$edges$edges
  left_frame <- task$neural_boundary$left_frame
  right_frame <- task$neural_boundary$right_frame
  stats::setNames(lapply(seq_len(nrow(edges)), function(index) {
    matrix(0,
      nrow(left_frame$legs[[edges$left[[index]]]]$operator),
      nrow(right_frame$legs[[edges$right[[index]]]]$operator)
    )
  }), paste0("edge_", sprintf("%06d", seq_len(nrow(edges)))))
}

.measurement_contract_projected <- function(task, left, right, route_plan,
                                            state) {
  h <- as.matrix(task$experimental_boundary$query$operator)
  partition_edges <- task$ordered_partition_products
  measurement_edges <- task$neural_boundary$edges$edges
  blocks <- .measurement_empty_blocks(task)
  route <- route_plan$route
  factorization <- if (route == "factorized_h") {
    .measurement_h_factorization(h, route_plan$rank_tolerance)
  } else {
    NULL
  }

  if (route == "scalar_stack") {
    left_nodes <- unique(measurement_edges$left)
    right_nodes <- unique(measurement_edges$right)
    left_index <- match(measurement_edges$left, left_nodes)
    right_index <- match(measurement_edges$right, right_nodes)
    values <- numeric(nrow(measurement_edges))
    for (partition_edge in seq_len(nrow(partition_edges))) {
      z_left <- do.call(cbind,
        left[[partition_edges$left[[partition_edge]]]][left_nodes])
      z_right <- do.call(cbind,
        right[[partition_edges$right[[partition_edge]]]][right_nodes])
      hz_right <- h %*% z_right
      edge_values <- colSums(
        z_left[, left_index, drop = FALSE] *
          hz_right[, right_index, drop = FALSE]
      )
      values <- values + partition_edges$weight[[partition_edge]] * edge_values
      state$completed_partition_edges <- partition_edge
    }
    for (edge in seq_along(blocks)) blocks[[edge]][1, 1] <- values[[edge]]
    return(blocks)
  }

  for (partition_edge in seq_len(nrow(partition_edges))) {
    left_partition <- left[[partition_edges$left[[partition_edge]]]]
    right_partition <- right[[partition_edges$right[[partition_edge]]]]
    weight <- partition_edges$weight[[partition_edge]]
    for (tile_start in seq.int(1L, nrow(measurement_edges),
      by = route_plan$edge_tile)) {
      tile <- tile_start:min(
        tile_start + route_plan$edge_tile - 1L,
        nrow(measurement_edges)
      )
      for (edge in tile) {
        z_left <- left_partition[[measurement_edges$left[[edge]]]]
        z_right <- right_partition[[measurement_edges$right[[edge]]]]
        value <- switch(route,
          forward_k = {
            out <- matrix(0, ncol(z_left), ncol(z_right))
            for (column in seq_len(ncol(z_right))) {
              h_right <- h %*% z_right[, column]
              for (row in seq_len(ncol(z_left))) {
                out[row, column] <- sum(z_left[, row] * h_right)
              }
            }
            out
          },
          pull_h = crossprod(z_left, h %*% z_right),
          factorized_h = {
            left_factor <- crossprod(factorization$left, z_left)
            right_factor <- crossprod(factorization$right, z_right)
            crossprod(left_factor, right_factor)
          },
          multivariate_blocks = crossprod(z_left, h %*% z_right),
          .invariant_error("Unknown measurement contraction route.")
        )
        blocks[[edge]] <- blocks[[edge]] + weight * value
      }
    }
    state$completed_partition_edges <- partition_edge
  }
  blocks
}

.measurement_source_bytes <- function(summary) {
  if (is.null(summary)) return(0)
  if (!is.null(summary$bytes_read)) return(sum(summary$bytes_read))
  sum(vapply(summary, .measurement_source_bytes, numeric(1)))
}

.measurement_failure_receipt <- function(task, route_plan, memory, error,
                                         completed_partition_edges,
                                         source_summary = NULL) {
  structure(list(
    scientific_plan_id = task$task_id,
    completion_status = "failed",
    route = if (is.null(route_plan)) NULL else route_plan$route,
    route_signature = if (is.null(route_plan)) NULL else route_plan$signature,
    memory = if (is.null(memory)) NULL else memory$plan,
    completed_partition_edges = as.integer(completed_partition_edges),
    source_revisions = task$source_uses[, c("side", "relation_id", "partition",
      "stable_revision"), drop = FALSE],
    source_summary = source_summary,
    error_class = class(error),
    message = conditionMessage(error)
  ), class = "effect_measurement_failure_receipt")
}

.measurement_kernel_condition <- function(receipt) {
  structure(list(
    message = paste0("Measurement contraction failed: ", receipt$message),
    call = NULL,
    receipt = receipt
  ), class = c("effect_measurement_kernel_error", "error", "condition"))
}

.run_measurement_contraction <- function(
    task, compute = compute_policy(),
    route = c("auto", "forward_k", "pull_h", "factorized_h",
              "scalar_stack", "multivariate_blocks"),
    edge_tile = 128L, rank_tolerance = 1e-12,
    storage = c("memory", "block"),
    open_descriptor = .open_source_descriptor) {
  task <- .validate_raw_measurement_task(task)
  storage <- match.arg(storage)
  compute_plan <- .execution_preflight(compute, function() {
    c(task$left_relation$capabilities,
      if (task$same_relation) list() else task$right_relation$capabilities)
  })
  feature_block <- compute_plan$compute$block_features
  if (is.null(feature_block)) {
    feature_block <- min(1024L,
      task$left_relation$n_features, task$right_relation$n_features)
  }
  route_plan <- .measurement_route_plan(
    task, route, feature_block, edge_tile, rank_tolerance
  )
  memory <- .measurement_kernel_memory_plan(
    task, route_plan, storage, compute_plan$compute$workspace_bytes
  )
  state <- new.env(parent = emptyenv())
  state$completed_partition_edges <- 0L
  session <- NULL
  source_summary <- NULL
  started <- proc.time()[["elapsed"]]

  result <- tryCatch({
    if (!is.na(memory$plan$fits_budget) && !memory$plan$fits_budget) {
      .input_error(
        "Measurement contraction exceeds the declared workspace budget."
      )
    }
    session <- .open_evidence_task_source_session(
      task, open_descriptor = open_descriptor
    )
    left_nodes <- .measurement_requested_nodes(task, "left")
    right_nodes <- .measurement_requested_nodes(task, "right")
    projected_left <- .project_measurement_family(
      task, session, "left", left_nodes, feature_block
    )
    projected_right <- if (task$same_relation &&
        identical(task$neural_boundary$left_frame$signature,
          task$neural_boundary$right_frame$signature) &&
        identical(left_nodes, right_nodes)) {
      projected_left
    } else {
      .project_measurement_family(
        task, session, "right", right_nodes, feature_block
      )
    }
    blocks <- .measurement_contract_projected(
      task, projected_left$values, projected_right$values,
      route_plan, state
    )
    source_summary <- session$summary()
    session$close()
    session <- NULL
    elapsed <- proc.time()[["elapsed"]] - started
    diagnostics <- list(
      route = route_plan$route,
      feature_block = feature_block,
      edge_tile = route_plan$edge_tile,
      partition_edges_completed = state$completed_partition_edges,
      measurement_edges = length(blocks),
      max_relation_block_bytes = max(
        projected_left$diagnostics$max_relation_block_bytes,
        projected_right$diagnostics$max_relation_block_bytes
      ),
      max_factor_slice_bytes = max(
        projected_left$diagnostics$max_factor_slice_bytes,
        projected_right$diagnostics$max_factor_slice_bytes
      ),
      projected_state_bytes =
        projected_left$diagnostics$projected_state_bytes +
        if (identical(projected_left, projected_right)) 0 else
          projected_right$diagnostics$projected_state_bytes,
      source_bytes_read = .measurement_source_bytes(source_summary),
      elapsed_seconds = elapsed,
      global_neural_operator_materialized = FALSE,
      kronecker_lift_materialized = FALSE
    )
    structure(list(
      blocks = blocks,
      edges = task$neural_boundary$edges,
      task_id = task$task_id,
      compute = compute_plan$compute,
      route_plan = route_plan,
      memory = memory,
      source_summary = source_summary,
      diagnostics = diagnostics
    ), class = "effect_measurement_contraction")
  }, error = function(error) {
    if (!is.null(session)) {
      source_summary <<- tryCatch(session$summary(), error = function(...) NULL)
      try(session$close(), silent = TRUE)
      session <<- NULL
    }
    receipt <- .measurement_failure_receipt(
      task, route_plan, memory, error,
      state$completed_partition_edges, source_summary
    )
    stop(.measurement_kernel_condition(receipt))
  })
  result
}

.measurement_form_dense_reference <- function(task, left_relations,
                                              right_relations = left_relations) {
  task <- .validate_raw_measurement_task(task)
  h <- as.matrix(task$experimental_boundary$query$operator)
  partition_edges <- task$ordered_partition_products
  measurement_edges <- task$neural_boundary$edges$edges
  blocks <- .measurement_empty_blocks(task)
  for (partition_edge in seq_len(nrow(partition_edges))) {
    left <- left_relations[[partition_edges$left[[partition_edge]]]]
    right <- right_relations[[partition_edges$right[[partition_edge]]]]
    if (!is.matrix(left) || !is.matrix(right)) {
      .input_error("Dense reference relations must be named matrix families.")
    }
    for (edge in seq_len(nrow(measurement_edges))) {
      left_leg <- task$neural_boundary$left_frame$legs[[
        measurement_edges$left[[edge]]
      ]]$operator
      right_leg <- task$neural_boundary$right_frame$legs[[
        measurement_edges$right[[edge]]
      ]]$operator
      value <- left_leg %*% t(left) %*% h %*% right %*% t(right_leg)
      blocks[[edge]] <- blocks[[edge]] +
        partition_edges$weight[[partition_edge]] * value
    }
  }
  blocks
}

.measurement_benchmark_record <- function(task, iterations = 3L,
                                          feature_block = 256L,
                                          edge_tile = 64L) {
  iterations <- .validate_tile_size(iterations, "iterations")
  routes <- c("auto", "pull_h", "forward_k", "multivariate_blocks")
  if (.measurement_route_plan(task, "auto", feature_block,
      edge_tile)$eligible[["scalar_stack"]]) {
    routes <- c(routes, "scalar_stack")
  }
  factor_plan <- .measurement_route_plan(task, "auto", feature_block, edge_tile)
  if (factor_plan$eligible[["factorized_h"]]) routes <- c(routes, "factorized_h")
  records <- lapply(unique(routes), function(route) {
    elapsed <- numeric(iterations)
    output_bytes <- numeric(iterations)
    for (iteration in seq_len(iterations)) {
      started <- proc.time()[["elapsed"]]
      run <- .run_measurement_contraction(
        task,
        compute = compute_policy(block_features = feature_block),
        route = route,
        edge_tile = edge_tile
      )
      elapsed[[iteration]] <- proc.time()[["elapsed"]] - started
      output_bytes[[iteration]] <- as.double(utils::object.size(run$blocks))
    }
    data.frame(
      route = route,
      iterations = iterations,
      median_elapsed_seconds = stats::median(elapsed),
      output_bytes = max(output_bytes),
      edge_count = nrow(task$neural_boundary$edges$edges),
      scalar_frame = all(vapply(c(
        task$neural_boundary$left_frame$legs,
        task$neural_boundary$right_frame$legs
      ), function(x) nrow(x$operator) == 1L, logical(1))),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, records)
}
