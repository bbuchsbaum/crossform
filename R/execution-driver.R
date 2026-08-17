# Geometry execution driver --------------------------------------------------
#
# Layer 4 (compiler / execution). This file is the runtime: it consumes an
# already-compiled `effect_geometry_execution_plan` and performs the reads,
# contractions, receipts, and cleanup that produce a result. Nothing here
# decides *what* to compute; `R/compiler.R` lowers a plan and hands it over.
#
# The executor used to be one 324-line function in `R/compiler.R` with seven
# inline closures and seven levels of nesting. It is now a sequence of named
# stages. Evaluation order, timing records, receipt fields, reporter events,
# and the guarded source cleanup are unchanged; only the nesting is gone.

.execution_measurement_index <- function(at) {
  index <- at$index
  if (is.null(index)) return(seq_len(nrow(at$weights)))
  if (is.data.frame(index) && "measurement" %in% names(index)) {
    return(index$measurement)
  }
  index
}

.execution_blas_record <- function() {
  vendor <- tryCatch(unname(extSoftVersion()[["BLAS"]]), error = function(e) NULL)
  if (is.null(vendor) || !.is_string(vendor)) vendor <- "unknown"
  list(vendor = vendor, requested_threads = 1L,
    observed_threads = NA_integer_)
}

.planned_execution_receipt <- function(plan, validate = TRUE) {
  if (isTRUE(validate)) .validate_geometry_execution_plan(plan)
  task <- plan$task
  sources <- task$left_relation$capabilities
  if (!isTRUE(task$same_relation)) {
    sources <- c(sources, task$right_relation$capabilities)
  }
  execution_receipt(
    scientific_plan_id = plan$scientific_plan_id,
    compute = plan$compute,
    sources = sources,
    memory = plan$memory,
    kernel_version = plan$kernel_version,
    task_partition_id = sprintf(
      "ascending-features-%d", plan$feature_block
    ),
    reduction_plan_id = task$stages$signature,
    numeric_contract = numerical_contract(),
    completion_status = "planned",
    task_count = plan$task_count,
    completed_task_count = 0L,
    blas = .execution_blas_record(),
    domain_signature = task$left_relation$domain$signature
  )
}

.execution_reporter <- function(reporter, final_receipt) {
  force(reporter)
  function(event) {
    if (event$type %in% c("complete", "failed", "interrupted")) {
      final_receipt$value <- event$receipt
    }
    if (!is.null(reporter)) reporter(event)
    invisible(NULL)
  }
}

.open_execution_storage <- function(storage, storage_path, dim) {
  storage <- match.arg(storage, c("memory", "block"))
  if (storage == "memory") {
    if (!is.null(storage_path)) {
      .input_error("`storage_path` is only valid for block-backed geometry.")
    }
    return(list(kind = storage, total = NULL, coherent = NULL,
      created = character(), created_directory = FALSE))
  }
  if (!.is_string(storage_path)) {
    .input_error("Block-backed geometry requires one nonempty `storage_path`.")
  }
  created_directory <- FALSE
  if (!dir.exists(storage_path)) {
    parent <- dirname(storage_path)
    if (!dir.exists(parent)) {
      .input_error("The parent of `storage_path` must already exist.")
    }
    if (!dir.create(storage_path, recursive = FALSE, showWarnings = FALSE)) {
      .input_error("Could not create `storage_path`.")
    }
    created_directory <- TRUE
  }
  total_path <- file.path(storage_path, "total.egm")
  coherent_path <- file.path(storage_path, "coherent.egm")
  if (file.exists(total_path) || file.exists(coherent_path)) {
    .input_error("Refusing to overwrite an existing geometry component.")
  }
  total <- .file_geometry_store(total_path, dim, create = TRUE)
  coherent <- tryCatch(
    .file_geometry_store(coherent_path, dim, create = TRUE),
    error = function(error) {
      unlink(total_path)
      if (created_directory) unlink(storage_path, recursive = TRUE)
      stop(error)
    }
  )
  list(kind = storage, total = total, coherent = coherent,
    created = c(total_path, coherent_path),
    created_directory = created_directory)
}

.execution_storage_cleanup <- function(storage, completed) {
  force(storage)
  force(completed)
  function() {
    if (!isTRUE(completed$value) && length(storage$created)) {
      unlink(storage$created)
      remaining <- storage$created[file.exists(storage$created)]
      if (length(remaining)) {
        .input_error("Failed to remove incomplete geometry components.")
      }
      if (storage$created_directory) {
        directory <- dirname(storage$created[[1L]])
        unlink(directory, recursive = TRUE)
        if (dir.exists(directory)) {
          .input_error("Failed to remove incomplete geometry directory.")
        }
      }
    }
    invisible(NULL)
  }
}

.ordered_pairing_marginals <- function(local, ordered_edges, mass) {
  if (!is.array(local) || length(dim(local)) != 3L ||
      !.is_finite_numeric(local)) {
    .invariant_error(
      "Compiled marginals require finite local relation first moments."
    )
  }
  partitions <- dimnames(local)[[3L]]
  effects <- dimnames(local)[[2L]]
  .validate_ordered_partition_edges(
    ordered_edges, partitions, partitions, TRUE
  )
  measurements <- dim(local)[[1L]]
  if (!.is_finite_numeric(mass) || !(length(mass) %in% c(1L, measurements)) ||
      any(mass <= 0)) {
    .invariant_error("Compiled marginal mass must be positive and finite.")
  }
  mass <- rep_len(mass, measurements)
  normalized <- function(partition) {
    matrix(local[, , match(partition, partitions), drop = FALSE],
      nrow = measurements, ncol = length(effects)) / mass
  }
  accumulate <- function(endpoint) {
    value <- matrix(0, measurements, length(effects))
    for (edge in seq_len(nrow(ordered_edges))) {
      value <- value + ordered_edges$weight[[edge]] *
        normalized(ordered_edges[[endpoint]][[edge]])
    }
    dimnames(value) <- dimnames(local)[1:2]
    value
  }
  if (identical(attr(ordered_edges, "expansion"),
      "self_adjoint_half_edges")) {
    value <- list(endpoint = accumulate("left"))
    semantics <- "undirected_endpoint"
  } else {
    value <- list(left = accumulate("left"), right = accumulate("right"))
    semantics <- "directed_roles"
  }
  attr(value, "semantics") <- semantics
  class(value) <- c("effect_marginals", "list")
  value
}

# Observation and event recording --------------------------------------------

.empty_geometry_execution_observations <- function(plan) {
  observed <- .empty_execution_observations()
  observed$task_counts[["planned"]] <- plan$task_count
  observed$tiles <- list(
    feature_block = plan$feature_block,
    row_tile = plan$row_tile,
    coordinate_tile = plan$coordinate_tile
  )
  observed
}

.record_execution_stage <- function(observed_state, stage, started) {
  observed_state$value$stage_seconds[[stage]] <-
    max(0, proc.time()[["elapsed"]] - started)
  invisible(NULL)
}

.report_execution_event <- function(observed_state) {
  force(observed_state)
  function(event, features) {
    observed_state$value$task_counts[[event]] <-
      observed_state$value$task_counts[[event]] + 1
    if (identical(event, "completed")) {
      observed_state$value$features_completed <-
        observed_state$value$features_completed + length(features)
    }
    invisible(NULL)
  }
}

# Sources --------------------------------------------------------------------

# The handle is returned with its own idempotent `close()` so the caller can
# both register it with `on.exit()` and close it explicitly at the end of a
# successful run without double-closing the session.
.open_execution_sources <- function(task) {
  source_session <- .open_effect_task_source_session(task, validate = FALSE)
  closed <- FALSE
  list(
    session = source_session,
    close = function() {
      if (!closed) {
        closed <<- TRUE
        source_session$close()
      }
      invisible(NULL)
    }
  )
}

# Two-sided sessions report per-side summaries; flatten them with
# side-qualified names so receipt observations keep one canonical shape on
# both self-form and rectangular executions.
.execution_session_summary <- function(source_session) {
  value <- source_session$summary()
  if (!is.null(value$access_mode)) return(value)
  side_named <- function(side, field) {
    entries <- value[[side]][[field]]
    stats::setNames(entries, paste0(side, ":", names(entries)))
  }
  list(
    access_mode = c(
      side_named("left", "access_mode"),
      side_named("right", "access_mode")
    ),
    bytes_read = c(
      side_named("left", "bytes_read"),
      side_named("right", "bytes_read")
    ),
    sides = value
  )
}

.execution_relation_reader <- function(relation, side, source_session,
                                       observed_state) {
  force(relation)
  force(side)
  force(source_session)
  force(observed_state)
  function(partition, features) {
    value <- .relation_block_with_reader(
      relation, partition, features,
      function(partition, features) {
        source_session$read(side, partition, features)
      },
      validate = FALSE
    )
    observed_state$value$bytes_read <-
      sum(.execution_session_summary(source_session)$bytes_read)
    value
  }
}

.execution_relation_readers <- function(task, source_session, observed_state) {
  left <- .execution_relation_reader(
    task$left_relation, "left", source_session, observed_state
  )
  list(
    left = left,
    right = if (isTRUE(task$same_relation)) {
      left
    } else {
      .execution_relation_reader(
        task$right_relation, "right", source_session, observed_state
      )
    }
  )
}

# Block-backed accumulation --------------------------------------------------

.execution_total_accumulator <- function(geometry_storage, storage) {
  if (!identical(storage, "block")) return(NULL)
  force(geometry_storage)
  function(rows, coordinates, increment) {
    existing <- .read_geometry_tile(geometry_storage$total, rows, coordinates)
    .write_geometry_tile(geometry_storage$total, rows, coordinates,
      existing + increment)
  }
}

.execution_coherent_writer <- function(geometry_storage, storage) {
  if (!identical(storage, "block")) return(NULL)
  force(geometry_storage)
  function(rows, coordinates, value) {
    .write_geometry_tile(geometry_storage$coherent, rows, coordinates, value)
  }
}

# Contraction stages ---------------------------------------------------------

.execution_support_streamed <- function(plan) {
  grepl("^support_streamed_metric_", plan$lowering)
}

.execute_node_block <- function(plan, task, readers, accumulate_tile,
                                task_observer) {
  at <- plan$frame
  requirements <- plan$requirements
  if (.execution_support_streamed(plan)) {
    return(.support_streamed_metric_contraction(
      frame = at,
      metric_schedule = plan$metric_schedule,
      read_relation = readers$left,
      partitions = task$left_relation$partitions,
      effects = task$left_space$coordinates,
      ordered_edges = task$ordered_edges,
      query = plan$query,
      form_total = requirements$total,
      form_coherent = requirements$coherent,
      retain_first_moments = requirements$marginals,
      task_observer = task_observer
    ))
  }
  kernel_frame <- if (grepl("^additive_metric_", plan$lowering)) {
    .metric_additive_frame(at, plan$metric_schedule)
  } else {
    at
  }
  same_relation <- isTRUE(task$same_relation)
  .streamed_effect_form_contraction(
    frame = kernel_frame,
    read_left = readers$left,
    read_right = readers$right,
    left_partitions = task$left_relation$partitions,
    right_partitions = task$right_relation$partitions,
    left_effects = task$left_space$coordinates,
    right_effects = task$right_space$coordinates,
    ordered_edges = task$ordered_edges,
    codec = if (same_relation) "symmetric_packed" else "rectangular",
    same_relation = same_relation,
    feature_block = plan$feature_block,
    row_tile = plan$row_tile,
    coordinate_tile = plan$coordinate_tile,
    accumulate_tile = accumulate_tile,
    retain_first_moments = requirements$coherent,
    query = plan$query,
    form_total = requirements$total,
    task_observer = task_observer
  )
}

.execute_coherent_block <- function(plan, task, streamed, write_tile) {
  if (!plan$requirements$coherent) return(NULL)
  if (.execution_support_streamed(plan)) {
    return(list(value = streamed$coherent, diagnostics = list(
      lowering = "metric_aware_rank_one_component",
      factorization_count = streamed$diagnostics$metric_factorizations,
      durable_output_bytes = 8 * length(streamed$coherent)
    )))
  }
  same_relation <- isTRUE(task$same_relation)
  .effect_form_coherent_from_first_moments(
    streamed$first_moments$left,
    streamed$first_moments$right,
    task$ordered_edges,
    streamed$mass,
    codec = if (same_relation) "symmetric_packed" else "rectangular",
    same_relation = same_relation,
    row_tile = plan$row_tile,
    write_tile = write_tile,
    query = plan$query
  )
}

.execute_marginal_block <- function(plan, task, streamed) {
  if (!plan$requirements$marginals) return(NULL)
  .ordered_pairing_marginals(
    streamed$first_moments$left,
    task$ordered_edges,
    streamed$mass
  )
}

.execution_component_values <- function(plan, streamed, coherent,
                                        geometry_storage) {
  requirements <- plan$requirements
  in_memory <- identical(plan$storage, "memory")
  list(
    total = if (!requirements$total) {
      NULL
    } else if (in_memory) {
      streamed$value
    } else {
      geometry_storage$total
    },
    coherent = if (!requirements$coherent) {
      NULL
    } else if (in_memory) {
      coherent$value
    } else {
      geometry_storage$coherent
    }
  )
}

# Result assembly ------------------------------------------------------------

.execution_total_diagnostics <- function(streamed) {
  diagnostics <- streamed$diagnostics
  first_moment_bytes <- diagnostics$durable_left_first_moment_bytes
  if (is.null(first_moment_bytes)) {
    first_moment_bytes <- diagnostics$durable_first_moment_bytes
  }
  if (!is.null(first_moment_bytes)) {
    diagnostics$durable_local_relation_bytes <- first_moment_bytes
  }
  diagnostics
}

.execution_metadata <- function(plan, task, streamed, coherent) {
  at <- plan$frame
  requirements <- plan$requirements
  list(
    frame = list(representation = at$representation,
      normalization = at$normalization, domain = at$domain),
    metric_schedule = list(
      kind = plan$metric_schedule$kind,
      signature = plan$metric_schedule$signature,
      feature_additive = plan$metric_schedule$feature_additive,
      support_dense = plan$metric_schedule$support_dense,
      scope = plan$metric_schedule$scope,
      frame_composition = plan$metric_schedule$frame_composition,
      materialization = plan$metric_schedule$materialization
    ),
    pairing_estimate = attr(task$ordered_edges, "source_estimate"),
    edge_operation = list(
      normalizer = task$normalizer,
      transform = task$transform,
      reducer = task$reducer,
      lowering = task$lowering,
      order = plan$task$stages$order
    ),
    measurement_bridge = list(
      kind = task$bridge$kind,
      signature = task$bridge$signature
    ),
    storage = plan$storage,
    tiles = list(
      feature_block = plan$feature_block,
      row_tile = plan$row_tile,
      coordinate_tile = plan$coordinate_tile
    ),
    requirements = requirements,
    diagnostics = list(
      total = if (requirements$total) {
        .execution_total_diagnostics(streamed)
      } else {
        NULL
      },
      coherent = if (requirements$coherent) coherent$diagnostics else NULL
    ),
    execution_plan = list(
      signature = plan$signature,
      parent = plan$parent_signature,
      lowering = plan$lowering,
      query_fused = plan$query_fused,
      kernel_version = plan$kernel_version,
      stage_signature = plan$task$stages$signature,
      materialization = requirements$materialization
    ),
    scientific_plan_id = plan$scientific_plan_id
  )
}

.execution_geometry_result <- function(plan, task, values, marginals, metadata,
                                       planned_receipt) {
  query <- plan$query
  component <- plan$component
  requirements <- plan$requirements
  same_relation <- isTRUE(task$same_relation)
  index <- .execution_measurement_index(plan$frame)
  if (is.null(query) && !same_relation) {
    return(effect_form(
      total = values$total,
      coherent = values$coherent,
      marginals = NULL,
      left_space = task$left_space,
      right_space = task$right_space,
      receipt = planned_receipt,
      index = index,
      metadata = metadata,
      codec = "rectangular",
      symmetric = FALSE,
      guaranteed_psd = FALSE
    ))
  }
  if (is.null(query)) {
    return(effect_geometry(values$total, values$coherent, marginals,
      effects = task$left_space, receipt = planned_receipt,
      index = index, metadata = metadata))
  }
  if (identical(component, "contrast")) {
    return(.new_effect_contrast_view(
      values$total, values$coherent, marginals, plan$signed_query,
      index, planned_receipt
    ))
  }
  total_matrix <- if (requirements$total) values$total else NULL
  coherent_matrix <- if (requirements$coherent) values$coherent else NULL
  view_values <- switch(component,
    total = total_matrix,
    coherent = coherent_matrix,
    configuration = total_matrix - coherent_matrix
  )
  colnames(view_values) <- .query_output_labels(query)
  effect_view(view_values, query, component, planned_receipt,
    effects = if (same_relation) task$left_space else NULL,
    left_space = if (same_relation) NULL else task$left_space,
    right_space = if (same_relation) NULL else task$right_space,
    index = index, metadata = metadata)
}

# The guarded body -----------------------------------------------------------

# `on.exit()` is registered here, in the one frame that owns the source
# session for the whole computation: the session closes exactly once, whether
# the run completes, errors, or is interrupted, and the explicit close before
# the result is assembled stays idempotent.
.execute_geometry_compute <- function(plan, task, geometry_storage,
                                      planned_receipt, observed_state,
                                      completed) {
  storage <- plan$storage
  admission_started <- proc.time()[["elapsed"]]
  sources <- .open_execution_sources(task)
  on.exit(sources$close(), add = TRUE)
  .record_execution_stage(observed_state, "source_admission", admission_started)
  observed_state$value$source_access <-
    .execution_session_summary(sources$session)$access_mode
  readers <- .execution_relation_readers(task, sources$session, observed_state)
  kernel_started <- proc.time()[["elapsed"]]
  streamed <- .execute_node_block(
    plan, task, readers,
    .execution_total_accumulator(geometry_storage, storage),
    .report_execution_event(observed_state)
  )
  .record_execution_stage(observed_state, "feature_tasks", kernel_started)
  coherent_started <- proc.time()[["elapsed"]]
  coherent <- .execute_coherent_block(
    plan, task, streamed,
    .execution_coherent_writer(geometry_storage, storage)
  )
  .record_execution_stage(observed_state, "coherent", coherent_started)
  values <- .execution_component_values(
    plan, streamed, coherent, geometry_storage
  )
  marginal_started <- proc.time()[["elapsed"]]
  marginals <- .execute_marginal_block(plan, task, streamed)
  .record_execution_stage(observed_state, "marginals", marginal_started)
  metadata <- .execution_metadata(plan, task, streamed, coherent)
  sources$close()
  source_summary <- .execution_session_summary(sources$session)
  metadata$source_session <- source_summary
  observed_state$value$bytes_read <- sum(source_summary$bytes_read)
  result <- .execution_geometry_result(
    plan, task, values, marginals, metadata, planned_receipt
  )
  completed$value <- TRUE
  result
}

.execute_geometry_plan <- function(plan, reporter = NULL) {
  .validate_geometry_execution_plan(plan)
  task <- .as_compiled_effect_task(plan$task, validate = FALSE)
  planned_receipt <- .planned_execution_receipt(plan, validate = FALSE)
  geometry_storage <- .open_execution_storage(
    plan$storage, plan$storage_path,
    c(nrow(plan$frame$weights), plan$output_width)
  )
  completed <- new.env(parent = emptyenv())
  completed$value <- FALSE
  final_receipt <- new.env(parent = emptyenv())
  final_receipt$value <- planned_receipt
  observed_state <- new.env(parent = emptyenv())
  observed_state$value <- .empty_geometry_execution_observations(plan)

  value <- .execute_guarded(
    compute = function() {
      .execute_geometry_compute(plan, task, geometry_storage, planned_receipt,
        observed_state, completed)
    },
    receipt = planned_receipt,
    reporter = .execution_reporter(reporter, final_receipt),
    cleanup = .execution_storage_cleanup(geometry_storage, completed),
    observations = function() observed_state$value,
    receipt_sink = function(receipt) final_receipt$value <- receipt
  )
  value$receipt <- final_receipt$value
  value
}

# The one entry that composes plan -> compile -> execute.
.run_geometry_compiler <- function(x, at = NULL, over = NULL,
                                   compute = NULL, storage = "memory",
                                   storage_path = NULL, query = NULL,
                                   component = NULL, reporter = NULL,
                                   signed_query = NULL) {
  if (inherits(x, "effect_geometry_plan")) {
    if (!is.null(at) || !is.null(over) || !is.null(compute)) {
      .input_error(
        "A compiled geometry plan cannot be combined with raw plan inputs."
      )
    }
    plan <- x
  } else {
    if (is.null(at) || is.null(over)) {
      .input_error(
        "A relation requires both `at` and `over` to compile geometry."
      )
    }
    if (is.null(compute)) compute <- compute_policy()
    plan <- plan_geometry(x, at, over, compute)
  }
  execution <- .compile_geometry_execution_plan(
    plan, query = query, component = component,
    storage = storage, storage_path = storage_path,
    signed_query = signed_query
  )
  .execute_geometry_plan(execution, reporter = reporter)
}
