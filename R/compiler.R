# Public relation-to-geometry compiler -------------------------------------

.compiler_capabilities <- function(x) {
  if (is.null(x$capabilities)) {
    stop(paste0(
      "Opaque relation sources require explicit `source_capabilities()` ",
      "before execution."
    ), call. = FALSE)
  }
  capabilities <- lapply(x$capabilities, .validate_source_capabilities)
  if (!all(vapply(capabilities, function(value) isTRUE(value$block_read),
    logical(1)))) {
    stop("All relation sources must support bounded block reads.", call. = FALSE)
  }
  capabilities
}

.effect_task_source_uses <- function(left, right, left_id, right_id) {
  side_uses <- function(relation, side, relation_id) {
    do.call(rbind, lapply(relation$partitions, function(partition) {
      descriptor <- relation$sources[[partition]]$descriptor
      handle_key <- if (is.null(descriptor) ||
          identical(descriptor$access, "coordinator")) {
        NA_character_
      } else {
        .source_descriptor_key(descriptor)
      }
      data.frame(
        side = side,
        relation_id = relation_id,
        partition = partition,
        stable_revision = relation$capabilities[[partition]]$stable_revision,
        handle_key = handle_key,
        stringsAsFactors = FALSE
      )
    }))
  }
  uses <- rbind(
    side_uses(left, "left", left_id),
    side_uses(right, "right", right_id)
  )
  rownames(uses) <- NULL
  list(
    uses = uses,
    distinct_handle_keys = unique(stats::na.omit(uses$handle_key))
  )
}

.effect_task_semantic <- function(left_id, right_id, left_space, right_space,
                                  edges, bridge, normalizer, transform, reducer,
                                  query_identity) {
  list(
    schema_version = 1L,
    left_relation = left_id,
    right_relation = right_id,
    left_space = left_space,
    right_space = right_space,
    ordered_edges = unclass(as.data.frame(edges)),
    edge_expansion = attr(edges, "expansion"),
    bridge = bridge$signature,
    normalizer = unclass(normalizer),
    transform = unclass(transform),
    reducer = unclass(reducer),
    query = query_identity
  )
}

.effect_task_id <- function(semantic) {
  paste0(
    "sha256:", digest::digest(semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L)
  )
}

# The identity the task would have carried without an embedded query. This is
# the estimand-bearing task id: a fused query is an execution strategy, so a
# queried task must still be verifiable against the complete-form identity of
# its parent plan.
.effect_task_base_id <- function(task) {
  if (identical(task$materialization$completeness, "complete_form")) {
    return(task$task_id)
  }
  semantic <- .effect_task_semantic(
    task$left_relation_id, task$right_relation_id,
    task$spaces$experimental_left, task$spaces$experimental_right,
    task$ordered_partition_products,
    task$neural_boundary$bridge,
    structure(task$stages$normalization$operation,
      class = "effect_edge_normalizer"),
    structure(task$stages$transform$operation,
      class = "effect_edge_transform"),
    structure(task$stages$reduction$operation,
      class = "effect_partition_reducer"),
    list(kind = "complete_form", query = NULL)
  )
  .effect_task_id(semantic)
}

.compile_effect_evidence_task <- function(
    left, over, right = NULL, query = NULL,
    normalizer = inner_product(), transform = NULL,
    reducer = .partition_reducer(), bridge = NULL) {
  .validate_relation(left)
  same_relation <- is.null(right)
  if (same_relation) right <- left else .validate_relation(right)
  left$capabilities <- .compiler_capabilities(left)
  right$capabilities <- if (same_relation) left$capabilities else
    .compiler_capabilities(right)
  normalizer <- .validate_edge_normalizer(normalizer)
  if (is.null(transform)) transform <- .identity_edge_transform()
  transform <- .validate_edge_transform(transform)
  reducer <- .validate_partition_reducer(reducer)
  operation <- .edge_operation_plan(normalizer, transform, reducer)
  bridge <- if (is.null(bridge)) {
    .identity_measurement_bridge(left, right)
  } else {
    .validate_measurement_bridge(bridge, left, right)
  }
  edges <- .ordered_partition_edges(
    over, left$partitions, right$partitions, same_relation
  )
  left_id <- .relation_family_identity(left)
  right_id <- if (same_relation) left_id else .relation_family_identity(right)

  query_identity <- if (is.null(query)) {
    list(kind = "complete_form", query = NULL)
  } else if (.is_pair_difference_query(query)) {
    .validate_pair_difference_for_task(
      query, left$effect_space$coordinates, right$effect_space$coordinates,
      same_relation
    )
    list(kind = "pair_difference_self_query", query = query)
  } else if (inherits(query, "effect_pair_query")) {
    .validate_query_for_compile(query)
    if (!.same_effect_space(query$left_space, left$effect_space) ||
        !.same_effect_space(query$right_space, right$effect_space)) {
      stop("Pair-query axes must match the ordered relation effect spaces.",
        call. = FALSE)
    }
    list(kind = "pair_query", query = query)
  } else if (inherits(query, "effect_query")) {
    .validate_query_for_compile(query)
    if (!identical(query$kind, "bilinear") || !isTRUE(query$fixed) ||
        !.same_effect_space(left$effect_space, right$effect_space) ||
        (!is.null(query$effect_space) &&
          !.same_effect_space(query$effect_space, left$effect_space))) {
      stop("A bilinear query requires matching self-form relation axes.",
        call. = FALSE)
    }
    list(kind = "bilinear_compatibility", query = query)
  } else if (is.matrix(query) && is.numeric(query) &&
      all(is.finite(query)) && nrow(query) > 0L && ncol(query) > 0L) {
    if (!same_relation) {
      stop("Raw physical queries are only compatible with self-form tasks.",
        call. = FALSE)
    }
    q <- length(left$effect_space$coordinates)
    if (nrow(query) != q * (q + 1L) / 2L) {
      stop("A raw self-form query must match the symmetric-packed width.",
        call. = FALSE)
    }
    list(kind = "physical_self_query", query = query)
  } else {
    stop("Task queries must be NULL, axis-bound queries, or finite matrices.",
      call. = FALSE)
  }
  semantic <- .effect_task_semantic(
    left_id, right_id, left$effect_space, right$effect_space,
    edges, bridge, normalizer, transform, reducer, query_identity
  )
  .new_effect_evidence_task(
    left, right, same_relation, edges, bridge, normalizer, transform, reducer,
    query_identity, semantic
  )
}

.compile_effect_task <- function(left, over, right = NULL, query = NULL,
                                 normalizer = inner_product(), transform = NULL,
                                 reducer = .partition_reducer(), bridge = NULL) {
  .as_compiled_effect_task(.compile_effect_evidence_task(
    left = left,
    over = over,
    right = right,
    query = query,
    normalizer = normalizer,
    transform = transform,
    reducer = reducer,
    bridge = bridge
  ))
}

.validate_compiled_effect_task <- function(task) {
  expected <- c(
    "left_relation", "right_relation", "same_relation", "left_relation_id",
    "right_relation_id", "left_space", "right_space", "ordered_edges",
    "bridge", "normalizer", "transform", "lowering", "reducer",
    "query_identity", "source_uses",
    "distinct_handle_keys", "task_id"
  )
  if (!inherits(task, "effect_compiled_task") || !is.list(task) ||
      !identical(names(task), expected) || !is.logical(task$same_relation) ||
      length(task$same_relation) != 1L || is.na(task$same_relation)) {
    stop("Compiled effect-task fields are missing or noncanonical.",
      call. = FALSE)
  }
  .validate_relation(task$left_relation)
  .validate_relation(task$right_relation)
  left_id <- .relation_family_identity(task$left_relation)
  right_id <- .relation_family_identity(task$right_relation)
  if (!identical(task$left_relation_id, left_id) ||
      !identical(task$right_relation_id, right_id) ||
      (task$same_relation && !identical(left_id, right_id)) ||
      !identical(task$left_space, task$left_relation$effect_space) ||
      !identical(task$right_space, task$right_relation$effect_space)) {
    stop("Compiled effect-task relation or axis identity is inconsistent.",
      call. = FALSE)
  }
  .validate_ordered_partition_edges(
    task$ordered_edges,
    task$left_relation$partitions,
    task$right_relation$partitions,
    task$same_relation
  )
  .validate_measurement_bridge(task$bridge, task$left_relation,
    task$right_relation)
  .validate_edge_normalizer(task$normalizer)
  .validate_edge_transform(task$transform)
  .validate_partition_reducer(task$reducer)
  operation <- .edge_operation_plan(
    task$normalizer, task$transform, task$reducer
  )
  if (!identical(task$lowering, operation$lowering)) {
    stop("Compiled effect-task lowering is inconsistent with its operations.",
      call. = FALSE)
  }
  sources <- .effect_task_source_uses(
    task$left_relation, task$right_relation, left_id, right_id
  )
  if (!identical(task$source_uses, sources$uses) ||
      !identical(task$distinct_handle_keys, sources$distinct_handle_keys)) {
    stop("Compiled effect-task source-use identity is inconsistent.",
      call. = FALSE)
  }
  semantic <- .effect_task_semantic(
    left_id, right_id, task$left_space, task$right_space,
    task$ordered_edges, task$bridge, task$normalizer, task$transform,
    task$reducer, task$query_identity
  )
  if (!identical(task$task_id, .effect_task_id(semantic))) {
    stop("Compiled effect-task identity is inconsistent with its semantics.",
      call. = FALSE)
  }
  invisible(task)
}

.compiler_query <- function(query, effect_space, right_space = NULL,
                            codec = "symmetric_packed") {
  effect_space <- .validate_effect_space(effect_space)
  q <- length(effect_space$coordinates)
  if (identical(codec, "rectangular")) {
    right_space <- .validate_effect_space(right_space)
    q_right <- length(right_space$coordinates)
    rectangular_width <- q * q_right
    if (inherits(query, "effect_pair_query")) {
      .validate_query_for_compile(query)
      if (!.same_effect_space(query$left_space, effect_space) ||
          !.same_effect_space(query$right_space, right_space)) {
        stop("Pair-query axes must match the ordered relation effect spaces.",
          call. = FALSE)
      }
      value <- matrix(as.vector(as.matrix(query$operator)), ncol = 1L)
      colnames(value) <- "view1"
      return(value)
    }
    stop(paste0(
      "Rectangular plans execute axis-bound `pair_query()` operators; ",
      "materialize with `materialize_geometry()` and use `query_geometry()` ",
      "for raw physical view matrices."
    ), call. = FALSE)
  }
  packed_width <- q * (q + 1L) / 2L
  if (.is_pair_difference_query(query)) {
    if (!identical(query$effects, effect_space$coordinates)) {
      stop("The query and relation effect spaces are incompatible.",
        call. = FALSE)
    }
    return(query)
  }
  if (inherits(query, "effect_query")) {
    .validate_query_for_compile(query)
    if (!identical(query$kind, "bilinear") || !isTRUE(query$fixed)) {
      stop("Direct execution requires a fixed bilinear query.", call. = FALSE)
    }
    if (nrow(query$operator) != q) {
      stop("The query operator dimension must equal the effect dimension.",
        call. = FALSE)
    }
    if (!is.null(query$effect_space) &&
        !.same_effect_space(query$effect_space, effect_space)) {
      stop("The query and relation effect spaces are incompatible.",
        call. = FALSE)
    }
    value <- matrix(.svec_symmetric(query$operator), ncol = 1L)
    colnames(value) <- "view1"
    return(value)
  }
  if (!is.matrix(query) || !is.numeric(query) || nrow(query) != packed_width ||
      ncol(query) < 1L || any(!is.finite(query))) {
    stop("`query` must be a finite packed-coordinate-by-view matrix.",
      call. = FALSE)
  }
  if (is.null(colnames(query))) colnames(query) <- paste0("view", seq_len(ncol(query)))
  query
}

.component_requirements <- function(query, component) {
  if (is.null(query)) {
    return(list(total = TRUE, coherent = TRUE, marginals = TRUE,
      materialization = "full_geometry"))
  }
  component <- match.arg(component, c("total", "coherent", "configuration"))
  list(
    total = component %in% c("total", "configuration"),
    coherent = component %in% c("coherent", "configuration"),
    marginals = FALSE,
    materialization = paste0("direct_", component)
  )
}

.validate_compiler_inputs <- function(x, at, over, right = NULL) {
  .validate_relation(x)
  if (!is.null(right)) .validate_relation(right)
  .validate_frame_for_compile(at)
  .validate_pairing(over)
  if (!identical(at$representation, "additive_diagonal")) {
    stop("crossform 0.1 executes only additive diagonal frames.", call. = FALSE)
  }
  if (!.same_domain_reference(at$domain, x$domain) ||
      (!is.null(right) && !.same_domain_reference(at$domain, right$domain))) {
    stop("Relation and frame exact neural-domain identities must agree.",
      call. = FALSE)
  }
  if (ncol(at$weights) != x$n_features ||
      (!is.null(right) && ncol(at$weights) != right$n_features)) {
    stop("The frame feature dimension must equal the relation feature dimension.",
      call. = FALSE)
  }
  right_partitions <- if (is.null(right)) x$partitions else right$partitions
  if (any(!over$left %in% x$partitions) ||
      any(!over$right %in% right_partitions)) {
    stop("Every pairing endpoint must identify a relation partition.",
      call. = FALSE)
  }
  invisible(TRUE)
}

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
    stop("Support-metric planning requires a fixed metric schedule.",
      call. = FALSE)
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

.compiler_plan_id <- function(x, at, over, materialization, query, component,
                              normalizer = inner_product(), transform = NULL,
                              reducer = .partition_reducer()) {
  task <- .compile_effect_task(
    x, over, query = query, normalizer = normalizer,
    transform = transform, reducer = reducer
  )
  semantic <- list(
    effect_task_id = task$task_id,
    left_space = task$left_space,
    right_space = task$right_space,
    ordered_edges = unclass(as.data.frame(task$ordered_edges)),
    edge_expansion = attr(task$ordered_edges, "expansion"),
    bridge = task$bridge$signature,
    normalizer = unclass(task$normalizer),
    transform = unclass(task$transform),
    lowering = task$lowering,
    reducer = unclass(task$reducer),
    frame = list(weights = at$weights, normalization = at$normalization,
      domain = at$domain),
    materialization = materialization,
    component = component
  )
  paste0("sha256:", digest::digest(semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L))
}

.compiler_index <- function(at) {
  index <- at$index
  if (is.null(index)) return(seq_len(nrow(at$weights)))
  if (is.data.frame(index) && "measurement" %in% names(index)) {
    return(index$measurement)
  }
  index
}

.compiler_blas <- function() {
  vendor <- tryCatch(unname(extSoftVersion()[["BLAS"]]), error = function(e) NULL)
  if (is.null(vendor) || !is.character(vendor) || length(vendor) != 1L ||
      is.na(vendor) || !nzchar(vendor)) vendor <- "unknown"
  list(vendor = vendor, requested_threads = 1L,
    observed_threads = NA_integer_)
}

.planned_compiler_receipt <- function(plan, validate = TRUE) {
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
    blas = .compiler_blas(),
    domain_signature = task$left_relation$domain$signature
  )
}

.compiler_reporter <- function(reporter, final_receipt) {
  force(reporter)
  function(event) {
    if (event$type %in% c("complete", "failed", "interrupted")) {
      final_receipt$value <- event$receipt
    }
    if (!is.null(reporter)) reporter(event)
    invisible(NULL)
  }
}

.compiler_storage <- function(storage, storage_path, dim) {
  storage <- match.arg(storage, c("memory", "block"))
  if (storage == "memory") {
    if (!is.null(storage_path)) {
      stop("`storage_path` is only valid for block-backed geometry.",
        call. = FALSE)
    }
    return(list(kind = storage, total = NULL, coherent = NULL,
      created = character(), created_directory = FALSE))
  }
  if (!is.character(storage_path) || length(storage_path) != 1L ||
      is.na(storage_path) || !nzchar(storage_path)) {
    stop("Block-backed geometry requires one nonempty `storage_path`.",
      call. = FALSE)
  }
  created_directory <- FALSE
  if (!dir.exists(storage_path)) {
    parent <- dirname(storage_path)
    if (!dir.exists(parent)) {
      stop("The parent of `storage_path` must already exist.", call. = FALSE)
    }
    if (!dir.create(storage_path, recursive = FALSE, showWarnings = FALSE)) {
      stop("Could not create `storage_path`.", call. = FALSE)
    }
    created_directory <- TRUE
  }
  total_path <- file.path(storage_path, "total.egm")
  coherent_path <- file.path(storage_path, "coherent.egm")
  if (file.exists(total_path) || file.exists(coherent_path)) {
    stop("Refusing to overwrite an existing geometry component.", call. = FALSE)
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

.compiler_cleanup <- function(storage, completed) {
  force(storage)
  force(completed)
  function() {
    if (!isTRUE(completed$value) && length(storage$created)) {
      unlink(storage$created)
      remaining <- storage$created[file.exists(storage$created)]
      if (length(remaining)) {
        stop("Failed to remove incomplete geometry components.", call. = FALSE)
      }
      if (storage$created_directory) {
        directory <- dirname(storage$created[[1L]])
        unlink(directory, recursive = TRUE)
        if (dir.exists(directory)) {
          stop("Failed to remove incomplete geometry directory.", call. = FALSE)
        }
      }
    }
    invisible(NULL)
  }
}

.ordered_pairing_marginals <- function(local, ordered_edges, mass) {
  if (!is.array(local) || length(dim(local)) != 3L ||
      !is.numeric(local) || any(!is.finite(local))) {
    stop("Compiled marginals require finite local relation first moments.",
      call. = FALSE)
  }
  partitions <- dimnames(local)[[3L]]
  effects <- dimnames(local)[[2L]]
  .validate_ordered_partition_edges(
    ordered_edges, partitions, partitions, TRUE
  )
  measurements <- dim(local)[[1L]]
  if (!is.numeric(mass) || !(length(mass) %in% c(1L, measurements)) ||
      any(!is.finite(mass)) || any(mass <= 0)) {
    stop("Compiled marginal mass must be positive and finite.", call. = FALSE)
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

.execute_geometry_plan <- function(plan, reporter = NULL) {
  .validate_geometry_execution_plan(plan)
  task <- .as_compiled_effect_task(plan$task, validate = FALSE)
  at <- plan$frame
  query <- plan$query
  component <- plan$component
  requirements <- plan$requirements
  storage <- plan$storage
  support_streamed <- grepl(
    "^support_streamed_metric_", plan$lowering
  )
  explicit_additive_metric <- grepl(
    "^additive_metric_", plan$lowering
  )
  planned_receipt <- .planned_compiler_receipt(plan, validate = FALSE)
  geometry_storage <- .compiler_storage(
    storage, plan$storage_path, c(nrow(at$weights), plan$output_width)
  )
  completed <- new.env(parent = emptyenv())
  completed$value <- FALSE
  final_receipt <- new.env(parent = emptyenv())
  final_receipt$value <- planned_receipt
  observed <- .empty_execution_observations()
  observed$task_counts[["planned"]] <- plan$task_count
  observed$tiles <- list(
    feature_block = plan$feature_block,
    row_tile = plan$row_tile,
    coordinate_tile = plan$coordinate_tile
  )
  observed_state <- new.env(parent = emptyenv())
  observed_state$value <- observed
  task_observer <- function(event, features) {
    observed_state$value$task_counts[[event]] <-
      observed_state$value$task_counts[[event]] + 1
    if (identical(event, "completed")) {
      observed_state$value$features_completed <-
        observed_state$value$features_completed + length(features)
    }
    invisible(NULL)
  }

  value <- .execute_guarded(
    compute = function() {
      admission_started <- proc.time()[["elapsed"]]
      source_session <- .open_effect_task_source_session(task, validate = FALSE)
      session_closed <- FALSE
      close_session <- function() {
        if (!session_closed) {
          session_closed <<- TRUE
          source_session$close()
        }
      }
      on.exit(close_session(), add = TRUE)
      # Two-sided sessions report per-side summaries; flatten them with
      # side-qualified names so receipt observations keep one canonical
      # shape on both self-form and rectangular executions.
      session_summary <- function() {
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
      observed_state$value$stage_seconds[["source_admission"]] <-
        max(0, proc.time()[["elapsed"]] - admission_started)
      observed_state$value$source_access <- session_summary()$access_mode
      total_accumulator <- NULL
      if (storage == "block") {
        total_accumulator <- function(rows, coordinates, increment) {
          existing <- .read_geometry_tile(geometry_storage$total, rows, coordinates)
          .write_geometry_tile(geometry_storage$total, rows, coordinates,
            existing + increment)
        }
      }
      kernel_started <- proc.time()[["elapsed"]]
      same_relation <- isTRUE(task$same_relation)
      codec <- if (same_relation) "symmetric_packed" else "rectangular"
      read_relation <- function(partition, features) {
          value <- .relation_block_with_reader(
            task$left_relation, partition, features,
            function(partition, features) {
              source_session$read("left", partition, features)
            },
            validate = FALSE
          )
          observed_state$value$bytes_read <-
            sum(session_summary()$bytes_read)
          value
      }
      read_right_relation <- if (same_relation) read_relation else
        function(partition, features) {
          value <- .relation_block_with_reader(
            task$right_relation, partition, features,
            function(partition, features) {
              source_session$read("right", partition, features)
            },
            validate = FALSE
          )
          observed_state$value$bytes_read <-
            sum(session_summary()$bytes_read)
          value
        }
      streamed <- if (support_streamed) {
        .support_streamed_metric_contraction(
          frame = at,
          metric_schedule = plan$metric_schedule,
          read_relation = read_relation,
          partitions = task$left_relation$partitions,
          effects = task$left_space$coordinates,
          ordered_edges = task$ordered_edges,
          query = query,
          form_total = requirements$total,
          form_coherent = requirements$coherent,
          retain_first_moments = requirements$marginals,
          task_observer = task_observer
        )
      } else {
        kernel_frame <- if (explicit_additive_metric) {
          .metric_additive_frame(at, plan$metric_schedule)
        } else {
          at
        }
        .streamed_effect_form_contraction(
          frame = kernel_frame,
          read_left = read_relation,
          read_right = read_right_relation,
          left_partitions = task$left_relation$partitions,
          right_partitions = task$right_relation$partitions,
          left_effects = task$left_space$coordinates,
          right_effects = task$right_space$coordinates,
          ordered_edges = task$ordered_edges,
          codec = codec,
          same_relation = same_relation,
          feature_block = plan$feature_block,
          row_tile = plan$row_tile,
          coordinate_tile = plan$coordinate_tile,
          accumulate_tile = total_accumulator,
          retain_first_moments = requirements$coherent,
          query = query,
          form_total = requirements$total,
          task_observer = task_observer
        )
      }
      observed_state$value$stage_seconds[["feature_tasks"]] <-
        max(0, proc.time()[["elapsed"]] - kernel_started)
      coherent_writer <- if (storage == "block") {
        function(rows, coordinates, value) {
          .write_geometry_tile(geometry_storage$coherent, rows, coordinates, value)
        }
      } else {
        NULL
      }
      coherent_started <- proc.time()[["elapsed"]]
      coherent <- if (requirements$coherent && support_streamed) {
        list(value = streamed$coherent, diagnostics = list(
          lowering = "metric_aware_rank_one_component",
          factorization_count = streamed$diagnostics$metric_factorizations,
          durable_output_bytes = 8 * length(streamed$coherent)
        ))
      } else if (requirements$coherent) {
        .effect_form_coherent_from_first_moments(
          streamed$first_moments$left,
          streamed$first_moments$right,
          task$ordered_edges,
          streamed$mass,
          codec = codec,
          same_relation = same_relation,
          row_tile = plan$row_tile,
          write_tile = coherent_writer,
          query = query
        )
      } else {
        NULL
      }
      observed_state$value$stage_seconds[["coherent"]] <-
        max(0, proc.time()[["elapsed"]] - coherent_started)
      total_value <- if (!requirements$total) NULL else if (storage == "memory") {
        streamed$value
      } else {
        geometry_storage$total
      }
      coherent_value <- if (!requirements$coherent) NULL else if (
        storage == "memory") {
        coherent$value
      } else {
        geometry_storage$coherent
      }
      marginal_started <- proc.time()[["elapsed"]]
      marginals <- if (requirements$marginals) {
        .ordered_pairing_marginals(
          streamed$first_moments$left,
          task$ordered_edges,
          streamed$mass
        )
      } else {
        NULL
      }
      observed_state$value$stage_seconds[["marginals"]] <-
        max(0, proc.time()[["elapsed"]] - marginal_started)
      metadata <- list(
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
        storage = storage,
        tiles = list(
          feature_block = plan$feature_block,
          row_tile = plan$row_tile,
          coordinate_tile = plan$coordinate_tile
        ),
        requirements = requirements,
        diagnostics = list(
          total = if (requirements$total) {
            diagnostics <- streamed$diagnostics
            first_moment_bytes <- diagnostics$durable_left_first_moment_bytes
            if (is.null(first_moment_bytes)) {
              first_moment_bytes <- diagnostics$durable_first_moment_bytes
            }
            if (!is.null(first_moment_bytes)) {
              diagnostics$durable_local_relation_bytes <- first_moment_bytes
            }
            diagnostics
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
      close_session()
      source_summary <- session_summary()
      metadata$source_session <- source_summary
      observed_state$value$bytes_read <- sum(source_summary$bytes_read)
      result <- if (is.null(query) && !same_relation) {
        rectangular_form <- effect_form(
          total = total_value,
          coherent = coherent_value,
          marginals = NULL,
          left_space = task$left_space,
          right_space = task$right_space,
          receipt = planned_receipt,
          index = .compiler_index(at),
          metadata = metadata,
          codec = "rectangular",
          symmetric = FALSE,
          guaranteed_psd = FALSE
        )
        rectangular_form
      } else if (is.null(query)) {
        effect_geometry(total_value, coherent_value, marginals,
          effects = task$left_space, receipt = planned_receipt,
          index = .compiler_index(at), metadata = metadata)
      } else if (identical(component, "contrast")) {
        .new_effect_contrast_view(
          total_value, coherent_value, marginals, plan$signed_query,
          .compiler_index(at), planned_receipt
        )
      } else {
        total_matrix <- if (requirements$total) total_value else NULL
        coherent_matrix <- if (requirements$coherent) coherent_value else NULL
        values <- switch(component,
          total = total_matrix,
          coherent = coherent_matrix,
          configuration = total_matrix - coherent_matrix
        )
        colnames(values) <- .query_output_labels(query)
        effect_view(values, query, component, planned_receipt,
          effects = if (same_relation) task$left_space else NULL,
          left_space = if (same_relation) NULL else task$left_space,
          right_space = if (same_relation) NULL else task$right_space,
          index = .compiler_index(at), metadata = metadata)
      }
      completed$value <- TRUE
      result
    },
    receipt = planned_receipt,
    reporter = .compiler_reporter(reporter, final_receipt),
    cleanup = .compiler_cleanup(geometry_storage, completed),
    observations = function() observed_state$value,
    receipt_sink = function(receipt) final_receipt$value <- receipt
  )
  value$receipt <- final_receipt$value
  value
}

.run_geometry_compiler <- function(x, at = NULL, over = NULL,
                                   compute = NULL, storage = "memory",
                                   storage_path = NULL, query = NULL,
                                   component = NULL, reporter = NULL,
                                   signed_query = NULL) {
  if (inherits(x, "effect_geometry_plan")) {
    if (!is.null(at) || !is.null(over) || !is.null(compute)) {
      stop("A compiled geometry plan cannot be combined with raw plan inputs.",
        call. = FALSE)
    }
    plan <- x
  } else {
    if (is.null(at) || is.null(over)) {
      stop("A relation requires both `at` and `over` to compile geometry.",
        call. = FALSE)
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

#' Materialize a complete cross-generalized geometry
#'
#' `materialize_geometry()` explicitly materializes complete total and
#' coherent packed geometry. Prefer [plan_geometry()] followed by
#' [evaluate_geometry()] when a fixed query is sufficient. The relation
#' compatibility form compiles that same plan first; it does not use a second
#' execution path.
#'
#' @param x An `effect_geometry_plan` or, for compatibility, an
#'   `effect_relation`.
#' @param at A compiled additive `effect_frame`, required only with a relation.
#' @param over An `effect_pairing`, required only with a relation.
#' @param storage Either `"memory"` or `"block"`.
#' @param storage_path Durable directory for block-backed geometry.
#' @param compute A sequential `compute_policy()` when `x` is a relation. A
#'   compiled plan already owns this choice.
#' @param reporter Optional nonsemantic coordinator-side event reporter.
#' @return A complete `effect_geometry`.
#' @export
materialize_geometry <- function(x, at = NULL, over = NULL,
                                 storage = c("memory", "block"),
                                 storage_path = NULL, compute = NULL,
                                 reporter = NULL) {
  storage <- match.arg(storage)
  .run_geometry_compiler(x, at, over, compute, storage, storage_path,
    query = NULL, component = NULL, reporter = reporter)
}

#' Evaluate a fixed query without materializing complete geometry
#'
#' @param x,at,over,compute,reporter As in [materialize_geometry()].
#' @param query A fixed `bilinear_query()` or packed-coordinate query matrix.
#' @param component One of `total`, `coherent`, or `configuration`.
#' @return A query-only `effect_view`.
#' @export
evaluate_geometry <- function(x, at = NULL, over = NULL, query = NULL,
                              component = c("total", "coherent", "configuration"),
                              compute = NULL, reporter = NULL) {
  if (inherits(x, "effect_geometry_plan") && is.null(query) &&
      !is.null(at) && is.null(over)) {
    query <- at
    at <- NULL
  }
  if (is.null(query)) {
    stop("Query-first execution requires an explicit fixed `query`.",
      call. = FALSE)
  }
  component <- match.arg(component)
  .run_geometry_compiler(x, at, over, compute, "memory", NULL,
    query = query, component = component, reporter = reporter)
}
