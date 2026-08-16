# Boundary-typed evidence-task intermediate representation -----------------

.evidence_boundary_signature <- function(semantic) {
  paste0("sha256:", digest::digest(
    semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
  ))
}

.open_experimental_boundary <- function(left_space, right_space) {
  left_space <- .validate_effect_space(left_space)
  right_space <- .validate_effect_space(right_space)
  semantic <- list(
    schema_version = 1L,
    boundary = "experimental",
    state = "open",
    left_space = left_space,
    right_space = right_space,
    query = NULL,
    role = NULL,
    sampling_axis = NULL
  )
  structure(c(semantic[-1L], list(
    signature = .evidence_boundary_signature(semantic)
  )), class = "effect_evidence_boundary")
}

.closed_experimental_boundary <- function(
    query, role = c("effect", "variation"), sampling_axis = NULL) {
  .validate_query_for_compile(query)
  if (!inherits(query, "effect_pair_query")) {
    stop("A closed experimental boundary requires an axis-bound pair query.",
      call. = FALSE)
  }
  role <- match.arg(role)
  if (role == "variation" &&
      (!is.character(sampling_axis) || length(sampling_axis) != 1L ||
       is.na(sampling_axis) || !nzchar(sampling_axis))) {
    stop("A variation boundary requires one declared sampling axis.",
      call. = FALSE)
  }
  if (role == "effect" && !is.null(sampling_axis)) {
    stop("An effect boundary cannot claim a repeated-sampling axis.",
      call. = FALSE)
  }
  semantic <- list(
    schema_version = 1L,
    boundary = "experimental",
    state = "closed",
    left_space = query$left_space,
    right_space = query$right_space,
    query = query,
    role = role,
    sampling_axis = sampling_axis
  )
  structure(c(semantic[-1L], list(
    signature = .evidence_boundary_signature(semantic)
  )), class = "effect_evidence_boundary")
}

.closed_neural_boundary <- function(bridge, left_relation, right_relation) {
  bridge <- .validate_measurement_bridge(
    bridge, left_relation, right_relation
  )
  semantic <- list(
    schema_version = 1L,
    boundary = "neural",
    state = "closed",
    left_space = left_relation$domain,
    right_space = right_relation$domain,
    closure_kind = "bridge",
    bridge = bridge,
    query = NULL,
    left_frame = NULL,
    right_frame = NULL,
    edges = NULL
  )
  structure(c(semantic[-1L], list(
    signature = .evidence_boundary_signature(list(
      schema_version = semantic$schema_version,
      boundary = semantic$boundary,
      state = semantic$state,
      left_space = semantic$left_space,
      right_space = semantic$right_space,
      closure_kind = semantic$closure_kind,
      bridge = bridge$signature
    ))
  )), class = "effect_evidence_boundary")
}

.neural_pair_query <- function(operator, left_domain, right_domain,
                               provenance = list()) {
  if (!is.matrix(operator) || !is.numeric(operator) ||
      any(dim(operator) < 1L) || any(!is.finite(operator))) {
    stop("A neural pair query must be a finite nonempty numeric matrix.",
      call. = FALSE)
  }
  operator <- unname(operator)
  storage.mode(operator) <- "double"
  left_domain <- .domain_reference(left_domain)
  right_domain <- .domain_reference(right_domain)
  if (nrow(operator) != left_domain$n_features ||
      ncol(operator) != right_domain$n_features) {
    stop("Neural pair-query dimensions do not match their identified domains.",
      call. = FALSE)
  }
  .validate_effect_provenance(provenance, "neural pair-query provenance")
  semantic <- list(
    schema_version = 1L,
    operator = operator,
    left_domain = left_domain,
    right_domain = right_domain,
    provenance = provenance
  )
  structure(c(semantic[-1L], list(
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    ))
  )), class = "effect_neural_pair_query")
}

.validate_neural_pair_query <- function(x) {
  expected <- c("operator", "left_domain", "right_domain", "provenance",
    "signature")
  if (!inherits(x, "effect_neural_pair_query") || !is.list(x) ||
      !identical(names(x), expected)) {
    stop("Neural pair-query fields are missing or noncanonical.",
      call. = FALSE)
  }
  rebuilt <- .neural_pair_query(
    x$operator, x$left_domain, x$right_domain, x$provenance
  )
  if (!identical(x, rebuilt)) {
    stop("Neural pair-query identity is inconsistent.", call. = FALSE)
  }
  rebuilt
}

.reverse_neural_pair_query <- function(x) {
  x <- .validate_neural_pair_query(x)
  .neural_pair_query(
    t(x$operator), x$right_domain, x$left_domain, x$provenance
  )
}

.closed_neural_query_boundary <- function(query) {
  query <- .validate_neural_pair_query(query)
  semantic <- list(
    schema_version = 1L,
    boundary = "neural",
    state = "closed",
    left_space = query$left_domain,
    right_space = query$right_domain,
    closure_kind = "query",
    bridge = NULL,
    query = query,
    left_frame = NULL,
    right_frame = NULL,
    edges = NULL
  )
  structure(c(semantic[-1L], list(
    signature = .evidence_boundary_signature(list(
      schema_version = semantic$schema_version,
      boundary = semantic$boundary,
      state = semantic$state,
      left_space = semantic$left_space,
      right_space = semantic$right_space,
      closure_kind = semantic$closure_kind,
      query = query$signature
    ))
  )), class = "effect_evidence_boundary")
}

.open_neural_boundary <- function(left_frame, right_frame = left_frame,
                                  edges) {
  left_frame <- .validate_measurement_frame(left_frame)
  right_frame <- .validate_measurement_frame(right_frame)
  edges <- .validate_measurement_edges(edges, left_frame, right_frame)
  semantic <- list(
    schema_version = 1L,
    boundary = "neural",
    state = "open",
    left_space = left_frame$source_domain,
    right_space = right_frame$source_domain,
    closure_kind = NULL,
    bridge = NULL,
    query = NULL,
    left_frame = left_frame,
    right_frame = right_frame,
    edges = edges
  )
  structure(c(semantic[-1L], list(
    signature = .evidence_boundary_signature(list(
      schema_version = semantic$schema_version,
      boundary = semantic$boundary,
      state = semantic$state,
      left_space = semantic$left_space,
      right_space = semantic$right_space,
      closure_kind = semantic$closure_kind,
      left_frame = left_frame$signature,
      right_frame = right_frame$signature,
      edges = edges$signature
    ))
  )), class = "effect_evidence_boundary")
}

.validate_evidence_boundary <- function(x) {
  if (!inherits(x, "effect_evidence_boundary") || !is.list(x) ||
      !identical(x$boundary %in% c("experimental", "neural"), TRUE) ||
      !identical(x$state %in% c("open", "closed"), TRUE)) {
    stop("Evidence-boundary fields are missing or noncanonical.",
      call. = FALSE)
  }
  if (x$boundary == "experimental") {
    expected <- c("boundary", "state", "left_space", "right_space", "query",
      "role", "sampling_axis", "signature")
    if (!identical(names(x), expected)) {
      stop("Experimental-boundary fields are missing or noncanonical.",
        call. = FALSE)
    }
    rebuilt <- if (x$state == "open") {
      if (!is.null(x$query) || !is.null(x$role) || !is.null(x$sampling_axis)) {
        stop("An open experimental boundary cannot carry a closure query.",
          call. = FALSE)
      }
      .open_experimental_boundary(x$left_space, x$right_space)
    } else {
      .closed_experimental_boundary(x$query, x$role, x$sampling_axis)
    }
  } else {
    expected <- c("boundary", "state", "left_space", "right_space",
      "closure_kind", "bridge", "query", "left_frame", "right_frame",
      "edges", "signature")
    if (!identical(names(x), expected)) {
      stop("Neural-boundary fields are missing or noncanonical.",
        call. = FALSE)
    }
    if (x$state == "closed") {
      if (!is.character(x$closure_kind) || length(x$closure_kind) != 1L ||
          is.na(x$closure_kind) ||
          !x$closure_kind %in% c("bridge", "query") ||
          !is.null(x$left_frame) ||
          !is.null(x$right_frame) || !is.null(x$edges)) {
        stop("A closed neural boundary requires one bridge or query closure.",
          call. = FALSE)
      }
      if (x$closure_kind == "bridge") {
        if (is.null(x$bridge) || !is.null(x$query) ||
            !inherits(x$bridge, "effect_measurement_bridge")) {
          stop("A bridge-closed neural boundary requires only its bridge.",
            call. = FALSE)
        }
        .validate_domain_reference(x$left_space)
        .validate_domain_reference(x$right_space)
        expected_signature <- .evidence_boundary_signature(list(
          schema_version = 1L,
          boundary = "neural",
          state = "closed",
          left_space = x$left_space,
          right_space = x$right_space,
          closure_kind = "bridge",
          bridge = x$bridge$signature
        ))
        if (!identical(x$signature, expected_signature)) {
          stop("Closed neural-boundary identity is inconsistent.",
            call. = FALSE)
        }
        rebuilt <- x
      } else {
        if (!is.null(x$bridge) || is.null(x$query)) {
          stop("A query-closed neural boundary requires only its query.",
            call. = FALSE)
        }
        rebuilt <- .closed_neural_query_boundary(x$query)
      }
    } else {
      if (!is.null(x$closure_kind) || !is.null(x$bridge) ||
          !is.null(x$query) || is.null(x$left_frame) ||
          is.null(x$right_frame) || is.null(x$edges)) {
        stop("An open neural boundary requires two frames and explicit edges.",
          call. = FALSE)
      }
      rebuilt <- .open_neural_boundary(
        x$left_frame, x$right_frame, x$edges
      )
    }
  }
  if (!identical(x, rebuilt)) {
    stop("Evidence-boundary identity is inconsistent.", call. = FALSE)
  }
  rebuilt
}

.evidence_materialization <- function(
    kind = c("effect_form", "measurement_form", "scalar_field"),
    completeness = c("complete_form", "query_only"), projection = NULL) {
  kind <- match.arg(kind)
  completeness <- match.arg(completeness)
  if (kind == "scalar_field" && completeness != "query_only") {
    stop("A scalar field is necessarily a query-only materialization.",
      call. = FALSE)
  }
  if (!is.null(projection) && !is.list(projection)) {
    stop("Materialization projection identity must be a list or NULL.",
      call. = FALSE)
  }
  semantic <- list(
    schema_version = 1L,
    kind = kind,
    completeness = completeness,
    projection = projection
  )
  structure(c(semantic[-1L], list(
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    ))
  )), class = "effect_evidence_materialization")
}

.validate_evidence_materialization <- function(x) {
  expected <- c("kind", "completeness", "projection", "signature")
  if (!inherits(x, "effect_evidence_materialization") || !is.list(x) ||
      !identical(names(x), expected)) {
    stop("Evidence-materialization fields are missing or noncanonical.",
      call. = FALSE)
  }
  rebuilt <- .evidence_materialization(
    x$kind, x$completeness, x$projection
  )
  if (!identical(x, rebuilt)) {
    stop("Evidence-materialization identity is inconsistent.",
      call. = FALSE)
  }
  rebuilt
}

.evidence_stage_plan <- function(normalizer = inner_product(), transform = NULL,
                                 reducer = reduce_partitions(),
                                 projection = NULL) {
  normalizer <- .validate_edge_normalizer(normalizer)
  if (is.null(transform)) transform <- .identity_edge_transform()
  transform <- .validate_edge_transform(transform)
  reducer <- .validate_partition_reducer(reducer)
  operation <- .edge_operation_plan(normalizer, transform, reducer)
  normalization_placement <- if (reducer$order == "edge_first") {
    "within_partition_pair"
  } else {
    "after_partition_aggregation"
  }
  order <- if (reducer$order == "edge_first") {
    c("ordered_partition_product", "neural_feature_normalization",
      "edge_transform", "partition_reduction", "experimental_projection")
  } else {
    c("ordered_partition_product", "partition_sufficient_reduction",
      "neural_feature_normalization", "edge_transform",
      "experimental_projection")
  }
  semantic <- list(
    schema_version = 1L,
    order = order,
    partition_product = list(
      operation = "ordered_partition_product",
      axis = "partition_pairs"
    ),
    normalization = list(
      operation = unclass(normalizer),
      axis = "neural_features",
      placement = normalization_placement
    ),
    transform = list(
      operation = unclass(transform),
      axis = "form_entries",
      placement = "before_partition_reduction"
    ),
    reduction = list(
      operation = unclass(reducer),
      axis = "partition_pairs"
    ),
    projection = list(
      operation = projection,
      axis = "experimental_coordinates"
    ),
    lowering = operation$lowering
  )
  structure(c(semantic[-1L], list(
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    ))
  )), class = "effect_evidence_stages")
}

.validate_evidence_stage_plan <- function(x) {
  expected <- c("order", "partition_product", "normalization", "transform",
    "reduction", "projection", "lowering", "signature")
  if (!inherits(x, "effect_evidence_stages") || !is.list(x) ||
      !identical(names(x), expected)) {
    stop("Evidence-stage fields are missing or noncanonical.",
      call. = FALSE)
  }
  normalizer <- structure(x$normalization$operation,
    class = "effect_edge_normalizer")
  transform <- structure(x$transform$operation,
    class = "effect_edge_transform")
  reducer <- structure(x$reduction$operation,
    class = "effect_partition_reducer")
  rebuilt <- .evidence_stage_plan(
    normalizer, transform, reducer, x$projection$operation
  )
  if (!identical(x, rebuilt)) {
    stop("Evidence-stage identity or operation order is inconsistent.",
      call. = FALSE)
  }
  rebuilt
}

.evidence_task_general_semantic <- function(
    left_id, right_id, spaces, ordered_edges, experimental_boundary,
    neural_boundary, stages, materialization) {
  list(
    schema_version = 1L,
    transport = "evidence_pairing",
    left_relation = left_id,
    right_relation = right_id,
    spaces = lapply(spaces, `[[`, "signature"),
    ordered_partition_products = unclass(as.data.frame(ordered_edges)),
    partition_expansion = attr(ordered_edges, "expansion"),
    experimental_boundary = experimental_boundary$signature,
    neural_boundary = neural_boundary$signature,
    stages = stages$signature,
    materialization = materialization$signature
  )
}

.validate_evidence_boundary_combination <- function(experimental_boundary,
                                                    neural_boundary,
                                                    materialization) {
  state <- paste(experimental_boundary$state, neural_boundary$state, sep = "/")
  valid <- switch(materialization$kind,
    effect_form = identical(state, "open/closed"),
    measurement_form = identical(state, "closed/open"),
    scalar_field = identical(state, "closed/closed"),
    FALSE
  )
  if (!isTRUE(valid)) {
    stop("Requested materialization is incompatible with its open boundaries.",
      call. = FALSE)
  }
  invisible(TRUE)
}

.new_evidence_task <- function(
    left_relation, right_relation, same_relation, ordered_edges,
    experimental_boundary, neural_boundary, stages, materialization,
    identity_schema = c("evidence-pairing-v1", "effect-form-v1"),
    compatibility_semantic = NULL) {
  .validate_relation(left_relation)
  .validate_relation(right_relation)
  left_relation$capabilities <- .compiler_capabilities(left_relation)
  right_relation$capabilities <- if (isTRUE(same_relation)) {
    left_relation$capabilities
  } else {
    .compiler_capabilities(right_relation)
  }
  if (!is.logical(same_relation) || length(same_relation) != 1L ||
      is.na(same_relation)) {
    stop("Evidence-task `same_relation` must be TRUE or FALSE.",
      call. = FALSE)
  }
  left_id <- .relation_family_identity(left_relation)
  right_id <- if (same_relation) left_id else
    .relation_family_identity(right_relation)
  .validate_ordered_partition_edges(
    ordered_edges, left_relation$partitions, right_relation$partitions,
    same_relation
  )
  experimental_boundary <- .validate_evidence_boundary(
    experimental_boundary
  )
  neural_boundary <- .validate_evidence_boundary(neural_boundary)
  if (experimental_boundary$boundary != "experimental" ||
      neural_boundary$boundary != "neural") {
    stop("Evidence-task boundaries occupy the wrong typed axes.",
      call. = FALSE)
  }
  if (!.same_effect_space(experimental_boundary$left_space,
        left_relation$effect_space) ||
      !.same_effect_space(experimental_boundary$right_space,
        right_relation$effect_space)) {
    stop("Experimental boundary axes differ from relation effect spaces.",
      call. = FALSE)
  }
  if (!.same_domain_reference(neural_boundary$left_space,
        left_relation$domain) ||
      !.same_domain_reference(neural_boundary$right_space,
        right_relation$domain)) {
    stop("Neural boundary axes differ from relation neural spaces.",
      call. = FALSE)
  }
  if (neural_boundary$state == "closed") {
    if (neural_boundary$closure_kind == "bridge") {
      .validate_measurement_bridge(
        neural_boundary$bridge, left_relation, right_relation
      )
    } else {
      .validate_neural_pair_query(neural_boundary$query)
    }
  } else {
    .validate_measurement_frame(neural_boundary$left_frame)
    .validate_measurement_frame(neural_boundary$right_frame)
    .validate_measurement_edges(
      neural_boundary$edges,
      neural_boundary$left_frame,
      neural_boundary$right_frame
    )
  }
  stages <- .validate_evidence_stage_plan(stages)
  materialization <- .validate_evidence_materialization(materialization)
  .validate_evidence_boundary_combination(
    experimental_boundary, neural_boundary, materialization
  )
  identity_schema <- match.arg(identity_schema)
  spaces <- list(
    experimental_left = left_relation$effect_space,
    experimental_right = right_relation$effect_space,
    neural_left = left_relation$domain,
    neural_right = right_relation$domain
  )
  if (identity_schema == "effect-form-v1") {
    if (is.null(compatibility_semantic) ||
        materialization$kind != "effect_form" ||
        neural_boundary$closure_kind != "bridge") {
      stop("Effect-form compatibility identity requires its legacy semantic.",
        call. = FALSE)
    }
    expected_compatibility <- .effect_task_semantic(
      left_id, right_id,
      left_relation$effect_space, right_relation$effect_space,
      ordered_edges, neural_boundary$bridge,
      structure(stages$normalization$operation,
        class = "effect_edge_normalizer"),
      structure(stages$transform$operation,
        class = "effect_edge_transform"),
      structure(stages$reduction$operation,
        class = "effect_partition_reducer"),
      materialization$projection
    )
    if (!identical(compatibility_semantic, expected_compatibility)) {
      stop("Effect-form compatibility semantic differs from certified fields.",
        call. = FALSE)
    }
    semantic <- expected_compatibility
  } else {
    if (!is.null(compatibility_semantic)) {
      stop("Native evidence tasks cannot carry a legacy compatibility semantic.",
        call. = FALSE)
    }
    semantic <- .evidence_task_general_semantic(
      left_id, right_id, spaces, ordered_edges, experimental_boundary,
      neural_boundary, stages, materialization
    )
  }
  sources <- .effect_task_source_uses(
    left_relation, right_relation, left_id, right_id
  )
  task_id <- .effect_task_id(semantic)
  structure(list(
    schema_version = 1L,
    left_relation = left_relation,
    right_relation = right_relation,
    same_relation = same_relation,
    left_relation_id = left_id,
    right_relation_id = right_id,
    spaces = spaces,
    ordered_partition_products = ordered_edges,
    experimental_boundary = experimental_boundary,
    neural_boundary = neural_boundary,
    stages = stages,
    materialization = materialization,
    source_uses = sources$uses,
    distinct_handle_keys = sources$distinct_handle_keys,
    identity_schema = identity_schema,
    semantic = semantic,
    task_id = task_id
  ), class = "effect_evidence_task")
}

.validate_evidence_task <- function(task) {
  expected <- c(
    "schema_version", "left_relation", "right_relation", "same_relation",
    "left_relation_id", "right_relation_id", "spaces",
    "ordered_partition_products", "experimental_boundary", "neural_boundary",
    "stages", "materialization", "source_uses", "distinct_handle_keys",
    "identity_schema", "semantic", "task_id"
  )
  if (!inherits(task, "effect_evidence_task") || !is.list(task) ||
      !identical(names(task), expected) ||
      !identical(task$schema_version, 1L)) {
    stop("Evidence-task fields are missing or noncanonical.", call. = FALSE)
  }
  rebuilt <- .new_evidence_task(
    task$left_relation, task$right_relation, task$same_relation,
    task$ordered_partition_products, task$experimental_boundary,
    task$neural_boundary, task$stages, task$materialization,
    task$identity_schema,
    if (task$identity_schema == "effect-form-v1") task$semantic else NULL
  )
  if (!identical(task, rebuilt)) {
    stop("Evidence-task identity is inconsistent with its typed semantics.",
      call. = FALSE)
  }
  invisible(task)
}

.reverse_query_identity <- function(x) {
  if (is.null(x) || identical(x$kind, "complete_form")) return(x)
  if (identical(x$kind, "pair_query")) {
    query <- x$query
    return(list(
      kind = "pair_query",
      query = pair_query(
        t(as.matrix(query$operator)), query$right_space, query$left_space,
        query$metadata
      )
    ))
  }
  if (x$kind %in% c("bilinear_compatibility", "physical_self_query",
      "pair_difference_self_query")) {
    return(x)
  }
  stop("Unknown effect-form projection identity.", call. = FALSE)
}

.reverse_evidence_task <- function(task) {
  .validate_evidence_task(task)
  edges <- task$ordered_partition_products
  reversed_edges <- edges
  if (!identical(attr(edges, "expansion"), "self_adjoint_half_edges")) {
    reversed_edges$left <- edges$right
    reversed_edges$right <- edges$left
  }
  attr(reversed_edges, "source_estimate") <- attr(edges, "source_estimate")
  attr(reversed_edges, "expansion") <- attr(edges, "expansion")

  experimental <- if (task$experimental_boundary$state == "open") {
    .open_experimental_boundary(
      task$experimental_boundary$right_space,
      task$experimental_boundary$left_space
    )
  } else {
    query <- task$experimental_boundary$query
    .closed_experimental_boundary(
      pair_query(
        t(as.matrix(query$operator)), query$right_space, query$left_space,
        query$metadata
      ),
      task$experimental_boundary$role,
      task$experimental_boundary$sampling_axis
    )
  }
  neural <- if (task$neural_boundary$state == "closed") {
    if (task$neural_boundary$closure_kind == "bridge") {
      bridge <- task$neural_boundary$bridge
      reversed_bridge <- if (bridge$kind == "identity") {
        .identity_measurement_bridge(task$right_relation, task$left_relation)
      } else {
        reverse_bridge(bridge)
      }
      .closed_neural_boundary(
        reversed_bridge, task$right_relation, task$left_relation
      )
    } else {
      .closed_neural_query_boundary(
        .reverse_neural_pair_query(task$neural_boundary$query)
      )
    }
  } else {
    .open_neural_boundary(
      task$neural_boundary$right_frame,
      task$neural_boundary$left_frame,
      .reverse_measurement_edges(
        task$neural_boundary$edges,
        task$neural_boundary$left_frame,
        task$neural_boundary$right_frame
      )
    )
  }
  projection <- .reverse_query_identity(task$materialization$projection)
  materialization <- .evidence_materialization(
    task$materialization$kind,
    task$materialization$completeness,
    projection
  )
  stages <- .evidence_stage_plan(
    structure(task$stages$normalization$operation,
      class = "effect_edge_normalizer"),
    structure(task$stages$transform$operation,
      class = "effect_edge_transform"),
    structure(task$stages$reduction$operation,
      class = "effect_partition_reducer"),
    projection
  )
  compatibility <- if (task$identity_schema == "effect-form-v1") {
    .effect_task_semantic(
      task$right_relation_id, task$left_relation_id,
      task$spaces$experimental_right, task$spaces$experimental_left,
      reversed_edges, neural$bridge,
      structure(task$stages$normalization$operation,
        class = "effect_edge_normalizer"),
      structure(task$stages$transform$operation,
        class = "effect_edge_transform"),
      structure(task$stages$reduction$operation,
        class = "effect_partition_reducer"),
      projection
    )
  } else {
    NULL
  }
  .new_evidence_task(
    task$right_relation, task$left_relation, task$same_relation,
    reversed_edges, experimental, neural, stages, materialization,
    task$identity_schema, compatibility
  )
}

.new_effect_evidence_task <- function(
    left, right, same_relation, edges, bridge, normalizer, transform, reducer,
    query_identity, semantic) {
  experimental <- .open_experimental_boundary(
    left$effect_space, right$effect_space
  )
  neural <- .closed_neural_boundary(bridge, left, right)
  completeness <- if (identical(query_identity$kind, "complete_form")) {
    "complete_form"
  } else {
    "query_only"
  }
  materialization <- .evidence_materialization(
    "effect_form", completeness, query_identity
  )
  stages <- .evidence_stage_plan(
    normalizer, transform, reducer, query_identity
  )
  .new_evidence_task(
    left, right, same_relation, edges, experimental, neural, stages,
    materialization, "effect-form-v1", semantic
  )
}

.as_compiled_effect_task <- function(task, validate = TRUE) {
  if (isTRUE(validate)) .validate_evidence_task(task)
  if (task$identity_schema != "effect-form-v1" ||
      task$materialization$kind != "effect_form" ||
      task$experimental_boundary$state != "open" ||
      task$neural_boundary$state != "closed") {
    stop("Only the certified effect-form specialization has a legacy adapter.",
      call. = FALSE)
  }
  structure(list(
    left_relation = task$left_relation,
    right_relation = task$right_relation,
    same_relation = task$same_relation,
    left_relation_id = task$left_relation_id,
    right_relation_id = task$right_relation_id,
    left_space = task$spaces$experimental_left,
    right_space = task$spaces$experimental_right,
    ordered_edges = task$ordered_partition_products,
    bridge = task$neural_boundary$bridge,
    normalizer = structure(task$stages$normalization$operation,
      class = "effect_edge_normalizer"),
    transform = structure(task$stages$transform$operation,
      class = "effect_edge_transform"),
    lowering = task$stages$lowering,
    reducer = structure(task$stages$reduction$operation,
      class = "effect_partition_reducer"),
    query_identity = task$materialization$projection,
    source_uses = task$source_uses,
    distinct_handle_keys = task$distinct_handle_keys,
    task_id = task$task_id
  ), class = "effect_compiled_task")
}
