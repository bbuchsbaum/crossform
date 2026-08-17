# Boundary-typed evidence-task intermediate representation -----------------

.evidence_boundary_signature <- function(semantic) {
  .sha256_signature(semantic)
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
    .input_error(
      "A closed experimental boundary requires an axis-bound pair query."
    )
  }
  role <- match.arg(role)
  if (role == "variation" &&
      (!is.character(sampling_axis) || length(sampling_axis) != 1L ||
       is.na(sampling_axis) || !nzchar(sampling_axis))) {
    .input_error("A variation boundary requires one declared sampling axis.")
  }
  if (role == "effect" && !is.null(sampling_axis)) {
    .input_error("An effect boundary cannot claim a repeated-sampling axis.")
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
  if (!.is_finite_matrix(operator) || any(dim(operator) < 1L)) {
    .input_error(
      "A neural pair query must be a finite nonempty numeric matrix."
    )
  }
  operator <- unname(operator)
  storage.mode(operator) <- "double"
  left_domain <- .domain_reference(left_domain)
  right_domain <- .domain_reference(right_domain)
  if (nrow(operator) != left_domain$n_features ||
      ncol(operator) != right_domain$n_features) {
    .contract_error(
      "Neural pair-query dimensions do not match their identified domains."
    )
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
    signature = .sha256_signature(semantic)
  )), class = "effect_neural_pair_query")
}

.validate_neural_pair_query <- function(x) {
  expected <- c("operator", "left_domain", "right_domain", "provenance",
    "signature")
  if (!.sealed_fields(x, "effect_neural_pair_query", expected)) {
    .input_error("Neural pair-query fields are missing or noncanonical.")
  }
  rebuilt <- .neural_pair_query(
    x$operator, x$left_domain, x$right_domain, x$provenance
  )
  if (!identical(x, rebuilt)) {
    .contract_error("Neural pair-query identity is inconsistent.")
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
    .input_error("Evidence-boundary fields are missing or noncanonical.")
  }
  if (x$boundary == "experimental") {
    expected <- c("boundary", "state", "left_space", "right_space", "query",
      "role", "sampling_axis", "signature")
    if (!identical(names(x), expected)) {
      .input_error("Experimental-boundary fields are missing or noncanonical.")
    }
    rebuilt <- if (x$state == "open") {
      if (!is.null(x$query) || !is.null(x$role) || !is.null(x$sampling_axis)) {
        .input_error(
          "An open experimental boundary cannot carry a closure query."
        )
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
      .input_error("Neural-boundary fields are missing or noncanonical.")
    }
    if (x$state == "closed") {
      if (!.is_string(x$closure_kind, allow_empty = TRUE) ||
          !x$closure_kind %in% c("bridge", "query") ||
          !is.null(x$left_frame) || !is.null(x$right_frame) ||
          !is.null(x$edges)) {
        .input_error(
          "A closed neural boundary requires one bridge or query closure."
        )
      }
      if (x$closure_kind == "bridge") {
        if (is.null(x$bridge) || !is.null(x$query) ||
            !inherits(x$bridge, "effect_measurement_bridge")) {
          .input_error(
            "A bridge-closed neural boundary requires only its bridge."
          )
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
        .check_signature(
          x$signature, expected_signature,
          "Closed neural-boundary identity is inconsistent."
        )
        rebuilt <- x
      } else {
        if (!is.null(x$bridge) || is.null(x$query)) {
          .input_error(
            "A query-closed neural boundary requires only its query."
          )
        }
        rebuilt <- .closed_neural_query_boundary(x$query)
      }
    } else {
      if (!is.null(x$closure_kind) || !is.null(x$bridge) ||
          !is.null(x$query) || is.null(x$left_frame) ||
          is.null(x$right_frame) || is.null(x$edges)) {
        .input_error(
          "An open neural boundary requires two frames and explicit edges."
        )
      }
      rebuilt <- .open_neural_boundary(
        x$left_frame, x$right_frame, x$edges
      )
    }
  }
  if (!identical(x, rebuilt)) {
    .contract_error("Evidence-boundary identity is inconsistent.")
  }
  rebuilt
}

.evidence_materialization <- function(
    kind = c("effect_form", "measurement_form", "scalar_field"),
    completeness = c("complete_form", "query_only"), projection = NULL) {
  kind <- match.arg(kind)
  completeness <- match.arg(completeness)
  if (kind == "scalar_field" && completeness != "query_only") {
    .input_error("A scalar field is necessarily a query-only materialization.")
  }
  if (!is.null(projection) && !is.list(projection)) {
    .input_error("Materialization projection identity must be a list or NULL.")
  }
  semantic <- list(
    schema_version = 1L,
    kind = kind,
    completeness = completeness,
    projection = projection
  )
  structure(c(semantic[-1L], list(
    signature = .sha256_signature(semantic)
  )), class = "effect_evidence_materialization")
}

.validate_evidence_materialization <- function(x) {
  expected <- c("kind", "completeness", "projection", "signature")
  if (!.sealed_fields(x, "effect_evidence_materialization", expected)) {
    .input_error("Evidence-materialization fields are missing or noncanonical.")
  }
  rebuilt <- .evidence_materialization(
    x$kind, x$completeness, x$projection
  )
  if (!identical(x, rebuilt)) {
    .contract_error("Evidence-materialization identity is inconsistent.")
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
    signature = .sha256_signature(semantic)
  )), class = "effect_evidence_stages")
}

.validate_evidence_stage_plan <- function(x) {
  expected <- c("order", "partition_product", "normalization", "transform",
    "reduction", "projection", "lowering", "signature")
  if (!.sealed_fields(x, "effect_evidence_stages", expected)) {
    .input_error("Evidence-stage fields are missing or noncanonical.")
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
    .contract_error(
      "Evidence-stage identity or operation order is inconsistent."
    )
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
    .contract_error(
      "Requested materialization is incompatible with its open boundaries."
    )
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
  left_relation$capabilities <- .relation_source_capabilities(left_relation)
  right_relation$capabilities <- if (isTRUE(same_relation)) {
    left_relation$capabilities
  } else {
    .relation_source_capabilities(right_relation)
  }
  if (!.is_flag(same_relation)) {
    .input_error("Evidence-task `same_relation` must be TRUE or FALSE.")
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
    .input_error("Evidence-task boundaries occupy the wrong typed axes.")
  }
  if (!.same_effect_space(experimental_boundary$left_space,
        left_relation$effect_space) ||
      !.same_effect_space(experimental_boundary$right_space,
        right_relation$effect_space)) {
    .contract_error(
      "Experimental boundary axes differ from relation effect spaces."
    )
  }
  if (!.same_domain_reference(neural_boundary$left_space,
        left_relation$domain) ||
      !.same_domain_reference(neural_boundary$right_space,
        right_relation$domain)) {
    .contract_error("Neural boundary axes differ from relation neural spaces.")
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
      .input_error(
        "Effect-form compatibility identity requires its legacy semantic."
      )
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
      .contract_error(
        "Effect-form compatibility semantic differs from certified fields."
      )
    }
    semantic <- expected_compatibility
  } else {
    if (!is.null(compatibility_semantic)) {
      .input_error(
        "Native evidence tasks cannot carry a legacy compatibility semantic."
      )
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
  if (!.sealed_fields(task, "effect_evidence_task", expected) ||
      !identical(task$schema_version, 1L)) {
    .input_error("Evidence-task fields are missing or noncanonical.")
  }
  rebuilt <- .new_evidence_task(
    task$left_relation, task$right_relation, task$same_relation,
    task$ordered_partition_products, task$experimental_boundary,
    task$neural_boundary, task$stages, task$materialization,
    task$identity_schema,
    if (task$identity_schema == "effect-form-v1") task$semantic else NULL
  )
  if (!identical(task, rebuilt)) {
    .contract_error(
      "Evidence-task identity is inconsistent with its typed semantics."
    )
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
  .input_error("Unknown effect-form projection identity.")
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
    .input_error(
      "Only the certified effect-form specialization has a legacy adapter."
    )
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

# Effect-task construction and identity --------------------------------------
#
# Building a task and naming it is plan work, not execution work. These
# functions lived in `R/compiler.R`, which forced every plan that wanted a
# task identity to call up into the executor's file and made the plan and
# compiler layers mutually recursive.

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
  .sha256_signature(semantic)
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
  left$capabilities <- .relation_source_capabilities(left)
  right$capabilities <- if (same_relation) left$capabilities else
    .relation_source_capabilities(right)
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
      .contract_error(
        "Pair-query axes must match the ordered relation effect spaces."
      )
    }
    list(kind = "pair_query", query = query)
  } else if (inherits(query, "effect_query")) {
    .validate_query_for_compile(query)
    if (!identical(query$kind, "bilinear") || !isTRUE(query$fixed) ||
        !.same_effect_space(left$effect_space, right$effect_space) ||
        (!is.null(query$effect_space) &&
          !.same_effect_space(query$effect_space, left$effect_space))) {
      .input_error(
        "A bilinear query requires matching self-form relation axes."
      )
    }
    list(kind = "bilinear_compatibility", query = query)
  } else if (is.matrix(query) && is.numeric(query) &&
      all(is.finite(query)) && nrow(query) > 0L && ncol(query) > 0L) {
    if (!same_relation) {
      .input_error(
        "Raw physical queries are only compatible with self-form tasks."
      )
    }
    q <- length(left$effect_space$coordinates)
    if (nrow(query) != q * (q + 1L) / 2L) {
      .contract_error(
        "A raw self-form query must match the symmetric-packed width."
      )
    }
    list(kind = "physical_self_query", query = query)
  } else {
    .input_error(
      "Task queries must be NULL, axis-bound queries, or finite matrices."
    )
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
  if (!.sealed_fields(task, "effect_compiled_task", expected) ||
      !.is_flag(task$same_relation)) {
    .input_error("Compiled effect-task fields are missing or noncanonical.")
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
    .contract_error(
      "Compiled effect-task relation or axis identity is inconsistent."
    )
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
    .contract_error(
      "Compiled effect-task lowering is inconsistent with its operations."
    )
  }
  sources <- .effect_task_source_uses(
    task$left_relation, task$right_relation, left_id, right_id
  )
  if (!identical(task$source_uses, sources$uses) ||
      !identical(task$distinct_handle_keys, sources$distinct_handle_keys)) {
    .contract_error("Compiled effect-task source-use identity is inconsistent.")
  }
  semantic <- .effect_task_semantic(
    left_id, right_id, task$left_space, task$right_space,
    task$ordered_edges, task$bridge, task$normalizer, task$transform,
    task$reducer, task$query_identity
  )
  if (!identical(task$task_id, .effect_task_id(semantic))) {
    .contract_error(
      "Compiled effect-task identity is inconsistent with its semantics."
    )
  }
  invisible(task)
}

.effect_task_plan_id <- function(x, at, over, materialization, query, component,
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
  .sha256_signature(semantic)
}
