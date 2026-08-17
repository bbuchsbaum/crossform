# Durable measurement-form semantics ---------------------------------------

.measurement_regularization <- function(
    kind = c("none", "ridge"), lambda_left = 0, lambda_right = lambda_left,
    applied = FALSE) {
  kind <- match.arg(kind)
  values <- c(lambda_left = lambda_left, lambda_right = lambda_right)
  if (!.is_finite_numeric(values) || anyNA(values) || any(values < 0) ||
      !.is_flag(applied)) {
    .input_error("Regularization values and application state are invalid.")
  }
  if ((kind == "none" && (any(values != 0) || applied)) ||
      (kind == "ridge" && !any(values > 0))) {
    .contract_error(
      "Regularization kind, values, and application state disagree."
    )
  }
  semantic <- list(
    schema_version = 1L,
    kind = kind,
    lambda_left = as.numeric(lambda_left),
    lambda_right = as.numeric(lambda_right),
    applied = applied
  )
  structure(c(semantic[-1L], list(
    signature = .sha256_signature(semantic)
  )), class = "effect_measurement_regularization")
}

.validate_measurement_regularization <- function(x) {
  expected <- c("kind", "lambda_left", "lambda_right", "applied", "signature")
  if (!.sealed_fields(x, "effect_measurement_regularization", expected)) {
    .input_error("Measurement regularization is missing or noncanonical.")
  }
  rebuilt <- .measurement_regularization(
    x$kind, x$lambda_left, x$lambda_right, x$applied
  )
  if (!identical(x, rebuilt)) {
    .contract_error("Measurement regularization identity is inconsistent.")
  }
  rebuilt
}

.variation_pair_query <- function(H, effect_space, sampling_axis,
                                  construction = c("psd_variation",
                                                   "joint_covariance"),
                                  provenance = list(), tolerance = 1e-12) {
  construction <- match.arg(construction)
  effect_space <- .validate_effect_space(effect_space)
  if (!.is_finite_matrix(H) || nrow(H) != ncol(H) ||
      nrow(H) != length(effect_space$coordinates) ||
      max(abs(H - t(H))) > tolerance) {
    .input_error(
      "A variation query must be finite, self-adjoint, and axis-bound."
    )
  }
  eigenvalues <- eigen((H + t(H)) / 2,
    symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) < -tolerance) {
    .input_error(
      "A variation query construction must be positive semidefinite."
    )
  }
  .measurement_identifier(sampling_axis, "Variation-query `sampling_axis`")
  .validate_effect_provenance(provenance, "variation-query provenance")
  proof <- list(
    schema_version = 1L,
    role = "variation",
    sampling_axis = sampling_axis,
    construction = construction,
    effect_space = effect_space$signature,
    operator = unname(H),
    provenance = provenance
  )
  metadata <- list(evidence_capability = c(proof[-1L], list(
    signature = .sha256_signature(proof)
  )))
  pair_query(H, effect_space, effect_space, metadata)
}

.validate_variation_pair_query <- function(query, sampling_axis,
                                           construction, tolerance = 1e-12) {
  .validate_query_for_compile(query)
  capability <- query$metadata$evidence_capability
  if (!is.list(capability) || !identical(names(capability), c(
      "role", "sampling_axis", "construction", "effect_space", "operator",
      "provenance", "signature")) ||
      !identical(capability$role, "variation") ||
      !identical(capability$sampling_axis, sampling_axis) ||
      !(identical(capability$construction, construction) ||
        (identical(construction, "psd_variation") &&
         identical(capability$construction, "joint_covariance"))) ||
      !identical(capability$effect_space, query$left_space$signature) ||
      !identical(capability$operator, unname(as.matrix(query$operator)))) {
    .contract_error(
      "Variation-query construction provenance is absent or inconsistent."
    )
  }
  proof <- list(
    schema_version = 1L,
    role = capability$role,
    sampling_axis = capability$sampling_axis,
    construction = capability$construction,
    effect_space = capability$effect_space,
    operator = capability$operator,
    provenance = capability$provenance
  )
  expected <- .sha256_signature(proof)
  .check_signature(
    capability$signature, expected,
    "Variation-query construction signature is inconsistent."
  )
  h <- as.matrix(query$operator)
  if (!.same_effect_space(query$left_space, query$right_space) ||
      max(abs(h - t(h))) > tolerance ||
      min(eigen((h + t(h)) / 2, symmetric = TRUE,
        only.values = TRUE)$values) < -tolerance) {
    .input_error(
      "Variation-query construction does not establish a positive self form."
    )
  }
  invisible(query)
}

.measurement_frame_manifest <- function(frame) {
  frame <- .validate_measurement_frame(frame)
  legs <- lapply(frame$legs, function(leg) {
    list(
      signature = leg$signature,
      output_space = leg$output_space,
      support = leg$support,
      estimation = leg$estimation,
      decomposition = leg$decomposition
    )
  })
  list(
    signature = frame$signature,
    source_domain = frame$source_domain,
    node_ids = frame$node_ids,
    legs = legs,
    coverage = frame$coverage,
    injectivity = frame$injectivity,
    dual = frame$dual
  )
}

.validate_measurement_frame_manifest <- function(x) {
  expected <- c("signature", "source_domain", "node_ids", "legs", "coverage",
    "injectivity", "dual")
  if (!is.list(x) || !identical(names(x), expected) ||
      !.strong_sha256(x$signature) ||
      !identical(names(x$legs), x$node_ids) ||
      any(!vapply(x$legs, function(leg) {
        is.list(leg) && identical(names(leg), c(
          "signature", "output_space", "support", "estimation", "decomposition"
        )) && inherits(leg$output_space, "effect_measurement_axis") &&
          (is.null(leg$decomposition) ||
           inherits(leg$decomposition, "effect_measurement_decomposition"))
      }, logical(1)))) {
    .input_error("Measurement-frame manifest is missing or noncanonical.")
  }
  .validate_domain_reference(x$source_domain)
  lapply(x$legs, function(leg) {
    .validate_measurement_axis(leg$output_space)
    if (!is.null(leg$decomposition)) {
      decomposition <- .validate_measurement_decomposition(leg$decomposition)
      if (!identical(decomposition$output_space, leg$output_space)) {
        .contract_error("Frame-manifest decomposition and output axes differ.")
      }
    }
    invisible(leg)
  })
  x
}

.measurement_edge_set_signature <- function(edges, left_frame, right_frame,
                                            weighted) {
  .sha256_signature(list(
    schema_version = 1L,
    left_frame = left_frame,
    right_frame = right_frame,
    weighted = weighted,
    edges = unclass(edges)
  ))
}

.measurement_edges_are_complete <- function(edges, left_nodes, right_nodes) {
  expected <- as.vector(outer(left_nodes, right_nodes, paste, sep = "\r"))
  observed <- paste(edges$left, edges$right, sep = "\r")
  length(observed) == length(expected) && setequal(observed, expected)
}

.measurement_edges_reverse_closed <- function(edges) {
  observed <- paste(edges$left, edges$right, sep = "\r")
  reversed <- paste(edges$right, edges$left, sep = "\r")
  setequal(observed, reversed)
}

.measurement_plan_semantic <- function(fields) {
  list(
    schema_version = 1L,
    origin_task_id = fields$origin_task_id,
    left_relation_id = fields$left_relation_id,
    right_relation_id = fields$right_relation_id,
    source_identity = unclass(fields$source_uses[, c(
      "side", "relation_id", "partition", "stable_revision"
    ), drop = FALSE]),
    same_relation = fields$same_relation,
    left_effect_space = fields$left_effect_space$signature,
    right_effect_space = fields$right_effect_space$signature,
    left_neural_space = fields$left_neural_space$signature,
    right_neural_space = fields$right_neural_space$signature,
    experimental_query = fields$experimental_query,
    query_role = fields$query_role,
    sampling_axis = fields$sampling_axis,
    query_construction = fields$query_construction,
    left_frame = fields$left_frame$signature,
    right_frame = fields$right_frame$signature,
    edge_set_signature = fields$edge_set_signature,
    block_index_signature = attr(fields$block_index, "signature", exact = TRUE),
    partition_products = unclass(fields$partition_products),
    partition_expansion = fields$partition_expansion,
    stages = fields$stages,
    regularization = fields$regularization,
    edge_scope = fields$edge_scope
  )
}

.new_measurement_plan_from_fields <- function(fields) {
  semantic <- .measurement_plan_semantic(fields)
  scientific_plan_id <- .sha256_signature(semantic, "measurement-sha256:")
  structure(c(fields, list(
    scientific_plan_id = scientific_plan_id,
    signature = scientific_plan_id
  )), class = "effect_measurement_plan")
}

.measurement_result_plan <- function(
    task, query_construction = c("arbitrary", "psd_variation",
                                 "joint_covariance"),
    edge_scope = c("requested", "frame_complete"),
    regularization = .measurement_regularization(), tolerance = 1e-12) {
  task <- .validate_raw_measurement_task(task)
  query_construction <- match.arg(query_construction)
  edge_scope <- match.arg(edge_scope)
  regularization <- .validate_measurement_regularization(regularization)
  if (regularization$applied) {
    .input_error("Raw measurement blocks cannot claim applied regularization.")
  }
  query <- task$experimental_boundary$query
  if (query_construction != "arbitrary") {
    if (!identical(task$experimental_boundary$role, "variation") ||
        is.null(task$experimental_boundary$sampling_axis)) {
      .input_error(
        "Positive variation capabilities require a variation boundary."
      )
    }
    .validate_variation_pair_query(
      query, task$experimental_boundary$sampling_axis,
      query_construction, tolerance
    )
  }
  if (query_construction == "joint_covariance" &&
      (!isTRUE(task$same_relation) ||
       any(task$ordered_partition_products$left !=
         task$ordered_partition_products$right))) {
    .input_error(paste0(
      "A joint-covariance claim requires same-relation, same-partition ",
      "self-products; crossvalidated partition forms can be indefinite."
    ))
  }
  left_frame <- task$neural_boundary$left_frame
  right_frame <- task$neural_boundary$right_frame
  edge_set <- task$neural_boundary$edges
  if (edge_scope == "frame_complete" &&
      !.measurement_edges_are_complete(
        edge_set$edges, left_frame$node_ids, right_frame$node_ids
      )) {
    .input_error(
      "A frame-complete claim requires every explicit ordered node pair."
    )
  }
  block_index <- .measurement_block_index(
    edge_set, left_frame, right_frame
  )
  fields <- list(
    origin_task_id = task$task_id,
    left_relation_id = task$left_relation_id,
    right_relation_id = task$right_relation_id,
    source_uses = task$source_uses,
    same_relation = task$same_relation,
    left_effect_space = task$spaces$experimental_left,
    right_effect_space = task$spaces$experimental_right,
    left_neural_space = task$spaces$neural_left,
    right_neural_space = task$spaces$neural_right,
    experimental_query = query,
    query_role = task$experimental_boundary$role,
    sampling_axis = task$experimental_boundary$sampling_axis,
    query_construction = query_construction,
    left_frame = .measurement_frame_manifest(left_frame),
    right_frame = .measurement_frame_manifest(right_frame),
    edges = edge_set$edges,
    edges_weighted = edge_set$weighted,
    edge_set_signature = edge_set$signature,
    block_index = block_index,
    partition_products = task$ordered_partition_products,
    partition_expansion = attr(task$ordered_partition_products, "expansion"),
    stages = task$stages,
    regularization = regularization,
    edge_scope = edge_scope,
    tolerance = tolerance
  )
  .new_measurement_plan_from_fields(fields)
}

.validate_measurement_plan <- function(plan) {
  expected <- c(
    "origin_task_id", "left_relation_id", "right_relation_id", "source_uses",
    "same_relation",
    "left_effect_space", "right_effect_space", "left_neural_space",
    "right_neural_space", "experimental_query", "query_role", "sampling_axis",
    "query_construction", "left_frame", "right_frame", "edges",
    "edges_weighted", "edge_set_signature", "block_index",
    "partition_products", "partition_expansion", "stages", "regularization",
    "edge_scope", "tolerance", "scientific_plan_id", "signature"
  )
  if (!.sealed_fields(plan, "effect_measurement_plan", expected) ||
      !.is_flag(plan$same_relation) ||
      !plan$query_construction %in% c("arbitrary", "psd_variation", "joint_covariance") ||
      !plan$edge_scope %in% c("requested", "frame_complete") ||
      !is.data.frame(plan$edges)) {
    .input_error("Measurement scientific plan is missing or noncanonical.")
  }
  expected_source_names <- c("side", "relation_id", "partition",
    "stable_revision", "handle_key")
  if (!is.data.frame(plan$source_uses) ||
      !identical(names(plan$source_uses), expected_source_names) ||
      nrow(plan$source_uses) < 2L ||
      anyNA(plan$source_uses[, c("side", "relation_id", "partition",
        "stable_revision"), drop = FALSE]) ||
      any(!plan$source_uses$side %in% c("left", "right")) ||
      any(!vapply(plan$source_uses$stable_revision, .strong_sha256,
        logical(1))) ||
      any(plan$source_uses$relation_id[plan$source_uses$side == "left"] !=
        plan$left_relation_id) ||
      any(plan$source_uses$relation_id[plan$source_uses$side == "right"] !=
        plan$right_relation_id)) {
    .input_error(
      "Measurement plan source-use ledger is missing or inconsistent."
    )
  }
  .validate_effect_space(plan$left_effect_space)
  .validate_effect_space(plan$right_effect_space)
  .validate_domain_reference(plan$left_neural_space)
  .validate_domain_reference(plan$right_neural_space)
  .validate_query_for_compile(plan$experimental_query)
  left_frame <- .validate_measurement_frame_manifest(plan$left_frame)
  right_frame <- .validate_measurement_frame_manifest(plan$right_frame)
  regularization <- .validate_measurement_regularization(plan$regularization)
  index <- .validate_measurement_block_index(plan$block_index)
  expected_edge_signature <- .measurement_edge_set_signature(
    plan$edges, left_frame$signature, right_frame$signature,
    plan$edges_weighted
  )
  if (!identical(plan$edge_set_signature, expected_edge_signature) ||
      !identical(attr(index, "edge_set", exact = TRUE),
        plan$edge_set_signature) ||
      !identical(attr(index, "left_frame", exact = TRUE),
        left_frame$signature) ||
      !identical(attr(index, "right_frame", exact = TRUE),
        right_frame$signature) ||
      !identical(index$left, plan$edges$left) ||
      !identical(index$right, plan$edges$right)) {
    .contract_error("Measurement plan edges, frames, and block index disagree.")
  }
  if (plan$edge_scope == "frame_complete" &&
      !.measurement_edges_are_complete(
        plan$edges, left_frame$node_ids, right_frame$node_ids
      )) {
    .input_error("Measurement plan falsely claims a complete frame edge set.")
  }
  if (plan$query_construction != "arbitrary") {
    .validate_variation_pair_query(
      plan$experimental_query, plan$sampling_axis,
      plan$query_construction, plan$tolerance
    )
  }
  fields <- plan[setdiff(names(plan), c("scientific_plan_id", "signature"))]
  rebuilt <- .new_measurement_plan_from_fields(fields)
  if (!identical(plan, rebuilt)) {
    .contract_error("Measurement scientific-plan identity is inconsistent.")
  }
  invisible(plan)
}

.measurement_form_capabilities <- function(plan) {
  .validate_measurement_plan(plan)
  h <- as.matrix(plan$experimental_query$operator)
  query_self <- .same_effect_space(
    plan$left_effect_space, plan$right_effect_space
  ) && nrow(h) == ncol(h)
  query_symmetric <- query_self && max(abs(h - t(h))) <= plan$tolerance
  frame_self <- identical(plan$left_frame$signature,
    plan$right_frame$signature)
  self_form <- isTRUE(plan$same_relation) && query_self && frame_self
  symmetric <- self_form && query_symmetric &&
    .measurement_edges_reverse_closed(plan$edges)
  repeated_variation <- plan$query_construction != "arbitrary"
  joint_covariance <- plan$query_construction == "joint_covariance" &&
    isTRUE(plan$same_relation) &&
    all(plan$partition_products$left == plan$partition_products$right)
  complete_edge_set <- plan$edge_scope == "frame_complete"
  guaranteed_psd <- symmetric && joint_covariance && complete_edge_set
  semantic <- list(
    schema_version = 1L,
    plan = plan$signature,
    self_form = self_form,
    symmetric = symmetric,
    guaranteed_psd = guaranteed_psd,
    repeated_variation = repeated_variation,
    joint_covariance = joint_covariance,
    complete_edge_set = complete_edge_set,
    sampling_axis = if (repeated_variation) plan$sampling_axis else NULL,
    construction = plan$query_construction
  )
  structure(c(semantic[-1L], list(
    signature = .sha256_signature(semantic)
  )), class = "effect_measurement_capabilities")
}

.validate_measurement_capabilities <- function(x, plan) {
  expected <- c("plan", "self_form", "symmetric", "guaranteed_psd",
    "repeated_variation", "joint_covariance", "complete_edge_set",
    "sampling_axis", "construction", "signature")
  if (!.sealed_fields(x, "effect_measurement_capabilities", expected)) {
    .input_error("Measurement-form capabilities are missing or noncanonical.")
  }
  rebuilt <- .measurement_form_capabilities(plan)
  if (!identical(x, rebuilt)) {
    .contract_error("Measurement-form capabilities cannot be forged.")
  }
  rebuilt
}

.measurement_block_diagnostics <- function(blocks, h, regularization,
                                           tolerance = 1e-10) {
  regularization <- .validate_measurement_regularization(regularization)
  h_singular <- svd(as.matrix(h), nu = 0L, nv = 0L)$d
  h_threshold <- if (length(h_singular)) max(h_singular) * tolerance else 0
  rows <- lapply(seq_along(blocks), function(position) {
    value <- blocks[[position]]
    singular <- svd(value, nu = 0L, nv = 0L)$d
    threshold <- if (length(singular)) max(singular) * tolerance else 0
    retained <- singular[singular > threshold]
    symmetric <- nrow(value) == ncol(value) &&
      max(abs(value - t(value))) <= tolerance
    eigenvalues <- if (symmetric) {
      eigen((value + t(value)) / 2, symmetric = TRUE,
        only.values = TRUE)$values
    } else {
      numeric()
    }
    data.frame(
      edge_id = names(blocks)[[position]],
      effective_rank = as.integer(length(retained)),
      condition_number = if (!length(retained)) Inf else
        max(retained) / min(retained),
      observed_min_eigenvalue = if (length(eigenvalues)) min(eigenvalues) else
        NA_real_,
      observed_max_eigenvalue = if (length(eigenvalues)) max(eigenvalues) else
        NA_real_,
      symmetric_observed = symmetric,
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  semantic <- list(
    schema_version = 1L,
    tolerance = tolerance,
    method = "svd-relative-v1",
    experimental_effective_rank = as.integer(sum(h_singular > h_threshold)),
    blocks = unclass(table),
    regularization = list(
      policy = regularization$signature,
      applied = regularization$applied,
      lambda_left = regularization$lambda_left,
      lambda_right = regularization$lambda_right
    )
  )
  structure(list(
    tolerance = semantic$tolerance,
    method = semantic$method,
    experimental_effective_rank = semantic$experimental_effective_rank,
    blocks = table,
    regularization = semantic$regularization,
    signature = .sha256_signature(semantic)
  ), class = "effect_measurement_diagnostics")
}

.validate_measurement_diagnostics <- function(x, plan) {
  expected <- c("tolerance", "method", "experimental_effective_rank", "blocks",
    "regularization", "signature")
  if (!.sealed_fields(x, "effect_measurement_diagnostics", expected) ||
      !identical(x$method, "svd-relative-v1") || !is.data.frame(x$blocks) ||
      !identical(x$blocks$edge_id, plan$block_index$edge_id)) {
    .input_error("Measurement diagnostics are missing or noncanonical.")
  }
  semantic <- list(
    schema_version = 1L,
    tolerance = x$tolerance,
    method = x$method,
    experimental_effective_rank = x$experimental_effective_rank,
    blocks = unclass(x$blocks),
    regularization = x$regularization
  )
  expected_signature <- .sha256_signature(semantic)
  .check_signature(
    x$signature, expected_signature,
    "Measurement diagnostics signature is inconsistent."
  )
  invisible(x)
}

.measurement_execution_receipt <- function(task, contraction, plan) {
  task <- .validate_raw_measurement_task(task)
  plan <- .validate_measurement_plan(plan)
  if (!inherits(contraction, "effect_measurement_contraction") ||
      !identical(contraction$task_id, task$task_id)) {
    .contract_error(
      "Measurement contraction and scientific task identities differ."
    )
  }
  sources <- c(task$left_relation$capabilities,
    if (task$same_relation) list() else task$right_relation$capabilities)
  observed <- .empty_execution_observations()
  count <- nrow(task$ordered_partition_products)
  observed$task_counts <- c(
    planned = count, started = count, completed = count, failed = 0, retried = 0
  )
  observed$features_completed <-
    task$left_relation$n_features + if (task$same_relation) 0 else
      task$right_relation$n_features
  observed$bytes_read <- contraction$diagnostics$source_bytes_read
  observed$tiles$feature_block <- contraction$diagnostics$feature_block
  observed$tiles$coordinate_tile <- contraction$diagnostics$edge_tile
  observed$stage_seconds <- c(
    measurement_contraction = contraction$diagnostics$elapsed_seconds
  )
  execution <- execution_receipt(
    scientific_plan_id = plan$scientific_plan_id,
    compute = contraction$compute,
    sources = sources,
    memory = contraction$memory$plan,
    kernel_version = "measurement-form-v1",
    task_partition_id = paste0("features-",
      contraction$diagnostics$feature_block),
    reduction_plan_id = task$stages$signature,
    completion_status = "complete",
    task_count = count,
    completed_task_count = count,
    elapsed_seconds = contraction$diagnostics$elapsed_seconds,
    observed = observed
  )
  semantic <- list(
    schema_version = 1L,
    scientific_plan_id = plan$scientific_plan_id,
    plan_signature = plan$signature,
    execution = execution,
    route = contraction$route_plan$signature,
    derivation = list(kind = "executed", parent = NULL)
  )
  structure(c(semantic[-1L], list(
    signature = .sha256_signature(semantic)
  )), class = "effect_measurement_receipt")
}

.validate_measurement_receipt <- function(x, plan) {
  expected <- c("scientific_plan_id", "plan_signature", "execution", "route",
    "derivation", "signature")
  if (!.sealed_fields(x, "effect_measurement_receipt", expected) ||
      !identical(x$scientific_plan_id, plan$scientific_plan_id) ||
      !identical(x$plan_signature, plan$signature) || !is.list(x$derivation) ||
      !identical(names(x$derivation), c("kind", "parent")) ||
      !x$derivation$kind %in% c("executed", "reversal")) {
    .contract_error(
      "Measurement receipt is missing or inconsistent with its plan."
    )
  }
  .validate_execution_receipt(x$execution)
  if (x$derivation$kind == "executed" &&
      (!is.null(x$derivation$parent) ||
       !identical(x$execution$scientific_plan_id,
         plan$scientific_plan_id))) {
    .contract_error(
      "Executed measurement receipt has contradictory provenance."
    )
  }
  if (x$derivation$kind == "reversal" &&
      (!is.character(x$derivation$parent) ||
       length(x$derivation$parent) != 1L)) {
    .contract_error("Reversed measurement receipt lacks its parent receipt.")
  }
  semantic <- list(
    schema_version = 1L,
    scientific_plan_id = x$scientific_plan_id,
    plan_signature = x$plan_signature,
    execution = x$execution,
    route = x$route,
    derivation = x$derivation
  )
  expected_signature <- .sha256_signature(semantic)
  .check_signature(
    x$signature, expected_signature,
    "Measurement receipt signature is inconsistent."
  )
  invisible(x)
}

.measurement_contract_signature <- function(plan, capabilities, completeness) {
  .sha256_signature(list(
    schema_version = 1L,
    type = "measurement_form",
    plan = plan$signature,
    capabilities = capabilities$signature,
    completeness = completeness
  ))
}

.new_measurement_form <- function(store, plan, capabilities, diagnostics,
                                  receipt) {
  plan <- .validate_measurement_plan(plan)
  capabilities <- .validate_measurement_capabilities(capabilities, plan)
  .validate_measurement_diagnostics(diagnostics, plan)
  .validate_measurement_receipt(receipt, plan)
  .validate_measurement_store(store, require_complete = TRUE, probe = TRUE)
  if (!identical(store$index, plan$block_index)) {
    .contract_error(
      "Measurement store and scientific plan block indices differ."
    )
  }
  completeness <- "complete"
  structure(list(
    store = store,
    block_index = plan$block_index,
    left_frame = plan$left_frame,
    right_frame = plan$right_frame,
    plan = plan,
    capabilities = capabilities,
    diagnostics = diagnostics,
    receipt = receipt,
    codec = store$codec,
    storage = store$representation,
    edge_completeness = if (capabilities$complete_edge_set) {
      "frame_complete"
    } else {
      "requested_complete"
    },
    contract_signature = .measurement_contract_signature(
      plan, capabilities, completeness
    ),
    # `completeness` reads "complete" here and "full" on an `effect_form`
    # (R/result.R). The word is hashed into `$contract_signature` just above,
    # so harmonizing the two would invalidate every recorded measurement-form
    # identity. `$result_capability` is the field with one vocabulary across
    # both result kinds; branch on that.
    result_capability = "complete_form",
    completeness = completeness
  ), class = "effect_measurement_form")
}

.measurement_form_from_contraction <- function(
    task, contraction, storage = c("memory", "block"), path = NULL,
    query_construction = c("arbitrary", "psd_variation", "joint_covariance"),
    edge_scope = c("requested", "frame_complete"),
    regularization = .measurement_regularization(), tolerance = 1e-10) {
  storage <- match.arg(storage)
  query_construction <- match.arg(query_construction)
  edge_scope <- match.arg(edge_scope)
  if (!inherits(contraction, "effect_measurement_contraction") ||
      !identical(contraction$task_id, task$task_id)) {
    .contract_error(
      "Measurement contraction does not belong to the requested task."
    )
  }
  plan <- .measurement_result_plan(
    task, query_construction, edge_scope, regularization, tolerance
  )
  index <- plan$block_index
  store <- if (storage == "memory") {
    .memory_measurement_store(contraction$blocks, index)
  } else {
    if (is.null(path)) {
      .input_error("Block-backed measurement results require an explicit path.")
    }
    value <- .file_measurement_store(path, index, create = TRUE)
    for (edge in seq_along(contraction$blocks)) {
      value <- .write_measurement_block(value, edge,
        contraction$blocks[[edge]])
    }
    value
  }
  diagnostics <- .measurement_block_diagnostics(
    contraction$blocks, task$experimental_boundary$query$operator,
    regularization, tolerance
  )
  capabilities <- .measurement_form_capabilities(plan)
  receipt <- .measurement_execution_receipt(task, contraction, plan)
  .new_measurement_form(store, plan, capabilities, diagnostics, receipt)
}

.validate_measurement_form <- function(x, probe = TRUE) {
  expected <- c("store", "block_index", "left_frame", "right_frame", "plan",
    "capabilities", "diagnostics", "receipt", "codec", "storage",
    "edge_completeness", "contract_signature", "result_capability",
    "completeness")
  if (!.sealed_fields(x, "effect_measurement_form", expected) ||
      inherits(x, "effect_form") ||
      !identical(x$result_capability, "complete_form") ||
      !identical(x$completeness, "complete")) {
    .input_error("`x` must be a canonical complete measurement_form.")
  }
  plan <- .validate_measurement_plan(x$plan)
  capabilities <- .validate_measurement_capabilities(x$capabilities, plan)
  .validate_measurement_diagnostics(x$diagnostics, plan)
  .validate_measurement_receipt(x$receipt, plan)
  if (!identical(x$block_index, plan$block_index) ||
      !identical(x$left_frame, plan$left_frame) ||
      !identical(x$right_frame, plan$right_frame) ||
      !identical(x$codec, "measurement-block-double-v1") ||
      !identical(x$storage, x$store$representation) ||
      !identical(x$edge_completeness,
        if (capabilities$complete_edge_set) "frame_complete" else
          "requested_complete") ||
      !identical(x$contract_signature,
        .measurement_contract_signature(plan, capabilities, "complete"))) {
    .contract_error(
      "Measurement-form axes, completeness, storage, or identity disagree."
    )
  }
  .validate_measurement_store(x$store, require_complete = TRUE, probe = probe)
  if (!identical(x$store$index, x$block_index)) {
    .contract_error(
      "Measurement-form store index is incompatible with its axes."
    )
  }
  invisible(x)
}

.measurement_block <- function(x, edge) {
  .validate_measurement_form(x, probe = FALSE)
  .read_measurement_store(x$store, edge)
}

.reverse_measurement_plan <- function(plan) {
  .validate_measurement_plan(plan)
  fields <- plan[setdiff(names(plan), c("scientific_plan_id", "signature"))]
  swap <- function(left, right) list(left = right, right = left)
  relation <- swap(fields$left_relation_id, fields$right_relation_id)
  effects <- swap(fields$left_effect_space, fields$right_effect_space)
  neural <- swap(fields$left_neural_space, fields$right_neural_space)
  frames <- swap(fields$left_frame, fields$right_frame)
  fields$left_relation_id <- relation$left
  fields$right_relation_id <- relation$right
  source_uses <- fields$source_uses
  left_uses <- source_uses[source_uses$side == "right", , drop = FALSE]
  right_uses <- source_uses[source_uses$side == "left", , drop = FALSE]
  left_uses$side <- "left"
  right_uses$side <- "right"
  fields$source_uses <- rbind(left_uses, right_uses)
  rownames(fields$source_uses) <- NULL
  fields$left_effect_space <- effects$left
  fields$right_effect_space <- effects$right
  fields$left_neural_space <- neural$left
  fields$right_neural_space <- neural$right
  fields$left_frame <- frames$left
  fields$right_frame <- frames$right
  query <- plan$experimental_query
  fields$experimental_query <- pair_query(
    t(as.matrix(query$operator)), query$right_space, query$left_space,
    query$metadata
  )
  reversed_edges <- fields$edges
  reversed_edges$left <- fields$edges$right
  reversed_edges$right <- fields$edges$left
  fields$edges <- reversed_edges
  fields$edge_set_signature <- .measurement_edge_set_signature(
    reversed_edges, fields$left_frame$signature,
    fields$right_frame$signature, fields$edges_weighted
  )
  index <- fields$block_index
  reversed_index <- data.frame(
    edge_id = index$edge_id,
    left = index$right,
    right = index$left,
    weight = index$weight,
    left_axis = index$right_axis,
    right_axis = index$left_axis,
    d_left = index$d_right,
    d_right = index$d_left,
    offset_elements = index$offset_elements,
    length_elements = index$length_elements,
    stringsAsFactors = FALSE
  )
  semantic <- list(
    schema_version = 1L,
    edge_set = fields$edge_set_signature,
    left_frame = fields$left_frame$signature,
    right_frame = fields$right_frame$signature,
    index = lapply(reversed_index, unname)
  )
  fields$block_index <- structure(reversed_index,
    edge_set = fields$edge_set_signature,
    left_frame = fields$left_frame$signature,
    right_frame = fields$right_frame$signature,
    signature = .sha256_signature(semantic),
    class = c("effect_measurement_block_index", "data.frame")
  )
  products <- fields$partition_products
  products$left <- fields$partition_products$right
  products$right <- fields$partition_products$left
  fields$partition_products <- products
  fields$regularization <- .measurement_regularization(
    fields$regularization$kind,
    fields$regularization$lambda_right,
    fields$regularization$lambda_left,
    fields$regularization$applied
  )
  .new_measurement_plan_from_fields(fields)
}

.reverse_measurement_form <- function(x) {
  .validate_measurement_form(x)
  reversed_plan <- .reverse_measurement_plan(x$plan)
  blocks <- stats::setNames(lapply(seq_len(nrow(x$block_index)), function(edge) {
    t(.measurement_block(x, edge))
  }), reversed_plan$block_index$edge_id)
  store <- .memory_measurement_store(blocks, reversed_plan$block_index)
  capabilities <- .measurement_form_capabilities(reversed_plan)
  diagnostics <- .measurement_block_diagnostics(
    blocks, reversed_plan$experimental_query$operator,
    reversed_plan$regularization, x$diagnostics$tolerance
  )
  receipt <- x$receipt
  semantic <- list(
    schema_version = 1L,
    scientific_plan_id = reversed_plan$scientific_plan_id,
    plan_signature = reversed_plan$signature,
    execution = receipt$execution,
    route = receipt$route,
    derivation = list(kind = "reversal", parent = receipt$signature)
  )
  reversed_receipt <- structure(c(semantic[-1L], list(
    signature = .sha256_signature(semantic)
  )), class = "effect_measurement_receipt")
  .new_measurement_form(
    store, reversed_plan, capabilities, diagnostics, reversed_receipt
  )
}

.measurement_view <- function(values, plan, receipt, view) {
  plan <- .validate_measurement_plan(plan)
  .validate_measurement_receipt(receipt, plan)
  if (!.is_finite_matrix(values) || nrow(values) != nrow(plan$block_index) ||
      ncol(values) < 1L || !is.list(view)) {
    .input_error("Measurement view values or query description are invalid.")
  }
  semantic <- list(
    schema_version = 1L,
    plan = plan$signature,
    receipt = receipt$signature,
    view = view,
    dimensions = dim(values)
  )
  structure(list(
    values = values,
    edge_index = plan$block_index,
    view = view,
    plan = plan,
    receipt = receipt,
    edge_completeness = if (plan$edge_scope == "frame_complete") {
      "frame_complete"
    } else {
      "requested_complete"
    },
    contract_signature = .sha256_signature(semantic),
    result_capability = "query_only",
    completeness = "query_only"
  ), class = "effect_measurement_view")
}

.validate_measurement_view <- function(x) {
  expected <- c("values", "edge_index", "view", "plan", "receipt",
    "edge_completeness", "contract_signature", "result_capability",
    "completeness")
  if (!.sealed_fields(x, "effect_measurement_view", expected) ||
      inherits(x, "effect_measurement_form") ||
      !identical(x$result_capability, "query_only") ||
      !identical(x$completeness, "query_only")) {
    .input_error("`x` must be a canonical query-only measurement view.")
  }
  plan <- .validate_measurement_plan(x$plan)
  .validate_measurement_receipt(x$receipt, plan)
  if (!identical(x$edge_index, plan$block_index) ||
      !identical(x$edge_completeness, if (plan$edge_scope == "frame_complete") "frame_complete" else "requested_complete") ||
      !.is_finite_matrix(x$values) || nrow(x$values) != nrow(x$edge_index)) {
    .input_error("Measurement view axes or values are inconsistent.")
  }
  semantic <- list(
    schema_version = 1L,
    plan = plan$signature,
    receipt = x$receipt$signature,
    view = x$view,
    dimensions = dim(x$values)
  )
  expected_signature <- .sha256_signature(semantic)
  .check_signature(
    x$contract_signature, expected_signature,
    "Measurement-view identity is inconsistent."
  )
  invisible(x)
}

.require_measurement_view_capabilities <- function(
    x, view = c("effect_coupling", "covariance_coupling",
                "canonical_coupling"),
    positive_self_blocks = FALSE) {
  view <- match.arg(view)
  .validate_measurement_form(x, probe = FALSE)
  .check_flag(positive_self_blocks, "positive_self_blocks")
  if (view == "effect_coupling") return(invisible(x))
  # Capability gating, not shape checking: these are contract-level refusals,
  # so callers branch on `$capability` rather than on the prose.
  if (!isTRUE(x$capabilities$repeated_variation)) {
    .capability_refusal(paste0(
      "Covariance coupling requires certified repeated variation: a ",
      "covariance is an average over repeats, and this form has not ",
      "established that its measurement axis repeats."
    ),
      capability = "certified_repeated_variation",
      namespace = "coupling_views",
      reasons = "repeated_variation_not_certified",
      remedies = paste0(
        "Build the form with a `variation_query()` that declares its ",
        "sampling axis, over a pairing whose partitions repeat that axis."
      )
    )
  }
  if (view == "canonical_coupling" &&
      (!isTRUE(x$capabilities$joint_covariance) ||
       !isTRUE(positive_self_blocks))) {
    .capability_refusal(paste0(
      "Canonical coupling requires a coherent joint covariance and ",
      "separately validated positive self-blocks."
    ),
      capability = "coherent_joint_covariance",
      namespace = "coupling_views",
      reasons = c(
        if (!isTRUE(x$capabilities$joint_covariance)) {
          "joint_covariance_not_certified"
        },
        if (!isTRUE(positive_self_blocks)) "self_blocks_not_validated"
      ),
      remedies = paste0(
        "Build the form over a pairing that certifies a joint covariance ",
        "across both measurement spaces, and read it through ",
        "`canonical_coupling()`, which validates the self-blocks it needs."
      )
    )
  }
  invisible(x)
}
