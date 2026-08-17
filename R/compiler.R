# Geometry compiler ----------------------------------------------------------
#
# Layer 4. The compiler consumes an `effect_geometry_plan` and lowers it to
# an `effect_geometry_execution_plan`: it validates and packs the query,
# resolves which components the requested view needs, selects tiles and a
# memory plan, names the kernel, and hashes the result back to the parent
# plan. It never builds a plan — `R/geometry-entry.R` does that for the
# relation compatibility forms — and it never executes one;
# `R/execution-driver.R` is the runtime.

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
        .contract_error(
          "Pair-query axes must match the ordered relation effect spaces."
        )
      }
      value <- matrix(as.vector(as.matrix(query$operator)), ncol = 1L)
      colnames(value) <- "view1"
      return(value)
    }
    .input_error(paste0(
      "Rectangular plans execute axis-bound `pair_query()` operators; ",
      "materialize with `materialize_geometry()` and use `query_geometry()` ",
      "for raw physical view matrices."
    ))
  }
  packed_width <- q * (q + 1L) / 2L
  if (.is_pair_difference_query(query)) {
    if (!identical(query$effects, effect_space$coordinates)) {
      .contract_error("The query and relation effect spaces are incompatible.")
    }
    return(query)
  }
  if (inherits(query, "effect_query")) {
    .validate_query_for_compile(query)
    if (!identical(query$kind, "bilinear") || !isTRUE(query$fixed)) {
      .input_error("Direct execution requires a fixed bilinear query.")
    }
    if (nrow(query$operator) != q) {
      .input_error(
        "The query operator dimension must equal the effect dimension."
      )
    }
    if (!is.null(query$effect_space) &&
        !.same_effect_space(query$effect_space, effect_space)) {
      .contract_error("The query and relation effect spaces are incompatible.")
    }
    value <- matrix(.svec_symmetric(query$operator), ncol = 1L)
    colnames(value) <- "view1"
    return(value)
  }
  if (!.is_finite_matrix(query) || nrow(query) != packed_width ||
      ncol(query) < 1L) {
    .input_error("`query` must be a finite packed-coordinate-by-view matrix.")
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

# Lowering a geometry plan to an execution plan ------------------------------
#
# A geometry plan is pure data: it names an estimand. This is where the
# compiler consumes one and decides how it will be executed — the physical
# query, the component requirements, the tile sizes and memory plan, the
# kernel version, and the execution identity that binds all of it back to the
# parent plan. It constructs no plan of its own.

.geometry_execution_signature <- function(fields) {
  .sha256_signature(list(
    schema_version = 1L,
    parent = fields$parent_signature,
    scientific_plan_id = fields$scientific_plan_id,
    task = fields$task$task_id,
    storage = fields$storage,
    storage_path = fields$storage_path,
    component = fields$component,
    signed_query = fields$signed_query,
    requirements = fields$requirements,
    output_width = fields$output_width,
    feature_block = fields$feature_block,
    row_tile = fields$row_tile,
    coordinate_tile = fields$coordinate_tile,
    memory = unclass(fields$memory),
    lowering = fields$lowering,
    kernel_version = fields$kernel_version
  ))
}

.geometry_kernel_version <- function(support_streamed, explicit_metric,
                                     query_fused, structured_query) {
  if (isTRUE(support_streamed)) {
    "support-streamed-metric-v1"
  } else if (isTRUE(explicit_metric)) {
    "additive-fixed-metric-v1"
  } else if (isTRUE(query_fused) && isTRUE(structured_query)) {
    "additive-query-fused-native-v1"
  } else if (isTRUE(query_fused)) {
    "additive-query-fused-v2"
  } else {
    "additive-packed-native-v1"
  }
}

.compile_geometry_execution_plan <- function(
    plan, query = NULL, component = NULL,
    storage = c("memory", "block"), storage_path = NULL,
    signed_query = NULL) {
  .validate_geometry_plan(plan)
  storage <- match.arg(storage)
  rectangular <- identical(plan$codec, "rectangular")
  if (!is.null(query)) {
    physical_query <- .compiler_query(
      query, plan$task$left_relation$effect_space,
      right_space = plan$task$right_relation$effect_space,
      codec = plan$codec
    )
    component <- match.arg(component,
      c("total", "coherent", "configuration", "contrast"))
    if (rectangular && identical(component, "contrast")) {
      .input_error("Signed contrast execution requires a self-form plan.")
    }
    if (storage != "memory" || !is.null(storage_path)) {
      .input_error(paste0(
        "Query-first execution returns an in-memory view and does not create ",
        "geometry stores."
      ))
    }
    task <- .compile_effect_evidence_task(
      plan$task$left_relation, plan$pairing,
      right = if (rectangular) plan$task$right_relation else NULL,
      query = query
    )
  } else {
    if (!is.null(signed_query)) {
      .input_error("Signed query weights require a fixed contrast query.")
    }
    physical_query <- NULL
    component <- "full"
    task <- plan$task
  }
  if (identical(component, "contrast")) {
    effects <- plan$task$left_relation$effect_space$coordinates
    signed_query <- .align_contrast(signed_query, effects)
    expected_query <- matrix(
      .svec_symmetric(tcrossprod(signed_query)), ncol = 1L
    )
    if (!identical(dim(physical_query), dim(expected_query)) ||
        !isTRUE(all.equal(
          unname(physical_query), expected_query, tolerance = 0
        ))) {
      .input_error(
        "A contrast plan's energy query must equal the signed outer product."
      )
    }
    requirements <- list(
      total = TRUE,
      coherent = TRUE,
      marginals = TRUE,
      materialization = "direct_contrast"
    )
  } else {
    if (!is.null(signed_query)) {
      .input_error(
        "Signed query weights are valid only for contrast execution."
      )
    }
    requirements <- .component_requirements(physical_query, component)
    if (rectangular) {
      # Signed endpoint marginals are a self-form concept; a rectangular
      # materialization carries total and coherent components only.
      requirements$marginals <- FALSE
    }
  }
  q <- plan$logical_shape[[1L]]
  output_width <- if (is.null(physical_query)) {
    plan$packed_width
  } else {
    .query_output_width(physical_query)
  }
  explicit_metric <- identical(
    plan$metric_schedule$kind, "fixed_metric_before_frame"
  )
  support_streamed <- identical(
    plan$metric_schedule$lowering, "support_streamed_pair_contraction"
  ) || (explicit_metric &&
    (isTRUE(requirements$coherent) || isTRUE(requirements$marginals)))
  if (support_streamed && storage == "block") {
    .input_error(paste0(
      "Support-streamed metric materialization currently requires in-memory ",
      "output; query-first execution remains bounded and is preferred."
    ))
  }
  selected <- if (support_streamed) {
    .support_metric_memory_plan(
      plan$task$left_relation, plan$frame, plan$metric_schedule,
      plan$compute, output_width, storage, requirements
    )
  } else {
    .select_compiler_memory_plan(
      plan$task$left_relation, plan$frame, plan$compute, output_width,
      storage, requirements, query = physical_query,
      right_relation = if (rectangular) plan$task$right_relation else NULL
    )
  }
  if (identical(selected$memory$fits_budget, FALSE)) {
    .input_error(sprintf(
      paste0("The conservative memory plan requires %.0f bytes, exceeding ",
        "the %.0f-byte budget."),
      selected$memory$planned_workspace_bytes,
      plan$compute$workspace_bytes
    ))
  }
  query_fused <- !is.null(physical_query)
  structured_query <- query_fused && .is_pair_difference_query(physical_query)
  lowering <- if (support_streamed && query_fused) {
    "support_streamed_metric_query_contraction"
  } else if (support_streamed) {
    "support_streamed_metric_form_contraction"
  } else if (explicit_metric && query_fused) {
    "additive_metric_query_fused_contraction"
  } else if (explicit_metric) {
    "additive_metric_form_contraction"
  } else if (query_fused) {
    "additive_query_fused_contraction"
  } else {
    "additive_form_contraction"
  }
  kernel_version <- .geometry_kernel_version(
    support_streamed, explicit_metric, query_fused, structured_query
  )
  scientific_plan_id <- if (is.null(physical_query)) {
    # Full materialization executes the plan's own estimand.
    plan$scientific_plan_id
  } else {
    # A fused query is an execution strategy for a view of the parent
    # estimand; its identity derives from the parent id plus the view
    # semantics, exactly as the projected route derives it.
    .geometry_view_scientific_id(
      plan$scientific_plan_id, component, physical_query, signed_query
    )
  }
  fields <- list(
    parent_signature = plan$signature,
    estimand_id = plan$scientific_plan_id,
    task = task,
    frame = plan$frame,
    metric_schedule = plan$metric_schedule,
    compute = plan$compute,
    storage = storage,
    storage_path = storage_path,
    query = physical_query,
    component = component,
    signed_query = signed_query,
    requirements = requirements,
    output_width = as.integer(output_width),
    feature_block = selected$feature_block,
    row_tile = selected$row_tile,
    coordinate_tile = selected$coordinate_tile,
    memory = selected$memory,
    task_count = if (support_streamed) {
      as.double(plan$measurements)
    } else {
      ceiling(plan$task$left_relation$n_features / selected$feature_block)
    },
    lowering = lowering,
    query_fused = query_fused,
    kernel_version = kernel_version,
    scientific_plan_id = scientific_plan_id
  )
  structure(c(fields, list(
    signature = .geometry_execution_signature(fields)
  )), class = "effect_geometry_execution_plan")
}

.validate_geometry_execution_plan <- function(x) {
  expected <- c("parent_signature", "estimand_id", "task", "frame",
    "metric_schedule",
    "compute", "storage", "storage_path", "query", "component",
    "signed_query", "requirements", "output_width", "feature_block", "row_tile",
    "coordinate_tile", "memory", "task_count", "lowering", "query_fused",
    "kernel_version", "scientific_plan_id", "signature")
  if (!.sealed_fields(x, "effect_geometry_execution_plan", expected) ||
      !.strong_sha256(x$parent_signature) ||
      !.strong_sha256(sub("^geometry-", "", x$estimand_id)) ||
      !x$storage %in% c("memory", "block") ||
      !x$component %in% c("full", "total", "coherent", "configuration", "contrast") ||
      !.is_flag(x$query_fused) || !.strong_sha256(x$signature) ||
      !.strong_sha256(sub("^geometry-", "", x$scientific_plan_id))) {
    .input_error("Geometry execution-plan fields are missing or noncanonical.")
  }
  .validate_evidence_task(x$task)
  .validate_frame_for_compile(x$frame)
  .validate_geometry_metric_schedule(x$metric_schedule)
  .validate_compute_policy(x$compute)
  .validate_memory_plan_for_receipt(x$memory)
  for (name in c("output_width", "feature_block", "row_tile",
      "coordinate_tile")) {
    value <- x[[name]]
    if (!is.integer(value) || length(value) != 1L || is.na(value) || value < 1L) {
      .input_error("Geometry execution-plan dimensions are invalid.")
    }
  }
  if (!.is_number(x$task_count) || x$task_count < 1 || x$task_count %% 1 != 0) {
    .input_error("Geometry execution-plan task count is invalid.")
  }
  expected_requirements <- if (identical(x$component, "contrast")) {
    list(
      total = TRUE,
      coherent = TRUE,
      marginals = TRUE,
      materialization = "direct_contrast"
    )
  } else {
    .component_requirements(x$query, x$component)
  }
  if (!isTRUE(x$task$same_relation)) {
    expected_requirements$marginals <- FALSE
  }
  if (identical(x$component, "contrast")) {
    expected_query <- matrix(
      .svec_symmetric(tcrossprod(x$signed_query)), ncol = 1L
    )
    if (!identical(dim(x$query), dim(expected_query)) ||
        !isTRUE(all.equal(unname(x$query), expected_query, tolerance = 0))) {
      .contract_error(
        "Contrast execution query and signed weights are inconsistent."
      )
    }
  }
  explicit_metric <- identical(
    x$metric_schedule$kind, "fixed_metric_before_frame"
  )
  support_streamed <- identical(
    x$metric_schedule$lowering, "support_streamed_pair_contraction"
  ) || (explicit_metric &&
    (isTRUE(x$requirements$coherent) || isTRUE(x$requirements$marginals)))
  expected_lowering <- if (support_streamed && x$query_fused) {
    "support_streamed_metric_query_contraction"
  } else if (support_streamed) {
    "support_streamed_metric_form_contraction"
  } else if (explicit_metric && x$query_fused) {
    "additive_metric_query_fused_contraction"
  } else if (explicit_metric) {
    "additive_metric_form_contraction"
  } else if (x$query_fused) {
    "additive_query_fused_contraction"
  } else {
    "additive_form_contraction"
  }
  expected_kernel <- .geometry_kernel_version(
    support_streamed, explicit_metric, x$query_fused,
    x$query_fused && .is_pair_difference_query(x$query)
  )
  expected_id <- if (is.null(x$query)) {
    x$estimand_id
  } else {
    .geometry_view_scientific_id(
      x$estimand_id, x$component, x$query, x$signed_query
    )
  }
  # The claimed estimand id is bound two ways: to the stored task through the
  # query-stripped base task id, and to the parent plan through the parent
  # signature, which digests exactly (estimand id, compute, payload bytes).
  expected_estimand_id <- .sha256_signature(list(
    schema_version = 1L,
    evidence_task = .effect_task_base_id(x$task),
    frame = .additive_frame_signature(x$frame),
    metric_schedule = x$metric_schedule$signature,
    component = "full",
    signed_query = NULL
  ), "geometry-sha256:")
  relation <- x$task$left_relation
  right_relation <- x$task$right_relation
  same <- isTRUE(x$task$same_relation)
  q <- length(relation$effect_space$coordinates)
  q_right <- length(right_relation$effect_space$coordinates)
  expected_parent_signature <- .geometry_plan_signature(
    x$estimand_id, x$compute,
    .geometry_dense_payload_bytes(
      nrow(x$frame$weights),
      if (same) q * (q + 1L) / 2L else q * q_right,
      if (same) q else q + q_right,
      length(unique(c(relation$partitions, if (same) NULL else
        right_relation$partitions)))
    )
  )
  fields <- x[names(x) != "signature"]
  if (!identical(x$requirements, expected_requirements) ||
      !identical(x$estimand_id, expected_estimand_id) ||
      !identical(x$parent_signature, expected_parent_signature) ||
      !identical(x$query_fused, !is.null(x$query)) ||
      !identical(x$lowering, expected_lowering) ||
      !identical(x$kernel_version, expected_kernel) ||
      !identical(x$scientific_plan_id, expected_id) ||
      !identical(x$signature, .geometry_execution_signature(fields))) {
    .contract_error(
      "Geometry execution-plan identity or lowering is inconsistent."
    )
  }
  invisible(x)
}
