# Query-first geometry plans -----------------------------------------------

.geometry_metric_schedule <- function(frame, metric = NULL) {
  .validate_frame_for_compile(frame)
  if (!is.null(metric)) {
    metric <- .validate_neural_metric(metric)
    if (!.same_domain_reference(metric$domain, frame$domain)) {
      stop("The neural metric and spatial frame must share one exact domain.",
        call. = FALSE)
    }
    full_support <- identical(metric$support, frame$domain$feature_ids)
    if (!full_support) {
      if (nrow(frame$weights) != 1L) {
        stop(paste0(
          "A support-local fixed metric can be used only with its one ",
          "matching frame node. Multi-node plans require a domain-wide ",
          "operator or an on-demand metric recipe."
        ), call. = FALSE)
      }
      node <- .frame_metric_node(frame, 1L)
      if (!identical(metric$support, node$support)) {
        stop("The support-local metric does not match the frame node.",
          call. = FALSE)
      }
    }
    scope <- if (full_support) "domain_operator" else "single_node"
    lowering <- .metric_lowering(metric)
    semantic <- list(
      schema_version = 1L,
      role = "same_space_metric_schedule",
      kind = "fixed_metric_before_frame",
      frame_composition = "sqrt_weight_congruence",
      feature_additive = metric$capabilities$feature_additive,
      support_dense = metric$capabilities$support_dense,
      materialization = "fixed_metric",
      scope = scope,
      lowering = lowering,
      metric_signature = metric$signature
    )
    return(structure(c(semantic[-1L], list(
      metric = metric,
      signature = paste0("sha256:", digest::digest(
        semantic, algo = "sha256", serialize = TRUE
      ))
    )), class = "effect_metric_schedule"))
  }
  semantic <- list(
    schema_version = 1L,
    role = "same_space_metric_schedule",
    kind = "implicit_identity_before_frame",
    frame_composition = "sqrt_weight_congruence",
    feature_additive = TRUE,
    support_dense = FALSE,
    materialization = "implicit",
    scope = "domain_operator",
    lowering = "additive_contraction",
    metric_signature = NULL
  )
  structure(c(semantic[-1L], list(
    metric = NULL,
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE
    ))
  )), class = "effect_metric_schedule")
}

.validate_geometry_metric_schedule <- function(x) {
  expected <- c("role", "kind", "frame_composition", "feature_additive",
    "support_dense", "materialization", "scope", "lowering",
    "metric_signature", "metric", "signature")
  if (!inherits(x, "effect_metric_schedule") || !is.list(x) ||
      !identical(names(x), expected) ||
      !identical(x$role, "same_space_metric_schedule") ||
      !x$kind %in% c(
        "implicit_identity_before_frame", "fixed_metric_before_frame"
      ) ||
      !identical(x$frame_composition, "sqrt_weight_congruence") ||
      !is.logical(x$feature_additive) || length(x$feature_additive) != 1L ||
      is.na(x$feature_additive) ||
      !is.logical(x$support_dense) || length(x$support_dense) != 1L ||
      is.na(x$support_dense) || identical(x$feature_additive, x$support_dense) ||
      !x$materialization %in% c("implicit", "fixed_metric") ||
      !x$scope %in% c("domain_operator", "single_node") ||
      !x$lowering %in% c(
        "additive_contraction", "support_streamed_pair_contraction"
      ) ||
      !.strong_sha256(x$signature)) {
    stop("Geometry metric-schedule fields are missing or noncanonical.",
      call. = FALSE)
  }
  if (identical(x$kind, "implicit_identity_before_frame")) {
    if (!identical(x$feature_additive, TRUE) ||
        !identical(x$support_dense, FALSE) ||
        !identical(x$materialization, "implicit") ||
        !identical(x$lowering, "additive_contraction") ||
        !is.null(x$metric_signature) || !is.null(x$metric)) {
      stop("The implicit identity metric schedule is inconsistent.",
        call. = FALSE)
    }
  } else {
    metric <- .validate_neural_metric(x$metric)
    if (!identical(x$materialization, "fixed_metric") ||
        !identical(x$metric_signature, metric$signature) ||
        !identical(x$feature_additive,
          metric$capabilities$feature_additive) ||
        !identical(x$support_dense, metric$capabilities$support_dense) ||
        !identical(x$lowering, .metric_lowering(metric))) {
      stop("The fixed neural metric schedule is inconsistent.", call. = FALSE)
    }
  }
  semantic <- c(list(schema_version = 1L), unclass(x[
    !names(x) %in% c("metric", "signature")
  ]))
  expected_signature <- paste0("sha256:", digest::digest(
    semantic, algo = "sha256", serialize = TRUE
  ))
  if (!identical(x$signature, expected_signature)) {
    stop("Geometry metric-schedule identity is inconsistent.", call. = FALSE)
  }
  x
}

.geometry_plan_scientific_id <- function(task, frame, metric_schedule,
                                         component = "full",
                                         signed_query = NULL,
                                         validate = TRUE) {
  if (isTRUE(validate)) {
    .validate_evidence_task(task)
    .validate_frame_for_compile(frame)
    metric_schedule <- .validate_geometry_metric_schedule(metric_schedule)
  }
  if (!is.character(component) || length(component) != 1L ||
      is.na(component) ||
      !component %in% c(
        "full", "total", "coherent", "configuration", "contrast"
      )) {
    stop("Geometry-plan component identity is invalid.", call. = FALSE)
  }
  if (identical(component, "contrast")) {
    effects <- task$left_relation$effect_space$coordinates
    if (!is.numeric(signed_query) || length(signed_query) != length(effects) ||
        any(!is.finite(signed_query)) ||
        !identical(names(signed_query), effects)) {
      stop("Contrast plans require one named signed weight per effect.",
        call. = FALSE)
    }
  } else if (!is.null(signed_query)) {
    stop("Signed query weights are valid only for contrast plans.",
      call. = FALSE)
  }
  paste0("geometry-sha256:", digest::digest(list(
    schema_version = 1L,
    evidence_task = task$task_id,
    frame = .additive_frame_signature(frame),
    metric_schedule = metric_schedule$signature,
    component = component,
    signed_query = signed_query
  ), algo = "sha256", serialize = TRUE))
}

.geometry_dense_payload_bytes <- function(measurements, packed_width,
                                          effects, partitions) {
  values <- c(measurements, packed_width, effects, partitions)
  if (any(!is.finite(values)) || any(values < 1)) {
    stop("Geometry logical dimensions must be positive and finite.",
      call. = FALSE)
  }
  # Two packed components plus one endpoint-equivalent marginal family.  The
  # number is a payload estimate, deliberately separate from R object overhead.
  bytes <- 8 * (2 * measurements * packed_width +
    measurements * effects * partitions)
  if (!is.finite(bytes) || bytes > 2^53) {
    stop("Dense geometry payload exceeds exact byte accounting.", call. = FALSE)
  }
  as.double(bytes)
}

.geometry_plan_signature <- function(scientific_plan_id, compute,
                                     dense_payload_bytes) {
  paste0("sha256:", digest::digest(list(
    schema_version = 1L,
    scientific_plan_id = scientific_plan_id,
    compute = unclass(compute),
    dense_payload_bytes = dense_payload_bytes
  ), algo = "sha256", serialize = TRUE))
}

#' Compile a reusable query-first geometry plan
#'
#' `plan_geometry()` validates the relation, spatial frame, pairing, source
#' capabilities, and compute policy without reading relation blocks. The plan
#' can answer fixed queries with [evaluate_geometry()] or be explicitly
#' materialized with [geometry()]. Complete packed geometry is therefore an
#' optional materialization, not the object that every analysis must allocate.
#'
#' @param x An `effect_relation`.
#' @param at A compiled additive `effect_frame`.
#' @param over An `effect_pairing`.
#' @param compute A sequential `compute_policy()`.
#' @param metric Optional fixed `neural_metric()`. A domain-wide metric is
#'   restricted to each frame support on demand; a support-local metric is
#'   accepted only for a matching one-node frame. Learned metric recipes are
#'   added by the statistical fitting layer rather than stored per node.
#' @return An immutable-by-convention `effect_geometry_plan`.
#' @examples
#' domain <- abstract_domain(3, id = "plan-example")
#' run1 <- rbind(a = c(1, 0, 2), b = c(0, 1, 1))
#' run2 <- rbind(a = c(1.1, 0.1, 1.9), b = c(0.1, 0.9, 1.2))
#' relation <- relation(list(run1 = run1, run2 = run2), domain = domain)
#' plan <- plan_geometry(
#'   relation,
#'   compile_frame(voxels(), domain),
#'   cross_partitions(relation)
#' )
#' result <- evaluate_geometry(plan, query = bilinear_query(diag(2)))
#' result
#' as.data.frame(result)
#' @export
plan_geometry <- function(x, at, over, compute = compute_policy(),
                          metric = NULL) {
  compute <- .validate_compute_policy(compute)
  .validate_compiler_inputs(x, at, over)
  capabilities <- .execution_preflight(compute, function() {
    .compiler_capabilities(x)
  })$source_capabilities
  x$capabilities <- capabilities
  task <- .compile_effect_evidence_task(x, over)
  metric_schedule <- .geometry_metric_schedule(at, metric)
  q <- length(x$effect_space$coordinates)
  packed_width <- q * (q + 1L) / 2L
  measurements <- nrow(at$weights)
  dense_payload_bytes <- .geometry_dense_payload_bytes(
    measurements, packed_width, q, length(x$partitions)
  )
  scientific_plan_id <- .geometry_plan_scientific_id(
    task, at, metric_schedule, "full", validate = FALSE
  )
  signature <- .geometry_plan_signature(
    scientific_plan_id, compute, dense_payload_bytes
  )
  value <- structure(list(
    task = task,
    frame = at,
    pairing = over,
    metric_schedule = metric_schedule,
    compute = compute,
    logical_shape = as.integer(c(q, q)),
    packed_width = as.integer(packed_width),
    measurements = as.integer(measurements),
    dense_payload_bytes = dense_payload_bytes,
    lowering = metric_schedule$lowering,
    scientific_plan_id = scientific_plan_id,
    signature = signature
  ), class = "effect_geometry_plan")
  .validate_geometry_plan(value, deep = FALSE)
  value
}

.validate_geometry_plan <- function(x, deep = TRUE) {
  expected <- c("task", "frame", "pairing", "metric_schedule", "compute",
    "logical_shape", "packed_width", "measurements",
    "dense_payload_bytes", "lowering", "scientific_plan_id", "signature")
  valid_lowering <- c(
    "additive_contraction",
    "support_streamed_pair_contraction"
  )
  if (!inherits(x, "effect_geometry_plan") || !is.list(x) ||
      !identical(names(x), expected) ||
      !is.character(x$lowering) || length(x$lowering) != 1L ||
      is.na(x$lowering) || !x$lowering %in% valid_lowering ||
      !is.integer(x$logical_shape) || length(x$logical_shape) != 2L ||
      !is.integer(x$packed_width) || length(x$packed_width) != 1L ||
      !is.integer(x$measurements) || length(x$measurements) != 1L ||
      !is.numeric(x$dense_payload_bytes) ||
      length(x$dense_payload_bytes) != 1L ||
      !is.finite(x$dense_payload_bytes) || x$dense_payload_bytes < 0 ||
      !.strong_sha256(sub("^geometry-", "", x$scientific_plan_id)) ||
      !.strong_sha256(x$signature)) {
    stop("Geometry-plan fields are missing or noncanonical.", call. = FALSE)
  }
  .validate_evidence_task(x$task)
  .validate_frame_for_compile(x$frame)
  .validate_pairing(x$pairing)
  metric_schedule <- .validate_geometry_metric_schedule(x$metric_schedule)
  compute <- .validate_compute_policy(x$compute)
  relation <- x$task$left_relation
  q <- length(relation$effect_space$coordinates)
  packed_width <- q * (q + 1L) / 2L
  measurements <- nrow(x$frame$weights)
  if (!isTRUE(x$task$same_relation) ||
      x$task$materialization$kind != "effect_form" ||
      x$task$materialization$completeness != "complete_form" ||
      x$task$experimental_boundary$state != "open" ||
      x$task$neural_boundary$state != "closed" ||
      x$task$neural_boundary$closure_kind != "bridge" ||
      !identical(x$logical_shape, as.integer(c(q, q))) ||
      !identical(x$packed_width, as.integer(packed_width)) ||
      !identical(x$measurements, as.integer(measurements)) ||
      !.same_domain_reference(x$frame$domain, relation$domain) ||
      !identical(x$lowering, metric_schedule$lowering) ||
      !identical(x$dense_payload_bytes, .geometry_dense_payload_bytes(
        measurements, packed_width, q, length(relation$partitions)
      ))) {
    stop("Geometry-plan axes, frame, or materialization are inconsistent.",
      call. = FALSE)
  }
  if (!is.null(metric_schedule$metric)) {
    metric <- metric_schedule$metric
    if (!.same_domain_reference(metric$domain, relation$domain) ||
        (metric_schedule$scope == "domain_operator" &&
         !identical(metric$support, relation$domain$feature_ids)) ||
        (metric_schedule$scope == "single_node" && measurements != 1L)) {
      stop("Geometry-plan metric support is inconsistent with its frame.",
        call. = FALSE)
    }
  }
  if (isTRUE(deep)) {
    expected_id <- .geometry_plan_scientific_id(
      x$task, x$frame, x$metric_schedule, "full", validate = FALSE
    )
    expected_signature <- .geometry_plan_signature(
      expected_id, compute, x$dense_payload_bytes
    )
    if (!identical(x$scientific_plan_id, expected_id) ||
        !identical(x$signature, expected_signature)) {
      stop("Geometry-plan identity is inconsistent.", call. = FALSE)
    }
  }
  invisible(x)
}

#' @export
print.effect_geometry_plan <- function(x, ...) {
  .validate_geometry_plan(x, deep = FALSE)
  payload <- structure(x$dense_payload_bytes, class = "object_size")
  cat("<effect_geometry_plan>\n", sep = "")
  cat("  effects:      ", x$logical_shape[[1L]], " x ",
    x$logical_shape[[2L]], "\n", sep = "")
  cat("  measurements: ", x$measurements, "\n", sep = "")
  cat("  features:     ", x$task$left_relation$n_features, "\n", sep = "")
  cat("  lowering:     ", x$lowering, "\n", sep = "")
  cat("  dense payload:", format(payload, units = "auto"), "\n")
  cat("  state:        query-first; call geometry(plan) to materialize\n")
  invisible(x)
}

.geometry_execution_signature <- function(fields) {
  paste0("sha256:", digest::digest(list(
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
  ), algo = "sha256", serialize = TRUE))
}

.compile_geometry_execution_plan <- function(
    plan, query = NULL, component = NULL,
    storage = c("memory", "block"), storage_path = NULL,
    signed_query = NULL) {
  .validate_geometry_plan(plan)
  storage <- match.arg(storage)
  if (!is.null(query)) {
    physical_query <- .compiler_query(
      query, plan$task$left_relation$effect_space
    )
    component <- match.arg(component,
      c("total", "coherent", "configuration", "contrast"))
    if (storage != "memory" || !is.null(storage_path)) {
      stop(paste0(
        "Query-first execution returns an in-memory view and does not create ",
        "geometry stores."
      ), call. = FALSE)
    }
    task <- .compile_effect_evidence_task(
      plan$task$left_relation, plan$pairing, query = query
    )
  } else {
    if (!is.null(signed_query)) {
      stop("Signed query weights require a fixed contrast query.",
        call. = FALSE)
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
      stop("A contrast plan's energy query must equal the signed outer product.",
        call. = FALSE)
    }
    requirements <- list(
      total = TRUE,
      coherent = TRUE,
      marginals = TRUE,
      materialization = "direct_contrast"
    )
  } else {
    if (!is.null(signed_query)) {
      stop("Signed query weights are valid only for contrast execution.",
        call. = FALSE)
    }
    requirements <- .component_requirements(physical_query, component)
  }
  q <- plan$logical_shape[[1L]]
  output_width <- if (is.null(physical_query)) {
    plan$packed_width
  } else {
    ncol(physical_query)
  }
  explicit_metric <- identical(
    plan$metric_schedule$kind, "fixed_metric_before_frame"
  )
  support_streamed <- identical(
    plan$metric_schedule$lowering, "support_streamed_pair_contraction"
  ) || (explicit_metric &&
    (isTRUE(requirements$coherent) || isTRUE(requirements$marginals)))
  if (support_streamed && storage == "block") {
    stop(paste0(
      "Support-streamed metric materialization currently requires in-memory ",
      "output; query-first execution remains bounded and is preferred."
    ), call. = FALSE)
  }
  selected <- if (support_streamed) {
    .support_metric_memory_plan(
      plan$task$left_relation, plan$frame, plan$metric_schedule,
      plan$compute, output_width, storage, requirements
    )
  } else {
    .select_compiler_memory_plan(
      plan$task$left_relation, plan$frame, plan$compute, output_width,
      storage, requirements
    )
  }
  if (identical(selected$memory$fits_budget, FALSE)) {
    stop(sprintf(
      paste0("The conservative memory plan requires %.0f bytes, exceeding ",
        "the %.0f-byte budget."),
      selected$memory$planned_workspace_bytes,
      plan$compute$workspace_bytes
    ), call. = FALSE)
  }
  query_fused <- !is.null(physical_query)
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
  kernel_version <- if (support_streamed) {
    "support-streamed-metric-v1"
  } else if (explicit_metric) {
    "additive-fixed-metric-v1"
  } else if (query_fused) {
    "additive-query-fused-v2"
  } else {
    "additive-effect-form-v2"
  }
  scientific_plan_id <- .geometry_plan_scientific_id(
    task, plan$frame, plan$metric_schedule, component, signed_query,
    validate = FALSE
  )
  fields <- list(
    parent_signature = plan$signature,
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
  expected <- c("parent_signature", "task", "frame", "metric_schedule",
    "compute", "storage", "storage_path", "query", "component",
    "signed_query", "requirements", "output_width", "feature_block", "row_tile",
    "coordinate_tile", "memory", "task_count", "lowering", "query_fused",
    "kernel_version", "scientific_plan_id", "signature")
  if (!inherits(x, "effect_geometry_execution_plan") || !is.list(x) ||
      !identical(names(x), expected) || !.strong_sha256(x$parent_signature) ||
      !x$storage %in% c("memory", "block") ||
      !x$component %in% c(
        "full", "total", "coherent", "configuration", "contrast"
      ) ||
      !is.logical(x$query_fused) || length(x$query_fused) != 1L ||
      is.na(x$query_fused) || !.strong_sha256(x$signature) ||
      !.strong_sha256(sub("^geometry-", "", x$scientific_plan_id))) {
    stop("Geometry execution-plan fields are missing or noncanonical.",
      call. = FALSE)
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
      stop("Geometry execution-plan dimensions are invalid.", call. = FALSE)
    }
  }
  if (!is.numeric(x$task_count) || length(x$task_count) != 1L ||
      is.na(x$task_count) || !is.finite(x$task_count) || x$task_count < 1 ||
      x$task_count %% 1 != 0) {
    stop("Geometry execution-plan task count is invalid.", call. = FALSE)
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
  if (identical(x$component, "contrast")) {
    expected_query <- matrix(
      .svec_symmetric(tcrossprod(x$signed_query)), ncol = 1L
    )
    if (!identical(dim(x$query), dim(expected_query)) ||
        !isTRUE(all.equal(unname(x$query), expected_query, tolerance = 0))) {
      stop("Contrast execution query and signed weights are inconsistent.",
        call. = FALSE)
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
  expected_kernel <- if (support_streamed) {
    "support-streamed-metric-v1"
  } else if (explicit_metric) {
    "additive-fixed-metric-v1"
  } else if (x$query_fused) {
    "additive-query-fused-v2"
  } else {
    "additive-effect-form-v2"
  }
  expected_id <- .geometry_plan_scientific_id(
    x$task, x$frame, x$metric_schedule, x$component, x$signed_query,
    validate = FALSE
  )
  fields <- x[names(x) != "signature"]
  if (!identical(x$requirements, expected_requirements) ||
      !identical(x$query_fused, !is.null(x$query)) ||
      !identical(x$lowering, expected_lowering) ||
      !identical(x$kernel_version, expected_kernel) ||
      !identical(x$scientific_plan_id, expected_id) ||
      !identical(x$signature, .geometry_execution_signature(fields))) {
    stop("Geometry execution-plan identity or lowering is inconsistent.",
      call. = FALSE)
  }
  invisible(x)
}
