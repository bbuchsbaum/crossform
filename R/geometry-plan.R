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
      signature = .sha256_signature(semantic)
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
    signature = .sha256_signature(semantic)
  )), class = "effect_metric_schedule")
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
  .sha256_signature(list(
    schema_version = 1L,
    evidence_task = task$task_id,
    frame = .additive_frame_signature(frame),
    metric_schedule = metric_schedule$signature,
    component = component,
    signed_query = signed_query
  ), "geometry-sha256:")
}

# A view's scientific identity is the parent geometry estimand plus the
# view's semantic descriptor. It must not depend on how the view was
# executed: the fused query-first route and projection from a materialized
# geometry are two executions of one estimand and must hash identically.
.geometry_view_scientific_id <- function(estimand_id, component, query,
                                         signed_query = NULL) {
  if (!is.character(estimand_id) || length(estimand_id) != 1L ||
      is.na(estimand_id) || !nzchar(estimand_id)) {
    stop("A view identity requires its parent geometry estimand id.",
      call. = FALSE)
  }
  .sha256_signature(list(
    schema_version = 1L,
    role = "geometry_view",
    parent = estimand_id,
    component = component,
    query = .query_identity_semantic(query),
    signed_query = if (is.null(signed_query)) NULL else unname(signed_query)
  ), "geometry-sha256:")
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
  .sha256_signature(list(
    schema_version = 1L,
    scientific_plan_id = scientific_plan_id,
    compute = unclass(compute),
    dense_payload_bytes = dense_payload_bytes
  ))
}

# Shared entry validation for the plan constructors: a relation, the frame
# it is measured at, and the pairing it may form must agree on the exact
# neural domain and on the feature and partition axes before anything is
# compiled.
.validate_geometry_plan_inputs <- function(x, at, over, right = NULL) {
  .validate_relation(x)
  if (!is.null(right)) .validate_relation(right)
  .validate_frame_for_compile(at)
  .validate_pairing(over)
  if (!identical(at$representation, "additive_diagonal")) {
    stop("crossform 0.1 executes only additive diagonal frames.", call. = FALSE)
  }
  mismatched <- if (!.same_domain_reference(at$domain, x$domain)) {
    x$domain
  } else if (!is.null(right) &&
      !.same_domain_reference(at$domain, right$domain)) {
    right$domain
  } else {
    NULL
  }
  if (!is.null(mismatched)) {
    side <- if (is.null(right) || identical(mismatched, x$domain)) {
      "relation"
    } else {
      "right relation"
    }
    stop(sprintf(paste0(
      "The %s and the frame were built on different neural domains, so a ",
      "feature index does not mean the same location on both sides. The %s ",
      "carries domain `%s` (%s, %s), the frame carries `%s` (%s, %s). ",
      "Compile the frame on the relation's own domain: ",
      "`compile_frame(<frame>, <relation>$domain)`."
    ),
      side, side,
      mismatched$id, .msg_count(mismatched$n_features, "feature"),
      .msg_signature(mismatched$signature),
      at$domain$id, .msg_count(at$domain$n_features, "feature"),
      .msg_signature(at$domain$signature)
    ), call. = FALSE)
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

#' Compile a reusable query-first geometry plan
#'
#' `plan_geometry()` validates the relation, spatial frame, pairing, source
#' capabilities, and compute policy without reading relation blocks. The plan
#' can answer fixed queries with [evaluate_geometry()] or be explicitly
#' materialized with [materialize_geometry()]. Complete packed geometry is
#' therefore an optional materialization, not the object that every analysis
#' must allocate.
#'
#' @param x An `effect_relation` supplying the left experimental axis.
#' @param at A compiled additive `effect_frame`.
#' @param over An `effect_pairing`. Rectangular plans require directed
#'   pairings whose left endpoints identify partitions of `x` and right
#'   endpoints identify partitions of `right`.
#' @param compute A sequential `compute_policy()`.
#' @param metric Optional fixed `neural_metric()`. A domain-wide metric is
#'   restricted to each frame support on demand; a support-local metric is
#'   accepted only for a matching one-node frame. Learned metric recipes are
#'   added by the statistical fitting layer rather than stored per node.
#'   Fixed metrics are not yet admitted on rectangular plans.
#' @param right Optional second `effect_relation` supplying a distinct right
#'   experimental axis. Supplying it compiles a rectangular cross-axis plan:
#'   the resulting form has one row axis per left effect and one column axis
#'   per right effect, is read with axis-bound [pair_query()]s through
#'   [evaluate_geometry()], and materializes to a rectangular effect form.
#'   Encoding-retrieval similarity is the canonical use.
#' @return An `effect_geometry_plan` recording the compiled `$task`,
#'   `$frame`, `$pairing`, `$metric_schedule`, and `$compute` policy, the
#'   `$logical_shape` and `$measurements` it will produce, and a
#'   `$scientific_plan_id` naming the estimand. Changing block size or storage
#'   changes the execution receipt, not this identity.
#' @section Structure:
#' A plan is a declaration, so every element describes what will be computed
#' rather than a result.
#'
#' - `$frame`: the compiled `effect_frame` the geometry is measured at.
#' - `$pairing`: the `effect_pairing` whose rows are the partition products
#'   the plan may form.
#' - `$measurements`: how many spatial measurements every view will return,
#'   one per row of `$frame$weights`.
#' - `$logical_shape`: the effect-axis extent as `c(left, right)`. The two
#'   entries are equal on a self-form plan.
#' - `$task$left_relation`: the relation supplying the left experimental
#'   axis. Its `$effects` is the order unnamed contrast weights are read in.
#' - `$compute`: the [compute_policy()] the plan was compiled under.
#' - `$scientific_plan_id`: the estimand identity. Keep it with the analysis
#'   record; block size, storage, and machine do not change it.
#'
#' The rest of `$task`, the `$metric_schedule`, and every other element not
#' listed here are the lowered form the executor consumes: internal, and
#' free to change.
#' @seealso [contrast_energy()], [rdm()], [rsa()], and [crossnobis()] for the
#'   named views of a plan; [evaluate_geometry()] for an arbitrary fixed
#'   query and [materialize_geometry()] for complete packed geometry;
#'   [coupling()] for the adjoint neural-side closure.
#' @family geometry plans and views
#' @examples
#' domain <- abstract_domain(3, id = "plan-example")
#' run1 <- rbind(a = c(1, 0, 2), b = c(0, 1, 1))
#' run2 <- rbind(a = c(1.1, 0.1, 1.9), b = c(0.1, 0.9, 1.2))
#' relation <- relation(list(run1 = run1, run2 = run2), domain = domain)
#'
#' # The plan validates relation, frame, pairing, and metric agreement
#' # without reading any neural value.
#' plan <- plan_geometry(
#'   relation,
#'   compile_frame(voxelwise(), domain),
#'   cross_partitions(relation, independence = "independent")
#' )
#' plan
#'
#' # One plan answers many fixed queries. Save its identity with the
#' # analysis record: it names the estimand, not the execution.
#' plan$scientific_plan_id
#' result <- evaluate_geometry(plan, query = bilinear_query(diag(2)))
#' as.data.frame(result)
#'
#' # The usual next step is a named view of the same plan.
#' contrast_energy(plan, c(a = 1, b = -1))$total
#' @export
plan_geometry <- function(x, at, over, compute = compute_policy(),
                          metric = NULL, right = NULL) {
  if (missing(at)) {
    stop(paste0(
      "`at` is required: pass a compiled frame from `compile_frame()`, for ",
      "example `compile_frame(searchlights(8), domain)`. It declares where ",
      "the geometry is measured."
    ), call. = FALSE)
  }
  if (missing(over)) {
    stop(paste0(
      "`over` is required: pass a pairing from `cross_partitions()` or ",
      "`pairing()`, for example ",
      "`cross_partitions(x, independence = \"independent\")`. It declares ",
      "which partition products the plan may form."
    ), call. = FALSE)
  }
  compute <- .validate_compute_policy(compute)
  .validate_geometry_plan_inputs(x, at, over, right = right)
  if (!is.null(right) && !is.null(metric)) {
    .capability_refusal(paste0(
      "Fixed neural metrics are not yet admitted on rectangular ",
      "cross-axis plans; crossform 0.1 compiles rectangular forms on ",
      "the implicit identity metric only."
    ),
      capability = "rectangular_fixed_metric",
      namespace = "geometry_plans",
      reasons = "rectangular_metric_law_not_compiled",
      remedies = paste0(
        "Omit `metric` for the rectangular plan, or use a self-form plan ",
        "for fixed-metric geometry."
      )
    )
  }
  capabilities <- .execution_preflight(compute, function() {
    .relation_source_capabilities(x)
  })$source_capabilities
  x$capabilities <- capabilities
  if (!is.null(right)) {
    right$capabilities <- .relation_source_capabilities(right)
  }
  task <- .compile_effect_evidence_task(x, over, right = right)
  metric_schedule <- .geometry_metric_schedule(at, metric)
  q_left <- length(x$effect_space$coordinates)
  q_right <- if (is.null(right)) {
    q_left
  } else {
    length(right$effect_space$coordinates)
  }
  codec <- if (is.null(right)) "symmetric_packed" else "rectangular"
  packed_width <- if (is.null(right)) {
    q_left * (q_left + 1L) / 2L
  } else {
    q_left * q_right
  }
  measurements <- nrow(at$weights)
  dense_payload_bytes <- .geometry_dense_payload_bytes(
    measurements, packed_width,
    if (is.null(right)) q_left else q_left + q_right,
    length(unique(c(x$partitions, if (is.null(right)) NULL else
      right$partitions)))
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
    codec = codec,
    logical_shape = as.integer(c(q_left, q_right)),
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
  if (!inherits(x, "effect_geometry_plan")) {
    stop(sprintf(paste0(
      "Expected an `effect_geometry_plan` from `plan_geometry()`; received %s."
    ), .msg_value(x)), call. = FALSE)
  }
  if (.validated_before(x, "geometry_plan", deep)) return(invisible(x))
  expected <- c("task", "frame", "pairing", "metric_schedule", "compute",
    "codec", "logical_shape", "packed_width", "measurements",
    "dense_payload_bytes", "lowering", "scientific_plan_id", "signature")
  valid_lowering <- c(
    "additive_contraction",
    "support_streamed_pair_contraction"
  )
  if (!inherits(x, "effect_geometry_plan") || !is.list(x) ||
      !identical(names(x), expected) ||
      !is.character(x$lowering) || length(x$lowering) != 1L ||
      is.na(x$lowering) || !x$lowering %in% valid_lowering ||
      !is.character(x$codec) || length(x$codec) != 1L ||
      !x$codec %in% c("symmetric_packed", "rectangular") ||
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
  metric_schedule <- .validate_geometry_metric_schedule(
    x$metric_schedule, deep = deep
  )
  compute <- .validate_compute_policy(x$compute)
  relation <- x$task$left_relation
  right_relation <- x$task$right_relation
  same <- isTRUE(x$task$same_relation)
  q_left <- length(relation$effect_space$coordinates)
  q_right <- length(right_relation$effect_space$coordinates)
  packed_width <- if (same) {
    q_left * (q_left + 1L) / 2L
  } else {
    q_left * q_right
  }
  measurements <- nrow(x$frame$weights)
  if (!identical(x$codec,
        if (same) "symmetric_packed" else "rectangular") ||
      (!same && !is.null(metric_schedule$metric)) ||
      x$task$materialization$kind != "effect_form" ||
      x$task$materialization$completeness != "complete_form" ||
      x$task$experimental_boundary$state != "open" ||
      x$task$neural_boundary$state != "closed" ||
      x$task$neural_boundary$closure_kind != "bridge" ||
      !identical(x$logical_shape, as.integer(c(q_left, q_right))) ||
      !identical(x$packed_width, as.integer(packed_width)) ||
      !identical(x$measurements, as.integer(measurements)) ||
      !.same_domain_reference(x$frame$domain, relation$domain) ||
      !.same_domain_reference(x$frame$domain, right_relation$domain) ||
      !identical(x$lowering, metric_schedule$lowering) ||
      !identical(x$dense_payload_bytes, .geometry_dense_payload_bytes(
        measurements, packed_width,
        if (same) q_left else q_left + q_right,
        length(unique(c(relation$partitions, if (same) NULL else
          right_relation$partitions)))
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
  .record_validated(x, "geometry_plan", deep)
  invisible(x)
}
