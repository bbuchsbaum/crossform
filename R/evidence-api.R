# Narrow public API for evidence-pairing materializations ------------------

.public_measurement_dense_limit_bytes <- 256 * 1024^2

.require_small_node_measurement_regime <- function(bytes, operation) {
  if (!is.numeric(bytes) || length(bytes) != 1L || is.na(bytes) ||
      !is.finite(bytes) || bytes < 0) {
    stop("Measurement resource preflight produced an invalid byte estimate.",
      call. = FALSE)
  }
  limit <- .public_measurement_dense_limit_bytes
  if (bytes > limit) {
    stop(sprintf(
      paste0(
        "%s requires approximately %.0f dense payload bytes, exceeding ",
        "the version 0.1 small-node limit of %.0f bytes. Use the ",
        "support-local geometry plan; brain-scale measurement tomography ",
        "is not yet exported."
      ),
      operation, bytes, limit
    ), call. = FALSE)
  }
  invisible(as.double(bytes))
}

.preflight_additive_measurement_frame <- function(frame, mode) {
  .validate_frame_for_compile(frame)
  p <- ncol(frame$weights)
  m <- nrow(frame$weights)
  nonzero <- if (inherits(frame$weights, "Matrix")) {
    Matrix::nnzero(frame$weights)
  } else {
    sum(frame$weights != 0)
  }
  leg_elements <- if (identical(mode, "coherent")) m * p else p * nonzero
  dense_payload <- 8 * (p * p + m * p + leg_elements)
  .require_small_node_measurement_regime(
    dense_payload, "The requested measurement frame"
  )
}

.preflight_custom_measurement_frame <- function(operators, domain) {
  .validate_domain(domain)
  if (!is.list(operators) || length(operators) < 1L) return(invisible(NULL))
  p <- domain$n_features
  leg_elements <- sum(vapply(operators, function(operator) {
    dimensions <- dim(operator)
    if (is.null(dimensions) || length(dimensions) != 2L ||
        any(!is.finite(dimensions)) || any(dimensions < 1)) {
      return(0)
    }
    prod(as.double(dimensions))
  }, numeric(1)))
  .require_small_node_measurement_regime(
    8 * (p * p + leg_elements), "The requested custom measurement frame"
  )
}

.public_measurement_frame_from_operators <- function(operators, domain, id,
                                                     units) {
  .validate_domain(domain)
  if (!is.list(operators) || length(operators) < 1L ||
      is.null(names(operators)) || anyNA(names(operators)) ||
      any(!nzchar(names(operators))) || anyDuplicated(names(operators))) {
    stop("Custom measurements require uniquely named operator matrices.",
      call. = FALSE)
  }
  .measurement_identifier(id, "Measurement-frame `id`")
  legs <- Map(function(operator, node) {
    if (inherits(operator, "Matrix")) operator <- as.matrix(operator)
    if (!is.matrix(operator) || !is.numeric(operator) ||
        nrow(operator) < 1L || ncol(operator) != domain$n_features ||
        any(!is.finite(operator))) {
      stop(paste0(
        "Every custom measurement must be a finite matrix with one column ",
        "per neural feature."
      ), call. = FALSE)
    }
    coordinates <- rownames(operator)
    if (is.null(coordinates) || length(coordinates) != nrow(operator) ||
        anyNA(coordinates) || any(!nzchar(coordinates)) ||
        anyDuplicated(coordinates)) {
      coordinates <- paste0("mode", seq_len(nrow(operator)))
    }
    axis_units <- if (length(units) == 1L) {
      rep(units, nrow(operator))
    } else {
      if (length(units) != nrow(operator)) {
        stop("Measurement units must be scalar or match every operator row.",
          call. = FALSE)
      }
      units
    }
    .measurement_leg(
      unname(operator), domain,
      .measurement_axis(
        coordinates,
        id = paste0(id, ":", node),
        basis_id = paste0(id, ":", node, ":fixed-basis"),
        units = axis_units,
        provenance = list(construction = "public_fixed_operator")
      ),
      estimation = "fixed",
      provenance = list(construction = "public_fixed_operator", node = node)
    )
  }, operators, names(operators))
  names(legs) <- names(operators)
  .measurement_frame(legs)
}

#' Declare identified neural measurements
#'
#' `measurement_frame()` preserves the oriented measurement legs needed for
#' cross-region forms. An additive frame can be adapted as total, coherent, or
#' coherent/configuration measurements. A named list of fixed matrices can be
#' used for prespecified measurements. Learned operators are intentionally not
#' accepted by this convenience constructor because they require frozen
#' training provenance.
#'
#' @param x A compiled additive `effect_frame`, or a uniquely named list of
#'   fixed measurement matrices.
#' @param domain Required for a list of matrices; the exact source
#'   `effect_domain`.
#' @param mode For an additive frame, one of `"total"`, `"coherent"`, or
#'   `"coherent_configuration"`.
#' @param id Stable identity used for custom measurement axes.
#' @param units One unit or one per row of each custom operator.
#' @details Version 0.1 exposes this as a deliberately small-node dense path.
#'   Construction fails before dense conversion when its estimated frame and
#'   leg payload exceeds 256 MiB. Brain-scale work should use the support-local
#'   geometry-plan path until matrix-free measurement frames are qualified.
#' @return An identified `effect_measurement_frame`.
#' @export
measurement_frame <- function(
    x, domain = NULL,
    mode = c("total", "coherent", "coherent_configuration"),
    id = "fixed-measurements:v1", units = "arbitrary") {
  mode <- match.arg(mode)
  if (inherits(x, "effect_frame")) {
    if (!is.null(domain)) {
      stop("An additive frame already identifies its neural domain.",
        call. = FALSE)
    }
    .preflight_additive_measurement_frame(x, mode)
    return(if (mode == "coherent_configuration") {
      .measurement_frame_from_additive_decomposition(x)
    } else {
      .measurement_frame_from_additive(x, mode)
    })
  }
  if (mode != "total") {
    stop(paste0(
      "Custom operator lists are already explicit measurements; `mode` must ",
      "remain `\"total\"`."
    ), call. = FALSE)
  }
  if (is.null(domain)) {
    stop("Custom measurement operators require their exact neural `domain`.",
      call. = FALSE)
  }
  .preflight_custom_measurement_frame(x, domain)
  .public_measurement_frame_from_operators(x, domain, id, units)
}

.edge_frame_signature <- function(fields) {
  paste0("sha256:", digest::digest(list(
    schema_version = 1L,
    from_frame = fields$from_frame$signature,
    to_frame = fields$to_frame$signature,
    edges = fields$edges$signature
  ), algo = "sha256", serialize = TRUE, serializeVersion = 2L))
}

.validate_edge_frame <- function(x) {
  expected <- c("from_frame", "to_frame", "edges", "signature")
  if (!inherits(x, "effect_edge_frame") || !is.list(x) ||
      !identical(names(x), expected)) {
    stop("`between` must be a canonical `edge_frame()`.", call. = FALSE)
  }
  from <- .validate_measurement_frame(x$from_frame)
  to <- .validate_measurement_frame(x$to_frame)
  .validate_measurement_edges(x$edges, from, to)
  fields <- x[setdiff(names(x), "signature")]
  if (!identical(x$signature, .edge_frame_signature(fields))) {
    stop("Edge-frame identity is inconsistent.", call. = FALSE)
  }
  x
}

#' Declare requested neural measurement edges
#'
#' Every edge is explicit. The constructor never creates all pairs on the
#' user's behalf. `from` and `to` name neural measurement nodes; the words
#' `left` and `right` are reserved for experimental relation sides in
#' [measurement_form()].
#'
#' @param from,to Equal-length vectors of node identifiers.
#' @param frame Measurement frame containing `from` nodes.
#' @param to_frame Measurement frame containing `to` nodes. Defaults to
#'   `frame`.
#' @param weight Optional finite edge weights recorded as part of edge identity.
#' @return An explicit `effect_edge_frame` for the `between` argument.
#' @export
edge_frame <- function(from, to, frame, to_frame = frame, weight = NULL) {
  frame <- .validate_measurement_frame(frame)
  to_frame <- .validate_measurement_frame(to_frame)
  fields <- list(
    from_frame = frame,
    to_frame = to_frame,
    edges = .measurement_edges(from, to, frame, to_frame, weight)
  )
  structure(c(fields, list(
    signature = .edge_frame_signature(fields)
  )), class = "effect_edge_frame")
}

#' Declare a repeated-variation experimental query
#'
#' A variation query binds a positive semidefinite operator to one identified
#' experimental space and records the sample axis that gives coupling its
#' repeated variation. A rank-one operator remains a valid effect query, but
#' normalized connectivity views reject it as degenerate.
#'
#' @param operator A finite symmetric positive-semidefinite matrix.
#' @param effects The exact `effect_space` on both query axes.
#' @param sampling_axis Identifier such as `"time"`, `"trial"`, or
#'   `"subject"`.
#' @param construction `"psd_variation"` for a positive variation form, or
#'   `"joint_covariance"` when the later partition construction is a coherent
#'   joint covariance.
#' @param provenance Named fixed-construction metadata.
#' @return An axis-bound `effect_pair_query` with variation capabilities.
#' @export
variation_query <- function(
    operator, effects, sampling_axis,
    construction = c("psd_variation", "joint_covariance"),
    provenance = list()) {
  .variation_pair_query(
    operator, effects, sampling_axis, match.arg(construction), provenance
  )
}

.public_query_capability <- function(by) {
  .validate_query_for_compile(by)
  capability <- by$metadata$evidence_capability
  if (is.null(capability)) {
    return(list(
      role = "effect", sampling_axis = NULL, construction = "arbitrary"
    ))
  }
  if (!is.list(capability) || !identical(capability$role, "variation") ||
      is.null(capability$sampling_axis) ||
      !capability$construction %in% c("psd_variation", "joint_covariance")) {
    stop("Experimental query capability metadata is invalid.", call. = FALSE)
  }
  .validate_variation_pair_query(
    by, capability$sampling_axis, capability$construction
  )
  list(
    role = "variation",
    sampling_axis = capability$sampling_axis,
    construction = capability$construction
  )
}

#' Materialize requested measurement-space forms
#'
#' `measurement_form()` closes the experimental boundary with `by` and leaves
#' only the explicitly requested neural edges in `between` open. `left` and
#' `right` always name experimental-neural relations; spatial endpoints are
#' contained inside `between`.
#'
#' Raw forms are bilinear. Covariance, correlation, CCA, geometry alignment,
#' and Gaussian information are separate capability-gated views. The default
#' reducer records that raw partition forms are aggregated before any later
#' nonlinear connectivity normalization.
#'
#' @param left The left `effect_relation`.
#' @param right The right `effect_relation`; defaults to `left`.
#' @param between An explicit `edge_frame()`.
#' @param by An axis-bound `pair_query()` or `variation_query()`.
#' @param over An explicit partition `pairing()`.
#' @param reducer A partition-order declaration. Raw forms are bilinear, but
#'   the choice remains part of scientific plan identity.
#' @param storage Either in-memory or block-backed storage.
#' @param storage_path Required for block-backed storage.
#' @param compute A `compute_policy()`.
#' @param route Internal contraction route. `"auto"` chooses an equivalent
#'   bounded route without changing the scientific plan.
#' @return A complete `effect_measurement_form` over the requested edge set.
#' @export
measurement_form <- function(
    left, between, by, over, right = left,
    reducer = aggregate_first(),
    storage = c("memory", "block"), storage_path = NULL,
    compute = compute_policy(),
    route = c("auto", "forward_k", "pull_h", "factorized_h",
              "scalar_stack", "multivariate_blocks")) {
  .validate_relation(left)
  .validate_relation(right)
  between <- .validate_edge_frame(between)
  .validate_pairing(over)
  reducer <- .validate_partition_reducer(reducer)
  storage <- match.arg(storage)
  route <- match.arg(route)
  same_relation <- identical(left, right)
  if (!.same_domain_reference(
      between$from_frame$source_domain, left$domain) ||
      !.same_domain_reference(
        between$to_frame$source_domain, right$domain)) {
    stop("Measurement frames do not belong to their experimental relation sides.",
      call. = FALSE)
  }
  capability <- .public_query_capability(by)
  if (!.same_effect_space(by$left_space, left$effect_space) ||
      !.same_effect_space(by$right_space, right$effect_space)) {
    stop("`by` does not match the left and right experimental spaces.",
      call. = FALSE)
  }
  partition_edges <- .ordered_partition_edges(
    over, left$partitions, right$partitions, same_relation
  )
  task <- .new_evidence_task(
    left, right, same_relation, partition_edges,
    .closed_experimental_boundary(
      by, capability$role, capability$sampling_axis
    ),
    .open_neural_boundary(
      between$from_frame, between$to_frame, between$edges
    ),
    .evidence_stage_plan(reducer = reducer),
    .evidence_materialization("measurement_form", "complete_form")
  )
  contraction <- .run_measurement_contraction(
    task, compute = compute, route = route, storage = storage
  )
  edges <- between$edges$edges
  edge_scope <- if (.measurement_edges_are_complete(
      edges, between$from_frame$node_ids, between$to_frame$node_ids)) {
    "frame_complete"
  } else {
    "requested"
  }
  .measurement_form_from_contraction(
    task, contraction,
    storage = storage, path = storage_path,
    query_construction = capability$construction,
    edge_scope = edge_scope
  )
}

#' Interpret completed measurement forms
#'
#' These functions are views over one completed measurement form, not separate
#' fitting engines. `effect_coupling()` makes no covariance claim.
#' `covariance_coupling()`, `canonical_coupling()`, and `connectivity()` require
#' certified repeated variation; normalized views additionally require
#' coherent positive self-covariances. `geometry_alignment()` is static linear
#' CKA/RV-like alignment and does not redefine established dynamic
#' informational-connectivity analyses.
#'
#' @param x A complete `effect_measurement_form`.
#' @param tolerance Positive numerical tolerance.
#' @return A typed `effect_coupling_result`.
#' @name coupling_views
NULL

#' @rdname coupling_views
#' @export
effect_coupling <- function(x) .effect_coupling(x)

#' @rdname coupling_views
#' @export
covariance_coupling <- function(x, tolerance = 1e-10) {
  .covariance_coupling(x, tolerance)
}

#' @rdname coupling_views
#' @param ridge,ridge_right Positive ridge values for the left and right
#'   covariance blocks.
#' @export
canonical_coupling <- function(x, ridge, ridge_right = ridge,
                               tolerance = 1e-10) {
  .canonical_coupling(
    x, .measurement_regularization("ridge", ridge, ridge_right), tolerance
  )
}

#' @rdname coupling_views
#' @export
geometry_alignment <- function(x, tolerance = 1e-10) {
  .geometry_alignment(x, tolerance)
}

#' Declare a joint Gaussian covariance interpretation
#'
#' This explicit declaration is required before canonical correlations are
#' transformed into Gaussian mutual information.
#'
#' @param provenance Named metadata describing the fixed model assumption.
#' @return An `effect_gaussian_covariance_model` declaration.
#' @export
gaussian_covariance_model <- function(provenance = list()) {
  .gaussian_covariance_model(provenance)
}

#' Request a validated connectivity view
#'
#' @param x A complete `effect_measurement_form`.
#' @param view One of signed scalar correlation, a canonical spectrum, static
#'   geometry alignment, or Gaussian mutual information.
#' @param ridge,ridge_right Explicit ridge values for canonical or Gaussian
#'   views.
#' @param model A [gaussian_covariance_model()] for Gaussian information.
#' @param units Information units, when applicable.
#' @param tolerance Positive numerical tolerance.
#' @return A typed `effect_coupling_result`.
#' @export
connectivity <- function(
    x, view = c("correlation", "canonical", "geometry_alignment",
                "gaussian_information"),
    ridge = NULL, ridge_right = ridge, model = NULL,
    units = c("nats", "bits"), tolerance = 1e-10) {
  view <- match.arg(view)
  regularization <- if (is.null(ridge)) NULL else
    .measurement_regularization("ridge", ridge, ridge_right)
  .connectivity_view(
    x, view, regularization, model, match.arg(units), tolerance
  )
}

#' Summarize a crossed node decomposition on one measurement edge
#'
#' The returned strengths are invariant to rotations of subspace-only
#' components. Raw configuration entries are not exposed as basis-free
#' scientific quantities.
#'
#' @param x A complete `effect_measurement_form` whose node measurements carry
#'   decompositions.
#' @param edge An edge number or edge identifier.
#' @return A data frame of crossed component dimensions, orientation status,
#'   Frobenius strengths, and strongest singular values.
#' @export
measurement_components <- function(x, edge) {
  lifted <- .lift_measurement_decomposition(x, edge)
  rows <- lapply(seq_len(nrow(lifted$index)), function(position) {
    value <- lifted$blocks[[position]]
    singular <- if (length(value)) {
      svd(value, nu = 0L, nv = 0L)$d
    } else {
      numeric()
    }
    data.frame(
      edge_id = lifted$edge_id,
      left_component = lifted$index$left_component[[position]],
      right_component = lifted$index$right_component[[position]],
      d_left = nrow(value),
      d_right = ncol(value),
      left_orientation = lifted$index$left_orientation[[position]],
      right_orientation = lifted$index$right_orientation[[position]],
      raw_entries_meaningful =
        lifted$index$left_orientation[[position]] == "oriented" &&
        lifted$index$right_orientation[[position]] == "oriented",
      frobenius_strength = sqrt(sum(value^2)),
      strongest_singular_value = if (length(singular)) max(singular) else 0,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  attr(result, "source_form") <- x$plan$signature
  attr(result, "decomposition") <- lifted$signature
  result
}

#' Reconstruct or project the global neural evidence operator
#'
#' Lossless reconstruction requires the complete node-edge block form and
#' full-column-rank identified frames. Rank-deficient frames return only an
#' explicitly projected operator. A workspace budget is checked before any
#' measurement block is read.
#'
#' @param x A frame-complete `effect_measurement_form`.
#' @param between The exact `edge_frame()` used to construct `x`.
#' @param method Reconstruction method; `"auto"` selects Parseval, canonical
#'   dual, or projected pseudoinverse according to frame diagnostics.
#' @param tolerance Positive relative singular-value and certification
#'   tolerance.
#' @param max_condition Maximum accepted retained condition number.
#' @param allow_projection Whether rank-deficient projected (lossy)
#'   reconstruction is admitted. It defaults to `FALSE` so a projection is an
#'   explicit choice, never a silent fallback: rank-deficient frames refuse
#'   until projection is explicitly admitted.
#' @param workspace_bytes Positive dense-workspace budget. The default is a
#'   hard 512 MiB ceiling for this explicitly small-node reconstruction path.
#' @param reference_operator Optional finite reference used only to certify the
#'   numerical reconstruction residual.
#' @return An `effect_tomography_result` distinguishing exact, certified, and
#'   projected reconstruction.
#' @export
reconstruct_evidence <- function(
    x, between,
    method = c("auto", "parseval", "canonical_dual",
               "projected_pseudoinverse"),
    tolerance = 1e-10, max_condition = 1e8,
    allow_projection = FALSE, workspace_bytes = 512 * 1024^2,
    reference_operator = NULL) {
  between <- .validate_edge_frame(between)
  .reconstruct_neural_evidence(
    x, between$from_frame, between$to_frame, match.arg(method),
    tolerance, max_condition, allow_projection, workspace_bytes,
    reference_operator
  )
}
