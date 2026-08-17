# Narrow public API for evidence-pairing materializations ------------------

.public_measurement_dense_limit_bytes <- 256 * 1024^2

.require_small_node_measurement_regime <- function(bytes, operation) {
  if (!.is_number(bytes) || bytes < 0) {
    .invariant_error(
      "Measurement resource preflight produced an invalid byte estimate."
    )
  }
  limit <- .public_measurement_dense_limit_bytes
  if (bytes > limit) {
    .input_error(sprintf(
      paste0(
        "%s requires approximately %.0f dense payload bytes, exceeding ",
        "the version 0.1 small-node limit of %.0f bytes. Use the ",
        "support-local geometry plan; brain-scale measurement tomography ",
        "is not yet exported."
      ),
      operation, bytes, limit
    ))
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
    .input_error(
      "Custom measurements require uniquely named operator matrices."
    )
  }
  .measurement_identifier(id, "Measurement-frame `id`")
  legs <- Map(function(operator, node) {
    if (inherits(operator, "Matrix")) operator <- as.matrix(operator)
    if (!.is_finite_matrix(operator) || nrow(operator) < 1L ||
        ncol(operator) != domain$n_features) {
      .input_error(paste0(
        "Every custom measurement must be a finite matrix with one column ",
        "per neural feature."
      ))
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
        .input_error(
          "Measurement units must be scalar or match every operator row."
        )
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
#' @return An `effect_measurement_frame` with one oriented `$legs` entry per
#'   node, the `$node_ids` naming them, the stacked `$frame_operator`, and
#'   `$coverage`, `$injectivity`, and `$dual` diagnostics used by
#'   [reconstruct_evidence()].
#' @seealso [edge_frame()] to request node pairs, [measurement_form()] to
#'   evaluate them, and [additive_frame()] or [compile_frame()] for the
#'   additive frames this can adapt.
#' @family neural domains and frames
#' @examples
#' # Two scalar regional measurements, each a fixed oriented row of weights.
#' native <- abstract_domain(4, id = "demo:native:v1")
#' nodes <- measurement_frame(
#'   list(
#'     anterior = matrix(c(1, 0, 0, 0), 1),
#'     posterior = matrix(c(0, 1, 0, 0), 1)
#'   ),
#'   domain = native, id = "demo:regional-means:v1"
#' )
#' nodes$node_ids
#' nodes$legs$anterior$operator
#'
#' # An additive frame can be adapted instead. In coherent/configuration
#' # mode each measurement keeps its weighted-mean direction and the
#' # orthogonal remainder as separate components.
#' additive <- additive_frame(
#'   matrix(c(1, 2, 1, 0, 0, 1, 2, 1), 2, 4, byrow = TRUE), domain = native
#' )
#' decomposed <- measurement_frame(additive, mode = "coherent_configuration")
#' decomposed$node_ids
#'
#' # An additive frame already names its domain, so supplying one is refused.
#' refused <- try(
#'   measurement_frame(additive, domain = native), silent = TRUE
#' )
#' conditionMessage(attr(refused, "condition"))
#' @export
measurement_frame <- function(
    x, domain = NULL,
    mode = c("total", "coherent", "coherent_configuration"),
    id = "fixed-measurements:v1", units = "arbitrary") {
  mode <- match.arg(mode)
  if (inherits(x, "effect_frame")) {
    if (!is.null(domain)) {
      .input_error("An additive frame already identifies its neural domain.")
    }
    .preflight_additive_measurement_frame(x, mode)
    return(if (mode == "coherent_configuration") {
      .measurement_frame_from_additive_decomposition(x)
    } else {
      .measurement_frame_from_additive(x, mode)
    })
  }
  if (mode != "total") {
    .input_error(paste0(
      "Custom operator lists are already explicit measurements; `mode` must ",
      "remain `\"total\"`."
    ))
  }
  if (is.null(domain)) {
    .input_error(
      "Custom measurement operators require their exact neural `domain`."
    )
  }
  .preflight_custom_measurement_frame(x, domain)
  .public_measurement_frame_from_operators(x, domain, id, units)
}

.edge_frame_signature <- function(fields) {
  .sha256_signature(list(
    schema_version = 1L,
    from_frame = fields$from_frame$signature,
    to_frame = fields$to_frame$signature,
    edges = fields$edges$signature
  ))
}

.validate_edge_frame <- function(x) {
  expected <- c("from_frame", "to_frame", "edges", "signature")
  if (!.sealed_fields(x, "effect_edge_frame", expected)) {
    .input_error("`between` must be a canonical `edge_frame()`.")
  }
  from <- .validate_measurement_frame(x$from_frame)
  to <- .validate_measurement_frame(x$to_frame)
  .validate_measurement_edges(x$edges, from, to)
  fields <- x[setdiff(names(x), "signature")]
  .check_signature(
    x$signature, .edge_frame_signature(fields),
    "Edge-frame identity is inconsistent."
  )
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
#' @return An `effect_edge_frame` for the `between` argument, holding the
#'   `$from_frame` and `$to_frame` it draws nodes from, the requested
#'   `$edges` table (`left`, `right`, `weight`), and a `$signature`.
#' @seealso [measurement_frame()] for the nodes, [measurement_form()] which
#'   consumes this edge set, and [reconstruct_evidence()], which needs every
#'   directed pair.
#' @family neural domains and frames
#' @examples
#' native <- abstract_domain(4, id = "demo:native:v1")
#' nodes <- measurement_frame(
#'   list(
#'     anterior = matrix(c(1, 0, 0, 0), 1),
#'     posterior = matrix(c(0, 1, 0, 0), 1)
#'   ),
#'   domain = native, id = "demo:regional-means:v1"
#' )
#'
#' # Request exactly the edges the question needs. A correlation view also
#' # needs the two self-pairs, because they supply its denominator.
#' pairs <- expand.grid(
#'   from = nodes$node_ids, to = nodes$node_ids, stringsAsFactors = FALSE
#' )
#' between <- edge_frame(pairs$from, pairs$to, nodes)
#' between$edges$edges
#'
#' # A seed-to-target edge set is just a shorter list; nothing is created on
#' # your behalf.
#' edge_frame("anterior", "posterior", nodes)$edges$edges
#'
#' # Node names are checked against the frame, so a typo cannot become a
#' # silently missing edge.
#' refused <- try(edge_frame("anterior", "postrior", nodes), silent = TRUE)
#' conditionMessage(attr(refused, "condition"))
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
#' @return An `effect_pair_query` carrying the `$operator`, its bound
#'   `$left_space` and `$right_space`, and `$metadata$evidence_capability`
#'   recording the `sampling_axis` and `construction` that later views check.
#' @seealso [measurement_form()] and [coupling()], which take this as `by`;
#'   [connectivity()], which refuses a rank-one variation query;
#'   [pair_query()] for a fixed query making no variation claim.
#' @family geometry plans and views
#' @examples
#' # Centering eight repeated time points and dividing by n - 1 is the
#' # within-session covariance operation, declared on the time axis.
#' times <- effect_space(paste0("time", 1:8), basis_id = "demo:time:v1")
#' center <- diag(8) - matrix(1 / 8, 8, 8)
#' sample_covariance <- variation_query(
#'   center / 7, times,
#'   sampling_axis = "time", construction = "joint_covariance",
#'   provenance = list(estimator = "centered within session")
#' )
#' sample_covariance$metadata$evidence_capability$sampling_axis
#' sample_covariance$metadata$evidence_capability$construction
#'
#' # Centering removes one direction, so eight time points leave rank seven:
#' # enough repeated variation for a normalized connectivity view.
#' qr(as.matrix(sample_covariance$operator))$rank
#'
#' # A contrast gives the rank-one query `c c'`. It is a valid effect query,
#' # but it retains no repeated variation for connectivity to normalize.
#' direction <- rep(c(-1, 1), each = 4)
#' rank_one <- variation_query(
#'   tcrossprod(direction), times, "time", "joint_covariance"
#' )
#' qr(as.matrix(rank_one$operator))$rank
#'
#' # The operator must be positive semidefinite to be called variation.
#' refused <- try(
#'   variation_query(diag(c(1, -1, rep(1, 6))), times, "time"), silent = TRUE
#' )
#' conditionMessage(attr(refused, "condition"))
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
    .input_error("Experimental query capability metadata is invalid.")
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
#' @return An `effect_measurement_form` over the requested edge set:
#'   `$block_index` names one row per edge (`edge_id`, `left`, `right`, block
#'   dimensions), `$diagnostics` reports the
#'   `experimental_effective_rank` that gates the normalized views,
#'   and `$capabilities`, `$plan`, and `$receipt` record what may be claimed.
#'   Read the blocks with [effect_coupling()] and the other views.
#' @seealso [effect_coupling()], [connectivity()],
#'   [measurement_components()], and [reconstruct_evidence()] for the views;
#'   [coupling()] to derive the same form from an existing geometry plan.
#' @family coupling and connectivity views
#' @examples
#' # Two sessions of six repeated time points over four native features.
#' native <- abstract_domain(4, id = "demo:native:v1")
#' times <- effect_space(paste0("time", 1:6), basis_id = "demo:time:v1")
#' trend <- seq(-1.5, 1.5, length.out = 6)
#' session <- function(shift) cbind(
#'   trend + shift,
#'   -0.8 * trend + c(0.2, -0.1, 0.1, 0, -0.1, -0.1),
#'   sin(seq(shift, pi, length.out = 6)),
#'   cos(seq(shift, pi, length.out = 6))
#' )
#' signals <- relation(
#'   list(session1 = session(0), session2 = session(0.1)),
#'   effects = times, domain = native
#' )
#'
#' # Two scalar node measurements and all four directed node pairs.
#' nodes <- measurement_frame(
#'   list(anterior = matrix(c(1, 0, 0, 0), 1),
#'        posterior = matrix(c(0, 1, 0, 0), 1)),
#'   domain = native, id = "demo:regional-means:v1"
#' )
#' pairs <- expand.grid(
#'   from = nodes$node_ids, to = nodes$node_ids, stringsAsFactors = FALSE
#' )
#'
#' # `by` closes the experimental axis (within-session covariance over time)
#' # and `over` says each session multiplies by itself, which is biased and
#' # must be declared as such.
#' covariance_over_time <- variation_query(
#'   (diag(6) - matrix(1 / 6, 6, 6)) / 5, times, "time", "joint_covariance"
#' )
#' within_session <- pairing(
#'   signals$partitions, signals$partitions, directed = TRUE,
#'   self_pairs = "allow_biased", independence = "not_independent"
#' )
#' form <- measurement_form(
#'   signals, edge_frame(pairs$from, pairs$to, nodes),
#'   covariance_over_time, within_session
#' )
#' form$block_index[, c("edge_id", "left", "right")]
#'
#' # Centering leaves five directions of repeated variation, which is what
#' # the normalized views require.
#' form$diagnostics$experimental_effective_rank
#' effect_coupling(form)$values[["edge_000002"]]
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
    .contract_error(
      "Measurement frames do not belong to their experimental relation sides."
    )
  }
  capability <- .public_query_capability(by)
  if (!.same_effect_space(by$left_space, left$effect_space) ||
      !.same_effect_space(by$right_space, right$effect_space)) {
    .contract_error(
      "`by` does not match the left and right experimental spaces."
    )
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
#' @return An `effect_coupling_result`. `$values` is one entry per edge (a
#'   matrix block for `effect_coupling()` and `covariance_coupling()`, a data
#'   frame of `canonical_correlation` per mode for `canonical_coupling()`,
#'   and a data frame of `geometry_alignment` for `geometry_alignment()`),
#'   with `$edge_index` naming the edges and `$kind`, `$terminology`,
#'   `$normalization_axis`, and `$regularization` recording what is claimed.
#' @seealso [measurement_form()] and [coupling()] to build `x`;
#'   [connectivity()] for the same views behind one capability-checked entry
#'   point; [measurement_components()] for decomposed nodes.
#'
#' @section Refusals:
#' `effect_coupling()` makes no covariance claim and so refuses nothing. Every
#' other view signals an `effect_capability_refusal` in namespace
#' `"coupling_views"`, so [catch_refusal()] can branch on the cause rather than
#' the prose:
#' \itemize{
#'   \item `"certified_repeated_variation"` (reason
#'     `"repeated_variation_not_certified"`) when the form has not established
#'     that its measurement axis repeats;
#'   \item `"nondegenerate_variation"` (reason `"rank_one_variation_axis"`)
#'     when the variation query has effective rank one;
#'   \item `"coherent_joint_covariance"`, carrying
#'     `"joint_covariance_not_certified"` and `"self_blocks_not_validated"` as
#'     applicable, when the form does not certify a joint covariance across
#'     both measurement spaces;
#'   \item `"nondegenerate_self_variance"` (reason
#'     `"self_variance_not_strictly_positive"`) when a scalar edge has no
#'     measured variation to divide by.
#' }
#' `canonical_coupling()` additionally refuses with
#' `"declared_regularization"` when `ridge` is absent.
#' @family coupling and connectivity views
#' @examples
#' # Two sessions of six repeated time points over four native features,
#' # measured at two multivariate populations.
#' native <- abstract_domain(4, id = "demo:native:v1")
#' times <- effect_space(paste0("time", 1:6), basis_id = "demo:time:v1")
#' trend <- seq(-1.5, 1.5, length.out = 6)
#' session <- function(shift) cbind(
#'   trend + shift,
#'   -0.8 * trend + c(0.2, -0.1, 0.1, 0, -0.1, -0.1),
#'   sin(seq(shift, pi, length.out = 6)),
#'   cos(seq(shift, pi, length.out = 6))
#' )
#' signals <- relation(
#'   list(session1 = session(0), session2 = session(0.1)),
#'   effects = times, domain = native
#' )
#' populations <- measurement_frame(
#'   list(anterior = diag(4)[1:2, , drop = FALSE],
#'        posterior = diag(4)[3:4, , drop = FALSE]),
#'   domain = native, id = "demo:populations:v1"
#' )
#' pairs <- expand.grid(
#'   from = populations$node_ids, to = populations$node_ids,
#'   stringsAsFactors = FALSE
#' )
#' form <- measurement_form(
#'   signals, edge_frame(pairs$from, pairs$to, populations),
#'   variation_query(
#'     (diag(6) - matrix(1 / 6, 6, 6)) / 5, times, "time", "joint_covariance"
#'   ),
#'   pairing(
#'     signals$partitions, signals$partitions, directed = TRUE,
#'     self_pairs = "allow_biased", independence = "not_independent"
#'   )
#' )
#' cross <- form$block_index$edge_id[
#'   form$block_index$left == "anterior" &
#'     form$block_index$right == "posterior"
#' ]
#'
#' # The raw block, with no covariance claim attached to it.
#' effect_coupling(form)$values[[cross]]
#'
#' # The same block, now certified as repeated-sample covariance.
#' covariance_coupling(form)$kind
#'
#' # Canonical correlations describe the shared modes in descending order.
#' # `ridge` is recorded because changing it changes the reported values.
#' canonical <- canonical_coupling(form, ridge = 0.05)
#' canonical$values[canonical$values$edge_id == cross, ]
#'
#' # Geometry alignment asks a different question: do the two populations
#' # induce similar geometry over the repeated observations?
#' alignment <- geometry_alignment(form)
#' alignment$values$geometry_alignment[alignment$values$edge_id == cross]
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
#' @return An `effect_gaussian_covariance_model` declaration recording its
#'   `$family`, that it is `$fixed`, the stated `$provenance`, and a
#'   `$signature` carried into the result identity. It performs no fitting
#'   and no goodness-of-fit test.
#' @seealso [connectivity()] with `view = "gaussian_information"`, the only
#'   place this declaration is used, and [canonical_coupling()] for the
#'   correlation spectrum the information is computed from.
#' @family neural metrics
#' @examples
#' # The declaration is required so that a Gaussian information number
#' # cannot be produced without someone stating the model it rests on.
#' model <- gaussian_covariance_model(
#'   list(assumption = "joint Gaussian time observations")
#' )
#' model$family
#' model$provenance$assumption
#'
#' # It records the assumption; it does not test it, so the signature is a
#' # provenance trail rather than evidence of fit.
#' model$fixed
#' @export
gaussian_covariance_model <- function(provenance = list()) {
  .gaussian_covariance_model(provenance)
}

#' Request a validated connectivity view
#'
#' `connectivity()` is the one entry point for the normalized views, and it
#' checks their preconditions before reporting a number: repeated variation of
#' effective rank above one, valid self-blocks, explicit regularization for
#' the canonical and Gaussian views, and an explicit model declaration for
#' Gaussian information.
#'
#' @param x A complete `effect_measurement_form`.
#' @param view One of signed scalar correlation, a canonical spectrum, static
#'   geometry alignment, or Gaussian mutual information.
#' @param ridge,ridge_right Explicit ridge values for canonical or Gaussian
#'   views.
#' @param model A [gaussian_covariance_model()] for Gaussian information.
#' @param units Information units, when applicable.
#' @param tolerance Positive numerical tolerance.
#' @return An `effect_coupling_result` whose `$values` data frame carries one
#'   row per edge with the view's column (`correlation`,
#'   `canonical_correlation` per `mode`, `geometry_alignment`, or
#'   `information` plus `units`), alongside `$kind`, `$normalization_axis`,
#'   `$regularization`, and `$terminology`.
#' @seealso [effect_coupling()] for the uninterpreted block,
#'   [gaussian_covariance_model()] for the declaration Gaussian information
#'   requires, and [measurement_form()] to build `x`.
#' @family coupling and connectivity views
#'
#' @section Refusals:
#' Each precondition signals an `effect_capability_refusal` in namespace
#' `"coupling_views"`, so [catch_refusal()] can branch on the cause instead of
#' the prose: capability `"certified_repeated_variation"` when the form has not
#' established repeated variation along a named sampling axis,
#' `"nondegenerate_variation"` when the variation query has effective rank one,
#' `"declared_regularization"` for a canonical view without `ridge`, and
#' `"declared_gaussian_model"` for Gaussian information without `model` or
#' `ridge`.
#' @examples
#' # Two scalar, oriented regional measurements over six repeated time
#' # points in two sessions.
#' native <- abstract_domain(4, id = "demo:native:v1")
#' times <- effect_space(paste0("time", 1:6), basis_id = "demo:time:v1")
#' trend <- seq(-1.5, 1.5, length.out = 6)
#' session <- function(shift) cbind(
#'   trend + shift,
#'   -0.8 * trend + c(0.2, -0.1, 0.1, 0, -0.1, -0.1),
#'   sin(seq(shift, pi, length.out = 6)),
#'   cos(seq(shift, pi, length.out = 6))
#' )
#' signals <- relation(
#'   list(session1 = session(0), session2 = session(0.1)),
#'   effects = times, domain = native
#' )
#' nodes <- measurement_frame(
#'   list(anterior = matrix(c(1, 0, 0, 0), 1),
#'        posterior = matrix(c(0, 1, 0, 0), 1)),
#'   domain = native, id = "demo:regional-means:v1"
#' )
#' pairs <- expand.grid(
#'   from = nodes$node_ids, to = nodes$node_ids, stringsAsFactors = FALSE
#' )
#' covariance_over_time <- variation_query(
#'   (diag(6) - matrix(1 / 6, 6, 6)) / 5, times, "time", "joint_covariance"
#' )
#' within_session <- pairing(
#'   signals$partitions, signals$partitions, directed = TRUE,
#'   self_pairs = "allow_biased", independence = "not_independent"
#' )
#' form <- measurement_form(
#'   signals, edge_frame(pairs$from, pairs$to, nodes),
#'   covariance_over_time, within_session
#' )
#'
#' # Both nodes are scalar and oriented, so the correlation view returns
#' # ordinary signed Pearson correlations. The generated data are negatively
#' # related, and the sign is recovered.
#' connectivity(form, view = "correlation")$values
#'
#' # Gaussian mutual information additionally needs an explicit model and
#' # information units, because it is a modeling claim, not a rescaling.
#' connectivity(
#'   form, view = "gaussian_information", ridge = 0.05,
#'   model = gaussian_covariance_model(
#'     list(assumption = "joint Gaussian time observations")
#'   ),
#'   units = "bits"
#' )$values
#'
#' # Omitting the model declaration is refused rather than defaulted, and the
#' # refusal names the missing capability.
#' refusal <- catch_refusal(
#'   connectivity(form, view = "gaussian_information", ridge = 0.05)
#' )
#' refusal$capability
#' refusal$reasons
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
#' @return A data frame with one row per crossed component, carrying
#'   `left_component`/`right_component`, the block dimensions `d_left` and
#'   `d_right`, `left_orientation`/`right_orientation`,
#'   `raw_entries_meaningful` (true only when both sides are oriented),
#'   `frobenius_strength`, and `strongest_singular_value`.
#' @seealso [measurement_frame()] with
#'   `mode = "coherent_configuration"`, which creates the decomposition, and
#'   [effect_coupling()] for the undecomposed block.
#' @family coupling and connectivity views
#' @examples
#' # Two overlapping additive measurements, each split into its
#' # weighted-mean (coherent) direction and the orthogonal remainder.
#' native <- abstract_domain(4, id = "demo:native:v1")
#' times <- effect_space(paste0("time", 1:6), basis_id = "demo:time:v1")
#' trend <- seq(-1.5, 1.5, length.out = 6)
#' session <- function(shift) cbind(
#'   trend + shift,
#'   -0.8 * trend + c(0.2, -0.1, 0.1, 0, -0.1, -0.1),
#'   sin(seq(shift, pi, length.out = 6)),
#'   cos(seq(shift, pi, length.out = 6))
#' )
#' signals <- relation(
#'   list(session1 = session(0), session2 = session(0.1)),
#'   effects = times, domain = native
#' )
#' decomposed <- measurement_frame(
#'   additive_frame(
#'     matrix(c(1, 2, 1, 0, 0, 1, 2, 1), 2, 4, byrow = TRUE), domain = native
#'   ),
#'   mode = "coherent_configuration"
#' )
#' form <- measurement_form(
#'   signals,
#'   edge_frame(
#'     decomposed$node_ids[1], decomposed$node_ids[2], decomposed
#'   ),
#'   variation_query(
#'     (diag(6) - matrix(1 / 6, 6, 6)) / 5, times, "time", "joint_covariance"
#'   ),
#'   pairing(
#'     signals$partitions, signals$partitions, directed = TRUE,
#'     self_pairs = "allow_biased", independence = "not_independent"
#'   )
#' )
#'
#' # Four crossed components on one edge. Only coherent-to-coherent has a
#' # fixed orientation on both sides, so only its raw entry is meaningful;
#' # the others are reported as rotation-invariant strengths.
#' components <- measurement_components(form, edge = 1)
#' components[, c(
#'   "left_component", "right_component", "raw_entries_meaningful",
#'   "frobenius_strength"
#' )]
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
#' @return An `effect_tomography_result` with the reconstructed `$operator`,
#'   the `$method` actually used (`"parseval"`, `"canonical_dual"`, or
#'   `"projected_pseudoinverse"`), a `$status` distinguishing exact,
#'   certified, and projected reconstruction, the `$lossless` and
#'   `$certified` flags, `$left_projection`/`$right_projection`, and frame
#'   `$diagnostics`.
#' @seealso [measurement_form()] and [edge_frame()], which must supply every
#'   directed node pair for a lossless claim.
#' @family sampling uncertainty
#'
#' @section Refusal:
#' A form that is not frame complete — a diagonal-only or requested-edge map —
#' signals an `effect_capability_refusal` carrying capability
#' `"complete_edge_set"` in namespace `"tomography"`, with reason
#' `"edge_set_is_not_frame_complete"`. Inspect it with [catch_refusal()].
#' @examples
#' # A Parseval frame: stacking the two node operators gives the identity, so
#' # the local blocks add back with no correction matrix.
#' native <- abstract_domain(4, id = "demo:native:v1")
#' times <- effect_space(paste0("time", 1:6), basis_id = "demo:time:v1")
#' trend <- seq(-1.5, 1.5, length.out = 6)
#' session <- function(shift) cbind(
#'   trend + shift,
#'   -0.8 * trend + c(0.2, -0.1, 0.1, 0, -0.1, -0.1),
#'   sin(seq(shift, pi, length.out = 6)),
#'   cos(seq(shift, pi, length.out = 6))
#' )
#' signals <- relation(
#'   list(session1 = session(0), session2 = session(0.1)),
#'   effects = times, domain = native
#' )
#' halves <- measurement_frame(
#'   list(first_half = diag(4)[1:2, , drop = FALSE],
#'        second_half = diag(4)[3:4, , drop = FALSE]),
#'   native, id = "demo:parseval:v1"
#' )
#' pairs <- expand.grid(
#'   from = halves$node_ids, to = halves$node_ids, stringsAsFactors = FALSE
#' )
#' between <- edge_frame(pairs$from, pairs$to, halves)
#' form <- measurement_form(
#'   signals, between,
#'   variation_query(
#'     (diag(6) - matrix(1 / 6, 6, 6)) / 5, times, "time", "joint_covariance"
#'   ),
#'   pairing(
#'     signals$partitions, signals$partitions, directed = TRUE,
#'     self_pairs = "allow_biased", independence = "not_independent"
#'   )
#' )
#'
#' # Every directed edge is present and the frame has full column rank, so
#' # the global 4-by-4 neural operator is recovered exactly.
#' reconstructed <- reconstruct_evidence(form, between)
#' c(method = reconstructed$method, status = reconstructed$status,
#'   lossless = reconstructed$lossless)
#' round(reconstructed$operator, 3)
#'
#' # Diagonal node blocks alone are not enough: dropping the cross-edges
#' # discards the between-node directions, and the lossless claim is refused.
#' self_only <- edge_frame(
#'   halves$node_ids, halves$node_ids, halves
#' )
#' partial <- measurement_form(
#'   signals, self_only,
#'   variation_query(
#'     (diag(6) - matrix(1 / 6, 6, 6)) / 5, times, "time", "joint_covariance"
#'   ),
#'   pairing(
#'     signals$partitions, signals$partitions, directed = TRUE,
#'     self_pairs = "allow_biased", independence = "not_independent"
#'   )
#' )
#' refusal <- catch_refusal(reconstruct_evidence(partial, self_only))
#' refusal$capability
#' refusal$remedies
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
