# Scientific views of complete geometry -----------------------------------

#' Read a contrast as signed, coherent, configuration, and total evidence
#'
#' This view never removes a mean pattern or reruns an analysis. It reads the
#' exact coherent/configuration decomposition already contained in `x`.
#'
#' @param x An `effect_geometry_plan` or complete `effect_geometry`.
#' @param weights One finite contrast weight per experimental effect. Named
#'   weights are reordered to the relation's declared effect order and must
#'   name every effect exactly once, so naming them is the safe form. Unnamed
#'   weights are accepted positionally, in the order given by
#'   `x$task$left_relation$effects`; the returned `$weights` always carries the
#'   effect names, so print it to confirm the alignment you intended.
#' @param remove_univariate Must be omitted or `FALSE`. Destructive demeaning
#'   is refused: the coherent/configuration/total decomposition already
#'   reports the common spatial mode and its orthogonal remainder as an
#'   exact, non-destructive partition.
#' @return An `effect_contrast_view` with one value per measurement in
#'   `$signed` (the signed contrast of the local weighted mean), `$coherent`,
#'   `$configuration`, and `$total = coherent + configuration`, plus
#'   `$coherence_fraction` (reported only where the raw cross-generalized
#'   components form a nonnegative partition, flagged by
#'   `$coherence_fraction_valid`), the aligned `$weights`, `$index`, and
#'   `$receipt`.
#' @section Structure:
#' Each value element holds one number per spatial measurement, in `$index`
#' order.
#'
#' - `$signed`: the signed contrast of the local weighted mean. It keeps its
#'   sign, so it says which way the effect goes; the energies below cannot.
#' - `$coherent`: the part of the energy carried by the measurement's own
#'   weighted common spatial mode.
#' - `$configuration`: the orthogonal remainder, the pattern part.
#' - `$total`: `$coherent + $configuration`, the crossvalidated energy of the
#'   contrast. Cross-generalized values may be negative.
#' - `$coherence_fraction`: `$coherent / $total`, and `NA` wherever the raw
#'   components do not form a nonnegative partition.
#' - `$coherence_fraction_valid`: `TRUE` exactly where that fraction was
#'   reported.
#' - `$weights`: the contrast, reordered to the relation's effect order and
#'   named. Print it to confirm the alignment.
#' - `$index`: the measurement identifiers, one per value, carried from the
#'   frame's `$index$measurement`.
#' - `$receipt`: the execution receipt for the run that produced the values.
#'
#' Any element not listed here is internal and may change.
#' @seealso [plan_geometry()] to build `x`, [rdm()] and [rsa()] for the other
#'   named views, and [crossnobis()] for the same total under a declared
#'   noise-precision metric.
#' @family geometry plans and views
#' @examples
#' # Where does an animate-versus-inanimate pattern reproduce across runs?
#' example <- example_fmri_effects()
#' plan <- plan_geometry(
#'   example$fit$relation, example$frame,
#'   cross_partitions(
#'     example$fit$relation,
#'     independence = "independent", generalizes_over = "run"
#'   )
#' )
#' effect <- contrast_energy(plan, example$contrast)
#'
#' # The strongest searchlight falls inside the planted signal, and its
#' # energy splits exactly into coherent and configuration parts.
#' peak <- which.max(effect$total)
#' c(
#'   signed = effect$signed[peak],
#'   coherent = effect$coherent[peak],
#'   configuration = effect$configuration[peak],
#'   total = effect$total[peak],
#'   planted = peak %in% example$truth$signal_measurements
#' )
#'
#' # Demeaning the univariate signal away is refused: the decomposition
#' # already separates the common spatial mode from its remainder.
#' refusal <- catch_refusal(
#'   contrast_energy(plan, example$contrast, remove_univariate = TRUE)
#' )
#' refusal$capability
#' refusal$remedies
#' @param ... Passed to the method. The generic dispatches on `x`: a geometry
#'   plan or a complete geometry is read by the default method documented
#'   here, and an `effect_population_result` from [estimate_population()] by
#'   the group-level method in [population_views].
#' @export
contrast_energy <- function(x, ...) UseMethod("contrast_energy")

#' @rdname contrast_energy
#' @export
contrast_energy.default <- function(x, weights, remove_univariate = FALSE,
                                    ...) {
  .check_no_extra_arguments("contrast_energy", ...)
  if (missing(weights)) {
    .input_error(paste0(
      "`weights` is required: pass one finite weight per experimental ",
      "effect, for example `contrast_energy(plan, c(face = 1, house = -1))`. ",
      "Unnamed weights are taken in the relation's declared effect order."
    ),
      arg = "weights", received = "no argument",
      expected = "one finite weight per experimental effect")
  }
  if (!isFALSE(remove_univariate)) {
    .capability_refusal(paste0(
      "`contrast_energy()` does not remove univariate signal: destructive ",
      "demeaning would change the estimand while leaving voxelwise mean ",
      "effects in the residual subspace. The returned view already reports ",
      "`coherent` (the weighted common spatial mode), `configuration` (its ",
      "orthogonal remainder), and `total` as an exact additive partition; ",
      "select the component you mean instead of deleting one."
    ),
      capability = "nondestructive_decomposition",
      namespace = "geometry_views",
      reasons = "univariate_removal_changes_estimand",
      remedies = paste0(
        "Read `$coherent`, `$configuration`, and `$total` from the ",
        "returned contrast view."
      )
    )
  }
  if (inherits(x, "effect_geometry_plan")) {
    .validate_geometry_plan(x)
    weights <- .align_contrast(
      weights, x$task$left_relation$effect_space$coordinates
    )
    return(.run_geometry_compiler(
      x,
      query = bilinear_query(tcrossprod(weights)),
      component = "contrast",
      signed_query = weights
    ))
  }
  if (!inherits(x, "effect_geometry")) {
    .input_error(sprintf(paste0(
      "`x` must be an `effect_geometry_plan` from `plan_geometry()` or a ",
      "complete `effect_geometry` from `materialize_geometry()`; received %s."
    ), .msg_value(x)),
      arg = "x", received = .msg_value(x),
      expected = "an `effect_geometry_plan` or a complete `effect_geometry`")
  }
  .validate_effect_geometry(x, probe = FALSE)
  weights <- .align_contrast(weights, x$effects)
  query <- bilinear_query(tcrossprod(weights))
  total <- drop(query_geometry(x, query, "total")$values)
  coherent <- drop(query_geometry(x, query, "coherent")$values)
  # The projected contrast carries the same route-stable view identity that
  # fused query-first execution derives for identical weights.
  packed_query <- matrix(.svec_symmetric(tcrossprod(weights)), ncol = 1L)
  receipt <- .projection_receipt(
    x$receipt,
    .geometry_view_scientific_id(
      x$receipt$scientific_plan_id, "contrast", packed_query, weights
    )
  )
  .new_effect_contrast_view(
    total, coherent, x$marginals, weights, x$index, receipt, x$metadata
  )
}

# `.new_effect_contrast_view()` is the `effect_contrast_view` record itself and
# now lives in R/result.R with the other result records, because the executor
# builds one directly on the query-first path and had to reach up into this
# view file to do it.

.self_geometry_source <- function(x, operation, complete = FALSE) {
  if (inherits(x, "effect_geometry_plan")) {
    if (isTRUE(complete)) {
      .capability_refusal(sprintf(paste0(
        "%s requires a complete effect form, and `x` is a query-first ",
        "`effect_geometry_plan`: a plan names an estimand but holds no ",
        "geometry to decompose."
      ), operation),
        capability = "complete_geometry",
        namespace = "geometry_views",
        reasons = "query_first_plan_has_no_materialized_geometry",
        remedies = "Call `materialize_geometry(x)` and pass the result."
      )
    }
    .validate_geometry_plan(x)
    if (identical(x$codec, "rectangular")) {
      .capability_refusal(sprintf(paste0(
        "%s requires a symmetric self form; this plan is rectangular (%d x ",
        "%d effects across two relations)."
      ), operation, x$logical_shape[[1L]], x$logical_shape[[2L]]),
        capability = "symmetric_self_form",
        namespace = "geometry_views",
        reasons = "rectangular_cross_axis_plan",
        remedies = paste0(
          "Read a rectangular plan with axis-bound `pair_query()`s through ",
          "`evaluate_geometry()`, or build a self-form plan by omitting ",
          "`right` in `plan_geometry()`."
        )
      )
    }
    space <- x$task$left_relation$effect_space
    return(list(
      kind = "plan",
      effects = space$coordinates,
      codec = "symmetric_packed",
      index = .execution_measurement_index(x$frame)
    ))
  }
  if (!inherits(x, "effect_form")) {
    .input_error(sprintf(paste0(
      "%s requires an `effect_geometry_plan` from `plan_geometry()` or a ",
      "complete effect form from `materialize_geometry()`; received %s."
    ), operation, .msg_value(x)),
      arg = "x", received = .msg_value(x),
      expected = "an `effect_geometry_plan` or a complete `effect_geometry`")
  }
  .validate_effect_form(x, probe = FALSE)
  if (!isTRUE(x$capabilities$self_form) ||
      !isTRUE(x$capabilities$symmetric)) {
    .capability_refusal(sprintf(paste0(
      "%s requires a symmetric self form; this form declares self_form = %s ",
      "and symmetric = %s."
    ), operation, isTRUE(x$capabilities$self_form),
      isTRUE(x$capabilities$symmetric)),
      capability = "symmetric_self_form",
      namespace = "geometry_views",
      reasons = "form_is_not_a_symmetric_self_form",
      remedies = paste0(
        "Materialize a self-form plan (one relation, no `right` argument), ",
        "or read a cross-axis form with `pair_query()` through ",
        "`evaluate_geometry()`."
      )
    )
  }
  list(
    kind = "form",
    effects = x$effects,
    codec = x$codec,
    index = x$index
  )
}

.require_self_form_view <- function(x, operation) {
  .self_geometry_source(x, operation, complete = TRUE)
  invisible(x)
}

#' Read squared experimental distances from geometry
#'
#' For effects `i` and `j`, this view applies the fixed contrast
#' `c = e_i - e_j` to the selected geometry component `G`:
#' \deqn{d_{ij} = c^T G c = G_{ii} + G_{jj} - 2G_{ij}.}{d_ij = c^T G c = G_ii + G_jj - 2 G_ij.}
#' With cross-partition geometry this is the signed crossvalidated squared
#' Euclidean distance, or squared Mahalanobis distance when the plan carries a
#' fixed neural metric. It is not `1 - Pearson correlation`: correlation
#' distance requires data-dependent diagonal normalization and is outside this
#' linear view.
#'
#' The returned values are point estimates. For equal-weight all-partition-
#' pairs geometry with a common fixed metric, [rdm_sampling_covariance()] can
#' construct a separate analytic within-measurement covariance law when the
#' plan was built from `lm_relation_fit()` and its residual error channel is
#' still available. Pair rows that share a partition are not independent
#' replicates for a spread-across-pairs standard error.
#'
#' @param x An `effect_geometry_plan` or a complete effect form carrying the
#'   symmetric self-form capability.
#' @param component One of `total`, `coherent`, or `configuration`.
#' @param pairs Optional two-column matrix selecting the effect pairs to
#'   report, by effect name or index. The default reports every unordered
#'   pair. Selected pairs execute without materializing the remaining
#'   geometry: the RDM is a view, not a mandatory intermediate object.
#' @param normalize Must be omitted. Correlation-style diagonal normalization
#'   of a signed cross-generalized form is refused: crossvalidated diagonals
#'   can be zero or negative, so `1 - r` here is not conventional Pearson
#'   distance. The boundary is documented in the correlation-distance policy.
#' @return An `effect_rdm_view`. `$values` has one row per spatial
#'   measurement and one column per requested experimental pair, `$pairs` is
#'   the `left`/`right` table naming those columns, and `$component`,
#'   `$index`, and `$receipt` record what was read. Cross-generalized
#'   distances may be negative.
#' @section Structure:
#' The distances are one measurement-by-pair matrix; the other elements name
#' its axes and record what was read.
#'
#' - `$values`: one row per spatial measurement, one column per requested
#'   pair, in `$pairs` row order.
#' - `$pairs`: a data frame whose `left` and `right` columns name the two
#'   effects behind each column of `$values`.
#' - `$component`: the geometry component the distances were taken from.
#' - `$index`: the measurement identifiers, one per row of `$values`,
#'   carried from the frame's `$index$measurement`.
#' - `$receipt`: the execution receipt for the run that produced the values.
#'
#' Any element not listed here is internal and may change.
#' @seealso [rsa()] to regress model RDMs on these distances,
#'   [contrast_energy()] for a single contrast, and
#'   [rdm_sampling_covariance()] for the admitted analytic uncertainty law.
#' @family geometry plans and views
#'
#' @section Refusals:
#' `normalize` signals an `effect_capability_refusal` with capability
#' `"guaranteed_psd"`, and a rectangular cross-axis plan or a non-symmetric
#' form signals capability `"symmetric_self_form"`, both in namespace
#' `"geometry_views"`. Branch on them with [catch_refusal()].
#' @examples
#' # Three conditions over two regions, generalizing across two runs.
#' domain <- abstract_domain(4, id = "rdm-example")
#' run1 <- rbind(
#'   face = c(1, 0.2, 0, 0), house = c(0, 1, 0.1, 0), tool = c(0, 0, 1, 0.3)
#' )
#' run2 <- rbind(
#'   face = c(0.9, 0.3, 0, 0), house = c(0.1, 0.9, 0, 0), tool = c(0, 0.1, 1.1, 0.2)
#' )
#' relation <- relation(list(run1 = run1, run2 = run2), domain = domain)
#' plan <- plan_geometry(
#'   relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
#'   cross_partitions(
#'     relation, independence = "independent", generalizes_over = "run"
#'   )
#' )
#'
#' # All three unordered pairs, one column each.
#' distances <- rdm(plan)
#' distances
#' distances$pairs
#'
#' # Selecting pairs is a narrower view, not a post-hoc subset: the remaining
#' # geometry is never computed.
#' rdm(plan, pairs = cbind("face", "house"))$values
#'
#' # Correlation-style normalization is refused, because crossvalidated
#' # diagonals can be zero or negative.
#' refusal <- catch_refusal(rdm(plan, normalize = "correlation"))
#' refusal$capability
#' @param ... Passed to the method. The generic dispatches on `x`: a geometry
#'   plan or a complete effect form is read by the default method documented
#'   here, and an `effect_population_result` from [estimate_population()] by
#'   the group-level method in [population_views].
#' @export
rdm <- function(x, ...) UseMethod("rdm")

#' @rdname rdm
#' @export
rdm.default <- function(x, component = c("total", "coherent", "configuration"),
                        pairs = NULL, normalize = NULL, ...) {
  .check_no_extra_arguments("rdm", ...)
  if (!is.null(normalize)) {
    .capability_refusal(paste0(
      "`rdm()` reports signed squared distances and will not apply ",
      "correlation-style normalization: crossvalidated diagonal estimates ",
      "can be zero or negative, so dividing by them is not conventional ",
      "`1 - Pearson` distance and would silently change the estimand. ",
      "Conventional correlation distance requires a guaranteed ",
      "positive-semidefinite self form and its own named view; see the ",
      "correlation-distance policy vignette."
    ),
      capability = "guaranteed_psd",
      namespace = "geometry_views",
      reasons = "signed_cross_generalized_diagonals",
      remedies = paste0(
        "Use the signed squared-distance RDM, or wait for the explicit ",
        "`correlation_rdm()` view with its own contract."
      )
    )
  }
  source <- .self_geometry_source(x, "RDM")
  component <- match.arg(component)
  query <- .pair_difference_query(source$effects, pairs = pairs)
  view <- if (source$kind == "plan") {
    evaluate_geometry(x, query = query, component = component)
  } else {
    query_geometry(x, query, component)
  }
  structure(
    list(
      values = view$values,
      pairs = data.frame(
        left = source$effects[query$pair_left],
        right = source$effects[query$pair_right],
        stringsAsFactors = FALSE
      ),
      component = component,
      index = view$index,
      receipt = view$receipt
    ),
    class = "effect_rdm_view"
  )
}

.validate_rdm_models <- function(models, effects, label) {
  q <- length(effects)
  if (is.null(models)) return(list())
  if (is.matrix(models)) models <- list(model = models)
  if (!is.list(models) || length(models) < 1L) {
    .input_error(sprintf(paste0(
      "`%s` must be one dissimilarity matrix or a nonempty named list of ",
      "them; received %s."
    ), label, .msg_value(models)),
      arg = label, received = .msg_value(models),
      expected = "one dissimilarity matrix, or a named list of them")
  }
  if (is.null(names(models)) || anyNA(names(models)) ||
      any(!nzchar(names(models))) || anyDuplicated(names(models))) {
    .input_error(sprintf(paste0(
      "`%s` must be a list with unique nonempty names; the names become the ",
      "coefficient columns of the fit."
    ), label),
      arg = label, received = .msg_names(names(models)),
      expected = "a list with unique nonempty names")
  }
  entries <- names(models)
  out <- lapply(entries, function(entry) {
    value <- models[[entry]]
    where <- sprintf("`%s` RDM `%s`", label, entry)
    if (!is.matrix(value) || !is.numeric(value)) {
      .input_error(sprintf("%s must be a numeric matrix; received %s.",
        where, .msg_value(value)),
        arg = label, received = .msg_value(value),
        expected = "a numeric matrix")
    }
    if (!identical(dim(value), c(q, q))) {
      .input_error(sprintf(paste0(
        "%s is %d x %d; the relation declares %s (%s), so every model RDM ",
        "must be %d x %d."
      ), where, nrow(value), ncol(value), .msg_count(q, "effect"),
        .msg_names(effects), q, q),
        arg = label,
        received = sprintf("%d x %d", nrow(value), ncol(value)),
        expected = sprintf("%d x %d", q, q))
    }
    if (any(!is.finite(value))) {
      .input_error(sprintf("%s contains %s non-finite %s.", where,
        sum(!is.finite(value)),
        if (sum(!is.finite(value)) == 1L) "entry" else "entries"),
        arg = label,
        received = sprintf("%d non-finite of %d", sum(!is.finite(value)),
          length(value)),
        expected = "all entries finite")
    }
    asymmetry <- max(abs(value - t(value)))
    if (asymmetry > 1e-12) {
      .input_error(sprintf(paste0(
        "%s is not symmetric; the largest difference between `m[i, j]` and ",
        "`m[j, i]` is %g. A dissimilarity between two effects has one value."
      ), where, asymmetry),
        arg = label, received = sprintf("asymmetry %g", asymmetry),
        expected = "a symmetric matrix")
    }
    if (max(abs(diag(value))) > 1e-12) {
      .input_error(sprintf(paste0(
        "%s has a nonzero diagonal (largest |m[i, i]| is %g). Pass a ",
        "dissimilarity matrix, not a similarity matrix: an effect is at ",
        "distance zero from itself."
      ), where, max(abs(diag(value)))),
        arg = label,
        received = sprintf("largest |m[i, i]| is %g", max(abs(diag(value)))),
        expected = "a zero diagonal")
    }
    row_ids <- rownames(value)
    column_ids <- colnames(value)
    if (!is.null(row_ids) || !is.null(column_ids)) {
      if (is.null(row_ids) || is.null(column_ids)) {
        .input_error(sprintf(paste0(
          "%s names only its %s; name both axes with the relation's effects ",
          "(%s), or neither."
        ), where, if (is.null(row_ids)) "columns" else "rows",
          .msg_names(effects)),
          arg = label,
          received = if (is.null(row_ids)) "column names only" else
            "row names only",
          expected = "both axes named with the relation's effects, or neither")
      }
      if (anyDuplicated(row_ids) || anyDuplicated(column_ids) ||
          anyNA(row_ids) || anyNA(column_ids) ||
          !setequal(row_ids, effects) || !setequal(column_ids, effects)) {
        unknown <- setdiff(unique(c(row_ids, column_ids)), effects)
        absent <- setdiff(effects, intersect(row_ids, column_ids))
        detail <- c(
          if (length(unknown)) sprintf("%s is not a declared effect",
            .msg_names(unknown)),
          if (length(absent)) sprintf("%s is missing from an axis",
            .msg_names(absent))
        )
        if (!length(detail)) detail <- "an axis repeats an effect"
        .input_error(sprintf(
          "%s axis names do not match the relation's effects (%s): %s.",
          where, .msg_names(effects), paste(detail, collapse = "; ")),
          arg = label, received = paste(detail, collapse = "; "),
          expected = .msg_names(effects))
      }
      value <- value[effects, effects, drop = FALSE]
    }
    value
  })
  names(out) <- entries
  out
}

# Name the columns that make an RSA design singular instead of asking the
# caller to guess. The pivoting QR already separates the retained basis from
# the dependent columns; regressing each dependent column on that basis says
# which retained terms reproduce it, which is the sentence a reader can act on.
.rsa_rank_deficiency_message <- function(design, qr_design, intercept) {
  labels <- colnames(design)
  rank <- qr_design$rank
  pivot <- qr_design$pivot
  if (rank < 1L) {
    return(paste0(
      "The RSA design is rank deficient: every column of the pair-space ",
      "design is zero, so no coefficient is identified. Supply at least one ",
      "model RDM with a nonzero off-diagonal entry."
    ))
  }
  retained <- labels[pivot[seq_len(rank)]]
  dependent <- labels[pivot[seq.int(rank + 1L, length(pivot))]]
  basis <- design[, retained, drop = FALSE]
  explanations <- vapply(dependent, function(term) {
    coefficients <- tryCatch(
      qr.solve(basis, design[, term], tol = 1e-10),
      error = function(condition) rep(NA_real_, length(retained))
    )
    contributors <- retained[is.finite(coefficients) &
        abs(coefficients) > 1e-8]
    if (!length(contributors)) {
      sprintf("`%s` is zero in pair space", term)
    } else {
      sprintf("`%s` is an exact linear combination of %s", term,
        .msg_names(contributors))
    }
  }, character(1))
  intercept_only <- isTRUE(intercept) &&
    all(vapply(dependent, function(term) {
      coefficients <- tryCatch(
        qr.solve(basis, design[, term], tol = 1e-10),
        error = function(condition) rep(NA_real_, length(retained))
      )
      contributors <- retained[is.finite(coefficients) &
          abs(coefficients) > 1e-8]
      length(contributors) > 0L && all(contributors == "(Intercept)")
    }, logical(1)))
  remedy <- if (intercept_only) {
    paste0(
      "The intercept column is added automatically, so a model RDM that is ",
      "constant off the diagonal duplicates it. Pass `intercept = FALSE` to ",
      "fit these RDMs without the constant column, or drop the redundant ",
      "model."
    )
  } else {
    paste0(
      "Remove one of the redundant RDMs, or combine the collinear models ",
      "into a single predictor."
    )
  }
  sprintf(
    "The RSA design is rank deficient (rank %d of %s): %s. %s",
    rank, .msg_count(ncol(design), "column"),
    paste(explanations, collapse = "; "), remedy
  )
}

.rdm_vector <- function(value) {
  pairs <- utils::combn(seq_len(nrow(value)), 2L)
  value[cbind(pairs[1L, ], pairs[2L, ])]
}

# Return X^+ without constructing an n-by-n identity as the response to
# qr.coef(). For X[, pivot] = Q R, the coefficient map in pivoted order is
# R^-1 Q^T; assigning those rows through qr$pivot restores model-term order.
.thin_qr_coefficient_map <- function(qr_design) {
  q <- qr.Q(qr_design, complete = FALSE)
  r <- qr.R(qr_design, complete = FALSE)
  pivoted <- backsolve(r, t(q))
  out <- matrix(0, nrow(pivoted), ncol(pivoted))
  out[qr_design$pivot, ] <- pivoted
  dimnames(out) <- list(colnames(qr_design$qr), NULL)
  out
}

#' Fit multiple-regression RSA as one compiled geometry query
#'
#' `rsa()` compiles the RDM transform and the least-squares coefficient map
#' into a single fixed query, so the regression is executed as one pass over
#' the plan rather than as a second analysis of a stored RDM. Model rows and
#' columns are aligned to the relation's effect names before any geometry is
#' read.
#'
#' @param x An `effect_geometry_plan` or a complete effect form carrying the
#'   symmetric self-form capability.
#' @param models Named model RDMs, each a finite symmetric zero-diagonal
#'   matrix over the experimental effects. Row and column names are optional;
#'   when supplied they must identify every effect exactly once and are
#'   reordered to the relation's effect order.
#' @param nuisance Optional named nuisance RDMs, in the same form.
#' @param intercept Whether to include an intercept in RDM space. It is `TRUE`
#'   by default, so a model RDM that is constant off the diagonal is collinear
#'   with it; the rank-deficiency message names the columns involved and
#'   `intercept = FALSE` fits the same models without the constant column.
#' @param component Geometry component to read.
#' @return An `effect_rsa_view`. `$coefficients` has one row per measurement
#'   and one named column per model, nuisance model, and the optional
#'   intercept; `$component`, `$index`, and `$receipt` record what was read.
#' @section Structure:
#' The fit is one measurement-by-term coefficient matrix; the other elements
#' name its axes and record what was read.
#'
#' - `$coefficients`: one row per spatial measurement, one named column per
#'   design term, in `$terms` row order.
#' - `$terms`: a data frame naming each column of `$coefficients` in `term`
#'   and labeling it `intercept`, `model`, or `nuisance` in `role`.
#' - `$component`: the geometry component the models were regressed on.
#' - `$index`: the measurement identifiers, one per row of `$coefficients`,
#'   carried from the frame's `$index$measurement`.
#' - `$receipt`: the execution receipt for the run that produced the fit.
#'
#' The compiled `$query` and any other element not listed here are internal
#' and may change.
#' @seealso [rdm()] for the distances the regression is fitted to, and
#'   [plan_geometry()] for the plan.
#' @family geometry plans and views
#'
#' @section Refusal:
#' A rectangular cross-axis plan or a non-symmetric form signals an
#' `effect_capability_refusal` with capability `"symmetric_self_form"` in
#' namespace `"geometry_views"`; see [catch_refusal()].
#' @examples
#' domain <- abstract_domain(4, id = "rsa-example")
#' run1 <- rbind(
#'   face = c(1, 0.2, 0, 0), house = c(0, 1, 0.1, 0), tool = c(0, 0, 1, 0.3)
#' )
#' run2 <- rbind(
#'   face = c(0.9, 0.3, 0, 0), house = c(0.1, 0.9, 0, 0), tool = c(0, 0.1, 1.1, 0.2)
#' )
#' relation <- relation(list(run1 = run1, run2 = run2), domain = domain)
#' plan <- plan_geometry(
#'   relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
#'   cross_partitions(relation, independence = "independent")
#' )
#'
#' # An animacy model: the two animate-inanimate pairs are far, the rest near.
#' conditions <- rownames(run1)
#' animacy <- matrix(
#'   c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3,
#'   dimnames = list(conditions, conditions)
#' )
#' fit <- rsa(plan, models = list(animacy = animacy))
#' round(fit$coefficients, 3)
#' as.data.frame(fit)
#'
#' # A model must be a dissimilarity matrix. Passing a similarity matrix,
#' # whose diagonal is nonzero, is rejected before geometry is read.
#' similarity <- 1 - animacy
#' wrong <- try(rsa(plan, models = list(animacy = similarity)), silent = TRUE)
#' conditionMessage(attr(wrong, "condition"))
#' @param ... Passed to the method. The generic dispatches on `x`: a geometry
#'   plan or a complete effect form is read by the default method documented
#'   here, and an `effect_population_result` from [estimate_population()] by
#'   the group-level method in [population_views].
#' @export
rsa <- function(x, ...) UseMethod("rsa")

#' @rdname rsa
#' @export
rsa.default <- function(x, models, nuisance = NULL, intercept = TRUE,
                        component = c("total", "coherent", "configuration"),
                        ...) {
  .check_no_extra_arguments("rsa", ...)
  source <- .self_geometry_source(x, "RSA")
  if (missing(models)) {
    .input_error(paste0(
      "`models` is required: pass one dissimilarity matrix over the ",
      "relation's effects, or a named list of them, for example ",
      "`rsa(plan, models = list(category = m))`."
    ),
      arg = "models", received = "no argument",
      expected = "one dissimilarity matrix, or a named list of them")
  }
  if (!.is_flag(intercept)) {
    .input_error(sprintf("`intercept` must be TRUE or FALSE; received %s.",
      .msg_value(intercept)),
      arg = "intercept", received = .msg_value(intercept),
      expected = "TRUE or FALSE")
  }
  component <- match.arg(component)
  q <- length(source$effects)
  models <- .validate_rdm_models(models, source$effects, "models")
  nuisance <- .validate_rdm_models(nuisance, source$effects, "nuisance")
  if (any(names(models) %in% names(nuisance))) {
    .input_error(sprintf(paste0(
      "Model and nuisance names must be distinct; %s appears in both. Each ",
      "name becomes one coefficient column."
    ), .msg_names(intersect(names(models), names(nuisance)))),
      arg = "nuisance",
      received = .msg_names(intersect(names(models), names(nuisance))),
      expected = "names distinct from `models`")
  }
  predictors <- c(models, nuisance)
  design <- do.call(cbind, lapply(predictors, .rdm_vector))
  colnames(design) <- names(predictors)
  roles <- c(rep("model", length(models)), rep("nuisance", length(nuisance)))
  if (intercept) {
    design <- cbind(`(Intercept)` = 1, design)
    roles <- c("intercept", roles)
  }
  qr_design <- qr(design, LAPACK = FALSE)
  if (qr_design$rank != ncol(design)) {
    .input_error(.rsa_rank_deficiency_message(design, qr_design, intercept),
      arg = "models", received = "a rank-deficient design",
      expected = "model and nuisance RDMs of full column rank")
  }
  # The OLS coefficient map is a fixed linear readout of pair space: the
  # compiled query stays in structured pair-difference form, so no packed
  # q(q+1)/2-by-pairs matrix is ever materialized.
  coefficient_map <- .thin_qr_coefficient_map(qr_design)
  rownames(coefficient_map) <- colnames(design)
  query <- .pair_difference_query(
    source$effects, coefficients = coefficient_map
  )
  view <- if (source$kind == "plan") {
    evaluate_geometry(x, query = query, component = component)
  } else {
    query_geometry(x, query, component)
  }

  structure(
    list(
      coefficients = view$values,
      terms = data.frame(term = colnames(design), role = roles,
        stringsAsFactors = FALSE),
      component = component,
      index = view$index,
      receipt = view$receipt,
      query = query
    ),
    class = "effect_rsa_view"
  )
}

#' Read the signed eigenvalue spectrum of cross-generalized geometry
#'
#' Unlike the linear views, eigenvalues are not additive across components:
#' the spectra of the `coherent` and `configuration` components do not sum
#' to the spectrum of `total`, even though the underlying matrices do.
#' Compare spectra across components only as separate decompositions.
#'
#' @param x A complete effect form carrying the symmetric self-form
#'   capability, as returned by [materialize_geometry()]. A query-only view
#'   has no geometry to decompose.
#' @param component Geometry component to decompose.
#' @param row_block Positive number of measurement rows read per block.
#' @return An `effect_spectrum_view`. `$values` has one row per measurement
#'   and one column per eigenvalue (`root1` largest), with `$component`,
#'   `$index`, `$receipt`, and `$indefinite_estimates_preserved = TRUE`
#'   recording that negative eigenvalues are never truncated at zero.
#' @section Structure:
#' One signed spectrum per spatial measurement, with the labels that say what
#' was decomposed.
#'
#' - `$values`: one row per measurement, one column per eigenvalue, ordered
#'   `root1` (largest) through `rootq`. Negative roots are retained.
#' - `$component`: the geometry component that was decomposed.
#' - `$index`: the measurement identifiers, one per row of `$values`, carried
#'   from the decomposed form.
#' - `$indefinite_estimates_preserved`: always `TRUE`, recording that no
#'   eigenvalue was truncated at zero.
#' - `$receipt`: the execution receipt of the form that was decomposed.
#'
#' Any element not listed here is internal and may change.
#' @seealso [materialize_geometry()], which produces the complete geometry
#'   this view requires, and [rdm()] for the linear distance view.
#' @family geometry plans and views
#'
#' @section Refusal:
#' Passing a query-first `effect_geometry_plan` signals an
#' `effect_capability_refusal` with capability `"complete_geometry"` in
#' namespace `"geometry_views"` and remedy `materialize_geometry(x)`; a
#' non-symmetric form signals capability `"symmetric_self_form"`. See
#' [catch_refusal()].
#' @examples
#' domain <- abstract_domain(4, id = "spectrum-example")
#' relation <- relation(
#'   list(run1 = rbind(a = c(1, 0, 2, 1), b = c(0, 1, 1, 0)),
#'        run2 = rbind(a = c(1.1, 0.1, 1.9, 0.8), b = c(0.1, 0.9, 1.2, 0.2))),
#'   domain = domain
#' )
#' geometry <- materialize_geometry(plan_geometry(
#'   relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
#'   cross_partitions(relation, independence = "independent")
#' ))
#'
#' # Signed eigenvalues per region, largest first. Small negative roots are
#' # retained: a cross-generalized form is not guaranteed positive.
#' spectrum <- geometry_spectrum(geometry)
#' spectrum
#' as.data.frame(spectrum)
#'
#' # Eigenvalues are not additive across components, so read each spectrum as
#' # its own decomposition rather than summing them.
#' geometry_spectrum(geometry, component = "coherent")$values
#' @export
geometry_spectrum <- function(x,
                              component = c("total", "coherent", "configuration"),
                              row_block = 1024L) {
  .require_self_form_view(x, "A geometry spectrum")
  component <- match.arg(component)
  row_block <- .validate_tile_size(row_block, "row_block")
  .validate_effect_form(x)
  .require_effect_form_component(x, component)
  q <- length(x$effects)
  values <- matrix(0, x$total$dim[[1L]], q,
    dimnames = list(NULL, paste0("root", seq_len(q))))
  for (start in .tile_starts(nrow(values), row_block)) {
    rows <- start:min(start + row_block - 1L, nrow(values))
    packed <- .geometry_component_validated(x, component, rows)
    for (position in seq_along(rows)) {
      form <- if (identical(x$codec, "symmetric_packed")) {
        .unsvec_symmetric(packed[position, ], q)
      } else {
        matrix(packed[position, ], q, q)
      }
      values[rows[[position]], ] <- eigen(
        form, symmetric = TRUE,
        only.values = TRUE
      )$values
    }
  }
  structure(
    list(values = values, component = component, index = x$index,
      indefinite_estimates_preserved = TRUE, receipt = x$receipt),
    class = "effect_spectrum_view"
  )
}

# `contribution()` -- additive aggregation over a territory -------------------
#
# The attribution reading of a conservative frame
# (`design/conservative-geometry-contract.md` section 1.2) says a node's value
# is its share of one fixed global budget. Reading a ledger means adding shares
# up over a territory, and this is the operation that does it. Four decisions,
# each forced by the contract rather than chosen for convenience:
#
#  1. Grouping is BY ROW. Every row of the view belongs to exactly one group,
#     so the group sums partition the budget exactly. Gap G4 of section 11.4
#     names the alternative -- splitting an overlapping node's mass across
#     several regions -- and it is deliberately not offered: a row is one node,
#     a node is not divisible, and any split needs a second partition of the
#     weights that can double-count. If overlap-splitting is ever wanted it has
#     to arrive as its own declared, separately certified reduction.
#  2. `total` is budget-exact: the group sums add back to the whole-domain
#     total (section 2, claim 2). The signed marginal is NOT, and this is the
#     one place the arithmetic is genuinely asymmetric. `$signed` is the local
#     weighted *mean* contrast -- `.ordered_pairing_marginals()` divides the
#     node's first moment by the node's own frame mass -- so it is intensive,
#     exactly the kind of quantity section 1.1 says must never be summed. A
#     mass-weighted mean would be the defensible aggregate and a view carries
#     no masses, so the aggregate reports `NA` and says why, on the same
#     masking discipline the coherence fraction uses.
#  3. `coherent` and `configuration` add up too, but their sums are
#     FRAME-RELATIVE (section 4, claim 4): `sum_x G_x^coh` is not a global
#     quantity, so a coherent budget is a share of that frame's coherent mass
#     and two frames give two incomparable denominators. The result records
#     `frame_relative` and the print says it.
#  4. Coherence fractions are recomputed, never averaged. A fraction of sums is
#     not a sum of fractions, so the group fraction comes from the aggregated
#     components under `.coherence_fraction()`, the same nonnegative-partition
#     mask a node's fraction gets.

.contribution_measurement_ids <- function(index) {
  if (is.data.frame(index)) {
    if ("measurement" %in% names(index)) {
      return(as.character(index$measurement))
    }
    return(as.character(seq_len(nrow(index))))
  }
  as.character(index)
}

# What `contribution()` accepts, and the reasoned refusal for everything else.
.contribution_source <- function(x) {
  if (inherits(x, "effect_contrast_view")) {
    return(list(kind = "contrast", measurements = length(x$total)))
  }
  if (inherits(x, "effect_view") && !inherits(x, "effect_form")) {
    .validate_effect_view(x)
    return(list(kind = "view", measurements = nrow(x$values)))
  }
  .contribution_unsupported(x)
}

.contribution_supported_sentence <- function() {
  paste0(
    "`contribution()` aggregates the two per-measurement views whose values ",
    "are additive contributions to one budget: an `effect_contrast_view` ",
    "from `contrast_energy()`, and a query-only `effect_view` from ",
    "`evaluate_geometry()` or `query_geometry()`."
  )
}

.contribution_unsupported <- function(x) {
  supported <- .contribution_supported_sentence()
  if (inherits(x, "effect_spectrum_view")) {
    .capability_refusal(paste0(
      "`contribution()` will not aggregate an `effect_spectrum_view`. ",
      "Eigenvalues are not additive across measurements: a territory's ",
      "spectrum is not the sum of its nodes' spectra, so the sums would not ",
      "be the spectrum of anything. ", supported
    ),
      capability = "additive_contribution",
      namespace = "geometry_views",
      reasons = "eigenvalues_are_not_additive",
      remedies = paste0(
        "Aggregate a linear read of the same geometry -- `contrast_energy()`, ",
        "or `evaluate_geometry(plan, query = ...)` -- and decompose the ",
        "aggregate afterwards if a spectrum of the territory is the question."
      )
    )
  }
  comparative <- c("effect_rdm_view", "effect_rsa_view",
    "effect_crossnobis_view")
  if (any(comparative %in% class(x))) {
    kind <- comparative[[which(comparative %in% class(x))[[1L]]]]
    .capability_refusal(sprintf(paste0(
      "`contribution()` does not aggregate an `%s`. Each row is read as a ",
      "dissimilarity between conditions or a regression coefficient on one ",
      "-- a comparison, not a share of a declared whole -- and summing those ",
      "over a territory would present a budget for an object nobody declared ",
      "to be one. %s"
    ), kind, supported),
      capability = "additive_contribution",
      namespace = "geometry_views",
      reasons = "comparative_readout_is_not_a_budget",
      remedies = paste0(
        "Run the same fixed query through `evaluate_geometry(plan, query = ",
        "...)`, which returns those numbers as an additive `effect_view`, ",
        "and aggregate that instead."
      )
    )
  }
  if (inherits(x, "effect_form") || inherits(x, "effect_geometry")) {
    .capability_refusal(paste0(
      "`contribution()` aggregates a readout, not packed geometry: a stored ",
      "form holds every effect pair at once and has no single per-measurement ",
      "value to add up. ", supported
    ),
      capability = "additive_contribution",
      namespace = "geometry_views",
      reasons = "packed_geometry_is_not_a_readout",
      remedies = paste0(
        "Read one component first -- `contrast_energy(x, weights)` or ",
        "`query_geometry(x, query, component)` -- then aggregate the view."
      )
    )
  }
  if (inherits(x, "effect_geometry_plan")) {
    .capability_refusal(paste0(
      "`contribution()` aggregates values, and an `effect_geometry_plan` ",
      "names an estimand without holding any. ", supported
    ),
      capability = "additive_contribution",
      namespace = "geometry_views",
      reasons = "plan_holds_no_values",
      remedies = paste0(
        "Evaluate the plan first -- `contrast_energy(plan, weights)` or ",
        "`evaluate_geometry(plan, query = ...)` -- then aggregate the view."
      )
    )
  }
  .input_error(sprintf("%s Received %s.", supported, .msg_value(x)),
    arg = "x", received = .msg_value(x),
    expected = "an `effect_contrast_view` or a query-only `effect_view`")
}

# The frame normalization the view was read through, preferring the declared
# one: a diagonal metric folds into the weights and leaves `normalization =
# "none"` on a frame that was declared conservative (contract section 5.3).
.contribution_frame_normalization <- function(x) {
  metadata <- x$metadata
  if (!is.list(metadata) || !is.list(metadata$frame)) return(NULL)
  value <- metadata$frame$declared_normalization
  if (is.null(value)) value <- metadata$frame$normalization
  if (.is_string(value)) value else NULL
}

.contribution_require_conservative <- function(x) {
  normalization <- .contribution_frame_normalization(x)
  if (identical(normalization, "conservative")) return(normalization)
  if (identical(normalization, "local")) {
    .capability_refusal(paste0(
      "`contribution()` refuses a locally normalized frame. A local frame is ",
      "a detection map: every node reports the mean evidence density inside ",
      "its own support, overlapping neighbourhoods double-count the features ",
      "they share, and adding those values up estimates nothing ",
      "(`design/conservative-geometry-contract.md` section 1.1). Only a ",
      "column-normalized conservative frame partitions one fixed global ",
      "budget, which is what makes a territory sum a ledger entry."
    ),
      capability = "conservative_frame",
      namespace = "geometry_views",
      reasons = "local_frame_is_a_detection_map",
      remedies = paste0(
        "Rebuild the frame with `normalization = \"conservative\"` -- for ",
        "example `regions(labels, \"conservative\")` or ",
        "`searchlights(radius, \"conservative\")` -- rerun the view, and ",
        "aggregate that. Which instrument to use is a decision about the ",
        "question, not about scaling: a detection map answers *where is ",
        "there evidence*, an attribution map answers *how is the total ",
        "distributed*."
      )
    )
  }
  .capability_refusal(sprintf(paste0(
    "`contribution()` needs a conservative frame, and this view %s. Without ",
    "column normalization there is no fixed global budget for the group sums ",
    "to divide, so the aggregate would not be an attribution ledger."
  ), if (is.null(normalization)) {
    "does not record the normalization of the frame it was read through"
  } else {
    sprintf("was read through a frame declaring `normalization = \"%s\"`",
      normalization)
  }),
    capability = "conservative_frame",
    namespace = "geometry_views",
    reasons = if (is.null(normalization)) {
      "frame_normalization_not_recorded"
    } else {
      "frame_is_not_column_normalized"
    },
    remedies = paste0(
      "Rerun the view through a frame built with `normalization = ",
      "\"conservative\"`, and confirm the budget with ",
      "`frame_conservation(frame)`."
    )
  )
}

# Align an external per-measurement metadata table to the view's rows. A frame
# family's `$index` is the intended input, so a `measurement` column is joined
# on rather than assumed to be in row order.
.contribution_using_table <- function(using, ids) {
  if (is.null(using)) return(NULL)
  if (!is.data.frame(using) || !ncol(using) || !nrow(using)) {
    .input_error(sprintf(paste0(
      "`using` must be a nonempty data frame of per-measurement metadata; a ",
      "frame family's `$index` is the intended one. Received %s."
    ), .msg_value(using)),
      arg = "using", received = .msg_value(using),
      expected = "a per-measurement metadata data frame")
  }
  if ("measurement" %in% names(using)) {
    keys <- as.character(using$measurement)
    if (anyDuplicated(keys)) {
      repeated <- unique(keys[duplicated(keys)])
      .input_error(sprintf(paste0(
        "`using$measurement` repeats %s. Each measurement of the view must ",
        "match exactly one metadata row, or a node would be counted twice."
      ), .msg_names(repeated)),
        arg = "using", received = sprintf("repeated %s", .msg_names(repeated)),
        expected = "one row per measurement")
    }
    position <- match(ids, keys)
    if (anyNA(position)) {
      absent <- ids[is.na(position)]
      .input_error(sprintf(paste0(
        "`using` has no row for %s of the view's %s (%s). The join is by ",
        "`measurement`, so every measurement in the view must appear in the ",
        "table exactly once."
      ), .msg_count(length(absent), "measurement"),
        .msg_count(length(ids), "measurement"), .msg_names(absent)),
        arg = "using",
        received = sprintf("no row for %s", .msg_names(absent)),
        expected = "one row per measurement of the view")
    }
    return(using[position, , drop = FALSE])
  }
  if (nrow(using) != length(ids)) {
    .input_error(sprintf(paste0(
      "`using` has %s but the view has %s, and the table carries no ",
      "`measurement` column to join on. Supply one row per measurement in ",
      "`$index` order, or add a `measurement` column."
    ), .msg_count(nrow(using), "row"),
      .msg_count(length(ids), "measurement")),
      arg = "using", received = .msg_count(nrow(using), "row"),
      expected = .msg_count(length(ids), "row"))
  }
  using
}

# Resolve `by` to one label per measurement: a column of the metadata, or the
# grouping itself.
.contribution_grouping <- function(by, label, table, x, n) {
  if (.is_string(by)) {
    sources <- list()
    if (!is.null(table)) sources[["using"]] <- table
    if (is.data.frame(x$index)) sources[["the view's `$index`"]] <- x$index
    for (name in names(sources)) {
      candidate <- sources[[name]]
      if (by %in% names(candidate)) {
        return(list(label = by, values = candidate[[by]]))
      }
    }
    available <- unique(unlist(lapply(sources, names), use.names = FALSE))
    .input_error(sprintf(paste0(
      "`by = \"%s\"` names no column of the grouping metadata. Available ",
      "columns: %s. Pass the per-measurement table as `using =` -- a frame ",
      "family's `$index` carries `measurement`, `family`, `node`, `scale`, ",
      "`center` and `alpha` -- or pass the grouping itself as a vector with ",
      "one entry per measurement."
    ), by, if (length(available)) {
      .msg_names(available)
    } else {
      paste0("none, because this view's `$index` is a plain identifier ",
        "vector and no `using` table was supplied")
    }),
      arg = "by", received = sprintf("\"%s\"", by),
      expected = if (length(available)) {
        .msg_names(available)
      } else {
        "a grouping vector, or a `using` table"
      })
  }
  if ((!is.atomic(by) && !is.factor(by)) || length(by) != n) {
    .input_error(sprintf(paste0(
      "`by` must be one column name, or one group label per measurement. ",
      "The view has %s and `by` has %s."
    ), .msg_count(n, "measurement"), .msg_count(length(by), "entry", "entries")),
      arg = "by", received = .msg_count(length(by), "entry", "entries"),
      expected = .msg_count(n, "entry", "entries"))
  }
  list(label = label, values = by)
}

# Group labels to an ordered factor with no empty and no missing group. A row
# outside every group would silently take part of the budget with it, so an
# `NA` label is refused rather than dropped.
.contribution_groups <- function(values, label) {
  if (anyNA(values)) {
    .input_error(sprintf(paste0(
      "The grouping `%s` is missing for %s. Every measurement must belong to ",
      "exactly one group, or the aggregate would drop part of the budget it ",
      "claims to partition. (A frame family's `scale` is `NA` for a member ",
      "that has no scale, such as a point or region member; group by ",
      "`family` instead, or supply a grouping that covers every row.)"
    ), label, .msg_count(sum(is.na(values)), "measurement")),
      arg = "by",
      received = sprintf("%s missing of %s",
        .msg_count(sum(is.na(values)), "label"),
        .msg_count(length(values), "measurement")),
      expected = "one group label for every measurement")
  }
  if (is.factor(values)) {
    groups <- droplevels(values)
    return(list(groups = groups, keys = levels(groups)))
  }
  levels <- as.character(sort(unique(values)))
  groups <- factor(as.character(values), levels = levels)
  list(groups = groups, keys = values[match(levels, as.character(values))])
}

.contribution_group_sums <- function(values, rows) {
  vapply(rows, function(subset) sum(values[subset]), numeric(1),
    USE.NAMES = FALSE)
}

.contribution_group_column_sums <- function(values, rows) {
  out <- matrix(0, length(rows), ncol(values),
    dimnames = list(NULL, colnames(values)))
  for (position in seq_along(rows)) {
    out[position, ] <- colSums(values[rows[[position]], , drop = FALSE])
  }
  out
}

# One row per group. The key goes under `measurement`, because a group *is* the
# aggregate's measurement identity and every reader of a crossform result index
# -- `as.data.frame()`, the plot highlight resolver, a `merge()` back onto
# frame metadata -- looks for that column. A second column under the grouping's
# own name is added only when `as.character()` would lose the key's type, which
# is what grouping a frame family by `scale` does; for character or factor keys
# it would be a verbatim copy, and the grouping's name is on the print line and
# in `$metadata$aggregation$aggregated_by` either way.
.contribution_index <- function(label, keys, rows) {
  index <- data.frame(measurement = as.character(keys),
    stringsAsFactors = FALSE, check.names = FALSE)
  typed <- !is.character(keys) && !is.factor(keys)
  if (typed && !identical(label, "measurement")) index[[label]] <- keys
  index[["n_rows"]] <- as.integer(lengths(rows, use.names = FALSE))
  index
}

.contribution_scientific_id <- function(parent, label, keys) {
  .sha256_signature(list(
    schema_version = 1L,
    role = "geometry_contribution",
    parent = parent,
    aggregated_by = label,
    groups = as.character(keys)
  ), "geometry-sha256:")
}

.contribution_provenance <- function(label, keys, rows, normalization,
                                     budget_exact, frame_relative_components,
                                     masked = character()) {
  list(
    aggregated_by = label,
    groups = length(rows),
    measurements = sum(lengths(rows)),
    frame_normalization = normalization,
    frame_relative = length(frame_relative_components) > 0L,
    budget_exact = budget_exact,
    frame_relative_components = frame_relative_components,
    masked = masked,
    overlap_split = FALSE,
    group_keys = as.character(keys)
  )
}

# `$signed` is the local weighted mean contrast, not a contribution: it is
# already divided by the node's own frame mass. Summing means over a territory
# is the error `contribution()` exists to prevent, and the mass-weighted mean
# that would be defensible needs masses no view carries. So the shape is kept
# and the values are masked, which is what the package does everywhere else a
# number is not earned.
.contribution_masked_signed <- function(signed, groups) {
  if (is.list(signed)) {
    return(lapply(signed, function(value) rep(NA_real_, groups)))
  }
  rep(NA_real_, groups)
}

#' Add a conservative attribution map up over a territory
#'
#' A conservative frame partitions one fixed global budget: each node's value
#' is its share, not its density, so the shares of a territory add up
#' (`design/conservative-geometry-contract.md` sections 1.2 and 2).
#' `contribution()` performs exactly that addition, and nothing else. It is the
#' ledger reading of an attribution map, and it refuses a detection map, whose
#' overlapping nodes double-count and whose sum estimates nothing.
#'
#' Grouping is **by row**: each row of `x` belongs to exactly one group, so the
#' group totals partition the whole-domain total exactly. Splitting an
#' overlapping node's mass across several territories is deliberately not
#' offered -- a row is one node, a node is not divisible, and any split needs a
#' second partition of the frame weights that can double-count.
#'
#' @param x A view whose per-measurement values are additive contributions: an
#'   `effect_contrast_view` from [contrast_energy()], or a query-only
#'   `effect_view` from [evaluate_geometry()] or [query_geometry()]. Every
#'   other result kind is refused with the reason; see *Refusals*.
#' @param by The grouping. Either one column name of the per-measurement
#'   metadata (searched in `using` first, then in `x$index` when that is a
#'   table), or the grouping itself as a vector or factor with one entry per
#'   measurement -- a region or network label, say. A length-one character
#'   value is always read as a column name. Every measurement must carry a
#'   group; an `NA` label is refused rather than dropped.
#' @param using Optional per-measurement metadata table naming the groups a
#'   view cannot carry itself. A [frame_family()]'s `$index` is the intended
#'   one: it is joined on its `measurement` column, so row order does not
#'   matter. A table without that column must have one row per measurement, in
#'   `$index` order.
#' @return The same class as `x`, with one row per group instead of one per
#'   measurement.
#'
#'   `$index` has one row per group: `measurement`, holding the group key, and
#'   `n_rows`, the number of measurements that went into it. A key that
#'   `as.character()` would flatten -- a numeric `scale`, say -- also gets a
#'   typed column under the grouping's own name. Groups appear in sorted key
#'   order, or in level order when `by` is a factor; an empty group is dropped.
#'   `$metadata$aggregation` records `aggregated_by`, the group count,
#'   `frame_relative`, and which components are budget-exact. `$receipt` is the
#'   parent receipt under a derived `scientific_plan_id`.
#' @section What conserves and what does not:
#' `total` is budget-exact: the group totals add back to the whole-domain
#' total, which is the value the same contrast takes under
#' `whole_brain("none")` (contract section 2, claim 2).
#'
#' `signed` is **not**, and the aggregate reports `NA` for it. A contrast
#' view's signed marginal is the local weighted *mean* contrast -- already
#' divided by the node's own frame mass -- so it is a density, and adding
#' densities over a territory is the error a conservative frame exists to
#' avoid (contract section 1.1). The defensible aggregate is a mass-weighted
#' mean, which needs the node masses a view does not carry, so the field is
#' masked rather than filled with a number nobody can interpret.
#'
#' `coherent` and `configuration` add up as arithmetic, but their sums are
#' **frame-relative** (contract section 4, claim 4): \eqn{\sum_x G_x^{coh}} is
#' not a global quantity, so a coherent budget is a share of *this frame's*
#' coherent mass and two frames give two incomparable denominators. The result
#' carries `frame_relative = TRUE` and the print says so. Never compare a
#' coherent budget across frames.
#'
#' `coherence_fraction` is recomputed from the aggregated components, not
#' averaged: a fraction of sums is not a sum of fractions. It is masked by the
#' same rule a node's fraction is masked by -- reported only where the
#' aggregated total is finite and positive and both components are nonnegative
#' -- so a group whose components do not form a nonnegative partition reports
#' `NA` rather than a clamped number.
#' @section Refusals:
#' A locally normalized frame signals an `effect_capability_refusal` with
#' capability `"conservative_frame"` in namespace `"geometry_views"`; so does a
#' view that does not record its frame normalization. A view whose values are
#' not additive contributions -- an `effect_spectrum_view`, `effect_rdm_view`,
#' `effect_rsa_view`, `effect_crossnobis_view` -- a stored `effect_geometry`,
#' and an unevaluated `effect_geometry_plan` each signal capability
#' `"additive_contribution"`. Branch on them with [catch_refusal()].
#' @seealso [contrast_energy()] and [evaluate_geometry()] for the views this
#'   aggregates, [frame_family()] for the `$index` that names scales and
#'   members, and [frame_conservation()] for the certificate that the frame
#'   conserves the budget being divided.
#' @family geometry plans and views
#' @examples
#' # Four regions over eight voxels, read as an attribution map.
#' domain <- abstract_domain(8, id = "contribution-example")
#' run1 <- rbind(
#'   face = c(1, 0.2, 0, 0, 0.5, 0.1, 0, 0),
#'   house = c(0, 1, 0.1, 0, 0, 0.6, 0.2, 0)
#' )
#' run2 <- rbind(
#'   face = c(0.9, 0.3, 0, 0, 0.4, 0.2, 0, 0),
#'   house = c(0.1, 0.9, 0, 0, 0.1, 0.5, 0.3, 0)
#' )
#' relation <- relation(list(run1 = run1, run2 = run2), domain = domain)
#' labels <- rep(c("r1", "r2", "r3", "r4"), each = 2)
#' plan <- plan_geometry(
#'   relation, compile_frame(regions(labels, "conservative"), domain),
#'   cross_partitions(relation, independence = "independent")
#' )
#' effect <- contrast_energy(plan, c(face = 1, house = -1))
#'
#' # Two networks, each two regions. The ledger adds up exactly.
#' network <- c("early", "early", "late", "late")
#' ledger <- contribution(effect, by = network)
#' ledger
#' sum(ledger$total) - sum(effect$total)
#'
#' # A detection map is refused: overlapping local nodes double-count, so
#' # their sum is not a share of anything.
#' detection <- contrast_energy(
#'   plan_geometry(
#'     relation, compile_frame(regions(labels), domain),
#'     cross_partitions(relation, independence = "independent")
#'   ),
#'   c(face = 1, house = -1)
#' )
#' refusal <- catch_refusal(contribution(detection, by = network))
#' refusal$capability
#' @param ... Passed to the method. The generic dispatches on `x`: a
#'   per-measurement view is aggregated by the default method documented here,
#'   and an `effect_population_result` or `effect_population_view` by the
#'   group-node method in [population_views].
#' @export
contribution <- function(x, ...) UseMethod("contribution")

#' @rdname contribution
#' @export
contribution.default <- function(x, by, using = NULL, ...) {
  .check_no_extra_arguments("contribution", ...)
  source <- .contribution_source(x)
  if (missing(by)) {
    .input_error(paste0(
      "`by` is required: name the grouping to add the contributions up over, ",
      "either as a column of the per-measurement metadata (for example ",
      "`contribution(view, by = \"family\", using = family$index)`) or as ",
      "one label per measurement (`contribution(view, by = network)`)."
    ),
      arg = "by", received = "no argument",
      expected = "a column name, or one group label per measurement")
  }
  label <- substitute(by)
  label <- if (is.symbol(label)) as.character(label) else "group"
  if (!is.null(using) && !.is_string(by)) {
    # Checked before the table is aligned, so the message names the real
    # mistake rather than reporting a join that was never going to be used.
    .input_error(paste0(
      "`using` is only meaningful when `by` names one of its columns. `by` ",
      "was given as a grouping vector, which already says which group every ",
      "measurement belongs to."
    ),
      arg = "using", received = "a table alongside a grouping vector",
      expected = "`by` as a column name, or `using = NULL`")
  }
  normalization <- .contribution_require_conservative(x)
  ids <- .contribution_measurement_ids(x$index)
  table <- .contribution_using_table(using, ids)
  resolved <- .contribution_grouping(by, label, table, x, source$measurements)
  grouped <- .contribution_groups(resolved$values, resolved$label)
  rows <- split(seq_len(source$measurements), grouped$groups)
  index <- .contribution_index(resolved$label, grouped$keys, rows)
  receipt <- .projection_receipt(
    x$receipt,
    .contribution_scientific_id(x$receipt$scientific_plan_id,
      resolved$label, grouped$keys)
  )
  metadata <- if (is.list(x$metadata)) x$metadata else list()
  metadata$scientific_plan_id <- NULL
  metadata$aggregated_from <- x$receipt$scientific_plan_id

  if (identical(source$kind, "contrast")) {
    metadata$aggregation <- .contribution_provenance(
      resolved$label, grouped$keys, rows, normalization,
      budget_exact = "total",
      frame_relative_components = c("coherent", "configuration"),
      masked = "signed"
    )
    signed <- .contribution_masked_signed(x$signed, length(rows))
    coherent <- .contribution_group_sums(x$coherent, rows)
    total <- .contribution_group_sums(x$total, rows)
    return(.effect_contrast_view_record(
      signed, coherent, total - coherent, total, x$weights, index, receipt,
      metadata
    ))
  }

  metadata$aggregation <- .contribution_provenance(
    resolved$label, grouped$keys, rows, normalization,
    budget_exact = if (identical(x$component, "total")) "values" else character(),
    frame_relative_components = if (identical(x$component, "total")) {
      character()
    } else {
      x$component
    }
  )
  effect_view(
    values = .contribution_group_column_sums(x$values, rows),
    query = x$query,
    component = x$component,
    receipt = receipt,
    index = index,
    metadata = metadata,
    left_space = x$left_space,
    right_space = x$right_space
  )
}

# `coherence_spectrum()` -- coherent share versus scale (and location) --------
#
# The reduction `design/conservative-geometry-contract.md` section 3.2 calls
# the scientific product of a conservative frame family. Its whole reason to
# exist is one asymmetry between the two columns it reports:
#
#   * Per-scale ENERGY is fixed by the family weighting. Section 3.1: every
#     member is column normalized on its own and the weights sum to one, so
#     the rows of scale `s` sum to `alpha_s * G_Omega` whatever the data say.
#     A panel of energy against scale is a plot of the analyst's own `alpha`
#     vector, and the contract makes it normative that no such panel may be
#     presented as evidence about spatial scale.
#   * The coherent SHARE is exactly invariant to it. Both components are
#     homogeneous of degree one under a row rescaling `w_x -> alpha * w_x`:
#     the total is linear in `w`, and the coherent part is
#     `K_coh = a a^T / (a^T K_x^-1 a)` with `a = w_x / sum(w_x)` invariant to
#     alpha while `K_x^-1` scales as `1/alpha`, so `K_coh` scales as alpha.
#     The ratio cancels it exactly (identity-metric form `(Bw)(Bw)^T/sum(w)`:
#     manifestly degree one). So the spectrum is a property of the family's
#     geometry rather than of the weights, and may be reported without
#     disclosing them -- the exact opposite of the energy panel.
#
# Three implementation decisions follow, and none of them is free:
#
#  1. The reduction is `contribution()` under a different grouping, so it is
#     built out of the same machinery: the conservative-frame refusal, the
#     metadata join, the row sums, the masked signed marginal, and above all
#     `.coherence_fraction()`. The share a scale reports is masked by exactly
#     the rule a node's fraction is masked by. What is new is the grouping --
#     `contribution()` takes one column, and the spectrum lives on the PAIR
#     (center, scale) whenever it is read location-wise (section 11.4, gap
#     G3) -- and the provenance saying which column is the finding.
#  2. `total` is summed and `configuration` is taken as `total - coherent`,
#     rather than summing `configuration` and adding. That keeps the
#     budget-exact column exact: `sum` of the scale's totals is `alpha_s`
#     times the whole-domain total to the tolerance section 9 states, and the
#     share is identically `sum_coh / (sum_coh + sum_cfg)` either way.
#  3. There is no location-wise collapse. Gap G3 says the coherent share is a
#     function of (location, scale) and not a number, so `by_location = TRUE`
#     returns that table and stops. A scale profile at one location is a
#     filter on it, which needs no named reduction; an alpha-weighted mean or
#     an argmax-scale over it would be one, and would have to arrive as its
#     own declared operation with its own certificate.

# Group keys as printable strings, one value at a time: a scale of 4 has to
# print as "4" and not "4.000000", and formatting the vector as a whole would
# pad every entry to the width of the longest.
.coherence_spectrum_key_strings <- function(values) {
  if (is.character(values)) return(values)
  if (is.factor(values)) return(as.character(values))
  if (is.numeric(values)) {
    return(vapply(values, function(value) {
      format(value, trim = TRUE, digits = 15)
    }, character(1), USE.NAMES = FALSE))
  }
  as.character(values)
}

# `by` and `by_location` are two spellings of one choice, so supplying both is
# an error rather than a precedence rule nobody can remember.
.coherence_spectrum_keys <- function(by, by_location) {
  if (!is.null(by)) {
    if (!isFALSE(by_location)) {
      .input_error(paste0(
        "`by` and `by_location` name the same choice twice. `by_location = ",
        "TRUE` is the shorthand for `by = c(\"scale\", \"center\")`; pass ",
        "one of them."
      ),
        arg = "by_location", received = "both `by` and `by_location`",
        expected = "`by`, or `by_location`, but not both")
    }
    if (!is.character(by) || !length(by) || anyNA(by) || !all(nzchar(by))) {
      .input_error(sprintf(paste0(
        "`by` must name one or more columns of the per-measurement metadata, ",
        "for example `by = \"scale\"` or `by = c(\"scale\", \"center\")`; ",
        "received %s."
      ), .msg_value(by)),
        arg = "by", received = .msg_value(by),
        expected = "one or more metadata column names")
    }
    if (anyDuplicated(by)) {
      .input_error(sprintf(paste0(
        "`by` repeats %s. Each grouping column may appear once."
      ), .msg_names(unique(by[duplicated(by)]))),
        arg = "by", received = sprintf("repeated %s",
          .msg_names(unique(by[duplicated(by)]))),
        expected = "distinct column names")
    }
    return(as.character(by))
  }
  if (!is.logical(by_location) || length(by_location) != 1L ||
      is.na(by_location)) {
    .input_error(sprintf(paste0(
      "`by_location` must be `TRUE` or `FALSE`; received %s. `TRUE` returns ",
      "the (scale, center) table, `FALSE` the per-scale spectrum."
    ), .msg_value(by_location)),
      arg = "by_location", received = .msg_value(by_location),
      expected = "TRUE or FALSE")
  }
  if (by_location) c("scale", "center") else "scale"
}

# What `coherence_spectrum()` reads, and the reasoned refusal for everything
# else. The share needs both components of one measurement at once, which is
# what rules out every readout that carries a single number per row.
.coherence_spectrum_needs <- function() {
  paste0(
    "`coherence_spectrum()` reads the coherent/configuration split of one ",
    "contrast over a frame family: pass a compiled `effect_geometry_plan` ",
    "with `weights`, or the `effect_contrast_view` that `contrast_energy()` ",
    "already returned for it."
  )
}

.coherence_spectrum_unsupported <- function(x) {
  needs <- .coherence_spectrum_needs()
  if (inherits(x, "effect_view") && !inherits(x, "effect_form")) {
    .capability_refusal(sprintf(paste0(
      "`coherence_spectrum()` cannot read a query-only `effect_view`. It ",
      "holds one component -- this one is `%s` -- and a coherent share is a ",
      "ratio of two, so there is nothing here to take the share of. %s"
    ), if (.is_string(x$component)) x$component else "one", needs),
      capability = "coherence_decomposition",
      namespace = "geometry_views",
      reasons = "single_component_view_has_no_share",
      remedies = paste0(
        "Read the same plan with `contrast_energy(plan, weights)`, which ",
        "returns `coherent`, `configuration` and `total` together, and pass ",
        "that view."
      )
    )
  }
  comparative <- c("effect_spectrum_view", "effect_rdm_view",
    "effect_rsa_view", "effect_crossnobis_view")
  if (any(comparative %in% class(x))) {
    kind <- comparative[[which(comparative %in% class(x))[[1L]]]]
    .capability_refusal(sprintf(paste0(
      "`coherence_spectrum()` does not read an `%s`. Its rows are ",
      "eigenvalues, dissimilarities or regression coefficients -- none of ",
      "them carries the coherent/configuration split a share is taken of. %s"
    ), kind, needs),
      capability = "coherence_decomposition",
      namespace = "geometry_views",
      reasons = "readout_carries_no_coherence_decomposition",
      remedies = paste0(
        "Run the same contrast through `contrast_energy(plan, weights)` and ",
        "pass that view, or pass the plan itself with `weights`."
      )
    )
  }
  if (inherits(x, "effect_form") || inherits(x, "effect_geometry")) {
    .capability_refusal(paste0(
      "`coherence_spectrum()` reads a contrast, not packed geometry: a ",
      "stored form holds every effect pair at once and names no single ",
      "contrast to take the share of. ", needs
    ),
      capability = "coherence_decomposition",
      namespace = "geometry_views",
      reasons = "packed_geometry_names_no_contrast",
      remedies = paste0(
        "Read one contrast first with `contrast_energy(x, weights)`, then ",
        "pass that view."
      )
    )
  }
  .input_error(sprintf("%s Received %s.", needs, .msg_value(x)),
    arg = "x", received = .msg_value(x),
    expected = "an `effect_geometry_plan` or an `effect_contrast_view`")
}

# A plan is evaluated here rather than by the caller, so that the frame's own
# `$index` is available as the grouping metadata: a family's per-row `family`
# / `scale` / `center` / `alpha` reaches a compiled frame but not a result,
# whose `$index` is the measurement identifier vector alone.
.coherence_spectrum_source <- function(x, weights, using) {
  if (inherits(x, "effect_geometry_plan")) {
    if (is.null(weights)) {
      .input_error(paste0(
        "`weights` is required with a plan: a coherence spectrum is the ",
        "share of one contrast, so name it, for example ",
        "`coherence_spectrum(plan, c(face = 1, house = -1))`. Unnamed ",
        "weights are taken in the relation's declared effect order."
      ),
        arg = "weights", received = "no argument",
        expected = "one finite weight per experimental effect")
    }
    .validate_geometry_plan(x)
    return(list(
      view = contrast_energy(x, weights),
      metadata = if (is.null(using)) x$frame$index else using,
      source = "plan"
    ))
  }
  if (inherits(x, "effect_contrast_view")) {
    if (!is.null(weights)) {
      .input_error(paste0(
        "`weights` applies only to a plan. This view was already read at one ",
        "contrast -- `x$weights` is it -- and a second contrast cannot be ",
        "applied to numbers that have been reduced. Pass the plan with the ",
        "weights you mean, or drop `weights`."
      ),
        arg = "weights", received = "weights alongside an evaluated view",
        expected = "`weights` with a plan, or no `weights` with a view")
    }
    metadata <- if (!is.null(using)) {
      using
    } else if (is.data.frame(x$index)) {
      x$index
    } else {
      NULL
    }
    return(list(view = x, metadata = metadata, source = "view"))
  }
  .coherence_spectrum_unsupported(x)
}

# The grouping columns, with the message a missing one has to carry: a view's
# `$index` is a bare identifier vector, so the family's `$index` almost always
# has to arrive as `using =`.
.coherence_spectrum_require_keys <- function(metadata, keys, source) {
  available <- if (is.data.frame(metadata)) names(metadata) else character()
  absent <- setdiff(keys, available)
  if (!length(absent)) return(invisible(NULL))
  .input_error(sprintf(paste0(
    "%s of the grouping metadata. Available columns: %s. A frame family's ",
    "`$index` is the intended table -- it carries `measurement`, `family`, ",
    "`node`, `scale`, `center` and `alpha` -- so build the frame with ",
    "`searchlights(c(...), \"conservative\")` or `frame_family()`, and %s"
  ), if (length(absent) == 1L) {
    sprintf("`by = \"%s\"` names no column", absent)
  } else {
    sprintf("%s name no column", .msg_names(absent))
  }, if (length(available)) {
    .msg_names(available)
  } else {
    "none, because no per-measurement metadata table was available"
  }, if (identical(source, "plan")) {
    "compile the plan on it, or pass the table as `using =`."
  } else {
    "pass its `$index` as `using =`."
  }),
    arg = if (length(available)) "by" else "using",
    received = if (length(available)) {
      .msg_names(absent)
    } else {
      "no per-measurement metadata"
    },
    expected = if (length(available)) {
      .msg_names(available)
    } else {
      "a frame family's `$index` as `using =`"
    })
}

# One group per distinct key tuple, ordered by the key columns in the order
# they were named. An `NA` key is refused rather than dropped, for the reason
# `contribution()` refuses one: a row outside every group takes part of the
# budget with it.
.coherence_spectrum_groups <- function(metadata, keys) {
  columns <- lapply(keys, function(key) metadata[[key]])
  for (position in seq_along(keys)) {
    values <- columns[[position]]
    if (!anyNA(values)) next
    .input_error(sprintf(paste0(
      "The grouping `%s` is missing for %s. Every measurement must belong to ",
      "exactly one group, or the spectrum would drop part of the budget it ",
      "splits. A frame family's `scale` is `NA` for a member that has no ",
      "scale, and its `center` is `NA` for a member whose rows are not ",
      "anchored at a feature -- a region or whole-brain member is both -- so ",
      "group by `family` instead, or build the family from members that all ",
      "carry the column you asked for."
    ), keys[[position]], .msg_count(sum(is.na(values)), "measurement")),
      arg = "by",
      received = sprintf("%s missing of %s",
        .msg_count(sum(is.na(values)), "label"),
        .msg_count(length(values), "measurement")),
      expected = "one group label for every measurement")
  }
  labels <- do.call(paste, c(
    lapply(columns, .coherence_spectrum_key_strings), list(sep = "::")
  ))
  ordering <- do.call(order, unname(columns))
  levels <- unique(labels[ordering])
  list(groups = factor(labels, levels = levels), keys = levels)
}

# One row per group: the composite key under `measurement` (which is what
# every reader of a crossform index looks for), the grouping columns
# themselves with their own types kept, and the row count.
#
# Any other metadata column that takes ONE value inside every group is carried
# through as well, which is how `family` and `alpha` reach the table without
# being named: the family weight belongs next to the energy it fixes. Row
# identities are excluded by name, because a group of one row would otherwise
# carry a second `measurement`-like column that means something else.
.coherence_spectrum_index <- function(metadata, keys, rows, group_keys) {
  first <- vapply(rows, function(subset) subset[[1L]], integer(1),
    USE.NAMES = FALSE)
  index <- data.frame(measurement = as.character(group_keys),
    stringsAsFactors = FALSE, check.names = FALSE)
  for (key in keys) index[[key]] <- metadata[[key]][first]
  carried <- setdiff(names(metadata), c(keys, "measurement", "node"))
  for (column in carried) {
    values <- metadata[[column]]
    constant <- all(vapply(rows, function(subset) {
      length(unique(values[subset])) == 1L
    }, logical(1), USE.NAMES = FALSE))
    if (constant) index[[column]] <- values[first]
  }
  index[["n_rows"]] <- as.integer(lengths(rows, use.names = FALSE))
  index
}

.coherence_spectrum_scientific_id <- function(parent, keys, group_keys) {
  .sha256_signature(list(
    schema_version = 1L,
    role = "geometry_coherence_spectrum",
    parent = parent,
    resolved_by = as.character(keys),
    groups = as.character(group_keys)
  ), "geometry-sha256:")
}

# `contribution()`'s provenance record, plus the three facts that are specific
# to this reduction: which grouping produced it, that the share and not the
# energy is the finding, and the family weights the energy column is fixed by.
.coherence_spectrum_provenance <- function(keys, group_keys, rows,
                                           normalization, index) {
  record <- .contribution_provenance(
    paste(keys, collapse = " x "), group_keys, rows, normalization,
    budget_exact = "total",
    frame_relative_components = c("coherent", "configuration"),
    masked = "signed"
  )
  record$reduction <- "coherence_spectrum"
  record$resolved_by <- as.character(keys)
  record$by_location <- identical(as.character(keys), c("scale", "center"))
  record$alpha_invariant <- "coherence_fraction"
  record$alpha_fixed <- "total"
  record$location_collapse <- "none"
  if ("alpha" %in% names(index)) {
    record$alpha <- stats::setNames(as.numeric(index$alpha),
      as.character(index$measurement))
  }
  if ("family" %in% names(index)) {
    record$family <- stats::setNames(as.character(index$family),
      as.character(index$measurement))
  }
  record
}

#' Read the coherent share of a conservative frame family against scale
#'
#' A multiscale conservative frame family divides one fixed global budget
#' between its scales, and it divides it by the weights the analyst chose:
#' the `total` component summed over the rows of scale \eqn{s} is exactly
#' \eqn{\alpha_s G_\Omega}{alpha_s * G_Omega} whatever the data say
#' (`design/conservative-geometry-contract.md` section 3.1). What the data do
#' decide is how each scale's fixed budget splits into a coherent and a
#' configuration part. `coherence_spectrum()` reports that split, per scale
#' and -- with `by_location = TRUE` -- per (scale, location).
#'
#' The reported share is
#' \deqn{\phi_s = \frac{\sum_{x \in s} \langle H, G^{\mathrm{coh}}_x\rangle}
#'                     {\sum_{x \in s} \langle H, G_x\rangle},}{phi_s = sum_coherent / sum_total,}
#' computed from the **aggregated** components. A fraction of sums is not a
#' sum of fractions, so averaging per-node fractions would answer a different
#' question.
#'
#' @param x A compiled `effect_geometry_plan` over a frame family, or the
#'   `effect_contrast_view` that [contrast_energy()] already returned for one.
#'   Every other result kind is refused with the reason; see *Refusals*.
#' @param weights One finite contrast weight per experimental effect, as in
#'   [contrast_energy()]. Required with a plan, and refused with an evaluated
#'   view, which was already read at one contrast.
#' @param by_location `FALSE` (the default) returns one row per scale;
#'   `TRUE` returns one row per (scale, center), the object gap G3 of
#'   `design/conservative-geometry-contract.md` section 11.4 requires. It is
#'   the shorthand for `by = c("scale", "center")`.
#' @param by The grouping columns, when they are not the two above: one or
#'   more column names of the per-measurement metadata, for example
#'   `by = "family"` for a family whose members have no numeric scale. Pass
#'   `by` or `by_location`, not both.
#' @param using Optional per-measurement metadata table, joined on its
#'   `measurement` column so row order does not matter. A [frame_family()]'s
#'   `$index` is the intended one, and it is required on the view route: a
#'   result's `$index` is the measurement identifier vector alone, so an
#'   evaluated view does not carry the scales it was read at. With a plan the
#'   frame's own `$index` is used unless `using` overrides it.
#' @return An `effect_contrast_view` with one row per group instead of one per
#'   measurement -- the same record `contribution()` returns, because the
#'   spectrum is that aggregation under a scale-resolved grouping.
#'
#'   `$coherent` is \eqn{E^{\mathrm{coh}}_s}, `$configuration` is
#'   \eqn{E^{\mathrm{cfg}}_s}, `$total` is their sum, and
#'   `$coherence_fraction` is \eqn{\phi_s}, masked by
#'   `$coherence_fraction_valid`. `$signed` is `NA`: a signed marginal is a
#'   local weighted mean, and means do not add over a territory.
#'
#'   `$index` has one row per group: `measurement` holds the composite key
#'   (`"1.01::v8"` for a (scale, center) row), the grouping columns keep their
#'   own types, `n_rows` counts the measurements behind the row, and any other
#'   metadata column that takes one value inside every group is carried
#'   through -- which is how `family` and `alpha` arrive next to the energy
#'   they fix. `$metadata$aggregation` records `reduction =
#'   "coherence_spectrum"`, `resolved_by`, `frame_relative = TRUE`,
#'   `alpha_invariant = "coherence_fraction"`, `alpha_fixed = "total"`,
#'   `location_collapse = "none"`, and the applied `alpha` per group.
#'   `$receipt` is the parent receipt under a derived `scientific_plan_id`.
#' @section Read the share, not the energy:
#' Every member of a conservative family is column normalized on its own and
#' the family weights sum to one, so each scale's rows carry exactly their
#' weight of the whole-domain total. **The `total` column is therefore the
#' `weights` vector times a constant, and a plot of energy against scale is a
#' plot of that vector.** The contract makes this normative: no such panel may
#' be presented as evidence about spatial scale.
#'
#' The share is the opposite. Both components are homogeneous of degree one
#' under a rescaling \eqn{w_x \mapsto \alpha w_x} of a row -- the total is
#' linear in \eqn{w}, and the coherent part
#' \eqn{K_{\mathrm{coh}} = a a^\top / (a^\top K_x^{-1} a)} has
#' \eqn{a = w_x / \sum w_x} invariant to \eqn{\alpha} while \eqn{K_x^{-1}}
#' scales as \eqn{1/\alpha} -- so the ratio cancels \eqn{\alpha} exactly. Two
#' families differing only in their weights give identical shares to machine
#' precision and different energies. The spectrum is a property of the
#' family's geometry, not of the weighting, and may be reported without
#' disclosing it (contract sections 3.1 and 3.2).
#' @section What is frame-relative, and what is masked:
#' `total` is budget-exact: the group totals add back to the whole-domain
#' total (contract section 2). `coherent` and `configuration` add up as
#' arithmetic, but \eqn{\sum_x G^{\mathrm{coh}}_x} is not a global quantity
#' (section 4, claim 4), so a coherent energy is a share of *this family's*
#' coherent mass and two families give two incomparable denominators. The
#' result carries `frame_relative = TRUE` and the print says so.
#'
#' `coherence_fraction` is masked, never clamped: it is reported only where
#' the aggregated total is finite and positive and both aggregated components
#' are nonnegative, and is `NA` elsewhere. Cross-generalized components are
#' signed, so a scale whose nodes' common modes anticorrelate across
#' partitions can carry a negative coherent energy; that scale reports `NA`
#' rather than a number outside \eqn{[0,1]}.
#'
#' One consequence is worth expecting. A singleton scale -- a point member, or
#' a radius covering one feature -- has exactly zero configuration
#' (`effect-form-v1` section 7), so its aggregated configuration sits within
#' one unit in the last place of zero and can land on either side of it. Where
#' it lands below, the mask fires and the share reads `NA` instead of `1`.
#' That is the mask working on an exactly-degenerate partition, and it is the
#' same behaviour a singleton node's own `coherence_fraction` already has.
#' @section No location-wise collapse:
#' The coherent share is a function of (location, scale) and not a number
#' (gap G3). `by_location = TRUE` returns that table and stops: the scale
#' profile at one location is a filter on it, which needs no named operation,
#' while an alpha-weighted mean over scales or an argmax-scale map is a
#' declared reduction that would need its own certificate.
#' `$metadata$aggregation$location_collapse` records `"none"`.
#' @section Refusals:
#' A locally normalized frame signals an `effect_capability_refusal` with
#' capability `"conservative_frame"` in namespace `"geometry_views"`; so does
#' a view that does not record its frame normalization. A query-only
#' `effect_view`, a comparative readout (`effect_spectrum_view`,
#' `effect_rdm_view`, `effect_rsa_view`, `effect_crossnobis_view`), and
#' packed geometry each signal capability `"coherence_decomposition"`, because
#' none of them carries the two components a share is taken of. Branch on them
#' with [catch_refusal()].
#' @seealso [searchlights()] and [frame_family()] for the families this reads,
#'   [contrast_energy()] for the view it aggregates, [contribution()] for the
#'   same arithmetic over a spatial territory, and [frame_conservation()] for
#'   the per-block certificate that each scale carries its own weight.
#' @family geometry plans and views
#' @examples
#' # A point effect on a line: one voxel carries the whole contrast.
#' n <- 15L
#' domain <- abstract_domain(
#'   n, coordinates = cbind(seq_len(n) - 1, 0),
#'   feature_ids = paste0("v", seq_len(n)), id = "spectrum-example"
#' )
#' pattern <- function(signal) {
#'   rbind(face = signal, house = rep(0, n))
#' }
#' signal <- replace(rep(0, n), 8L, 1)
#' relation <- relation(
#'   list(run1 = pattern(signal), run2 = pattern(signal)), domain = domain
#' )
#' family <- compile_frame(
#'   searchlights(c(0.5, 1.01, 2.01), "conservative"), domain
#' )
#' plan <- plan_geometry(
#'   relation, family,
#'   cross_partitions(relation, independence = "independent")
#' )
#'
#' # The share falls as the neighborhood grows past the signal: one voxel of
#' # evidence looks entirely coherent at a scale that sees only it.
#' spectrum <- coherence_spectrum(plan, c(face = 1, house = -1))
#' as.data.frame(spectrum)[, c("scale", "alpha", "total", "coherence_fraction")]
#'
#' # The energy column is the family weighting and nothing else, so it is the
#' # same at every scale here; the share is what varies.
#' spectrum$total
#'
#' # One row per (scale, center) instead. A location's scale profile is a
#' # filter on that table, not a separate reduction.
#' located <- coherence_spectrum(plan, c(face = 1, house = -1),
#'   by_location = TRUE)
#' profile <- as.data.frame(located)
#' profile[profile$center == "v8", c("scale", "center", "coherence_fraction")]
#' @export
coherence_spectrum <- function(x, weights = NULL, by_location = FALSE,
                               by = NULL, using = NULL) {
  keys <- .coherence_spectrum_keys(by, by_location)
  source <- .coherence_spectrum_source(x, weights, using)
  view <- source$view
  normalization <- .contribution_require_conservative(view)
  ids <- .contribution_measurement_ids(view$index)
  # A frame that is not a family carries an `$index` of measurement labels and
  # nothing else, and a declared frame may carry none at all. Neither is a
  # `using` mistake, so both fall through to the message that names the family
  # constructors rather than to the join's.
  supplied <- if (is.data.frame(source$metadata)) source$metadata else NULL
  metadata <- .contribution_using_table(supplied, ids)
  .coherence_spectrum_require_keys(metadata, keys, source$source)

  grouped <- .coherence_spectrum_groups(metadata, keys)
  rows <- split(seq_along(view$total), grouped$groups)
  index <- .coherence_spectrum_index(metadata, keys, rows, grouped$keys)

  coherent <- .contribution_group_sums(view$coherent, rows)
  total <- .contribution_group_sums(view$total, rows)
  signed <- .contribution_masked_signed(view$signed, length(rows))

  record <- if (is.list(view$metadata)) view$metadata else list()
  record$scientific_plan_id <- NULL
  record$aggregated_from <- view$receipt$scientific_plan_id
  record$aggregation <- .coherence_spectrum_provenance(
    keys, grouped$keys, rows, normalization, index
  )
  receipt <- .projection_receipt(
    view$receipt,
    .coherence_spectrum_scientific_id(
      view$receipt$scientific_plan_id, keys, grouped$keys
    )
  )
  .effect_contrast_view_record(signed, coherent, total - coherent, total,
    view$weights, index, receipt, record)
}
