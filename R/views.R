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
#' @export
contrast_energy <- function(x, weights, remove_univariate = FALSE) {
  if (missing(weights)) {
    .input_error(paste0(
      "`weights` is required: pass one finite weight per experimental ",
      "effect, for example `contrast_energy(plan, c(face = 1, house = -1))`. ",
      "Unnamed weights are taken in the relation's declared effect order."
    ))
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
    ), .msg_value(x)))
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
    total, coherent, x$marginals, weights, x$index, receipt
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
    ), operation, .msg_value(x)))
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
#' \deqn{d_{ij} = c^T G c = G_{ii} + G_{jj} - 2G_{ij}.}
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
#' @export
rdm <- function(x, component = c("total", "coherent", "configuration"),
                pairs = NULL, normalize = NULL) {
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
    ), label, .msg_value(models)))
  }
  if (is.null(names(models)) || anyNA(names(models)) ||
      any(!nzchar(names(models))) || anyDuplicated(names(models))) {
    .input_error(sprintf(paste0(
      "`%s` must be a list with unique nonempty names; the names become the ",
      "coefficient columns of the fit."
    ), label))
  }
  entries <- names(models)
  out <- lapply(entries, function(entry) {
    value <- models[[entry]]
    where <- sprintf("`%s` RDM `%s`", label, entry)
    if (!is.matrix(value) || !is.numeric(value)) {
      .input_error(sprintf("%s must be a numeric matrix; received %s.",
        where, .msg_value(value)))
    }
    if (!identical(dim(value), c(q, q))) {
      .input_error(sprintf(paste0(
        "%s is %d x %d; the relation declares %s (%s), so every model RDM ",
        "must be %d x %d."
      ), where, nrow(value), ncol(value), .msg_count(q, "effect"),
        .msg_names(effects), q, q))
    }
    if (any(!is.finite(value))) {
      .input_error(sprintf("%s contains %s non-finite %s.", where,
        sum(!is.finite(value)),
        if (sum(!is.finite(value)) == 1L) "entry" else "entries"))
    }
    asymmetry <- max(abs(value - t(value)))
    if (asymmetry > 1e-12) {
      .input_error(sprintf(paste0(
        "%s is not symmetric; the largest difference between `m[i, j]` and ",
        "`m[j, i]` is %g. A dissimilarity between two effects has one value."
      ), where, asymmetry))
    }
    if (max(abs(diag(value))) > 1e-12) {
      .input_error(sprintf(paste0(
        "%s has a nonzero diagonal (largest |m[i, i]| is %g). Pass a ",
        "dissimilarity matrix, not a similarity matrix: an effect is at ",
        "distance zero from itself."
      ), where, max(abs(diag(value)))))
    }
    row_ids <- rownames(value)
    column_ids <- colnames(value)
    if (!is.null(row_ids) || !is.null(column_ids)) {
      if (is.null(row_ids) || is.null(column_ids)) {
        .input_error(sprintf(paste0(
          "%s names only its %s; name both axes with the relation's effects ",
          "(%s), or neither."
        ), where, if (is.null(row_ids)) "columns" else "rows",
          .msg_names(effects)))
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
          where, .msg_names(effects), paste(detail, collapse = "; ")))
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
#' @export
rsa <- function(x, models, nuisance = NULL, intercept = TRUE,
                component = c("total", "coherent", "configuration")) {
  source <- .self_geometry_source(x, "RSA")
  if (missing(models)) {
    .input_error(paste0(
      "`models` is required: pass one dissimilarity matrix over the ",
      "relation's effects, or a named list of them, for example ",
      "`rsa(plan, models = list(category = m))`."
    ))
  }
  if (!.is_flag(intercept)) {
    .input_error(sprintf("`intercept` must be TRUE or FALSE; received %s.",
      .msg_value(intercept)))
  }
  component <- match.arg(component)
  q <- length(source$effects)
  models <- .validate_rdm_models(models, source$effects, "models")
  nuisance <- .validate_rdm_models(nuisance, source$effects, "nuisance")
  if (any(names(models) %in% names(nuisance))) {
    .input_error(sprintf(paste0(
      "Model and nuisance names must be distinct; %s appears in both. Each ",
      "name becomes one coefficient column."
    ), .msg_names(intersect(names(models), names(nuisance)))))
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
    .input_error(.rsa_rank_deficiency_message(design, qr_design, intercept))
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
