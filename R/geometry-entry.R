# Public geometry entry points -----------------------------------------------
#
# Layer 5. These are the two exported ways to run a geometry plan. They also
# carry the relation compatibility form, `f(relation, at, over)`, which builds
# the plan first: plan construction belongs to the caller's side of the
# boundary, so the compiler below only ever consumes a plan it is handed.

#' Materialize a complete cross-generalized geometry
#'
#' `materialize_geometry()` explicitly materializes complete total and
#' coherent packed geometry. Prefer [plan_geometry()] followed by
#' [evaluate_geometry()] when a fixed query is sufficient. The relation
#' compatibility form compiles that same plan first; it does not use a second
#' execution path.
#'
#' @param x An `effect_geometry_plan` or, for compatibility, an
#'   `effect_relation`.
#' @param at A compiled additive `effect_frame`, required only with a relation.
#' @param over An `effect_pairing`, required only with a relation.
#' @param storage Either `"memory"` or `"block"`.
#' @param storage_path Durable directory for block-backed geometry.
#' @param compute A sequential `compute_policy()` when `x` is a relation. A
#'   compiled plan already owns this choice.
#' @param reporter Optional nonsemantic coordinator-side event reporter.
#' @return An `effect_geometry` holding the packed `$total` and `$coherent`
#'   stores, the signed endpoint `$marginals`, the `$effects` names, the
#'   measurement `$index`, and the execution `$receipt`. Read it with
#'   [geometry_component()], [query_geometry()], or the named views.
#'
#'   `$result_capability` is `"complete_form"` here and `"query_only"` on a
#'   view; it is the one field whose vocabulary is the same across every
#'   result kind, so branch on it rather than on the displayed
#'   `$completeness`, which reads `"full"` on an `effect_geometry` and
#'   `"complete"` on an `effect_measurement_form`. Those two words are frozen:
#'   `completeness` is hashed into the measurement contract signature and into
#'   evidence-task identity, so reconciling them would invalidate recorded
#'   identities.
#' @seealso [plan_geometry()] then [evaluate_geometry()] for the query-first
#'   route; [geometry_spectrum()], which requires a complete geometry.
#' @family geometry plans and views
#' @examples
#' domain <- abstract_domain(4, id = "materialize-example")
#' relation <- relation(
#'   list(run1 = rbind(a = c(1, 0, 2, 1), b = c(0, 1, 1, 0)),
#'        run2 = rbind(a = c(1.1, 0.1, 1.9, 0.8), b = c(0.1, 0.9, 1.2, 0.2))),
#'   domain = domain
#' )
#' plan <- plan_geometry(
#'   relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
#'   cross_partitions(relation, independence = "independent")
#' )
#'
#' # Materialize only when several different readouts will share the same
#' # packed geometry; a single fixed query does not need this step.
#' geometry <- materialize_geometry(plan)
#' geometry$effects
#' geometry_component(geometry, "total")
#'
#' # The complete object unlocks readouts that a query-only view cannot give,
#' # such as the signed eigenvalue spectrum.
#' geometry_spectrum(geometry)$values
#' @export
materialize_geometry <- function(x, at = NULL, over = NULL,
                                 storage = c("memory", "block"),
                                 storage_path = NULL, compute = NULL,
                                 reporter = NULL) {
  storage <- match.arg(storage)
  .run_geometry_compiler(x, at, over, compute, storage, storage_path,
    query = NULL, component = NULL, reporter = reporter)
}

#' Evaluate a fixed query without materializing complete geometry
#'
#' `evaluate_geometry()` is the query-first route: the requested view is
#' contracted during execution, so complete packed geometry is never
#' allocated. Use it whenever a fixed query answers the question; use
#' [materialize_geometry()] only when several readouts must share one stored
#' geometry.
#'
#' @param x,at,over,compute,reporter As in [materialize_geometry()].
#' @param query A fixed `bilinear_query()`, an axis-bound [pair_query()] for
#'   rectangular cross-axis plans, or a packed-coordinate query matrix.
#' @param component One of `total`, `coherent`, or `configuration`.
#' @return A query-only `effect_view`: `$values` has one row per measurement
#'   and one column per query column, with `$component`, `$index`, and the
#'   execution `$receipt`. It carries no packed geometry, so
#'   [geometry_component()] does not apply to it.
#' @seealso [plan_geometry()] to build the plan, [contrast_energy()],
#'   [rdm()], and [rsa()] for named views of the same plan.
#' @family geometry plans and views
#' @examples
#' domain <- abstract_domain(4, id = "evaluate-example")
#' relation <- relation(
#'   list(run1 = rbind(a = c(1, 0, 2, 1), b = c(0, 1, 1, 0)),
#'        run2 = rbind(a = c(1.1, 0.1, 1.9, 0.8), b = c(0.1, 0.9, 1.2, 0.2))),
#'   domain = domain
#' )
#' plan <- plan_geometry(
#'   relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
#'   cross_partitions(relation, independence = "independent")
#' )
#'
#' # A fixed contrast read straight off the plan, with no packed geometry.
#' view <- evaluate_geometry(plan, query = bilinear_query(tcrossprod(c(1, -1))))
#' view
#' as.data.frame(view)
#'
#' # A query is mandatory: query-first execution will not guess one.
#' missing_query <- try(evaluate_geometry(plan), silent = TRUE)
#' conditionMessage(attr(missing_query, "condition"))
#' @export
evaluate_geometry <- function(x, at = NULL, over = NULL, query = NULL,
                              component = c("total", "coherent", "configuration"),
                              compute = NULL, reporter = NULL) {
  if (inherits(x, "effect_geometry_plan") && is.null(query) &&
      !is.null(at) && is.null(over)) {
    query <- at
    at <- NULL
  }
  if (is.null(query)) {
    .input_error("Query-first execution requires an explicit fixed `query`.")
  }
  component <- match.arg(component)
  .run_geometry_compiler(x, at, over, compute, "memory", NULL,
    query = query, component = component, reporter = reporter)
}
