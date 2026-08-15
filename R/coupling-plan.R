# Adjoint-side closure from the plan vocabulary ------------------------------

#' Take the adjoint-side coupling closure of a geometry plan
#'
#' A geometry plan closes the neural boundary and leaves the experimental
#' axes open; `coupling()` asks the adjoint question from the same plan
#' vocabulary: close the experimental boundary with a declared query and
#' leave two neural measurements open. The plan's spatial frame supplies the
#' measurement legs and the plan's pairing supplies the partition products,
#' so no second frame or pairing declaration is required.
#'
#' The measurement pipeline this compiles into is deliberately small-node:
#' frames whose dense payload exceeds the admitted ceiling refuse before
#' reading data, so searchlight-resolved coupling remains explicitly
#' unadmitted rather than silently truncated.
#'
#' @param x A self-form `effect_geometry_plan`.
#' @param between A two-column matrix or data frame of frame node pairs,
#'   given as node names from the plan's frame index or as node indices.
#'   Each row becomes one ordered measurement edge.
#' @param by The experimental closure: a [variation_query()].
#' @param over Optional `effect_pairing`; defaults to the plan's own pairing.
#' @param mode Measurement mode passed to [measurement_frame()].
#' @param reducer Partition reducer for the measurement pipeline.
#' @param compute A `compute_policy()`.
#' @param route Measurement execution route; `"auto"` selects one.
#' @return An effect measurement form; read it with [effect_coupling()],
#'   [connectivity()], [measurement_components()], or the other coupling
#'   views.
#' @examples
#' domain <- abstract_domain(4, id = "coupling-plan-example")
#' relation <- relation(
#'   list(
#'     run1 = matrix(rnorm(8), 2, 4, dimnames = list(c("a", "b"), NULL)),
#'     run2 = matrix(rnorm(8), 2, 4, dimnames = list(c("a", "b"), NULL))
#'   ),
#'   effects = effect_space(c("a", "b")), domain = domain
#' )
#' plan <- plan_geometry(
#'   relation, compile_frame(regions(c("r1", "r1", "r2", "r2")), domain),
#'   cross_partitions(relation)
#' )
#' form <- coupling(
#'   plan, cbind("r1", "r2"),
#'   by = variation_query(
#'     diag(2) - 0.5, relation$effect_space, "trial", "psd_variation"
#'   )
#' )
#' effect_coupling(form)
#' @export
coupling <- function(x, between, by, over = NULL,
                     mode = c("total", "coherent", "coherent_configuration"),
                     reducer = aggregate_first(),
                     compute = compute_policy(),
                     route = "auto") {
  .validate_geometry_plan(x)
  if (identical(x$codec, "rectangular")) {
    stop(paste0(
      "Coupling closures currently require a self-form plan; rectangular ",
      "cross-axis coupling is not yet compiled."
    ), call. = FALSE)
  }
  mode <- match.arg(mode)
  frame <- measurement_frame(x$frame, mode = mode)
  node_ids <- frame$node_ids
  resolve <- function(value, side) {
    if (is.numeric(value)) {
      if (anyNA(value) || any(value < 1) || any(value > length(node_ids)) ||
          any(value %% 1 != 0)) {
        stop(sprintf(
          "`between` %s indices must identify frame measurements.", side
        ), call. = FALSE)
      }
      return(node_ids[as.integer(value)])
    }
    value <- as.character(value)
    if (anyNA(value) || any(!value %in% node_ids)) {
      stop(sprintf(paste0(
        "`between` %s names must identify measurements of the plan's ",
        "frame."
      ), side), call. = FALSE)
    }
    value
  }
  if (is.data.frame(between)) between <- as.matrix(between)
  if (!is.matrix(between) || ncol(between) != 2L || nrow(between) < 1L) {
    stop("`between` must be a two-column matrix of frame node pairs.",
      call. = FALSE)
  }
  edges <- edge_frame(
    resolve(between[, 1L], "from"),
    resolve(between[, 2L], "to"),
    frame
  )
  if (is.null(over)) over <- x$pairing
  measurement_form(
    x$task$left_relation, edges, by, over,
    reducer = reducer, compute = compute, route = route
  )
}
