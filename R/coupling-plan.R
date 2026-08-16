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
#' @return An `effect_measurement_form` over the requested node pairs, with
#'   one block per edge in `$values` once read. Read it with
#'   [effect_coupling()], [connectivity()], [measurement_components()], or
#'   the other coupling views.
#' @seealso [plan_geometry()] for the plan this closes the other way,
#'   [measurement_form()] for the general constructor, and
#'   [variation_query()] for the experimental closure.
#' @family coupling and connectivity views
#'
#' @section Refusal:
#' A rectangular cross-axis plan signals an `effect_capability_refusal` with
#' capability `"self_form_coupling"` in namespace `"coupling_plans"`; branch on
#' it with [catch_refusal()].
#' @examples
#' # A plan already fixes the relation, the spatial frame, and the pairing.
#' # `coupling()` reuses all three and asks the adjoint question: keep two
#' # neural measurements open and close the experimental axis instead.
#' set.seed(4)
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
#'   cross_partitions(relation, independence = "independent")
#' )
#'
#' # `between` names frame nodes, so no second frame or pairing is declared.
#' form <- coupling(
#'   plan, cbind("r1", "r2"),
#'   by = variation_query(
#'     diag(2) - 0.5, relation$effect_space, "trial", "psd_variation"
#'   )
#' )
#' effect_coupling(form)$values[[1]]
#'
#' # Node names are resolved against the plan's own frame index.
#' unknown <- try(coupling(plan, cbind("r1", "r9"), by = variation_query(
#'   diag(2) - 0.5, relation$effect_space, "trial", "psd_variation"
#' )), silent = TRUE)
#' conditionMessage(attr(unknown, "condition"))
#' @export
coupling <- function(x, between, by, over = NULL,
                     mode = c("total", "coherent", "coherent_configuration"),
                     reducer = aggregate_first(),
                     compute = compute_policy(),
                     route = "auto") {
  .validate_geometry_plan(x)
  if (identical(x$codec, "rectangular")) {
    .capability_refusal(paste0(
      "Coupling closures require a self-form plan; this plan is rectangular ",
      "(two distinct experimental axes), and crossform 0.1 compiles no ",
      "cross-axis coupling closure."
    ),
      capability = "self_form_coupling",
      namespace = "coupling_plans",
      reasons = "rectangular_cross_axis_plan",
      remedies = paste0(
        "Build the plan from one relation by omitting `right` in ",
        "`plan_geometry()`."
      )
    )
  }
  if (missing(between)) {
    stop(paste0(
      "`between` is required: pass a two-column matrix of frame node pairs, ",
      "for example `cbind(from_ids, to_ids)` naming measurements of the ",
      "plan's frame."
    ), call. = FALSE)
  }
  if (missing(by)) {
    stop(paste0(
      "`by` is required: pass the `variation_query()` that says which ",
      "experimental variation the coupling is read over."
    ), call. = FALSE)
  }
  mode <- match.arg(mode)
  frame <- measurement_frame(x$frame, mode = mode)
  node_ids <- frame$node_ids
  resolve <- function(value, side) {
    if (is.numeric(value)) {
      if (anyNA(value) || any(value < 1) || any(value > length(node_ids)) ||
          any(value %% 1 != 0)) {
        stop(sprintf(paste0(
          "`between` %s indices must be whole numbers in 1:%d, one per frame ",
          "measurement; received %s."
        ), side, length(node_ids), .msg_names(format(value))), call. = FALSE)
      }
      return(node_ids[as.integer(value)])
    }
    value <- as.character(value)
    if (anyNA(value) || any(!value %in% node_ids)) {
      unknown <- unique(value[is.na(value) | !value %in% node_ids])
      stop(sprintf(paste0(
        "`between` %s names %s, which %s a measurement of the plan's frame. ",
        "The frame has %s: %s."
      ), side, .msg_names(unknown),
        if (length(unknown) == 1L) "is not" else "are not",
        .msg_count(length(node_ids), "measurement"), .msg_names(node_ids)),
        call. = FALSE)
    }
    value
  }
  if (is.data.frame(between)) between <- as.matrix(between)
  if (!is.matrix(between) || ncol(between) != 2L || nrow(between) < 1L) {
    stop(sprintf(paste0(
      "`between` must be a nonempty two-column matrix of frame node pairs ",
      "(from, to); received %s."
    ), .msg_value(between)), call. = FALSE)
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
