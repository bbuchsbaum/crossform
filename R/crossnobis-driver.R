# The crossnobis runtime entry ------------------------------------------------
#
# `R/crossnobis.R` is the plan layer: `noise_precision()`, `plan_crossnobis()`,
# and the metric- and pairing-role checks `crossnobis()` applies. What lives
# here is the public entry that reads a contrast off a compiled plan.
#
# There is one runtime. Until B3 this file also carried a second one -- a
# planned receipt, a private driver, and a dispatch on a retired
# `effect_crossnobis_plan` class -- for the learned-metric route. That route is
# now an ordinary compiler lowering (B2), so `crossnobis()` validates the
# metric role and hands the plan to `.run_geometry_compiler()` exactly as the
# fixed-metric route always did. The plan constructs nothing here and this file
# constructs no plan; the direction is plan -> runtime, one way.

#' Evaluate a signed local crossnobis contrast
#'
#' The result is the query-first evidence
#' `tr(c c^T B_a K B_b^T)` aggregated over the plan's independent partition
#' edges. This function is a validating view over the ordinary geometry
#' compiler; it does not introduce a second numerical engine. Negative finite
#' estimates are retained.
#'
#' @param x An `effect_geometry_plan` carrying either an explicit fixed
#'   `noise_precision()` metric, or a provenance-frozen learned metric
#'   schedule compiled from a recipe by [plan_crossnobis()] or
#'   `plan_geometry(metric = )`.
#' @param weights One finite contrast weight per experimental effect.
#' @return An `effect_crossnobis_view` whose `$values` holds one signed
#'   crossvalidated squared Mahalanobis value per spatial measurement, with
#'   the aligned `$contrast`, the named `$estimand`, the `$metric` and
#'   `$pairing` identities, `$index`, and the executed `$receipt`.
#' @section Structure:
#' One signed value per spatial measurement, alongside the declarations that
#' make it a Mahalanobis reading.
#'
#' - `$values`: the signed crossvalidated squared Mahalanobis contrast, one
#'   per measurement, in `$index` order. Negative estimates are retained.
#' - `$contrast`: the weights, reordered to the relation's effect order and
#'   named.
#' - `$estimand`: the named estimand,
#'   `"crossvalidated_squared_mahalanobis_contrast"`.
#' - `$metric`: the signature of the noise-precision metric the values were
#'   read under.
#' - `$pairing`: the identity of the independent partition edges they were
#'   generalized over.
#' - `$index`: the measurement identifiers, one per value, carried from the
#'   frame's `$index$measurement`.
#' - `$receipt`: the execution receipt for the run that produced the values.
#'
#' The `$metadata` block and any other element not listed here are internal
#' and may change.
#' @seealso [noise_precision()] for the fixed metric a geometry plan needs,
#'   [plan_crossnobis()] for the learned-metric route, and
#'   [contrast_energy()] for the decomposed reading of the same estimand.
#' @family geometry plans and views
#'
#' @section Refusal:
#' A plan carrying the implicit identity metric, or a fixed metric that was not
#' built by [noise_precision()], signals an `effect_capability_refusal` with
#' capability `"declared_noise_metric"` in namespace `"geometry_views"`.
#' Inspect it with [catch_refusal()].
#'
#' @section One estimand, two views:
#' On a fixed-metric geometry plan, `crossnobis(x, weights)` is the named
#' Mahalanobis reading of exactly `contrast_energy(x, weights)$total`: the same
#' compiled estimand, exposed as a single signed value. When the analysis
#' also needs the signed endpoint marginals or the exact
#' coherent/configuration decomposition of the same quantity, call
#' [contrast_energy()] on the same plan; every component it returns inherits the
#' plan's fixed metric.
#' @examples
#' domain <- abstract_domain(3, id = "crossnobis-example")
#' run1 <- rbind(a = c(1, 0, 0), b = c(0, 1, 0))
#' run2 <- rbind(a = c(1.1, 0, 0), b = c(0, 0.9, 0))
#' relation <- relation(list(run1 = run1, run2 = run2), domain = domain)
#' plan <- plan_geometry(
#'   relation,
#'   compile_frame(whole_brain(), domain),
#'   cross_partitions(relation, independence = "independent"),
#'   metric = noise_precision(diag(3), domain)
#' )
#' result <- crossnobis(plan, c(a = 1, b = -1))
#' result
#' as.data.frame(result)
#' @export
crossnobis <- function(x, weights) {
  if (missing(weights)) {
    .input_error(paste0(
      "`weights` is required: pass one finite weight per experimental ",
      "effect, for example `crossnobis(plan, c(face = 1, house = -1))`."
    ))
  }
  if (!inherits(x, "effect_geometry_plan")) {
    .input_error(sprintf(paste0(
      "`x` must be an `effect_geometry_plan` from `plan_geometry()` or ",
      "`plan_crossnobis()`; received %s."
    ), .msg_value(x)))
  }
  # One argument type, validated rather than dispatched on. Both admitted
  # schedules make the same claim -- that the operator between the two
  # endpoints is a noise precision -- and the fixed and learned routes differ
  # only in where that claim is recorded.
  learned <- identical(
    x$metric_schedule$kind, "learned_local_before_frame"
  )
  metric_identity <- if (learned) {
    .crossnobis_learned_metric(x)$signature
  } else {
    .crossnobis_plan_metric(x)$signature
  }
  .require_crossnobis_pairing(x$pairing)
  weights <- .align_contrast(
    weights, x$task$left_relation$effect_space$coordinates
  )
  query <- bilinear_query(tcrossprod(weights))
  # The internal runner rather than `evaluate_geometry()`: this is the plan
  # layer, and the public entry point sits above the views it feeds. The
  # signed weights are an executor hint on the learned route, where the
  # scheduled kernel contracts `c B_a` against `c B_b` instead of applying
  # the packed operator; they do not move the estimand identity.
  evaluated <- .run_geometry_compiler(
    x, query = query, component = "total",
    signed_query = if (learned) weights else NULL
  )
  values <- drop(evaluated$values)
  structure(list(
    values = values,
    contrast = weights,
    estimand = "crossvalidated_squared_mahalanobis_contrast",
    metric = metric_identity,
    pairing = .metric_pairing_identity(x$pairing),
    index = evaluated$index,
    receipt = evaluated$receipt,
    metadata = evaluated$metadata
  ), class = "effect_crossnobis_view")
}
