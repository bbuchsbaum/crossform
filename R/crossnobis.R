# Crossvalidated Mahalanobis contrast views -------------------------------

#' Construct a fixed neural noise-precision metric
#'
#' `noise_precision()` is the semantically specific constructor for a known or
#' externally estimated precision that is fixed before effect evaluation. It
#' uses the same `neural_metric()` representation while recording that the
#' operator has inverse-noise-covariance meaning.
#'
#' @param value A finite symmetric positive-semidefinite precision matrix.
#' @param domain Exact neural feature domain.
#' @param support Ordered support identities for a local matrix.
#' @param covariance Optional retained inverse covariance.
#' @param tolerance Positive numerical validation tolerance.
#' @param provenance Additional compact provenance. Reserved semantic fields
#'   cannot be replaced.
#' @return An `effect_neural_metric` whose `$provenance$metric_role` is
#'   `"noise_precision"`, carrying the canonicalized `$value`, its `$domain`
#'   and `$support`, `$capabilities`, and a `$signature`. The recorded role is
#'   what lets [crossnobis()] name the result a Mahalanobis distance.
#' @seealso [crossnobis()] and [plan_geometry()] (its `metric` argument);
#'   [shrinkage_precision()] when the precision must instead be learned from
#'   residuals.
#' @family neural metrics
#' @examples
#' # A precision estimated outside crossform and fixed before evaluation.
#' domain <- abstract_domain(3, id = "noise-precision-example")
#' covariance <- matrix(c(1, 0.4, 0.1, 0.4, 1, 0.3, 0.1, 0.3, 1), 3, 3)
#' metric <- noise_precision(
#'   solve(covariance), domain, covariance = covariance,
#'   provenance = list(source = "resting-run residual covariance")
#' )
#' metric$provenance$metric_role
#' metric$capabilities$positive_definite
#'
#' # The role is what a crossnobis plan checks for, so a plain
#' # `neural_metric()` with the same numbers is deliberately not equivalent.
#' identical(
#'   neural_metric(solve(covariance), domain)$provenance$metric_role,
#'   metric$provenance$metric_role
#' )
#'
#' # The semantic provenance fields are reserved and cannot be overwritten.
#' refused <- try(
#'   noise_precision(diag(3), domain, provenance = list(metric_role = "other")),
#'   silent = TRUE
#' )
#' conditionMessage(attr(refused, "condition"))
#' @export
noise_precision <- function(value, domain, support = NULL,
                            covariance = NULL, tolerance = 1e-10,
                            provenance = list()) {
  .validate_effect_provenance(provenance, "noise-precision provenance")
  reserved <- c("metric_role", "estimator_status")
  if (!is.null(names(provenance)) && any(names(provenance) %in% reserved)) {
    .input_error("Noise-precision semantic provenance fields are reserved.")
  }
  neural_metric(
    value, domain, support = support, inverse = covariance,
    estimation = "fixed", tolerance = tolerance,
    provenance = c(list(
      metric_role = "noise_precision",
      estimator_status = "fixed_before_effect_evaluation"
    ), provenance)
  )
}

.crossnobis_plan_metric <- function(x) {
  .validate_geometry_plan(x)
  schedule <- x$metric_schedule
  if (!identical(schedule$kind, "fixed_metric_before_frame")) {
    .capability_refusal(paste0(
      "Crossnobis is a Mahalanobis reading, so it requires an explicit ",
      "noise-precision metric; this plan carries the implicit identity ",
      "metric, which is Euclidean and declares nothing about the noise."
    ),
      capability = "declared_noise_metric",
      namespace = "geometry_views",
      reasons = "implicit_identity_metric_is_not_a_noise_model",
      remedies = paste0(
        "Compile the plan with `plan_geometry(..., metric = ",
        "noise_precision(...))`, use `plan_crossnobis()` to learn the ",
        "precision from residuals, or read the same estimand without a ",
        "noise claim through `contrast_energy()`."
      )
    )
  }
  metric <- .validate_neural_metric(schedule$metric, deep = FALSE)
  if (!identical(metric$provenance$metric_role, "noise_precision")) {
    .capability_refusal(paste0(
      "Crossnobis interpretation requires a metric declared as a noise ",
      "precision; this plan's fixed metric was not constructed by ",
      "`noise_precision()` and carries no such role, so calling its output a ",
      "Mahalanobis distance would overstate what the metric means."
    ),
      capability = "declared_noise_metric",
      namespace = "geometry_views",
      reasons = "metric_role_is_not_noise_precision",
      remedies = paste0(
        "Build the metric with `noise_precision()`, or read the plan with ",
        "`contrast_energy()`, which makes no noise-model claim."
      )
    )
  }
  metric
}

# The learned counterpart of `.crossnobis_plan_metric()`. A recipe earns the
# Mahalanobis reading by being a residual-derived precision; the identity
# recipe is a legitimate schedule but declares nothing about the noise, so it
# gets the same refusal a plan with the implicit identity metric gets.
.crossnobis_learned_metric <- function(x) {
  schedule <- x$metric_schedule$schedule
  if (identical(schedule$recipe$kind, "identity")) {
    .capability_refusal(paste0(
      "Crossnobis is a Mahalanobis reading, so it requires a residual-",
      "derived precision; this plan's recipe is `identity_metric()`, which ",
      "is Euclidean and declares nothing about the noise."
    ),
      capability = "declared_noise_metric",
      namespace = "geometry_views",
      reasons = "identity_recipe_is_not_a_noise_model",
      remedies = paste0(
        "Compile the plan with `shrinkage_precision()` or ",
        "`diagonal_precision()`, or read the same estimand without a noise ",
        "claim through `contrast_energy()`."
      )
    )
  }
  schedule
}

.require_crossnobis_pairing <- function(over) {
  .validate_pairing(over)
  if (!identical(attr(over, "independence", exact = TRUE), "independent") ||
      !identical(attr(over, "estimate", exact = TRUE),
        "cross_generalized") || any(over$left == over$right)) {
    .input_error(sprintf(paste0(
      "Crossnobis requires cross-partition edges declared independent and ",
      "containing no self-products; this pairing declares independence `%s` ",
      "with estimate `%s`%s. The crossvalidated distance is unbiased only ",
      "because the two endpoints carry independent noise."
    ), attr(over, "independence", exact = TRUE),
      attr(over, "estimate", exact = TRUE),
      if (any(over$left == over$right)) {
        sprintf(" and contains %s",
          .msg_count(sum(over$left == over$right), "self-product"))
      } else {
        ""
      }))
  }
  invisible(TRUE)
}

#' Compile an on-demand learned-metric crossnobis plan
#'
#' The plan freezes the residual-statistics identity, training policy,
#' regularization recipe, evaluation pairing, and spatial support graph. It
#' does not retain one covariance or precision matrix per node and fold; each
#' local solve handle is derived on demand from canonical pair statistics.
#' The learned metric is itself random. Version 0.1 therefore returns signed
#' crossnobis point estimates but does not apply the fixed-metric analytic
#' covariance law in [rdm_sampling_covariance()] or claim calibrated intervals
#' for this plan; doing so requires propagation of metric-estimation
#' uncertainty (for example, an admitted LD-t specialization).
#'
#' `plan_crossnobis()` is a convenience wrapper around [plan_geometry()], not
#' a deprecated alias: it names the intent, and compiles exactly
#' `plan_geometry(x, at, over, metric = <recipe>, training = )`. The plan it
#' returns is an ordinary `effect_geometry_plan` whose metric schedule has
#' kind `learned_local_before_frame`, lowered by the geometry compiler like
#' any other. The pairing contract crossnobis requires -- independent,
#' cross-partition, self-product-free -- is enforced by [crossnobis()] when
#' the plan is read, which is where the fixed-metric route enforces it.
#'
#' @param x An `effect_relation_fit` with residual-block capability.
#' @param at A support-index-backed compiled spatial frame.
#' @param over Independent cross-partition evaluation edges.
#' @param metric An on-demand metric recipe such as `shrinkage_precision()`.
#' @param training A `metric_training_policy()`.
#' @param compute A sequential `compute_policy()`.
#' @param residual_workspace_bytes Positive budget used while accumulating
#'   canonical residual pair sufficient statistics. It changes cache capacity,
#'   never the canonical numerical tile shape. Defaults to
#'   `compute$workspace_bytes`, or 512 MiB when the policy declares none.
#' @return An `effect_geometry_plan` reusable across fixed contrasts. It
#'   carries the learned `$metric_schedule` wrapping the frozen schedule, the
#'   `$frame`, `$pairing`, and `$execution_hints`, and a
#'   `$scientific_plan_id` that identifies the estimand independently of
#'   execution choices.
#' @seealso [plan_geometry()], which this wraps; [crossnobis()] to read a
#'   contrast from this plan,
#'   [shrinkage_precision()] and [metric_training_policy()] for the metric
#'   declarations it freezes, and [residual_pair_statistics()] for the
#'   sufficient statistics it compiles.
#' @family geometry plans and views
#' @examples
#' set.seed(7)
#' domain <- abstract_domain(
#'   3, coordinates = cbind(x = 0:2, y = 0),
#'   id = "learned-crossnobis-example"
#' )
#' design <- cbind(intercept = 1, condition = rep(c(-0.5, 0.5), 4))
#' coefficients <- rbind(
#'   intercept = c(0, 0, 0),
#'   condition = c(0.6, -0.4, 0.2)
#' )
#' responses <- setNames(lapply(seq_len(3), function(run) {
#'   design %*% coefficients + matrix(rnorm(8 * 3, sd = 0.25), 8, 3)
#' }), paste0("run", seq_len(3)))
#' fit <- lm_relation_fit(
#'   responses, design, rbind(condition = c(0, 1)), domain = domain
#' )
#' plan <- plan_crossnobis(
#'   fit,
#'   compile_frame(searchlights(1.01), domain),
#'   pairing("run1", "run2", independence = "independent"),
#'   metric = shrinkage_precision(0.2)
#' )
#' result <- crossnobis(plan, c(condition = 1))
#' result
#' as.data.frame(result)
#' @export
plan_crossnobis <- function(
    x, at, over, metric = shrinkage_precision(),
    training = metric_training_policy("exclude_evaluation"),
    compute = compute_policy(),
    residual_workspace_bytes = NULL) {
  if (missing(at)) {
    .input_error(paste0(
      "`at` is required: pass a compiled frame from `compile_frame()`, for ",
      "example `compile_frame(searchlights(8), domain)`."
    ))
  }
  if (missing(over)) {
    .input_error(paste0(
      "`over` is required: pass a pairing from `cross_partitions()` or ",
      "`pairing()` declaring which partition products may be formed."
    ))
  }
  if (inherits(metric, "effect_metric_recipe") &&
      identical(.validate_metric_recipe(metric)$kind, "identity")) {
    .input_error(paste0(
      "A learned crossnobis plan requires a residual-derived precision ",
      "recipe; use `noise_precision(diag(...))` for a fixed identity metric."
    ))
  }
  plan_geometry(
    x, at, over, compute = compute, metric = metric, training = training,
    residual_workspace_bytes = residual_workspace_bytes
  )
}
