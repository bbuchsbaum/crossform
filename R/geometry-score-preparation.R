# Independent evaluation admission (layer 5) -------------------------------
# All checks here use metadata. No training/evaluation values are opened.

.geometry_prediction_match <- function(training_target, evaluation_target) {
  .geometry_prediction_target_check(training_target)
  .geometry_prediction_target_check(evaluation_target)
  refuse <- function(reason, message) .capability_refusal(message,
    capability = "matching_prediction_target", namespace = "predictive_geometry",
    reasons = reason, remedies = "Evaluate the same effect coordinates, fixed metric and spatial measurements as the fitted target.")
  fields <- c("effect_space", "domain", "metric_signature", "normalization", "component", "generalizes_over")
  for (field in fields) {
    if (!identical(training_target[[field]], evaluation_target[[field]])) {
      refuse(paste0("prediction_", field, "_mismatch"),
        sprintf("The evaluation %s differs from the frozen prediction target.", field))
    }
  }
  mapping <- match(training_target$measurement_ids, evaluation_target$measurement_ids)
  if (length(mapping) != length(evaluation_target$measurement_ids) || anyNA(mapping)) {
    refuse("prediction_measurements_mismatch", "Evaluation must identify exactly the fitted spatial measurements.")
  }
  if (!identical(unname(training_target$measurement_signatures),
      unname(evaluation_target$measurement_signatures[mapping]))) {
    refuse("prediction_measurement_support_mismatch", "A spatial measurement has different support or weights at evaluation.")
  }
  as.integer(mapping)
}

.geometry_score_prepare <- function(fit, plan, component = fit$parameters$component) {
  fit <- .validate_geometry_fit(fit)
  .check_class(plan, "effect_geometry_plan", "plan", from = "plan_geometry()")
  .self_geometry_source(plan, "Predictive evaluation")
  target <- .geometry_prediction_target(plan, component)
  mapping <- .geometry_prediction_match(fit$target, target)
  support <- .geometry_training_support(plan)
  .geometry_support_independence(fit$training, support)
  prepared <- .geometry_fit_prepare(plan, fit$model, component)
  list(prepared = prepared, mapping = mapping,
    validation = list(contract = "predictive-geometry-v1",
      assumption = "Conditional unbiasedness of evaluation geometry given all fitting and selection dependencies.",
      training = fit$training$signature, evaluation = support$signature,
      origin_scope = support$scope_signature, independent = TRUE,
      sampling_assumption = support$assumption, mapping = mapping))
}
