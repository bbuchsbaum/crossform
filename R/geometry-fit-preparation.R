# Plan metadata collectors and reduced execution preparation (layer 5).
# Shared by fitting and scoring; this file never calls either public verb.

.geometry_training_support <- function(plan, upstream = list()) {
  .check_class(plan, "effect_geometry_plan", "plan", from = "plan_geometry()")
  relation <- .validate_relation(plan$task$left_relation)
  over <- .validate_pairing(plan$pairing)
  # Zero-weight edges are declared but are not numerical dependencies.
  products <- as.data.frame(over[over$weight > 0, c("left", "right", "weight"), drop = FALSE])
  attributes(products) <- list(names = c("left", "right", "weight"),
    row.names = seq_len(nrow(products)), class = "data.frame")
  used <- sort(unique(c(products$left, products$right)), method = "radix")
  capabilities <- .relation_source_capabilities(relation)
  endpoints <- stats::setNames(lapply(used, function(p) {
    list(origin = relation$sources[[p]]$origin,
      source_revision = capabilities[[p]]$stable_revision,
      extractor_signature = .sha256_signature(relation$extractors[[p]], "extractor-sha256:"))
  }), used)
  metadata <- list(directed = attr(over, "directed"), self_pairs = attr(over, "self_pairs"),
    independence = attr(over, "independence"), generalizes_over = attr(over, "generalizes_over"))
  .new_geometry_training_support(relation$provenance$observation_origins,
    endpoints, products, metadata, upstream)
}

.geometry_prediction_target <- function(plan, component) {
  index <- .execution_measurement_index(plan$frame)
  if ((!is.character(index) && !is.numeric(index)) || length(index) != plan$measurements ||
      anyNA(index) || anyDuplicated(index) || any(!nzchar(as.character(index)))) {
    .capability_refusal("Prediction requires one unique stable identifier per spatial measurement.",
      capability = "matching_prediction_target", namespace = "predictive_geometry",
      reasons = "ambiguous_measurement_ids", remedies = "Declare unique measurement IDs on the spatial frame.")
  }
  node_at <- .frame_metric_node_accessor(plan$frame)
  node_signatures <- vapply(seq_len(plan$measurements), function(i) {
    node <- node_at(i)
    .sha256_signature(list(id = as.character(index[[i]]), support = node$support,
      weight = node$weight), "prediction-node-sha256:")
  }, "")
  list(effect_space = plan$task$left_relation$effect_space,
    domain = plan$frame$domain, metric_signature = plan$metric_schedule$signature,
    normalization = plan$frame$normalization, component = component,
    generalizes_over = attr(plan$pairing, "generalizes_over"),
    measurement_ids = as.character(index), measurement_signatures = node_signatures,
    index = index)
}

.geometry_fit_prepare <- function(plan, basis, component = "total", upstream = list()) {
  .check_class(plan, "effect_geometry_plan", "plan", from = "plan_geometry()")
  .self_geometry_source(plan, "Predictive geometry")
  basis <- .validate_model_basis(basis)
  component <- match.arg(component, c("total", "coherent", "configuration"))
  if (isTRUE(attr(plan$pairing, "directed"))) {
    .capability_refusal("The symmetric packed executor requires an undirected self-form pairing.",
      capability = "predictive_symmetric_pairing", namespace = "predictive_geometry",
      reasons = "directed_full_form_not_admitted",
      remedies = "Declare the symmetric self-form estimand with an undirected pairing.")
  }
  if (!plan$metric_schedule$kind %in% c("implicit_identity_before_frame", "fixed_metric_before_frame")) {
    .capability_refusal("Predictive geometry currently requires a fixed native neural metric.",
      capability = "predictive_fixed_metric", namespace = "predictive_geometry",
      reasons = "metric_schedule_not_lowerable",
      remedies = "Use the implicit identity or a declared fixed neural_metric() with native composition.")
  }
  if (!is.null(plan$metric_schedule$metric) &&
      !identical(plan$metric_schedule$metric$estimation, "fixed")) {
    .capability_refusal("A learned neural metric needs its own admitted training schedule for prediction.",
      capability = "predictive_fixed_metric", namespace = "predictive_geometry",
      reasons = "unadmitted_learned_neural_metric",
      remedies = "Use an externally fixed neural metric in this version.")
  }
  if (!is.null(plan$metric_schedule$metric) &&
      !isTRUE(plan$metric_schedule$metric$capabilities$positive_definite)) {
    .capability_refusal("The existing complete-form executor requires an SPD metric for its coherent decomposition.",
      capability = "predictive_complete_form_metric", namespace = "predictive_geometry",
      reasons = "complete_form_metric_requires_spd",
      remedies = "Use a fixed positive definite neural metric or the implicit identity.")
  }
  rel <- plan$task$left_relation
  space <- rel$effect_space
  if (!identical(space$coordinates, basis$conditions) ||
      !identical(unique(unname(space$units)), basis$centering$units) ||
      !identical(unique(unname(space$scale)), basis$centering$scale) ||
      (identical(basis$centering$declaration$kind, "effect_space") &&
       !identical(space$signature, basis$centering$declaration$signature))) {
    .capability_refusal("The model declaration does not match the plan's effect meaning, order, units or scale.",
      capability = "matching_model_coordinates", namespace = "predictive_geometry",
      reasons = "model_effect_space_mismatch",
      remedies = "Build model_basis(conditions = relation$effect_space) for this target.")
  }
  if (.is_model_coordinate_space(space)) {
    .capability_refusal("Pass the plan on the original effects; this plan is already lowered.",
      capability = "original_prediction_coordinates", namespace = "predictive_geometry",
      reasons = "already_lowered_prediction_plan", remedies = "Use the original condition-effect plan.")
  }
  training_support <- .geometry_training_support(plan, upstream)
  used <- names(training_support$endpoints)
  reduced_relation <- .relation_lowered(.relation_subset(rel, used), basis)
  products <- training_support$products
  over <- pairing(products$left, products$right, weight = products$weight,
    directed = attr(plan$pairing, "directed"), self_pairs = attr(plan$pairing, "self_pairs"),
    independence = if (identical(attr(plan$pairing, "independence"), "undeclared")) NULL else
      attr(plan$pairing, "independence"), generalizes_over = attr(plan$pairing, "generalizes_over"))
  reduced_plan <- plan_geometry(reduced_relation, plan$frame, over,
    compute = plan$compute, metric = plan$metric_schedule$metric)
  if (!identical(reduced_plan$measurements, plan$measurements) ||
      !identical(reduced_plan$metric_schedule$signature, plan$metric_schedule$signature) ||
      !identical(reduced_plan$compute, plan$compute) ||
      reduced_plan$packed_width != basis$dimension * (basis$dimension + 1L)/2L) {
    .invariant_error("Predictive lowering changed the target beyond its effect coordinates and used partitions.")
  }
  list(plan = reduced_plan, basis = basis, support = training_support,
    target = .geometry_prediction_target(plan, component), component = component)
}

.geometry_fit_materialize <- function(prepared, storage = "memory", storage_path = NULL) {
  materialize_geometry(prepared$plan, storage = storage, storage_path = storage_path)
}
