# Declared edge cross-fitting (layer 5, internal orchestration) ------------
# This estimates the fitted procedure at the fold training sizes, not the
# risk of an all-data refit. Every fold reuses the public fixed fit/score.

.geometry_prediction_pair_plan <- function(plan, edges, compute = plan$compute) {
  used <- unique(c(edges$left, edges$right))
  rel <- .relation_subset(plan$task$left_relation, used)
  independence <- attr(plan$pairing, "independence")
  if (identical(independence, "undeclared")) independence <- NULL
  plan_geometry(rel, plan$frame, pairing(edges$left, edges$right, weight = edges$weight,
    directed = FALSE, independence = independence,
    generalizes_over = attr(plan$pairing, "generalizes_over")),
    metric = plan$metric_schedule$metric, compute = compute)
}

.geometry_crossfit_schedule <- function(plan, basis, component = "total") {
  .geometry_fit_prepare(plan, basis, component) # metadata-only route admission
  edges <- as.data.frame(plan$pairing)
  edges <- edges[edges$weight > 0, c("left", "right", "weight"), drop = FALSE]
  rownames(edges) <- NULL
  if (any(edges$left == edges$right)) {
    .capability_refusal("Predictive cross-fitting requires independent cross-partition edges.",
      capability = "predictive_crossfit_schedule", namespace = "predictive_geometry",
      reasons = "crossfit_self_products_not_admitted", remedies = "Declare independent cross-partition products.")
  }
  folds <- lapply(seq_len(nrow(edges)), function(i) {
    held_out <- unlist(edges[i, c("left", "right")], use.names = FALSE)
    disjoint <- !edges$left %in% held_out & !edges$right %in% held_out
    if (!any(disjoint)) {
      .capability_refusal(sprintf("Evaluation edge %s / %s has no disjoint training edge in the declared pairing.",
        held_out[[1L]], held_out[[2L]]), capability = "predictive_crossfit_schedule",
        namespace = "predictive_geometry", reasons = "insufficient_disjoint_edges",
        remedies = "Declare a pairing in which every evaluation edge has disjoint training products.")
    }
    training <- .geometry_prediction_pair_plan(plan, edges[disjoint, , drop = FALSE])
    evaluation <- .geometry_prediction_pair_plan(plan, edges[i, , drop = FALSE])
    train_support <- .geometry_training_support(training)
    eval_support <- .geometry_training_support(evaluation)
    .geometry_support_independence(train_support, eval_support)
    list(id = sprintf("edge%04d", i), weight = edges$weight[[i]],
      training = training, evaluation = evaluation,
      training_support = train_support, evaluation_support = eval_support)
  })
  list(edges = edges, folds = folds, target = .geometry_prediction_target(plan, component),
    estimand = "Expected predictive gain of the declared fitting procedure at these fold training sizes.")
}

.geometry_crossfit <- function(plan, basis, weights = NULL, rank, penalty,
                               component = "total", storage = c("memory", "block"),
                               storage_path = NULL, row_block = 128L) {
  storage <- match.arg(storage)
  basis <- .validate_model_basis(basis)
  pool <- .model_basis_pool(basis, weights)
  rank <- .check_count(rank, "rank", min = 0L, max = .Machine$integer.max)
  .check_number(penalty, "penalty", nonnegative = TRUE)
  row_block <- .check_count(row_block, "row_block", max = .Machine$integer.max)
  schedule <- .geometry_crossfit_schedule(plan, basis, component)
  nf <- length(schedule$folds); n <- plan$measurements
  # Retained score tables, rank paths and metadata coexist across folds.
  # Reserve them before either fitter/scorer performs its own admission.
  k <- min(rank, pool$dimension)
  factor_width <- basis$dimension * k + k + 2 * basis$dimension
  retained_bound <- 4 * 512 * as.double(n) * nf +
    2 * as.double(utils::object.size(schedule)) +
    if (storage == "memory") 2 * 8 * as.double(n) * nf * factor_width else 0
  remaining <- if (is.null(plan$compute$workspace_bytes)) NULL else plan$compute$workspace_bytes - retained_bound
  if (!is.null(remaining) && remaining <= 0) {
    .capability_refusal("The cross-fit score tables and schedule exceed the declared workspace.",
      capability = "predictive_crossfit_workspace", namespace = "predictive_geometry",
      reasons = "crossfit_workspace_budget_exceeded", remedies = "Increase the declared workspace or reduce the number of measurements/folds.")
  }
  compute <- compute_policy(block_features = plan$compute$block_features, workspace_bytes = remaining)
  created <- FALSE; success <- FALSE
  if (storage == "block") {
    .check_string(storage_path, "storage_path", what = "a new durable cross-fit directory")
    if (file.exists(storage_path)) .input_error("Refusing to overwrite an existing cross-fit directory.")
    if (!dir.create(storage_path, recursive = TRUE)) .input_error("Cannot create the cross-fit directory.")
    storage_path <- normalizePath(storage_path, mustWork = TRUE); created <- TRUE
  } else if (!is.null(storage_path)) .input_error("`storage_path` is used only with block storage.")
  on.exit(if (created && !success) unlink(storage_path, recursive = TRUE), add = TRUE)
  scores <- vector("list", nf)
  fits <- vector("list", nf)
  for (i in seq_len(nf)) {
    fold <- schedule$folds[[i]]
    train <- .geometry_prediction_pair_plan(fold$training, as.data.frame(fold$training$pairing), compute)
    test <- .geometry_prediction_pair_plan(fold$evaluation, as.data.frame(fold$evaluation$pairing), compute)
    fit_path <- if (storage == "block") file.path(storage_path, paste0(fold$id, "-fit")) else NULL
    score_path <- if (storage == "block") file.path(storage_path, paste0(fold$id, "-score")) else NULL
    fit <- fit_geometry(train, basis, weights = pool$weights, rank = rank, penalty = penalty,
      component = component, storage = storage, storage_path = fit_path, row_block = row_block)
    scores[[i]] <- score_geometry(fit, test, storage = storage, storage_path = score_path, row_block = row_block)
    fits[[i]] <- fit
  }
  names(scores) <- vapply(schedule$folds, `[[`, "", "id")
  names(fits) <- names(scores)
  edge_weights <- schedule$edges$weight
  reduce <- function(field) as.numeric(do.call(cbind, lapply(scores, function(x) x$table[[field]])) %*% edge_weights)
  value <- list(table = data.frame(measurement = schedule$target$measurement_ids,
      gain = reduce("gain"), inner_product = reduce("inner_product"),
      prediction_norm_sq = reduce("prediction_norm_sq")),
    rank_by_fold = do.call(cbind, lapply(scores, function(x) x$table$rank)),
    edges = schedule$edges, scores = scores, fits = fits, target = schedule$target,
    estimand = schedule$estimand, model_signature = basis$signature,
    parameters = list(rank = rank, penalty = penalty, weights = pool$weights, component = component),
    workspace = list(retained_bound = retained_bound, fold_budget = remaining,
      budget = plan$compute$workspace_bytes), complete = TRUE)
  value$signature <- .sha256_signature(.geometry_crossfit_semantic(value), "geometry-crossfit-sha256:")
  result <- structure(value, class = "effect_geometry_crossfit")
  .validate_geometry_crossfit(result)
  if (storage == "block") {
    pending <- file.path(storage_path, "crossfit.pending.rds")
    saveRDS(result, pending)
    if (!file.rename(pending, file.path(storage_path, "crossfit.rds"))) .input_error("Cannot finalize the cross-fit record.")
  }
  success <- TRUE
  result
}

.geometry_crossfit_semantic <- function(x) {
  value <- unclass(x[setdiff(names(x), c("signature", "fits", "scores"))])
  value$fits <- vapply(x$fits, `[[`, "", "signature")
  value$scores <- vapply(x$scores, `[[`, "", "signature")
  value
}

.validate_geometry_crossfit <- function(x) {
  fields <- c("table", "rank_by_fold", "edges", "scores", "fits", "target", "estimand",
    "model_signature", "parameters", "workspace", "complete", "signature")
  if (!.sealed_fields(x, "effect_geometry_crossfit", fields) || !isTRUE(x$complete) ||
      !is.list(x$scores) || !length(x$scores) || length(x$scores) != nrow(x$edges) ||
      any(!is.finite(x$edges$weight)) || any(x$edges$weight <= 0) || abs(sum(x$edges$weight) - 1) > 1e-12) {
    .contract_error("A cross-fit record requires all declared fold scores and their weights.")
  }
  lapply(x$scores, .validate_geometry_score)
  if (length(x$fits) != length(x$scores)) .contract_error("Every cross-fit score requires its frozen fitted prediction.")
  lapply(x$fits, .validate_geometry_fit)
  if (!identical(vapply(x$fits, `[[`, "", "prediction_id"),
      vapply(x$scores, `[[`, "", "prediction_id"))) .contract_error("Cross-fit scores and frozen predictions disagree.")
  .geometry_prediction_target_check(x$target)
  if (!identical(x$table$measurement, x$target$measurement_ids) ||
      !identical(x$rank_by_fold, do.call(cbind, lapply(x$scores, function(score) score$table$rank)))) {
    .contract_error("Cross-fit measurements and fold ranks are inconsistent.")
  }
  for (field in c("gain", "inner_product", "prediction_norm_sq")) {
    expected <- as.numeric(do.call(cbind, lapply(x$scores, function(score) score$table[[field]])) %*% x$edges$weight)
    if (!isTRUE(all.equal(x$table[[field]], expected, tolerance = 1e-12))) {
      .contract_error("Cross-fit evidence must aggregate each frozen fold's gain and squared norm separately.")
    }
  }
  .check_signature(x$signature, .sha256_signature(.geometry_crossfit_semantic(x), "geometry-crossfit-sha256:"),
    "Cross-fit record identity is inconsistent.")
  x
}
