# Nested training-only selection (layer 5, internal orchestration) ----------

.geometry_candidates <- function(candidates) {
  if (!is.list(candidates) || !length(candidates)) .input_error("Declare a nonempty list of candidate tuples.")
  result <- lapply(candidates, function(candidate) {
    if (!is.list(candidate) || !all(c("basis", "rank", "penalty") %in% names(candidate)) ||
        anyDuplicated(names(candidate)) || any(!names(candidate) %in% c("basis", "rank", "penalty", "weights"))) {
      .input_error("Every candidate declares basis, rank, penalty and optional named model weights.")
    }
    basis <- .validate_model_basis(candidate$basis)
    pool <- .model_basis_pool(basis, candidate$weights)
    rank <- .check_count(candidate$rank, "candidate rank", min = 0L, max = .Machine$integer.max)
    .check_number(candidate$penalty, "candidate penalty", nonnegative = TRUE)
    fields <- list(basis = basis$signature, weights = pool$weights, rank = rank, penalty = unname(candidate$penalty))
    list(basis = basis, weights = pool$weights, rank = rank, penalty = unname(candidate$penalty),
      id = .sha256_signature(fields, "geometry-candidate-sha256:"))
  })
  ids <- vapply(result, `[[`, "", "id")
  if (anyDuplicated(ids)) .input_error("Candidate tuples must be distinct; duplicate recipes add no choice.")
  result <- result[order(ids, method = "radix")]
  stats::setNames(result, vapply(result, `[[`, "", "id"))
}

.geometry_selection_weights <- function(value, ids, label) {
  if (is.null(value)) return(stats::setNames(rep(1/length(ids), length(ids)), ids))
  if (!.is_finite_numeric(value) || length(value) != length(ids) ||
      !.is_strings(names(value), unique = TRUE) || !setequal(names(value), ids) ||
      any(value < 0) || abs(sum(value) - 1) > 1e-12) {
    .input_error(sprintf("`%s` must name nonnegative weights summing to one for the declared IDs.", label))
  }
  value[ids]
}

.geometry_select_indices <- function(values, candidates, scope, measurement_weights = NULL) {
  if (!.is_finite_matrix(values) || ncol(values) != length(candidates) || !nrow(values)) {
    .contract_error("Selection requires finite inner-validation gains for each candidate and measurement.")
  }
  rank <- vapply(candidates, `[[`, 0L, "rank")
  penalty <- vapply(candidates, `[[`, 0, "penalty")
  ids <- vapply(candidates, `[[`, "", "id")
  select <- function(gains) order(-gains, rank, -penalty, ids, method = "radix")[[1L]]
  if (identical(scope, "global")) {
    if (!.is_finite_numeric(measurement_weights) || length(measurement_weights) != nrow(values)) {
      .contract_error("Global selection requires declared measurement weights.")
    }
    rep.int(select(as.numeric(crossprod(measurement_weights, values))), nrow(values))
  } else if (identical(scope, "per_measurement")) {
    as.integer(vapply(seq_len(nrow(values)), function(i) select(values[i, ]), 0L))
  } else .input_error("Selection scope must be global or per_measurement.")
}

.geometry_inner_schedule <- function(train_plan, test_plan, inner, candidates, component) {
  .check_class(train_plan, "effect_geometry_plan", "train_plan", from = "plan_geometry()")
  .check_class(test_plan, "effect_geometry_plan", "test_plan", from = "plan_geometry()")
  train_support <- .geometry_training_support(train_plan)
  test_support <- .geometry_training_support(test_plan)
  .geometry_support_independence(train_support, test_support)
  .geometry_prediction_match(.geometry_prediction_target(train_plan, component),
    .geometry_prediction_target(test_plan, component))
  for (candidate in candidates) {
    .geometry_fit_prepare(train_plan, candidate$basis, component)
    .geometry_fit_prepare(test_plan, candidate$basis, component)
  }
  if (!is.list(inner) || !length(inner)) .input_error("Declare at least one feasible inner split and its weight.")
  weights <- vapply(inner, function(fold) {
    if (!is.list(fold) || !setequal(names(fold), c("train", "validate", "weight")) ||
        anyDuplicated(names(fold)) || !.is_number(fold$weight) || fold$weight <= 0) {
      .input_error("Each inner split declares train/validate partition sets and a positive weight.")
    }
    fold$weight
  }, 0)
  if (abs(sum(weights) - 1) > 1e-12) .input_error("Declared inner-fold weights must sum to one.")
  used <- names(train_support$endpoints)
  declared <- train_support$products
  folds <- lapply(seq_along(inner), function(i) {
    fold <- inner[[i]]
    if (!.is_strings(fold$train, unique = TRUE) || !.is_strings(fold$validate, unique = TRUE) ||
        length(intersect(fold$train, fold$validate)) ||
        !all(c(fold$train, fold$validate) %in% used)) {
      .capability_refusal("Inner training and validation must be disjoint subsets of actual outer-training partitions.",
        capability = "nested_geometry_selection", namespace = "predictive_geometry",
        reasons = "invalid_nested_partition_split", remedies = "Declare feasible inner folds entirely within outer training.")
    }
    subset_plan <- function(parts) {
      edges <- declared[declared$left %in% parts & declared$right %in% parts, , drop = FALSE]
      if (!nrow(edges)) {
        .capability_refusal("An inner subset contains no positive-weight product declared by outer training.",
          capability = "nested_geometry_selection", namespace = "predictive_geometry",
          reasons = "inner_split_has_no_declared_products", remedies = "Supply inner subsets with their own declared independent neural products.")
      }
      .geometry_prediction_pair_plan(train_plan, edges)
    }
    training <- subset_plan(fold$train); validation <- subset_plan(fold$validate)
    a <- .geometry_training_support(training); b <- .geometry_training_support(validation)
    .geometry_support_independence(a, b)
    list(training = training, validation = validation, training_support = a, validation_support = b,
      weight = weights[[i]], id = .sha256_signature(list(train = a$signature, validation = b$signature,
        weight = weights[[i]]), "geometry-inner-fold-sha256:"))
  })
  ids <- vapply(folds, `[[`, "", "id")
  if (anyDuplicated(ids)) .input_error("Duplicate inner splits add no validation information.")
  folds <- folds[order(ids, method = "radix")]
  list(folds = folds, training_support = train_support, evaluation_support = test_support,
    target = .geometry_prediction_target(train_plan, component))
}

.geometry_fit_candidate <- function(plan, candidate, component, row_block, upstream = list(), cache = NULL) {
  pool <- if (is.null(cache)) .model_basis_pool(candidate$basis, candidate$weights) else
    .geometry_cache_pool(cache, candidate$basis, candidate$weights)
  prepared <- .geometry_fit_prepare(plan, candidate$basis, component, upstream = upstream)
  parameters <- list(rank = candidate$rank, penalty = candidate$penalty, weights = candidate$weights,
    component = component, loss = "frobenius")
  .geometry_fit_execute(prepared, pool, parameters, "memory", NULL, row_block, FALSE, cache = cache)
}

.geometry_score_candidate <- function(fit, plan, row_block, cache = NULL) {
  if (is.null(cache)) return(score_geometry(fit, plan, row_block = row_block))
  prepared <- .geometry_score_prepare(fit, plan)
  .geometry_score_execute(fit, prepared, "memory", NULL, row_block,
    modes = FALSE, components = FALSE, cache = cache)
}

.geometry_selection_caches <- function(train, evaluate, candidates, reuse) {
  if (!reuse) return(list(training = NULL, evaluation = NULL))
  train_bytes <- .geometry_cache_bytes(train, candidates, "training")
  eval_bytes <- .geometry_cache_bytes(evaluate, candidates, "evaluation")
  list(training = .new_geometry_training_cache(candidates, "training", partner_bytes = eval_bytes),
    evaluation = .new_geometry_training_cache(candidates, "evaluation", partner_bytes = train_bytes))
}

.geometry_selection_cache_summary <- function(caches) {
  lapply(caches, function(cache) if (is.null(cache)) NULL else .geometry_cache_summary(cache))
}

.geometry_selection_semantic <- function(x) {
  value <- unclass(x[setdiff(names(x), c("signature", "candidates"))])
  value$candidates <- vapply(x$candidates, `[[`, "", "id")
  value
}

.validate_geometry_selection <- function(x) {
  fields <- c("candidates", "scope", "measurement_weights", "fold_weights", "inner_gains",
    "mean_inner_gain", "selected", "dependencies", "target", "tie_rule", "signature")
  if (!is.list(x) || !identical(names(x), fields) ||
      !identical(dim(x$inner_gains), c(length(x$target$measurement_ids), length(x$candidates), length(x$fold_weights)))) {
    .contract_error("A selection record requires all declared candidates, inner folds and target measurements.")
  }
  .geometry_prediction_target_check(x$target)
  lapply(x$dependencies, .validate_geometry_training_support)
  for (candidate in x$candidates) {
    rebuilt <- .geometry_candidates(list(candidate[setdiff(names(candidate), "id")]))[[1]]
    if (!identical(candidate, rebuilt)) .contract_error("A selected candidate recipe was changed.")
  }
  expected <- matrix(0, length(x$target$measurement_ids), length(x$candidates))
  for (i in seq_along(x$fold_weights)) expected <- expected + x$fold_weights[[i]] *
    matrix(x$inner_gains[, , i], nrow(expected), ncol(expected))
  selected <- .geometry_select_indices(expected, x$candidates, x$scope, x$measurement_weights)
  if (!isTRUE(all.equal(x$mean_inner_gain, expected, tolerance = 0, check.attributes = FALSE)) ||
      !identical(x$selected, vapply(x$candidates, `[[`, "", "id")[selected])) {
    .contract_error("Selected recipes must follow only the declared inner gains and deterministic tie rule.")
  }
  .check_signature(x$signature, .sha256_signature(.geometry_selection_semantic(x), "geometry-selection-sha256:"),
    "Nested selection identity is inconsistent.")
  x
}

.geometry_nested <- function(train_plan, test_plan, candidates, inner,
                              scope = c("global", "per_measurement"), measurement_weights = NULL,
                              component = "total", row_block = 128L, reuse = TRUE) {
  scope <- match.arg(scope)
  .check_flag(reuse, "reuse")
  row_block <- .check_count(row_block, "row_block", max = .Machine$integer.max)
  candidates <- .geometry_candidates(candidates)
  schedule <- .geometry_inner_schedule(train_plan, test_plan, inner, candidates, component)
  ids <- schedule$target$measurement_ids
  if (scope == "global") measurement_weights <- .geometry_selection_weights(measurement_weights, ids, "measurement_weights") else {
    if (!is.null(measurement_weights)) .input_error("Per-measurement selection does not pool measurements; omit measurement_weights.")
  }
  n <- length(ids); nc <- length(candidates); nf <- length(schedule$folds)
  # This initial orchestration retains memory results. Its explicit reserve
  # covers candidate/fold summaries and every possible selected final factor.
  factor_widths <- vapply(candidates, function(x) {
    k <- min(x$rank, x$basis$dimension)
    x$basis$dimension * k + k + 2 * x$basis$dimension
  }, 0)
  retained_bound <- 2 * (as.double(utils::object.size(schedule)) + as.double(utils::object.size(candidates)) +
    512 * as.double(n) * nc * (nf + 1) + 8 * n * sum(factor_widths))
  budget_plan <- function(plan) {
    remaining <- if (is.null(plan$compute$workspace_bytes)) NULL else plan$compute$workspace_bytes - retained_bound
    if (!is.null(remaining) && remaining <= 0) {
      .capability_refusal("Nested selection records exceed the declared workspace.",
        capability = "predictive_selection_workspace", namespace = "predictive_geometry",
        reasons = "selection_workspace_budget_exceeded", remedies = "Increase workspace or reduce the declared candidate/fold/measurement set.")
    }
    .geometry_prediction_pair_plan(plan, as.data.frame(plan$pairing),
      compute_policy(block_features = plan$compute$block_features, workspace_bytes = remaining))
  }
  train_plan <- budget_plan(train_plan); test_plan <- budget_plan(test_plan)
  gains <- array(0, c(n, nc, nf), dimnames = list(ids, names(candidates),
    vapply(schedule$folds, `[[`, "", "id")))
  dependencies <- list()
  cache_statistics <- vector("list", nf + 1L)
  for (i in seq_len(nf)) {
    fold <- schedule$folds[[i]]
    train <- budget_plan(fold$training); validate <- budget_plan(fold$validation)
    caches <- .geometry_selection_caches(train, validate, candidates, reuse)
    dependencies <- c(dependencies, list(fold$training_support, fold$validation_support))
    for (j in seq_len(nc)) {
      fit <- .geometry_fit_candidate(train, candidates[[j]], component, row_block, cache = caches$training)
      evidence <- .geometry_score_candidate(fit, validate, row_block, cache = caches$evaluation)
      gains[, j, i] <- evidence$table$gain
    }
    cache_statistics[[i]] <- .geometry_selection_cache_summary(caches)
  }
  fold_weights <- vapply(schedule$folds, `[[`, 0, "weight")
  mean_gain <- matrix(0, n, nc)
  for (i in seq_len(nf)) mean_gain <- mean_gain + fold_weights[[i]] * matrix(gains[, , i], n, nc)
  chosen <- .geometry_select_indices(mean_gain, candidates, scope, measurement_weights)
  selection <- list(candidates = candidates, scope = scope, measurement_weights = measurement_weights,
    fold_weights = fold_weights, inner_gains = gains, mean_inner_gain = mean_gain,
    selected = vapply(candidates, `[[`, "", "id")[chosen], dependencies = dependencies,
    target = schedule$target, tie_rule = "lower rank, larger penalty, canonical candidate ID")
  selection$signature <- .sha256_signature(.geometry_selection_semantic(selection), "geometry-selection-sha256:")
  .validate_geometry_selection(selection)
  # The selection record is complete before any outer evaluation values are
  # read. Refit chosen recipes using outer training, retaining inner ancestry.
  chosen_ids <- sort(unique(selection$selected), method = "radix")
  caches <- .geometry_selection_caches(train_plan, test_plan, candidates[chosen_ids], reuse)
  fits <- stats::setNames(lapply(chosen_ids, function(id) .geometry_fit_candidate(train_plan,
    candidates[[id]], component, row_block, upstream = dependencies, cache = caches$training)), chosen_ids)
  scores <- lapply(fits, function(fit) .geometry_score_candidate(fit, test_plan, row_block, caches$evaluation))
  cache_statistics[[nf + 1L]] <- .geometry_selection_cache_summary(caches)
  table <- scores[[1]]$table
  for (i in seq_len(n)) table[i, ] <- scores[[selection$selected[[i]]]]$table[i, ]
  result <- list(table = table, selection = selection, fits = fits, scores = scores,
    estimand = "Predictive gain after inner selection and refitting on the declared outer-training data.",
    retained_workspace_bound = retained_bound, cache_statistics = cache_statistics, complete = TRUE)
  result$signature <- .sha256_signature(.geometry_nested_semantic(result), "geometry-nested-sha256:")
  .validate_geometry_nested(structure(result, class = "effect_geometry_nested"))
}

.geometry_nested_semantic <- function(x) {
  list(table = x$table, selection = x$selection$signature,
    fits = vapply(x$fits, `[[`, "", "signature"), scores = vapply(x$scores, `[[`, "", "signature"),
    estimand = x$estimand, retained_workspace_bound = x$retained_workspace_bound,
    cache_statistics = x$cache_statistics, complete = x$complete)
}

.validate_geometry_nested <- function(x) {
  fields <- c("table", "selection", "fits", "scores", "estimand", "retained_workspace_bound", "cache_statistics", "complete", "signature")
  if (!.sealed_fields(x, "effect_geometry_nested", fields) || !isTRUE(x$complete)) {
    .contract_error("Nested predictive evidence requires a complete selection, refit and outer score.")
  }
  .validate_geometry_selection(x$selection)
  chosen <- sort(unique(x$selection$selected), method = "radix")
  if (!identical(names(x$fits), chosen) || !identical(names(x$scores), chosen) ||
      !identical(x$table$measurement, x$selection$target$measurement_ids)) {
    .contract_error("Nested evidence must retain exactly the chosen recipes on the declared measurements.")
  }
  for (id in chosen) {
    fit <- .validate_geometry_fit(x$fits[[id]])
    score <- .validate_geometry_score(x$scores[[id]])
    candidate <- x$selection$candidates[[id]]
    if (!identical(fit$prediction_id, score$prediction_id) ||
        !identical(fit$model$signature, candidate$basis$signature) ||
        !identical(fit$parameters$rank, candidate$rank) ||
        !identical(fit$parameters$penalty, candidate$penalty) ||
        !identical(fit$parameters$weights, candidate$weights) ||
        !identical(fit$training$upstream, x$selection$dependencies)) {
      .contract_error("An outer fit must use the selected recipe and retain every inner training/validation dependency.")
    }
  }
  for (i in seq_len(nrow(x$table))) {
    expected <- x$scores[[x$selection$selected[[i]]]]$table[i, , drop = FALSE]
    if (!isTRUE(all.equal(x$table[i, , drop = FALSE], expected, tolerance = 0))) {
      .contract_error("Outer score rows must follow the frozen per-measurement selection.")
    }
  }
  .check_signature(x$signature, .sha256_signature(.geometry_nested_semantic(x), "geometry-nested-sha256:"),
    "Nested evidence identity is inconsistent.")
  x
}
