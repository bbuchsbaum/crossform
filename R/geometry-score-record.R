# Frozen predictive evidence records (layer 4) -----------------------------

.geometry_score_fields <- c("prediction_id", "model_signature", "parameters", "target",
  "evaluation_target", "training", "evaluation", "validation", "table", "components",
  "mode_values", "mode_budget", "row_signatures", "row_block", "receipt", "workspace", "complete", "signature")

.geometry_score_semantic <- function(x) {
  result <- unclass(x[setdiff(.geometry_score_fields, c("signature", "mode_values"))])
  result$storage <- if (is.null(x$mode_values)) NULL else list(dim = x$mode_values$dim,
    representation = x$mode_values$representation, manifest = x$mode_values$manifest,
    path = x$mode_values$path)
  result
}

.validate_geometry_score <- function(x) {
  if (!.sealed_fields(x, "effect_geometry_score", .geometry_score_fields) || !isTRUE(x$complete)) {
    .contract_error("Predictive evidence must be a complete, canonical score record.")
  }
  .geometry_prediction_target_check(x$target)
  .geometry_prediction_target_check(x$evaluation_target)
  .geometry_support_independence(x$training, x$evaluation)
  .validate_execution_receipt(x$receipt)
  .validate_prediction_workspace(x$workspace, x$receipt)
  mapping <- match(x$target$measurement_ids, x$evaluation_target$measurement_ids)
  if (!identical(.geometry_prediction_target_semantic(x$target),
      .geometry_prediction_target_semantic(x$evaluation_target)) ||
      !identical(x$validation$mapping, mapping) || !isTRUE(x$validation$independent) ||
      !identical(x$validation$training, x$training$signature) ||
      !identical(x$validation$evaluation, x$evaluation$signature) ||
      !identical(x$validation$origin_scope, x$training$scope_signature)) {
    .contract_error("Predictive validation does not match its target and observation dependencies.")
  }
  tab <- x$table
  fields <- c("measurement", "gain", "rank", "inner_product", "prediction_norm_sq")
  if (!is.data.frame(tab) || !identical(names(tab), fields) ||
      !identical(tab$measurement, x$target$measurement_ids) ||
      !.is_finite_matrix(as.matrix(tab[-1L])) || any(tab$rank < 0 | tab$rank %% 1 != 0) ||
      any(tab$prediction_norm_sq < 0) ||
      any(abs(tab$gain - (2 * tab$inner_product - tab$prediction_norm_sq)) >
        1e-10 * pmax(1, abs(tab$gain), abs(tab$inner_product), tab$prediction_norm_sq))) {
    .contract_error("Predictive gain, cost or measurement rows are inconsistent.")
  }
  if (!is.null(x$components)) {
    cmp <- x$components
    if (!identical(x$parameters$component, "total") || !is.data.frame(cmp) ||
        !identical(names(cmp), c("measurement", "coherent_inner_product", "configuration_inner_product")) ||
        !identical(cmp$measurement, tab$measurement) || !.is_finite_matrix(as.matrix(cmp[-1L])) ||
        any(abs(cmp$coherent_inner_product + cmp$configuration_inner_product - tab$inner_product) >
          1e-10 * pmax(1, abs(tab$inner_product), abs(cmp$coherent_inner_product), abs(cmp$configuration_inner_product)))) {
      .contract_error("Frozen component evidence must add before the prediction cost is subtracted once.")
    }
  }
  if (!.is_count(x$mode_budget, min = 0L) || !.is_count(x$row_block) ||
      any(tab$rank > x$mode_budget)) .contract_error("Predictive mode bounds are inconsistent.")
  if (!is.null(x$mode_values)) {
    .validate_geometry_store(x$mode_values, "score modes", "rectangular", probe = FALSE)
    if (!identical(x$mode_values$dim, c(as.integer(nrow(tab)), as.integer(4L * x$mode_budget))) ||
        length(x$row_signatures) != nrow(tab) ||
        !all(vapply(x$row_signatures, .strong_sha256, FALSE))) {
      .contract_error("Predictive mode storage does not match its declared rows.")
    }
  } else if (length(x$row_signatures)) .contract_error("Mode identities require a mode store.")
  .check_signature(x$signature, .sha256_signature(.geometry_score_semantic(x), "geometry-score-sha256:"),
    "Predictive score record identity is inconsistent.")
  x
}

.geometry_score_mode_rows <- function(x, rows) {
  if (is.null(x$mode_values) || !length(rows)) {
    if (is.null(x$mode_values) && x$mode_budget > 0L) .input_error("Mode evidence was not retained; request modes = TRUE when scoring.")
    return(data.frame(measurement = character(), mode = integer(), group = integer(),
      amplitude = numeric(), evidence = numeric(), prediction_cost = numeric(), gain = numeric()))
  }
  blocks <- split(rows, ceiling(seq_along(rows) / x$row_block))
  result <- lapply(blocks, function(at) {
    values <- .read_geometry_store(x$mode_values, at)
    hashes <- vapply(seq_len(nrow(values)), function(i) .sha256_signature(as.numeric(values[i, ])), "")
    if (!identical(unname(hashes), unname(x$row_signatures[at]))) {
      .contract_error("Stored score modes no longer match their frozen row identities.")
    }
    do.call(rbind, lapply(seq_along(at), function(i) {
      k <- x$mode_budget
      amp <- values[i, seq_len(k)]
      evidence <- values[i, k + seq_len(k)]
      gain <- values[i, 2L * k + seq_len(k)]
      group <- values[i, 3L * k + seq_len(k)]
      active <- which(amp > 0)
      expected <- 2 * amp * evidence - amp^2
      tab <- x$table[at[[i]], ]
      if (any(!is.finite(c(amp, evidence, gain, group))) || any(amp < 0) ||
          length(active) != tab$rank || any(group[active] < 1 | group[active] %% 1 != 0) ||
          max(c(abs(gain - expected), abs(sum(gain) - tab$gain),
            abs(sum(amp^2) - tab$prediction_norm_sq), 0)) >
              1e-10 * max(1, abs(expected), abs(tab$gain), tab$prediction_norm_sq)) {
        .contract_error("Predictive modes disagree with frozen amplitudes, costs or total gain.")
      }
      data.frame(measurement = rep(tab$measurement, length(active)), mode = active,
        group = as.integer(group[active]), amplitude = amp[active], evidence = evidence[active],
        prediction_cost = amp[active]^2, gain = gain[active])
    }))
  })
  do.call(rbind, result)
}
