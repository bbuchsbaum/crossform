# Frozen predictive values and bounded readers (layer 4) --------------------
# No fitting/scoring calls live here. The record is shared by those readers
# without introducing a fit/score cycle or raw-geometry capabilities.

.geometry_prediction_layout <- function(dimension, support_dimension, budget) {
  widths <- c(factor = dimension * budget, amplitudes = budget,
    signed = dimension, shifted = support_dimension)
  ends <- cumsum(widths)
  columns <- Map(function(end, width) if (width == 0) integer() else
    seq.int(end - width + 1L, end), ends, widths)
  list(dimension = as.integer(dimension), support_dimension = as.integer(support_dimension),
    budget = as.integer(budget), width = as.integer(sum(widths)), columns = columns)
}

.geometry_prediction_target_check <- function(target) {
  fields <- c("effect_space", "domain", "metric_signature", "normalization", "component",
    "generalizes_over", "measurement_ids", "measurement_signatures", "index")
  if (!is.list(target) || !identical(names(target), fields) ||
      !.is_strings(target$measurement_ids, unique = TRUE) || !length(target$measurement_ids) ||
      !identical(as.character(target$index), target$measurement_ids) ||
      length(target$measurement_signatures) != length(target$measurement_ids) ||
      !.strong_sha256(target$metric_signature) ||
      !all(vapply(sub("^prediction-node-", "", target$measurement_signatures), .strong_sha256, FALSE))) {
    .contract_error("Prediction target identities are missing or inconsistent.")
  }
  .validate_effect_space(target$effect_space)
  .validate_domain_reference(target$domain)
  invisible(target)
}

.geometry_prediction_target_semantic <- function(target) {
  ordered <- order(target$measurement_ids, method = "radix")
  result <- target
  result$measurement_ids <- target$measurement_ids[ordered]
  result$measurement_signatures <- target$measurement_signatures[ordered]
  result$index <- NULL
  result
}

.geometry_prediction_interpretation <- function(basis, pool, penalty) {
  full_span <- pool$dimension == basis$n_conditions - 1L
  isotropic <- max(pool$values) - min(pool$values) <= pool$tolerance * max(pool$values)
  list(layer = "learned_prediction", unbiased = FALSE,
    preference = if (penalty == 0) "span_only" else if (isotropic) "isotropic_shrinkage" else "kernel_regularized",
    full_centered_span = full_span,
    generic_full_span_baseline = full_span && penalty == 0,
    loss = "frobenius_geometry", score_units = "squared geometry units")
}

.geometry_fit_fields <- c("model", "pool", "parameters", "target", "training",
  "layout", "values", "row_signatures", "diagnostics", "interpretation",
  "source", "receipt", "workspace", "prediction_id", "complete", "signature")

.geometry_fit_prediction_id <- function(x) {
  .sha256_signature(list(contract = "predictive-geometry-v1", model = x$model$signature,
    pool = x$pool$signature, parameters = x$parameters,
    target = .geometry_prediction_target_semantic(x$target), training = x$training$signature),
    "geometry-prediction-sha256:")
}

.geometry_fit_semantic <- function(x) {
  list(prediction_id = x$prediction_id, layout = x$layout,
    row_signatures = x$row_signatures, diagnostics = x$diagnostics,
    interpretation = x$interpretation, complete = x$complete,
    target = x$target, receipt = x$receipt, workspace = x$workspace,
    storage = list(dim = x$values$dim, representation = x$values$representation,
      manifest = x$values$manifest, path = x$values$path),
    source = if (is.null(x$source)) NULL else x$source$contract_signature)
}

.new_geometry_fit <- function(model, pool, parameters, target, training, layout,
                              values, row_signatures, diagnostics, source, receipt, workspace) {
  value <- list(model = model, pool = pool, parameters = parameters, target = target,
    training = training, layout = layout, values = values, row_signatures = row_signatures,
    diagnostics = diagnostics,
    interpretation = .geometry_prediction_interpretation(model, pool, parameters$penalty),
    source = source, receipt = receipt, workspace = workspace, prediction_id = NULL, complete = TRUE)
  value$prediction_id <- .geometry_fit_prediction_id(value)
  value$signature <- .sha256_signature(.geometry_fit_semantic(value), "geometry-fit-sha256:")
  structure(value, class = "effect_geometry_fit")
}

.validate_geometry_fit <- function(x) {
  if (!.sealed_fields(x, "effect_geometry_fit", .geometry_fit_fields) || !isTRUE(x$complete)) {
    .contract_error("A geometry fit must be a complete, canonical frozen prediction.")
  }
  model <- .validate_model_basis(x$model)
  .geometry_prediction_target_check(x$target)
  .validate_geometry_training_support(x$training)
  .validate_execution_receipt(x$receipt)
  .validate_prediction_workspace(x$workspace, x$receipt)
  pars <- x$parameters
  if (!is.list(pars) || !identical(names(pars), c("rank", "penalty", "weights", "component", "loss")) ||
      !.is_count(pars$rank, min = 0L) || pars$rank > .Machine$integer.max ||
      !.is_number(pars$penalty) || pars$penalty < 0 ||
      !identical(pars$component, x$target$component) || !identical(pars$loss, "frobenius")) {
    .contract_error("Fitted prediction parameters are invalid or inconsistent.")
  }
  pool <- .model_basis_pool(model, pars$weights)
  layout <- .geometry_prediction_layout(model$dimension, pool$dimension, min(pars$rank, pool$dimension))
  if (!identical(x$pool, pool) || !identical(x$layout, layout) ||
      !identical(x$interpretation, .geometry_prediction_interpretation(model, pool, pars$penalty))) {
    .contract_error("Fitted model support, layout or interpretation is inconsistent.")
  }
  .validate_geometry_store(x$values, "prediction values", "rectangular", probe = FALSE)
  n <- length(x$target$measurement_ids)
  if (!identical(x$values$dim, c(as.integer(n), layout$width)) ||
      length(x$row_signatures) != n ||
      !all(vapply(x$row_signatures, .strong_sha256, FALSE)) ||
      !is.data.frame(x$diagnostics) || nrow(x$diagnostics) != n ||
      !identical(x$diagnostics$measurement, x$target$measurement_ids)) {
    .contract_error("Fitted measurements and stored prediction rows disagree.")
  }
  if (!is.null(x$source)) {
    .validate_effect_form(x$source, probe = FALSE)
    if (!identical(x$source$effect_space, model$effect_space) ||
        !identical(x$source$index, x$target$index)) {
      .contract_error("The retained signed form does not match the fitted coordinates or measurements.")
    }
  }
  .check_signature(x$prediction_id, .geometry_fit_prediction_id(x),
    "Fitted prediction identity is inconsistent with its training and target.")
  .check_signature(x$signature, .sha256_signature(.geometry_fit_semantic(x), "geometry-fit-sha256:"),
    "Fitted prediction record identity is inconsistent.")
  x
}

.geometry_prediction_row <- function(value, layout) {
  list(factor = matrix(value[layout$columns$factor], layout$dimension, layout$budget),
    amplitudes = value[layout$columns$amplitudes], signed_spectrum = value[layout$columns$signed],
    shifted_spectrum = value[layout$columns$shifted])
}

.geometry_prediction_rows <- function(x, rows, validate = TRUE) {
  if (validate) x <- .validate_geometry_fit(x)
  if (!is.numeric(rows) || !length(rows) || anyNA(rows) || any(rows %% 1 != 0) ||
      any(rows < 1) || any(rows > x$values$dim[[1L]])) {
    .input_error("Prediction rows must identify existing measurements.")
  }
  values <- .read_geometry_store(x$values, rows)
  signatures <- vapply(seq_len(nrow(values)), function(i) .sha256_signature(as.numeric(values[i, ])), "")
  if (!identical(unname(signatures), unname(x$row_signatures[rows]))) {
    .contract_error("Stored prediction values no longer match their frozen row identities.")
  }
  lapply(seq_len(nrow(values)), function(i) {
    value <- .geometry_prediction_row(values[i, ], x$layout)
    amp <- value$amplitudes
    shifted <- value$shifted_spectrum
    expected_amp <- pmax(shifted[seq_len(x$layout$budget)], 0)
    norm_sq <- sum(amp^2)
    penalty_trace <- sum(rowSums(crossprod(x$pool$vectors, value$factor)^2) / x$pool$values)
    diag_row <- x$diagnostics[rows[[i]], , drop = FALSE]
    scale <- max(1, abs(shifted), abs(value$signed_spectrum))
    error <- max(c(abs(amp - expected_amp),
      abs(crossprod(value$factor) - diag(amp, length(amp))),
      abs(value$factor - x$pool$vectors %*% crossprod(x$pool$vectors, value$factor)), 0))
    if (error > 1e-8 * scale || diag_row$rank != sum(amp > 0) ||
        abs(diag_row$prediction_norm_sq - norm_sq) > 1e-8 * max(1, norm_sq) ||
        abs(diag_row$penalty_trace - penalty_trace) > 1e-8 * max(1, penalty_trace) ||
        abs(diag_row$objective - .5 * (sum(value$signed_spectrum^2) - norm_sq)) >
          1e-8 * max(1, sum(value$signed_spectrum^2), norm_sq)) {
      .contract_error("Stored prediction factors or derived fit diagnostics are inconsistent.")
    }
    value
  })
}
