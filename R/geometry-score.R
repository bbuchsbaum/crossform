# Independent readout of frozen representational predictions (layer 5) -----

.geometry_score_row <- function(prediction, S, tolerance) {
  S <- .geometry_fit_symmetric(S)
  amp <- prediction$amplitudes
  Z <- prediction$factor
  weighted_evidence <- colSums(Z * (S %*% Z))
  active <- which(amp > 0)
  evidence <- numeric(length(amp))
  evidence[active] <- weighted_evidence[active] / amp[active]
  gains <- 2 * weighted_evidence - amp^2
  groups <- integer(length(amp))
  groups[active] <- .geometry_fit_mode_groups(amp[active],
    max(abs(prediction$shifted_spectrum)), tolerance)
  result <- list(inner_product = sum(weighted_evidence), prediction_norm_sq = sum(amp^2),
    gain = sum(gains), rank = as.integer(length(active)),
    amplitudes = amp, evidence = evidence, mode_gain = gains, groups = groups)
  if (any(!is.finite(unlist(result)))) .input_error("Predictive gain overflows; rescale geometry before fitting.")
  result
}

.geometry_score_from_form <- function(fit, prepared, source, storage, storage_path,
                                      row_block, modes, components) {
  n <- length(fit$target$measurement_ids)
  .validate_effect_form(source, probe = FALSE)
  .require_effect_form_component(source, fit$parameters$component)
  k <- fit$layout$budget
  table <- data.frame(measurement = fit$target$measurement_ids, gain = numeric(n),
    rank = integer(n), inner_product = numeric(n), prediction_norm_sq = numeric(n))
  cmp <- if (components) data.frame(measurement = table$measurement,
    coherent_inner_product = numeric(n), configuration_inner_product = numeric(n)) else NULL
  keep_modes <- modes && k > 0L
  hashes <- if (keep_modes) character(n) else character()
  output <- if (keep_modes && storage == "memory") matrix(0, n, 4L * k) else NULL
  mode_store <- if (keep_modes && storage == "block") .file_geometry_store(
    file.path(storage_path, "modes.bin"), c(n, 4L * k), create = TRUE,
    codec = "rectangular") else NULL
  for (start in seq.int(1L, n, by = row_block)) {
    rows <- seq.int(start, min(n, start + row_block - 1L))
    mapped <- prepared$mapping[rows]
    predictions <- .geometry_prediction_rows(fit, rows, validate = FALSE)
    packed <- .geometry_component_validated(source, fit$parameters$component, rows = mapped)
    coherent <- if (components) .geometry_component_validated(source, "coherent", rows = mapped) else NULL
    configuration <- if (components) .geometry_component_validated(source, "configuration", rows = mapped) else NULL
    mode_values <- if (keep_modes) matrix(0, length(rows), 4L * k) else NULL
    for (i in seq_along(rows)) {
      prediction <- predictions[[i]]
      result <- .geometry_score_row(prediction, .unsvec_symmetric(packed[i, ], fit$model$dimension),
        fit$model$tolerance)
      for (field in names(table)[-1L]) table[rows[[i]], field] <- result[[field]]
      if (components) {
        for (component in c("coherent", "configuration")) {
          small <- .unsvec_symmetric(if (component == "coherent") coherent[i, ] else configuration[i, ],
            fit$model$dimension)
          cmp[rows[[i]], paste0(component, "_inner_product")] <-
            sum(prediction$factor * (small %*% prediction$factor))
        }
      }
      if (keep_modes) {
        mode_values[i, ] <- c(result$amplitudes, result$evidence, result$mode_gain, result$groups)
        hashes[rows[[i]]] <- .sha256_signature(as.numeric(mode_values[i, ]))
      }
    }
    if (keep_modes) {
      if (storage == "memory") output[rows, ] <- mode_values else
        .write_geometry_tile(mode_store, rows, seq_len(4L * k), mode_values)
    }
  }
  if (keep_modes && storage == "memory") mode_store <- .memory_geometry_store(output, "rectangular")
  value <- list(prediction_id = fit$prediction_id, model_signature = fit$model$signature,
    parameters = fit$parameters, target = fit$target, evaluation_target = prepared$prepared$target,
    training = fit$training, evaluation = prepared$prepared$support, validation = prepared$validation,
    table = table, components = cmp, mode_values = mode_store, mode_budget = k,
    row_signatures = hashes, row_block = row_block, receipt = source$receipt,
    workspace = prepared$prepared$workspace, complete = TRUE)
  value$signature <- .sha256_signature(.geometry_score_semantic(value), "geometry-score-sha256:")
  .validate_geometry_score(structure(value, class = "effect_geometry_score"))
}

#' Score a frozen representational prediction on independent geometry
#'
#' Measure gain over the zero-geometry predictor: twice the signed held-out
#' inner product minus the fitted prediction's squared Frobenius norm.
#' Positive gain improves squared geometry prediction error in expectation
#' under the declared independence assumptions. Gain and mode evidence may
#' be negative; no fitting or selection occurs during scoring.
#'
#' @param fit A frozen [fit_geometry()] result.
#' @param plan Geometry plan on disjoint observation origins from the same
#'   declared parent manifest, with the same effects, fixed neural metric,
#'   generalization axis and spatial measurements. Measurement IDs are matched
#'   explicitly, so row order may differ. Both neural cross-product independence
#'   and fitting/evaluation independence are required.
#' @param storage `"memory"` or `"block"` for evaluation forms and mode evidence.
#' @param storage_path A new durable directory for block storage. The complete
#'   score is saved as `score.rds`; existing directories are never overwritten.
#' @param row_block Maximum number of measurement rows processed together.
#' @param modes Retain signed per-mode evidence and gain. Tied modes have
#'   arbitrary individual orientations; use the grouped view for invariant
#'   evidence over a fully retained tied eigenspace.
#' @param components For a total-geometry fit, retain coherent and configuration
#'   inner products read with the same frozen prediction. They add before the
#'   prediction cost is subtracted once. Separately fitted latent components
#'   do not in general conserve the signed decomposition.
#' @return An `effect_geometry_score` with per-measurement gain, rank, inner
#'   product and prediction norm, actual training/evaluation dependencies,
#'   validation assumptions and an execution receipt. `as.data.frame()` returns
#'   the summary; `view = "modes"` or `"groups"` reads retained details. Gain has
#'   squared geometry units. There is no explained fraction, generic standard
#'   error or raw crossvalidated geometry inheritance.
#' @details Conditional on all training and selection data, the mean gain is
#'   `||G*||^2 - ||G* - F||^2` when the independent test form is unbiased for
#'   `G*`. In particular, under a zero signal the expected gain is `-||F||^2`,
#'   not zero. Rank, penalty, model weights and preprocessing must be fixed
#'   without consulting these evaluation outcomes. Same-condition partition
#'   prediction does not establish generalization to new conditions.
#'
#'   The summary view retains measurement ID, gain, effective rank, held-out
#'   inner product and prediction cost. The mode view separates fitted
#'   amplitude, signed test evidence and gain. The group view sums these over
#'   a fully retained tied eigenspace. `rows` selects measurement rows by
#'   position, preserving order and repetitions; empty selections retain the
#'   corresponding data-frame schema. Mode/group views require `modes = TRUE`
#'   at scoring unless the rank budget is zero.
#'
#'   See `vignette("predictive-geometry")` for a complete fixed-split workflow.
#' @examples
#' d <- abstract_domain(2, id = "score-example")
#' B <- matrix(c(1, -1, 0, 0, 0, 0, 1, -1), 4,
#'   dimnames = list(c("face", "body", "house", "tool"), NULL))
#' origins <- list(id = "four-runs", partitions = list(a = "run1", b = "run2",
#'   c = "run3", d = "run4"), independence = "independent",
#'   assumption = "Independent runs with externally fixed preprocessing.")
#' rel <- relation(list(a = B, b = B, c = B, d = B), domain = d,
#'   provenance = list(observation_origins = origins))
#' at <- compile_frame(whole_brain(normalization = "none"), d)
#' train <- plan_geometry(rel, at, pairing("a", "b",
#'   independence = "independent", generalizes_over = "run"))
#' test <- plan_geometry(rel, at, pairing("c", "d",
#'   independence = "independent", generalizes_over = "run"))
#' basis <- model_basis(features = B, conditions = rel$effect_space, normalize = "trace")
#' fit <- fit_geometry(train, basis, rank = 2, penalty = 0.1)
#' evidence <- score_geometry(fit, test, modes = TRUE)
#' as.data.frame(evidence)
#' as.data.frame(evidence, view = "groups")
#' @seealso [fit_geometry()], [model_basis()]
#' @export
score_geometry <- function(fit, plan, storage = c("memory", "block"),
                            storage_path = NULL, row_block = 128L,
                            modes = FALSE, components = FALSE) {
  storage <- match.arg(storage)
  row_block <- .check_count(row_block, "row_block", max = .Machine$integer.max)
  .check_flag(modes, "modes"); .check_flag(components, "components")
  fit <- .validate_geometry_fit(fit)
  if (components && fit$parameters$component != "total") {
    .input_error("Component decomposition is a readout of a frozen total-geometry prediction.")
  }
  prepared <- .geometry_score_prepare(fit, plan)
  .geometry_score_execute(fit, prepared, storage, storage_path, row_block, modes, components)
}

.geometry_score_execute <- function(fit, prepared, storage, storage_path, row_block,
                                    modes, components, cache = NULL) {
  cache_bytes <- 0; frozen_bound <- NULL
  budget <- fit$layout$budget
  if (!is.null(cache)) {
    .geometry_cache_check(cache, "evaluation")
    if (storage != "memory") .input_error("The local evaluation cache currently retains memory geometry only.")
    cache_bytes <- .geometry_cache_bind(cache, prepared$prepared, row_block)
    budget <- fit$model$dimension
    n <- length(fit$target$measurement_ids); r <- as.double(fit$model$dimension)
    frozen_bound <- 512 * n + 8 * n * (r^2 + 3 * r)
  }
  prepared$prepared <- .geometry_prediction_workspace(prepared$prepared, budget,
    storage, row_block, "score", fit, modes, components, cache_bytes = cache_bytes,
    frozen_prediction_bytes = frozen_bound)
  row_block <- prepared$prepared$workspace$row_block
  created <- FALSE; success <- FALSE
  if (storage == "block") {
    .check_string(storage_path, "storage_path", what = "a new durable score directory")
    if (file.exists(storage_path)) .input_error("Refusing to overwrite an existing score directory.")
    if (!dir.create(storage_path, recursive = TRUE)) .input_error("Cannot create the score directory.")
    storage_path <- normalizePath(storage_path, mustWork = TRUE); created <- TRUE
  } else if (!is.null(storage_path)) .input_error("`storage_path` is used only with block storage.")
  on.exit(if (created && !success) unlink(storage_path, recursive = TRUE), add = TRUE)
  source_path <- if (storage == "block") file.path(storage_path, "evaluation") else NULL
  source <- if (is.null(cache)) .geometry_fit_materialize(prepared$prepared, storage, source_path) else
    .geometry_cache_source(cache, prepared$prepared)
  result <- .geometry_score_from_form(fit, prepared, source, storage, storage_path, row_block, modes, components)
  if (storage == "block") {
    unlink(source_path, recursive = TRUE)
    pending <- file.path(storage_path, "score.pending.rds")
    saveRDS(result, pending)
    if (!file.rename(pending, file.path(storage_path, "score.rds"))) .input_error("Cannot finalize the score record.")
  }
  success <- TRUE
  result
}

#' @export
as.data.frame.effect_geometry_score <- function(x, row.names = NULL, optional = FALSE,
                                               ..., view = c("summary", "modes", "groups"), rows = NULL) {
  .check_no_extra_arguments("as.data.frame.effect_geometry_score", ...)
  x <- .validate_geometry_score(x)
  view <- match.arg(view)
  if (is.null(rows)) rows <- seq_len(nrow(x$table))
  if (!is.numeric(rows) || anyNA(rows) || any(rows %% 1 != 0) || any(rows < 1) || any(rows > nrow(x$table))) {
    .input_error("Score rows must identify existing measurements.")
  }
  value <- if (view == "summary") x$table[rows, , drop = FALSE] else .geometry_score_mode_rows(x, rows)
  if (view == "groups") {
    if (nrow(value)) {
      # Preserve requested row order, including repeated measurements. Each
      # new mode sequence identifies one occurrence of a measurement row.
      occurrence <- cumsum(c(TRUE, diff(value$mode) <= 0L |
        utils::head(value$measurement, -1L) != utils::tail(value$measurement, -1L)))
      key <- paste(occurrence, value$group, sep = ":")
      groups <- split(seq_len(nrow(value)), factor(key, levels = unique(key)))
      value <- do.call(rbind, lapply(groups, function(i) data.frame(
        measurement = value$measurement[i[[1L]]], group = value$group[i[[1L]]], size = length(i),
        amplitude_sum = sum(value$amplitude[i]), evidence_sum = sum(value$evidence[i]),
        prediction_cost = sum(value$prediction_cost[i]), gain = sum(value$gain[i]))))
    } else value <- data.frame(measurement = character(), group = integer(), size = integer(),
      amplitude_sum = numeric(), evidence_sum = numeric(), prediction_cost = numeric(), gain = numeric())
  }
  rownames(value) <- row.names
  value
}

#' @export
format.effect_geometry_score <- function(x, ...) {
  .check_no_extra_arguments("format.effect_geometry_score", ...)
  x <- .validate_geometry_score(x)
  sprintf("geometry_score<%d measurements; independent %s prediction; squared geometry units>",
    nrow(x$table), x$parameters$component)
}

#' @export
print.effect_geometry_score <- function(x, ...) {
  .check_no_extra_arguments("print.effect_geometry_score", ...)
  cat(format(x), "\n", sep = "")
  cat("Gain = 2 * held-out inner product - prediction norm squared.\n")
  cat(sprintf("Training/evaluation: %d/%d used partitions; rank and penalty were frozen.\n",
    length(x$training$endpoints), length(x$evaluation$endpoints)))
  print(utils::head(x$table, 6L), row.names = FALSE)
  invisible(x)
}
