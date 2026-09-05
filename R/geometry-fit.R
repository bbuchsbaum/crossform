# Fitting a frozen representational prediction (layer 5) --------------------

.geometry_fit_diagnostic_names <- c("measurement", "rank", "prediction_norm_sq",
  "penalty_trace", "compressed_residual_sq", "objective", "source_negative_mass",
  "shifted_negative_mass", "shifted_truncated_positive_mass")

.geometry_fit_from_form <- function(source, prepared, pool, parameters,
                                    storage, storage_path, row_block, retain_source, cache = NULL) {
  n <- prepared$plan$measurements
  .validate_effect_form(source, probe = FALSE)
  .require_effect_form_component(source, prepared$component)
  r <- prepared$basis$dimension
  budget <- min(parameters$rank, pool$dimension)
  layout <- .geometry_prediction_layout(r, pool$dimension, budget)
  output <- if (storage == "memory") matrix(0, n, layout$width) else NULL
  store <- if (storage == "block") .file_geometry_store(
    file.path(storage_path, "prediction.bin"), c(n, layout$width), create = TRUE,
    codec = "rectangular") else NULL
  diagnostics <- data.frame(measurement = prepared$target$measurement_ids,
    rank = integer(n), prediction_norm_sq = numeric(n), penalty_trace = numeric(n),
    compressed_residual_sq = numeric(n), objective = numeric(n),
    source_negative_mass = numeric(n), shifted_negative_mass = numeric(n),
    shifted_truncated_positive_mass = numeric(n))
  hashes <- character(n)
  for (start in seq.int(1L, n, by = row_block)) {
    rows <- seq.int(start, min(n, start + row_block - 1L))
    packed <- .geometry_component_validated(source, prepared$component, rows = rows)
    values <- matrix(0, length(rows), layout$width)
    for (i in seq_along(rows)) {
      S <- .unsvec_symmetric(packed[i, ], r)
      path <- if (is.null(cache)) .geometry_fit_path(S, pool, parameters$penalty) else
        .geometry_cache_path(cache, rows[[i]], S, pool, parameters$penalty)
      fit <- .geometry_fit_rank(path, parameters$rank)
      factor <- matrix(0, r, budget)
      amplitudes <- numeric(budget)
      if (fit$rank_effective > 0L) {
        active <- seq_len(fit$rank_effective)
        factor[, active] <- fit$factor
        amplitudes[active] <- fit$amplitudes
      }
      values[i, ] <- c(as.numeric(factor), amplitudes,
        path$signed_spectrum, path$shifted_spectrum)
      hashes[rows[[i]]] <- .sha256_signature(as.numeric(values[i, ]))
      diagnostics$rank[rows[[i]]] <- fit$rank_effective
      for (field in .geometry_fit_diagnostic_names[-c(1, 2)]) {
        diagnostics[rows[[i]], field] <- fit[[field]]
      }
    }
    if (storage == "memory") output[rows, ] <- values else {
      .write_geometry_tile(store, rows, seq_len(layout$width), values)
    }
  }
  if (storage == "memory") store <- .memory_geometry_store(output, codec = "rectangular")
  .new_geometry_fit(prepared$basis, pool, parameters, prepared$target, prepared$support,
    layout, store, hashes, diagnostics, if (retain_source) source else NULL, source$receipt,
    prepared$workspace)
}

.geometry_fit_execute <- function(prepared, pool, parameters, storage, storage_path,
                                  row_block, retain_source, cache = NULL) {
  cache_bytes <- 0
  if (!is.null(cache)) {
    .geometry_cache_check(cache, "training")
    .geometry_cache_require_recipe(cache, pool, parameters$penalty)
    if (storage != "memory") .input_error("The local candidate cache currently retains memory geometry only.")
    cache_bytes <- .geometry_cache_bind(cache, prepared, row_block)
  }
  budget <- if (is.null(cache)) min(parameters$rank, pool$dimension) else prepared$basis$dimension
  prepared <- .geometry_prediction_workspace(prepared, budget,
    storage, row_block, "fit", cache_bytes = cache_bytes)
  row_block <- prepared$workspace$row_block
  created <- FALSE
  success <- FALSE
  if (storage == "block") {
    .check_string(storage_path, "storage_path", what = "a new durable prediction directory")
    if (file.exists(storage_path)) .input_error("Refusing to overwrite an existing prediction directory.")
    if (!dir.create(storage_path, recursive = TRUE)) .input_error("Cannot create the prediction directory.")
    storage_path <- normalizePath(storage_path, mustWork = TRUE)
    created <- TRUE
  } else if (!is.null(storage_path)) {
    .input_error("`storage_path` is used only with block storage.")
  }
  on.exit(if (created && !success) unlink(storage_path, recursive = TRUE), add = TRUE)
  source_path <- if (storage == "block") file.path(storage_path, "training") else NULL
  source <- if (is.null(cache)) .geometry_fit_materialize(prepared, storage, source_path) else
    .geometry_cache_source(cache, prepared)
  result <- .geometry_fit_from_form(source, prepared, pool, parameters,
    storage, storage_path, row_block, retain_source, cache = cache)
  .validate_geometry_fit(result)
  if (storage == "block") {
    if (!retain_source) unlink(source_path, recursive = TRUE)
    pending <- file.path(storage_path, "fit.pending.rds")
    saveRDS(result, pending)
    if (!file.rename(pending, file.path(storage_path, "fit.rds"))) {
      .input_error("Cannot finalize the frozen prediction record.")
    }
  }
  success <- TRUE
  result
}

#' Fit a regularized representational form
#'
#' Learn a positive semidefinite prediction from signed cross-partition
#' geometry in a declared model family. Model eigenvalues determine the
#' inverse-kernel penalty; overlapping models share one union space.
#' Fitting is biased and adaptive. Independent predictive evidence is a
#' separate operation on a frozen fit and disjoint observation origins.
#'
#' @param plan A geometry plan over the training pairing. Only positive-weight
#'   pairing endpoints are read. The current complete-form executor admits
#'   undirected self forms with implicit identity or a fixed SPD neural metric.
#' @param basis A data-independent [model_basis()], preferably declaring
#'   `conditions = relation$effect_space` and `normalize = "trace"`.
#' @param weights Named nonnegative model weights summing to one. Required
#'   for multiple models; a single model defaults to weight one. Zero weights
#'   remove that model's exclusive directions from admissible support.
#' @param rank Required nonnegative response rank budget. Zero gives the zero
#'   predictor; an oversized budget clamps to positive pooled support. A cut
#'   through a positive tied eigenspace refuses rather than choosing an
#'   arbitrary prediction. This rank is separate from each model's input rank.
#' @param penalty Required finite nonnegative inverse-kernel penalty. Declare
#'   it before evaluation or choose it using training-only validation.
#' @param component The signed target: `"total"`, `"coherent"`, or
#'   `"configuration"`. It remains fixed when the prediction is evaluated.
#' @param storage `"memory"` or `"block"` for signed computation and retained
#'   prediction factors/spectra.
#' @param storage_path New durable directory for block storage. A completed
#'   record is saved as `fit.rds`; existing directories are never overwritten.
#' @param row_block Positive number of measurements processed together.
#' @param retain_source Retain the signed compressed training form as `$source`.
#'   Otherwise only its signed spectrum and fit diagnostics are retained.
#' @param loss The implemented loss is `"frobenius"` on effect geometry.
#'   RDM or covariance-weighted loss requires a different numerical solver.
#' @return A frozen `effect_geometry_fit` with model/support specifications,
#'   parameters, exact target identities, actual training dependencies,
#'   reduced prediction factors and spectra, diagnostics, and execution
#'   receipt. It has no raw-geometry unbiasedness or inference capability.
#'   `as.data.frame()` reads its per-measurement diagnostics without reading
#'   neural sources. Block records reopen with `readRDS("<path>/fit.rds")`.
#'
#' @details For pooled kernel `K = Q U D U' Q'` and signed reduced form
#'   `S = Q' G Q`, the estimator is `F = Q A Q'`, where
#'   `A = U [U' S U - penalty * D^-1]_{+,rank} U'`. The bracket keeps only
#'   the largest positive roots within the rank budget. The source is never
#'   clipped before compression or penalty shifting. Diagnostics distinguish
#'   negative source mass, negative shifted mass and truncated shifted mass.
#'
#'   The reported residual and objective are in the model union coordinates.
#'   Compression does not identify the whole-geometry residual or an explained
#'   fraction. A zero-penalty full-span model is labelled a generic baseline.
#'   Positive anisotropic penalties retain model-specific preferences.
#'   With a trace-normalized model kernel the penalty has geometry units;
#'   changing the neural frame normalization changes its appropriate scale.
#'   The source and any retained source spectrum stay signed. See
#'   `vignette("predictive-geometry")` for independent scoring, model-rank
#'   choices and the distinction from new-condition prediction.
#' @seealso [model_basis()], [score_geometry()], [model_geometry()]
#' @examples
#' d <- abstract_domain(3, id = "fit-example")
#' effect_names <- c("face", "body", "house", "tool")
#' B <- matrix(seq_len(12)/10, 4, dimnames = list(effect_names, NULL))
#' rel <- relation(list(run1 = B, run2 = B + 0.05), domain = d)
#' model <- model_basis(features = matrix(c(1, 1, -1, -1), 4),
#'   conditions = rel$effect_space, normalize = "trace")
#' plan <- plan_geometry(rel, compile_frame(whole_brain(), d),
#'   cross_partitions(rel, independence = "independent", generalizes_over = "run"))
#' fit <- fit_geometry(plan, model, rank = 1, penalty = 0.01)
#' fit
#' as.data.frame(fit)
#' @export
fit_geometry <- function(plan, basis, weights = NULL, rank, penalty,
                          component = c("total", "coherent", "configuration"),
                          storage = c("memory", "block"), storage_path = NULL,
                          row_block = 128L, retain_source = FALSE, loss = "frobenius") {
  if (missing(rank) || missing(penalty)) {
    .input_error("Declare both `rank` and `penalty` before fitting the geometry.")
  }
  if (!identical(loss, "frobenius")) {
    .capability_refusal("The spectral estimator implements Frobenius geometry loss only.",
      capability = "predictive_geometry_loss", namespace = "predictive_geometry",
      reasons = "predictive_loss_not_implemented", remedies = "Use loss = 'frobenius' for this estimator.")
  }
  rank <- .check_count(rank, "rank", min = 0L, max = .Machine$integer.max)
  .check_number(penalty, "penalty", nonnegative = TRUE)
  row_block <- .check_count(row_block, "row_block", max = .Machine$integer.max)
  .check_flag(retain_source, "retain_source")
  component <- match.arg(component)
  storage <- match.arg(storage)
  basis <- .validate_model_basis(basis)
  pool <- .model_basis_pool(basis, weights)
  prepared <- .geometry_fit_prepare(plan, basis, component)
  parameters <- list(rank = rank, penalty = unname(penalty), weights = pool$weights,
    component = component, loss = loss)
  .geometry_fit_execute(prepared, pool, parameters, storage, storage_path, row_block, retain_source)
}

#' @export
as.data.frame.effect_geometry_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
  .check_no_extra_arguments("as.data.frame.effect_geometry_fit", ...)
  x <- .validate_geometry_fit(x)
  value <- x$diagnostics
  if (!is.null(row.names)) rownames(value) <- row.names
  value
}

#' @export
format.effect_geometry_fit <- function(x, ...) {
  .check_no_extra_arguments("format.effect_geometry_fit", ...)
  x <- .validate_geometry_fit(x)
  sprintf("geometry_fit<%d measurements; model support %d; rank <= %d; penalty %g; %s>",
    length(x$target$measurement_ids), x$pool$dimension, x$layout$budget,
    x$parameters$penalty, x$interpretation$preference)
}

#' @export
print.effect_geometry_fit <- function(x, ...) {
  .check_no_extra_arguments("print.effect_geometry_fit", ...)
  cat(format(x), "\n", sep = "")
  cat("Learned PSD prediction; fitting is biased. Source spectra remain signed.\n")
  cat(sprintf("Training: %d used partitions; component: %s.\n",
    length(x$training$endpoints), x$parameters$component))
  if (!x$training$known) cat("Training ancestry is undeclared; independent scoring is unavailable.\n")
  if (x$interpretation$generic_full_span_baseline) cat("Generic unregularized full-span baseline.\n")
  print(utils::head(x$diagnostics[, c("measurement", "rank", "prediction_norm_sq",
    "compressed_residual_sq"), drop = FALSE], 6L), row.names = FALSE)
  invisible(x)
}
