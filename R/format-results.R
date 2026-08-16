# Compact printing and row-wise coercion ------------------------------------

.result_index_data_frame <- function(index) {
  if (is.data.frame(index)) {
    value <- index
    rownames(value) <- NULL
    return(value)
  }
  data.frame(measurement = index, check.names = FALSE,
    stringsAsFactors = FALSE)
}

.bind_result_values <- function(index, values) {
  if (is.null(dim(values))) values <- matrix(values, ncol = 1L)
  values <- as.matrix(values)
  if (nrow(values) != NROW(index)) {
    stop("Result values and measurement index have different lengths.",
      call. = FALSE)
  }
  if (is.null(colnames(values))) {
    colnames(values) <- paste0("value", seq_len(ncol(values)))
  }
  data.frame(.result_index_data_frame(index), values,
    check.names = FALSE, row.names = NULL)
}

.format_result_preview <- function(x, title, fields = list(), ...) {
  data <- as.data.frame(x)
  # `.pf_emit()` aligns the header block; `measurements` is the longest key,
  # so adding a field never shifts the existing columns.
  .pf_emit(title, c(list(measurements = nrow(data)), fields))
  # Four significant digits. The default seven prints the last bits of a
  # LAPACK result, which differ across BLAS implementations and platforms
  # and would make printed output non-reproducible.
  arguments <- list(utils::head(data), row.names = FALSE, ...)
  if (!"digits" %in% names(arguments)) {
    arguments$digits <- 4L
  }
  do.call(print, arguments)
  if (nrow(data) > 6L) {
    cat("  ... ", nrow(data) - 6L, " more measurements\n", sep = "")
  }
  invisible(x)
}

.format_counted_result <- function(class, measurements, detail = NULL) {
  suffix <- if (is.null(detail) || !nzchar(detail)) "" else paste0(", ", detail)
  sprintf("<%s: %d measurements%s>", class, as.integer(measurements), suffix)
}

#' @export
format.effect_geometry_plan <- function(x, detail = FALSE, ...) {
  .validate_geometry_plan(x, deep = FALSE)
  if (isTRUE(detail)) {
    return(.geometry_plan_lines(x, detail = TRUE))
  }
  .format_counted_result(
    "effect_geometry_plan", x$measurements,
    paste0(x$logical_shape[[1L]], " effects, ",
      nrow(x$pairing), " partition pairs")
  )
}

#' @export
format.effect_crossnobis_plan <- function(x, ...) {
  .validate_crossnobis_plan(x, deep = FALSE)
  .format_counted_result(
    "effect_crossnobis_plan", nrow(x$frame$weights),
    paste0(x$metric_schedule$recipe$kind, ", ",
      x$metric_schedule$training_policy$kind)
  )
}

#' @export
format.effect_geometry <- function(x, ...) {
  .validate_effect_geometry(x, probe = FALSE)
  .format_counted_result(
    "effect_geometry", x$total$dim[[1L]],
    paste0(x$logical_shape[[1L]], " effects, ",
      paste(x$storage, collapse = "+"), " storage")
  )
}

#' @export
format.effect_relation_fit <- function(x, ...) {
  .validate_relation_fit(x, deep = FALSE)
  .format_counted_result(
    "effect_relation_fit", length(x$relation$partitions),
    paste0(length(x$relation$effects), " effects, ",
      x$relation$n_features, " features")
  )
}

#' @export
format.effect_rdm_sampling_covariance <- function(x, ...) {
  .validate_sampling_covariance(x, deep = FALSE)
  .format_counted_result(
    "effect_rdm_sampling_covariance", x$dimension,
    paste0(x$plan$target$policy, ", factorized")
  )
}

#' @export
format.effect_view <- function(x, ...) {
  .validate_effect_view(x)
  .format_counted_result("effect_view", nrow(x$values),
    paste0(ncol(x$values), " queries"))
}

#' @export
format.effect_crossnobis_view <- function(x, ...) {
  .format_counted_result("effect_crossnobis_view", length(x$values),
    "signed estimates")
}

#' @export
format.effect_contrast_view <- function(x, ...) {
  .format_counted_result("effect_contrast_view", length(x$total),
    "signed + energy decomposition")
}

#' @export
format.effect_rdm_view <- function(x, ...) {
  .format_counted_result("effect_rdm_view", nrow(x$values),
    paste0(ncol(x$values), " distances"))
}

#' @export
format.effect_rsa_view <- function(x, ...) {
  .format_counted_result("effect_rsa_view", nrow(x$coefficients),
    paste0(ncol(x$coefficients), " coefficients"))
}

#' @export
format.effect_spectrum_view <- function(x, ...) {
  .format_counted_result("effect_spectrum_view", nrow(x$values),
    paste0(ncol(x$values), " roots"))
}

#' @export
print.effect_geometry <- function(x, ...) {
  .validate_effect_geometry(x, probe = FALSE)
  cat("<effect_geometry>\n", sep = "")
  cat("  effects:      ", x$logical_shape[[1L]], " x ",
    x$logical_shape[[2L]], "\n", sep = "")
  cat("  measurements: ", x$total$dim[[1L]], "\n", sep = "")
  cat("  components:   total + coherent + configuration\n", sep = "")
  cat("  storage:      ", paste(x$storage, collapse = ", "), "\n", sep = "")
  cat("  estimate:     signed cross-generalized form; PSD not assumed\n",
    sep = "")
  invisible(x)
}

#' @export
print.effect_relation_fit <- function(x, ...) {
  .validate_relation_fit(x, deep = FALSE)
  capability <- relation_fit_capabilities(x)
  cat("<effect_relation_fit>\n", sep = "")
  cat("  effects:      ", length(x$relation$effects), "\n", sep = "")
  cat("  features:     ", x$relation$n_features, "\n", sep = "")
  cat("  partitions:   ", length(x$relation$partitions), "\n", sep = "")
  cat("  residuals:    ", sum(capability$residual_blocks), "/",
    nrow(capability), " partitions\n", sep = "")
  cat("  covariance:   ", sum(capability$effect_covariance), "/",
    nrow(capability), " partitions\n", sep = "")
  invisible(x)
}

#' @export
print.effect_rdm_sampling_covariance <- function(x, ...) {
  .validate_sampling_covariance(x, deep = FALSE)
  cat("<effect_rdm_sampling_covariance>\n", sep = "")
  cat("  distances:    ", x$dimension, "\n", sep = "")
  cat("  measurement:  ", x$source$node, "\n", sep = "")
  cat("  partitions:   ", x$partitions, " (dependent pair products)\n",
    sep = "")
  cat("  target:       ", x$plan$target$target, " / ",
    x$plan$target$policy, "\n", sep = "")
  cat("  metric:       fixed\n", sep = "")
  cat("  storage:      exact factorized covariance\n", sep = "")
  cat("  spatial law:  local marginal only\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.effect_view <- function(x, row.names = NULL, optional = FALSE,
                                      ...) {
  .validate_effect_view(x)
  .bind_result_values(x$index, x$values)
}

#' @export
print.effect_view <- function(x, ...) {
  .validate_effect_view(x)
  .format_result_preview(x, "effect_view", ...)
}

#' @export
as.data.frame.effect_crossnobis_view <- function(x, row.names = NULL,
                                                 optional = FALSE, ...) {
  if (!is.numeric(x$values) || any(!is.finite(x$values)) ||
      length(x$values) != NROW(x$index)) {
    stop("Crossnobis result values or measurement index are invalid.",
      call. = FALSE)
  }
  .bind_result_values(x$index,
    matrix(x$values, ncol = 1L, dimnames = list(NULL, "crossnobis")))
}

#' @export
print.effect_crossnobis_view <- function(x, ...) {
  .format_result_preview(x, "effect_crossnobis_view",
    fields = list(contrast = .pf_weights(x$contrast)), ...)
}

#' @export
as.data.frame.effect_contrast_view <- function(x, row.names = NULL,
                                               optional = FALSE, ...) {
  signed <- if (is.list(x$signed)) {
    do.call(cbind, x$signed)
  } else {
    matrix(x$signed, ncol = 1L, dimnames = list(NULL, "signed"))
  }
  if (is.list(x$signed)) {
    colnames(signed) <- paste0("signed_", names(x$signed))
  }
  values <- cbind(
    signed,
    coherent = x$coherent,
    configuration = x$configuration,
    total = x$total,
    coherence_fraction = x$coherence_fraction
  )
  .bind_result_values(x$index, values)
}

#' @export
print.effect_contrast_view <- function(x, ...) {
  # The contrast is shown because unnamed weights are accepted positionally:
  # a reader who mis-ordered them sees the alignment that was actually used.
  .format_result_preview(x, "effect_contrast_view",
    fields = list(contrast = .pf_weights(x$weights)), ...)
  # A column of NA in the preview is otherwise unexplained, and the reason is
  # a real property of the estimate rather than missing data.
  valid <- x$coherence_fraction_valid
  if (is.null(valid)) valid <- !is.na(x$coherence_fraction)
  cat(strwrap(sprintf(paste0(
    "coherence_fraction: %d of %d valid; NA where coherent and ",
    "configuration are not a nonnegative partition"
  ), sum(valid), length(valid)),
    width = .pf_line_width, prefix = "    ", initial = "  "), sep = "\n")
  invisible(x)
}

#' @export
as.data.frame.effect_rdm_view <- function(x, row.names = NULL,
                                          optional = FALSE, ...) {
  .bind_result_values(x$index, x$values)
}

#' @export
print.effect_rdm_view <- function(x, ...) {
  .format_result_preview(x, "effect_rdm_view", ...)
}

#' @export
as.data.frame.effect_rsa_view <- function(x, row.names = NULL,
                                          optional = FALSE, ...) {
  .bind_result_values(x$index, x$coefficients)
}

#' @export
print.effect_rsa_view <- function(x, ...) {
  .format_result_preview(x, "effect_rsa_view", ...)
}

#' @export
as.data.frame.effect_spectrum_view <- function(x, row.names = NULL,
                                               optional = FALSE, ...) {
  .bind_result_values(x$index, x$values)
}

#' @export
print.effect_spectrum_view <- function(x, ...) {
  .format_result_preview(x, "effect_spectrum_view", ...)
}
