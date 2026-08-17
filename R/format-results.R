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
    .input_error("Result values and measurement index have different lengths.")
  }
  if (is.null(colnames(values))) {
    colnames(values) <- paste0("value", seq_len(ncol(values)))
  }
  data.frame(.result_index_data_frame(index), values,
    check.names = FALSE, row.names = NULL)
}

# A value that must be shown whole. `.pf_emit()` truncates a line past
# `.pf_line_width`, which is right for a signature or a route but wrong for a
# contrast: `?contrast_energy` tells the reader to check exactly that line
# against their intended alignment, and a forty-effect contrast was arriving
# as `... wo...`. A wrapped value is emitted across as many lines as it needs,
# each continuation indented to the value column.
.pf_wrap <- function(items) structure(list(items = as.character(items)),
  class = "pf_wrap")

# Greedily pack already-formatted items into `, `-separated lines of at most
# `width` characters. Never breaks inside an item, so `house -0.5` stays whole
# however narrow the column is.
.pf_fill <- function(items, width) {
  if (!length(items)) {
    return("none")
  }
  lines <- character()
  current <- ""
  for (i in seq_along(items)) {
    piece <- if (i < length(items)) paste0(items[[i]], ",") else items[[i]]
    candidate <- if (nzchar(current)) paste(current, piece) else piece
    if (nchar(candidate) > width && nzchar(current)) {
      lines <- c(lines, current)
      current <- piece
    } else {
      current <- candidate
    }
  }
  c(lines, current)
}

# `.pf_emit()` with support for `.pf_wrap()` values. Plain values keep the
# existing behaviour exactly, so every view that does not wrap prints
# byte-identically to before.
.pf_emit_block <- function(class, fields) {
  keep <- !vapply(fields, is.null, logical(1))
  fields <- fields[keep]
  if (!any(vapply(fields, inherits, logical(1), "pf_wrap"))) {
    return(.pf_emit(class, fields))
  }
  cat("<", class, ">\n", sep = "")
  keys <- names(fields)
  width <- max(nchar(keys)) + 2L
  indent <- strrep(" ", width + 2L)
  for (i in seq_along(fields)) {
    value <- fields[[i]]
    head <- sprintf("  %-*s", width, paste0(keys[[i]], ":"))
    if (inherits(value, "pf_wrap")) {
      wrapped <- .pf_fill(value$items, .pf_line_width - width - 2L)
      cat(paste0(c(head, rep(indent, length(wrapped) - 1L)), wrapped,
        collapse = "\n"), "\n", sep = "")
    } else {
      line <- paste0(head, paste0(as.character(value), collapse = ""))
      if (nchar(line) > .pf_line_width) {
        line <- paste0(substr(line, 1L, .pf_line_width - 3L), "...")
      }
      cat(line, "\n", sep = "")
    }
  }
  invisible(NULL)
}

.format_result_preview <- function(x, title, fields = list(), ...) {
  data <- as.data.frame(x)
  # `.pf_emit_block()` aligns the header block; `measurements` is the longest
  # key, so adding a field never shifts the existing columns.
  .pf_emit_block(title, c(list(measurements = nrow(data)), fields))
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
  cat("  residual:     ", .pf_residual_noise(x$source), "\n", sep = "")
  cat("  storage:      exact factorized covariance\n", sep = "")
  cat("  spatial law:  local marginal only\n", sep = "")
  if (!is.null(x$source$execution)) {
    cat("  execution:    ", x$source$execution$route, " / ",
      x$source$execution$residual_strategy, "\n", sep = "")
  }
  invisible(x)
}

#' @export
print.effect_rdm_sampling_covariance_batch <- function(x, ...) {
  if (!length(x) || !inherits(x[[1L]], "effect_rdm_sampling_covariance")) {
    .input_error("`x` must come from `rdm_sampling_covariance()`.")
  }
  first <- x[[1L]]
  execution <- first$source$execution
  cat("<effect_rdm_sampling_covariance_batch>\n", sep = "")
  cat("  measurements: ", length(x), "\n", sep = "")
  cat("  distances:    ", first$dimension, "\n", sep = "")
  if (!is.null(execution)) {
    cat("  execution:    ", execution$route, " / ",
      execution$residual_strategy, "\n", sep = "")
    cat("  shared residual statistics: ",
      if (isTRUE(execution$shared_residual_statistics)) "yes" else "no",
      "\n", sep = "")
  }
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
  if (!.is_finite_numeric(x$values) || length(x$values) != NROW(x$index)) {
    .input_error("Crossnobis result values or measurement index are invalid.")
  }
  .bind_result_values(x$index,
    matrix(x$values, ncol = 1L, dimnames = list(NULL, "crossnobis")))
}

# The contrast in relation order, one `name value` item per effect. Named so
# that a positionally supplied contrast shows the alignment that was actually
# used. Capped at `max` effects, because past that a print method is no longer
# a summary; `x$weights` is the whole vector.
.pf_contrast_items <- function(weights, max = 12L) {
  if (is.null(weights) || !length(weights)) {
    return("none")
  }
  labels <- names(weights)
  if (is.null(labels) || anyNA(labels) || !all(nzchar(labels))) {
    labels <- paste0("[", seq_along(weights), "]")
  }
  items <- paste0(labels, " ",
    format(signif(as.numeric(weights), 4L), trim = TRUE))
  if (length(items) <= max) {
    return(items)
  }
  c(items[seq_len(max)], sprintf("(+%d more)", length(items) - max))
}

#' @export
print.effect_crossnobis_view <- function(x, ...) {
  .format_result_preview(x, "effect_crossnobis_view",
    fields = list(contrast = .pf_wrap(.pf_contrast_items(x$contrast))), ...)
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
    fields = list(contrast = .pf_wrap(.pf_contrast_items(x$weights))), ...)
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
  # Without a header the columns of the coefficient table are the only clue
  # to what was fitted, and `(Intercept)` next to a model name does not say
  # which names are models, which are nuisance, or over how many pairs the
  # regression ran.
  terms <- x$terms
  named <- function(kind) {
    if (is.data.frame(terms)) terms$term[terms$role == kind] else character()
  }
  pairs <- x$query$pair_labels
  effects <- x$query$effects
  .format_result_preview(x, "effect_rsa_view",
    fields = list(
      models = .pf_wrap(named("model")),
      nuisance = if (length(named("nuisance"))) {
        .pf_wrap(named("nuisance"))
      },
      intercept = if (length(named("intercept"))) "fitted" else "omitted",
      pairs = if (!is.null(pairs)) {
        sprintf("%s over %s", .msg_count(length(pairs), "distance"),
          .msg_count(length(effects), "effect"))
      },
      estimate = sprintf("OLS coefficients on the %s component", x$component)
    ), ...)
  cat(sprintf("  %-14s%s\n", "next:",
    "as.data.frame(x), plot(x, terms = ...)"))
  invisible(x)
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
