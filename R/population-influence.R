# Population leave-one-subject influence -------------------------------------

#' Diagnose leave-one-subject population influence
#'
#' Computes descriptive leave-one-subject coefficient changes from the exact
#' retained transported values in a [population_decomposition()]. It also binds
#' optional cross-fitted heterogeneity output and transport/coverage provenance.
#'
#' @param x An `effect_population_decomposition`.
#' @param heterogeneity Optional `effect_population_heterogeneity` from the same
#'   planned subjects, normally with `estimator = "cross_fit"`.
#' @param mode `"bounded"` refuses work beyond `max_subjects` or `max_cells`;
#'   `"deep"` is the explicit opt-in for a larger run.
#' @param max_subjects,max_cells Positive bounds for the default mode.
#' @return An `effect_population_influence` with a long influence table,
#'   subject provenance, optional heterogeneity summary, and descriptive scope.
#' @export
population_influence <- function(
    x, heterogeneity = NULL, mode = c("bounded", "deep"),
    max_subjects = 50L, max_cells = 250000L) {
  if (!inherits(x, "effect_population_decomposition")) {
    .input_error("`x` must be returned by `population_decomposition()`.",
      arg = "x", received = .msg_value(x),
      expected = "an `effect_population_decomposition`")
  }
  mode <- match.arg(mode)
  max_subjects <- .check_count(max_subjects, "max_subjects", min = 1L,
    what = "a positive subject bound")
  max_cells <- .check_count(max_cells, "max_cells", min = 1L,
    what = "a positive output-row bound")
  subjects <- x$coverage$planned_subjects
  dimensions <- dim(x$coefficient_gap)
  requested_cells <- length(subjects) * prod(dimensions) * 3L
  if (identical(mode, "bounded") &&
      (length(subjects) > max_subjects || requested_cells > max_cells)) {
    .capability_refusal(paste0(
      "The requested leave-one-subject court exceeds the bounded default: ",
      length(subjects), " subjects and ", requested_cells, " rows."
    ), capability = "population_influence_deep_mode",
      namespace = "population_influence", reasons = "bounded_default_exceeded",
      remedies = "Review the size, then rerun with `mode = \"deep\"` explicitly.")
  }
  components <- c("total", "coherent", "configuration")
  terms <- dimnames(x$coefficient_gap)[[3L]]
  nodes <- dimnames(x$coefficient_gap)[[1L]]
  queries <- dimnames(x$coefficient_gap)[[2L]]
  rows <- list()
  for (component in components) {
    result <- x$results[[component]]
    design <- result$uncertainty$between$design
    for (node in seq_along(nodes)) for (query in seq_along(queries)) {
      available <- as.logical(result$coverage$availability[node, query, ])
      full <- result$coefficients[node, query, ]
      for (subject in seq_along(subjects)) {
        keep <- available
        keep[[subject]] <- FALSE
        X <- design[keep, , drop = FALSE]
        y <- result$values[node, query, keep]
        estimable <- nrow(X) >= ncol(X) && qr(X)$rank == ncol(X)
        loo <- if (estimable) {
          stats::setNames(drop(qr.solve(X, y)), colnames(X))
        } else stats::setNames(rep(NA_real_, length(terms)), terms)
        rows[[length(rows) + 1L]] <- data.frame(
          subject = subjects[[subject]], node = nodes[[node]],
          query = queries[[query]], term = terms, component = component,
          contributed_to_primary = available[[subject]],
          primary_subject_set_id =
            result$coverage$subject_set_id[node, query],
          primary_estimate = as.numeric(full),
          leave_one_out_estimate = as.numeric(loo[terms]),
          delta = as.numeric(loo[terms] - full),
          abs_delta = abs(as.numeric(loo[terms] - full)),
          leave_one_out_estimable = estimable,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  influence <- do.call(rbind, rows)
  rownames(influence) <- NULL
  index <- x$results$total$receipt$subjects
  index <- index[match(subjects, index$subject), , drop = FALSE]
  split_sample <- function(value) {
    if (is.na(value) || !nzchar(value)) character() else
      trimws(strsplit(value, ";", fixed = TRUE)[[1L]])
  }
  self_fit <- vapply(seq_along(subjects), function(position) {
    subjects[[position]] %in% split_sample(index$fitting_sample[[position]])
  }, logical(1))
  if (any(self_fit)) {
    .contract_error("A subject appears in its own declared transport fitting sample.")
  }
  provenance <- data.frame(
    subject = subjects, plan_id = index$plan_id,
    transport_signature = index$transport_signature,
    transport_source = index$transport_source,
    transport_status = index$transport_status,
    fitting_sample = index$fitting_sample, cross_fit = index$cross_fit,
    transport_self_fit = self_fit,
    sink_territory = index$sink_territory,
    retained_territory = 1 - index$sink_territory,
    stringsAsFactors = FALSE
  )
  at <- match(influence$subject, provenance$subject)
  influence$transport_signature <- provenance$transport_signature[at]
  influence$cross_fit <- provenance$cross_fit[at]
  influence$sink_territory <- provenance$sink_territory[at]
  influence$scale <- if ("scale" %in% names(x$index)) {
    x$index$scale[match(influence$node, x$index$node)]
  } else NA_real_
  heterogeneity_summary <- NULL
  if (!is.null(heterogeneity)) {
    if (!inherits(heterogeneity, "effect_population_heterogeneity")) {
      .input_error("`heterogeneity` must be returned by `heterogeneity()`.",
        arg = "heterogeneity", received = .msg_value(heterogeneity),
        expected = "an `effect_population_heterogeneity`")
    }
    loading_subjects <- rownames(heterogeneity$loadings)
    if (!setequal(loading_subjects, subjects)) {
      .population_decomposition_refusal("heterogeneity_subjects_differ")
    }
    heterogeneity_summary <- list(
      estimator = heterogeneity$estimator,
      spectrum = heterogeneity$spectrum,
      leading_loading = stats::setNames(
        heterogeneity$loadings[subjects, 1L], subjects
      ),
      cross_fit = heterogeneity$receipt$cross_fit,
      transport_partition_overlap =
        heterogeneity$receipt$cross_fit$transport_partition_overlap
    )
  }
  structure(list(
    contract_version = "population-influence-v1",
    parent_population_plan_id = x$parent_population_plan_id,
    mode = mode, bounds = list(max_subjects = max_subjects,
      max_cells = max_cells, requested_cells = requested_cells),
    influence = influence, subject_provenance = provenance,
    heterogeneity = heterogeneity_summary,
    interpretation = paste0(
      "Leave-one-subject changes and heterogeneity loadings are descriptive ",
      "influence diagnostics, not formal exclusion or inference rules."
    )
  ), class = "effect_population_influence")
}

#' @export
print.effect_population_influence <- function(x, ...) {
  cat("<effect_population_influence>\n")
  cat("  subjects: ", nrow(x$subject_provenance), "; rows: ",
      nrow(x$influence), "; mode: ", x$mode, "\n", sep = "")
  if (!is.null(x$heterogeneity)) {
    cat("  heterogeneity: ", x$heterogeneity$estimator, "\n", sep = "")
  }
  cat("  reading: ", x$interpretation, "\n", sep = "")
  invisible(x)
}
