# Equal-axis population component views --------------------------------------

#' Build an equal-axis population component view
#'
#' Converts a [population_decomposition()] object into directly inspectable
#' plotted data for total, coherent, and configuration coefficients. One
#' symmetric axis limit is computed across every selected panel, so equal
#' signed magnitudes always receive equal visual magnitude.
#'
#' @param x An `effect_population_decomposition`.
#' @param term One population model coefficient name.
#' @param query Optional query names. The default keeps every query.
#' @return An `effect_population_component_view` with `data`, `coverage`,
#'   one `axis` contract, and a plotting receipt.
#' @export
population_component_view <- function(x, term, query = NULL) {
  if (!inherits(x, "effect_population_decomposition")) {
    .input_error("`x` must be returned by `population_decomposition()`.",
      arg = "x", received = .msg_value(x),
      expected = "an `effect_population_decomposition`")
  }
  terms <- dimnames(x$coefficient_gap)[[3L]]
  if (!.is_string(term) || !term %in% terms) {
    .input_error(sprintf("`term` must name one coefficient: %s.",
      .msg_names(terms)), arg = "term", received = .msg_value(term),
      expected = paste0("one of ", .msg_names(terms)))
  }
  queries <- dimnames(x$coefficient_gap)[[2L]]
  if (is.null(query)) query <- queries
  if (!.is_strings(query, unique = TRUE) || !all(query %in% queries)) {
    .input_error(sprintf("`query` must select from %s.", .msg_names(queries)),
      arg = "query", received = .msg_value(query),
      expected = paste0("names from ", .msg_names(queries)))
  }
  components <- c("total", "coherent", "configuration")
  labels <- c(total = "Total", coherent = "Coherent component",
              configuration = "Configurational component")
  rows <- lapply(components, function(component) {
    uncertainty <- x$uncertainty[[component]]$between
    estimate <- uncertainty$estimate[, query, term, drop = FALSE]
    se <- uncertainty$se[, query, term, drop = FALSE]
    lower <- uncertainty$lower[, query, term, drop = FALSE]
    upper <- uncertainty$upper[, query, term, drop = FALSE]
    grid <- expand.grid(
      node = dimnames(estimate)[[1L]], query = dimnames(estimate)[[2L]],
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    flatten <- function(value) as.numeric(aperm(value, c(1L, 2L, 3L)))
    data.frame(
      component = component, component_label = unname(labels[[component]]),
      node = grid$node, query = grid$query, term = term,
      estimate = flatten(estimate), se = flatten(se),
      lower = flatten(lower), upper = flatten(upper),
      sign = ifelse(flatten(estimate) < 0, "negative",
                    ifelse(flatten(estimate) > 0, "positive", "zero")),
      units = "signed transported evidence coefficient",
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, rows)
  rownames(data) <- NULL
  finite <- unlist(data[c("estimate", "lower", "upper")], use.names = FALSE)
  finite <- finite[is.finite(finite)]
  limit <- if (length(finite)) max(abs(finite)) else 1
  if (!is.finite(limit) || limit <= 0) limit <- 1
  data$visual_magnitude <- abs(data$estimate) / limit
  node_position <- match(data$node, as.character(x$index$node))
  data$scale <- if ("scale" %in% names(x$index)) {
    as.numeric(x$index$scale[node_position])
  } else NA_real_
  coverage_rows <- expand.grid(
    node = dimnames(x$coverage$n)[[1L]], query = query,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  at <- cbind(match(coverage_rows$node, dimnames(x$coverage$n)[[1L]]),
              match(coverage_rows$query, dimnames(x$coverage$n)[[2L]]))
  coverage <- data.frame(
    coverage_rows,
    coverage_policy = x$coverage$policy,
    subject_set_id = x$coverage$subject_set_id[at],
    n = x$coverage$n[at], fraction = x$coverage$fraction[at],
    n_eff = x$coverage$n_eff[at], mass_n_eff =
      x$coverage$mass_n_eff[match(coverage_rows$node,
                                 names(x$coverage$mass_n_eff))],
    stringsAsFactors = FALSE
  )
  structure(list(
    contract_version = "population-component-view-v1",
    parent_population_plan_id = x$parent_population_plan_id,
    term = term, queries = query, components = components,
    data = data, coverage = coverage,
    axis = list(type = "shared_symmetric", limits = c(-limit, limit),
      zero = 0, sign = "positive_up_negative_down",
      units = "signed transported evidence coefficient"),
    receipt = list(
      aggregation = "population_coefficient_not_participant_ratio_average",
      uncertainty = x$estimator,
      comparability_checked = TRUE,
      coverage_preserved = TRUE
    )
  ), class = "effect_population_component_view")
}

#' @export
print.effect_population_component_view <- function(x, ...) {
  cat("<effect_population_component_view>\n")
  cat("  term: ", x$term, "; queries: ", paste(x$queries, collapse = ", "),
      "\n", sep = "")
  cat("  shared axis: [", paste(format(x$axis$limits), collapse = ", "),
      "]\n", sep = "")
  cat("  rows: ", nrow(x$data), "; coverage rows: ", nrow(x$coverage),
      "\n", sep = "")
  invisible(x)
}

#' @export
plot.effect_population_component_view <- function(x, ...) {
  components <- x$components
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(length(x$queries), length(components)),
                mar = c(5, 4, 3, 1))
  for (query in x$queries) {
    for (component in components) {
      data <- x$data[x$data$query == query & x$data$component == component, ]
      position <- seq_len(nrow(data))
      graphics::plot(position, data$estimate, ylim = x$axis$limits,
        xaxt = "n", xlab = "group node / scale", ylab = x$axis$units,
        main = paste(data$component_label[[1L]], query, sep = "\n"),
        pch = 19, ...)
      graphics::axis(1, at = position, labels = data$node, las = 2)
      graphics::abline(h = 0, lty = 2, col = "grey50")
      graphics::segments(position, data$lower, position, data$upper)
    }
  }
  invisible(x)
}

#' Build a population scale profile with explicit pointwise uncertainty
#'
#' @param x An `effect_population_decomposition`.
#' @param term One population coefficient name.
#' @param query Optional query names.
#' @param interval Explicit interval method: `"HC3"`, `"classical"`, or
#'   `"wild_bootstrap"`.
#' @param bootstrap For `interval = "wild_bootstrap"`, a named list containing
#'   total, coherent, and configuration bootstrap objects made from the three
#'   component results and the same coefficient contrast.
#' @return An `effect_population_scale_profile`. Bands are pointwise; no
#'   simultaneous or maxT interpretation is supplied.
#' @export
population_scale_profile <- function(
    x, term, query = NULL,
    interval = c("HC3", "classical", "wild_bootstrap"),
    bootstrap = NULL) {
  if (!inherits(x, "effect_population_decomposition")) {
    .input_error("`x` must be returned by `population_decomposition()`.",
      arg = "x", received = .msg_value(x),
      expected = "an `effect_population_decomposition`")
  }
  interval <- match.arg(interval)
  temporary <- x
  if (!identical(interval, "wild_bootstrap")) {
    temporary$uncertainty <- lapply(x$results, population_uncertainty,
                                    estimator = interval)
    temporary$estimator <- interval
  }
  view <- population_component_view(temporary, term, query)
  if (identical(interval, "wild_bootstrap")) {
    components <- c("total", "coherent", "configuration")
    if (!is.list(bootstrap) || !identical(names(bootstrap), components) ||
        !all(vapply(bootstrap, inherits, logical(1),
                    "effect_population_wild_bootstrap"))) {
      .input_error(paste0(
        "`bootstrap` must be named total, coherent, and configuration and ",
        "contain three population wild-bootstrap objects."
      ), arg = "bootstrap", received = .msg_value(bootstrap),
        expected = "three named `effect_population_wild_bootstrap` objects")
    }
    terms <- dimnames(x$coefficient_gap)[[3L]]
    expected <- stats::setNames(numeric(length(terms)), terms)
    expected[[term]] <- 1
    for (component in components) {
      boot <- bootstrap[[component]]
      if (!identical(boot$parent_result_id,
                     x$results[[component]]$scientific_plan_id) ||
          !identical(boot$contrast, expected)) {
        .population_decomposition_refusal(
          "bootstrap_parent_or_coefficient_contrast_differs"
        )
      }
      rows <- view$data$component == component
      at <- cbind(match(view$data$node[rows], rownames(boot$observed_estimate)),
                  match(view$data$query[rows], colnames(boot$observed_estimate)))
      estimate <- boot$observed_estimate[at]
      se <- boot$observed_se[at]
      critical <- boot$critical_value[at]
      view$data$estimate[rows] <- estimate
      view$data$se[rows] <- se
      view$data$lower[rows] <- estimate - critical * se
      view$data$upper[rows] <- estimate + critical * se
    }
    view$receipt$uncertainty <- "null_imposed_wild_bootstrap_HC3_t"
  }
  finite <- unlist(view$data[c("estimate", "lower", "upper")], use.names = FALSE)
  finite <- finite[is.finite(finite)]
  limit <- if (length(finite)) max(abs(finite)) else 1
  if (!is.finite(limit) || limit <= 0) limit <- 1
  view$axis$limits <- c(-limit, limit)
  view$data$visual_magnitude <- abs(view$data$estimate) / limit
  at <- match(paste(view$data$node, view$data$query),
              paste(view$coverage$node, view$coverage$query))
  view$data$n <- view$coverage$n[at]
  view$data$fraction <- view$coverage$fraction[at]
  view$data$n_eff <- view$coverage$n_eff[at]
  view$data$subject_set_id <- view$coverage$subject_set_id[at]
  view$data$gap <- !is.finite(view$data$estimate) |
    !is.finite(view$data$lower) | !is.finite(view$data$upper)
  view$interval <- list(
    method = interval, semantics = "pointwise",
    simultaneous_coverage = "not_available_unimplemented_uncalibrated",
    calibration_scope = paste0(
      "Matched simulation regimes in ",
      "inst/extdata/certification/population-calibration-results.csv; ",
      "no marginal claim under informative coverage or transport estimation."
    )
  )
  view$receipt$interpolation <- "none_gaps_and_subject_set_changes_retained"
  view$contract_version <- "population-scale-profile-v1"
  class(view) <- c("effect_population_scale_profile",
                   "effect_population_component_view")
  view
}

#' @export
print.effect_population_scale_profile <- function(x, ...) {
  cat("<effect_population_scale_profile>\n")
  cat("  interval: ", x$interval$method, " (", x$interval$semantics, ")\n",
      sep = "")
  cat("  simultaneous coverage: ", x$interval$simultaneous_coverage, "\n",
      sep = "")
  cat("  rows: ", nrow(x$data), "; gaps: ", sum(x$data$gap), "\n", sep = "")
  invisible(x)
}
