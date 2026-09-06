# Population coverage and transport sensitivity diagnostics -----------------

.population_diagnostic_correlation <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || stats::sd(x[keep]) == 0 || stats::sd(y[keep]) == 0) {
    return(NA_real_)
  }
  stats::cor(x[keep], y[keep])
}

.population_diagnostic_smd <- function(value, available) {
  keep <- is.finite(value) & !is.na(available)
  value <- value[keep]
  available <- available[keep]
  if (!any(available) || all(available) || length(value) < 3L ||
      stats::sd(value) == 0) return(NA_real_)
  (mean(value[available]) - mean(value[!available])) / stats::sd(value)
}

#' Diagnose population coverage, sink exposure, and transport sensitivity
#'
#' Produces descriptive diagnostics for a fitted query-bank population result.
#' It never refits or mutates the primary estimand. Thresholded rows are
#' explicitly labelled sensitivity summaries, and associations are descriptive
#' warnings rather than causal corrections.
#'
#' @param x An `effect_population_result` returned by [estimate_population()].
#' @param minimum_coverage Minimum planned-subject fraction for a warning.
#' @param minimum_transport_quality Minimum retained-territory fraction for a
#'   subject to enter the sensitivity summary.
#' @param material_change Fraction of contributors removed at which a
#'   composition-change warning is raised.
#' @return An `effect_population_diagnostics` object with subject provenance,
#'   cell summaries, descriptive associations, sensitivity summaries, and
#'   warnings.
#' @export
population_diagnostics <- function(x, minimum_coverage = 0.8,
                                   minimum_transport_quality = 0.7,
                                   material_change = 0.2) {
  .validate_population_result(x)
  if (!identical(x$basis, "query_bank") || is.null(x$values)) {
    .capability_refusal(
      "Population diagnostics need retained subject-level query values.",
      capability = "population_diagnostic_subject_values",
      namespace = "population_diagnostics",
      reasons = "complete_form_retains_no_subject_values",
      remedies = "Use `estimate_population(plan, queries)` for diagnostics."
    )
  }
  for (argument in c("minimum_coverage", "minimum_transport_quality",
                     "material_change")) {
    value <- get(argument)
    if (!.is_number(value) || value < 0 || value > 1) {
      .input_error(paste0("`", argument, "` must be one finite number in [0, 1]."),
        arg = argument, received = .msg_value(value), expected = "a number in [0, 1]")
    }
  }

  subjects <- x$coverage$planned_subjects
  nodes <- dimnames(x$values)[[1L]]
  queries <- dimnames(x$values)[[2L]]
  index <- x$receipt$subjects
  index <- index[match(subjects, index$subject), , drop = FALSE]
  quality <- 1 - index$sink_territory
  design <- x$uncertainty$between$design
  rownames(design) <- subjects
  covariate_columns <- setdiff(colnames(design), "(Intercept)")
  subject_provenance <- data.frame(
    subject = subjects,
    plan_id = index$plan_id,
    transport_signature = index$transport_signature,
    transport_source = index$transport_source,
    transport_status = index$transport_status,
    cross_fit = index$cross_fit,
    sink_territory = index$sink_territory,
    all_sink_rows = index$all_sink_rows,
    retained_territory = quality,
    stringsAsFactors = FALSE
  )

  cells <- associations <- sensitivity <- list()
  for (node in seq_along(nodes)) {
    for (query in seq_along(queries)) {
      available <- as.logical(x$coverage$availability[node, query, ])
      response <- as.numeric(x$values[node, query, ])
      primary_subjects <- subjects[available]
      keep <- available & quality >= minimum_transport_quality
      sensitivity_subjects <- subjects[keep]
      removed_fraction <- if (length(primary_subjects)) {
        1 - length(sensitivity_subjects) / length(primary_subjects)
      } else 0
      cells[[length(cells) + 1L]] <- data.frame(
        node = nodes[[node]], query = queries[[query]],
        coverage_policy = x$coverage$policy,
        subject_set_id = x$coverage$subject_set_id[node, query],
        planned_n = length(subjects), n = sum(available),
        coverage_fraction = mean(available),
        mass_n_eff = x$coverage$mass_n_eff[[node]],
        mean_sink_territory = if (any(available))
          mean(index$sink_territory[available]) else NA_real_,
        mean_transport_quality = if (any(available))
          mean(quality[available]) else NA_real_,
        stringsAsFactors = FALSE
      )
      sensitivity[[length(sensitivity) + 1L]] <- data.frame(
        node = nodes[[node]], query = queries[[query]],
        primary_subject_set_id = x$coverage$subject_set_id[node, query],
        primary_n = sum(available), sensitivity_n = sum(keep),
        removed_fraction = removed_fraction,
        primary_mean = if (any(available)) mean(response[available]) else NA_real_,
        sensitivity_mean = if (any(keep)) mean(response[keep]) else NA_real_,
        mean_difference = if (any(available) && any(keep))
          mean(response[keep]) - mean(response[available]) else NA_real_,
        minimum_coverage = minimum_coverage,
        minimum_transport_quality = minimum_transport_quality,
        target_status = "sensitivity_descriptive_not_primary",
        primary_subjects = paste(primary_subjects, collapse = ","),
        sensitivity_subjects = paste(sensitivity_subjects, collapse = ","),
        stringsAsFactors = FALSE
      )
      diagnostic_values <- cbind(
        sink_territory = index$sink_territory,
        retained_territory = quality,
        design[, covariate_columns, drop = FALSE]
      )
      for (name in colnames(diagnostic_values)) {
        associations[[length(associations) + 1L]] <- data.frame(
          node = nodes[[node]], query = queries[[query]], variable = name,
          availability_smd = .population_diagnostic_smd(
            diagnostic_values[, name], available
          ),
          outcome_correlation = .population_diagnostic_correlation(
            diagnostic_values[, name], response
          ),
          interpretation = "descriptive_association_not_causal_correction",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  cells <- do.call(rbind, cells)
  associations <- do.call(rbind, associations)
  sensitivity <- do.call(rbind, sensitivity)
  low_coverage <- cells$coverage_fraction < minimum_coverage
  low_quality <- cells$mean_transport_quality < minimum_transport_quality
  changed <- sensitivity$removed_fraction >= material_change
  warnings <- data.frame(
    node = c(cells$node[low_coverage], cells$node[low_quality],
             sensitivity$node[changed]),
    query = c(cells$query[low_coverage], cells$query[low_quality],
              sensitivity$query[changed]),
    warning = c(rep("coverage_below_declared_minimum", sum(low_coverage)),
                rep("sink_territory_above_declared_maximum", sum(low_quality)),
                rep("group_composition_changes_materially", sum(changed))),
    stringsAsFactors = FALSE
  )
  value <- structure(list(
    contract_version = "population-diagnostics-v1",
    parent_result_id = x$scientific_plan_id,
    parent_population_plan_id = x$receipt$population_plan_id,
    conditioning = x$coverage$conditioning,
    thresholds = list(minimum_coverage = minimum_coverage,
      minimum_transport_quality = minimum_transport_quality,
      material_change = material_change),
    subject_provenance = subject_provenance,
    cells = cells, associations = associations,
    sensitivity = sensitivity, warnings = warnings,
    interpretation = paste0(
      "Descriptive diagnostics conditional on realized transport and coverage; ",
      "threshold sensitivity does not replace the primary estimand."
    )
  ), class = "effect_population_diagnostics")
  value
}

#' @export
print.effect_population_diagnostics <- function(x, ...) {
  cat("<effect_population_diagnostics>\n")
  cat("  cells: ", nrow(x$cells), "; subjects: ",
      nrow(x$subject_provenance), "\n", sep = "")
  cat("  warnings: ", nrow(x$warnings), "\n", sep = "")
  cat("  reading: ", x$interpretation, "\n", sep = "")
  invisible(x)
}
