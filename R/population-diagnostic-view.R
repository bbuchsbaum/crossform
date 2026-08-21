# Bound effect and population-support views ----------------------------------

#' Bind population effects to coverage and transport diagnostics
#'
#' @param effect_view An `effect_population_component_view` or
#'   `effect_population_scale_profile`.
#' @param diagnostics An `effect_population_diagnostics` for the total result
#'   underlying the same decomposition.
#' @param node,query Optional synchronized selections.
#' @return An `effect_population_diagnostic_view` with exact-key effect and
#'   support rows, warnings, filter provenance, and separate diagnostic axes.
#' @export
population_diagnostic_view <- function(effect_view, diagnostics,
                                       node = NULL, query = NULL) {
  if (!inherits(effect_view, "effect_population_component_view")) {
    .input_error("`effect_view` must be a population component or scale view.",
      arg = "effect_view", received = .msg_value(effect_view),
      expected = "an `effect_population_component_view`")
  }
  if (!inherits(diagnostics, "effect_population_diagnostics")) {
    .input_error("`diagnostics` must be returned by `population_diagnostics()`.",
      arg = "diagnostics", received = .msg_value(diagnostics),
      expected = "an `effect_population_diagnostics`")
  }
  if (!identical(effect_view$parent_population_plan_id,
                 diagnostics$parent_population_plan_id)) {
    .population_decomposition_refusal("effect_and_diagnostic_plans_differ")
  }
  available_nodes <- unique(effect_view$data$node)
  available_queries <- unique(effect_view$data$query)
  if (is.null(node)) node <- available_nodes
  if (is.null(query)) query <- available_queries
  if (!.is_strings(node, unique = TRUE) || !all(node %in% available_nodes)) {
    .input_error("`node` contains an unavailable group node.",
      arg = "node", received = .msg_value(node),
      expected = paste0("names from ", .msg_names(available_nodes)))
  }
  if (!.is_strings(query, unique = TRUE) ||
      !all(query %in% available_queries)) {
    .input_error("`query` contains an unavailable query.",
      arg = "query", received = .msg_value(query),
      expected = paste0("names from ", .msg_names(available_queries)))
  }
  select_key <- function(data) data$node %in% node & data$query %in% query
  effect <- effect_view$data[select_key(effect_view$data), , drop = FALSE]
  support <- diagnostics$cells[select_key(diagnostics$cells), , drop = FALSE]
  warnings <- diagnostics$warnings[select_key(diagnostics$warnings), ,
                                   drop = FALSE]
  sensitivity <- diagnostics$sensitivity[
    select_key(diagnostics$sensitivity), , drop = FALSE
  ]
  key <- paste(support$node, support$query, sep = "::")
  effect_key <- paste(effect$node, effect$query, sep = "::")
  at <- match(effect_key, key)
  if (anyNA(at)) {
    .contract_error("Every displayed effect must have one exact support row.")
  }
  effect$coverage_fraction <- support$coverage_fraction[at]
  effect$n <- support$n[at]
  effect$mass_n_eff <- support$mass_n_eff[at]
  effect$mean_sink_territory <- support$mean_sink_territory[at]
  effect$mean_transport_quality <- support$mean_transport_quality[at]
  effect$support_key <- key[at]
  panel <- function(metric, label, units, limits) list(
    metric = metric, label = label, units = units, limits = limits,
    scale_identity = paste0("separate_", metric, "_axis")
  )
  structure(list(
    contract_version = "population-diagnostic-view-v1",
    parent_population_plan_id = effect_view$parent_population_plan_id,
    effect = effect, support = support, sensitivity = sensitivity,
    warnings = warnings,
    panels = list(
      coverage = panel("coverage_fraction", "Coverage", "planned-subject proportion",
                       c(0, 1)),
      effective_n = panel("mass_n_eff", "Mass effective N", "subjects", c(0, NA)),
      sink = panel("mean_sink_territory", "Sink territory",
                   "native-territory proportion", c(0, 1)),
      transport_quality = panel("mean_transport_quality", "Retained territory",
                                "native-territory proportion", c(0, 1))
    ),
    filters = list(
      node = node, query = query,
      effect_rows_before = nrow(effect_view$data),
      effect_rows_after = nrow(effect),
      support_rows_before = nrow(diagnostics$cells),
      support_rows_after = nrow(support),
      operation = "synchronized_exact_key_selection"
    ),
    thresholds = diagnostics$thresholds,
    interpretation = paste0(
      "Effect and support panels share node/query identity. Diagnostic axes ",
      "have separate units and scales; filters do not redefine the primary fit."
    )
  ), class = "effect_population_diagnostic_view")
}

#' @export
print.effect_population_diagnostic_view <- function(x, ...) {
  cat("<effect_population_diagnostic_view>\n")
  cat("  effect rows: ", nrow(x$effect), "; support rows: ",
      nrow(x$support), "\n", sep = "")
  cat("  warnings: ", nrow(x$warnings), "\n", sep = "")
  cat("  filter: ", x$filters$operation, "\n", sep = "")
  cat("  reading: ", x$interpretation, "\n", sep = "")
  invisible(x)
}
