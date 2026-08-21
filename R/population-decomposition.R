# Population coefficient decomposition --------------------------------------

.population_decomposition_refusal <- function(reason) {
  .capability_refusal(paste0(
    "Population component coefficients are comparable only when total, ",
    "coherent, and configuration results use the same population plan, ",
    "query bank, model, subject set at every cell, and additive normalization."
  ), capability = "population_component_decomposition",
    namespace = "population_decomposition", reasons = reason,
    remedies = paste0(
      "Estimate all three components from one `effect_population_plan` with ",
      "the same queries and `normalization = \"none\"`."
    ))
}

.population_decomposition_cross_covariance <- function(coherent,
                                                        configuration,
                                                        estimator) {
  nodes <- dim(coherent$coefficients)[[1L]]
  queries <- dim(coherent$coefficients)[[2L]]
  terms <- dimnames(coherent$coefficients)[[3L]]
  components <- c("total", "coherent", "configuration")
  labels <- list(node = dimnames(coherent$coefficients)[[1L]],
    query = dimnames(coherent$coefficients)[[2L]], term = terms,
    component = components, term2 = terms, component2 = components)
  output <- array(NA_real_, c(nodes, queries, length(terms), 3L,
                              length(terms), 3L), dimnames = labels)
  design <- coherent$uncertainty$between$design
  for (node in seq_len(nodes)) {
    for (query in seq_len(queries)) {
      available <- coherent$coverage$availability[node, query, ]
      X <- design[available, , drop = FALSE]
      if (nrow(X) <= ncol(X) || qr(X)$rank < ncol(X)) next
      bread <- solve(crossprod(X))
      ec <- coherent$residuals[node, query, available]
      eg <- configuration$residuals[node, query, available]
      if (identical(estimator, "HC3")) {
        leverage <- rowSums((X %*% bread) * X)
        if (any(1 - leverage <= sqrt(.Machine$double.eps))) next
        ec <- ec / (1 - leverage)
        eg <- eg / (1 - leverage)
        covariance <- function(left, right) {
          bread %*% crossprod(X, X * (left * right)) %*% bread
        }
      } else {
        covariance <- function(left, right) {
          sum(left * right) / (nrow(X) - ncol(X)) * bread
        }
      }
      cc <- covariance(ec, ec)
      gg <- covariance(eg, eg)
      cg <- covariance(ec, eg)
      blocks <- list(
        total = list(total = cc + gg + 2 * cg,
          coherent = cc + cg, configuration = gg + cg),
        coherent = list(total = cc + cg, coherent = cc,
          configuration = cg),
        configuration = list(total = gg + cg, coherent = cg,
          configuration = gg)
      )
      for (left in components) for (right in components) {
        output[node, query, , left, , right] <- blocks[[left]][[right]]
      }
    }
  }
  output
}

#' Enforce the population coefficient decomposition law
#'
#' Checks `beta_total = beta_coherent + beta_configuration` for every
#' estimable population coefficient, node, and query. The three results must
#' come from the same plan, model, query bank, and cellwise subject sets.
#'
#' @param total,coherent,configuration Comparable `effect_population_result`
#'   objects for the named components.
#' @param estimator `"HC3"` or `"classical"` for the joint component
#'   coefficient covariance.
#' @param tolerance Absolute numerical tolerance for the conservation law.
#' @return An `effect_population_decomposition` carrying the coefficient gap,
#'   joint component covariance, uncertainty objects, and comparability receipt.
#' @export
population_decomposition <- function(total, coherent, configuration,
                                     estimator = c("HC3", "classical"),
                                     tolerance = 1e-10) {
  inputs <- list(total = total, coherent = coherent,
                 configuration = configuration)
  for (input in inputs) .validate_population_result(input)
  estimator <- match.arg(estimator)
  tolerance <- .check_number(tolerance, "tolerance", nonnegative = TRUE,
    what = "one nonnegative absolute conservation tolerance")
  if (!identical(unname(vapply(inputs, `[[`, character(1), "component")),
                 names(inputs))) {
    .population_decomposition_refusal("component_roles_do_not_match_arguments")
  }
  parents <- vapply(inputs, function(value)
    value$receipt$population_plan_id, character(1))
  if (length(unique(parents)) != 1L) {
    .population_decomposition_refusal("population_plans_differ")
  }
  if (!all(vapply(inputs, function(value)
      identical(value$normalization, "none"), logical(1)))) {
    .population_decomposition_refusal("componentwise_normalization_is_nonadditive")
  }
  reference <- total
  comparable <- vapply(inputs[-1L], function(value) {
    identical(value$queries, reference$queries) &&
      identical(value$coverage$subject_set_id,
                reference$coverage$subject_set_id) &&
      identical(value$coverage$coefficient_estimable,
                reference$coverage$coefficient_estimable) &&
      identical(value$uncertainty$between$design,
                reference$uncertainty$between$design) &&
      identical(dimnames(value$coefficients),
                dimnames(reference$coefficients))
  }, logical(1))
  if (!all(comparable)) {
    .population_decomposition_refusal("queries_models_or_subject_sets_differ")
  }
  coefficient_gap <- total$coefficients - coherent$coefficients -
    configuration$coefficients
  value_gap <- total$values - coherent$values - configuration$values
  estimable <- reference$coverage$coefficient_estimable
  max_coefficient_gap <- max(abs(coefficient_gap[estimable]), na.rm = TRUE)
  max_value_gap <- max(abs(value_gap), na.rm = TRUE)
  if (!is.finite(max_coefficient_gap) || max_coefficient_gap > tolerance ||
      !is.finite(max_value_gap) || max_value_gap > tolerance) {
    .contract_error(sprintf(paste0(
      "Population decomposition failed: coefficient gap %.3g and subject-value ",
      "gap %.3g exceed tolerance %.3g."
    ), max_coefficient_gap, max_value_gap, tolerance))
  }
  uncertainty <- lapply(inputs, population_uncertainty,
                        estimator = estimator)
  joint <- .population_decomposition_cross_covariance(
    coherent, configuration, estimator
  )
  direct_total <- uncertainty$total$between$covariance
  derived_total <- joint[, , , "total", , "total", drop = FALSE]
  dim(derived_total) <- dim(direct_total)
  dimnames(derived_total) <- dimnames(direct_total)
  covariance_gap <- max(abs(direct_total - derived_total), na.rm = TRUE)
  value <- structure(list(
    contract_version = "population-decomposition-v1",
    parent_population_plan_id = parents[[1L]], estimator = estimator,
    tolerance = tolerance, results = inputs,
    coefficients = lapply(inputs, `[[`, "coefficients"),
    coefficient_gap = coefficient_gap, max_coefficient_gap = max_coefficient_gap,
    max_subject_value_gap = max_value_gap,
    component_covariance = joint, uncertainty = uncertainty,
    direct_derived_total_covariance_gap = covariance_gap,
    coverage = reference$coverage, index = reference$index,
    queries = reference$queries, normalization = reference$normalization,
    law = "beta_total = beta_coherent + beta_configuration",
    interpretation = paste0(
      "Components are additive estimand ledgers under one model and subject set; ",
      "they are not separate biological mechanisms."
    )
  ), class = "effect_population_decomposition")
  value
}

#' @export
print.effect_population_decomposition <- function(x, ...) {
  cat("<effect_population_decomposition>\n")
  cat("  law: ", x$law, "\n", sep = "")
  cat("  max coefficient gap: ", format(x$max_coefficient_gap), "\n", sep = "")
  cat("  estimator: ", x$estimator, "\n", sep = "")
  cat("  reading: ", x$interpretation, "\n", sep = "")
  invisible(x)
}
