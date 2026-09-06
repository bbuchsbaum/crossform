# Population wild-bootstrap inference ---------------------------------------
#
# The bootstrap unit is the participant. One subject-by-replicate weight
# matrix is generated before any cell is visited and the same row is reused at
# every node and query in that replicate. That is the structural guarantee:
# within-subject dependence across the jointly analysed result is preserved by
# construction, including when D1 gives different cells different subject
# subsets.

.population_bootstrap_distributions <- c("rademacher", "mammen")
.population_bootstrap_min_subjects <- 6L

.population_bootstrap_contrast <- function(x, contrast) {
  terms <- dimnames(x$coefficients)[[3L]]
  if (.is_string(contrast)) {
    if (!contrast %in% terms) {
      .input_error(sprintf(
        "`contrast` must name a population coefficient: %s.",
        .msg_names(terms)
      ), arg = "contrast", received = .msg_value(contrast),
        expected = paste0("one of ", .msg_names(terms)))
    }
    value <- stats::setNames(numeric(length(terms)), terms)
    value[[contrast]] <- 1
    return(value)
  }
  if (!is.numeric(contrast) || length(contrast) != length(terms) ||
      any(!is.finite(contrast))) {
    .input_error(sprintf(paste0(
      "`contrast` must be one coefficient name or %d finite weights over ",
      "the population model terms."
    ), length(terms)), arg = "contrast", received = .msg_value(contrast),
      expected = paste0("one name, or weights over ", .msg_names(terms)))
  }
  if (!is.null(names(contrast))) {
    if (!.is_strings(names(contrast), unique = TRUE) ||
        !setequal(names(contrast), terms)) {
      .input_error(paste0(
        "Named `contrast` weights must identify every population model term ",
        "exactly once."
      ), arg = "contrast", received = .msg_value(contrast),
        expected = paste0("weights named ", .msg_names(terms)))
    }
    contrast <- contrast[terms]
  }
  if (!any(contrast != 0)) {
    .input_error("`contrast` must contain at least one nonzero weight.",
      arg = "contrast", received = .msg_value(contrast),
      expected = "a nonzero coefficient contrast")
  }
  stats::setNames(as.numeric(contrast), terms)
}

.population_bootstrap_weights <- function(subjects, replicates, distribution,
                                          seed) {
  global <- .GlobalEnv
  had_seed <- exists(".Random.seed", envir = global, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = global, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = global)
    } else if (exists(".Random.seed", envir = global, inherits = FALSE)) {
      rm(".Random.seed", envir = global)
    }
  }, add = TRUE)
  set.seed(seed)
  count <- length(subjects) * replicates
  draw <- if (identical(distribution, "rademacher")) {
    ifelse(stats::runif(count) < 0.5, -1, 1)
  } else {
    root <- sqrt(5)
    lower <- -(root - 1) / 2
    upper <- (root + 1) / 2
    probability_lower <- (root + 1) / (2 * root)
    ifelse(stats::runif(count) < probability_lower, lower, upper)
  }
  matrix(draw, length(subjects), replicates,
    dimnames = list(subject = subjects,
      replicate = paste0("r", seq_len(replicates))))
}

.population_bootstrap_hc3 <- function(design, response, contrast,
                                      leverage_tolerance) {
  factorization <- qr(design)
  if (factorization$rank < ncol(design)) {
    return(list(ok = FALSE, reason = "rank_deficient_design"))
  }
  bread <- .population_bread(factorization, ncol(design))
  coefficient <- as.numeric(bread %*% crossprod(design, response))
  fitted <- as.numeric(design %*% coefficient)
  residual <- response - fitted
  leverage <- rowSums((design %*% bread) * design)
  leverage[leverage < 0 & leverage > -leverage_tolerance] <- 0
  leverage[leverage > 1 & leverage < 1 + leverage_tolerance] <- 1
  if (any(1 - leverage <= leverage_tolerance)) {
    return(list(ok = FALSE, reason = "leverage_near_one"))
  }
  adjusted <- residual / (1 - leverage)
  meat <- crossprod(design, design * adjusted^2)
  covariance <- bread %*% meat %*% bread
  covariance <- (covariance + t(covariance)) / 2
  variance <- as.numeric(crossprod(contrast, covariance %*% contrast))
  scale <- sqrt(max(variance, 0))
  if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps) *
      max(1, abs(as.numeric(crossprod(contrast, coefficient))))) {
    return(list(
      ok = FALSE, reason = "nonpositive_standard_error",
      coefficient = coefficient, fitted = fitted, residual = residual,
      leverage = leverage, adjusted_residual = adjusted,
      bread = bread, covariance = covariance, se = scale
    ))
  }
  list(
    ok = TRUE, reason = NA_character_, coefficient = coefficient,
    estimate = as.numeric(crossprod(contrast, coefficient)),
    fitted = fitted, residual = residual, leverage = leverage,
    adjusted_residual = adjusted, bread = bread, covariance = covariance,
    se = scale
  )
}

.population_bootstrap_id <- function(parent, contrast, null, replicates,
                                     distribution, seed, weights,
                                     leverage_tolerance, conditioning) {
  .sha256_signature(list(
    schema_version = 1L,
    contract_version = "population-estimand-v1",
    role = "population_null_imposed_wild_bootstrap",
    parent = parent,
    contrast = contrast,
    null = null,
    replicates = replicates,
    weight_distribution = distribution,
    seed = seed,
    weight_signature = .sha256_signature(weights),
    leverage_tolerance = leverage_tolerance,
    conditioning = conditioning
  ), "population-sha256:")
}

.population_bootstrap_fields <- c(
  "estimator", "studentization", "weight_distribution", "contrast", "null",
  "alternative", "replicates", "seed", "rng", "leverage_tolerance",
  "weights", "weight_signature",
  "observed_estimate", "observed_se", "observed_t", "replicate_t",
  "replicate_success", "replicate_failure_reason", "successful_replicates",
  "failed_replicates", "p_value", "monte_carlo_se", "critical_value",
  "reject", "level", "status", "reason", "coverage", "index", "queries",
  "conditioning", "receipt", "parent_result_id", "scientific_plan_id"
)

.validate_population_bootstrap <- function(x) {
  if (!inherits(x, "effect_population_wild_bootstrap") ||
      !.sealed_fields(x, "effect_population_wild_bootstrap",
        .population_bootstrap_fields) ||
      !identical(x$estimator, "null_imposed_wild_bootstrap") ||
      !identical(x$studentization, "HC3") ||
      !.is_string(x$weight_distribution) ||
      !x$weight_distribution %in% .population_bootstrap_distributions ||
      !.is_string(x$alternative) || !identical(x$alternative, "two_sided") ||
      !is.matrix(x$weights) || !is.array(x$replicate_t) ||
      !identical(dim(x$replicate_t), dim(x$replicate_success)) ||
      !identical(dim(x$replicate_t), dim(x$replicate_failure_reason)) ||
      !identical(dim(x$replicate_t)[[3L]], x$replicates) ||
      !identical(nrow(x$weights), length(x$coverage$planned_subjects)) ||
      !identical(ncol(x$weights), x$replicates) ||
      !identical(rownames(x$weights), x$coverage$planned_subjects) ||
      !identical(x$weight_signature, .sha256_signature(x$weights)) ||
      !identical(dim(x$observed_estimate), dim(x$p_value)) ||
      !identical(dim(x$p_value), dim(x$status)) ||
      !identical(dim(x$p_value), dim(x$reason)) ||
      !.is_number(x$leverage_tolerance) ||
      x$leverage_tolerance <= 0 || x$leverage_tolerance >= 1 ||
      !.is_string(x$parent_result_id) ||
      !.strong_sha256(sub("^population-", "", x$scientific_plan_id))) {
    .input_error("Population wild-bootstrap fields are missing or noncanonical.")
  }
  .validate_population_conditioning(
    x$conditioning, x$coverage$planned_subjects
  )
  if (!identical(x$conditioning, x$coverage$conditioning)) {
    .contract_error(paste0(
      "Wild-bootstrap transport conditioning must be inherited unchanged ",
      "from the population result."
    ))
  }
  counted <- apply(x$replicate_success, c(1L, 2L), sum)
  storage.mode(counted) <- "integer"
  if (!identical(counted, x$successful_replicates) ||
      !identical(x$failed_replicates, x$replicates - counted) ||
      any(x$replicate_success & !is.finite(x$replicate_t)) ||
      any(!x$replicate_success & is.finite(x$replicate_t)) ||
      any(x$replicate_success & !is.na(x$replicate_failure_reason)) ||
      any(!x$replicate_success & is.na(x$replicate_failure_reason)) ||
      any(x$status == "estimated" & !is.finite(x$p_value)) ||
      any(x$status != "estimated" & is.finite(x$p_value)) ||
      any(x$p_value[is.finite(x$p_value)] <= 0 |
        x$p_value[is.finite(x$p_value)] > 1)) {
    .contract_error("Population wild-bootstrap replicate accounting is inconsistent.")
  }
  expected_mcse <- sqrt(x$p_value * (1 - x$p_value) /
    (1 + x$successful_replicates))
  finite_mcse <- is.finite(expected_mcse)
  if (!identical(is.na(x$monte_carlo_se), is.na(expected_mcse)) ||
      any(abs(x$monte_carlo_se[finite_mcse] -
        expected_mcse[finite_mcse]) > 0)) {
    .contract_error("Population wild-bootstrap Monte Carlo error is inconsistent.")
  }
  expected_id <- .population_bootstrap_id(
    x$parent_result_id, x$contrast, x$null, x$replicates,
    x$weight_distribution, x$seed, x$weights, x$leverage_tolerance,
    x$conditioning
  )
  if (!identical(x$scientific_plan_id, expected_id)) {
    .contract_error("Population wild-bootstrap identity is inconsistent.")
  }
  invisible(x)
}

#' Null-imposed subject-level wild bootstrap for population coefficients
#'
#' `population_wild_bootstrap()` tests one explicit linear contrast of the
#' population-model coefficients against one explicit null at every admitted
#' node and query. It uses a null-imposed, HC3-studentized wild bootstrap. One
#' random weight is drawn per planned subject and replicate and reused across
#' **all** nodes and queries, preserving the within-subject dependence of the
#' jointly analysed result.
#'
#' @section Algorithm and null:
#' For coefficient contrast \eqn{c}{c} and null \eqn{c'\beta = b_0}{c' beta =
#' b0}, each cell is first fit under the linear restriction. Bootstrap response
#' \eqn{y^* = X\hat\beta_R + w_i e_{Ri}/(1-h_i)}{y* = X beta_R + w_i e_Ri /
#' (1-h_i)} uses the restricted residual and its HC3 leverage adjustment. Each
#' replicate is refit without the restriction and studentized by its own HC3
#' covariance. The two-sided p-value uses the plus-one rule, so it is never
#' zero; `$monte_carlo_se` reports its finite-replication uncertainty.
#'
#' `"rademacher"` weights are symmetric \eqn{-1,+1}{-1,+1} draws with equal
#' probability. `"mammen"` uses the two-point, mean-zero, variance-one Mammen
#' distribution, whose third moment also equals one; it is provided for
#' asymmetric-error sensitivity analyses, not as an automatic improvement.
#'
#' @section Coverage, failures, and conditioning:
#' The exact subject set from `$coverage` is used separately at every cell,
#' while the weights come from one shared planned-subject matrix. Cells with
#' fewer than six available subjects, rank deficiency, no residual df,
#' leverage within `leverage_tolerance` of one, or a nonpositive observed HC3
#' standard error are refused with a status and reason. Replicate failures are
#' retained in `$replicate_failure_reason`; a cell is reported only when at
#' least 90 percent and at least 20 replicates succeed. All inference is
#' conditional on the realized transport and coverage. It does not propagate
#' uncertainty from estimating either. Cross-fitting can limit circular reuse
#' of responses, but does not make this bootstrap marginal over transport
#' estimation or fold assignment.
#'
#' The returned `$weights`, `$seed`, `$rng`, `$weight_signature`, replicate
#' statistics, success flags, failure reasons, contrast, null and parent
#' receipt provide full resampling provenance. The caller's random-number state
#' is restored on exit.
#'
#' @param x An `effect_population_result` from [estimate_population()].
#' @param contrast One population coefficient name, or a numeric vector with
#'   one weight per model term. Named weights are aligned by term name.
#' @param null Finite scalar value of the tested coefficient contrast.
#' @param replicates Number of bootstrap replicates, at least 99.
#' @param seed Required nonnegative integer random seed.
#' @param weights Wild-weight distribution: `"rademacher"` or `"mammen"`.
#' @param level Nominal two-sided test level used for the absolute bootstrap-t
#'   critical value and `$reject` diagnostic.
#' @param leverage_tolerance Positive tolerance below which \eqn{1-h_i}{1-h_i}
#'   is treated as numerically zero.
#' @return An `effect_population_wild_bootstrap` carrying the explicit null,
#'   observed and replicate HC3-t statistics, shared subject weights,
#'   replicate-level failures, plus-one p-values, critical values, rejection
#'   flags, Monte Carlo SEs, exact coverage, conditioning and provenance.
#' @family population transports
#' @seealso [population_uncertainty()] for analytic classical and HC3
#'   covariance without resampling.
#' @export
population_wild_bootstrap <- function(
    x, contrast, null = 0, replicates = 999L, seed = NULL,
    weights = c("rademacher", "mammen"), level = 0.95,
    leverage_tolerance = 1e-8) {
  .validate_population_result(x)
  if (!identical(x$basis, "query_bank")) {
    .population_uncertainty_refuse_basis()
  }
  contrast <- .population_bootstrap_contrast(x, contrast)
  null <- .check_number(null, "null", what = "one finite tested null")
  replicates <- .check_count(replicates, "replicates", min = 99L,
    what = "an integer number of bootstrap replicates, at least 99")
  if (is.null(seed) || !.is_number(seed) || seed %% 1 != 0 || seed < 0 ||
      seed > .Machine$integer.max) {
    .input_error("`seed` must be one explicit nonnegative integer.",
      arg = "seed", received = .msg_value(seed),
      expected = "one nonnegative integer")
  }
  seed <- as.integer(seed)
  weights <- match.arg(weights)
  level <- .check_number(level, "level",
    what = "one nominal two-sided test level in (0, 1)")
  if (level <= 0 || level >= 1) {
    .input_error("`level` must lie strictly between zero and one.", arg = "level")
  }
  leverage_tolerance <- .check_number(
    leverage_tolerance, "leverage_tolerance",
    what = "one positive finite leverage-complement tolerance"
  )
  if (leverage_tolerance <= 0 || leverage_tolerance >= 1) {
    .input_error(
      "`leverage_tolerance` must lie strictly between zero and one.",
      arg = "leverage_tolerance")
  }

  subjects <- x$coverage$planned_subjects
  weight_matrix <- .population_bootstrap_weights(
    subjects, replicates, weights, seed
  )
  nodes <- dim(x$values)[[1L]]
  queries <- dim(x$values)[[2L]]
  cell_names <- dimnames(x$values)[1:2]
  cell_matrix <- function(mode = "double", fill = NA) {
    matrix(fill, nodes, queries, dimnames = cell_names)
  }
  observed_estimate <- observed_se <- observed_t <- cell_matrix()
  p_value <- monte_carlo_se <- critical_value <- cell_matrix()
  reject <- cell_matrix("logical", NA)
  status <- cell_matrix("character", "refused")
  reason <- cell_matrix("character", NA_character_)
  replicate_t <- array(NA_real_, c(nodes, queries, replicates),
    dimnames = c(cell_names, list(replicate = colnames(weight_matrix))))
  replicate_success <- array(FALSE, dim(replicate_t), dimnames = dimnames(replicate_t))
  replicate_failure_reason <- array(NA_character_, dim(replicate_t),
    dimnames = dimnames(replicate_t))
  design_all <- x$uncertainty$between$design

  for (node in seq_len(nodes)) {
    for (query in seq_len(queries)) {
      available <- x$coverage$availability[node, query, ]
      n <- sum(available)
      if (n < .population_bootstrap_min_subjects) {
        reason[node, query] <- "insufficient_subjects_for_wild_bootstrap"
        replicate_failure_reason[node, query, ] <- reason[node, query]
        next
      }
      design <- design_all[available, , drop = FALSE]
      response <- x$values[node, query, available]
      observed <- .population_bootstrap_hc3(
        design, response, contrast, leverage_tolerance
      )
      if (!isTRUE(observed$ok)) {
        reason[node, query] <- observed$reason
        replicate_failure_reason[node, query, ] <- observed$reason
        next
      }
      observed_estimate[node, query] <- observed$estimate
      observed_se[node, query] <- observed$se
      observed_t[node, query] <- (observed$estimate - null) / observed$se

      denominator <- as.numeric(crossprod(
        contrast, observed$bread %*% contrast
      ))
      if (!is.finite(denominator) || denominator <= 0) {
        reason[node, query] <- "contrast_not_estimable"
        replicate_failure_reason[node, query, ] <- reason[node, query]
        next
      }
      restricted_coefficient <- observed$coefficient -
        as.numeric(observed$bread %*% contrast) *
        (observed$estimate - null) / denominator
      restricted_fitted <- as.numeric(design %*% restricted_coefficient)
      restricted_residual <- response - restricted_fitted
      restricted_adjusted <- restricted_residual / (1 - observed$leverage)
      local_weights <- weight_matrix[available, , drop = FALSE]

      for (replicate in seq_len(replicates)) {
        bootstrap_response <- restricted_fitted +
          restricted_adjusted * local_weights[, replicate]
        bootstrap <- .population_bootstrap_hc3(
          design, bootstrap_response, contrast, leverage_tolerance
        )
        if (!isTRUE(bootstrap$ok)) {
          replicate_failure_reason[node, query, replicate] <- bootstrap$reason
          next
        }
        replicate_t[node, query, replicate] <-
          (bootstrap$estimate - null) / bootstrap$se
        replicate_success[node, query, replicate] <- TRUE
      }
      successful <- sum(replicate_success[node, query, ])
      required <- max(20L, ceiling(0.9 * replicates))
      if (successful < required) {
        reason[node, query] <- "insufficient_successful_replicates"
        next
      }
      draws <- replicate_t[node, query, replicate_success[node, query, ]]
      exceed <- sum(abs(draws) >= abs(observed_t[node, query]))
      p_value[node, query] <- (1 + exceed) / (1 + successful)
      monte_carlo_se[node, query] <- sqrt(
        p_value[node, query] * (1 - p_value[node, query]) /
          (1 + successful)
      )
      critical_value[node, query] <- as.numeric(stats::quantile(
        abs(draws), probs = level, names = FALSE, type = 8
      ))
      reject[node, query] <- abs(observed_t[node, query]) >
        critical_value[node, query]
      status[node, query] <- "estimated"
    }
  }
  successful_replicates <- apply(replicate_success, c(1L, 2L), sum)
  storage.mode(successful_replicates) <- "integer"
  weight_signature <- .sha256_signature(weight_matrix)
  value <- structure(list(
    estimator = "null_imposed_wild_bootstrap",
    studentization = "HC3",
    weight_distribution = weights,
    contrast = contrast,
    null = as.numeric(null),
    alternative = "two_sided",
    replicates = replicates,
    seed = seed,
    rng = list(kind = RNGkind(), r_version = as.character(getRversion())),
    leverage_tolerance = as.numeric(leverage_tolerance),
    weights = weight_matrix,
    weight_signature = weight_signature,
    observed_estimate = observed_estimate,
    observed_se = observed_se,
    observed_t = observed_t,
    replicate_t = replicate_t,
    replicate_success = replicate_success,
    replicate_failure_reason = replicate_failure_reason,
    successful_replicates = successful_replicates,
    failed_replicates = replicates - successful_replicates,
    p_value = p_value,
    monte_carlo_se = monte_carlo_se,
    critical_value = critical_value,
    reject = reject,
    level = as.numeric(level),
    status = status,
    reason = reason,
    coverage = x$coverage,
    index = x$index,
    queries = x$queries,
    conditioning = x$coverage$conditioning,
    receipt = x$receipt,
    parent_result_id = x$scientific_plan_id,
    scientific_plan_id = .population_bootstrap_id(
      x$scientific_plan_id, contrast, null, replicates, weights, seed,
      weight_matrix, leverage_tolerance, x$coverage$conditioning
    )
  ), class = "effect_population_wild_bootstrap")
  .validate_population_bootstrap(value)
  value
}
