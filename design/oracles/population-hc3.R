# Independent explicit-matrix oracle for population HC3 covariance.
#
# This file is deliberately outside R/ and uses no crossform implementation
# helper. Production uses a pivoted QR and a cellwise execution record. The
# oracle uses singular values to admit the design, normal equations solved by
# base LAPACK, and the HC3 definition written out as dense matrices.

population_hc3_oracle <- function(design, response, terms = NULL,
                                  rank_tolerance = sqrt(.Machine$double.eps),
                                  leverage_tolerance = 1e-8) {
  if (!is.matrix(design) || !is.numeric(design) || any(!is.finite(design)) ||
      !is.numeric(response) || length(response) != nrow(design) ||
      any(!is.finite(response))) {
    stop("HC3 oracle requires a finite numeric design matrix and response.",
      call. = FALSE)
  }
  if (!nrow(design) || !ncol(design)) {
    stop("HC3 oracle requires a nonempty design.", call. = FALSE)
  }
  if (nrow(design) <= ncol(design)) {
    stop("HC3 oracle requires positive residual degrees of freedom.",
      call. = FALSE)
  }
  singular <- svd(design, nu = 0L, nv = 0L)$d
  cutoff <- rank_tolerance * max(singular)
  if (min(singular) <= cutoff) {
    stop("HC3 oracle refuses a rank-deficient or numerically singular design.",
      call. = FALSE)
  }
  columns <- colnames(design)
  if (is.null(columns)) columns <- paste0("term", seq_len(ncol(design)))
  if (is.null(terms)) terms <- columns
  if (!is.character(terms) || !length(terms) || anyDuplicated(terms) ||
      !all(terms %in% columns)) {
    stop("HC3 oracle terms must uniquely select design columns.",
      call. = FALSE)
  }

  crossproduct <- crossprod(design)
  bread <- solve(crossproduct, diag(ncol(design)))
  dimnames(bread) <- list(columns, columns)
  coefficient <- as.numeric(bread %*% crossprod(design, response))
  names(coefficient) <- columns
  fitted <- as.numeric(design %*% coefficient)
  residual <- response - fitted
  hat_matrix <- design %*% bread %*% t(design)
  leverage <- diag(hat_matrix)
  if (any(1 - leverage <= leverage_tolerance)) {
    stop("HC3 oracle refuses leverage whose complement is numerically zero.",
      call. = FALSE)
  }
  adjusted_residual <- residual / (1 - leverage)
  meat <- crossprod(design, design * adjusted_residual^2)
  covariance <- bread %*% meat %*% bread
  covariance <- (covariance + t(covariance)) / 2
  dimnames(covariance) <- list(columns, columns)
  selected <- match(terms, columns)

  list(
    estimator = "HC3_explicit_matrix_oracle_v1",
    assumptions = c(
      "independent_subjects", "fixed_group_design",
      "heteroskedasticity_robust_sandwich",
      "finite_sample_leverage_adjustment"
    ),
    coefficient = coefficient[terms],
    fitted = fitted,
    residual = residual,
    leverage = leverage,
    adjusted_residual = adjusted_residual,
    bread = bread,
    meat = meat,
    covariance = covariance[selected, selected, drop = FALSE],
    standard_error = sqrt(pmax(diag(
      covariance[selected, selected, drop = FALSE]
    ), 0)),
    n = nrow(design),
    rank = length(singular),
    residual_df = nrow(design) - ncol(design),
    condition_number = max(singular) / min(singular),
    leverage_tolerance = leverage_tolerance
  )
}
if (sys.nframe() == 0L) {
  x <- c(-1, 0, 1, 2, 4)
  design <- cbind("(Intercept)" = 1, x = x)
  response <- c(1, 2, 2, 4, 9)
  value <- population_hc3_oracle(design, response)
  stopifnot(
    identical(value$estimator, "HC3_explicit_matrix_oracle_v1"),
    max(abs(value$leverage - diag(
      design %*% solve(crossprod(design)) %*% t(design)
    ))) < 1e-14,
    all(is.finite(value$covariance)),
    inherits(try(population_hc3_oracle(
      cbind(1, x, 2 * x), response
    ), silent = TRUE), "try-error")
  )
  message("population-hc3-oracle PASS")
}
