# Numerical reproducibility contract ----------------------------------------

#' Define effectagram's numerical reproducibility contract
#'
#' The contract separates three claims. Reordering completion of fixed tasks is
#' bitwise invariant because reduction order is canonical. Changing feature
#' blocks or numerical platform is guaranteed only within the declared combined
#' absolute and relative tolerance. Bitwise equality across block partitions or
#' platforms is deliberately not promised.
#'
#' @param atol,rtol Nonnegative finite absolute and relative tolerances.
#' @return A declarative numerical contract.
#' @export
numerical_contract <- function(atol = 1e-12, rtol = 1e-10) {
  if (!is.numeric(atol) || length(atol) != 1L || is.na(atol) ||
      !is.finite(atol) || atol < 0) {
    stop("`atol` must be one nonnegative finite number.", call. = FALSE)
  }
  if (!is.numeric(rtol) || length(rtol) != 1L || is.na(rtol) ||
      !is.finite(rtol) || rtol < 0) {
    stop("`rtol` must be one nonnegative finite number.", call. = FALSE)
  }
  structure(
    list(
      atol = atol,
      rtol = rtol,
      scheduling = list(
        guarantee = "bitwise",
        condition = "same tasks and canonical reduction order"
      ),
      block_partition = list(
        guarantee = "tolerance",
        condition = "same estimand and precision"
      ),
      cross_platform = list(
        guarantee = "tolerance",
        condition = "same estimand, precision, and supported kernels"
      ),
      bitwise_across_blocking = FALSE,
      bitwise_across_platforms = FALSE
    ),
    class = "effect_numerical_contract"
  )
}

#' Assess results under a named numerical guarantee
#'
#' Tolerance agreement uses
#' `abs(x - y) <= atol + rtol * max(abs(x), abs(y))` elementwise.
#'
#' @param x,y Numeric objects with identical dimensions.
#' @param guarantee One of `scheduling`, `block_partition`, or `cross_platform`.
#' @param contract An `effect_numerical_contract`.
#' @return A diagnostic `effect_numeric_agreement` value.
#' @export
numerical_agreement <- function(x, y,
                                guarantee = c("scheduling", "block_partition",
                                              "cross_platform"),
                                contract = numerical_contract()) {
  guarantee <- match.arg(guarantee)
  if (!inherits(contract, "effect_numerical_contract")) {
    stop("`contract` must be an effectagram numerical contract.", call. = FALSE)
  }
  if (!is.numeric(x) || !is.numeric(y) || !identical(dim(x), dim(y)) ||
      length(x) != length(y) || any(!is.finite(x)) || any(!is.finite(y))) {
    stop("`x` and `y` must be finite numeric objects with identical dimensions.",
      call. = FALSE)
  }

  difference <- abs(x - y)
  scale <- pmax(abs(x), abs(y))
  allowed <- contract$atol + contract$rtol * scale
  passed <- if (guarantee == "scheduling") identical(x, y) else all(difference <= allowed)

  structure(
    list(
      passed = passed,
      guarantee = guarantee,
      comparison = if (guarantee == "scheduling") "bitwise" else "tolerance",
      max_absolute_error = if (length(difference)) max(difference) else 0,
      max_allowed_error = if (length(allowed)) max(allowed) else contract$atol,
      atol = contract$atol,
      rtol = contract$rtol
    ),
    class = "effect_numeric_agreement"
  )
}

.canonical_reduce <- function(values, task_id) {
  if (!is.list(values) || length(values) < 1L || length(values) != length(task_id)) {
    stop("`values` and `task_id` must describe at least one result per task.",
      call. = FALSE)
  }
  if (anyNA(task_id) || anyDuplicated(task_id)) {
    stop("`task_id` values must be unique and non-missing.", call. = FALSE)
  }
  first_dim <- dim(values[[1L]])
  valid <- vapply(values, function(x) {
    is.numeric(x) && identical(dim(x), first_dim) && all(is.finite(x))
  }, logical(1))
  if (!all(valid)) {
    stop("Every task result must be finite numeric data with identical dimensions.",
      call. = FALSE)
  }

  ordered <- values[order(task_id, method = "radix")]
  Reduce(`+`, ordered)
}
