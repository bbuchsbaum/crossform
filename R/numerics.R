# Numerical reproducibility contract ----------------------------------------

#' Define crossform's numerical reproducibility contract
#'
#' The contract separates three claims. Reordering completion of fixed tasks is
#' bitwise invariant because reduction order is canonical. Changing feature
#' blocks or numerical platform is guaranteed only within the declared combined
#' absolute and relative tolerance. Bitwise equality across block partitions or
#' platforms is deliberately not promised.
#'
#' @param atol,rtol Nonnegative finite absolute and relative tolerances.
#' @return An `effect_numerical_contract`: the declared `$atol` and `$rtol`
#'   plus one guarantee record per claim (`$scheduling`, `$block_partition`,
#'   `$cross_platform`) and the two explicit non-promises
#'   `$bitwise_across_blocking` and `$bitwise_across_platforms`.
#' @seealso [numerical_agreement()], which tests two results against one of
#'   these named guarantees.
#' @family numerical contracts and receipts
#' @examples
#' # The default contract. Read the guarantee you intend to rely on rather
#' # than assuming bitwise equality everywhere.
#' contract <- numerical_contract()
#' contract$scheduling
#' contract$block_partition
#'
#' # Bitwise reproducibility across block partitions is deliberately not
#' # promised, so a tolerance is the only portable claim.
#' contract$bitwise_across_blocking
#'
#' # Tighten the tolerance when a downstream check needs it.
#' numerical_contract(atol = 1e-14, rtol = 1e-12)$rtol
#' @export
numerical_contract <- function(atol = 1e-12, rtol = 1e-10) {
  .check_number(atol, "atol", nonnegative = TRUE)
  .check_number(rtol, "rtol", nonnegative = TRUE)
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
#'   `scheduling` is checked bitwise; the other two are checked against the
#'   contract tolerance.
#' @param contract An `effect_numerical_contract`.
#' @return An `effect_numeric_agreement` with `$passed`, the `$guarantee` and
#'   `$comparison` used (`"bitwise"` or `"tolerance"`),
#'   `$max_absolute_error`, `$max_allowed_error`, and the `$atol`/`$rtol` in
#'   force.
#' @seealso [numerical_contract()] for the guarantees themselves.
#' @family numerical contracts and receipts
#' @examples
#' # Re-running the same fixed plan under a different feature-block size is a
#' # `block_partition` change: crossform promises tolerance, not bitwise
#' # equality.
#' baseline <- c(1, 2, 3) / 7
#' reblocked <- baseline + c(0, 1e-15, -2e-15)
#' agreement <- numerical_agreement(baseline, reblocked, "block_partition")
#' agreement$passed
#' agreement$max_absolute_error
#'
#' # The same pair fails the strict bitwise scheduling guarantee, which is
#' # reserved for reordering completion of identical tasks.
#' numerical_agreement(baseline, reblocked, "scheduling")$passed
#' numerical_agreement(baseline, baseline, "scheduling")$passed
#' @export
numerical_agreement <- function(x, y,
                                guarantee = c("scheduling", "block_partition",
                                              "cross_platform"),
                                contract = numerical_contract()) {
  guarantee <- match.arg(guarantee)
  contract <- .validate_numerical_contract(contract)
  if (!.is_finite_numeric(x) || !.is_finite_numeric(y) ||
      !identical(dim(x), dim(y)) || length(x) != length(y)) {
    .input_error(
      "`x` and `y` must be finite numeric objects with identical dimensions."
    )
  }

  difference <- abs(x - y)
  scale <- pmax(abs(x), abs(y))
  allowed <- contract$atol + contract$rtol * scale
  # `num.eq = FALSE` makes the comparison genuinely bitwise: it distinguishes
  # 0 from -0, which `==` (and identical()'s default) treats as equal.
  passed <- if (guarantee == "scheduling") {
    identical(x, y, num.eq = FALSE)
  } else {
    all(difference <= allowed)
  }

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

.validate_numerical_contract <- function(contract) {
  if (!inherits(contract, "effect_numerical_contract") || !is.list(contract)) {
    .input_error("`contract` must be a crossform numerical contract.")
  }
  rebuilt <- numerical_contract(contract$atol, contract$rtol)
  if (!identical(contract, rebuilt)) {
    .input_error("Numerical contract fields are inconsistent or noncanonical.")
  }
  rebuilt
}

.canonical_reduce <- function(values, task_id) {
  if (!is.list(values) || length(values) < 1L || length(values) != length(task_id)) {
    .input_error(
      "`values` and `task_id` must describe at least one result per task."
    )
  }
  if (anyNA(task_id) || anyDuplicated(task_id)) {
    .input_error("`task_id` values must be unique and non-missing.")
  }
  first_dim <- dim(values[[1L]])
  valid <- vapply(values, function(x) {
    is.numeric(x) && identical(dim(x), first_dim) && all(is.finite(x))
  }, logical(1))
  if (!all(valid)) {
    .input_error(
      "Every task result must be finite numeric data with identical dimensions."
    )
  }

  ordered <- values[order(task_id, method = "radix")]
  Reduce(`+`, ordered)
}
