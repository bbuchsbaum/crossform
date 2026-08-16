# Classed capability refusals ------------------------------------------------
#
# A capability refusal is the package's contract-level "no": the requested
# interpretation cannot be earned from the supplied objects. Refusals must be
# machine-distinguishable from shape errors and internal-invariant failures,
# so they are raised as classed conditions carrying the missing capability,
# the namespace it belongs to, every unmet reason (not only the first), and
# concrete remedies.

.capability_refusal <- function(message, capability, namespace,
                                reasons = character(),
                                remedies = character()) {
  stopifnot(
    is.character(message), length(message) == 1L, nzchar(message),
    is.character(capability), length(capability) == 1L, nzchar(capability),
    is.character(namespace), length(namespace) == 1L, nzchar(namespace),
    is.character(reasons), is.character(remedies)
  )
  condition <- structure(list(
    message = message,
    call = NULL,
    capability = capability,
    namespace = namespace,
    reasons = reasons,
    remedies = remedies
  ), class = c("effect_capability_refusal", "error", "condition"))
  stop(condition)
}

#' Inspect a capability refusal
#'
#' Every contract-level refusal in crossform signals a condition of class
#' `effect_capability_refusal`. The condition carries the missing
#' `$capability`, its `$namespace`, all unmet `$reasons`, and suggested
#' `$remedies`, so callers can branch on the *cause* rather than matching
#' message prose.
#'
#' @param expr An expression expected to refuse.
#' @return The captured `effect_capability_refusal` condition, or `NULL` if
#'   `expr` succeeded.
#' @examples
#' domain <- abstract_domain(2, id = "refusal-example")
#' relation <- relation(
#'   list(a = matrix(1:4, 2), b = matrix(2:5, 2)),
#'   effects = c("x", "y"), domain = domain
#' )
#' refusal <- catch_refusal(
#'   rdm_sampling_covariance(
#'     plan_geometry(relation, compile_frame(whole_brain(), domain),
#'       cross_partitions(relation, independence = "independent")),
#'     relation, target = "null", at = 1L
#'   )
#' )
#' refusal$capability
#' refusal$reasons
#' @export
catch_refusal <- function(expr) {
  tryCatch({
    force(expr)
    NULL
  }, effect_capability_refusal = function(condition) condition)
}
