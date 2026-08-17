# The condition taxonomy -----------------------------------------------------
#
# Every failure this package raises answers one of four questions, and the
# class says which, so a caller can branch on the cause instead of matching
# message prose:
#
#   effect_input_error       an argument has the wrong type, shape, or value.
#                            The caller fixes it by passing something else.
#   effect_contract_error    two objects, or an object and its own recorded
#                            identity, disagree. Nothing about either object is
#                            malformed; they simply do not belong together.
#   effect_invariant_error   the package computed something impossible. This is
#                            a crossform bug, and the message says so.
#   effect_capability_refusal  the requested interpretation cannot be earned
#                            from the supplied objects. Not an error in the
#                            input: an honest refusal, with reasons and
#                            remedies attached (see `catch_refusal()`).
#
# The first three share the parent class `effect_error` and carry the same
# three optional fields --- `arg`, `received`, `expected` --- so one print
# method serves all of them. Refusals keep their own richer shape.
#
# `call. = FALSE` is accepted, and defaults to FALSE, purely so that these
# constructors are drop-in replacements for the `stop(message, call. = FALSE)`
# they were converted from. A condition with a NULL `call` prints its message
# without an "Error in f():" prefix, which is what every message in this
# package was written to assume.

.effect_condition <- function(class, message, arg = NULL, received = NULL,
                              expected = NULL, call. = FALSE) {
  if (!is.character(message) || length(message) != 1L || is.na(message)) {
    message <- paste(as.character(message), collapse = "")
  }
  structure(list(
    message = message,
    call = if (isTRUE(call.) && sys.nframe() > 2L) sys.call(-2L) else NULL,
    arg = arg,
    received = received,
    expected = expected
  ), class = c(class, "effect_error", "error", "condition"))
}

# An argument is wrong: bad type, bad shape, bad value. User fixable.
.input_error <- function(message, arg = NULL, received = NULL,
                         expected = NULL, call. = FALSE) {
  stop(.effect_condition("effect_input_error", message, arg, received,
    expected, call.))
}

# Two objects disagree: an identity, receipt, or signature does not match the
# thing it is supposed to describe. Both objects may be individually valid.
.contract_error <- function(message, arg = NULL, received = NULL,
                            expected = NULL, call. = FALSE) {
  stop(.effect_condition("effect_contract_error", message, arg, received,
    expected, call.))
}

# The package violated its own invariant. Nothing the caller passed can
# explain this, so the message asks for a report rather than a fix.
.invariant_error <- function(message, arg = NULL, received = NULL,
                             expected = NULL, call. = FALSE) {
  stop(.effect_condition("effect_invariant_error",
    paste0(message, " This is a crossform bug; please report it at ",
      "<https://github.com/bbuchsbaum/crossform/issues>."),
    arg, received, expected, call.))
}

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
#' @seealso [crossform_conditions] for the other condition classes crossform
#'   raises, and how to branch on them.
#' @family conditions
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

#' Conditions raised by crossform
#'
#' Every failure crossform raises is a classed condition, so a caller can
#' branch on the *cause* rather than matching message prose. There are four
#' classes, and each answers a different question about what went wrong.
#'
#' @section The classes:
#'
#' \describe{
#'   \item{`effect_input_error`}{An argument has the wrong type, shape, or
#'     value --- a character vector where one string was expected, a matrix
#'     holding `NA`, a negative count. The caller fixes it by passing
#'     something else. May carry `$arg` (the argument name), `$received` (a
#'     short description of what arrived), and `$expected` (what was wanted);
#'     see "The three optional fields" below.}
#'   \item{`effect_contract_error`}{Two objects disagree, or an object
#'     disagrees with its own recorded identity: a result whose receipt names
#'     a different scientific plan, a frame compiled against a different
#'     neural domain, a signature that no longer matches the fields it
#'     summarizes. Each object may be individually well formed; they simply do
#'     not belong together. The fix is to pass objects that were built from
#'     the same declarations, not to repair either one.}
#'   \item{`effect_invariant_error`}{The package computed something
#'     impossible. Nothing the caller passed can explain it, so the message
#'     asks for a bug report rather than a change of input.}
#'   \item{`effect_capability_refusal`}{Not an error in the input: the
#'     requested interpretation cannot be earned from the objects supplied ---
#'     for example, analytic standard errors from a fit that did not retain
#'     the residual channel they require. Refusals carry `$capability`,
#'     `$namespace`, all unmet `$reasons`, and concrete `$remedies`. See
#'     [catch_refusal()].}
#' }
#'
#' The first three inherit from `effect_error`, so
#' `tryCatch(expr, effect_error = ...)` catches any of them while letting a
#' refusal through. All four inherit from `error`, so ordinary `try()` and
#' `tryCatch(expr, error = ...)` behave exactly as before.
#'
#' @section The three optional fields:
#'
#' `$arg`, `$received`, and `$expected` are a convenience, not a contract.
#' They are filled in by the shared argument guards (the `.check_*` family
#' behind most exported entry points) and by the hand-written checks at the
#' entry points that know the three values, which covers the argument errors
#' a caller provokes in ordinary use. Every other site --- internal checks
#' whose failure is about a whole record rather than one argument, and checks
#' whose "expected" is a paragraph rather than a phrase --- leaves them
#' `NULL` and says everything in `conditionMessage()`.
#'
#' So branch on the **class** first, which is always there, and treat the
#' fields as something to report when present:
#'
#' ```r
#' if (!is.null(condition$arg)) {
#'   message("bad argument `", condition$arg, "`: ", condition$received)
#' }
#' ```
#'
#' The same three fields exist on `effect_contract_error` and
#' `effect_invariant_error` under the same rule; a contract error that
#' compares two identities typically fills `$received` and `$expected` with
#' the two shortened signatures.
#'
#' @return These are condition classes, not functions; there is nothing to
#'   call. [catch_refusal()] captures a refusal as a value.
#' @family conditions
#' @seealso [catch_refusal()]
#' @examples
#' domain <- abstract_domain(2, id = "conditions-example")
#' relation <- relation(
#'   list(a = matrix(1:4, 2), b = matrix(2:5, 2)),
#'   effects = c("x", "y"), domain = domain
#' )
#' plan <- plan_geometry(
#'   relation, compile_frame(whole_brain(), domain),
#'   cross_partitions(relation, independence = "independent")
#' )
#'
#' # Branch on the class. Three weights for two effects is an input error.
#' tryCatch(
#'   contrast_energy(plan, c(1, -1, 0)),
#'   effect_input_error = function(e) conditionMessage(e),
#'   effect_contract_error = function(e) "these objects do not belong together"
#' )
#'
#' # A guard-raised error additionally names the argument and the value.
#' guarded <- tryCatch(
#'   effect_space(c("x", "y"), basis_id = 42),
#'   effect_input_error = function(e) e
#' )
#' guarded$arg
#' guarded$received
#' guarded$expected
#'
#' # Not every site fills them, so test before you use them. Here the whole
#' # explanation is in the message.
#' bare <- tryCatch(abstract_domain(2, id = ""), effect_input_error = identity)
#' is.null(bare$arg)
#' conditionMessage(bare)
#'
#' # `effect_error` catches input, contract, and invariant failures alike,
#' # while letting a capability refusal through.
#' tryCatch(abstract_domain(-1), effect_error = function(e) "caught")
#' @name crossform_conditions
NULL
