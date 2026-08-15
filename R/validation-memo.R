# Structural-equality validation memo ---------------------------------------
#
# Every identity boundary in this package re-checks the object it is handed,
# and that is the right default: a compiled plan, relation fit, or covariance
# form must never be trusted because of where it came from. The cost is that a
# searchlight sweep hands the *same* immutable plan, fit, and metric to the
# same validators tens of thousands of times. On a 577-centre frame the
# analytic sampling-covariance path spent essentially all of its time
# revalidating objects it had already validated microseconds earlier.
#
# Re-admitting an object because it is `identical()` to one already validated
# in full is not a weaker check. `identical()` is total structural equality, so
# a hit is a proof that the object would pass the same validation byte for
# byte; R's copy-on-modify semantics mean a mutated object is never identical
# to the one recorded. Nothing is skipped on a miss, and a recorded depth
# satisfies a later request only if it is at least as deep.
#
# This is an execution optimization in the sense of effect-form-contract.md
# section 9: it removes administrative work without touching any semantic
# input, so no scientific plan identity changes.
#
# The store holds one entry per literal key, and the keys are a closed set
# written in this package's sources, so it cannot grow with the length of a
# sweep. Each entry retains one STRONG reference for the process lifetime:
# after the caller drops its own copy, the last-validated plan, fit, or
# metric under each key stays reachable here (including any in-memory
# relation sources it embeds) until another object replaces it or
# `.reset_validation_memo()` runs. The retained set is bounded by the key
# count times the largest single object, not by sweep length;
# `.validation_memo_bytes()` reports the current retention.

.validation_memo <- new.env(parent = emptyenv())

# TRUE when `object` is structurally identical to the object last validated
# under `key` at a depth at least as deep as the one requested.
.validated_before <- function(object, key, deep = TRUE) {
  entry <- .validation_memo[[key]]
  !is.null(entry) && (entry$deep || !isTRUE(deep)) &&
    identical(entry$object, object)
}

.record_validated <- function(object, key, deep = TRUE) {
  assign(key, list(object = object, deep = isTRUE(deep)),
    envir = .validation_memo)
  invisible(object)
}

# Drops every retained object. Validation behaviour is unchanged by this; only
# the next validation of each object does full work again. Used by tests that
# assert the memo is an optimization rather than a semantic layer, and
# available to memory-sensitive callers that have dropped large fitted
# objects and want the memo's strong references released too.
.reset_validation_memo <- function() {
  rm(list = ls(.validation_memo, all.names = TRUE), envir = .validation_memo)
  invisible(NULL)
}

# Bytes currently retained by the memo's strong references, per key.
.validation_memo_bytes <- function() {
  keys <- ls(.validation_memo, all.names = TRUE)
  stats::setNames(vapply(keys, function(key) {
    as.numeric(utils::object.size(.validation_memo[[key]]$object))
  }, numeric(1)), keys)
}
