# Argument guards -------------------------------------------------------------
#
# Layer 1 (see design/architecture.md). This file calls only
# `message-helpers.R` and `conditions.R`, so every layer above is free to use
# it.
#
# Before this file the package hand-rolled the same five shape tests at
# hundreds of sites: `!is.character(x) || length(x) != 1L || is.na(x) ||
# !nzchar(x)` appeared verbatim dozens of times, each followed by its own
# `stop()` and its own idea of how to phrase the failure. The tests were
# consistent; the messages were not, and half of them never said what had
# actually arrived.
#
# There are two kinds of helper here, because validators need both. The
# `.is_*` predicates answer; the `.check_*` guards raise. Below them are the
# two halves of the sealed-record validator prologue.

# Shape predicates -------------------------------------------------------------
#
# The `.check_*` helpers raise; these only answer. A validator that checks a
# dozen fields at once and reports them with one message cannot use a raiser
# for each field, but it can say what it means:
#
#   if (!.is_string(x$kind) || !.is_flag(x$complete) ||
#       !.strong_sha256(x$signature)) {
#     .input_error("Record fields are missing or noncanonical.")
#   }
#
# in place of eleven inlined clauses. Each predicate tests its own
# precondition first, in the order the hand-written chains did --- `is.numeric`
# before `is.finite`, `is.character` before `nzchar` --- so it is total on any
# input and can be evaluated in any position of a condition.

# One nonempty character string (NA is never a string).
.is_string <- function(x, allow_empty = FALSE) {
  is.character(x) && length(x) == 1L && !is.na(x) && (allow_empty || nzchar(x))
}

# A character vector, every entry present and nonempty.
.is_strings <- function(x, unique = FALSE) {
  is.character(x) && !anyNA(x) && all(nzchar(x)) &&
    (!unique || !anyDuplicated(x))
}

# One TRUE or FALSE.
.is_flag <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}

# One number, finite unless `finite = FALSE`.
.is_number <- function(x, finite = TRUE) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && (!finite || is.finite(x))
}

# One whole number at least `min`.
.is_count <- function(x, min = 1L) {
  .is_number(x) && x %% 1 == 0 && x >= min
}

# A numeric vector with no missing or infinite entries.
.is_finite_numeric <- function(x) {
  is.numeric(x) && all(is.finite(x))
}

# A numeric matrix with no missing or infinite entries.
.is_finite_matrix <- function(x) {
  is.matrix(x) && is.numeric(x) && all(is.finite(x))
}

# Raising guards ---------------------------------------------------------------
#
# Each takes `arg`, the name to show in backticks, and returns the value
# (coerced where the coercion is part of the contract, as with counts), so a
# guard may be written as an assignment:
#
#   id <- .check_string(id, "id")
#   workers <- .check_count(workers, "workers")
#
# `what` overrides the expectation clause. It exists so a call site can keep
# wording that was chosen deliberately --- "one nonempty identifier" rather
# than the generic "one nonempty character string" --- while still gaining
# the observed value and the condition class.

# "`arg` must be <expectation>; received <value>."
.check_failed <- function(arg, what, x, expected = what) {
  received <- .msg_received(x)
  .input_error(
    sprintf("`%s` must be %s; received %s.", arg, what, received),
    arg = arg, received = received, expected = expected
  )
}

# One character string. `allow_empty = TRUE` admits "" but never NA.
.check_string <- function(x, arg, allow_empty = FALSE, what = NULL) {
  if (is.null(what)) {
    what <- if (allow_empty) {
      "one character string"
    } else {
      "one nonempty character string"
    }
  }
  if (!.is_string(x, allow_empty)) .check_failed(arg, what, x)
  x
}

# One nonempty identifier, unnamed. Identifiers reach the package as elements
# of named vectors often enough that the `unname()` is part of the contract:
# an id that carried its own name would put that name into a signature.
# It lived in `R/effect-map.R` because a map was the first record with ids in
# it, which made `study-facts.R`, `study.R`, `design-model.R` and
# `observation-model.R` all appear to depend on effect maps to spell a string.
.validate_nonempty_id <- function(value, name) {
  .check_string(value, name, what = "one nonempty identifier")
  unname(value)
}

# One TRUE or FALSE. NA is not a flag: a guard that admitted it would push the
# three-valued logic downstream into an `if` that errors far from the cause.
.check_flag <- function(x, arg, what = "TRUE or FALSE") {
  if (!.is_flag(x)) .check_failed(arg, what, x)
  x
}

# One whole number at least `min`. Returns an integer, because a count that
# arrived as a double and is used as an index or a dimension must not carry a
# fractional part into arithmetic downstream.
.check_count <- function(x, arg, min = 1L, max = NULL, what = NULL) {
  if (is.null(what)) {
    what <- if (min == 1L) {
      "one positive whole number"
    } else {
      sprintf("one whole number of at least %s", format(min))
    }
    if (!is.null(max)) {
      what <- paste0(what, sprintf(" and at most %s", format(max)))
    }
  }
  if (!.is_count(x, min) || (!is.null(max) && x > max)) {
    .check_failed(arg, what, x)
  }
  as.integer(x)
}

# One number. `finite` rejects NA, NaN, and infinities; `positive` requires
# strictly greater than zero, `nonnegative` requires at least zero.
.check_number <- function(x, arg, finite = TRUE, positive = FALSE,
                          nonnegative = FALSE, what = NULL) {
  if (is.null(what)) {
    what <- paste0(
      "one ",
      if (positive) "positive " else if (nonnegative) "nonnegative " else "",
      if (finite) "finite " else "",
      "number"
    )
  }
  if (!.is_number(x, finite) || (positive && x <= 0) ||
      (nonnegative && x < 0)) {
    .check_failed(arg, what, x)
  }
  x
}

# One matrix, optionally of a declared shape. The dimension clause is part of
# the expectation rather than a second guard, so a caller who passed a
# transposed matrix reads both numbers in one sentence.
.check_matrix <- function(x, arg, finite = TRUE, nrow = NULL, ncol = NULL,
                          numeric = TRUE, what = NULL) {
  shape <- if (is.null(nrow) && is.null(ncol)) {
    ""
  } else {
    sprintf(" with %s and %s",
      if (is.null(nrow)) "any number of rows" else .msg_count(nrow, "row"),
      if (is.null(ncol)) {
        "any number of columns"
      } else {
        .msg_count(ncol, "column")
      })
  }
  # Finiteness is only meaningful, and only checked, for a numeric matrix.
  finite <- finite && numeric
  if (is.null(what)) {
    what <- paste0("a ", if (numeric) "numeric " else "",
      if (finite) "matrix with no missing or infinite entries" else "matrix",
      shape)
  }
  ok <- if (finite) {
    .is_finite_matrix(x)
  } else {
    is.matrix(x) && (!numeric || is.numeric(x))
  }
  if (!ok || (!is.null(nrow) && nrow(x) != nrow) ||
      (!is.null(ncol) && ncol(x) != ncol)) {
    .check_failed(arg, what, x)
  }
  x
}

# One object of a declared class. `from` names the constructor that builds
# one, because "must be an effect_geometry_plan" is only actionable to a
# reader who already knows which call returns one.
.check_class <- function(x, class, arg, from = NULL, what = NULL) {
  if (!inherits(x, class)) {
    if (is.null(what)) {
      what <- sprintf("an object of class `%s`", class[[1L]])
    }
    if (!is.null(from)) {
      what <- sprintf("%s from `%s`", what, from)
    }
    .check_failed(arg, what, x)
  }
  x
}

# Sealed records ---------------------------------------------------------------
#
# Every value the package hands back is sealed: a classed list whose field
# names are fixed at construction and whose last field is a content hash over
# the rest. Its validator therefore opens the same way every time --- right
# class, is a list, exactly these names in exactly this order --- and closes
# by rebuilding the hash and comparing.
#
# `.sealed_fields()` is the opening test as a predicate, so a validator can
# write it as one clause of the compound condition it already has rather than
# three. `.check_signature()` is the closing test, and it raises
# `effect_contract_error`: the fields are all individually well formed, and
# what failed is the agreement between the record and its own identity.

.sealed_fields <- function(x, class, expected) {
  inherits(x, class) && is.list(x) && identical(names(x), expected)
}

.check_signature <- function(actual, expected, message, arg = NULL) {
  if (!identical(actual, expected)) {
    .contract_error(message, arg = arg,
      received = .msg_signature(actual), expected = .msg_signature(expected))
  }
  invisible(TRUE)
}
