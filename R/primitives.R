# Leaf primitives ------------------------------------------------------------
#
# This file is the bottom of the package. Nothing here may call another
# crossform internal outside the primitive layer, so every other file is free
# to depend on it without creating a cycle.
#
# Layering (a file may only call downward; see design/architecture.md):
#
#   1. primitives   R/primitives.R, R/message-helpers.R, R/conditions.R,
#                   R/check.R
#   2. values       R/domain.R, R/frame.R, R/pairing.R, R/relation.R,
#                   R/effect-space.R, R/effect-map.R, R/metric.R, R/source.R,
#                   R/study.R, R/receipt.R, R/measurement.R, R/bridge.R,
#                   R/capabilities.R, R/compute-policy.R, R/memory-plan.R, ...
#   3. plans        R/geometry-plan.R, R/relation-plan.R, R/crossnobis.R,
#                   R/evidence-*.R, R/coupling-plan.R, ...
#   4. compiler     R/compiler.R, R/execution-driver.R, R/kernel.R, R/task.R,
#      /execution   R/storage.R, R/measurement-kernel.R
#   5. results      R/result.R, R/views.R, R/geometry-entry.R,
#      /views       R/coupling-views.R, R/tomography.R,
#                   R/measurement-result.R, R/format-results.R,
#                   R/print-methods.R, R/plot-methods.R
#   6. adapters     R/adapter-*.R, R/neuroim2-adapter.R, R/benchmark.R,
#      /facade      R/example-data.R, R/evidence-api.R
#
# The rule is checked by tests/testthat/test-architecture.R, which carries the
# full layer map and an explicit register of the upward edges that remain,
# each with its follow-up. The package has no `Collate:` field on purpose:
# every definition resolves at call time, so file order carries no meaning and
# a second statement of dependency would only drift.

# Content addressing ---------------------------------------------------------
#
# Every identity in the package is a SHA-256 over an R value. `serialize` and
# `serializeVersion` are pinned here, in one place, because a change to either
# would silently invalidate every recorded signature in the package and in
# users' saved objects.

.sha256_signature <- function(x, prefix = "sha256:") {
  paste0(prefix, digest::digest(
    x, algo = "sha256", serialize = TRUE, serializeVersion = 2L
  ))
}

# For values that are already a character representation, where serializing the
# R object would fold R's internal encoding into the identity.
.sha256_string <- function(x, prefix = "sha256:") {
  paste0(prefix, digest::digest(x, algo = "sha256", serialize = FALSE))
}

.sha256_file <- function(path) {
  paste0("sha256:", digest::digest(file = path, algo = "sha256"))
}

# Hot path: this runs on every field of every validated plan, so it avoids the
# regular-expression engine for the two cheap structural checks first. A
# TRE-compiled `[[:xdigit:]]{64}` costs about 55us per call; the length and
# prefix tests cost about 0.6us and reject almost everything. The `nchar()`
# guard also pins the semantics of the regex it replaces, whose `$` would
# otherwise admit a trailing newline under `perl = TRUE`.
#
# Equivalent to `grepl("^sha256:[[:xdigit:]]{64}$", x)` on every valid input,
# including NA, zero-length, uppercase hex, and embedded or trailing newlines.
# One boundary differs: a string carrying invalid multibyte bytes raises an
# encoding error here, where the regex returned FALSE with a warning. Both
# fail closed, no value the package produces can reach it, and the sole caller
# is `.strong_sha256()`.
.is_sha256_signature <- function(x) {
  nchar(x) == 71L & startsWith(x, "sha256:") &
    grepl("^[[:xdigit:]]{64}$", substr(x, 8L, 71L), perl = TRUE)
}

# One scalar signature, not a vector of them.
.strong_sha256 <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    .is_sha256_signature(value)
}

# Tiling ---------------------------------------------------------------------

.tile_starts <- function(n, size) seq.int(1L, n, by = size)

.validate_tile_size <- function(x, name) {
  .check_count(x, name, what = "one positive integer")
  as.integer(x)
}

# Symmetric packing ----------------------------------------------------------
#
# svec packs the lower triangle of a symmetric matrix with off-diagonal
# entries scaled by sqrt(2), so the Euclidean inner product of two packed
# vectors equals the Frobenius inner product of the matrices. unsvec inverts
# it exactly.

.svec_symmetric <- function(x) {
  q <- nrow(x)
  out <- numeric(q * (q + 1L) / 2L)
  k <- 0L
  for (column in seq_len(q)) {
    for (row in column:q) {
      k <- k + 1L
      out[[k]] <- x[row, column] * if (row == column) 1 else sqrt(2)
    }
  }
  out
}

# One query column back to the operator it names. Kernels that must stay
# memory bounded decode a single column with `.physical_query_operator()`;
# `.physical_query_operators()` is the eager list form, for callers that
# already hold every operator at once.
.physical_query_operator <- function(column, q_left, q_right, codec) {
  if (identical(codec, "rectangular")) {
    matrix(column, q_left, q_right)
  } else {
    .unsvec_symmetric(column, q_left)
  }
}

.physical_query_operators <- function(query, q_left, q_right, codec) {
  lapply(seq_len(ncol(query)), function(view) {
    .physical_query_operator(query[, view], q_left, q_right, codec)
  })
}

.unsvec_symmetric <- function(value, q) {
  if (!.is_finite_numeric(value) || length(value) != q * (q + 1L)/2L) {
    .input_error("Packed geometry has the wrong width or non-finite values.")
  }
  out <- matrix(0, q, q)
  coordinate <- 0L
  for (column in seq_len(q)) {
    for (row in column:q) {
      coordinate <- coordinate + 1L
      entry <- value[[coordinate]] / if (row == column) 1 else sqrt(2)
      out[row, column] <- entry
      out[column, row] <- entry
    }
  }
  out
}

# Contrast alignment ---------------------------------------------------------

.align_contrast <- function(value, effects, label = "weights") {
  if (!is.numeric(value) || is.matrix(value)) {
    .input_error(sprintf(paste0(
      "`%s` must be a numeric contrast vector with one weight per effect; ",
      "received %s. The relation declares %s: %s."
    ), label, .msg_value(value), .msg_count(length(effects), "effect"),
      .msg_names(effects)))
  }
  if (length(value) != length(effects)) {
    .input_error(sprintf(paste0(
      "`%s` has %s but the relation declares %s (%s). Supply one weight per ",
      "effect, or name the weights to have them aligned for you."
    ), label, .msg_count(length(value), "value"),
      .msg_count(length(effects), "effect"), .msg_names(effects)))
  }
  if (any(!is.finite(value))) {
    .input_error(sprintf(paste0(
      "`%s` must be finite, but the weight%s for %s %s NA, NaN, or Inf. A ",
      "contrast weight of zero excludes an effect."
    ), label, if (sum(!is.finite(value)) == 1L) "" else "s",
      .msg_positions(!is.finite(value),
        if (is.null(names(value))) effects else names(value)),
      if (sum(!is.finite(value)) == 1L) "is" else "are"))
  }
  if (!is.null(names(value))) {
    supplied <- names(value)
    if (anyNA(supplied) || any(!nzchar(supplied))) {
      .input_error(sprintf(paste0(
        "`%s` is partially named: every weight must be named, or none. The ",
        "relation declares %s."
      ), label, .msg_names(effects)))
    }
    if (anyDuplicated(supplied)) {
      .input_error(sprintf(
        "`%s` names %s more than once; each effect takes exactly one weight.",
        label, .msg_names(unique(supplied[duplicated(supplied)]))))
    }
    if (!setequal(supplied, effects)) {
      unknown <- setdiff(supplied, effects)
      absent <- setdiff(effects, supplied)
      detail <- c(
        if (length(unknown)) sprintf("%s is not a declared effect",
          .msg_names(unknown)),
        if (length(absent)) sprintf("%s has no weight", .msg_names(absent))
      )
      .input_error(sprintf(
        "`%s` does not match the relation's effects (%s): %s.",
        label, .msg_names(effects), paste(detail, collapse = "; ")))
    }
    value <- value[effects]
  }
  stats::setNames(as.numeric(value), effects)
}
