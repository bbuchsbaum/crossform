# Compact printing for user-facing objects -----------------------------------
#
# Every class a user can hold gets a one-screen summary: a `<class>` header,
# a short block of aligned `key: value` lines, and -- where the object is lazy
# or has an obvious successor -- a closing `next:` or `state:` hint naming the
# exact call to make. Printing never materializes, never dereferences a source,
# and never emits closures, environments, whole matrices, or whole hashes.

# One descriptor, not two methods --------------------------------------------
#
# A sealed value record used to carry a `print.<class>` whose body was
# `.pf_emit("<class>", list(...))` and a `format.<class>` whose body was
# `.pf_inline("<class>", ...)`. The class name was written six times --- twice
# in a method name, twice in an emitted header, twice in NAMESPACE --- and 163
# such methods differed only in their field list. `.pf_records` below keeps
# only that difference; `.onLoad()` registers one shared printer and one
# shared formatter per entry. The header is now emitted by the engine from the
# name the method was registered under, so `<class>` and the class it
# dispatches on can no longer drift apart.
#
# The alternative was a shared parent class appended to every record's class
# vector, so that one `print.crossform_record` caught them all by inheritance.
# It was rejected on two counts. The class vector is set at some two hundred
# constructor sites across forty files, none of them a printing file. And
# every identity in this package is `.sha256_signature()` over a *serialized*
# R value, where the class attribute is part of what `serialize()` writes:
# each of the 143 signature sites would have had to be audited for whether it
# digests a record or only that record's signature field before a tail class
# could be called safe, and a wrong answer at any one of them silently moves a
# published estimand id. Registering methods changes no class vector, no
# dispatch, and no printed byte.
#
# Records whose printer is not a header over aligned key/value lines --- the
# `cat()`-based ones, the result previews, the `.format_counted_result()`
# formats --- remain ordinary methods, defined below and named in NAMESPACE.

.pf_max_lines <- 12L
.pf_line_width <- 76L

# Emit a `<class>` header followed by aligned key/value lines. `fields` is a
# named list; NULL entries are dropped so callers can omit lines conditionally.
.pf_emit <- function(class, fields) {
  keep <- !vapply(fields, is.null, logical(1))
  fields <- fields[keep]
  cat("<", class, ">\n", sep = "")
  if (!length(fields)) {
    return(invisible(NULL))
  }
  keys <- names(fields)
  values <- vapply(fields, function(value) {
    paste0(as.character(value), collapse = "")
  }, character(1))
  width <- max(nchar(keys)) + 2L
  lines <- sprintf("  %-*s%s", width, paste0(keys, ":"), values)
  long <- nchar(lines) > .pf_line_width
  lines[long] <- paste0(substr(lines[long], 1L, .pf_line_width - 3L), "...")
  cat(lines, sep = "\n")
  invisible(NULL)
}

# One-line `<class: detail>` summary used by the format() methods.
.pf_inline <- function(class, ...) {
  detail <- paste0(c(...), collapse = ", ")
  if (!nzchar(detail)) {
    return(sprintf("<%s>", class))
  }
  sprintf("<%s: %s>", class, detail)
}

# Show a content hash as its first 12 hex digits. Never print a whole digest.
# The same shortening errors use, so a printed identity and a quoted one in a
# refusal are comparable by eye.
.pf_sig <- function(x, chars = 12L) .msg_signature(x, chars)

# Show a set as `a, b, c (+k more)`.
.pf_set <- function(x, max = 3L, empty = "none") {
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (!length(x)) {
    return(empty)
  }
  if (length(x) <= max) {
    return(paste(x, collapse = ", "))
  }
  paste0(paste(x[seq_len(max)], collapse = ", "), " (+", length(x) - max,
    " more)")
}

# Four significant digits, never a 15-digit float.
.pf_num <- function(x, digits = 4L, empty = "none") {
  if (is.null(x) || !length(x)) {
    return(empty)
  }
  value <- suppressWarnings(as.numeric(x))
  # A printer must never warn. If the input was not numeric at all, show it
  # as text rather than emitting a coercion warning and a row of NA.
  if (anyNA(value) && !anyNA(x)) {
    return(.pf_set(as.character(x), max = 4L, empty = empty))
  }
  paste(format(signif(value, digits), trim = TRUE), collapse = ", ")
}

# A regularization record carries `kind` and `applied` alongside its lambdas,
# so the fields are read by name; unlisting it would coerce them all to text.
.pf_regularization <- function(regularization) {
  if (is.null(regularization)) {
    return("none")
  }
  kind <- .pf_or(regularization$kind, "none")
  lambdas <- c(regularization$lambda_left, regularization$lambda_right)
  if (!length(lambdas)) {
    return(kind)
  }
  if (length(unique(lambdas)) == 1L) {
    return(paste0(kind, " ", .pf_num(lambdas[[1L]])))
  }
  paste0(kind, " (left ", .pf_num(regularization$lambda_left), ", right ",
    .pf_num(regularization$lambda_right), ")")
}

# Dimensions as `n x m`, for objects we refuse to print in full.
.pf_dim <- function(x) {
  dimensions <- if (is.null(dim(x))) length(x) else dim(x)
  paste(dimensions, collapse = " x ")
}

# Same rendering for a field that already holds an extent vector.
.pf_shape <- function(x) {
  if (is.null(x) || !length(x)) {
    return("unknown")
  }
  paste(as.integer(x), collapse = " x ")
}

.pf_yn <- function(x) if (isTRUE(x)) "yes" else "no"

.pf_or <- function(x, empty = "none") {
  if (is.null(x) || !length(x)) empty else as.character(x)[[1L]]
}

# A wrapped caveat under the key/value block. It lives here with the other
# printing primitives rather than beside its first caller: `format-results.R`
# and this file are both layer 5, so either may call the other, but only one
# direction at a time --- a helper defined there and used here closes the
# cycle `benchmarks/call-graph-scc.R` ratchets at one component per file.
.pf_note <- function(text) {
  cat(strwrap(text, width = .pf_line_width, prefix = "    ", initial = "  "),
    sep = "\n")
}

# Number of stored (structurally nonzero) entries in a sparse or dense matrix.
.pf_stored <- function(x) {
  if (methods::is(x, "sparseMatrix")) {
    slots <- methods::slotNames(x)
    if ("x" %in% slots) {
      return(length(methods::slot(x, "x")))
    }
    if ("i" %in% slots) {
      return(length(methods::slot(x, "i")))
    }
  }
  length(x)
}

# Byte counts in human units. Deterministic given the count, so safe to print.
.pf_bytes <- function(n) {
  if (is.null(n) || !length(n) || !is.finite(n[[1L]])) {
    return("unknown")
  }
  format(structure(as.numeric(n[[1L]]), class = "object_size"), units = "auto")
}

# One compiled pipeline stage as `<operation> over <axis>`.
.pf_stage <- function(stage) {
  operation <- stage$operation
  kind <- if (is.list(operation)) operation$kind else operation
  paste0(.pf_or(kind, "none"), " over ", .pf_or(stage$axis, "nothing"))
}

# Wrap a character vector into `    - text` bullets that respect the line
# budget, continuing over-long entries with a hanging indent.
.pf_bullets <- function(x) {
  unlist(lapply(as.character(x), function(entry) {
    strwrap(entry, width = .pf_line_width, prefix = "      ",
      initial = "    - ")
  }), use.names = FALSE)
}

# Capability records ---------------------------------------------------------
#
# Capability objects come in two shapes: a flat named record of logical or
# character flags, and an adjudicated record carrying `$available`, a named
# `$capabilities` record, and a `$reasons` data frame. One formatter covers
# both so every `*_capabilities` class reads the same way.

.pf_capability_flags <- function(x, drop = character()) {
  drop <- c(drop, "signature", "schema_version", "plan")
  flags <- x[!names(x) %in% drop]
  flags <- flags[vapply(flags, function(value) {
    length(value) == 1L && (is.logical(value) || is.character(value))
  }, logical(1))]
  vapply(flags, function(value) {
    if (is.logical(value)) .pf_yn(value) else .pf_sig(as.character(value))
  }, character(1))
}

# The single shared capability renderer. `x` may be a flat flag record or an
# adjudicated `$available`/`$capabilities`/`$reasons` record.
.pf_capabilities_emit <- function(class, x, max = .pf_max_lines) {
  adjudicated <- is.list(x) && !is.null(x$capabilities) && !is.null(x$reasons)
  reasons <- if (adjudicated) x$reasons else NULL
  flags <- if (adjudicated) {
    .pf_capability_flags(x$capabilities)
  } else {
    .pf_capability_flags(x)
  }
  first_reason <- NULL
  if (is.data.frame(reasons) && nrow(reasons) > 0L) {
    column <- if ("why" %in% names(reasons)) "why" else names(reasons)[[1L]]
    first_reason <- paste0(reasons[[1L]][[1L]], " - ",
      as.character(reasons[[column]][[1L]]))
  } else if (is.character(reasons) && length(reasons)) {
    first_reason <- reasons[[1L]]
  }
  budget <- max - (!is.null(first_reason)) - adjudicated
  shown <- flags
  hidden <- 0L
  if (length(shown) > budget) {
    keep <- max(0L, budget - 1L)
    hidden <- length(shown) - keep
    shown <- shown[seq_len(keep)]
  }
  fields <- list()
  if (adjudicated) {
    fields[["available"]] <- .pf_yn(x$available)
  }
  fields <- c(fields, as.list(shown))
  if (hidden > 0L) {
    fields[["..."]] <- paste0("+", hidden, " more capabilities")
  }
  if (!is.null(first_reason)) {
    fields[["first reason"]] <- first_reason
  }
  .pf_emit(class, fields)
}

.pf_capabilities_inline <- function(class, x) {
  adjudicated <- is.list(x) && !is.null(x$capabilities) && !is.null(x$reasons)
  flags <- if (adjudicated) {
    .pf_capability_flags(x$capabilities)
  } else {
    .pf_capability_flags(x)
  }
  logicals <- flags[flags %in% c("yes", "no")]
  granted <- sum(logicals == "yes")
  detail <- if (length(logicals)) {
    paste0(granted, "/", length(logicals), " granted")
  } else {
    paste0(length(flags), " properties")
  }
  if (adjudicated) {
    detail <- c(if (isTRUE(x$available)) "available" else "unavailable", detail)
  }
  .pf_inline(class, detail)
}

# Capability refusals --------------------------------------------------------

#' @export
format.effect_capability_refusal <- function(x, ...) {
  header <- c(
    "<effect_capability_refusal>",
    sprintf("  %-13s%s", "capability:", .pf_or(x$capability, "unknown")),
    sprintf("  %-13s%s", "namespace:", .pf_or(x$namespace, "unknown"))
  )
  reasons <- if (length(x$reasons)) {
    c("  reasons:", .pf_bullets(x$reasons))
  } else {
    "  reasons:     none recorded"
  }
  remedies <- if (length(x$remedies)) {
    c("  remedies:", .pf_bullets(x$remedies))
  } else {
    "  remedies:    none recorded"
  }
  c(header, reasons, remedies,
    "  state:       refused; no partial result was produced")
}

#' @export
print.effect_capability_refusal <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}

# Input, contract, and invariant failures -------------------------------------
#
# The three classes under `effect_error` share one shape --- a message plus an
# optional argument name, observed value, and expectation --- so they share
# one printer, in the refusal's layout. A caught condition prints as a record
# of what was asked and what was wrong with it, not as a wall of prose.

#' @export
format.effect_error <- function(x, ...) {
  lines <- c(
    sprintf("<%s>", class(x)[[1L]]),
    strwrap(conditionMessage(x), width = .pf_line_width, prefix = "  ")
  )
  field <- function(label, value) {
    if (is.null(value) || !length(value)) return(character())
    sprintf("  %-11s%s", paste0(label, ":"), as.character(value)[[1L]])
  }
  c(lines,
    field("arg", x$arg),
    field("received", x$received),
    field("expected", x$expected))
}

#' @export
print.effect_error <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}

# "point (alpha 0.25), narrow (alpha 0.75, scale 1.01)" -- the one line that
# says what a stacked family is made of and with what weight.
.pf_frame_family_members <- function(specification) {
  paste0(vapply(specification$members, function(member) {
    paste0(member$family, " (alpha ", .pf_num(member$alpha),
      if (is.na(member$scale)) "" else paste0(", scale ", .pf_num(member$scale)),
      ")")
  }, character(1)), collapse = ", ")
}

# The RDM block keeps the `cat()` rendering it has always had rather than
# moving to `.pf_emit()`, because `.pf_emit()` clips at 76 columns and the
# residual line -- degrees of freedom, effective dimension, and which
# tr(Sigma^2) estimator produced them -- runs past that. Truncating the one
# line that says how trustworthy the noise term is would be the wrong saving,
# so this block is wider than the compact budget on purpose. Everything else
# is the same aligned key/value form `.pf_emit()` produces.
.print_rdm_sampling_covariance <- function(x) {
  cat("<effect_sampling_covariance>\n", sep = "")
  cat("  basis:        rdm\n", sep = "")
  cat("  distances:    ", x$dimension, "\n", sep = "")
  cat("  measurement:  ", x$source$node, "\n", sep = "")
  cat("  partitions:   ", x$partitions, " (dependent pair products)\n",
    sep = "")
  cat("  target:       ", x$plan$target$target, " / ",
    x$plan$target$policy, "\n", sep = "")
  cat("  metric:       fixed\n", sep = "")
  cat("  residual:     ", .pf_residual_noise(x$source), "\n", sep = "")
  cat("  storage:      exact factorized covariance\n", sep = "")
  cat("  spatial law:  local marginal only\n", sep = "")
  if (!is.null(x$source$execution)) {
    cat("  execution:    ", x$source$execution$route, " / ",
      x$source$execution$residual_strategy, "\n", sep = "")
  }
  invisible(x)
}

# The one printer that chooses between two whole layouts rather than between
# two field lists, so it cannot be a descriptor: a descriptor's `fields`
# returns the list `.pf_emit()` renders, and there is no field list that means
# "print this other block instead". The RDM basis keeps the `cat()` form above;
# every other basis prints the aligned block here. `format()` does branch on
# the basis and *is* a descriptor, because all three branches are inline
# details.

#' @export
print.effect_sampling_covariance <- function(x, ...) {
  .validate_sampling_covariance(x, deep = FALSE)
  if (identical(x$basis, "rdm")) {
    return(.print_rdm_sampling_covariance(x))
  }
  fields <- list(
    basis = x$basis,
    coordinates = x$dimension,
    labels = .pf_set(x$labels, max = 3L),
    partitions = x$partitions,
    target = paste0(.pf_or(x$plan$target$target), " / ",
      .pf_or(x$plan$target$policy)),
    factors = paste0("signal ", .pf_dim(x$signal_factor), ", xi ",
      .pf_dim(x$xi_factor)),
    "residual noise" = .pf_residual_noise(x$source),
    storage = "exact factorized covariance",
    signature = .pf_sig(x$signature)
  )
  if (identical(x$basis, "query_bank")) {
    # A bank of contrast energies is read at one measurement, and the block is
    # a marginal. Both facts are printed, because a K-by-K matrix invites the
    # reader to treat it as a joint covariance over everything in sight.
    fields <- append(fields, list(
      queries = paste0(x$dimension, " contrast energies, lowered onto ",
        length(x$source$distance_labels), " distances"),
      measurement = .pf_or(as.character(x$source$node)),
      "spatial law" = "local marginal only; no cross-measurement covariance"
    ), after = 3L)
  }
  .pf_emit("effect_sampling_covariance", fields)
  invisible(x)
}

#' @export
print.effect_sampling_covariance_batch <- function(x, ...) {
  if (!length(x) || !inherits(x[[1L]], "effect_sampling_covariance")) {
    .input_error("`x` must come from `rdm_sampling_covariance()`.")
  }
  first <- x[[1L]]
  basis <- .pf_or(attr(x, "basis"), first$basis)
  execution <- first$source$execution
  cat("<effect_sampling_covariance_batch>\n", sep = "")
  cat("  basis:        ", basis, "\n", sep = "")
  cat("  measurements: ", length(x), "\n", sep = "")
  cat("  ", if (identical(basis, "rdm")) "distances:   " else "coordinates: ",
    " ", first$dimension, "\n", sep = "")
  if (!is.null(execution)) {
    cat("  execution:    ", execution$route, " / ",
      execution$residual_strategy, "\n", sep = "")
    cat("  shared residual statistics: ",
      if (isTRUE(execution$shared_residual_statistics)) "yes" else "no",
      "\n", sep = "")
  }
  invisible(x)
}

# The quadratic noise term is only as good as the residual covariance behind
# it, so the two numbers that decide that -- the residual degrees of freedom
# and the effective dimension the support spends its variance on -- are
# printed rather than left in `$source` for the curious.
.pf_residual_noise <- function(source) {
  estimator <- source$noise_trace_estimator
  if (is.null(estimator)) return("unrecorded")
  if (!identical(estimator, "wishart_unbiased_quadratic")) {
    return("known residual covariance; no plug-in correction")
  }
  paste0(
    source$residual_df, " residual df, ",
    .pf_num(source$residual_effective_dimension), " effective dimensions",
    " (Wishart-corrected tr(Sigma^2))"
  )
}

# Query-first geometry plans -------------------------------------------------

# The plan is the object a user is told to keep with an analysis record, so the
# default block says what was declared, in the words the declaration used. The
# compiler's own vocabulary -- how the metric is lowered, how large a full
# materialization would be -- is real but is not what the record is for, so it
# is shown only on request.
.geometry_plan_lines <- function(x, detail = FALSE) {
  metric <- x$metric_schedule$metric
  learned <- identical(
    x$metric_schedule$kind, "learned_local_before_frame"
  )
  metric_status <- if (learned) {
    sprintf("learned %s, trained %s (derived per support)",
      x$metric_schedule$schedule$recipe$kind,
      x$metric_schedule$schedule$training_policy$kind)
  } else if (is.null(metric)) {
    "implicit identity"
  } else if (isTRUE(metric$capabilities$learned_frozen)) {
    "fixed (learned, frozen before evaluation)"
  } else {
    "fixed"
  }
  axis <- attr(x$pairing, "generalizes_over", exact = TRUE)
  independence <- attr(x$pairing, "independence", exact = TRUE)
  lines <- c(
    "<effect_geometry_plan>",
    sprintf("  effects:      %d x %d", x$logical_shape[[1L]],
      x$logical_shape[[2L]]),
    sprintf("  measurements: %d", x$measurements),
    sprintf("  features:     %d", x$task$left_relation$n_features),
    sprintf("  metric:       %s", metric_status),
    sprintf("  generalizes:  %d partition pairs%s, endpoints %s",
      nrow(x$pairing),
      if (is.null(axis)) " (axis undeclared)" else paste0(" over ", axis),
      independence),
    "  execution:    query-first, in memory"
  )
  if (isTRUE(detail)) {
    payload <- structure(x$dense_payload_bytes, class = "object_size")
    lines <- c(lines,
      sprintf("  lowering:     %s", x$lowering),
      sprintf("  dense payload: %s", format(payload, units = "auto")),
      sprintf("  codec:        %s", x$codec),
      sprintf("  plan id:      %s", .pf_sig(x$scientific_plan_id)),
      sprintf("  signature:    %s", .pf_sig(x$signature))
    )
  }
  if (learned) {
    # The one plan kind that reads before it returns. Say so on the object
    # rather than only in the receipt.
    return(c(lines,
      "  state:        residual statistics accumulated; no effect read",
      "  next:         crossnobis(plan, weights)"))
  }
  c(lines,
    "  state:        nothing computed yet",
    "  next:         contrast_energy(plan, weights), rdm(plan), rsa(plan)")
}

#' @export
print.effect_geometry_plan <- function(x, detail = FALSE, ...) {
  .validate_geometry_plan(x, deep = FALSE)
  cat(paste0(.geometry_plan_lines(x, detail = detail), "\n"), sep = "")
  invisible(x)
}

# Contrast weights in relation order, as `face 1, house -1, body 0, tool 0`.
# Named so that a positionally supplied contrast shows the alignment used.
.pf_weights <- function(weights, max = 6L) {
  if (is.null(weights) || !length(weights)) {
    return("none")
  }
  labels <- names(weights)
  if (is.null(labels) || anyNA(labels) || !all(nzchar(labels))) {
    labels <- paste0("[", seq_along(weights), "]")
  }
  .pf_set(paste0(labels, " ",
    format(signif(as.numeric(weights), 4L), trim = TRUE)), max = max)
}

# Location transports ---------------------------------------------------------
#
# A transport is estimand-bearing, so its one-screen summary has to carry the
# three things that change what a group number means -- the semantics, the
# unmapped territory, and how the operator was built -- alongside the shape.
# Density additionally names its denominator, because the declared row mass is
# what the ratio is per unit of.

.pf_transport_shape <- function(x) {
  paste0(nrow(x$matrix), " native -> ", nrow(x$group_index),
    " group + sink")
}

.pf_transport_provenance <- function(provenance) {
  cross_fit <- if (is.null(provenance$cross_fit)) {
    "none"
  } else {
    .pf_set(provenance$cross_fit, max = 3L)
  }
  paste0(.pf_or(provenance$method, "undeclared"), " (cross-fit: ", cross_fit,
    ")")
}

# Population plans -------------------------------------------------------------
#
# The plan's one-screen summary has to carry every choice that decides what a
# group number means, because each of them is a plan-identity field: the
# transport semantics, the budget normalization, the fit and its evaluation
# order. Coverage joins them -- a transport can raise a consensus share simply
# by sinking the nodes that disagree, so unmapped territory is printed beside
# the shape rather than left to a diagnostic call the reader may not make.

.pf_population_subjects <- function(x) {
  paste0(length(x$subjects), " (", .pf_set(names(x$subjects), max = 4L), ")")
}

.pf_population_sink <- function(x) {
  share <- x$subject_index$sink_territory
  reached <- sum(share > 0)
  if (!reached) {
    return("empty in every subject (full native coverage)")
  }
  paste0("present in ", reached, " of ", nrow(x$subject_index),
    " subjects, worst ", sprintf("%.1f%%", 100 * max(share)),
    " of territory")
}

.pf_population_model <- function(x) {
  paste0(x$model$formula_text, " -> ",
    .msg_count(length(x$model$columns), "column"), ", rank ", x$model$rank)
}

.pf_population_fit <- function(x) {
  paste0(toupper(x$fit$kind), " (", gsub("_", "-", x$fit$weights),
    " weights), ", gsub("_", " ", x$fit$evaluation_order))
}

# Population results -----------------------------------------------------------
#
# `population-form-v1` section 8.1 makes one line of this print method
# normative: a transported component view must state the native frame family
# it is a ledger of and the transport identity that carried it, so a reader
# never has to infer which frame a "coherent" number belongs to. The closing
# note is the other half of the same obligation --- the additive identity does
# survive transport, and the failure is only in reading the coherent ledger as
# a group-node common mode.

.pf_population_ledger <- function(x) {
  paste0(x$ledger, " (component \"", x$component, "\")")
}

.pf_population_frame <- function(x) {
  frame <- x$receipt$frame
  family <- if (all(is.na(frame$family))) "undeclared" else
    .pf_set(frame$family, max = 3L)
  scale <- if (all(is.na(frame$scale))) NULL else
    paste0(", scale ", .pf_num(frame$scale))
  paste0(family, scale, ", ", .pf_set(frame$normalization, max = 2L))
}

.pf_population_carrier <- function(x) {
  audit <- x$receipt$subjects
  cross_fit <- audit$cross_fit[!is.na(audit$cross_fit)]
  paste0(x$semantics, ", ", .pf_set(unique(audit$provenance), max = 3L),
    ", cross-fit ", if (length(cross_fit)) {
      .pf_set(unique(cross_fit), max = 2L)
    } else {
      "not declared"
    })
}

.pf_population_budget <- function(x) {
  budget <- x$receipt$budget
  if (!isTRUE(budget$asserted)) {
    return("not asserted (density conserves no total)")
  }
  paste0("preserved, worst relative deviation ",
    format(signif(budget$max_relative_deviation, 3L)), " against ",
    format(budget$tolerance))
}

.pf_population_normalization <- function(x) {
  record <- x$receipt$normalization
  if (identical(record$mode, "none")) {
    return("none (mean subject ledger, native evidence units)")
  }
  excluded <- sum(!record$admitted)
  paste0(record$mode, " (", gsub("_", " ", record$budget_estimate),
    ", floor ", format(record$budget_floor), ")",
    if (excluded) paste0(", ", excluded, " subject-query cells refused")
    else "")
}

# Section 7 requires the two layers to be reported apart, and so does the
# printed block: they are two keys, never one line with a semicolon in it,
# because a reader who sees them joined will read the second as a refinement of
# the first. The between-subject layer is always available on a query-bank
# result --- the participants are its replicates and it needs no error channel
# --- which is why this line no longer disappears when `uncertainty =` was not
# supplied. `uncalibrated` is on it unconditionally: a reader who sees a
# standard error offered without that word will supply a p-value from memory.
.pf_population_uncertainty <- function(x) {
  if (is.null(x$uncertainty)) {
    return(NULL)
  }
  if (x$residual_df < 1L) {
    return("between-subject SE not estimable (saturated group model)")
  }
  paste0("between-subject SE, df ", x$residual_df, " (uncalibrated)")
}

# The within-subject line, said as what it is rather than as a status word. A
# refusal names the thing that is missing; an admission names how much of the
# transport it covers, because "exact" over two of forty columns is a different
# report from "exact" over all forty.
.pf_population_within <- function(x) {
  within <- x$uncertainty$within
  if (is.null(within)) {
    return(NULL)
  }
  if (is.null(within$admitted)) {
    reason <- within$refusal$reasons[[1L]]
    return(paste0("refused: ", if (grepl("^same_data_ratio", reason)) {
      "unit_budget divisor has no standard error"
    } else if (grepl("^native_node_labels_unaligned", reason)) {
      "transport and covariance name nodes differently"
    } else {
      gsub("_", " ", reason)
    }))
  }
  if (!within$admitted_columns) {
    return("refused: no group column is fed by a single native row")
  }
  paste0("exact at ", within$admitted_columns, " of ", within$columns,
    " node-participant columns")
}

# The readout axis, said the way each basis reads it: a query bank is `K`
# declared contrasts, a complete form is the whole packed geometry and declares
# none. Keeping the two apart in the printed line is the same discipline
# section 8.1 imposes on the ledger name --- a reader must never have to guess
# whether a number came from a contrast someone chose.
.pf_population_readout <- function(x) {
  if (identical(x$basis, "complete_form")) {
    q <- length(x$effects)
    return(paste0("complete ", q, "x", q, " form, ",
      nrow(x$coordinates), " packed coordinates (", .pf_set(
        x$coordinates$coordinate, max = 3L), ")"))
  }
  paste0(nrow(x$queries), " (", .pf_set(rownames(x$queries), max = 3L), ")")
}

# The streamed route's bound, printed because it is the acceptance claim: the
# dense stack the route refuses to build is quoted beside the one it built.
.pf_population_streaming <- function(x) {
  streaming <- x$receipt$streaming
  if (is.null(streaming)) return(NULL)
  paste0("coordinate tile ", streaming$coordinate_tile, " of ",
    streaming$packed_width, ", ",
    .msg_count(streaming$passes_per_subject, "pass", "passes"),
    "/participant; peak group stack ",
    format(streaming$group_stack_doubles, big.mark = ",", scientific = FALSE),
    " doubles vs ",
    format(streaming$refused_dense_doubles, big.mark = ",",
      scientific = FALSE), " dense")
}

#' @export
print.effect_population_result <- function(x, ...) {
  .validate_population_result(x)
  audit <- x$receipt$subjects
  form <- identical(x$basis, "complete_form")
  # The readout line is labelled by what it holds: a query bank names the
  # contrasts someone chose, a complete form declares none. One label for both
  # would let a reader of a form believe a contrast was picked for them.
  fields <- list(
    ledger = .pf_population_ledger(x),
    subjects = paste0(nrow(audit), " (", .pf_set(audit$subject, max = 4L),
      "), residual df ", x$residual_df),
    `group nodes` = paste0(nrow(x$index) - 1L, " + sink")
  )
  fields[[if (form) "readout" else "queries"]] <- .pf_population_readout(x)
  .pf_emit("effect_population_result", c(fields, list(
    streaming = .pf_population_streaming(x),
    frame = .pf_population_frame(x),
    transport = .pf_population_carrier(x),
    normalization = .pf_population_normalization(x),
    fit = paste0(toupper(x$receipt$fit$kind), " (",
      gsub("_", "-", x$receipt$fit$weights), " weights), ",
      gsub("_", " ", x$receipt$evaluation_order)),
    budget = .pf_population_budget(x),
    uncertainty = .pf_population_uncertainty(x),
    within = .pf_population_within(x),
    unresolved = if (x$receipt$unresolved_columns) {
      paste0(x$receipt$unresolved_columns, " node-",
        if (form) "coordinate" else "query", " cells not estimated")
    },
    estimand = .pf_sig(x$scientific_plan_id)
  )))
  if (!identical(x$component, "total")) {
    .pf_note(paste0(
      "native_coherent_ledger + native_configuration_ledger = ",
      "transported_total holds exactly; this is native-node ", x$component,
      " evidence carried to a group location, not a group-node common mode."
    ))
  }
  if (form) {
    .pf_note(paste0(
      "No participant-level arrays: indexed by participant, group node and ",
      "packed coordinate they are the dense stack this route streams in ",
      "order not to build. Read them through estimate_population() with the ",
      "contrasts you want."
    ))
  }
  cat(sprintf("  %-15s%s\n", "next:", if (form) {
    "as.data.frame(x), x$coefficient_forms[node, term, ]"
  } else {
    "as.data.frame(x), x$coefficients[, , term], population_uncertainty(x)"
  }))
  invisible(x)
}

# Population uncertainty -------------------------------------------------------
#
# The printed record's job is to make three things unmissable: that there are
# two layers, that they are not added together, and that the `t` is
# uncalibrated for real data whatever the recorded simulation measured. The
# summary numbers are ranges rather than tables --- a printer never emits a
# whole array -- and the closing note is the separation itself.

.pf_uncertainty_range <- function(values) {
  finite <- values[is.finite(values)]
  if (!length(finite)) {
    return("all NA")
  }
  paste0(.pf_num(min(finite)), " to ", .pf_num(max(finite)),
    if (length(finite) < length(values)) {
      paste0(" (", length(values) - length(finite), " NA)")
    })
}

.pf_uncertainty_max <- function(values) {
  finite <- values[is.finite(values)]
  if (!length(finite)) {
    return("all NA")
  }
  paste0(.pf_num(max(finite)),
    if (length(finite) < length(values)) {
      paste0(" (", length(values) - length(finite), " NA)")
    })
}

.pf_uncertainty_within <- function(x) {
  within <- x$within
  if (is.null(within)) {
    return(paste0("absent (estimate_population() was given no `uncertainty`)"))
  }
  if (is.null(within$admitted)) {
    return(paste0("refused: ", .pf_set(within$refusal$reasons, max = 2L)))
  }
  paste0(within$admitted_columns, " of ", within$columns,
    " node-participant columns exact, variance ",
    .pf_uncertainty_range(within$variance[within$admitted[
      , rep(seq_len(dim(within$variance)[[3L]]),
        each = dim(within$variance)[[2L]]), drop = FALSE]]))
}

#' @export
print.effect_population_uncertainty <- function(x, ...) {
  .validate_population_uncertainty(x)
  .pf_emit("effect_population_uncertainty", list(
    ledger = paste0(x$ledger, " (component \"", x$receipt$component, "\")"),
    `group nodes` = paste0(nrow(x$index) - 1L, " + sink"),
    queries = paste0(nrow(x$queries), " (",
      .pf_set(rownames(x$queries), max = 3L), ")"),
    terms = .pf_set(x$term, max = 4L),
    `between-subject` = paste0("SE ",
      .pf_uncertainty_range(x$between$se), ", df ",
      x$between$residual_df, ", |t| up to ",
      .pf_uncertainty_max(abs(x$between$t))),
    interval = paste0(format(100 * x$between$level), "% nominal, ",
      x$between$calibration),
    `within-subject` = .pf_uncertainty_within(x),
    normalization = x$normalization,
    estimand = .pf_sig(x$scientific_plan_id)
  ))
  .pf_note(paste0(
    "The two layers are reported separately and are never pooled: a ",
    "within-subject sampling variance and a between-subject residual ",
    "variance answer different questions, and their sum answers neither."
  ))
  .pf_note(paste0(
    "The t is UNCALIBRATED for real data. Measured against t_df, the nominal ",
    "95% interval covers 0.948 to 0.952 of the time under a correctly ",
    "specified group model, and 0.885 to 0.923 when the participants' noise ",
    "is linked to the group covariates -- a nominal 5% test rejecting a true ",
    "null 11.5% of the time at N = 24. See ",
    "benchmarks/POPULATION-NULL-COVERAGE.md."
  ))
  cat(sprintf("  %-15s%s\n", "next:",
    "as.data.frame(x, layer = \"between\"), x$between$se[, , term]"))
  invisible(x)
}

# Population views -------------------------------------------------------------
#
# `population-form-v1` section 8.1 makes two of these lines normative and one
# of them forbidden. Required: the native frame family the ledger belongs to
# (family, scale, normalization) and the transport that carried it (semantics,
# provenance, cross-fit status), so a reader never has to infer which frame a
# transported number is a ledger of. Forbidden: the bare words `coherent` and
# `configuration` anywhere in the output. The closing note therefore names the
# ledger rather than the component -- the additive identity does survive
# transport, and the only false reading is of the ledger as a group-node
# common mode.
#
# The basis line is the other half of what a reader needs: it says whether the
# number was *selected* out of the estimated bank, *recombined* from several
# of its columns, or read off an assembled coefficient form, which is the
# difference between three routes that are exact for three different reasons.

.pf_population_view_names <- c(
  contrast_energy = "contrast energy",
  rdm = "squared distances",
  rsa = "RSA coefficients",
  ledger = "estimated ledger",
  contribution = "aggregated ledger"
)

.pf_population_view_kind <- function(x) {
  position <- match(x$view, names(.pf_population_view_names))
  if (is.na(position)) x$view else unname(.pf_population_view_names[[position]])
}

.pf_population_view_basis <- function(x) {
  basis <- x$receipt$basis
  route <- gsub("_", " ", basis$route)
  if (identical(basis$kind, "complete_form")) {
    return(paste0(route, ", ", basis$width, " packed coordinates"))
  }
  paste0(route, " over ", .msg_count(basis$width, "estimated query",
    "estimated queries"), " (", .pf_set(basis$queries, max = 3L), ")")
}

.pf_population_view_rows <- function(x) {
  if (identical(x$view, "contribution")) {
    # The sink is a row of the aggregate and not one of the group nodes the
    # territories were built from: section 3.3 is that it is never a value at
    # a location, and counting it as one here would put it back.
    return(paste0(nrow(x$index) - 1L, " territories + sink, from ",
      .msg_count(sum(x$index$n_nodes[!x$index$sink]), "group node")))
  }
  paste0(nrow(x$index) - 1L, " + sink")
}

#' @export
print.effect_population_view <- function(x, ...) {
  .validate_population_view(x)
  .pf_emit("effect_population_view", list(
    view = .pf_population_view_kind(x),
    ledger = x$ledger,
    term = x$term,
    `group nodes` = .pf_population_view_rows(x),
    columns = paste0(ncol(x$values), " (",
      .pf_set(as.character(x$columns$column), max = 4L), ")"),
    frame = .pf_population_frame(x),
    transport = .pf_population_carrier(x),
    normalization = .pf_population_normalization(x),
    basis = .pf_population_view_basis(x),
    budget = .pf_population_budget(x),
    aggregation = if (identical(x$view, "contribution")) {
      paste0("by ", x$metadata$aggregation$aggregated_by, ", ",
        if (x$native_ledger) "frame-relative, " else "",
        "budget-exact; the sink is its own row")
    },
    estimand = .pf_sig(x$scientific_plan_id)
  ))
  if (x$native_ledger) {
    .pf_note(paste0(
      "native_coherent_ledger + native_configuration_ledger = ",
      "transported_total holds exactly; a ", x$ledger, " is native-node ",
      "evidence carried to a group location, not a group-node common mode."
    ))
  }
  cat(sprintf("  %-15s%s\n", "next:", if (identical(x$view, "contribution")) {
    "as.data.frame(x), x$values"
  } else {
    "as.data.frame(x), contribution(x, by = ...)"
  }))
  invisible(x)
}

# The record registry --------------------------------------------------------
#
# One entry per sealed value class, in the order the classes are built. A
# `fields` function returns the named list `.pf_emit()` renders under a
# `<class>` header; an `inline` function returns the pieces `.pf_inline()`
# joins into a one-line `<class: a, b>`. Either may instead be the string
# "capabilities", routing to the two shared capability renderers above.

.pf_records <- list(
  effect_source_capabilities = list(
    inline = "capabilities",
    fields = function(x) list(
      block_read = .pf_yn(x$block_read),
      reopenable = .pf_yn(x$reopenable),
      thread_safe = .pf_yn(x$thread_safe),
      revision = .pf_sig(x$stable_revision)
    )
  ),
  effect_error_capabilities = list(
    inline = "capabilities",
    emit = "capabilities"
  ),
  effect_metric_capabilities = list(
    inline = "capabilities",
    fields = function(x) {
      granted <- .pf_capability_flags(x, drop = "role")
      list(
        role = .pf_or(x$role),
        granted = .pf_set(names(granted)[granted == "yes"], max = 4L,
          empty = "none"),
        withheld = .pf_set(names(granted)[granted == "no"], max = 4L,
          empty = "none")
      )
    }
  ),
  effect_measurement_capabilities = list(
    inline = "capabilities",
    emit = "capabilities"
  ),

  # Neural domains -------------------------------------------------------------
  effect_domain = list(
    inline = function(x) c(x$id,
      paste0(x$n_features, " ", x$kind, " features")),
    fields = function(x) {
      .validate_domain(x)
      coordinates <- if (is.null(x$coordinates)) {
        "none (abstract)"
      } else {
        paste0(.pf_dim(x$coordinates), " (",
          .pf_set(x$coordinate_units, max = 3L), ")")
      }
      list(
        id = x$id,
        kind = x$kind,
        features = x$n_features,
        "feature ids" = .pf_set(x$feature_ids, max = 4L),
        coordinates = coordinates,
        geometry = .pf_sig(x$geometry_signature),
        `next` = "compile_frame(whole_brain(), domain) to place a frame on it"
      )
    }
  ),
  effect_domain_reference = list(
    inline = function(x) c(x$id,
      paste0(x$n_features, " features")),
    fields = function(x) {
      .validate_domain_reference(x)
      list(
        id = x$id,
        features = x$n_features,
        units = .pf_set(x$coordinate_units, max = 3L),
        geometry = .pf_sig(x$geometry_signature),
        signature = .pf_sig(x$signature),
        state = "identity only; the referenced domain carries the coordinates"
      )
    }
  ),

  # Frames ---------------------------------------------------------------------
  effect_frame_spec = list(
    inline = function(x) {
      detail <- if (is.null(x$radius)) {
        NULL
      } else {
        paste0("radius ", .pf_num(x$radius))
      }
      c(x$kind, detail,
        paste0(x$normalization, " normalization"))
    },
    fields = function(x) list(
      kind = x$kind,
      radius = if (is.null(x$radius)) NULL else .pf_num(x$radius),
      normalization = x$normalization,
      state = "unplaced; call compile_frame(spec, domain) to bind a domain"
    )
  ),
  effect_frame_family_spec = list(
    inline = function(x) c(
      paste0(length(x$members), " members"),
      paste0(x$normalization, " normalization")),
    fields = function(x) list(
      kind = x$kind,
      normalization = x$normalization,
      members = .pf_frame_family_members(x),
      alpha = paste0("sums to ", .pf_num(sum(x$alpha)))
    )
  ),
  effect_frame = list(
    inline = function(x) c(paste0(nrow(x$index), " nodes"),
      x$representation, x$specification$kind),
    fields = function(x) {
      weights <- x$weights
      specification <- x$specification
      placement <- paste0(specification$kind,
        if (is.null(specification$radius)) {
          ""
        } else {
          paste0(", radius ", .pf_num(specification$radius), " ",
            .pf_or(x$domain$coordinate_units, ""))
        })
      list(
        representation = x$representation,
        nodes = nrow(x$index),
        features = x$domain$n_features,
        specification = placement,
        members = if (identical(specification$kind, "frame_family")) {
          .pf_frame_family_members(specification)
        } else {
          NULL
        },
        normalization = x$normalization,
        weights = paste0(.pf_dim(weights), ", ", .pf_stored(weights),
          " stored"),
        domain = paste0(x$domain_id, " (", x$domain_kind, ")"),
        fixed = paste0(.pf_yn(x$fixed), "; locally estimated ",
          .pf_yn(x$locally_estimated)),
        `next` = "plan_geometry(relation, frame, pairing)"
      )
    }
  ),
  effect_support_index = list(
    inline = function(x) c(paste0(length(x$node_ids), " nodes"),
      paste0(length(x$members), " memberships"), x$construction$kind),
    fields = function(x) {
      .validate_support_index(x, deep = FALSE)
      cost <- x$cost
      sizes <- cost$support_size
      construction <- x$construction
      list(
        nodes = length(x$node_ids),
        construction = paste0(construction$kind, " via ",
          construction$provider),
        radius = if (is.null(construction$radius)) {
          NULL
        } else {
          paste0(.pf_num(construction$radius), " ",
            .pf_or(construction$coordinate_units, ""))
        },
        memberships = length(x$members),
        "support size" = paste0("min ", .pf_num(sizes[["min"]]), ", median ",
          .pf_num(sizes[["median"]]), ", max ", .pf_num(sizes[["max"]])),
        "overlapping pairs" = if (is.na(cost$pair_pattern_stored_nnz)) {
          "not materialized"
        } else {
          cost$pair_pattern_stored_nnz
        },
        domain = x$domain$id,
        signature = .pf_sig(x$signature)
      )
    }
  ),
  effect_support_cost = list(
    inline = function(x) c(paste0(x$nodes, " nodes"),
      paste0(x$support_memberships, " memberships")),
    fields = function(x) list(
      nodes = .pf_num(x$nodes),
      features = .pf_num(x$features),
      memberships = .pf_num(x$support_memberships),
      "pair pattern" = if (is.na(x$pair_pattern_stored_nnz)) {
        "not materialized"
      } else {
        paste0(.pf_num(x$pair_pattern_stored_nnz),
          " stored of ", .pf_num(x$pair_pattern_nnz))
      },
      "diagonal metric" = paste0(.pf_num(x$diagonal_metric_entries),
        " entries"),
      "dense metric" = paste0(.pf_num(x$dense_metric_entries), " entries, ",
        .pf_bytes(x$materialized_dense_metric_bytes)),
      "largest local" = .pf_bytes(x$largest_local_metric_bytes),
      structural = .pf_bytes(x$estimated_structural_bytes),
      state = "cost projection only; nothing is materialized"
    )
  ),

  # Relations and their parts --------------------------------------------------
  effect_relation = list(
    inline = function(x) c(paste0(length(x$effects), " effects"),
      paste0(length(x$partitions), " partitions"),
      paste0(x$n_features, " features")),
    fields = function(x) {
      .validate_relation(x, deep = FALSE)
      kinds <- vapply(x$sources,
        function(source) .pf_or(source$kind, "unknown"), character(1))
      estimated <- is.null(x$extractors) ||
        all(vapply(x$extractors, is.null, logical(1)))
      list(
        effects = paste0(length(x$effects), " (",
          .pf_set(x$effects, max = 4L), ")"),
        partitions = paste0(length(x$partitions), " (",
          .pf_set(x$partitions, max = 4L), ")"),
        features = x$n_features,
        domain = x$domain_id,
        basis = .pf_or(x$effect_space$basis_id, "unspecified"),
        sources = paste0(length(x$sources), " ",
          .pf_set(unique(kinds), max = 2L), " (unread)"),
        extraction = if (estimated) {
          "none; sources are already effect-by-feature estimates"
        } else {
          paste0(length(x$extractors), " extractors, ",
            .pf_or(x$extractors[[1L]]$estimator, "unknown"))
        },
        state = "lazy; sources are read only by neural feature block",
        `next` = "plan_geometry(relation, frame, pairing)"
      )
    }
  ),
  effect_space = list(
    inline = function(x) c(paste0(length(x$coordinates), " coordinates"),
      x$basis_id),
    fields = function(x) {
      .validate_effect_space(x)
      unit <- unique(as.character(x$units))
      scale <- unique(as.numeric(x$scale))
      list(
        coordinates = paste0(length(x$coordinates), " (",
          .pf_set(x$coordinates, max = 4L), ")"),
        basis = .pf_or(x$basis_id, "unspecified"),
        units = if (length(unit) == 1L) unit else .pf_set(unit, max = 3L),
        scale = if (length(scale) == 1L) .pf_num(scale) else
          .pf_num(range(scale)),
        signature = .pf_sig(x$signature)
      )
    }
  ),
  effect_extractor = list(
    inline = function(x) c(paste0(length(x$effects), " effects"),
      paste0("from ", x$n_observations, " observations"), x$estimator),
    fields = function(x) {
      .validate_effect_extractor(x)
      diagnostics <- x$diagnostics
      list(
        effects = paste0(length(x$effects), " (",
          .pf_set(x$effects, max = 4L), ")"),
        observations = x$n_observations,
        map = .pf_dim(x$map),
        estimator = x$estimator,
        solver = paste0(.pf_or(diagnostics$solver, "unknown"), ", rank ",
          .pf_or(diagnostics$rank, "unknown"),
          if (isTRUE(diagnostics$rank_deficient)) " (rank deficient)" else ""),
        estimability = paste0("max error ",
          .pf_num(max(c(0, diagnostics$estimability_error)))),
        whitener = .pf_or(diagnostics$observation_whitener$kind, "none")
      )
    }
  ),
  effect_response_source = list(
    inline = function(x) c(x$kind, .pf_shape(x$dim)),
    fields = function(x) list(
      kind = x$kind,
      dim = paste0(.pf_shape(x$dim), " (observations x features)"),
      access = .pf_or(x$descriptor$access, "coordinator"),
      revision = .pf_sig(x$descriptor$stable_revision),
      state = "unread; blocks are pulled on demand by feature index"
    )
  ),
  effect_source_descriptor = list(
    inline = function(x) c(x$kind, .pf_shape(x$dim)),
    fields = function(x) list(
      kind = x$kind,
      dim = .pf_shape(x$dim),
      access = .pf_or(x$access, "coordinator"),
      revision = .pf_sig(x$stable_revision),
      spec = .pf_set(names(x$spec), max = 4L, empty = "none")
    )
  ),
  effect_error_model = list(
    inline = function(x) c(x$kind,
      paste0(x$residual_df, " residual df"), x$sampling_unit),
    fields = function(x) {
      .validate_error_model(x, deep = FALSE)
      list(
        kind = x$kind,
        assumptions = .pf_set(x$assumptions, max = 1L),
        "sampling unit" = x$sampling_unit,
        "residual df" = x$residual_df,
        "effect covariance" = .pf_dim(x$effect_covariance),
        scale = x$scale_convention,
        whitener = .pf_or(x$observation_whitener$kind, "none"),
        residuals = if (is.null(x$residual_source)) {
          "not retained"
        } else {
          paste0("available, ", .pf_shape(x$residual_source$dim), " (unread)")
        },
        signature = .pf_sig(x$signature)
      )
    }
  ),

  # Pairings -------------------------------------------------------------------
  #
  # `effect_pairing` subclasses `data.frame`, so only print() is defined here:
  # a format() method would shadow `format.data.frame` for internal callers.
  effect_pairing = list(
    fields = function(x) list(
      pairs = nrow(x),
      left = .pf_set(unique(x$left), max = 4L),
      right = .pf_set(unique(x$right), max = 4L),
      weights = if (length(unique(x$weight)) == 1L) {
        paste0("equal (", .pf_num(x$weight[[1L]]), ")")
      } else {
        paste0("unequal, ", .pf_num(range(x$weight)))
      },
      directed = .pf_yn(attr(x, "directed", exact = TRUE)),
      "self pairs" = .pf_or(attr(x, "self_pairs", exact = TRUE), "forbid"),
      independence = .pf_or(attr(x, "independence", exact = TRUE),
        "undeclared"),
      generalizes = .pf_or(attr(x, "generalizes_over", exact = TRUE),
        "axis undeclared"),
      estimate = .pf_or(attr(x, "estimate", exact = TRUE), "unspecified")
    )
  ),

  # Metrics --------------------------------------------------------------------
  effect_neural_metric = list(
    inline = function(x) c(x$role,
      paste0(length(x$support), " supported features"), x$estimation),
    fields = function(x) {
      .validate_neural_metric(x, deep = FALSE)
      capabilities <- x$capabilities
      list(
        role = x$role,
        domain = x$domain$id,
        support = paste0(length(x$support), " of ", x$domain$n_features,
          " features"),
        value = paste0(.pf_dim(x$value), " (not shown)"),
        estimation = x$estimation,
        definite = paste0("psd ", .pf_yn(capabilities$positive_semidefinite),
          ", pd ", .pf_yn(capabilities$positive_definite)),
        inverse = .pf_or(x$inverse_representation$kind, "none"),
        tolerance = .pf_num(x$tolerance),
        signature = .pf_sig(x$signature)
      )
    }
  ),
  effect_metric_recipe = list(
    inline = function(x) c(x$kind, x$role, "unfitted"),
    fields = function(x) {
      .validate_metric_recipe(x)
      list(
        kind = x$kind,
        role = x$role,
        domain = x$domain$id,
        support = if (isTRUE(x$support_local)) {
          "local (per node)"
        } else {
          "whole domain"
        },
        hyperparameters = .pf_set(names(x$hyperparameters), max = 4L),
        signature = .pf_sig(x$signature),
        state = "recipe only; the metric is fit when a plan schedules it"
      )
    }
  ),
  effect_metric_schedule = list(
    inline = function(x) c(x$kind, x$scope, x$materialization),
    fields = function(x) list(
      role = x$role,
      kind = x$kind,
      scope = x$scope,
      composition = x$frame_composition,
      lowering = x$lowering,
      materialization = x$materialization,
      properties = paste0("feature additive ", .pf_yn(x$feature_additive),
        ", dense support ", .pf_yn(x$support_dense)),
      metric = if (identical(x$kind, "learned_local_before_frame")) {
        paste0("derived on demand from ", x$schedule$recipe$kind, " (frozen ",
          .pf_sig(x$schedule$signature), ")")
      } else if (is.null(x$metric)) {
        "implicit identity"
      } else {
        .pf_sig(x$metric_signature)
      },
      signature = .pf_sig(x$signature)
    )
  ),
  effect_frozen_metric_schedule = list(
    inline = function(x) c(x$recipe$kind,
      x$training_policy$kind, "frozen before evaluation"),
    fields = function(x) {
      .validate_frozen_metric_schedule(x, deep = FALSE)
      list(
        role = .pf_or(x$role),
        recipe = .pf_or(x$recipe$kind, "none"),
        training = .pf_or(x$training_policy$kind, "none"),
        support = paste0(length(x$support_index$node_ids), " nodes"),
        pairing = paste0(nrow(x$pairing), " partition pairs"),
        statistics = if (is.null(x$statistics)) {
          "none"
        } else {
          paste0(length(x$statistics$partitions), " residual partitions")
        },
        signature = .pf_sig(x$signature),
        state = "frozen; training excludes the evaluated partitions"
      )
    }
  ),
  effect_metric_training_policy = list(
    inline = function(x) c(x$kind),
    fields = function(x) {
      .validate_metric_training_policy(x)
      list(
        kind = x$kind,
        "uses evaluation residuals" =
          .pf_yn(x$includes_evaluation_residuals),
        assumption = .pf_or(x$assumption, "none"),
        justification = .pf_or(x$justification, "none recorded"),
        signature = .pf_sig(x$signature)
      )
    }
  ),

  # Policies and contracts -----------------------------------------------------
  effect_compute_policy = list(
    inline = function(x) c(paste0(x$workers, " workers"),
      x$process_backend),
    fields = function(x) list(
      workers = x$workers,
      backend = .pf_or(x$process_backend, "sequential"),
      "feature block" = if (is.null(x$block_features)) {
        "chosen by the kernel"
      } else {
        as.character(x$block_features)
      },
      workspace = if (is.null(x$workspace_bytes)) {
        "unbounded"
      } else {
        .pf_bytes(x$workspace_bytes)
      }
    )
  ),
  effect_memory_plan = list(
    inline = function(x) c(
      paste0(.pf_bytes(x$planned_workspace_bytes), " planned"),
      paste0(x$workers, " workers")),
    fields = function(x) {
      categories <- x$categories
      ranked <- sort(categories[categories > 0], decreasing = TRUE)
      list(
        workers = paste0(x$workers, " (", x$n_active, " active)"),
        persistent = .pf_bytes(x$persistent_workspace_bytes),
        "per active task" = .pf_bytes(x$task_workspace_per_active_bytes),
        modeled = .pf_bytes(x$modeled_workspace_bytes),
        planned = paste0(.pf_bytes(x$planned_workspace_bytes), " (safety ",
          .pf_num(x$safety_factor), "x)"),
        budget = if (is.null(x$budget_bytes)) {
          "none declared"
        } else {
          paste0(.pf_bytes(x$budget_bytes), ", fits ", .pf_yn(x$fits_budget))
        },
        largest = .pf_set(names(ranked), max = 3L),
        measured = if (is.null(x$measured_workspace_bytes)) {
          "not yet executed"
        } else {
          .pf_bytes(x$measured_workspace_bytes)
        },
        kind = x$prediction_kind
      )
    }
  ),
  effect_numerical_contract = list(
    inline = function(x) c(paste0("atol ", .pf_num(x$atol)),
      paste0("rtol ", .pf_num(x$rtol))),
    fields = function(x) {
      .validate_numerical_contract(x)
      list(
        atol = .pf_num(x$atol),
        rtol = .pf_num(x$rtol),
        scheduling = .pf_or(x$scheduling$guarantee, "unspecified"),
        "block partition" = .pf_or(x$block_partition$guarantee, "unspecified"),
        "cross platform" = .pf_or(x$cross_platform$guarantee, "unspecified"),
        bitwise = paste0("across blocking ",
          .pf_yn(x$bitwise_across_blocking), ", across platforms ",
          .pf_yn(x$bitwise_across_platforms))
      )
    }
  ),
  effect_numeric_agreement = list(
    inline = function(x) c(
      if (isTRUE(x$passed)) "passed" else "failed", x$guarantee),
    fields = function(x) list(
      result = if (isTRUE(x$passed)) "passed" else "failed",
      guarantee = x$guarantee,
      comparison = x$comparison,
      "max error" = .pf_num(x$max_absolute_error),
      "max allowed" = .pf_num(x$max_allowed_error),
      tolerance = paste0("atol ", .pf_num(x$atol), ", rtol ", .pf_num(x$rtol))
    )
  ),
  effect_partition_reducer = list(
    inline = function(x) c(x$kind, x$weight_convention),
    fields = function(x) {
      .validate_partition_reducer(x)
      list(
        kind = x$kind,
        weights = x$weight_convention,
        order = x$order
      )
    }
  ),
  effect_gaussian_covariance_model = list(
    inline = function(x) c(x$family,
      if (isTRUE(x$fixed)) "fixed" else "estimated"),
    fields = function(x) {
      .validate_gaussian_covariance_model(x)
      list(
        family = x$family,
        fixed = .pf_yn(x$fixed),
        signature = .pf_sig(x$signature),
        state = "joint Gaussian declared; information units are claimable"
      )
    }
  ),

  # Compiled tasks, sampling plans, and receipts -------------------------------
  effect_evidence_task = list(
    inline = function(x) c(
      paste0(nrow(x$ordered_partition_products), " partition products"),
      if (isTRUE(x$same_relation)) "one relation" else "two relations",
      x$stages$lowering),
    fields = function(x) {
      spaces <- x$spaces
      list(
        relations = if (isTRUE(x$same_relation)) {
          "one (left and right are the same)"
        } else {
          "two (cross-relation)"
        },
        "left space" = paste0(
          length(spaces$experimental_left$coordinates), " coordinates, ",
          .pf_or(spaces$experimental_left$basis_id, "unspecified")),
        "neural space" = .pf_or(spaces$neural_left$id, "unspecified"),
        "partition products" = nrow(x$ordered_partition_products),
        stages = .pf_set(x$stages$order, max = 3L),
        lowering = x$stages$lowering,
        reducer = .pf_or(x$semantic$reducer$kind, "none"),
        materialization = .pf_or(x$materialization$kind, "query_only"),
        "task id" = .pf_sig(x$task_id),
        state = "compiled contract; no neural data has been read"
      )
    }
  ),
  effect_sampling_record = list(
    inline = function(x) c(x$kind),
    fields = function(x) {
      flags <- .pf_capability_flags(x, drop = c("kind", "signature"))
      fields <- c(list(kind = x$kind), as.list(utils::head(flags, 8L)))
      fields[["signature"]] <- .pf_sig(x$signature)
      fields
    }
  ),
  effect_evidence_sampling_plan = list(
    inline = function(x) c(x$contract, x$sampling_axis,
      x$target$target),
    fields = function(x) {
      .validate_evidence_sampling_plan(x, deep = FALSE)
      unmet <- x$unavailable_reasons
      list(
        contract = x$contract,
        "sampling axis" = x$sampling_axis,
        target = paste0(.pf_or(x$target$target), " / ",
          .pf_or(x$target$policy)),
        "error channel" = .pf_or(x$error_channel$kind, "none"),
        metric = .pf_or(x$metric$status, "unknown"),
        partitions = .pf_or(x$partition$model, x$partition$kind),
        "spatial scope" = x$spatial_scope,
        operation = .pf_or(x$operation$operation, x$operation$kind),
        unmet = if (length(unmet)) .pf_set(unmet, max = 2L) else "none",
        signature = .pf_sig(x$signature)
      )
    }
  ),

  # A sampling covariance is one class over two coordinate bases, so there is
  # one print method and it reports the basis it is reading. The RDM basis
  # keeps its own block: a reader of crossvalidated distances wants the
  # measurement, the metric status, and the reminder that the law is local,
  # none of which the general evidence basis has to say.
  effect_sampling_covariance = list(
    inline = function(x) {
      if (identical(x$basis, "rdm")) {
        return(c(
          paste0(x$dimension, " distances"), x$plan$target$policy,
          "factorized"))
      }
      if (identical(x$basis, "query_bank")) {
        return(c(
          paste0(x$dimension, " queries"), x$plan$target$policy, "factorized"))
      }
      c(
        paste0(x$dimension, " coordinates"), x$plan$target$policy)
    }
  ),
  effect_execution_receipt = list(
    inline = function(x) c(x$completion_status,
      paste0(x$completed_task_count, "/", x$task_count, " tasks")),
    fields = function(x) {
      .validate_execution_receipt(x)
      list(
        status = paste0(x$completion_status, " (", x$completed_task_count, "/",
          x$task_count, " tasks)"),
        "plan id" = .pf_sig(x$scientific_plan_id),
        kernel = x$kernel_version,
        partitioning = x$task_partition_id,
        precision = x$precision,
        compute = paste0(x$compute$workers, " workers, ",
          .pf_or(x$compute$process_backend, "sequential")),
        memory = paste0(.pf_bytes(x$memory$planned_workspace_bytes),
          " planned"),
        tolerance = paste0("atol ", .pf_num(x$numeric_contract$atol), ", rtol ",
          .pf_num(x$numeric_contract$rtol)),
        blas = paste0(.pf_or(x$blas$requested_threads, "unrecorded"),
          " thread(s) requested"),
        domain = .pf_sig(x$domain_signature)
      )
    }
  ),
  effect_geometry_store = list(
    inline = function(x) c(.pf_shape(x$dim), x$representation),
    fields = function(x) list(
      dim = paste0(.pf_shape(x$dim), " (measurements x packed values)"),
      representation = x$representation,
      format = .pf_or(x$manifest$format, "unspecified"),
      `next` = "geometry_component(geometry, \"total\")"
    )
  ),
  effect_residual_pair_statistics = list(
    inline = function(x) c(
      paste0(length(x$pair_i), " feature pairs"),
      paste0(length(x$partitions), " partitions")),
    fields = function(x) {
      .validate_residual_pair_statistics(x, deep = FALSE)
      list(
        "feature pairs" = length(x$pair_i),
        partitions = paste0(length(x$partitions), " (",
          .pf_set(x$partitions, max = 3L), ")"),
        domain = x$domain$id,
        support = .pf_sig(x$support_index),
        "relation fit" = .pf_sig(x$relation_fit),
        signature = .pf_sig(x$signature),
        state = "sufficient statistics only; residuals are not retained"
      )
    }
  ),

  # Study facts ----------------------------------------------------------------
  effect_observations = list(
    inline = function(x) c(paste0(length(x$partitions), " partitions"),
      paste0(x$n_features, " features")),
    fields = function(x) {
      .validate_observations(x, deep = FALSE)
      counts <- vapply(x$indexes, function(index) {
        length(index$observation_id)
      }, integer(1))
      list(
        partitions = paste0(length(x$partitions), " (",
          .pf_set(x$partitions, max = 3L), ")"),
        observations = paste0(sum(counts), " total (",
          .pf_set(as.character(counts), max = 3L), " per partition)"),
        features = x$n_features,
        domain = x$domain$id,
        timing = if (any(vapply(x$indexes, function(index) {
          isTRUE(index$timing)
        }, logical(1)))) "resolved" else "none declared",
        id = .pf_sig(x$observations_id),
        state = "lazy; sources are unread until a relation is fitted",
        `next` = "study(observations, events, confounds, hierarchy)"
      )
    }
  ),
  effect_observation_index = list(
    inline = function(x) c(x$partition,
      paste0(length(x$observation_id), " observations")),
    fields = function(x) {
      .validate_observation_index(x)
      list(
        partition = x$partition,
        observations = length(x$observation_id),
        ids = .pf_set(x$observation_id, max = 4L),
        time = if (isTRUE(x$timing)) {
          paste0(.pf_num(min(x$time)), " to ", .pf_num(max(x$time)), " ",
            x$units)
        } else {
          "none declared"
        },
        signature = .pf_sig(x$signature)
      )
    }
  ),
  effect_events = list(
    inline = function(x) c(paste0(nrow(x$data), " events"),
      paste0(length(unique(x$data[[x$partition_column]])), " partitions")),
    fields = function(x) {
      .validate_events(x)
      columns <- setdiff(names(x$data), c(x$partition_column, x$event_id_column,
        x$onset_column, x$duration_column))
      list(
        events = nrow(x$data),
        partitions = .pf_set(unique(x$data[[x$partition_column]]), max = 3L),
        timing = if (isTRUE(x$timing)) {
          paste0("onset + duration in ", x$units)
        } else {
          "none declared"
        },
        attributes = .pf_set(columns, max = 4L),
        id = .pf_sig(x$events_id),
        `next` = "study_axis(study, \"<attribute>\") to model an event attribute"
      )
    }
  ),
  effect_observation_confounds = list(
    inline = function(x) c(paste0(nrow(x$data), " rows"),
      paste0(length(x$schema) - 2L, " regressors")),
    fields = function(x) {
      .validate_observation_confounds(x)
      regressors <- setdiff(names(x$data), c(x$partition_column,
        x$observation_id_column, x$censor_column))
      censored <- if (is.null(x$censor_column)) {
        "none declared"
      } else {
        paste0(sum(!x$data[[x$censor_column]]), " of ", nrow(x$data),
          " observations dropped by `", x$censor_column, "`")
      }
      list(
        rows = nrow(x$data),
        partitions = .pf_set(unique(x$data[[x$partition_column]]), max = 3L),
        regressors = paste0(length(regressors), " (",
          .pf_set(regressors, max = 3L), ")"),
        censor = censored,
        id = .pf_sig(x$confounds_id)
      )
    }
  ),
  effect_partition_hierarchy = list(
    inline = function(x) c(paste0(nrow(x$data), " partitions"),
      .pf_set(x$axes, max = 3L)),
    fields = function(x) {
      .validate_partition_hierarchy(x)
      widths <- vapply(x$levels, length, integer(1))
      list(
        leaf = x$leaf,
        partitions = nrow(x$data),
        axes = paste0(length(x$axes), " (", .pf_set(x$axes, max = 4L), ")"),
        levels = .pf_set(paste0(names(widths), " ", widths), max = 4L),
        signature = .pf_sig(x$signature),
        `next` = "study_axis(study, \"<axis>\") to resample over an axis"
      )
    }
  ),
  effect_observation_model = list(
    inline = function(x) c(.pf_or(x$kind),
      .pf_or(x$sampling_unit)),
    fields = function(x) list(
      kind = .pf_or(x$kind),
      "sampling unit" = .pf_or(x$sampling_unit),
      independence = .pf_or(x$independence, "undeclared"),
      whitener = .pf_or(x$whitener$kind, "identity"),
      assumptions = .pf_set(x$assumptions, max = 2L),
      training = .pf_sig(x$training_revision),
      id = .pf_sig(x$observation_model_id)
    )
  ),

  # Measurement spaces, frames, and edges --------------------------------------
  effect_measurement_space = list(
    inline = function(x) c(x$id,
      paste0(x$n_measurements, " measurements")),
    fields = function(x) {
      .validate_measurement_space(x)
      list(
        id = x$id,
        measurements = x$n_measurements,
        signature = .pf_sig(x$signature)
      )
    }
  ),
  effect_measurement_frame = list(
    inline = function(x) c(paste0(length(x$node_ids), " nodes"),
      paste0("from ", x$source_domain$id)),
    fields = function(x) {
      .validate_measurement_frame(x)
      widths <- vapply(x$legs, function(leg) nrow(leg$operator), integer(1))
      list(
        nodes = paste0(length(x$node_ids), " (", .pf_set(x$node_ids, max = 3L),
          ")"),
        "node widths" = if (length(unique(widths)) == 1L) {
          paste0(widths[[1L]], " each")
        } else {
          paste0(min(widths), " to ", max(widths))
        },
        domain = x$source_domain$id,
        operator = .pf_dim(x$frame_operator),
        coverage = paste0(sum(x$coverage$count > 0L), " of ",
          length(x$coverage$feature_ids), " domain features",
          if (isTRUE(x$coverage$complete)) " (complete)" else " (partial)"),
        injectivity = x$injectivity$construction,
        signature = .pf_sig(x$signature)
      )
    }
  ),
  effect_edge_frame = list(
    inline = function(x) c(paste0(nrow(x$edges$edges), " edges"),
      paste0(length(x$from_frame$node_ids), " x ",
        length(x$to_frame$node_ids), " nodes")),
    fields = function(x) {
      .validate_edge_frame(x)
      list(
        edges = paste0(NROW(x$edges$edges), " of ",
          length(x$from_frame$node_ids) * length(x$to_frame$node_ids),
          " possible"),
        from = paste0(length(x$from_frame$node_ids), " nodes (",
          .pf_set(x$from_frame$node_ids, max = 3L), ")"),
        to = paste0(length(x$to_frame$node_ids), " nodes (",
          .pf_set(x$to_frame$node_ids, max = 3L), ")"),
        weighted = .pf_yn(x$edges$weighted),
        signature = .pf_sig(x$signature),
        state = "only the listed pairs are computed"
      )
    }
  ),
  effect_pair_query = list(
    inline = function(x) c(x$kind, .pf_dim(x$operator),
      if (isTRUE(x$fixed)) "fixed" else "estimated"),
    fields = function(x) list(
      kind = x$kind,
      operator = .pf_dim(x$operator),
      left = paste0(length(x$left_space$coordinates), " coordinates, ",
        .pf_or(x$left_space$basis_id, "unspecified")),
      right = paste0(length(x$right_space$coordinates), " coordinates, ",
        .pf_or(x$right_space$basis_id, "unspecified")),
      fixed = .pf_yn(x$fixed),
      role = .pf_or(x$metadata$evidence_capability$role, "effect"),
      "sampling axis" =
        .pf_or(x$metadata$evidence_capability$sampling_axis, "not declared"),
      construction =
        .pf_or(x$metadata$evidence_capability$construction, "arbitrary")
    )
  ),

  # Measurement forms and coupling views ---------------------------------------
  effect_measurement_form = list(
    inline = function(x) c(paste0(nrow(x$block_index), " blocks"),
      x$storage, x$completeness),
    fields = function(x) {
      .validate_measurement_form(x, probe = FALSE)
      capabilities <- x$capabilities
      list(
        blocks = nrow(x$block_index),
        nodes = paste0(length(x$left_frame$node_ids), " x ",
          length(x$right_frame$node_ids)),
        storage = paste0(x$storage, " (", x$codec, ")"),
        completeness = paste0(x$completeness, ", edges ", x$edge_completeness),
        construction = .pf_or(capabilities$construction, "arbitrary"),
        claims = paste0("symmetric ", .pf_yn(capabilities$symmetric),
          ", psd ", .pf_yn(capabilities$guaranteed_psd),
          ", repeated variation ", .pf_yn(capabilities$repeated_variation)),
        capability = x$result_capability,
        `next` = "effect_coupling(form), covariance_coupling(form)"
      )
    }
  ),
  effect_coupling_result = list(
    inline = function(x) c(x$kind,
      paste0(nrow(x$edge_index), " edges"), x$terminology),
    fields = function(x) {
      .validate_coupling_result(x)
      list(
        kind = x$kind,
        terminology = x$terminology,
        edges = paste0(nrow(x$edge_index), ", ", x$edge_completeness),
        # `$values` is a named list of blocks for two kinds and a table for the
        # other five; `length()` on a data frame would report its column count.
        values = if (identical(.coupling_value_shape(x), "edge_blocks")) {
          paste0(length(x$values), " blocks")
        } else {
          paste0(nrow(x$values), " rows x ", ncol(x$values), " columns")
        },
        normalization = .pf_or(x$normalization_axis, "none"),
        summary = .pf_or(x$summary_axis, "none"),
        regularization = .pf_regularization(x$regularization),
        units = .pf_or(x$units, "not claimed"),
        stages = .pf_set(x$stage_order, max = 3L),
        signature = .pf_sig(x$signature)
      )
    }
  ),

  # Design models and their receipts -------------------------------------------
  effect_design_model = list(
    inline = function(x) c(paste0(length(x$designs), " partitions"),
      paste0(length(x$condition_space$coordinates), " conditions")),
    fields = function(x) {
      .validate_design_model(x)
      rows <- vapply(x$designs, nrow, integer(1))
      columns <- vapply(x$designs, ncol, integer(1))
      granted <- .pf_capability_flags(x$capabilities)
      list(
        partitions = paste0(length(x$partitions), " (",
          .pf_set(x$partitions, max = 3L), ")"),
        conditions = paste0(length(x$condition_space$coordinates), " (",
          .pf_set(x$condition_space$coordinates, max = 3L), ")"),
        designs = paste0(.pf_set(as.character(rows), max = 2L), " rows x ",
          .pf_set(as.character(unique(columns)), max = 2L), " columns"),
        solver = .pf_set(unique(as.character(x$solver)), max = 2L),
        guarantees = .pf_set(names(granted)[granted == "yes"], max = 3L),
        id = .pf_sig(x$design_model_id),
        `next` = "plan_relation(study, model, effects, observation_model)"
      )
    }
  ),
  effect_raw_design_model = list(
    inline = function(x) c(
      paste0(length(x$designs), " partitions"), "no symbolic model"),
    fields = function(x) {
      rows <- vapply(x$designs, nrow, integer(1))
      columns <- vapply(x$designs, ncol, integer(1))
      list(
        partitions = paste0(length(x$designs), " (",
          .pf_set(names(x$designs), max = 3L), ")"),
        designs = paste0(.pf_set(as.character(rows), max = 2L), " rows x ",
          .pf_set(as.character(unique(columns)), max = 2L), " columns"),
        solver = .pf_set(unique(as.character(x$solver)), max = 2L),
        state = "raw columns; coding invariance and row lineage are not claimed"
      )
    }
  ),
  effect_design_receipt = list(
    inline = function(x) c(x$partition,
      paste0("rank ", x$rank), paste0(x$residual_df, " residual df")),
    fields = function(x) {
      .validate_design_receipt(x)
      censoring <- x$censoring
      list(
        partition = x$partition,
        design = .pf_dim(x$design),
        coefficients = .pf_set(x$coefficient_axis, max = 4L),
        effects = paste0(nrow(x$lowered_target), " lowered from ",
          ncol(x$lowered_target), " coefficients"),
        solver = paste0(x$solver, " (policy ", x$solver_policy, "), rank ",
          x$rank),
        aliases = if (length(x$aliases)) {
          .pf_set(x$aliases, max = 3L)
        } else {
          "none"
        },
        "residual df" = x$residual_df,
        censoring = if (!isTRUE(censoring$declared)) {
          "none declared"
        } else {
          paste0(length(censoring$retained_rows), " of ", censoring$input_rows,
            " rows retained by `", censoring$column, "`")
        },
        whitener = .pf_or(x$observation_whitener$kind, "identity"),
        id = .pf_sig(x$design_receipt_id)
      )
    }
  ),

  # Generic effect forms -------------------------------------------------------
  effect_form = list(
    inline = function(x) c(.pf_shape(x$logical_shape),
      .pf_or(x$completeness, "unknown")),
    fields = function(x) list(
      effects = .pf_shape(x$logical_shape),
      measurements = .pf_or(x$total$dim[[1L]], "unknown"),
      completeness = .pf_or(x$completeness, "unknown"),
      capability = .pf_or(x$result_capability, "unspecified"),
      storage = .pf_set(x$storage, max = 3L),
      estimate = "signed cross-generalized form; PSD not assumed"
    )
  ),

  # The worked example bundle --------------------------------------------------
  effect_example_effects = list(
    fields = function(x) {
      truth <- x$truth
      list(
        fit = format(x$fit),
        domain = format(x$domain),
        frame = format(x$frame),
        contrast = paste0(length(x$contrast), " weights (",
          .pf_num(x$contrast), ")"),
        model_rdm = paste0(.pf_dim(x$model_rdm), " model dissimilarities"),
        planted = paste0(length(truth$planted_features), " pattern voxels + ",
          length(truth$mean_features), " mean-shift voxels"),
        truth = paste0(length(truth$signal_measurements),
          " signal measurements, seed ", truth$seed),
        `next` = "plan_geometry(example$fit$relation, example$frame, pairing)"
      )
    }
  ),

  # Measurement bridges, stores, and diagnostics -------------------------------
  effect_measurement_bridge = list(
    inline = function(x) c(
      paste0(x$common_space$n_measurements, " common measurements")),
    fields = function(x) list(
      "common space" = paste0(.pf_or(x$common_space$id, "unnamed"), ", ",
        .pf_or(x$common_space$n_measurements, "unknown"), " measurements"),
      left = paste0(.pf_or(x$left_domain$id, "unknown"), " via ",
        .pf_dim(x$left_leg)),
      right = paste0(.pf_or(x$right_domain$id, "unknown"), " via ",
        .pf_dim(x$right_leg)),
      signature = .pf_sig(x$signature)
    )
  ),
  effect_measurement_store = list(
    inline = function(x) c(paste0(nrow(x$index), " blocks"),
      x$representation),
    fields = function(x) list(
      blocks = nrow(x$index),
      representation = x$representation,
      codec = x$codec,
      elements = .pf_num(sum(x$index$length_elements)),
      path = if (is.null(x$path)) "in memory" else "on disk",
      signature = .pf_sig(x$signature),
      `next` = "effect_coupling(form) or rdm(form) to read blocks"
    )
  ),
  effect_measurement_receipt = list(
    inline = function(x) c(.pf_sig(x$route),
      .pf_or(x$execution$completion_status, "unknown")),
    fields = function(x) list(
      route = .pf_sig(x$route),
      status = .pf_or(x$execution$completion_status, "unknown"),
      "plan id" = .pf_sig(x$scientific_plan_id),
      derivation = .pf_set(names(x$derivation), max = 3L),
      signature = .pf_sig(x$signature)
    )
  ),
  effect_measurement_diagnostics = list(
    inline = function(x) c(x$method,
      paste0("effective rank ", x$experimental_effective_rank)),
    fields = function(x) list(
      method = x$method,
      "effective rank" = x$experimental_effective_rank,
      blocks = paste0(nrow(x$blocks), " (",
        .pf_set(names(x$blocks), max = 3L), ")"),
      regularization = .pf_or(x$regularization$kind, "none"),
      tolerance = .pf_num(x$tolerance),
      signature = .pf_sig(x$signature)
    )
  ),

  # Data-frame-backed indices --------------------------------------------------
  #
  # These subclass `data.frame`; only print() is defined so that
  # `format.data.frame` keeps working for internal callers.
  effect_measurement_block_index = list(
    fields = function(x) list(
      blocks = nrow(x),
      edges = .pf_set(x$edge_id, max = 3L),
      left = .pf_set(unique(x$left), max = 3L),
      right = .pf_set(unique(x$right), max = 3L),
      widths = paste0(.pf_set(as.character(unique(x$d_left)), max = 2L), " x ",
        .pf_set(as.character(unique(x$d_right)), max = 2L)),
      elements = .pf_num(sum(x$length_elements)),
      columns = .pf_set(names(x), max = 4L)
    )
  ),
  effect_ordered_edges = list(
    fields = function(x) list(
      products = nrow(x),
      left = .pf_set(unique(x$left), max = 3L),
      right = .pf_set(unique(x$right), max = 3L),
      orientation = .pf_set(unique(x$orientation), max = 3L),
      weights = if (length(unique(x$weight)) == 1L) {
        paste0("equal (", .pf_num(x$weight[[1L]]), ")")
      } else {
        paste0("unequal, ", .pf_num(range(x$weight)))
      },
      expansion = .pf_or(attr(x, "expansion", exact = TRUE), "none"),
      estimate = .pf_or(attr(x, "source_estimate", exact = TRUE), "unspecified")
    )
  ),

  # Task fragments -------------------------------------------------------------
  effect_evidence_stages = list(
    inline = function(x) c(paste0(length(x$order), " stages"),
      x$lowering),
    fields = function(x) list(
      order = .pf_set(x$order, max = 3L),
      lowering = x$lowering,
      normalization = .pf_stage(x$normalization),
      transform = .pf_stage(x$transform),
      reduction = .pf_stage(x$reduction),
      projection = .pf_stage(x$projection),
      signature = .pf_sig(x$signature)
    )
  ),
  effect_frame_conservation = list(
    inline = function(x) c(
      if (isTRUE(x$conserved)) "conserved" else "not conserved",
      x$component, paste0("max deviation ", .pf_num(x$max_deviation))),
    fields = function(x) {
      mass <- x$feature_mass
      list(
        conserved = .pf_yn(x$conserved),
        component = x$component,
        normalization = x$normalization,
        "max deviation" = .pf_num(x$max_deviation),
        tolerance = .pf_num(x$tolerance),
        "feature mass" = paste0(length(mass), " features, ",
          .pf_num(min(mass)), " to ", .pf_num(max(mass))),
        members = if (is.null(x$members)) {
          NULL
        } else {
          paste0(nrow(x$members), " family blocks, ",
            if (all(x$members$conserved)) "each" else "not each",
            " carrying its own alpha (max deviation ",
            .pf_num(max(x$members$max_deviation)), ")")
        }
      )
    }
  ),

  # Remaining query and metric records -----------------------------------------
  effect_query = list(
    inline = function(x) c(.pf_or(x$kind, "query"),
      .pf_dim(x$operator)),
    fields = function(x) list(
      kind = .pf_or(x$kind, "query"),
      operator = .pf_dim(x$operator),
      "effect space" = .pf_or(x$effect_space$basis_id, "unspecified"),
      fixed = .pf_yn(x$fixed)
    )
  ),
  effect_neural_pair_query = list(
    inline = function(x) c(.pf_dim(x$operator),
      x$left_domain$id),
    fields = function(x) {
      .validate_neural_pair_query(x)
      list(
        operator = .pf_dim(x$operator),
        left = paste0(x$left_domain$id, " (", x$left_domain$n_features,
          " features)"),
        right = paste0(x$right_domain$id, " (", x$right_domain$n_features,
          " features)"),
        signature = .pf_sig(x$signature)
      )
    }
  ),
  effect_pair_difference_query = list(
    inline = function(x) c(x$kind,
      paste0(length(x$pair_labels), " pairs")),
    fields = function(x) list(
      kind = x$kind,
      effects = paste0(length(x$effects), " (", .pf_set(x$effects, max = 4L),
        ")"),
      pairs = paste0(length(x$pair_labels), " (",
        .pf_set(x$pair_labels, max = 3L), ")"),
      coefficients = .pf_dim(x$coefficients)
    )
  ),
  effect_metric_inverse_representation = list(
    inline = function(x) c(x$kind),
    fields = function(x) list(
      kind = x$kind,
      signature = .pf_sig(x$signature),
      state = if (identical(x$kind, "none")) {
        "no inverse supplied; quadratic forms must be solved"
      } else {
        "inverse supplied; quadratic forms use it directly"
      }
    )
  ),
  effect_metric_handle = list(
    inline = function(x) c(.pf_or(x$inverse_mode, "unknown"),
      paste0(length(x$factorizations), " factorizations")),
    fields = function(x) list(
      role = .pf_or(x$metric$role, "unknown"),
      "inverse mode" = .pf_or(x$inverse_mode, "unknown"),
      shortcut = .pf_yn(x$shortcut),
      factorizations = length(x$factorizations),
      state = "resolved handle; the metric is fixed for this evaluation"
    )
  ),
  effect_metric_components = list(
    inline = function(x) c(
      paste0("coherent rank ", .pf_or(x$coherent_rank, "unknown")),
      .pf_or(x$inverse_quadratic_mode, "unknown")),
    fields = function(x) list(
      "coherent rank" = .pf_or(x$coherent_rank, "unknown"),
      "configuration psd" = .pf_yn(x$configuration_psd),
      denominator = .pf_or(x$denominator, "unspecified"),
      "inverse mode" = .pf_or(x$inverse_quadratic_mode, "unknown"),
      factorizations = .pf_or(x$factorization_count, "0"),
      signature = .pf_sig(x$signature)
    )
  ),
  effect_lowering = list(
    inline = function(x) c(x$kind,
      if (isTRUE(x$collapsed)) "collapsed" else "not collapsed"),
    fields = function(x) list(
      kind = x$kind,
      collapsed = .pf_yn(x$collapsed),
      reason = .pf_or(x$reason, "none recorded")
    )
  ),
  effect_residual_source = list(
    inline = function(x) c(.pf_shape(x$dim), "unread"),
    fields = function(x) list(
      dim = paste0(.pf_shape(x$dim), " (observations x features)"),
      revision = .pf_sig(x$stable_revision),
      state = "unread; residual blocks are pulled by feature index"
    )
  ),
  effect_coupling_partition_policy = list(
    inline = function(x) c(.pf_or(x$placement),
      paste0(length(x$weights), " partition weights")),
    fields = function(x) {
      .validate_coupling_partition_policy(x)
      weights <- x$weights
      list(
        placement = .pf_or(x$placement),
        weights = paste0(length(weights), ", ",
          if (length(unique(weights)) == 1L) {
            paste0("equal (", .pf_num(weights[[1L]]), ")")
          } else {
            paste0("unequal, ", .pf_num(range(weights)))
          }),
        transform = .pf_or(x$transform$kind, "none"),
        signature = .pf_sig(x$signature)
      )
    }
  ),

  # Measurement views, decompositions, and tomography --------------------------
  effect_measurement_view = list(
    inline = function(x) c(.pf_or(x$view),
      paste0(nrow(x$edge_index), " edges")),
    fields = function(x) {
      .validate_measurement_view(x)
      list(
        view = .pf_or(x$view),
        edges = paste0(nrow(x$edge_index), ", ", x$edge_completeness),
        dimensions = .pf_shape(x$dimensions),
        values = paste0(length(x$values), " blocks"),
        capability = .pf_or(x$result_capability, "unspecified"),
        completeness = .pf_or(x$completeness, "unknown")
      )
    }
  ),
  effect_measurement_decomposition = list(
    inline = function(x) c(x$id,
      paste0(length(x$components), " components")),
    fields = function(x) {
      .validate_measurement_decomposition(x)
      list(
        id = x$id,
        components = paste0(length(x$components), " (",
          .pf_set(names(x$components), max = 3L), ")"),
        "output space" = .pf_or(x$output_space$id, "unnamed"),
        orientation = .pf_or(x$orientation, "unspecified"),
        signature = .pf_sig(x$signature)
      )
    }
  ),
  effect_tomography_result = list(
    inline = function(x) c(x$method, x$status,
      if (isTRUE(x$lossless)) "lossless" else "lossy"),
    fields = function(x) {
      .validate_tomography_result(x)
      list(
        method = x$method,
        status = x$status,
        lossless = .pf_yn(x$lossless),
        certified = .pf_yn(x$certified),
        projections = paste0("left ", .pf_dim(x$left_projection), ", right ",
          .pf_dim(x$right_projection)),
        reference = .pf_sig(x$reference_signature),
        signature = .pf_sig(x$signature)
      )
    }
  ),

  # Measurement plan fragments -------------------------------------------------
  effect_measurement_plan = list(
    inline = function(x) c(paste0(nrow(x$block_index), " blocks"),
      x$query_construction, x$edge_scope),
    fields = function(x) list(
      blocks = nrow(x$block_index),
      relations = if (isTRUE(x$same_relation)) "one" else
        "two (cross-relation)",
      query = paste0(x$query_role, ", ", x$query_construction),
      "sampling axis" = .pf_or(x$sampling_axis, "not declared"),
      edges = paste0(nrow(x$edges), " (", x$edge_scope, ")"),
      "partition products" = paste0(nrow(x$partition_products), ", ",
        x$partition_expansion),
      regularization = .pf_or(x$regularization$kind, "none"),
      tolerance = .pf_num(x$tolerance),
      "plan id" = .pf_sig(x$scientific_plan_id),
      state = "query-first; no measurement block has been computed"
    )
  ),
  effect_measurement_leg = list(
    inline = function(x) c(x$kind, .pf_dim(x$operator)),
    fields = function(x) list(
      kind = x$kind,
      operator = .pf_dim(x$operator),
      domain = .pf_or(x$source_domain$id, "unknown"),
      "output space" = .pf_or(x$output_space$id, "unnamed"),
      support = paste0(length(x$support), " features"),
      estimation = .pf_or(x$estimation, "fixed"),
      signature = .pf_sig(x$signature)
    )
  ),
  effect_measurement_axis = list(
    inline = function(x) c(x$id,
      paste0(length(x$coordinates), " coordinates")),
    fields = function(x) {
      .validate_measurement_axis(x)
      units <- unique(as.character(x$units))
      list(
        id = x$id,
        coordinates = paste0(length(x$coordinates), " (",
          .pf_set(x$coordinates, max = 4L), ")"),
        basis = .pf_or(x$basis_id, "unspecified"),
        units = if (length(units) == 1L) units else .pf_set(units, max = 3L),
        signature = .pf_sig(x$signature)
      )
    }
  ),
  effect_measurement_edges = list(
    inline = function(x) c(paste0(nrow(x$edges), " edges"),
      if (isTRUE(x$weighted)) "weighted" else "unweighted"),
    fields = function(x) list(
      edges = nrow(x$edges),
      from = .pf_set(unique(x$edges[[1L]]), max = 3L),
      to = .pf_set(unique(x$edges[[2L]]), max = 3L),
      weighted = .pf_yn(x$weighted),
      signature = .pf_sig(x$signature)
    )
  ),
  effect_measurement_regularization = list(
    inline = function(x) c(x$kind,
      if (isTRUE(x$applied)) "applied" else "not applied"),
    fields = function(x) {
      .validate_measurement_regularization(x)
      list(
        kind = x$kind,
        applied = .pf_yn(x$applied),
        "lambda left" = .pf_num(x$lambda_left, empty = "none"),
        "lambda right" = .pf_num(x$lambda_right, empty = "none"),
        signature = .pf_sig(x$signature)
      )
    }
  ),
  effect_observation_whitener = list(
    inline = function(x) c(.pf_or(x$kind, "identity"),
      .pf_shape(x$dim)),
    fields = function(x) list(
      kind = .pf_or(x$kind, "identity"),
      dim = .pf_shape(x$dim),
      signature = .pf_sig(x$signature),
      state = "treated as fixed and known when the relation was fitted"
    )
  ),

  # Matched and control pair couplings -----------------------------------------
  effect_pair_coupling = list(
    inline = function(x) c(x$kind,
      paste0(sum(x$eligible), " eligible pairs")),
    fields = function(x) {
      .validate_pair_coupling(x)
      list(
        kind = x$kind,
        shape = .pf_dim(x$value),
        left = paste0(length(x$left_space$coordinates), " coordinates, ",
          .pf_or(x$left_space$basis_id, "unspecified")),
        right = paste0(length(x$right_space$coordinates), " coordinates, ",
          .pf_or(x$right_space$basis_id, "unspecified")),
        eligible = paste0(sum(x$eligible), " of ", length(x$eligible),
          " pairs"),
        mass = paste0("total ", .pf_num(sum(x$value)), ", max ",
          .pf_num(max(x$value))),
        `next` = if (identical(x$kind, "match")) {
          "control_coupling(matches) then coupling_contrast(matches, controls)"
        } else {
          "coupling_contrast(matches, controls)"
        }
      )
    }
  ),
  effect_location_transport = list(
    inline = function(x) c(.pf_transport_shape(x),
      x$semantics, .pf_or(x$provenance$method, "undeclared")),
    fields = function(x) {
      .validate_location_transport(x)
      sink <- .transport_sink_territory(x)
      declared <- !all(x$row_mass == 1)
      list(
        nodes = .pf_transport_shape(x),
        semantics = x$semantics,
        `row mass` = if (identical(x$semantics, "density") || declared) {
          if (declared) {
            paste0("declared, total ", .pf_num(sum(x$row_mass)))
          } else {
            "unit (one per native node)"
          }
        },
        sink = paste0("mass ", .pf_num(sink$mass), " of ", nrow(x$matrix),
          " rows, ", sprintf("%.1f%%", 100 * sink$share), " of territory"),
        provenance = .pf_transport_provenance(x$provenance),
        built = .pf_or(x$provenance$details, "undeclared"),
        signature = .pf_sig(x$signature)
      )
    }
  ),
  effect_population_plan = list(
    inline = function(x) c(
      .msg_count(length(x$subjects), "subject"),
      paste0(nrow(x$group_index), " group nodes"),
      x$model$formula_text, x$normalization),
    fields = function(x) {
      .validate_population_plan(x, deep = FALSE)
      list(
        subjects = .pf_population_subjects(x),
        `group nodes` = paste0(nrow(x$group_index), " + sink"),
        sink = .pf_population_sink(x),
        transport = paste0(x$semantics, ", ",
          .pf_set(unique(x$subject_index$provenance), max = 3L)),
        model = .pf_population_model(x),
        normalization = x$normalization,
        fit = .pf_population_fit(x),
        frames = if (isTRUE(x$allow_nonconservative)) {
          "non-conservative admitted (declared)"
        },
        estimand = .pf_sig(x$scientific_plan_id),
        signature = .pf_sig(x$signature)
      )
    }
  ),
  effect_population_result = list(
    inline = function(x) c(x$ledger,
      .msg_count(nrow(x$receipt$subjects), "subject"),
      paste0(nrow(x$index) - 1L, " group nodes + sink"),
      if (identical(x$basis, "complete_form")) {
        paste0("complete ", length(x$effects), "x", length(x$effects), " form")
      } else {
        .msg_count(nrow(x$queries), "query", "queries")
      })
  ),
  effect_population_uncertainty = list(
    inline = function(x) c(x$ledger,
      paste0(nrow(x$index) - 1L, " group nodes + sink"),
      .msg_count(length(x$term), "term"), "uncalibrated")
  ),
  effect_population_view = list(
    inline = function(x) c(.pf_population_view_kind(x), x$ledger,
      paste0("term ", x$term), .pf_population_view_rows(x),
      .msg_count(ncol(x$values), "column"))
  )
)

# Registering the records ----------------------------------------------------
#
# One printer and one formatter, each closed over the class it is registered
# for, stand in for every entry in `.pf_records`. Adding a sealed record costs
# a descriptor rather than two methods and two hand-maintained NAMESPACE
# lines. An entry with no `inline` gets no `format()` method, which is how
# `effect_pairing` keeps `format.data.frame` for its internal callers.

.pf_record_printer <- function(class, descriptor) {
  force(class)
  if (identical(descriptor$emit, "capabilities")) {
    return(function(x, ...) {
      .pf_capabilities_emit(class, x)
      invisible(x)
    })
  }
  fields <- descriptor$fields
  function(x, ...) {
    .pf_emit(class, fields(x))
    invisible(x)
  }
}

.pf_record_formatter <- function(class, descriptor) {
  force(class)
  if (identical(descriptor$inline, "capabilities")) {
    return(function(x, ...) .pf_capabilities_inline(class, x))
  }
  inline <- descriptor$inline
  function(x, ...) .pf_inline(class, inline(x))
}

.onLoad <- function(libname, pkgname) {
  for (class in names(.pf_records)) {
    descriptor <- .pf_records[[class]]
    if (!is.null(descriptor$fields) ||
        identical(descriptor$emit, "capabilities")) {
      .S3method("print", class, .pf_record_printer(class, descriptor))
    }
    if (!is.null(descriptor$inline)) {
      .S3method("format", class, .pf_record_formatter(class, descriptor))
    }
  }
  invisible(NULL)
}
