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
#                   R/evidence-*.R, ...
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

# Contiguous float64 run I/O, shared by the package's two block stores.
#
# `R/storage.R` (`effect_geometry_store`) and `R/measurement-storage.R`
# (`effect_measurement_store`) do *not* share an on-disk format and are not
# interchangeable. A geometry store is one dense column-major
# `dim[1] x dim[2]` matrix addressed by a column stride, with no sidecar. A
# measurement store is a concatenation of variable-shaped blocks addressed by
# an offset table, and it requires a `.manifest.rds` sidecar that records
# which blocks have been written and refuses a second write to any of them.
# Neither layout can be read with the other's offset arithmetic, so unifying
# them would change bytes on disk and invalidate recorded store signatures.
#
# What the two genuinely share is the transfer primitive underneath both:
# a headerless little-endian float64 payload, opened once, seeked to an
# element offset, transferred as one run, and closed once. That is this
# function, and it is the only place in the package that names the element
# size, the endianness, and the byte arithmetic `offset * 8`.
#
# `offsets` is in elements, not bytes. On read, `n` is the run length (one
# value, or one per offset) and the result is a list of double vectors; on
# write, `values` supplies one double vector per offset.
.tile_io <- function(path, offsets, mode = c("read", "write"), n = NULL,
                     values = NULL) {
  mode <- match.arg(mode)
  reading <- identical(mode, "read")
  connection <- file(path, open = if (reading) "rb" else "r+b")
  on.exit(close(connection), add = TRUE)
  if (reading) {
    n <- rep_len(n, length(offsets))
    lapply(seq_along(offsets), function(run) {
      seek(connection, where = offsets[[run]] * 8, origin = "start",
        rw = "read")
      readBin(connection, "double", n = n[[run]], size = 8,
        endian = "little")
    })
  } else {
    for (run in seq_along(offsets)) {
      seek(connection, where = offsets[[run]] * 8, origin = "start",
        rw = "write")
      writeBin(as.double(values[[run]]), connection, size = 8,
        endian = "little")
    }
    invisible(NULL)
  }
}

# The storage-format tag a store records in its manifest for a given effect-
# form codec. Both stores stamp it -- the in-memory one in R/result.R and the
# file-backed one in R/storage.R -- and a manifest is compared against it on
# read, so it is one constant in one place rather than a string a result file
# owns and the executor's storage layer borrows.
.effect_form_codec_format <- function(codec) {
  switch(codec,
    symmetric_packed = "packed-double-v1",
    rectangular = "rectangular-double-v1",
    .input_error("Unknown effect-form storage codec.")
  )
}

# Allocate a store payload: `n_values` float64 zeros, written in bounded
# chunks so that creating a large store never materializes it in memory.
.tile_zero_fill <- function(path, n_values) {
  connection <- file(path, open = "w+b")
  on.exit(if (!is.null(connection)) try(close(connection), silent = TRUE),
    add = TRUE)
  remaining <- n_values
  zero <- numeric(min(8192, remaining))
  while (remaining > 0) {
    count <- min(remaining, length(zero))
    writeBin(zero[seq_len(count)], connection, size = 8, endian = "little")
    remaining <- remaining - count
  }
  close(connection)
  connection <- NULL
  invisible(NULL)
}

# Workspace-plan record ------------------------------------------------------
#
# The byte accounting behind an `effect_memory_plan`: the category vocabulary,
# the validation of every scalar the plan carries, the arithmetic that turns
# categories into derived totals, and the record itself.
#
# It sits at the bottom of the package rather than beside `memory_plan()` in
# R/memory-plan.R because two files must be able to build this record. One is
# `memory_plan()`, the named-argument constructor. The other is
# `.validate_memory_plan_for_receipt()` in R/receipt.R, which rebuilds a plan
# from its own categories in order to prove that its derived fields were not
# edited after the fact. A receipt records canonical values; it must not be a
# dependency of the vocabulary that produces them, so it cannot call the
# constructor and pull itself into the layer-2 value cycle. Shared arithmetic
# in the primitive layer is reachable from both without joining them.

.workspace_plan_category_names <- function() {
  c("frame", "resident_source", "source_handles", "source_block",
    "relation_block", "atom_block", "local_state", "output", "contraction",
    "replacement_copy", "serialization_overlap", "reorder_buffer",
    "checkpoint_buffer")
}

.workspace_plan_record <- function(categories, workers, n_active,
                                   safety_factor, budget_bytes = NULL,
                                   measured_workspace_bytes = NULL,
                                   baseline_rss_bytes = NULL,
                                   peak_rss_bytes = NULL) {
  if (!identical(names(categories), .workspace_plan_category_names())) {
    .invariant_error("Workspace categories are missing or noncanonical.")
  }
  max_exact <- 2^53
  whole <- function(value, name, positive = FALSE) {
    lower <- if (positive) 1 else 0
    if (!.is_number(value) || value < lower || value %% 1 != 0 ||
        value > max_exact) {
      .input_error(sprintf(
        "`%s` must be one %sfinite whole scalar no greater than 2^53.",
        name, if (positive) "positive " else "nonnegative "))
    }
    value
  }
  for (name in names(categories)) whole(categories[[name]], name)
  workers <- whole(workers, "workers", positive = TRUE)
  n_active <- whole(n_active, "n_active", positive = TRUE)
  if (n_active > workers) .input_error("`n_active` cannot exceed `workers`.")
  if (!.is_number(safety_factor) || safety_factor < 1) {
    .input_error(
      "`safety_factor` must be one finite number greater than or equal to one."
    )
  }
  optional <- function(value, name) {
    if (!is.null(value) && (!is.numeric(value) || length(value) != 1L ||
        is.na(value) || !is.finite(value) || value < 0)) {
      .input_error(sprintf(
        "`%s` must be NULL or one nonnegative finite byte count.",
        name))
    }
    value
  }
  optional(budget_bytes, "budget_bytes")
  if (!is.null(budget_bytes) && budget_bytes == 0) {
    .input_error("`budget_bytes` must be NULL or positive.")
  }
  optional(measured_workspace_bytes, "measured_workspace_bytes")
  optional(baseline_rss_bytes, "baseline_rss_bytes")
  optional(peak_rss_bytes, "peak_rss_bytes")
  if (!is.null(baseline_rss_bytes) && !is.null(peak_rss_bytes) &&
      peak_rss_bytes < baseline_rss_bytes) {
    .input_error("Peak RSS cannot be smaller than baseline RSS.")
  }

  persistent_names <- c("frame", "resident_source", "source_handles",
    "local_state", "output")
  task_names <- c("source_block", "relation_block", "atom_block",
    "contraction", "replacement_copy", "serialization_overlap",
    "reorder_buffer", "checkpoint_buffer")
  persistent <- sum(categories[persistent_names])
  per_active <- sum(categories[task_names])
  active <- n_active * per_active
  if (!is.finite(active) || active > max_exact ||
      !is.finite(persistent + active) || persistent + active > max_exact) {
    .invariant_error(
      "Workspace byte accounting overflows exact representation."
    )
  }
  modeled <- persistent + active
  conservative <- ceiling(modeled * safety_factor)
  incremental_rss <- if (is.null(baseline_rss_bytes) || is.null(peak_rss_bytes)) {
    NULL
  } else {
    peak_rss_bytes - baseline_rss_bytes
  }

  structure(list(
    categories = categories,
    workers = as.integer(workers),
    n_active = as.integer(n_active),
    persistent_workspace_bytes = persistent,
    task_workspace_per_active_bytes = per_active,
    active_task_workspace_bytes = active,
    modeled_workspace_bytes = modeled,
    safety_factor = safety_factor,
    planned_workspace_bytes = conservative,
    budget_bytes = budget_bytes,
    fits_budget = if (is.null(budget_bytes)) NA else conservative <= budget_bytes,
    measured_workspace_bytes = measured_workspace_bytes,
    measured_workspace_within_plan = if (is.null(measured_workspace_bytes)) NA else
      measured_workspace_bytes <= conservative,
    baseline_rss_bytes = baseline_rss_bytes,
    incremental_peak_rss_bytes = incremental_rss,
    absolute_peak_rss_bytes = peak_rss_bytes,
    prediction_kind = "crossform_owned_workspace_upper_bound"
  ), class = "effect_memory_plan")
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
