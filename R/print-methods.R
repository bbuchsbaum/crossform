# Compact printing for user-facing objects -----------------------------------
#
# Every class a user can hold gets a one-screen summary: a `<class>` header,
# a short block of aligned `key: value` lines, and -- where the object is lazy
# or has an obvious successor -- a closing `next:` or `state:` hint naming the
# exact call to make. Printing never materializes, never dereferences a source,
# and never emits closures, environments, whole matrices, or whole hashes.

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
  paste(format(signif(as.numeric(x), digits), trim = TRUE), collapse = ", ")
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

#' @export
format.effect_source_capabilities <- function(x, ...) {
  .pf_capabilities_inline("effect_source_capabilities", x)
}

#' @export
print.effect_source_capabilities <- function(x, ...) {
  .pf_emit("effect_source_capabilities", list(
    block_read = .pf_yn(x$block_read),
    reopenable = .pf_yn(x$reopenable),
    thread_safe = .pf_yn(x$thread_safe),
    revision = .pf_sig(x$stable_revision)
  ))
  invisible(x)
}

#' @export
format.effect_error_capabilities <- function(x, ...) {
  .pf_capabilities_inline("effect_error_capabilities", x)
}

#' @export
print.effect_error_capabilities <- function(x, ...) {
  .pf_capabilities_emit("effect_error_capabilities", x)
  invisible(x)
}

#' @export
format.effect_metric_capabilities <- function(x, ...) {
  .pf_capabilities_inline("effect_metric_capabilities", x)
}

#' @export
print.effect_metric_capabilities <- function(x, ...) {
  granted <- .pf_capability_flags(x, drop = "role")
  .pf_emit("effect_metric_capabilities", list(
    role = .pf_or(x$role),
    granted = .pf_set(names(granted)[granted == "yes"], max = 4L,
      empty = "none"),
    withheld = .pf_set(names(granted)[granted == "no"], max = 4L,
      empty = "none")
  ))
  invisible(x)
}

#' @export
format.effect_measurement_capabilities <- function(x, ...) {
  .pf_capabilities_inline("effect_measurement_capabilities", x)
}

#' @export
print.effect_measurement_capabilities <- function(x, ...) {
  .pf_capabilities_emit("effect_measurement_capabilities", x)
  invisible(x)
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

# Neural domains -------------------------------------------------------------

#' @export
format.effect_domain <- function(x, ...) {
  .pf_inline("effect_domain", x$id,
    paste0(x$n_features, " ", x$kind, " features"))
}

#' @export
print.effect_domain <- function(x, ...) {
  .validate_domain(x)
  coordinates <- if (is.null(x$coordinates)) {
    "none (abstract)"
  } else {
    paste0(.pf_dim(x$coordinates), " (",
      .pf_set(x$coordinate_units, max = 3L), ")")
  }
  .pf_emit("effect_domain", list(
    id = x$id,
    kind = x$kind,
    features = x$n_features,
    "feature ids" = .pf_set(x$feature_ids, max = 4L),
    coordinates = coordinates,
    geometry = .pf_sig(x$geometry_signature),
    `next` = "compile_frame(whole_brain(), domain) to place a frame on it"
  ))
  invisible(x)
}

#' @export
format.effect_domain_reference <- function(x, ...) {
  .pf_inline("effect_domain_reference", x$id,
    paste0(x$n_features, " features"))
}

#' @export
print.effect_domain_reference <- function(x, ...) {
  .validate_domain_reference(x)
  .pf_emit("effect_domain_reference", list(
    id = x$id,
    features = x$n_features,
    units = .pf_set(x$coordinate_units, max = 3L),
    geometry = .pf_sig(x$geometry_signature),
    signature = .pf_sig(x$signature),
    state = "identity only; the referenced domain carries the coordinates"
  ))
  invisible(x)
}

# Frames ---------------------------------------------------------------------

#' @export
format.effect_frame_spec <- function(x, ...) {
  detail <- if (is.null(x$radius)) {
    NULL
  } else {
    paste0("radius ", .pf_num(x$radius))
  }
  .pf_inline("effect_frame_spec", x$kind, detail,
    paste0(x$normalization, " normalization"))
}

#' @export
print.effect_frame_spec <- function(x, ...) {
  .pf_emit("effect_frame_spec", list(
    kind = x$kind,
    radius = if (is.null(x$radius)) NULL else .pf_num(x$radius),
    normalization = x$normalization,
    state = "unplaced; call compile_frame(spec, domain) to bind a domain"
  ))
  invisible(x)
}

#' @export
format.effect_frame <- function(x, ...) {
  .pf_inline("effect_frame", paste0(nrow(x$index), " nodes"),
    x$representation, x$specification$kind)
}

#' @export
print.effect_frame <- function(x, ...) {
  weights <- x$weights
  specification <- x$specification
  placement <- paste0(specification$kind,
    if (is.null(specification$radius)) {
      ""
    } else {
      paste0(", radius ", .pf_num(specification$radius), " ",
        .pf_or(x$domain$coordinate_units, ""))
    })
  .pf_emit("effect_frame", list(
    representation = x$representation,
    nodes = nrow(x$index),
    features = x$domain$n_features,
    specification = placement,
    normalization = x$normalization,
    weights = paste0(.pf_dim(weights), ", ", .pf_stored(weights), " stored"),
    domain = paste0(x$domain_id, " (", x$domain_kind, ")"),
    fixed = paste0(.pf_yn(x$fixed), "; locally estimated ",
      .pf_yn(x$locally_estimated)),
    `next` = "plan_geometry(relation, frame, pairing)"
  ))
  invisible(x)
}

#' @export
format.effect_support_index <- function(x, ...) {
  .pf_inline("effect_support_index", paste0(length(x$node_ids), " nodes"),
    paste0(length(x$members), " memberships"), x$construction$kind)
}

#' @export
print.effect_support_index <- function(x, ...) {
  .validate_support_index(x, deep = FALSE)
  cost <- x$cost
  sizes <- cost$support_size
  construction <- x$construction
  .pf_emit("effect_support_index", list(
    nodes = length(x$node_ids),
    construction = paste0(construction$kind, " via ", construction$provider),
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
  ))
  invisible(x)
}

#' @export
format.effect_support_cost <- function(x, ...) {
  .pf_inline("effect_support_cost", paste0(x$nodes, " nodes"),
    paste0(x$support_memberships, " memberships"))
}

#' @export
print.effect_support_cost <- function(x, ...) {
  .pf_emit("effect_support_cost", list(
    nodes = .pf_num(x$nodes),
    features = .pf_num(x$features),
    memberships = .pf_num(x$support_memberships),
    "pair pattern" = if (is.na(x$pair_pattern_stored_nnz)) {
      "not materialized"
    } else {
      paste0(.pf_num(x$pair_pattern_stored_nnz),
        " stored of ", .pf_num(x$pair_pattern_nnz))
    },
    "diagonal metric" = paste0(.pf_num(x$diagonal_metric_entries), " entries"),
    "dense metric" = paste0(.pf_num(x$dense_metric_entries), " entries, ",
      .pf_bytes(x$materialized_dense_metric_bytes)),
    "largest local" = .pf_bytes(x$largest_local_metric_bytes),
    structural = .pf_bytes(x$estimated_structural_bytes),
    state = "cost projection only; nothing is materialized"
  ))
  invisible(x)
}

# Relations and their parts --------------------------------------------------

#' @export
format.effect_relation <- function(x, ...) {
  .pf_inline("effect_relation", paste0(length(x$effects), " effects"),
    paste0(length(x$partitions), " partitions"),
    paste0(x$n_features, " features"))
}

#' @export
print.effect_relation <- function(x, ...) {
  .validate_relation(x, deep = FALSE)
  kinds <- vapply(x$sources, function(source) .pf_or(source$kind, "unknown"),
    character(1))
  estimated <- is.null(x$extractors) ||
    all(vapply(x$extractors, is.null, logical(1)))
  .pf_emit("effect_relation", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_space <- function(x, ...) {
  .pf_inline("effect_space", paste0(length(x$coordinates), " coordinates"),
    x$basis_id)
}

#' @export
print.effect_space <- function(x, ...) {
  .validate_effect_space(x)
  unit <- unique(as.character(x$units))
  scale <- unique(as.numeric(x$scale))
  .pf_emit("effect_space", list(
    coordinates = paste0(length(x$coordinates), " (",
      .pf_set(x$coordinates, max = 4L), ")"),
    basis = .pf_or(x$basis_id, "unspecified"),
    units = if (length(unit) == 1L) unit else .pf_set(unit, max = 3L),
    scale = if (length(scale) == 1L) .pf_num(scale) else .pf_num(range(scale)),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

#' @export
format.effect_extractor <- function(x, ...) {
  .pf_inline("effect_extractor", paste0(length(x$effects), " effects"),
    paste0("from ", x$n_observations, " observations"), x$estimator)
}

#' @export
print.effect_extractor <- function(x, ...) {
  .validate_effect_extractor(x)
  diagnostics <- x$diagnostics
  .pf_emit("effect_extractor", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_response_source <- function(x, ...) {
  .pf_inline("effect_response_source", x$kind, .pf_shape(x$dim))
}

#' @export
print.effect_response_source <- function(x, ...) {
  .pf_emit("effect_response_source", list(
    kind = x$kind,
    dim = paste0(.pf_shape(x$dim), " (observations x features)"),
    access = .pf_or(x$descriptor$access, "coordinator"),
    revision = .pf_sig(x$descriptor$stable_revision),
    state = "unread; blocks are pulled on demand by feature index"
  ))
  invisible(x)
}

#' @export
format.effect_source_descriptor <- function(x, ...) {
  .pf_inline("effect_source_descriptor", x$kind, .pf_shape(x$dim))
}

#' @export
print.effect_source_descriptor <- function(x, ...) {
  .pf_emit("effect_source_descriptor", list(
    kind = x$kind,
    dim = .pf_shape(x$dim),
    access = .pf_or(x$access, "coordinator"),
    revision = .pf_sig(x$stable_revision),
    spec = .pf_set(names(x$spec), max = 4L, empty = "none")
  ))
  invisible(x)
}

#' @export
format.effect_error_model <- function(x, ...) {
  .pf_inline("effect_error_model", x$kind,
    paste0(x$residual_df, " residual df"), x$sampling_unit)
}

#' @export
print.effect_error_model <- function(x, ...) {
  .validate_error_model(x, deep = FALSE)
  .pf_emit("effect_error_model", list(
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
  ))
  invisible(x)
}

# Pairings -------------------------------------------------------------------
#
# `effect_pairing` subclasses `data.frame`, so only print() is defined here:
# a format() method would shadow `format.data.frame` for internal callers.

#' @export
print.effect_pairing <- function(x, ...) {
  .pf_emit("effect_pairing", list(
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
    independence = .pf_or(attr(x, "independence", exact = TRUE), "undeclared"),
    generalizes = .pf_or(attr(x, "generalizes_over", exact = TRUE),
      "axis undeclared"),
    estimate = .pf_or(attr(x, "estimate", exact = TRUE), "unspecified")
  ))
  invisible(x)
}

# Metrics --------------------------------------------------------------------

#' @export
format.effect_neural_metric <- function(x, ...) {
  .pf_inline("effect_neural_metric", x$role,
    paste0(length(x$support), " supported features"), x$estimation)
}

#' @export
print.effect_neural_metric <- function(x, ...) {
  .validate_neural_metric(x, deep = FALSE)
  capabilities <- x$capabilities
  .pf_emit("effect_neural_metric", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_metric_recipe <- function(x, ...) {
  .pf_inline("effect_metric_recipe", x$kind, x$role, "unfitted")
}

#' @export
print.effect_metric_recipe <- function(x, ...) {
  .validate_metric_recipe(x)
  .pf_emit("effect_metric_recipe", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_metric_schedule <- function(x, ...) {
  .pf_inline("effect_metric_schedule", x$kind, x$scope, x$materialization)
}

#' @export
print.effect_metric_schedule <- function(x, ...) {
  .pf_emit("effect_metric_schedule", list(
    role = x$role,
    kind = x$kind,
    scope = x$scope,
    composition = x$frame_composition,
    lowering = x$lowering,
    materialization = x$materialization,
    properties = paste0("feature additive ", .pf_yn(x$feature_additive),
      ", dense support ", .pf_yn(x$support_dense)),
    metric = if (is.null(x$metric)) {
      "implicit identity"
    } else {
      .pf_sig(x$metric_signature)
    },
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

#' @export
format.effect_frozen_metric_schedule <- function(x, ...) {
  .pf_inline("effect_frozen_metric_schedule", x$recipe$kind,
    x$training_policy$kind, "frozen before evaluation")
}

#' @export
print.effect_frozen_metric_schedule <- function(x, ...) {
  .validate_frozen_metric_schedule(x, deep = FALSE)
  .pf_emit("effect_frozen_metric_schedule", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_metric_training_policy <- function(x, ...) {
  .pf_inline("effect_metric_training_policy", x$kind)
}

#' @export
print.effect_metric_training_policy <- function(x, ...) {
  .validate_metric_training_policy(x)
  .pf_emit("effect_metric_training_policy", list(
    kind = x$kind,
    "uses evaluation residuals" =
      .pf_yn(x$includes_evaluation_residuals),
    assumption = .pf_or(x$assumption, "none"),
    justification = .pf_or(x$justification, "none recorded"),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

# Policies and contracts -----------------------------------------------------

#' @export
format.effect_compute_policy <- function(x, ...) {
  .pf_inline("effect_compute_policy", paste0(x$workers, " workers"),
    x$process_backend)
}

#' @export
print.effect_compute_policy <- function(x, ...) {
  .pf_emit("effect_compute_policy", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_memory_plan <- function(x, ...) {
  .pf_inline("effect_memory_plan",
    paste0(.pf_bytes(x$planned_workspace_bytes), " planned"),
    paste0(x$workers, " workers"))
}

#' @export
print.effect_memory_plan <- function(x, ...) {
  categories <- x$categories
  ranked <- sort(categories[categories > 0], decreasing = TRUE)
  .pf_emit("effect_memory_plan", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_numerical_contract <- function(x, ...) {
  .pf_inline("effect_numerical_contract", paste0("atol ", .pf_num(x$atol)),
    paste0("rtol ", .pf_num(x$rtol)))
}

#' @export
print.effect_numerical_contract <- function(x, ...) {
  .validate_numerical_contract(x)
  .pf_emit("effect_numerical_contract", list(
    atol = .pf_num(x$atol),
    rtol = .pf_num(x$rtol),
    scheduling = .pf_or(x$scheduling$guarantee, "unspecified"),
    "block partition" = .pf_or(x$block_partition$guarantee, "unspecified"),
    "cross platform" = .pf_or(x$cross_platform$guarantee, "unspecified"),
    bitwise = paste0("across blocking ",
      .pf_yn(x$bitwise_across_blocking), ", across platforms ",
      .pf_yn(x$bitwise_across_platforms))
  ))
  invisible(x)
}

#' @export
format.effect_numeric_agreement <- function(x, ...) {
  .pf_inline("effect_numeric_agreement",
    if (isTRUE(x$passed)) "passed" else "failed", x$guarantee)
}

#' @export
print.effect_numeric_agreement <- function(x, ...) {
  .pf_emit("effect_numeric_agreement", list(
    result = if (isTRUE(x$passed)) "passed" else "failed",
    guarantee = x$guarantee,
    comparison = x$comparison,
    "max error" = .pf_num(x$max_absolute_error),
    "max allowed" = .pf_num(x$max_allowed_error),
    tolerance = paste0("atol ", .pf_num(x$atol), ", rtol ", .pf_num(x$rtol))
  ))
  invisible(x)
}

#' @export
format.effect_partition_reducer <- function(x, ...) {
  .pf_inline("effect_partition_reducer", x$kind, x$weight_convention)
}

#' @export
print.effect_partition_reducer <- function(x, ...) {
  .validate_partition_reducer(x)
  .pf_emit("effect_partition_reducer", list(
    kind = x$kind,
    weights = x$weight_convention,
    order = x$order
  ))
  invisible(x)
}

#' @export
format.effect_gaussian_covariance_model <- function(x, ...) {
  .pf_inline("effect_gaussian_covariance_model", x$family,
    if (isTRUE(x$fixed)) "fixed" else "estimated")
}

#' @export
print.effect_gaussian_covariance_model <- function(x, ...) {
  .validate_gaussian_covariance_model(x)
  .pf_emit("effect_gaussian_covariance_model", list(
    family = x$family,
    fixed = .pf_yn(x$fixed),
    signature = .pf_sig(x$signature),
    state = "joint Gaussian declared; information units are claimable"
  ))
  invisible(x)
}

# Compiled tasks, sampling plans, and receipts -------------------------------

#' @export
format.effect_evidence_task <- function(x, ...) {
  .pf_inline("effect_evidence_task",
    paste0(nrow(x$ordered_partition_products), " partition products"),
    if (isTRUE(x$same_relation)) "one relation" else "two relations",
    x$stages$lowering)
}

#' @export
print.effect_evidence_task <- function(x, ...) {
  spaces <- x$spaces
  .pf_emit("effect_evidence_task", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_sampling_record <- function(x, ...) {
  .pf_inline("effect_sampling_record", x$kind)
}

#' @export
print.effect_sampling_record <- function(x, ...) {
  flags <- .pf_capability_flags(x, drop = c("kind", "signature"))
  fields <- c(list(kind = x$kind), as.list(utils::head(flags, 8L)))
  fields[["signature"]] <- .pf_sig(x$signature)
  .pf_emit("effect_sampling_record", fields)
  invisible(x)
}

#' @export
format.effect_evidence_sampling_plan <- function(x, ...) {
  .pf_inline("effect_evidence_sampling_plan", x$contract, x$sampling_axis,
    x$target$target)
}

#' @export
print.effect_evidence_sampling_plan <- function(x, ...) {
  .validate_evidence_sampling_plan(x, deep = FALSE)
  unmet <- x$unavailable_reasons
  .pf_emit("effect_evidence_sampling_plan", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_sampling_covariance <- function(x, ...) {
  .pf_inline("effect_sampling_covariance",
    paste0(x$dimension, " coordinates"), x$plan$target$policy)
}

#' @export
print.effect_sampling_covariance <- function(x, ...) {
  .validate_sampling_covariance(x, deep = FALSE)
  .pf_emit("effect_sampling_covariance", list(
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
  ))
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

#' @export
format.effect_execution_receipt <- function(x, ...) {
  .pf_inline("effect_execution_receipt", x$completion_status,
    paste0(x$completed_task_count, "/", x$task_count, " tasks"))
}

#' @export
print.effect_execution_receipt <- function(x, ...) {
  .validate_execution_receipt(x)
  .pf_emit("effect_execution_receipt", list(
    status = paste0(x$completion_status, " (", x$completed_task_count, "/",
      x$task_count, " tasks)"),
    "plan id" = .pf_sig(x$scientific_plan_id),
    kernel = x$kernel_version,
    partitioning = x$task_partition_id,
    precision = x$precision,
    compute = paste0(x$compute$workers, " workers, ",
      .pf_or(x$compute$process_backend, "sequential")),
    memory = paste0(.pf_bytes(x$memory$planned_workspace_bytes), " planned"),
    tolerance = paste0("atol ", .pf_num(x$numeric_contract$atol), ", rtol ",
      .pf_num(x$numeric_contract$rtol)),
    blas = paste0(.pf_or(x$blas$requested_threads, "unrecorded"),
      " thread(s) requested"),
    domain = .pf_sig(x$domain_signature)
  ))
  invisible(x)
}

#' @export
format.effect_geometry_store <- function(x, ...) {
  .pf_inline("effect_geometry_store", .pf_shape(x$dim), x$representation)
}

#' @export
print.effect_geometry_store <- function(x, ...) {
  .pf_emit("effect_geometry_store", list(
    dim = paste0(.pf_shape(x$dim), " (measurements x packed values)"),
    representation = x$representation,
    format = .pf_or(x$manifest$format, "unspecified"),
    `next` = "geometry_component(geometry, \"total\")"
  ))
  invisible(x)
}

#' @export
format.effect_residual_pair_statistics <- function(x, ...) {
  .pf_inline("effect_residual_pair_statistics",
    paste0(length(x$pair_i), " feature pairs"),
    paste0(length(x$partitions), " partitions"))
}

#' @export
print.effect_residual_pair_statistics <- function(x, ...) {
  .validate_residual_pair_statistics(x, deep = FALSE)
  .pf_emit("effect_residual_pair_statistics", list(
    "feature pairs" = length(x$pair_i),
    partitions = paste0(length(x$partitions), " (",
      .pf_set(x$partitions, max = 3L), ")"),
    domain = x$domain$id,
    support = .pf_sig(x$support_index),
    "relation fit" = .pf_sig(x$relation_fit),
    signature = .pf_sig(x$signature),
    state = "sufficient statistics only; residuals are not retained"
  ))
  invisible(x)
}

# Study facts ----------------------------------------------------------------

#' @export
format.effect_observations <- function(x, ...) {
  .pf_inline("effect_observations", paste0(length(x$partitions), " partitions"),
    paste0(x$n_features, " features"))
}

#' @export
print.effect_observations <- function(x, ...) {
  .validate_observations(x, deep = FALSE)
  counts <- vapply(x$indexes, function(index) {
    length(index$observation_id)
  }, integer(1))
  .pf_emit("effect_observations", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_observation_index <- function(x, ...) {
  .pf_inline("effect_observation_index", x$partition,
    paste0(length(x$observation_id), " observations"))
}

#' @export
print.effect_observation_index <- function(x, ...) {
  .validate_observation_index(x)
  .pf_emit("effect_observation_index", list(
    partition = x$partition,
    observations = length(x$observation_id),
    ids = .pf_set(x$observation_id, max = 4L),
    time = if (isTRUE(x$timing)) {
      paste0(.pf_num(min(x$time)), " to ", .pf_num(max(x$time)), " ", x$units)
    } else {
      "none declared"
    },
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

#' @export
format.effect_events <- function(x, ...) {
  .pf_inline("effect_events", paste0(nrow(x$data), " events"),
    paste0(length(unique(x$data[[x$partition_column]])), " partitions"))
}

#' @export
print.effect_events <- function(x, ...) {
  .validate_events(x)
  columns <- setdiff(names(x$data), c(x$partition_column, x$event_id_column,
    x$onset_column, x$duration_column))
  .pf_emit("effect_events", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_observation_confounds <- function(x, ...) {
  .pf_inline("effect_observation_confounds", paste0(nrow(x$data), " rows"),
    paste0(length(x$schema) - 2L, " regressors"))
}

#' @export
print.effect_observation_confounds <- function(x, ...) {
  .validate_observation_confounds(x)
  regressors <- setdiff(names(x$data), c(x$partition_column,
    x$observation_id_column, x$censor_column))
  censored <- if (is.null(x$censor_column)) {
    "none declared"
  } else {
    paste0(sum(!x$data[[x$censor_column]]), " of ", nrow(x$data),
      " observations dropped by `", x$censor_column, "`")
  }
  .pf_emit("effect_observation_confounds", list(
    rows = nrow(x$data),
    partitions = .pf_set(unique(x$data[[x$partition_column]]), max = 3L),
    regressors = paste0(length(regressors), " (",
      .pf_set(regressors, max = 3L), ")"),
    censor = censored,
    id = .pf_sig(x$confounds_id)
  ))
  invisible(x)
}

#' @export
format.effect_partition_hierarchy <- function(x, ...) {
  .pf_inline("effect_partition_hierarchy", paste0(nrow(x$data), " partitions"),
    .pf_set(x$axes, max = 3L))
}

#' @export
print.effect_partition_hierarchy <- function(x, ...) {
  .validate_partition_hierarchy(x)
  widths <- vapply(x$levels, length, integer(1))
  .pf_emit("effect_partition_hierarchy", list(
    leaf = x$leaf,
    partitions = nrow(x$data),
    axes = paste0(length(x$axes), " (", .pf_set(x$axes, max = 4L), ")"),
    levels = .pf_set(paste0(names(widths), " ", widths), max = 4L),
    signature = .pf_sig(x$signature),
    `next` = "study_axis(study, \"<axis>\") to resample over an axis"
  ))
  invisible(x)
}

#' @export
format.effect_observation_model <- function(x, ...) {
  .pf_inline("effect_observation_model", .pf_or(x$kind),
    .pf_or(x$sampling_unit))
}

#' @export
print.effect_observation_model <- function(x, ...) {
  .pf_emit("effect_observation_model", list(
    kind = .pf_or(x$kind),
    "sampling unit" = .pf_or(x$sampling_unit),
    independence = .pf_or(x$independence, "undeclared"),
    whitener = .pf_or(x$whitener$kind, "identity"),
    assumptions = .pf_set(x$assumptions, max = 2L),
    training = .pf_sig(x$training_revision),
    id = .pf_sig(x$observation_model_id)
  ))
  invisible(x)
}

# Measurement spaces, frames, and edges --------------------------------------

#' @export
format.effect_measurement_space <- function(x, ...) {
  .pf_inline("effect_measurement_space", x$id,
    paste0(x$n_measurements, " measurements"))
}

#' @export
print.effect_measurement_space <- function(x, ...) {
  .validate_measurement_space(x)
  .pf_emit("effect_measurement_space", list(
    id = x$id,
    measurements = x$n_measurements,
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

#' @export
format.effect_measurement_frame <- function(x, ...) {
  .pf_inline("effect_measurement_frame", paste0(length(x$node_ids), " nodes"),
    paste0("from ", x$source_domain$id))
}

#' @export
print.effect_measurement_frame <- function(x, ...) {
  .validate_measurement_frame(x)
  widths <- vapply(x$legs, function(leg) nrow(leg$operator), integer(1))
  .pf_emit("effect_measurement_frame", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_edge_frame <- function(x, ...) {
  .pf_inline("effect_edge_frame", paste0(nrow(x$edges$edges), " edges"),
    paste0(length(x$from_frame$node_ids), " x ",
      length(x$to_frame$node_ids), " nodes"))
}

#' @export
print.effect_edge_frame <- function(x, ...) {
  .validate_edge_frame(x)
  .pf_emit("effect_edge_frame", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_pair_query <- function(x, ...) {
  .pf_inline("effect_pair_query", x$kind, .pf_dim(x$operator),
    if (isTRUE(x$fixed)) "fixed" else "estimated")
}

#' @export
print.effect_pair_query <- function(x, ...) {
  .pf_emit("effect_pair_query", list(
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
  ))
  invisible(x)
}

# Measurement forms and coupling views ---------------------------------------

#' @export
format.effect_measurement_form <- function(x, ...) {
  .pf_inline("effect_measurement_form", paste0(nrow(x$block_index), " blocks"),
    x$storage, x$completeness)
}

#' @export
print.effect_measurement_form <- function(x, ...) {
  .validate_measurement_form(x, probe = FALSE)
  capabilities <- x$capabilities
  .pf_emit("effect_measurement_form", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_coupling_result <- function(x, ...) {
  .pf_inline("effect_coupling_result", x$kind,
    paste0(nrow(x$edge_index), " edges"), x$terminology)
}

#' @export
print.effect_coupling_result <- function(x, ...) {
  .validate_coupling_result(x)
  .pf_emit("effect_coupling_result", list(
    kind = x$kind,
    terminology = x$terminology,
    edges = paste0(nrow(x$edge_index), ", ", x$edge_completeness),
    values = paste0(length(x$values), " blocks"),
    normalization = .pf_or(x$normalization_axis, "none"),
    summary = .pf_or(x$summary_axis, "none"),
    regularization = if (is.null(x$regularization)) {
      "none"
    } else {
      .pf_num(unlist(x$regularization, use.names = FALSE))
    },
    units = .pf_or(x$units, "not claimed"),
    stages = .pf_set(x$stage_order, max = 3L),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

# Design models and their receipts -------------------------------------------

#' @export
format.effect_design_model <- function(x, ...) {
  .pf_inline("effect_design_model", paste0(length(x$designs), " partitions"),
    paste0(length(x$condition_space$coordinates), " conditions"))
}

#' @export
print.effect_design_model <- function(x, ...) {
  .validate_design_model(x)
  rows <- vapply(x$designs, nrow, integer(1))
  columns <- vapply(x$designs, ncol, integer(1))
  granted <- .pf_capability_flags(x$capabilities)
  .pf_emit("effect_design_model", list(
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
  ))
  invisible(x)
}

#' @export
format.effect_raw_design_model <- function(x, ...) {
  .pf_inline("effect_raw_design_model",
    paste0(length(x$designs), " partitions"), "no symbolic model")
}

#' @export
print.effect_raw_design_model <- function(x, ...) {
  rows <- vapply(x$designs, nrow, integer(1))
  columns <- vapply(x$designs, ncol, integer(1))
  .pf_emit("effect_raw_design_model", list(
    partitions = paste0(length(x$designs), " (",
      .pf_set(names(x$designs), max = 3L), ")"),
    designs = paste0(.pf_set(as.character(rows), max = 2L), " rows x ",
      .pf_set(as.character(unique(columns)), max = 2L), " columns"),
    solver = .pf_set(unique(as.character(x$solver)), max = 2L),
    state = "raw columns; coding invariance and row lineage are not claimed"
  ))
  invisible(x)
}

#' @export
format.effect_design_receipt <- function(x, ...) {
  .pf_inline("effect_design_receipt", x$partition,
    paste0("rank ", x$rank), paste0(x$residual_df, " residual df"))
}

#' @export
print.effect_design_receipt <- function(x, ...) {
  .validate_design_receipt(x)
  censoring <- x$censoring
  .pf_emit("effect_design_receipt", list(
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
  ))
  invisible(x)
}

# Generic effect forms -------------------------------------------------------

#' @export
format.effect_form <- function(x, ...) {
  .pf_inline("effect_form", .pf_shape(x$logical_shape),
    .pf_or(x$completeness, "unknown"))
}

#' @export
print.effect_form <- function(x, ...) {
  .pf_emit("effect_form", list(
    effects = .pf_shape(x$logical_shape),
    measurements = .pf_or(x$total$dim[[1L]], "unknown"),
    completeness = .pf_or(x$completeness, "unknown"),
    capability = .pf_or(x$result_capability, "unspecified"),
    storage = .pf_set(x$storage, max = 3L),
    estimate = "signed cross-generalized form; PSD not assumed"
  ))
  invisible(x)
}

# Query-first geometry plans -------------------------------------------------

# The plan is the object a user is told to keep with an analysis record, so the
# default block says what was declared, in the words the declaration used. The
# compiler's own vocabulary -- how the metric is lowered, how large a full
# materialization would be -- is real but is not what the record is for, so it
# is shown only on request.
.geometry_plan_lines <- function(x, detail = FALSE) {
  metric <- x$metric_schedule$metric
  metric_status <- if (is.null(metric)) {
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

# The worked example bundle --------------------------------------------------

#' @export
print.effect_example_effects <- function(x, ...) {
  truth <- x$truth
  .pf_emit("effect_example_effects", list(
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
  ))
  invisible(x)
}

# Measurement bridges, stores, and diagnostics -------------------------------

#' @export
format.effect_measurement_bridge <- function(x, ...) {
  .pf_inline("effect_measurement_bridge",
    paste0(x$common_space$n_measurements, " common measurements"))
}

#' @export
print.effect_measurement_bridge <- function(x, ...) {
  .pf_emit("effect_measurement_bridge", list(
    "common space" = paste0(.pf_or(x$common_space$id, "unnamed"), ", ",
      .pf_or(x$common_space$n_measurements, "unknown"), " measurements"),
    left = paste0(.pf_or(x$left_domain$id, "unknown"), " via ",
      .pf_dim(x$left_leg)),
    right = paste0(.pf_or(x$right_domain$id, "unknown"), " via ",
      .pf_dim(x$right_leg)),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

#' @export
format.effect_measurement_store <- function(x, ...) {
  .pf_inline("effect_measurement_store", paste0(nrow(x$index), " blocks"),
    x$representation)
}

#' @export
print.effect_measurement_store <- function(x, ...) {
  .pf_emit("effect_measurement_store", list(
    blocks = nrow(x$index),
    representation = x$representation,
    codec = x$codec,
    elements = .pf_num(sum(x$index$length_elements)),
    path = if (is.null(x$path)) "in memory" else "on disk",
    signature = .pf_sig(x$signature),
    `next` = "effect_coupling(form) or rdm(form) to read blocks"
  ))
  invisible(x)
}

#' @export
format.effect_measurement_receipt <- function(x, ...) {
  .pf_inline("effect_measurement_receipt", .pf_sig(x$route),
    .pf_or(x$execution$completion_status, "unknown"))
}

#' @export
print.effect_measurement_receipt <- function(x, ...) {
  .pf_emit("effect_measurement_receipt", list(
    route = .pf_sig(x$route),
    status = .pf_or(x$execution$completion_status, "unknown"),
    "plan id" = .pf_sig(x$scientific_plan_id),
    derivation = .pf_set(names(x$derivation), max = 3L),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

#' @export
format.effect_measurement_diagnostics <- function(x, ...) {
  .pf_inline("effect_measurement_diagnostics", x$method,
    paste0("effective rank ", x$experimental_effective_rank))
}

#' @export
print.effect_measurement_diagnostics <- function(x, ...) {
  .pf_emit("effect_measurement_diagnostics", list(
    method = x$method,
    "effective rank" = x$experimental_effective_rank,
    blocks = paste0(nrow(x$blocks), " (",
      .pf_set(names(x$blocks), max = 3L), ")"),
    regularization = .pf_or(x$regularization$kind, "none"),
    tolerance = .pf_num(x$tolerance),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

# Data-frame-backed indices --------------------------------------------------
#
# These subclass `data.frame`; only print() is defined so that
# `format.data.frame` keeps working for internal callers.

#' @export
print.effect_measurement_block_index <- function(x, ...) {
  .pf_emit("effect_measurement_block_index", list(
    blocks = nrow(x),
    edges = .pf_set(x$edge_id, max = 3L),
    left = .pf_set(unique(x$left), max = 3L),
    right = .pf_set(unique(x$right), max = 3L),
    widths = paste0(.pf_set(as.character(unique(x$d_left)), max = 2L), " x ",
      .pf_set(as.character(unique(x$d_right)), max = 2L)),
    elements = .pf_num(sum(x$length_elements)),
    columns = .pf_set(names(x), max = 4L)
  ))
  invisible(x)
}

#' @export
print.effect_ordered_edges <- function(x, ...) {
  .pf_emit("effect_ordered_edges", list(
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
  ))
  invisible(x)
}

# Task fragments -------------------------------------------------------------

#' @export
format.effect_evidence_stages <- function(x, ...) {
  .pf_inline("effect_evidence_stages", paste0(length(x$order), " stages"),
    x$lowering)
}

#' @export
print.effect_evidence_stages <- function(x, ...) {
  .pf_emit("effect_evidence_stages", list(
    order = .pf_set(x$order, max = 3L),
    lowering = x$lowering,
    normalization = .pf_stage(x$normalization),
    transform = .pf_stage(x$transform),
    reduction = .pf_stage(x$reduction),
    projection = .pf_stage(x$projection),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

#' @export
format.effect_frame_conservation <- function(x, ...) {
  .pf_inline("effect_frame_conservation",
    if (isTRUE(x$conserved)) "conserved" else "not conserved",
    x$component, paste0("max deviation ", .pf_num(x$max_deviation)))
}

#' @export
print.effect_frame_conservation <- function(x, ...) {
  mass <- x$feature_mass
  .pf_emit("effect_frame_conservation", list(
    conserved = .pf_yn(x$conserved),
    component = x$component,
    normalization = x$normalization,
    "max deviation" = .pf_num(x$max_deviation),
    tolerance = .pf_num(x$tolerance),
    "feature mass" = paste0(length(mass), " features, ",
      .pf_num(min(mass)), " to ", .pf_num(max(mass)))
  ))
  invisible(x)
}

# Remaining query and metric records -----------------------------------------

#' @export
format.effect_query <- function(x, ...) {
  .pf_inline("effect_query", .pf_or(x$kind, "query"),
    .pf_dim(x$operator))
}

#' @export
print.effect_query <- function(x, ...) {
  .pf_emit("effect_query", list(
    kind = .pf_or(x$kind, "query"),
    operator = .pf_dim(x$operator),
    "effect space" = .pf_or(x$effect_space$basis_id, "unspecified"),
    fixed = .pf_yn(x$fixed)
  ))
  invisible(x)
}

#' @export
format.effect_neural_pair_query <- function(x, ...) {
  .pf_inline("effect_neural_pair_query", .pf_dim(x$operator),
    x$left_domain$id)
}

#' @export
print.effect_neural_pair_query <- function(x, ...) {
  .validate_neural_pair_query(x)
  .pf_emit("effect_neural_pair_query", list(
    operator = .pf_dim(x$operator),
    left = paste0(x$left_domain$id, " (", x$left_domain$n_features,
      " features)"),
    right = paste0(x$right_domain$id, " (", x$right_domain$n_features,
      " features)"),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

#' @export
format.effect_pair_difference_query <- function(x, ...) {
  .pf_inline("effect_pair_difference_query", x$kind,
    paste0(length(x$pair_labels), " pairs"))
}

#' @export
print.effect_pair_difference_query <- function(x, ...) {
  .pf_emit("effect_pair_difference_query", list(
    kind = x$kind,
    effects = paste0(length(x$effects), " (", .pf_set(x$effects, max = 4L),
      ")"),
    pairs = paste0(length(x$pair_labels), " (",
      .pf_set(x$pair_labels, max = 3L), ")"),
    coefficients = .pf_dim(x$coefficients)
  ))
  invisible(x)
}

#' @export
format.effect_metric_inverse_representation <- function(x, ...) {
  .pf_inline("effect_metric_inverse_representation", x$kind)
}

#' @export
print.effect_metric_inverse_representation <- function(x, ...) {
  .pf_emit("effect_metric_inverse_representation", list(
    kind = x$kind,
    signature = .pf_sig(x$signature),
    state = if (identical(x$kind, "none")) {
      "no inverse supplied; quadratic forms must be solved"
    } else {
      "inverse supplied; quadratic forms use it directly"
    }
  ))
  invisible(x)
}

#' @export
format.effect_metric_handle <- function(x, ...) {
  .pf_inline("effect_metric_handle", .pf_or(x$inverse_mode, "unknown"),
    paste0(length(x$factorizations), " factorizations"))
}

#' @export
print.effect_metric_handle <- function(x, ...) {
  .pf_emit("effect_metric_handle", list(
    role = .pf_or(x$metric$role, "unknown"),
    "inverse mode" = .pf_or(x$inverse_mode, "unknown"),
    shortcut = .pf_yn(x$shortcut),
    factorizations = length(x$factorizations),
    state = "resolved handle; the metric is fixed for this evaluation"
  ))
  invisible(x)
}

#' @export
format.effect_metric_components <- function(x, ...) {
  .pf_inline("effect_metric_components",
    paste0("coherent rank ", .pf_or(x$coherent_rank, "unknown")),
    .pf_or(x$inverse_quadratic_mode, "unknown"))
}

#' @export
print.effect_metric_components <- function(x, ...) {
  .pf_emit("effect_metric_components", list(
    "coherent rank" = .pf_or(x$coherent_rank, "unknown"),
    "configuration psd" = .pf_yn(x$configuration_psd),
    denominator = .pf_or(x$denominator, "unspecified"),
    "inverse mode" = .pf_or(x$inverse_quadratic_mode, "unknown"),
    factorizations = .pf_or(x$factorization_count, "0"),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

#' @export
format.effect_lowering <- function(x, ...) {
  .pf_inline("effect_lowering", x$kind,
    if (isTRUE(x$collapsed)) "collapsed" else "not collapsed")
}

#' @export
print.effect_lowering <- function(x, ...) {
  .pf_emit("effect_lowering", list(
    kind = x$kind,
    collapsed = .pf_yn(x$collapsed),
    reason = .pf_or(x$reason, "none recorded")
  ))
  invisible(x)
}

#' @export
format.effect_residual_source <- function(x, ...) {
  .pf_inline("effect_residual_source", .pf_shape(x$dim), "unread")
}

#' @export
print.effect_residual_source <- function(x, ...) {
  .pf_emit("effect_residual_source", list(
    dim = paste0(.pf_shape(x$dim), " (observations x features)"),
    revision = .pf_sig(x$stable_revision),
    state = "unread; residual blocks are pulled by feature index"
  ))
  invisible(x)
}

#' @export
format.effect_coupling_partition_policy <- function(x, ...) {
  .pf_inline("effect_coupling_partition_policy", .pf_or(x$placement),
    paste0(length(x$weights), " partition weights"))
}

#' @export
print.effect_coupling_partition_policy <- function(x, ...) {
  .validate_coupling_partition_policy(x)
  weights <- x$weights
  .pf_emit("effect_coupling_partition_policy", list(
    placement = .pf_or(x$placement),
    weights = paste0(length(weights), ", ",
      if (length(unique(weights)) == 1L) {
        paste0("equal (", .pf_num(weights[[1L]]), ")")
      } else {
        paste0("unequal, ", .pf_num(range(weights)))
      }),
    transform = .pf_or(x$transform$kind, "none"),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

# Measurement views, decompositions, and tomography --------------------------

#' @export
format.effect_measurement_view <- function(x, ...) {
  .pf_inline("effect_measurement_view", .pf_or(x$view),
    paste0(nrow(x$edge_index), " edges"))
}

#' @export
print.effect_measurement_view <- function(x, ...) {
  .validate_measurement_view(x)
  .pf_emit("effect_measurement_view", list(
    view = .pf_or(x$view),
    edges = paste0(nrow(x$edge_index), ", ", x$edge_completeness),
    dimensions = .pf_shape(x$dimensions),
    values = paste0(length(x$values), " blocks"),
    capability = .pf_or(x$result_capability, "unspecified"),
    completeness = .pf_or(x$completeness, "unknown")
  ))
  invisible(x)
}

#' @export
format.effect_measurement_decomposition <- function(x, ...) {
  .pf_inline("effect_measurement_decomposition", x$id,
    paste0(length(x$components), " components"))
}

#' @export
print.effect_measurement_decomposition <- function(x, ...) {
  .validate_measurement_decomposition(x)
  .pf_emit("effect_measurement_decomposition", list(
    id = x$id,
    components = paste0(length(x$components), " (",
      .pf_set(names(x$components), max = 3L), ")"),
    "output space" = .pf_or(x$output_space$id, "unnamed"),
    orientation = .pf_or(x$orientation, "unspecified"),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

#' @export
format.effect_tomography_result <- function(x, ...) {
  .pf_inline("effect_tomography_result", x$method, x$status,
    if (isTRUE(x$lossless)) "lossless" else "lossy")
}

#' @export
print.effect_tomography_result <- function(x, ...) {
  .validate_tomography_result(x)
  .pf_emit("effect_tomography_result", list(
    method = x$method,
    status = x$status,
    lossless = .pf_yn(x$lossless),
    certified = .pf_yn(x$certified),
    projections = paste0("left ", .pf_dim(x$left_projection), ", right ",
      .pf_dim(x$right_projection)),
    reference = .pf_sig(x$reference_signature),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

# Measurement plan fragments -------------------------------------------------

#' @export
format.effect_measurement_plan <- function(x, ...) {
  .pf_inline("effect_measurement_plan", paste0(nrow(x$block_index), " blocks"),
    x$query_construction, x$edge_scope)
}

#' @export
print.effect_measurement_plan <- function(x, ...) {
  .pf_emit("effect_measurement_plan", list(
    blocks = nrow(x$block_index),
    relations = if (isTRUE(x$same_relation)) "one" else "two (cross-relation)",
    query = paste0(x$query_role, ", ", x$query_construction),
    "sampling axis" = .pf_or(x$sampling_axis, "not declared"),
    edges = paste0(nrow(x$edges), " (", x$edge_scope, ")"),
    "partition products" = paste0(nrow(x$partition_products), ", ",
      x$partition_expansion),
    regularization = .pf_or(x$regularization$kind, "none"),
    tolerance = .pf_num(x$tolerance),
    "plan id" = .pf_sig(x$scientific_plan_id),
    state = "query-first; no measurement block has been computed"
  ))
  invisible(x)
}

#' @export
format.effect_measurement_leg <- function(x, ...) {
  .pf_inline("effect_measurement_leg", x$kind, .pf_dim(x$operator))
}

#' @export
print.effect_measurement_leg <- function(x, ...) {
  .pf_emit("effect_measurement_leg", list(
    kind = x$kind,
    operator = .pf_dim(x$operator),
    domain = .pf_or(x$source_domain$id, "unknown"),
    "output space" = .pf_or(x$output_space$id, "unnamed"),
    support = paste0(length(x$support), " features"),
    estimation = .pf_or(x$estimation, "fixed"),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

#' @export
format.effect_measurement_axis <- function(x, ...) {
  .pf_inline("effect_measurement_axis", x$id,
    paste0(length(x$coordinates), " coordinates"))
}

#' @export
print.effect_measurement_axis <- function(x, ...) {
  .validate_measurement_axis(x)
  units <- unique(as.character(x$units))
  .pf_emit("effect_measurement_axis", list(
    id = x$id,
    coordinates = paste0(length(x$coordinates), " (",
      .pf_set(x$coordinates, max = 4L), ")"),
    basis = .pf_or(x$basis_id, "unspecified"),
    units = if (length(units) == 1L) units else .pf_set(units, max = 3L),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

#' @export
format.effect_measurement_edges <- function(x, ...) {
  .pf_inline("effect_measurement_edges", paste0(nrow(x$edges), " edges"),
    if (isTRUE(x$weighted)) "weighted" else "unweighted")
}

#' @export
print.effect_measurement_edges <- function(x, ...) {
  .pf_emit("effect_measurement_edges", list(
    edges = nrow(x$edges),
    from = .pf_set(unique(x$edges[[1L]]), max = 3L),
    to = .pf_set(unique(x$edges[[2L]]), max = 3L),
    weighted = .pf_yn(x$weighted),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

#' @export
format.effect_measurement_regularization <- function(x, ...) {
  .pf_inline("effect_measurement_regularization", x$kind,
    if (isTRUE(x$applied)) "applied" else "not applied")
}

#' @export
print.effect_measurement_regularization <- function(x, ...) {
  .validate_measurement_regularization(x)
  .pf_emit("effect_measurement_regularization", list(
    kind = x$kind,
    applied = .pf_yn(x$applied),
    "lambda left" = .pf_num(x$lambda_left, empty = "none"),
    "lambda right" = .pf_num(x$lambda_right, empty = "none"),
    signature = .pf_sig(x$signature)
  ))
  invisible(x)
}

#' @export
format.effect_observation_whitener <- function(x, ...) {
  .pf_inline("effect_observation_whitener", .pf_or(x$kind, "identity"),
    .pf_shape(x$dim))
}

#' @export
print.effect_observation_whitener <- function(x, ...) {
  .pf_emit("effect_observation_whitener", list(
    kind = .pf_or(x$kind, "identity"),
    dim = .pf_shape(x$dim),
    signature = .pf_sig(x$signature),
    state = "treated as fixed and known when the relation was fitted"
  ))
  invisible(x)
}

# Matched and control pair couplings -----------------------------------------

#' @export
format.effect_pair_coupling <- function(x, ...) {
  .pf_inline("effect_pair_coupling", x$kind,
    paste0(sum(x$eligible), " eligible pairs"))
}

#' @export
print.effect_pair_coupling <- function(x, ...) {
  .validate_pair_coupling(x)
  .pf_emit("effect_pair_coupling", list(
    kind = x$kind,
    shape = .pf_dim(x$value),
    left = paste0(length(x$left_space$coordinates), " coordinates, ",
      .pf_or(x$left_space$basis_id, "unspecified")),
    right = paste0(length(x$right_space$coordinates), " coordinates, ",
      .pf_or(x$right_space$basis_id, "unspecified")),
    eligible = paste0(sum(x$eligible), " of ", length(x$eligible), " pairs"),
    mass = paste0("total ", .pf_num(sum(x$value)), ", max ",
      .pf_num(max(x$value))),
    `next` = if (identical(x$kind, "match")) {
      "control_coupling(matches) then coupling_contrast(matches, controls)"
    } else {
      "coupling_contrast(matches, controls)"
    }
  ))
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
