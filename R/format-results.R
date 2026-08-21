# Compact printing and row-wise coercion ------------------------------------

.result_index_data_frame <- function(index) {
  if (is.data.frame(index)) {
    value <- index
    rownames(value) <- NULL
    return(value)
  }
  data.frame(measurement = index, check.names = FALSE,
    stringsAsFactors = FALSE)
}

.bind_result_values <- function(index, values) {
  if (is.null(dim(values))) values <- matrix(values, ncol = 1L)
  values <- as.matrix(values)
  if (nrow(values) != NROW(index)) {
    .input_error("Result values and measurement index have different lengths.")
  }
  if (is.null(colnames(values))) {
    colnames(values) <- paste0("value", seq_len(ncol(values)))
  }
  data.frame(.result_index_data_frame(index), values,
    check.names = FALSE, row.names = NULL)
}

# A value that must be shown whole. `.pf_emit()` truncates a line past
# `.pf_line_width`, which is right for a signature or a route but wrong for a
# contrast: `?contrast_energy` tells the reader to check exactly that line
# against their intended alignment, and a forty-effect contrast was arriving
# as `... wo...`. A wrapped value is emitted across as many lines as it needs,
# each continuation indented to the value column.
.pf_wrap <- function(items) structure(list(items = as.character(items)),
  class = "pf_wrap")

# Greedily pack already-formatted items into `, `-separated lines of at most
# `width` characters. Never breaks inside an item, so `house -0.5` stays whole
# however narrow the column is.
.pf_fill <- function(items, width) {
  if (!length(items)) {
    return("none")
  }
  lines <- character()
  current <- ""
  for (i in seq_along(items)) {
    piece <- if (i < length(items)) paste0(items[[i]], ",") else items[[i]]
    candidate <- if (nzchar(current)) paste(current, piece) else piece
    if (nchar(candidate) > width && nzchar(current)) {
      lines <- c(lines, current)
      current <- piece
    } else {
      current <- candidate
    }
  }
  c(lines, current)
}

# `.pf_emit()` with support for `.pf_wrap()` values. Plain values keep the
# existing behaviour exactly, so every view that does not wrap prints
# byte-identically to before.
.pf_emit_block <- function(class, fields) {
  keep <- !vapply(fields, is.null, logical(1))
  fields <- fields[keep]
  if (!any(vapply(fields, inherits, logical(1), "pf_wrap"))) {
    return(.pf_emit(class, fields))
  }
  cat("<", class, ">\n", sep = "")
  keys <- names(fields)
  width <- max(nchar(keys)) + 2L
  indent <- strrep(" ", width + 2L)
  for (i in seq_along(fields)) {
    value <- fields[[i]]
    head <- sprintf("  %-*s", width, paste0(keys[[i]], ":"))
    if (inherits(value, "pf_wrap")) {
      wrapped <- .pf_fill(value$items, .pf_line_width - width - 2L)
      cat(paste0(c(head, rep(indent, length(wrapped) - 1L)), wrapped,
        collapse = "\n"), "\n", sep = "")
    } else {
      line <- paste0(head, paste0(as.character(value), collapse = ""))
      if (nchar(line) > .pf_line_width) {
        line <- paste0(substr(line, 1L, .pf_line_width - 3L), "...")
      }
      cat(line, "\n", sep = "")
    }
  }
  invisible(NULL)
}

.format_result_preview <- function(x, title, fields = list(), ...) {
  data <- as.data.frame(x)
  # `.pf_emit_block()` aligns the header block; `measurements` is the longest
  # key, so adding a field never shifts the existing columns.
  .pf_emit_block(title, c(list(measurements = nrow(data)), fields))
  # Four significant digits. The default seven prints the last bits of a
  # LAPACK result, which differ across BLAS implementations and platforms
  # and would make printed output non-reproducible.
  arguments <- list(utils::head(data), row.names = FALSE, ...)
  if (!"digits" %in% names(arguments)) {
    arguments$digits <- 4L
  }
  do.call(print, arguments)
  if (nrow(data) > 6L) {
    cat("  ... ", nrow(data) - 6L, " more measurements\n", sep = "")
  }
  invisible(x)
}

.format_counted_result <- function(class, measurements, detail = NULL) {
  suffix <- if (is.null(detail) || !nzchar(detail)) "" else paste0(", ", detail)
  sprintf("<%s: %d measurements%s>", class, as.integer(measurements), suffix)
}

#' @export
format.effect_geometry_plan <- function(x, detail = FALSE, ...) {
  .validate_geometry_plan(x, deep = FALSE)
  if (isTRUE(detail)) {
    return(.geometry_plan_lines(x, detail = TRUE))
  }
  .format_counted_result(
    "effect_geometry_plan", x$measurements,
    paste0(x$logical_shape[[1L]], " effects, ",
      nrow(x$pairing), " partition pairs")
  )
}

#' @export
format.effect_geometry <- function(x, ...) {
  .validate_effect_geometry(x, probe = FALSE)
  .format_counted_result(
    "effect_geometry", x$total$dim[[1L]],
    paste0(x$logical_shape[[1L]], " effects, ",
      paste(x$storage, collapse = "+"), " storage")
  )
}

#' @export
format.effect_relation_fit <- function(x, ...) {
  .validate_relation_fit(x, deep = FALSE)
  .format_counted_result(
    "effect_relation_fit", length(x$relation$partitions),
    paste0(length(x$relation$effects), " effects, ",
      x$relation$n_features, " features")
  )
}

#' @export
format.effect_view <- function(x, ...) {
  .validate_effect_view(x)
  .format_counted_result("effect_view", nrow(x$values),
    paste0(ncol(x$values), " queries"))
}

#' @export
format.effect_crossnobis_view <- function(x, ...) {
  .format_counted_result("effect_crossnobis_view", length(x$values),
    "signed estimates")
}

#' @export
format.effect_contrast_view <- function(x, ...) {
  .format_counted_result("effect_contrast_view", length(x$total),
    "signed + energy decomposition")
}

#' @export
format.effect_rdm_view <- function(x, ...) {
  .format_counted_result("effect_rdm_view", nrow(x$values),
    paste0(ncol(x$values), " distances"))
}

#' @export
format.effect_rsa_view <- function(x, ...) {
  .format_counted_result("effect_rsa_view", nrow(x$coefficients),
    paste0(ncol(x$coefficients), " coefficients"))
}

#' @export
format.effect_spectrum_view <- function(x, ...) {
  .format_counted_result("effect_spectrum_view", nrow(x$values),
    paste0(ncol(x$values), " roots"))
}

#' @export
print.effect_geometry <- function(x, ...) {
  .validate_effect_geometry(x, probe = FALSE)
  cat("<effect_geometry>\n", sep = "")
  cat("  effects:      ", x$logical_shape[[1L]], " x ",
    x$logical_shape[[2L]], "\n", sep = "")
  cat("  measurements: ", x$total$dim[[1L]], "\n", sep = "")
  cat("  components:   total + coherent + configuration\n", sep = "")
  cat("  storage:      ", paste(x$storage, collapse = ", "), "\n", sep = "")
  cat("  estimate:     signed cross-generalized form; PSD not assumed\n",
    sep = "")
  invisible(x)
}

#' @export
print.effect_relation_fit <- function(x, ...) {
  .validate_relation_fit(x, deep = FALSE)
  capability <- relation_fit_capabilities(x)
  cat("<effect_relation_fit>\n", sep = "")
  cat("  effects:      ", length(x$relation$effects), "\n", sep = "")
  cat("  features:     ", x$relation$n_features, "\n", sep = "")
  cat("  partitions:   ", length(x$relation$partitions), "\n", sep = "")
  cat("  residuals:    ", sum(capability$residual_blocks), "/",
    nrow(capability), " partitions\n", sep = "")
  cat("  covariance:   ", sum(capability$effect_covariance), "/",
    nrow(capability), " partitions\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.effect_view <- function(x, row.names = NULL, optional = FALSE,
                                      ...) {
  .validate_effect_view(x)
  .bind_result_values(x$index, x$values)
}

#' @export
print.effect_view <- function(x, ...) {
  .validate_effect_view(x)
  .format_result_preview(x, "effect_view", ...)
  .pf_aggregation_note(x)
  invisible(x)
}

# An aggregated view prints as one row per territory, and nothing in that table
# says the rows are group sums rather than measurements -- or that a coherent
# budget belongs to one frame and cannot be compared across frames
# (`design/conservative-geometry-contract.md` section 4). `contribution()` and
# `coherence_spectrum()` are the only things that set `$metadata$aggregation`,
# so the note appears exactly on the objects it describes and no other print
# changes at all.
.pf_aggregation_note <- function(x) {
  record <- x$metadata$aggregation
  if (!is.list(record)) return(invisible(NULL))
  .pf_note(sprintf(
    "aggregated_by: `%s`; %s over %s, summed by row",
    record$aggregated_by, .msg_count(record$groups, "group"),
    .msg_count(record$measurements, "measurement")
  ))
  # A spectrum's two numeric columns are read in opposite directions, and the
  # table alone does not say which is which: `total` per group is alpha times
  # the whole-domain total by construction (section 3.1), so it reports the
  # family weighting, while the share cancels alpha exactly (section 3.2) and
  # is the scale-resolved finding. Printing the table without saying so is how
  # the panel the contract forbids gets read off it.
  if (identical(record$reduction, "coherence_spectrum")) {
    .pf_note(paste0(
      "coherence_spectrum: `total` per group is alpha times the whole-domain ",
      "total by construction, so that column is the family weighting and not ",
      "a finding; `coherence_fraction` is exactly invariant to alpha and is ",
      "the scale-resolved quantity to read"
    ))
  }
  if (isTRUE(record$frame_relative)) {
    .pf_note(sprintf(paste0(
      "frame_relative: TRUE -- %s %s a share of this frame's own mass, not of ",
      "a global one, and %s not comparable across frames"
    ), .msg_names(record$frame_relative_components),
      if (length(record$frame_relative_components) == 1L) "is" else "are",
      if (length(record$frame_relative_components) == 1L) "is" else "are"))
  }
  if (length(record$masked)) {
    .pf_note(sprintf(paste0(
      "masked: %s NA over groups; a local weighted mean is a density, and ",
      "densities do not add over a territory"
    ), .msg_names(record$masked)))
  }
  invisible(NULL)
}

#' @export
as.data.frame.effect_crossnobis_view <- function(x, row.names = NULL,
                                                 optional = FALSE, ...) {
  if (!.is_finite_numeric(x$values) || length(x$values) != NROW(x$index)) {
    .input_error("Crossnobis result values or measurement index are invalid.")
  }
  .bind_result_values(x$index,
    matrix(x$values, ncol = 1L, dimnames = list(NULL, "crossnobis")))
}

# The contrast in relation order, one `name value` item per effect. Named so
# that a positionally supplied contrast shows the alignment that was actually
# used. Capped at `max` effects, because past that a print method is no longer
# a summary; `x$weights` is the whole vector.
.pf_contrast_items <- function(weights, max = 12L) {
  if (is.null(weights) || !length(weights)) {
    return("none")
  }
  labels <- names(weights)
  if (is.null(labels) || anyNA(labels) || !all(nzchar(labels))) {
    labels <- paste0("[", seq_along(weights), "]")
  }
  items <- paste0(labels, " ",
    format(signif(as.numeric(weights), 4L), trim = TRUE))
  if (length(items) <= max) {
    return(items)
  }
  c(items[seq_len(max)], sprintf("(+%d more)", length(items) - max))
}

#' @export
print.effect_crossnobis_view <- function(x, ...) {
  .format_result_preview(x, "effect_crossnobis_view",
    fields = list(contrast = .pf_wrap(.pf_contrast_items(x$contrast))), ...)
}

#' @export
as.data.frame.effect_contrast_view <- function(x, row.names = NULL,
                                               optional = FALSE, ...) {
  signed <- if (is.list(x$signed)) {
    do.call(cbind, x$signed)
  } else {
    matrix(x$signed, ncol = 1L, dimnames = list(NULL, "signed"))
  }
  if (is.list(x$signed)) {
    colnames(signed) <- paste0("signed_", names(x$signed))
  }
  values <- cbind(
    signed,
    coherent = x$coherent,
    configuration = x$configuration,
    total = x$total,
    coherence_fraction = x$coherence_fraction
  )
  .bind_result_values(x$index, values)
}

#' @export
print.effect_contrast_view <- function(x, ...) {
  # The contrast is shown because unnamed weights are accepted positionally:
  # a reader who mis-ordered them sees the alignment that was actually used.
  .format_result_preview(x, "effect_contrast_view",
    fields = list(contrast = .pf_wrap(.pf_contrast_items(x$weights))), ...)
  # A column of NA in the preview is otherwise unexplained, and the reason is
  # a real property of the estimate rather than missing data.
  valid <- x$coherence_fraction_valid
  if (is.null(valid)) valid <- !is.na(x$coherence_fraction)
  cat(strwrap(sprintf(paste0(
    "coherence_fraction: %d of %d valid; NA where coherent and ",
    "configuration are not a nonnegative partition"
  ), sum(valid), length(valid)),
    width = .pf_line_width, prefix = "    ", initial = "  "), sep = "\n")
  .pf_aggregation_note(x)
  invisible(x)
}

#' @export
as.data.frame.effect_rdm_view <- function(x, row.names = NULL,
                                          optional = FALSE, ...) {
  .bind_result_values(x$index, x$values)
}

#' @export
print.effect_rdm_view <- function(x, ...) {
  .format_result_preview(x, "effect_rdm_view", ...)
}

#' @export
as.data.frame.effect_rsa_view <- function(x, row.names = NULL,
                                          optional = FALSE, ...) {
  .bind_result_values(x$index, x$coefficients)
}

#' @export
print.effect_rsa_view <- function(x, ...) {
  # Without a header the columns of the coefficient table are the only clue
  # to what was fitted, and `(Intercept)` next to a model name does not say
  # which names are models, which are nuisance, or over how many pairs the
  # regression ran.
  terms <- x$terms
  named <- function(kind) {
    if (is.data.frame(terms)) terms$term[terms$role == kind] else character()
  }
  pairs <- x$query$pair_labels
  effects <- x$query$effects
  .format_result_preview(x, "effect_rsa_view",
    fields = list(
      models = .pf_wrap(named("model")),
      nuisance = if (length(named("nuisance"))) {
        .pf_wrap(named("nuisance"))
      },
      intercept = if (length(named("intercept"))) "fitted" else "omitted",
      pairs = if (!is.null(pairs)) {
        sprintf("%s over %s", .msg_count(length(pairs), "distance"),
          .msg_count(length(effects), "effect"))
      },
      estimate = sprintf("OLS coefficients on the %s component", x$component)
    ), ...)
  cat(sprintf("  %-14s%s\n", "next:",
    "as.data.frame(x), plot(x, terms = ...)"))
  invisible(x)
}

#' @export
as.data.frame.effect_spectrum_view <- function(x, row.names = NULL,
                                               optional = FALSE, ...) {
  .bind_result_values(x$index, x$values)
}

#' @export
print.effect_spectrum_view <- function(x, ...) {
  .format_result_preview(x, "effect_spectrum_view", ...)
}

# Population results -----------------------------------------------------------
#
# The coefficient table in long form: one row per group node, query and term,
# with the group node index --- the sink included, marked, and carrying the
# units its column is reported in --- leading it exactly as a measurement index
# leads every other result table. Long rather than wide because the result has
# three axes and a wide table has to fold two of them into column names, which
# stops being readable at the second query.
#
# A complete-form result is the same table with the readout axis renamed: one
# row per group node, term and packed coordinate, and the coordinate's own
# `row`, `column` and `scale` carried beside it so a reader can rebuild the
# symmetric form --- or check the `sqrt(2)` --- without consulting the codec.
# The axis order differs from the query-bank case because the arrays do
# (`node x term x coordinate` against `node x query x term`), and the long
# table follows the array rather than imposing an order on it.

#' @export
as.data.frame.effect_population_result <- function(x, row.names = NULL,
                                                   optional = FALSE, ...) {
  .validate_population_result(x)
  if (identical(x$basis, "complete_form")) {
    shape <- dim(x$coefficient_forms)
    labels <- dimnames(x$coefficient_forms)
    nodes <- shape[[1L]]
    coordinates <- x$coordinates[
      rep(seq_len(shape[[3L]]), each = nodes * shape[[2L]]), , drop = FALSE
    ]
    node_position <- rep(seq_len(nodes), times = shape[[2L]] * shape[[3L]])
    term_position <- rep(rep(seq_len(shape[[2L]]), each = nodes),
      times = shape[[3L]])
    readout_position <- rep(seq_len(shape[[3L]]),
      each = nodes * shape[[2L]])
    return(data.frame(
      .result_index_data_frame(x$index)[
        rep(seq_len(nodes), times = shape[[2L]] * shape[[3L]]), , drop = FALSE
      ],
      term = rep(rep(labels[[2L]], each = nodes), times = shape[[3L]]),
      coordinate = coordinates$coordinate,
      row = coordinates$row,
      column = coordinates$column,
      scale = coordinates$scale,
      .population_result_coverage_frame(x$coverage, node_position,
        readout_position, term_position),
      estimate = as.numeric(x$coefficient_forms),
      check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE
    ))
  }
  shape <- dim(x$coefficients)
  labels <- dimnames(x$coefficients)
  nodes <- shape[[1L]]
  node_position <- rep(seq_len(nodes), times = shape[[2L]] * shape[[3L]])
  readout_position <- rep(rep(seq_len(shape[[2L]]), each = nodes),
    times = shape[[3L]])
  term_position <- rep(seq_len(shape[[3L]]), each = nodes * shape[[2L]])
  data.frame(
    .result_index_data_frame(x$index)[
      rep(seq_len(nodes), times = shape[[2L]] * shape[[3L]]), ,
      drop = FALSE
    ],
    query = rep(rep(labels[[2L]], each = nodes), times = shape[[3L]]),
    term = rep(labels[[3L]], each = nodes * shape[[2L]]),
    .population_result_coverage_frame(x$coverage, node_position,
      readout_position, term_position),
    estimate = as.numeric(x$coefficients),
    check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE
  )
}

.population_subject_list <- function(coverage, node, readout, include = TRUE) {
  selected <- coverage$availability[node, readout, ]
  if (!include) selected <- !selected
  paste(coverage$planned_subjects[selected], collapse = ",")
}

.population_result_coverage_frame <- function(coverage, node, readout, term) {
  cell <- cbind(node, readout)
  coefficient <- cbind(node, readout, term)
  data.frame(
    coverage_policy = coverage$policy,
    planned_n = length(coverage$planned_subjects),
    n = coverage$n[cell],
    fraction = coverage$fraction[cell],
    n_eff = coverage$n_eff[cell],
    mass_n_eff = coverage$mass_n_eff[node],
    design_rank = coverage$design_rank[cell],
    residual_df = coverage$residual_df[cell],
    coverage_status = coverage$status[cell],
    coefficient_estimable = coverage$coefficient_estimable[coefficient],
    exclusion_reason = coverage$exclusion_reason[coefficient],
    subject_set_id = coverage$subject_set_id[cell],
    available_subjects = mapply(.population_subject_list,
      MoreArgs = list(coverage = coverage, include = TRUE), node, readout,
      USE.NAMES = FALSE),
    excluded_subjects = mapply(.population_subject_list,
      MoreArgs = list(coverage = coverage, include = FALSE), node, readout,
      USE.NAMES = FALSE),
    check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE
  )
}

# Population views ------------------------------------------------------------
#
# One row per group node and view column, with the group node index -- the
# sink included, marked, and carrying its units -- leading it exactly as it
# leads the population result's own table. The transported component is named
# by its `population-form-v1` section 8.1 ledger name and by nothing else:
# there is no `component` column and no cell reading `coherent`, because a
# transported coherent ledger is not the group node's coherent part and a
# column header is where that confusion would start.

#' @export
as.data.frame.effect_population_view <- function(x, row.names = NULL,
                                                 optional = FALSE, ...) {
  .validate_population_view(x)
  nodes <- nrow(x$values)
  columns <- x$columns[rep(seq_len(nrow(x$columns)), each = nodes), ,
    drop = FALSE]
  names(columns)[names(columns) == "column"] <- "view_column"
  node_position <- rep(seq_len(nodes), times = ncol(x$values))
  column_position <- rep(seq_len(ncol(x$values)), each = nodes)
  cell <- cbind(node_position, column_position)
  data.frame(
    .result_index_data_frame(x$index)[
      rep(seq_len(nodes), times = ncol(x$values)), , drop = FALSE
    ],
    view = x$view,
    ledger = x$ledger,
    term = x$term,
    columns,
    coverage_policy = x$coverage$policy,
    planned_n = length(x$coverage$planned_subjects),
    n = x$coverage$n[cell],
    fraction = x$coverage$fraction[cell],
    n_eff = x$coverage$n_eff[cell],
    design_rank = x$coverage$design_rank[cell],
    residual_df = x$coverage$residual_df[cell],
    coverage_status = x$coverage$status[cell],
    coefficient_estimable = x$coverage$estimable[cell],
    exclusion_reason = x$coverage$exclusion_reason[cell],
    subject_set_id = x$coverage$subject_set_id[cell],
    available_subjects = mapply(.population_subject_list,
      MoreArgs = list(coverage = x$coverage, include = TRUE),
      node_position, column_position, USE.NAMES = FALSE),
    excluded_subjects = mapply(.population_subject_list,
      MoreArgs = list(coverage = x$coverage, include = FALSE),
      node_position, column_position, USE.NAMES = FALSE),
    estimate = as.numeric(x$values),
    check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE
  )
}

# Population uncertainty -------------------------------------------------------
#
# One layer per call, never both. The between-subject block is indexed by group
# node, query and model term; the within-subject block is indexed by group
# node, query and *participant* and has no term axis at all. There is no join
# key between them and no arithmetic that would produce one, so a single table
# would have to invent a column to hold both --- and a reader who found them in
# one frame would add them. `layer` is therefore an argument and not a
# convenience.

.population_uncertainty_between_frame <- function(x) {
  shape <- dim(x$between$estimate)
  labels <- dimnames(x$between$estimate)
  nodes <- shape[[1L]]
  data.frame(
    .result_index_data_frame(x$index)[
      rep(seq_len(nodes), times = shape[[2L]] * shape[[3L]]), , drop = FALSE
    ],
    layer = "between_subject",
    ledger = x$ledger,
    query = rep(rep(labels[[2L]], each = nodes), times = shape[[3L]]),
    term = rep(labels[[3L]], each = nodes * shape[[2L]]),
    estimator = x$between$estimator,
    assumptions = paste(x$between$assumptions, collapse = ";"),
    estimate = as.numeric(x$between$estimate),
    se = as.numeric(x$between$se),
    n = as.numeric(x$between$n),
    design_rank = as.numeric(x$between$design_rank),
    residual_df = as.numeric(x$between$residual_df),
    max_leverage = as.numeric(x$between$max_leverage),
    max_abs_adjusted_residual = as.numeric(apply(
      abs(x$between$adjusted_residual), c(1L, 2L), function(value) {
        finite <- value[is.finite(value)]
        if (length(finite)) max(finite) else NA_real_
      }
    )),
    uncertainty_status = as.character(x$between$status),
    uncertainty_reason = as.character(x$between$reason),
    transport_conditioning = x$between$conditioning$transport,
    coverage_conditioning = x$between$conditioning$coverage,
    t = as.numeric(x$between$t),
    level = x$between$level,
    lower = as.numeric(x$between$lower),
    upper = as.numeric(x$between$upper),
    calibration = x$between$calibration,
    check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE
  )
}

.population_uncertainty_within_frame <- function(x) {
  within <- x$within
  if (is.null(within) || is.null(within$admitted)) {
    refusal <- if (is.null(within)) NULL else within$refusal
    .capability_refusal(paste0(
      "This result carries no within-subject layer. Either ",
      "`estimate_population()` was not given `uncertainty`, or the layer was ",
      "refused before any group column was examined; ",
      "`population_uncertainty(x)$within$refusal$reasons` names which."
    ),
      capability = if (is.null(refusal)) {
        "transported_sampling_covariance"
      } else {
        refusal$capability
      },
      namespace = "population_uncertainty",
      reasons = if (is.null(refusal)) {
        "within_layer_absent"
      } else {
        refusal$reasons
      },
      remedies = if (is.null(refusal)) {
        paste0("Pass `uncertainty = ` to `estimate_population()`, one ",
          "`rdm_sampling_covariance()` per participant.")
      } else {
        refusal$remedies
      })
  }
  shape <- dim(within$variance)
  labels <- dimnames(within$variance)
  nodes <- shape[[1L]]
  admitted <- within$admitted[, rep(seq_len(shape[[3L]]),
    each = shape[[2L]]), drop = FALSE]
  data.frame(
    .result_index_data_frame(x$index)[
      rep(seq_len(nodes), times = shape[[2L]] * shape[[3L]]), , drop = FALSE
    ],
    layer = "within_subject",
    ledger = x$ledger,
    query = rep(rep(labels[[2L]], each = nodes), times = shape[[3L]]),
    subject = rep(labels[[3L]], each = nodes * shape[[2L]]),
    admitted = as.logical(admitted),
    source_node = as.character(within$source_node[, rep(
      seq_len(shape[[3L]]), each = shape[[2L]]), drop = FALSE]),
    coefficient = as.numeric(within$coefficient[, rep(
      seq_len(shape[[3L]]), each = shape[[2L]]), drop = FALSE]),
    variance = as.numeric(within$variance),
    se = sqrt(as.numeric(within$variance)),
    calibration = "uncalibrated",
    check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE
  )
}

#' @export
as.data.frame.effect_population_uncertainty <- function(
    x, row.names = NULL, optional = FALSE, ...,
    layer = c("between", "within")) {
  .validate_population_uncertainty(x)
  layer <- match.arg(layer)
  if (identical(layer, "between")) {
    return(.population_uncertainty_between_frame(x))
  }
  .population_uncertainty_within_frame(x)
}

# Population wild bootstrap --------------------------------------------------

#' @export
as.data.frame.effect_population_wild_bootstrap <- function(
    x, row.names = NULL, optional = FALSE, ...) {
  .validate_population_bootstrap(x)
  nodes <- nrow(x$index)
  queries <- nrow(x$queries)
  node_position <- rep(seq_len(nodes), times = queries)
  query_position <- rep(seq_len(queries), each = nodes)
  cell <- cbind(node_position, query_position)
  contrast <- paste0(
    names(x$contrast), "=", format(x$contrast, trim = TRUE),
    collapse = ";"
  )
  data.frame(
    .result_index_data_frame(x$index)[node_position, , drop = FALSE],
    query = rep(rownames(x$queries), each = nodes),
    estimator = x$estimator,
    studentization = x$studentization,
    weight_distribution = x$weight_distribution,
    contrast = contrast,
    null = x$null,
    alternative = x$alternative,
    observed_estimate = x$observed_estimate[cell],
    observed_se = x$observed_se[cell],
    observed_t = x$observed_t[cell],
    successful_replicates = x$successful_replicates[cell],
    failed_replicates = x$failed_replicates[cell],
    p_value = x$p_value[cell],
    monte_carlo_se = x$monte_carlo_se[cell],
    critical_value = x$critical_value[cell],
    reject = x$reject[cell],
    level = x$level,
    status = x$status[cell],
    reason = x$reason[cell],
    n = x$coverage$n[cell],
    design_rank = x$coverage$design_rank[cell],
    residual_df = x$coverage$residual_df[cell],
    subject_set_id = x$coverage$subject_set_id[cell],
    available_subjects = mapply(.population_subject_list,
      MoreArgs = list(coverage = x$coverage, include = TRUE),
      node_position, query_position, USE.NAMES = FALSE),
    transport_conditioning = x$conditioning$transport,
    coverage_conditioning = x$conditioning$coverage,
    seed = x$seed,
    replicates = x$replicates,
    weight_signature = x$weight_signature,
    check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE
  )
}
