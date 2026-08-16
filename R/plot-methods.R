# Base-graphics pictures of compiled geometry views -------------------------
#
# These methods draw what the view objects already contain. They never
# recompute geometry, never smooth, and never clip signed estimates at zero:
# a picture of a cross-generalized quantity must show its negative half.

# Okabe-Ito, the standard colorblind-safe categorical set. Held as literal
# hex so the palette does not depend on a grDevices version or on the
# session's palette state.
.geometry_colors <- c(
  black = "#000000", orange = "#E69F00", sky = "#56B4E9",
  green = "#009E73", yellow = "#F0E442", blue = "#0072B2",
  vermillion = "#D55E00", purple = "#CC79A7",
  neutral = "#8C8C8C", faint = "#C8C8C8"
)

# A symmetric blue-white-red ramp. Diverging, colorblind-safe, and centered
# on white so that zero is visually neutral for signed distances.
.diverging_colors <- function(n = 128L) {
  grDevices::colorRampPalette(c(
    "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
    "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
  ))(n)
}

# User `...` wins over the method's defaults, so every call remains a single
# `do.call()` without duplicate-argument errors.
.merge_graphics_args <- function(defaults, dots) {
  if (!length(dots)) return(defaults)
  named <- if (is.null(names(dots))) {
    rep(FALSE, length(dots))
  } else {
    nzchar(names(dots))
  }
  defaults[names(dots)[named]] <- dots[named]
  c(defaults, dots[!named])
}

.measurement_ids <- function(index) {
  if (is.data.frame(index)) {
    if ("measurement" %in% names(index)) {
      return(as.character(index$measurement))
    }
    return(as.character(seq_len(nrow(index))))
  }
  as.character(index)
}

# Accepts positions, a logical mask, or measurement identifiers, and always
# returns positions. An out-of-range request is an error, not a silent drop:
# a highlight set that quietly loses members would misreport ground truth.
.resolve_highlight <- function(highlight, index, n) {
  if (is.null(highlight)) return(integer())
  if (is.logical(highlight)) {
    if (length(highlight) != n || anyNA(highlight)) {
      stop(paste0(
        "A logical `highlight` must have one non-missing value per ",
        "measurement."
      ), call. = FALSE)
    }
    return(which(highlight))
  }
  if (is.character(highlight)) {
    positions <- match(highlight, .measurement_ids(index))
    if (anyNA(positions)) {
      stop("`highlight` names measurements that are not in this view.",
        call. = FALSE)
    }
    return(sort(unique(positions)))
  }
  if (!is.numeric(highlight) || anyNA(highlight) ||
      any(!is.finite(highlight)) || any(highlight %% 1 != 0) ||
      any(highlight < 1L) || any(highlight > n)) {
    stop(sprintf(
      "`highlight` must select measurement positions between 1 and %d.", n
    ), call. = FALSE)
  }
  sort(unique(as.integer(highlight)))
}

.resolve_measurement <- function(measurement, index, n) {
  if (is.character(measurement)) {
    if (length(measurement) != 1L || is.na(measurement)) {
      stop("`measurement` must be one measurement position or identifier.",
        call. = FALSE)
    }
    position <- match(measurement, .measurement_ids(index))
    if (is.na(position)) {
      stop("`measurement` does not identify a measurement in this view.",
        call. = FALSE)
    }
    return(position)
  }
  if (!is.numeric(measurement) || length(measurement) != 1L ||
      is.na(measurement) || !is.finite(measurement) ||
      measurement %% 1 != 0 || measurement < 1L || measurement > n) {
    stop(sprintf(
      "`measurement` must be one measurement position between 1 and %d.", n
    ), call. = FALSE)
  }
  as.integer(measurement)
}

.validate_highlight_label <- function(highlight_label) {
  if (!is.character(highlight_label) || length(highlight_label) != 1L ||
      is.na(highlight_label) || !nzchar(highlight_label)) {
    stop("`highlight_label` must be one non-empty character string.",
      call. = FALSE)
  }
  highlight_label
}

.finite_range <- function(...) {
  value <- unlist(list(...), use.names = FALSE)
  value <- value[is.finite(value)]
  if (!length(value)) return(c(0, 1))
  out <- range(value)
  if (out[[1L]] == out[[2L]]) out <- out + c(-0.5, 0.5)
  out
}

.main_or <- function(main, default) {
  if (is.null(main)) return(default)
  if (!is.character(main) || !length(main) || anyNA(main)) {
    stop("`main` must be one or more non-missing character titles.",
      call. = FALSE)
  }
  main
}

.panel_main <- function(main, panel, default) {
  if (is.null(main)) return(default)
  main <- .main_or(main, default)
  main[[((panel - 1L) %% length(main)) + 1L]]
}

# One measurement-index panel, shared by the contrast profile, the RSA
# coefficients, and the crossnobis statistic.
.plot_measurement_profile <- function(values, highlight, xlab, ylab, main,
                                      color, dots = list(),
                                      reference = 0, headroom = 0,
                                      legend_label = NULL) {
  positions <- seq_along(values)
  ylim <- .finite_range(values, reference)
  ylim[[2L]] <- ylim[[2L]] + headroom * diff(ylim)
  defaults <- list(
    x = positions, y = values, type = "h",
    col = .geometry_colors[["faint"]], cex.main = 1,
    xlab = xlab, ylab = ylab, main = main, ylim = ylim
  )
  do.call(graphics::plot.default, .merge_graphics_args(defaults, dots))
  if (is.finite(reference)) {
    graphics::abline(h = reference, lty = 3, col = .geometry_colors[["black"]])
  }
  graphics::points(positions, values, pch = 20, cex = 0.5,
    col = .geometry_colors[["neutral"]])
  if (length(highlight)) {
    graphics::segments(highlight, reference, highlight, values[highlight],
      col = color, lwd = 2)
    graphics::points(highlight, values[highlight], pch = 19, cex = 0.9,
      col = color)
    if (!is.null(legend_label)) {
      graphics::legend("topleft", bty = "n", cex = 0.8, pch = 19, col = color,
        legend = sprintf("%s (%d measurements)", legend_label,
          length(highlight)))
    }
  }
  invisible(NULL)
}

# `signed` is a plain vector for a self view and a named list of endpoint
# marginals for a coupling view; both are legal, so name which one is shown.
.contrast_signed <- function(x) {
  signed <- x$signed
  if (is.list(signed)) {
    if (!length(signed)) {
      return(list(values = rep(NA_real_, length(x$total)), label = "signed"))
    }
    name <- names(signed)[[1L]]
    if (is.null(name) || !nzchar(name)) name <- "endpoint 1"
    return(list(values = as.numeric(signed[[1L]]),
      label = paste0("signed (", name, ")")))
  }
  list(values = as.numeric(signed), label = "signed")
}

#' Plot geometry views
#'
#' Base-graphics pictures of the objects returned by [contrast_energy()],
#' [rdm()], [rsa()], [crossnobis()], and [rdm_sampling_covariance()]. Each
#' method draws only what the view already holds: nothing is recomputed,
#' nothing is smoothed, and negative crossvalidated estimates are drawn where
#' they fall rather than clipped at zero.
#'
#' `plot()` on an `effect_contrast_view` is the fastest way to see the
#' central claim of the package. The left panel places every measurement in
#' the coherent/configuration plane, where the dashed guides are lines of
#' constant `total` because `coherent + configuration = total` exactly. The
#' right panel walks the same total along the measurement index. Passing
#' `highlight` (for the generated fixture,
#' `example_fmri_effects()$truth$signal_measurements`) fills the
#' measurements that overlap the planted signal, so a correct analysis is
#' visible rather than asserted.
#'
#' @param x A view object: `effect_contrast_view`, `effect_rdm_view`,
#'   `effect_rsa_view`, `effect_crossnobis_view`, or
#'   `effect_rdm_sampling_covariance`.
#' @param which Which contrast panels to draw: `"both"` (the default
#'   two-panel figure), `"decomposition"` for the coherent/configuration
#'   plane alone, or `"profile"` for the total-energy index plot alone.
#' @param highlight Optional measurements to mark: positions, a logical mask
#'   with one value per measurement, or measurement identifiers. Measurements
#'   outside the view are an error rather than a silent drop.
#' @param highlight_label What the highlighted measurements are, named in the
#'   legend alongside their count, as in `"planted signal (57 measurements)"`.
#'   Say what the set means rather than that it was selected.
#' @param measurement Which measurement to draw. For an `effect_rdm_view`,
#'   one position or identifier, or the string `"mean"` for the mean
#'   dissimilarity matrix over measurements. The default is the measurement
#'   with the largest mean distance.
#' @param annotate Whether to print the numeric distance inside each cell of
#'   an RDM heatmap. The default prints them when there are at most eight
#'   experimental effects.
#' @param terms Optional character vector selecting which RSA coefficient
#'   columns to draw. The default draws every column, including the
#'   intercept and any nuisance terms.
#' @param estimate Optional point estimates to center the sampling intervals
#'   on: a numeric vector over the distance labels carried by `x`, or an
#'   `effect_rdm_view`, from which the row for the measurement `x` describes
#'   is taken. The default centers the intervals on zero, which shows the
#'   interval half-width implied by the declared calibration target.
#' @param level Two-sided normal coverage for the sampling interval. This is
#'   a width computed from the analytic sampling covariance, not a
#'   calibrated confidence interval or a group-level test.
#' @param sort Whether to order the distance rows by their plotted center.
#' @param main Optional title. For a multi-panel figure a character vector is
#'   recycled across panels.
#' @param ... Further arguments passed to the underlying base-graphics call
#'   for each panel, such as `pch`, `cex`, `xlim`, or `col`.
#' @return The view `x`, invisibly. Called for the picture it draws.
#' @seealso [contrast_energy()], [rdm()], [rsa()], [crossnobis()], and
#'   [rdm_sampling_covariance()] for the objects these methods draw, and
#'   [as.data.frame()] methods for the same values as a table.
#' @family geometry plans and views
#' @name plot_views
#' @importFrom graphics abline axis box image layout lcm legend mtext par
#'   plot.default points rect segments text
#' @importFrom grDevices adjustcolor colorRampPalette
#' @examples
#' example <- example_fmri_effects()
#' plan <- plan_geometry(
#'   example$fit$relation, example$frame,
#'   cross_partitions(
#'     example$fit$relation,
#'     independence = "independent", generalizes_over = "run"
#'   )
#' )
#' planted <- example$truth$signal_measurements
#'
#' # Where does the animate-versus-inanimate contrast reproduce, and how does
#' # its energy split? The planted measurements fill in.
#' energy <- contrast_energy(plan, example$contrast)
#' plot(energy, highlight = planted, highlight_label = "planted signal")
#'
#' # The signed squared distances at the strongest measurement.
#' distances <- rdm(plan)
#' plot(distances, measurement = which.max(energy$total))
#'
#' # One panel per RSA coefficient, over measurements.
#' plot(rsa(plan, models = list(category = example$model_rdm)),
#'      highlight = planted)
#'
#' # The same estimand under a declared fixed noise metric.
#' metric <- noise_precision(
#'   diag(example$domain$n_features), example$domain
#' )
#' mahalanobis_plan <- plan_geometry(
#'   example$fit$relation, example$frame,
#'   cross_partitions(example$fit$relation, independence = "independent"),
#'   metric = metric
#' )
#' plot(crossnobis(mahalanobis_plan, example$contrast), highlight = planted)
#'
#' # Analytic within-measurement uncertainty for those distances.
#' uncertainty <- rdm_sampling_covariance(
#'   plan, example$fit, target = "null", at = which.max(energy$total)
#' )
#' plot(uncertainty, estimate = distances)
#' @rdname plot_views
#' @export
plot.effect_contrast_view <- function(x,
                                      which = c("both", "decomposition",
                                        "profile"),
                                      highlight = NULL,
                                      highlight_label = "highlighted",
                                      main = NULL, ...) {
  which <- match.arg(which)
  dots <- list(...)
  total <- as.numeric(x$total)
  coherent <- as.numeric(x$coherent)
  configuration <- as.numeric(x$configuration)
  signed <- .contrast_signed(x)
  highlight <- .resolve_highlight(highlight, x$index, length(total))
  highlight_label <- .validate_highlight_label(highlight_label)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  if (identical(which, "both")) {
    graphics::par(mfrow = c(1L, 2L))
  }
  graphics::par(mar = c(4.4, 4.4, 3.2, 1.1))

  if (which %in% c("both", "decomposition")) {
    sign_color <- ifelse(
      !is.finite(signed$values) | signed$values == 0,
      .geometry_colors[["neutral"]],
      ifelse(signed$values > 0, .geometry_colors[["blue"]],
        .geometry_colors[["vermillion"]])
    )
    marked <- logical(length(total))
    marked[highlight] <- TRUE
    # Every point keeps its sign color. Highlighting is carried by weight and
    # a ring, not by a color of its own: unhighlighted measurements, which
    # pile up near the origin, are simply drawn lighter and smaller.
    fill <- sign_color
    sign_color[!marked] <- grDevices::adjustcolor(
      sign_color[!marked], alpha.f = if (length(highlight)) 0.35 else 0.7
    )
    sign_color[marked] <- .geometry_colors[["black"]]
    xlim <- .finite_range(coherent)
    ylim <- .finite_range(configuration)
    ylim[[2L]] <- ylim[[2L]] + 0.2 * diff(ylim)
    defaults <- list(
      x = coherent, y = configuration,
      col = sign_color, bg = fill, pch = ifelse(marked, 21L, 19L),
      cex = ifelse(marked, 1.05, 0.65),
      xlim = xlim, ylim = ylim, cex.main = 1,
      xlab = "Coherent energy (region-average contrast)",
      ylab = "Configuration energy (pattern beyond that mean)",
      main = .panel_main(main, 1L, "Coherent versus configuration energy")
    )
    do.call(graphics::plot.default, .merge_graphics_args(defaults, dots))
    graphics::abline(h = 0, v = 0, lty = 3, col = .geometry_colors[["black"]])

    # coherent + configuration = total exactly, so a constant total is a
    # line of slope -1. Two guides make the partition readable.
    guides <- unique(signif(stats::quantile(
      total[is.finite(total)], c(0.9, 0.99), names = FALSE
    ), 2))
    for (guide in guides[is.finite(guides)]) {
      graphics::abline(a = guide, b = -1, lty = 2,
        col = .geometry_colors[["neutral"]])
      # Label the guide where it enters the panel from the upper left, which
      # is the sparse corner of a coherent/configuration cloud.
      entry <- max(xlim[[1L]], guide - ylim[[2L]])
      exit <- min(xlim[[2L]], guide - ylim[[1L]])
      if (entry >= exit) next
      label_x <- entry + 0.05 * (exit - entry)
      graphics::text(label_x, guide - label_x, labels = format(guide),
        adj = c(-0.1, -0.4), cex = 0.7, col = .geometry_colors[["neutral"]])
    }
    graphics::legend(
      "topright", bty = "n", cex = 0.8,
      legend = c(
        paste0(signed$label, " > 0"), paste0(signed$label, " < 0"),
        if (length(highlight)) {
          sprintf("%s (%d measurements)", highlight_label, length(highlight))
        }
      ),
      col = c(.geometry_colors[["blue"]], .geometry_colors[["vermillion"]],
        if (length(highlight)) .geometry_colors[["black"]]),
      pt.bg = .geometry_colors[["neutral"]],
      pch = c(19L, 19L, if (length(highlight)) 21L),
      pt.cex = c(0.9, 0.9, if (length(highlight)) 1.05)
    )
    # The dashed guides no longer need a legend row of their own.
    graphics::mtext("dashed lines: constant total energy", side = 3,
      line = 0.2, cex = 0.75, col = .geometry_colors[["neutral"]])
  }

  if (which %in% c("both", "profile")) {
    .plot_measurement_profile(
      total, highlight,
      xlab = "Measurement (searchlight, region, or voxel)",
      ylab = "Total contrast energy reproducing across partitions",
      main = .panel_main(main, if (identical(which, "both")) 2L else 1L,
        "Total energy across measurements"),
      color = .geometry_colors[["vermillion"]], dots = dots,
      headroom = if (length(highlight)) 0.12 else 0,
      legend_label = highlight_label
    )
  }
  invisible(x)
}

# Rebuild the symmetric condition-by-condition matrix from the packed pair
# columns the view actually carries. Pairs the view did not request stay NA
# and are simply not drawn.
.rdm_matrix <- function(x, values) {
  pairs <- x$pairs
  effects <- unique(c(pairs$left, pairs$right))
  q <- length(effects)
  out <- matrix(NA_real_, q, q, dimnames = list(effects, effects))
  diag(out) <- 0
  rows <- match(pairs$left, effects)
  columns <- match(pairs$right, effects)
  out[cbind(rows, columns)] <- values
  out[cbind(columns, rows)] <- values
  out
}

#' @rdname plot_views
#' @export
plot.effect_rdm_view <- function(x, measurement = NULL, annotate = NULL,
                                 main = NULL, ...) {
  values <- as.matrix(x$values)
  dots <- list(...)
  if (is.character(measurement) && length(measurement) == 1L &&
      !is.na(measurement) && identical(measurement, "mean")) {
    row <- colMeans(values, na.rm = TRUE)
    label <- sprintf("mean of %d measurements", nrow(values))
  } else {
    if (is.null(measurement)) {
      means <- rowMeans(values, na.rm = TRUE)
      means[!is.finite(means)] <- -Inf
      measurement <- which.max(means)
    }
    position <- .resolve_measurement(measurement, x$index, nrow(values))
    row <- values[position, ]
    label <- sprintf("measurement %s",
      .measurement_ids(x$index)[[position]])
  }

  matrix_values <- .rdm_matrix(x, row)
  effects <- rownames(matrix_values)
  q <- length(effects)
  if (is.null(annotate)) annotate <- q <= 8L
  if (!is.logical(annotate) || length(annotate) != 1L || is.na(annotate)) {
    stop("`annotate` must be TRUE or FALSE.", call. = FALSE)
  }
  extreme <- max(abs(matrix_values[is.finite(matrix_values)]), 0)
  if (!is.finite(extreme) || extreme <= 0) extreme <- 1
  limits <- c(-extreme, extreme)
  colors <- .diverging_colors()

  op <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(1L)
    graphics::par(op)
  }, add = TRUE)
  # `asp = 1` keeps the cells square, so the heatmap is limited by whichever
  # panel dimension is smaller. Give the key a fixed physical width, just
  # enough for its strip, axis, and title, and the margins only the room the
  # effect names need, so the matrix itself takes everything that is left.
  graphics::layout(matrix(c(1L, 2L), nrow = 1L),
    widths = c(1, graphics::lcm(2.6)))
  left <- min(10, 0.55 * max(nchar(effects)) + 1.4)
  graphics::par(mar = c(left, left, 3.2, 0.6))

  # image() reads z[i, j] at (x[i], y[j]); transposing and reversing the
  # columns puts effect 1 in the top-left cell, as the matrix reads.
  drawn <- t(matrix_values)[, rev(seq_len(q)), drop = FALSE]
  defaults <- list(
    x = seq_len(q), y = seq_len(q), z = drawn, zlim = limits,
    col = colors, axes = FALSE, asp = 1, xlab = "", ylab = "",
    main = .panel_main(main, 1L,
      sprintf("Signed squared distances at %s", label))
  )
  do.call(graphics::image, .merge_graphics_args(defaults, dots))
  # `asp = 1` letterboxes the image inside the plot region, so the frame, the
  # separators, and the effect names must hug the cells, not the region.
  graphics::text(seq_len(q), 0.5 - 0.15, labels = effects, srt = 90,
    adj = c(1, 0.5), cex = 0.9, xpd = NA)
  graphics::text(0.5 - 0.15, rev(seq_len(q)), labels = effects,
    adj = c(1, 0.5), cex = 0.9, xpd = NA)
  edges <- seq(0.5, q + 0.5, by = 1)
  graphics::segments(edges, 0.5, edges, q + 0.5, col = "#FFFFFF")
  graphics::segments(0.5, edges, q + 0.5, edges, col = "#FFFFFF")
  graphics::rect(0.5, 0.5, q + 0.5, q + 0.5,
    border = .geometry_colors[["neutral"]])
  graphics::mtext(
    sprintf("%s component; distances are signed and may be negative",
      x$component),
    side = 3, line = 0.2, cex = 0.8, col = .geometry_colors[["neutral"]]
  )
  if (annotate) {
    for (column in seq_len(q)) {
      for (row_position in seq_len(q)) {
        value <- drawn[column, row_position]
        if (!is.finite(value)) next
        graphics::text(column, row_position, labels = format(value,
          digits = 2, trim = TRUE), cex = 0.8,
          col = if (abs(value) > 0.65 * extreme) "#FFFFFF" else "#000000")
      }
    }
  }

  key <- seq(limits[[1L]], limits[[2L]], length.out = length(colors))
  graphics::par(mar = c(left, 0.4, 3.2, 3.4))
  graphics::image(x = 1, y = key, z = matrix(key, nrow = 1L), zlim = limits,
    col = colors, axes = FALSE, xlab = "", ylab = "")
  graphics::axis(4, las = 1, cex.axis = 0.8)
  graphics::box(col = .geometry_colors[["neutral"]])
  graphics::mtext("Signed squared distance", side = 4, line = 2.3, cex = 0.8)
  invisible(x)
}

#' @rdname plot_views
#' @export
plot.effect_rsa_view <- function(x, terms = NULL, highlight = NULL,
                                 highlight_label = "highlighted",
                                 main = NULL, ...) {
  coefficients <- as.matrix(x$coefficients)
  labels <- colnames(coefficients)
  if (is.null(labels)) {
    labels <- paste0("coefficient", seq_len(ncol(coefficients)))
  }
  if (!is.null(terms)) {
    if (!is.character(terms) || !length(terms) || anyNA(terms) ||
        !all(terms %in% labels)) {
      stop("`terms` must name RSA coefficient columns of this view.",
        call. = FALSE)
    }
    coefficients <- coefficients[, terms, drop = FALSE]
    labels <- terms
  }
  highlight <- .resolve_highlight(highlight, x$index, nrow(coefficients))
  highlight_label <- .validate_highlight_label(highlight_label)
  panels <- ncol(coefficients)
  roles <- stats::setNames(x$terms$role, x$terms$term)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  if (panels > 1L) {
    columns <- ceiling(sqrt(panels))
    graphics::par(mfrow = c(ceiling(panels / columns), columns))
  }
  graphics::par(mar = c(4.4, 4.4, 3.2, 1.1))
  for (panel in seq_len(panels)) {
    term <- labels[[panel]]
    role <- if (term %in% names(roles)) roles[[term]] else "coefficient"
    default_main <- switch(role,
      intercept = term,
      model = sprintf("%s (model RDM)", term),
      nuisance = sprintf("%s (nuisance RDM)", term),
      term
    )
    .plot_measurement_profile(
      as.numeric(coefficients[, panel]), highlight,
      xlab = "Measurement (searchlight, region, or voxel)",
      ylab = "Fitted RSA coefficient (distance per model unit)",
      main = .panel_main(main, panel, default_main),
      color = .geometry_colors[["vermillion"]], dots = list(...),
      headroom = if (length(highlight) && panel == 1L) 0.12 else 0,
      # One legend for the figure: repeating it on every model panel would
      # say the same thing four times.
      legend_label = if (panel == 1L) highlight_label else NULL
    )
  }
  invisible(x)
}

#' @rdname plot_views
#' @export
plot.effect_crossnobis_view <- function(x, highlight = NULL,
                                        highlight_label = "highlighted",
                                        main = NULL, ...) {
  values <- as.numeric(x$values)
  highlight <- .resolve_highlight(highlight, x$index, length(values))
  highlight_label <- .validate_highlight_label(highlight_label)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(4.4, 4.4, 3.2, 1.1))
  .plot_measurement_profile(
    values, highlight,
    xlab = "Measurement (searchlight, region, or voxel)",
    ylab = "Crossvalidated squared Mahalanobis contrast",
    main = .main_or(main,
      "Crossnobis contrast along the measurement index")[[1L]],
    color = .geometry_colors[["vermillion"]], dots = list(...),
    headroom = if (length(highlight)) 0.12 else 0,
    legend_label = highlight_label
  )
  graphics::mtext(
    "Signed estimate; about half fall below zero where there is no effect",
    side = 3, line = 0.2, cex = 0.8, col = .geometry_colors[["neutral"]]
  )
  invisible(x)
}

# The covariance object carries uncertainty for one measurement only; it does
# not carry the distances themselves. Centering on a point estimate is
# therefore an explicit argument, and the estimate must come from the same
# measurement.
.sampling_center <- function(x, estimate) {
  labels <- x$labels
  if (is.null(estimate)) {
    return(list(values = rep(0, length(labels)), centered = FALSE))
  }
  if (inherits(estimate, "effect_rdm_view")) {
    values <- as.matrix(estimate$values)
    columns <- colnames(values)
    if (is.null(columns) || !setequal(columns, labels)) {
      stop(paste0(
        "The `estimate` view reports different distance pairs than this ",
        "sampling covariance."
      ), call. = FALSE)
    }
    position <- match(x$source$node, .measurement_ids(estimate$index))
    if (is.na(position)) {
      stop(sprintf(paste0(
        "The `estimate` view does not contain measurement %s, which this ",
        "sampling covariance describes."
      ), x$source$node), call. = FALSE)
    }
    return(list(values = as.numeric(values[position, labels]),
      centered = TRUE))
  }
  if (!is.numeric(estimate) || length(estimate) != length(labels) ||
      anyNA(estimate) || any(!is.finite(estimate))) {
    stop(sprintf(
      "`estimate` must be %d finite distances or an `effect_rdm_view`.",
      length(labels)
    ), call. = FALSE)
  }
  if (!is.null(names(estimate))) {
    if (!setequal(names(estimate), labels)) {
      stop("Named `estimate` values must name every distance exactly once.",
        call. = FALSE)
    }
    estimate <- estimate[labels]
  }
  list(values = as.numeric(estimate), centered = TRUE)
}

#' @rdname plot_views
#' @export
plot.effect_rdm_sampling_covariance <- function(x, estimate = NULL,
                                                level = 0.95, sort = TRUE,
                                                main = NULL, ...) {
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be one number strictly between 0 and 1.",
      call. = FALSE)
  }
  if (!is.logical(sort) || length(sort) != 1L || is.na(sort)) {
    stop("`sort` must be TRUE or FALSE.", call. = FALSE)
  }
  standard_errors <- sqrt(sampling_covariance(x, "diagonal"))
  center <- .sampling_center(x, estimate)
  labels <- x$labels
  half <- stats::qnorm(1 - (1 - level) / 2) * as.numeric(standard_errors)
  order_by <- if (sort) {
    order(if (center$centered) center$values else half)
  } else {
    rev(seq_along(labels))
  }
  values <- center$values[order_by]
  half <- half[order_by]
  labels <- labels[order_by]
  positions <- seq_along(labels)

  policy <- if (identical(x$plan$target$target, "null")) "null" else "plugin"
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  left <- min(20, 0.55 * max(nchar(labels)) + 3)
  graphics::par(mar = c(4.6, left, 3.6, 1.6))

  defaults <- list(
    x = values, y = positions, pch = 19, cex = 0.9,
    col = .geometry_colors[["blue"]],
    xlim = .finite_range(values - half, values + half, 0),
    ylim = c(0.5, length(positions) + 0.5), yaxt = "n",
    xlab = if (center$centered) {
      sprintf("Signed squared distance with %g%% sampling interval",
        100 * level)
    } else {
      sprintf("%g%% sampling interval (no point estimate supplied)",
        100 * level)
    },
    ylab = "",
    main = .main_or(main, sprintf(
      "Distance uncertainty at measurement %s", x$source$node
    ))[[1L]]
  )
  do.call(graphics::plot.default, .merge_graphics_args(defaults, list(...)))
  graphics::abline(v = 0, lty = 3, col = .geometry_colors[["black"]])
  graphics::segments(values - half, positions, values + half, positions,
    col = .geometry_colors[["blue"]], lwd = 1.6)
  cap <- 0.16
  graphics::segments(values - half, positions - cap, values - half,
    positions + cap, col = .geometry_colors[["blue"]])
  graphics::segments(values + half, positions - cap, values + half,
    positions + cap, col = .geometry_colors[["blue"]])
  graphics::points(values, positions, pch = 19, cex = 0.9,
    col = .geometry_colors[["blue"]])
  graphics::axis(2, at = positions, labels = labels, las = 1, tick = FALSE,
    line = -0.5, cex.axis = 0.9)
  graphics::mtext(
    sprintf(
      "%s target, %d partitions; within-measurement only, no spatial law",
      policy, x$partitions
    ),
    side = 3, line = 0.2, cex = 0.8, col = .geometry_colors[["neutral"]]
  )
  invisible(x)
}
