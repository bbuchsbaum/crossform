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

# Above this many measurements one glyph per measurement stops being a
# picture of the data: the profile panel fills in as a solid slab and the
# decomposition panel as an undifferentiated blob. Beyond it the panels
# summarize what they cannot draw rather than drawing it badly.
.dense_measurements <- 5000L

# Below this many measurements nothing is rescaled at all, so the small
# figures the package documents are exactly what they were.
.sparse_measurements <- 500L

# Glyph size and opacity as functions of how many glyphs share one panel. The
# cube root is the compromise a two-dimensional panel wants: total ink grows
# slowly enough that a few thousand measurements stay separable, and the
# floors keep a crowded panel from fading to nothing. At or below
# `.sparse_measurements` both factors are exactly one, so the scaling is a
# no-op for the sizes the vignettes and README draw.
.measurement_scale <- function(n) {
  n <- max(as.integer(n), 1L)
  if (n <= .sparse_measurements) return(list(cex = 1, alpha = 1))
  ratio <- (.sparse_measurements / n)^(1 / 3)
  list(cex = max(0.45, ratio), alpha = max(0.25, ratio))
}

# Apply an opacity to one or many colors. `adjustcolor()` takes only a scalar
# `alpha.f`, and it would rewrite an opaque request as "#C8C8C8FF" rather
# than leaving it alone; returning the color untouched at alpha 1 keeps the
# small-measurement path byte-identical to what it drew before.
.fade <- function(color, alpha) {
  if (length(alpha) == 1L && alpha >= 1) return(color)
  channels <- grDevices::col2rgb(color)
  grDevices::rgb(
    channels[1L, ], channels[2L, ], channels[3L, ],
    alpha = round(255 * pmax(0, pmin(1, alpha))), maxColorValue = 255
  )
}

# Blue where the signed marginal is positive, vermillion where it is
# negative, neutral where it is zero or missing. NA in the condition of
# `ifelse()` would propagate, so the missing case is tested first.
.signed_colors <- function(signed) {
  ifelse(
    !is.finite(signed) | signed == 0,
    .geometry_colors[["neutral"]],
    ifelse(signed > 0, .geometry_colors[["blue"]],
      .geometry_colors[["vermillion"]])
  )
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
      .input_error(paste0(
        "A logical `highlight` must have one non-missing value per ",
        "measurement."
      ))
    }
    return(which(highlight))
  }
  if (is.character(highlight)) {
    positions <- match(highlight, .measurement_ids(index))
    if (anyNA(positions)) {
      .input_error("`highlight` names measurements that are not in this view.")
    }
    return(sort(unique(positions)))
  }
  if (!.is_finite_numeric(highlight) || anyNA(highlight) ||
      any(highlight %% 1 != 0) || any(highlight < 1L) || any(highlight > n)) {
    .input_error(sprintf(
      "`highlight` must select measurement positions between 1 and %d.", n
    ))
  }
  sort(unique(as.integer(highlight)))
}

.resolve_measurement <- function(measurement, index, n) {
  if (is.character(measurement)) {
    if (length(measurement) != 1L || is.na(measurement)) {
      .input_error(
        "`measurement` must be one measurement position or identifier."
      )
    }
    position <- match(measurement, .measurement_ids(index))
    if (is.na(position)) {
      .input_error(
        "`measurement` does not identify a measurement in this view."
      )
    }
    return(position)
  }
  if (!.is_number(measurement) || measurement %% 1 != 0 || measurement < 1L ||
      measurement > n) {
    .input_error(sprintf(
      "`measurement` must be one measurement position between 1 and %d.", n
    ))
  }
  as.integer(measurement)
}

.validate_highlight_label <- function(highlight_label) {
  .check_string(
    highlight_label, "highlight_label", what = "one non-empty character string"
  )
  highlight_label
}

# How many measurements keep a needle of their own once a profile panel is
# too crowded to draw them all. Checked on every call, not only the crowded
# ones, so a typo surfaces on the small view the reader is developing against
# rather than on the large view they finally run.
.validate_top <- function(top) {
  if (!.is_number(top)) {
    .input_error(sprintf(paste0(
      "`top` must be one positive whole number of measurements to draw ",
      "individually; received %s."
    ), .msg_value(top)))
  }
  if (top %% 1 != 0 || top < 1) {
    .input_error(sprintf(paste0(
      "`top` = %s must be a positive whole number of measurements."
    ), format(top, scientific = FALSE)))
  }
  as.integer(top)
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
    .input_error("`main` must be one or more non-missing character titles.")
  }
  main
}

.panel_main <- function(main, panel, default) {
  if (is.null(main)) return(default)
  main <- .main_or(main, default)
  main[[((panel - 1L) %% length(main)) + 1L]]
}

# The bulk of a crowded profile, drawn once instead of once per measurement:
# a shaded polygon spanning the central 95% of the values with the median
# ruled across it, and needles kept only for the `top` measurements furthest
# from the reference line. Extremity is measured in both directions on
# purpose. These are crossvalidated estimates, about half of which fall below
# zero where there is no effect, and a summary that kept only the largest
# positive values would quietly flatter the analysis.
#
# Returns the positions that kept a needle, so a caller can say how many
# measurements it drew.
.plot_profile_band <- function(values, reference, top) {
  n <- length(values)
  baseline <- if (is.finite(reference)) reference else 0
  finite <- values[is.finite(values)]
  if (length(finite)) {
    band <- stats::quantile(finite, c(0.025, 0.975), names = FALSE)
    graphics::polygon(
      c(1, n, n, 1),
      c(band[[1L]], band[[1L]], band[[2L]], band[[2L]]),
      col = .fade(.geometry_colors[["faint"]], 0.55), border = NA
    )
    center <- stats::median(finite)
    graphics::segments(1, center, n, center,
      col = .geometry_colors[["neutral"]])
  }
  ranked <- order(abs(values - baseline), decreasing = TRUE, na.last = NA)
  keep <- ranked[seq_len(min(as.integer(top), length(ranked)))]
  if (length(keep)) {
    scale <- .measurement_scale(length(keep))
    graphics::segments(keep, baseline, keep, values[keep],
      col = .fade(.geometry_colors[["faint"]], scale$alpha))
    graphics::points(keep, values[keep], pch = 20, cex = 0.5 * scale$cex,
      col = .fade(.geometry_colors[["neutral"]], scale$alpha))
  }
  keep
}

# One measurement-index panel, shared by the contrast profile, the RSA
# coefficients, and the crossnobis statistic. Returns, invisibly, whether the
# panel summarized the bulk and which measurements it drew individually.
.plot_measurement_profile <- function(values, highlight, xlab, ylab, main,
                                      color, dots = list(),
                                      reference = 0, headroom = 0,
                                      legend_label = NULL, top = 200L,
                                      highlight_fill = NULL) {
  positions <- seq_along(values)
  dense <- length(values) > .dense_measurements
  scale <- .measurement_scale(length(values))
  ylim <- .finite_range(values, reference)
  ylim[[2L]] <- ylim[[2L]] + headroom * diff(ylim)
  defaults <- list(
    x = positions, y = values, type = if (dense) "n" else "h",
    col = .fade(.geometry_colors[["faint"]], scale$alpha), cex.main = 1,
    xlab = xlab, ylab = ylab, main = main, ylim = ylim
  )
  do.call(graphics::plot.default, .merge_graphics_args(defaults, dots))
  if (is.finite(reference)) {
    graphics::abline(h = reference, lty = 3, col = .geometry_colors[["black"]])
  }
  if (dense) {
    drawn <- .plot_profile_band(values, reference, top)
  } else {
    drawn <- positions
    graphics::points(positions, values, pch = 20, cex = 0.5 * scale$cex,
      col = .fade(.geometry_colors[["neutral"]], scale$alpha))
  }
  # `highlight_fill` is how a caller that already has a per-measurement color
  # -- the contrast view's signed marginal -- keeps one meaning per color
  # across a two-panel figure. With it, highlighting is carried by a black
  # ring and the fill still means the sign, exactly as in the decomposition
  # panel. Without it the panel keeps its single accent color.
  ringed <- !is.null(highlight_fill)
  if (length(highlight)) {
    # Every highlighted measurement is drawn whatever the panel's size: the
    # highlight set is the reader's question, so it is never summarized away.
    marked <- .measurement_scale(length(highlight))
    needle <- if (ringed) highlight_fill else color
    graphics::segments(highlight, reference, highlight, values[highlight],
      col = .fade(needle, marked$alpha), lwd = 2)
    if (ringed) {
      graphics::points(highlight, values[highlight], pch = 21L,
        cex = 0.95 * marked$cex, col = .geometry_colors[["black"]],
        bg = .fade(highlight_fill, marked$alpha))
    } else {
      graphics::points(highlight, values[highlight], pch = 19,
        cex = 0.9 * marked$cex, col = .fade(color, marked$alpha))
    }
  }
  named <- length(highlight) && !is.null(legend_label)
  # The highlight swatch is drawn as the encoding itself: an unfilled black
  # ring, because the fill of a ringed point carries the sign and a gray
  # filled swatch would claim gray means "highlighted" while the same gray
  # means "not highlighted" two rows above it.
  ring_pch <- if (ringed) 21L else 19L
  ring_col <- if (ringed) .geometry_colors[["black"]] else color
  if (dense) {
    graphics::legend(
      "topleft", bty = "n", cex = 0.8,
      pch = c(15L, 20L, if (named) ring_pch),
      col = c(.fade(.geometry_colors[["faint"]], 0.55),
        .geometry_colors[["neutral"]], if (named) ring_col),
      pt.bg = NA,
      legend = c(
        sprintf("central 95 %% of %s", .msg_count(length(values),
          "measurement")),
        sprintf("%s furthest from the reference",
          .msg_count(length(drawn), "measurement")),
        if (named) {
          sprintf("%s (%d measurements)", legend_label, length(highlight))
        }
      )
    )
  } else if (named && ringed) {
    graphics::legend(
      "topleft", bty = "n", cex = 0.8,
      pch = c(20L, 21L),
      col = c(.geometry_colors[["neutral"]], .geometry_colors[["black"]]),
      pt.bg = NA, pt.cex = c(1, 0.95),
      legend = c("not highlighted",
        sprintf("%s (%d measurements)", legend_label, length(highlight)))
    )
  } else if (named) {
    graphics::legend("topleft", bty = "n", cex = 0.8, pch = 19, col = color,
      legend = sprintf("%s (%d measurements)", legend_label,
        length(highlight)))
  }
  invisible(list(dense = dense, drawn = drawn))
}

# Density shading for a panel holding more measurements than it has room for
# glyphs. Plain rectangular binning over the panel's own limits, drawn with
# `rect()`: no kernel, no smoothing parameter, and no new dependency, so the
# picture is a count of the measurements that fall in each cell and nothing
# else. Opacity follows the log count, which keeps one crowded cell near the
# origin from erasing the sparse tail, and each cell takes the color of the
# sign that dominates it so the signed reading of the panel survives.
#
# Returns the number of cells that held at least one measurement.
.plot_density_bins <- function(x, y, signed, xlim, ylim, bins = 64L) {
  usable <- is.finite(x) & is.finite(y)
  if (!any(usable)) return(invisible(0L))
  if (length(signed) != length(x)) signed <- rep(NA_real_, length(x))
  x <- x[usable]
  y <- y[usable]
  signed <- signed[usable]
  x_edges <- seq(xlim[[1L]], xlim[[2L]], length.out = bins + 1L)
  y_edges <- seq(ylim[[1L]], ylim[[2L]], length.out = bins + 1L)
  column <- findInterval(x, x_edges, rightmost.closed = TRUE,
    all.inside = TRUE)
  row <- findInterval(y, y_edges, rightmost.closed = TRUE, all.inside = TRUE)
  cell <- (row - 1L) * bins + column
  cells <- bins * bins
  count <- tabulate(cell, cells)
  positive <- tabulate(cell[is.finite(signed) & signed > 0], cells)
  negative <- tabulate(cell[is.finite(signed) & signed < 0], cells)
  filled <- which(count > 0L)
  color <- ifelse(
    positive[filled] > negative[filled], .geometry_colors[["blue"]],
    ifelse(negative[filled] > positive[filled],
      .geometry_colors[["vermillion"]], .geometry_colors[["neutral"]])
  )
  intensity <- log1p(count[filled]) / log1p(max(count))
  filled_column <- ((filled - 1L) %% bins) + 1L
  filled_row <- ((filled - 1L) %/% bins) + 1L
  graphics::rect(
    x_edges[filled_column], y_edges[filled_row],
    x_edges[filled_column + 1L], y_edges[filled_row + 1L],
    col = .fade(color, 0.15 + 0.85 * intensity), border = NA
  )
  invisible(length(filled))
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
#' `example_fmri_effects()$truth$signal_measurements`) marks the
#' measurements that overlap the planted signal, so a correct analysis is
#' visible rather than asserted.
#'
#' One color means one thing across both panels. Color is the sign of the
#' signed marginal, blue above zero and orange-red below; a black ring is
#' "highlighted", and the fill inside it still carries the sign; plain gray is
#' "not highlighted". When no `highlight` is given there is nothing for gray
#' to mean, and every measurement keeps its sign color.
#'
#' Glyph size and opacity ease down as a panel fills, so a view of a few
#' thousand measurements stays readable; at or below five hundred
#' measurements nothing is rescaled. Above five thousand measurements one
#' glyph per measurement is no longer a picture of the data, and the contrast
#' panels change what they draw rather than overplotting. The profile panel
#' shades the central 95% of the values as a band, rules the median across
#' it, and keeps needles only for the `top` measurements furthest from the
#' reference line in either direction, so the negative half of a
#' crossvalidated estimate is still shown. The decomposition panel bins the
#' coherent/configuration plane into rectangular cells whose opacity follows
#' the log count and whose color is the signed marginal that dominates the
#' cell; the binning is a plain count, with no kernel and no smoothing.
#' Highlighted measurements are always drawn individually, on top of either
#' summary, however many of them there are.
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
#'   legend alongside their count, as in `"planted signal (104 measurements)"`.
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
#' @param top How many measurements keep a needle of their own once a profile
#'   panel holds more than five thousand of them, taken as the `top`
#'   measurements furthest from the reference line in either direction. The
#'   default is 200, capped at the number of measurements the view has. The
#'   argument is validated on every call, including views small enough to
#'   draw every measurement, where it is otherwise unused.
#' @param ... Further arguments passed to the underlying base-graphics call
#'   for each panel, such as `pch`, `cex`, `xlim`, or `col`.
#' @return The view `x`, invisibly. Called for the picture it draws.
#' @seealso [contrast_energy()], [rdm()], [rsa()], [crossnobis()], and
#'   [rdm_sampling_covariance()] for the objects these methods draw, and
#'   [as.data.frame()] methods for the same values as a table.
#' @family geometry plans and views
#' @name plot_views
#' @importFrom graphics abline axis box image layout lcm legend mtext par
#'   plot.default points polygon rect segments text
#' @importFrom grDevices adjustcolor col2rgb colorRampPalette rgb
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
#'      highlight = planted, highlight_label = "planted signal")
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
                                      main = NULL, top = 200L, ...) {
  which <- match.arg(which)
  dots <- list(...)
  total <- as.numeric(x$total)
  coherent <- as.numeric(x$coherent)
  configuration <- as.numeric(x$configuration)
  signed <- .contrast_signed(x)
  highlight <- .resolve_highlight(highlight, x$index, length(total))
  highlight_label <- .validate_highlight_label(highlight_label)
  top <- .validate_top(top)
  dense <- length(total) > .dense_measurements

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  if (identical(which, "both")) {
    graphics::par(mfrow = c(1L, 2L))
  }
  graphics::par(mar = c(4.4, 4.4, 3.2, 1.1))

  if (which %in% c("both", "decomposition")) {
    marked <- logical(length(total))
    marked[highlight] <- TRUE
    # A line of constant total has slope -1 only when the two axes share a
    # scale, so the panel is drawn square with one common range and `asp = 1`.
    # The guides are then visually true rather than merely labeled as such.
    limits <- .finite_range(coherent, configuration)
    limits[[2L]] <- limits[[2L]] + 0.2 * diff(limits)
    xlim <- limits
    ylim <- limits
    defaults <- list(
      x = coherent, y = configuration,
      xlim = xlim, ylim = ylim, asp = 1, cex.main = 1,
      xlab = "Coherent energy (region-average contrast)",
      ylab = "Configuration energy (pattern beyond that mean)",
      main = .panel_main(main, 1L, "Coherent versus configuration energy")
    )
    if (dense) {
      # The cloud is drawn as counts below, so the panel opens empty.
      defaults$type <- "n"
    } else {
      scale <- .measurement_scale(length(total))
      # One meaning per color across the figure. Highlighting is carried by a
      # black ring, and the fill inside it carries the sign. When a highlight
      # set is given, the measurements outside it go plain gray -- the same
      # gray the profile panel already uses for "not highlighted" -- so no
      # color says two things. With no highlight set there is nothing for
      # gray to mean, and every point keeps its sign color.
      fill <- .signed_colors(signed$values)
      point_color <- if (length(highlight)) {
        ifelse(marked, .geometry_colors[["black"]],
          .fade(.geometry_colors[["neutral"]], 0.55 * scale$alpha))
      } else {
        .fade(fill, 0.7 * scale$alpha)
      }
      defaults$col <- point_color
      defaults$bg <- fill
      defaults$pch <- ifelse(marked, 21L, 19L)
      defaults$cex <- ifelse(marked, 1.05, 0.65) * scale$cex
    }
    merged <- .merge_graphics_args(defaults, dots)
    do.call(graphics::plot.default, merged)
    if (dense) {
      # Bin over the limits the panel was actually drawn with, so a caller's
      # `xlim` moves the cells with the points they stand for.
      .plot_density_bins(coherent, configuration, signed$values,
        .finite_range(merged$xlim), .finite_range(merged$ylim))
    }
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
    if (dense && length(highlight)) {
      # The density cells stand for the cloud; the highlight set is the
      # reader's question and is drawn measurement by measurement over it.
      ringed <- .measurement_scale(length(highlight))
      graphics::points(coherent[highlight], configuration[highlight],
        pch = 21L, cex = 1.05 * ringed$cex,
        col = .geometry_colors[["black"]],
        bg = .fade(.signed_colors(signed$values[highlight]), ringed$alpha))
    }
    # Four glyphs, four rows, and each color appears once. The highlight row
    # is drawn as the encoding itself -- an unfilled black ring -- because
    # the fill of a ringed point is its sign, and a gray filled swatch here
    # would collide with the gray that means "not highlighted".
    plain <- if (dense) 15L else 19L
    graphics::legend(
      "topright", bty = "n", cex = 0.75,
      legend = c(
        paste0(signed$label, " > 0"), paste0(signed$label, " < 0"),
        if (length(highlight) && !dense) "not highlighted",
        if (length(highlight)) {
          sprintf("%s (%d measurements)", highlight_label, length(highlight))
        }
      ),
      col = c(.geometry_colors[["blue"]], .geometry_colors[["vermillion"]],
        if (length(highlight) && !dense) .geometry_colors[["neutral"]],
        if (length(highlight)) .geometry_colors[["black"]]),
      pt.bg = NA,
      pch = c(plain, plain,
        if (length(highlight) && !dense) plain,
        if (length(highlight)) 21L),
      pt.cex = c(0.9, 0.9,
        if (length(highlight) && !dense) 0.65,
        if (length(highlight)) 1.05)
    )
    # The dashed guides no longer need a legend row of their own.
    graphics::mtext(
      if (dense) {
        sprintf("dashed lines: constant total energy; %s as density bins",
          .msg_count(length(total), "measurement"))
      } else {
        "dashed lines: constant total energy"
      },
      side = 3, line = 0.2, cex = 0.75, col = .geometry_colors[["neutral"]]
    )
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
      legend_label = highlight_label, top = top,
      # The same glyph as the decomposition panel: sign-colored fill, black
      # ring. Without it the profile's accent color would be vermillion,
      # which the other panel spends on "signed < 0".
      highlight_fill = if (length(highlight)) {
        .signed_colors(signed$values[highlight])
      }
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
  .check_flag(annotate, "annotate")
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
    col = colors, axes = FALSE, asp = 1, xlab = "", ylab = "", cex.main = 1,
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
                                 main = NULL, top = 200L, ...) {
  coefficients <- as.matrix(x$coefficients)
  labels <- colnames(coefficients)
  if (is.null(labels)) {
    labels <- paste0("coefficient", seq_len(ncol(coefficients)))
  }
  if (!is.null(terms)) {
    if (!is.character(terms) || !length(terms) || anyNA(terms) ||
        !all(terms %in% labels)) {
      .input_error("`terms` must name RSA coefficient columns of this view.")
    }
    coefficients <- coefficients[, terms, drop = FALSE]
    labels <- terms
  }
  highlight <- .resolve_highlight(highlight, x$index, nrow(coefficients))
  highlight_label <- .validate_highlight_label(highlight_label)
  top <- .validate_top(top)
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
      legend_label = if (panel == 1L) highlight_label else NULL,
      top = top
    )
  }
  invisible(x)
}

#' @rdname plot_views
#' @export
plot.effect_crossnobis_view <- function(x, highlight = NULL,
                                        highlight_label = "highlighted",
                                        main = NULL, top = 200L, ...) {
  values <- as.numeric(x$values)
  highlight <- .resolve_highlight(highlight, x$index, length(values))
  highlight_label <- .validate_highlight_label(highlight_label)
  top <- .validate_top(top)
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
    legend_label = highlight_label, top = top
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
      .contract_error(paste0(
        "The `estimate` view reports different distance pairs than this ",
        "sampling covariance."
      ))
    }
    position <- match(x$source$node, .measurement_ids(estimate$index))
    if (is.na(position)) {
      .contract_error(sprintf(paste0(
        "The `estimate` view does not contain measurement %s, which this ",
        "sampling covariance describes."
      ), x$source$node))
    }
    return(list(values = as.numeric(values[position, labels]),
      centered = TRUE))
  }
  if (!.is_finite_numeric(estimate) || length(estimate) != length(labels) ||
      anyNA(estimate)) {
    .input_error(sprintf(
      "`estimate` must be %d finite distances or an `effect_rdm_view`.",
      length(labels)
    ))
  }
  if (!is.null(names(estimate))) {
    if (!setequal(names(estimate), labels)) {
      .input_error(
        "Named `estimate` values must name every distance exactly once."
      )
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
  if (!.is_number(level) || level <= 0 || level >= 1) {
    .input_error("`level` must be one number strictly between 0 and 1.")
  }
  .check_flag(sort, "sort")
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
      sprintf("Signed squared distance with %g %% sampling interval",
        100 * level)
    } else {
      sprintf("%g %% sampling interval (no point estimate supplied)",
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
