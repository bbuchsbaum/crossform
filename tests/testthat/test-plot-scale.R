# A contrast view over tens of thousands of measurements must still be a
# picture. Above `.dense_measurements` the profile panel summarizes its bulk
# as a band and the decomposition panel as counted cells; at and below the
# sizes the README and the vignettes draw, nothing changes at all.

# `pkgload::load_all()` registers S3 methods from NAMESPACE, so before
# `roxygen2::roxygenise()` regenerates it these methods exist but do not
# dispatch. Registering them here keeps the tests exercising `plot()` itself.
local({
  registered <- get(".__S3MethodsTable__.", envir = asNamespace("base"))
  if (!exists("plot.effect_contrast_view", envir = registered,
      inherits = FALSE)) {
    registerS3method("plot", "effect_contrast_view",
      utils::getFromNamespace("plot.effect_contrast_view", "crossform"),
      envir = asNamespace("base"))
  }
})

with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

# A synthetic view rather than a compiled one: these branches are chosen by
# the number of measurements and by nothing else, and 20000 real searchlights
# would cost minutes to build for a picture that reads the same random
# numbers. `.new_effect_contrast_view()` is the constructor the executor uses
# (R/views.R), so the object under test is a real `effect_contrast_view` with
# every field the plot methods read; only `$receipt` is empty, and no plot
# method reads it.
synthetic_contrast_view <- function(n, planted = integer(), seed = 20000L) {
  set.seed(seed)
  effects <- c("face", "house", "body", "tool")
  weights <- stats::setNames(c(1, -1, 0, 0), effects)
  marginals <- structure(
    list(endpoint = matrix(stats::rnorm(n * length(effects), 0, 0.3), n,
      length(effects), dimnames = list(NULL, effects))),
    semantics = "undirected_endpoint"
  )
  coherent <- stats::rnorm(n, 0.05, 0.12)
  configuration <- stats::rnorm(n, 0.08, 0.20)
  # A planted set the summary must not swallow, and a signed estimate that
  # falls below zero, which the band and the needles both have to survive.
  coherent[planted] <- coherent[planted] + 0.5
  configuration[planted] <- configuration[planted] + 0.9
  constructor <- utils::getFromNamespace(".new_effect_contrast_view",
    "crossform")
  constructor(coherent + configuration, coherent, marginals, weights,
    as.character(seq_len(n)), NULL)
}

test_that("glyph scaling is a no-op at the sizes the figures draw", {
  scale <- utils::getFromNamespace(".measurement_scale", "crossform")
  fade <- utils::getFromNamespace(".fade", "crossform")
  # The README fixture is 280 measurements; nothing about it may move.
  for (n in c(1L, 100L, 280L, 500L)) {
    expect_identical(scale(n), list(cex = 1, alpha = 1))
  }
  expect_identical(fade("#C8C8C8", 1), "#C8C8C8")
  # Beyond that the mid-size panels thin out, monotonically.
  sizes <- c(1000L, 2000L, 4000L, 5000L)
  factors <- vapply(sizes, function(n) scale(n)$cex, numeric(1))
  expect_true(all(diff(factors) < 0))
  expect_true(all(factors < 1 & factors > 0.45))
  alphas <- vapply(sizes, function(n) scale(n)$alpha, numeric(1))
  expect_true(all(diff(alphas) < 0))
  expect_true(all(alphas < 1 & alphas > 0.25))
  # Both factors have a floor, so a very crowded panel still shows something.
  expect_identical(scale(500000L), list(cex = 0.45, alpha = 0.25))
  expect_gte(scale(.Machine$integer.max)$alpha, 0.25)
  # A faded color is the color `adjustcolor()` would have produced, so the
  # scaled panels stay on the declared palette.
  expect_identical(fade("#C8C8C8", 0.55),
    grDevices::adjustcolor("#C8C8C8", alpha.f = 0.55))
  expect_length(fade(c("#0072B2", "#D55E00"), c(0.2, 0.9)), 2L)
})

test_that("a 20000-measurement contrast view summarises both panels", {
  planted <- seq(1L, 20000L, by = 200L)
  view <- synthetic_contrast_view(20000L, planted)
  expect_length(view$total, 20000L)

  bins <- 0L
  original <- crossform:::.plot_density_bins
  testthat::local_mocked_bindings(
    .plot_density_bins = function(...) {
      bins <<- bins + 1L
      original(...)
    },
    .package = "crossform"
  )
  with_null_device({
    drawn <- withVisible(plot(view, highlight = planted,
      highlight_label = "planted signal"))
    expect_false(drawn$visible)
    expect_identical(drawn$value, view)
    # The decomposition panel went through the binning path, once.
    expect_identical(bins, 1L)
    # ... and only for the panel that draws it.
    expect_silent(plot(view, which = "profile"))
    expect_identical(bins, 1L)
    expect_silent(plot(view, which = "decomposition"))
    expect_identical(bins, 2L)
  })
})

test_that("the crowded profile keeps a band and only the top needles", {
  profile <- utils::getFromNamespace(".plot_measurement_profile", "crossform")
  view <- synthetic_contrast_view(20000L)
  draw <- function(values, top = 200L) {
    profile(values, integer(), xlab = "", ylab = "", main = "",
      color = "#D55E00", top = top)
  }
  with_null_device({
    summary <- draw(view$total)
    expect_true(summary$dense)
    expect_length(summary$drawn, 200L)
    # `top` chooses how many, and caps at the number of measurements.
    expect_length(draw(view$total, top = 25L)$drawn, 25L)
    expect_length(draw(view$total, top = 50000L)$drawn, 20000L)
    # Below the threshold every measurement is still drawn individually.
    small <- draw(view$total[seq_len(280L)])
    expect_false(small$dense)
    expect_identical(small$drawn, seq_len(280L))
  })
})

test_that("the top needles are the extremes on both sides of zero", {
  band <- utils::getFromNamespace(".plot_profile_band", "crossform")
  # Two large positives, one large negative, and 5000 measurements of noise
  # in between. A summary that ranked by value alone would lose the negative
  # one, and these estimates are crossvalidated: the negative half is data.
  values <- c(rep(c(0.01, -0.01), length.out = 5000L), 9, -8, 7)
  with_null_device({
    graphics::plot.default(seq_along(values), values, type = "n")
    kept <- band(values, reference = 0, top = 3L)
  })
  expect_identical(sort(kept), c(5001L, 5002L, 5003L))
  expect_identical(values[kept], c(9, -8, 7))
})

test_that("the density panel bins rather than drawing every measurement", {
  bins <- utils::getFromNamespace(".plot_density_bins", "crossform")
  view <- synthetic_contrast_view(20000L)
  with_null_device({
    graphics::plot.default(view$coherent, view$configuration, type = "n")
    filled <- bins(view$coherent, view$configuration, view$signed,
      range(view$coherent), range(view$configuration))
  })
  # Cells, not measurements: far fewer marks than there are values, and never
  # more than the grid holds.
  expect_gt(filled, 0L)
  expect_lt(filled, 64L * 64L)
  expect_lt(filled, length(view$total) / 4)
})

test_that("`top` is refused rather than guessed", {
  small <- synthetic_contrast_view(120L)
  large <- synthetic_contrast_view(20000L)
  with_null_device({
    for (view in list(small, large)) {
      for (bad in list(0L, -5L, 2.5, NA_integer_, Inf, "200", NULL,
          c(10L, 20L))) {
        expect_error(plot(view, top = bad), "`top`",
          class = "effect_input_error")
      }
      # The check does not wait for a view large enough to use the argument.
      expect_silent(plot(view, which = "profile", top = 12L))
    }
  })
})
