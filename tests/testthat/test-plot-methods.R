# Every picture must draw from a single call with no arguments, return its
# input invisibly, and refuse a measurement or highlight it cannot honour.

# `pkgload::load_all()` registers S3 methods from NAMESPACE, so before
# `roxygen2::roxygenise()` regenerates it these methods exist but do not
# dispatch. Registering them here keeps the tests exercising `plot()` itself,
# which is the user's path, rather than the unexported functions.
local({
  registered <- get(".__S3MethodsTable__.", envir = asNamespace("base"))
  classes <- c("effect_contrast_view", "effect_rdm_view", "effect_rsa_view",
    "effect_crossnobis_view", "effect_rdm_sampling_covariance")
  for (class in classes) {
    method <- paste0("plot.", class)
    if (exists(method, envir = registered, inherits = FALSE)) next
    registerS3method("plot", class,
      utils::getFromNamespace(method, "crossform"),
      envir = asNamespace("base"))
  }
})

plot_fixture <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    # A smaller volume than the default fixture: 100 searchlights still carry
    # planted and unplanted measurements, negative estimates, and missing
    # coherence fractions, which is everything these pictures must handle.
    example <- example_fmri_effects(dimensions = c(5L, 5L, 4L))
    plan <- plan_geometry(
      example$fit$relation, example$frame,
      cross_partitions(
        example$fit$relation,
        independence = "independent", generalizes_over = "run"
      )
    )
    energy <- contrast_energy(plan, example$contrast)
    distances <- rdm(plan)
    peak <- which.max(energy$total)
    cached <<- list(
      example = example,
      plan = plan,
      energy = energy,
      distances = distances,
      coefficients = rsa(plan, models = list(category = example$model_rdm)),
      peak = peak,
      planted = example$truth$signal_measurements,
      uncertainty = rdm_sampling_covariance(
        plan, example$fit, target = "null", at = peak
      )
    )
    cached
  }
})

with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

expect_drawn <- function(object, view) {
  drawn <- expect_invisible(object)
  expect_identical(drawn, view)
}

test_that("the contrast view draws its decomposition and profile", {
  fixture <- plot_fixture()
  view <- fixture$energy
  with_null_device({
    expect_drawn(plot(view), view)
    expect_drawn(plot(view, highlight = fixture$planted), view)
    expect_drawn(plot(view, which = "decomposition"), view)
    expect_drawn(plot(view, which = "profile",
      highlight = fixture$planted), view)
    # A logical mask and measurement identifiers select the same set.
    mask <- logical(length(view$total))
    mask[fixture$planted] <- TRUE
    expect_silent(plot(view, which = "profile", highlight = mask))
    expect_silent(plot(view, which = "profile",
      highlight = as.character(fixture$planted)))
    # `...` reaches the underlying call.
    expect_silent(plot(view, which = "decomposition", pch = 3L,
      col = "#0072B2", main = "custom"))
  })
})

test_that("a highlight set can be named in the legend", {
  fixture <- plot_fixture()
  with_null_device({
    expect_drawn(
      plot(fixture$energy, highlight = fixture$planted,
        highlight_label = "planted signal"),
      fixture$energy
    )
    expect_silent(plot(fixture$energy, which = "profile",
      highlight = fixture$planted, highlight_label = "planted signal"))
    expect_drawn(
      plot(fixture$coefficients, highlight = fixture$planted,
        highlight_label = "planted signal"),
      fixture$coefficients
    )
    # The label is checked even without a highlight set, so a typo surfaces
    # rather than sitting unused until a highlight is passed.
    expect_silent(plot(fixture$energy, highlight_label = "planted signal"))
    for (bad in list(NA_character_, "", character(), c("a", "b"), 1)) {
      expect_error(
        plot(fixture$energy, highlight = fixture$planted,
          highlight_label = bad),
        "one non-empty character string"
      )
      expect_error(
        plot(fixture$coefficients, highlight_label = bad),
        "one non-empty character string"
      )
    }
  })
})

test_that("the contrast view refuses a highlight it cannot honour", {
  fixture <- plot_fixture()
  view <- fixture$energy
  with_null_device({
    expect_error(plot(view, highlight = 0L), "measurement positions")
    expect_error(plot(view, highlight = length(view$total) + 1L),
      "measurement positions")
    expect_error(plot(view, highlight = c(TRUE, FALSE)), "logical")
    expect_error(plot(view, highlight = "no-such-measurement"),
      "not in this view")
    expect_error(plot(view, which = "spectrum"), "should be one of")
  })
})

test_that("the contrast decomposition survives non-finite coherence", {
  fixture <- plot_fixture()
  view <- fixture$energy
  # `coherence_fraction` is NA wherever the raw components are not a
  # nonnegative partition, which is the common case; the picture must not
  # depend on it.
  expect_true(anyNA(view$coherence_fraction))
  broken <- view
  broken$coherence_fraction[] <- NA_real_
  with_null_device(expect_silent(plot(broken)))
})

test_that("the contrast view draws directed endpoint marginals", {
  fixture <- plot_fixture()
  # A coupling plan reports `signed` as one marginal per endpoint role rather
  # than a single vector; the decomposition panel colors by the first of them
  # and names which one it drew.
  directed <- fixture$energy
  directed$signed <- list(
    left = fixture$energy$signed, right = -fixture$energy$signed
  )
  with_null_device({
    expect_drawn(plot(directed), directed)
    expect_drawn(plot(directed, which = "decomposition"), directed)
  })
})

test_that("the RDM view draws a symmetric heatmap at a chosen measurement", {
  fixture <- plot_fixture()
  view <- fixture$distances
  with_null_device({
    expect_drawn(plot(view), view)
    expect_drawn(plot(view, measurement = fixture$peak), view)
    expect_drawn(plot(view, measurement = "mean"), view)
    expect_silent(plot(view, measurement = as.character(fixture$peak)))
    expect_silent(plot(view, measurement = fixture$peak, annotate = FALSE))
    expect_silent(plot(view, measurement = 1L, main = "custom"))
  })
})

test_that("the RDM view refuses a measurement outside the view", {
  fixture <- plot_fixture()
  view <- fixture$distances
  with_null_device({
    expect_error(plot(view, measurement = nrow(view$values) + 1L),
      "measurement position between")
    expect_error(plot(view, measurement = 0L), "measurement position between")
    expect_error(plot(view, measurement = 1.5), "measurement position between")
    expect_error(plot(view, measurement = "not-a-measurement"),
      "does not identify a measurement")
    expect_error(plot(view, measurement = c(1L, 2L)),
      "measurement position between")
    expect_error(plot(view, annotate = NA), "TRUE or FALSE")
  })
})

test_that("the RDM heatmap reconstructs the packed pair columns", {
  fixture <- plot_fixture()
  view <- fixture$distances
  rebuild <- utils::getFromNamespace(".rdm_matrix", "crossform")
  matrix_values <- rebuild(view, view$values[fixture$peak, ])
  expect_identical(dim(matrix_values), c(4L, 4L))
  expect_identical(rownames(matrix_values), colnames(matrix_values))
  expect_equal(diag(matrix_values), c(face = 0, body = 0, house = 0, tool = 0))
  expect_equal(matrix_values, t(matrix_values))
  expect_equal(
    matrix_values[view$pairs$left[[2L]], view$pairs$right[[2L]]],
    unname(view$values[fixture$peak, 2L])
  )
})

test_that("the RSA view draws one panel per coefficient", {
  fixture <- plot_fixture()
  view <- fixture$coefficients
  with_null_device({
    expect_drawn(plot(view), view)
    expect_drawn(plot(view, highlight = fixture$planted), view)
    expect_silent(plot(view, terms = "category"))
    expect_silent(plot(view, terms = colnames(view$coefficients),
      main = c("one", "two")))
    expect_error(plot(view, terms = "no-such-model"),
      "RSA coefficient columns")
    expect_error(plot(view, terms = character()), "RSA coefficient columns")
  })
})

test_that("the crossnobis view draws a signed index plot", {
  example <- plot_fixture()$example
  metric <- noise_precision(diag(example$domain$n_features), example$domain)
  plan <- plan_geometry(
    example$fit$relation, example$frame,
    cross_partitions(example$fit$relation, independence = "independent"),
    metric = metric
  )
  view <- crossnobis(plan, example$contrast)
  with_null_device({
    expect_drawn(plot(view), view)
    expect_drawn(plot(view, highlight = example$truth$signal_measurements),
      view)
    expect_silent(plot(view, main = "custom", col = "#009E73"))
    expect_drawn(
      plot(view, highlight = example$truth$signal_measurements,
        highlight_label = "planted signal"),
      view
    )
    expect_error(plot(view, highlight = length(view$values) + 1L),
      "measurement positions")
    expect_error(plot(view, highlight_label = NA_character_),
      "one non-empty character string")
  })
  # Signed estimates are drawn where they fall, so the panel must span them.
  expect_true(any(view$values < 0))
})

test_that("the sampling covariance draws intervals for one measurement", {
  fixture <- plot_fixture()
  view <- fixture$uncertainty
  with_null_device({
    expect_drawn(plot(view), view)
    expect_drawn(plot(view, estimate = fixture$distances), view)
    expect_silent(plot(view, estimate = fixture$distances, sort = FALSE))
    expect_silent(plot(view, level = 0.9, main = "custom"))
    named <- stats::setNames(seq_along(view$labels), view$labels)
    expect_silent(plot(view, estimate = rev(named)))
    expect_silent(plot(view, estimate = as.numeric(named)))
  })
})

test_that("the sampling covariance refuses mismatched inputs", {
  fixture <- plot_fixture()
  view <- fixture$uncertainty
  with_null_device({
    expect_error(plot(view, level = 0), "strictly between 0 and 1")
    expect_error(plot(view, level = 1), "strictly between 0 and 1")
    expect_error(plot(view, sort = NA), "TRUE or FALSE")
    expect_error(plot(view, estimate = 1), "finite distances")
    expect_error(
      plot(view, estimate = stats::setNames(
        rep(0, length(view$labels)), rev(seq_along(view$labels))
      )),
      "name every distance exactly once"
    )
    # The estimate must come from the measurement this covariance describes,
    # not merely from a view whose pair labels happen to agree.
    at_peak <- fixture$distances
    at_peak$index <- at_peak$index[fixture$peak]
    at_peak$values <- at_peak$values[fixture$peak, , drop = FALSE]
    expect_silent(plot(view, estimate = at_peak))
    elsewhere <- fixture$distances
    keep <- setdiff(seq_len(nrow(fixture$distances$values)), fixture$peak)
    elsewhere$index <- elsewhere$index[keep]
    elsewhere$values <- elsewhere$values[keep, , drop = FALSE]
    expect_error(plot(view, estimate = elsewhere),
      "does not contain measurement")
    mislabelled <- at_peak
    colnames(mislabelled$values) <- paste0("q", seq_along(view$labels))
    expect_error(plot(view, estimate = mislabelled),
      "different distance pairs")
  })
})

test_that("plotting restores the graphical parameters it changed", {
  fixture <- plot_fixture()
  with_null_device({
    before <- graphics::par(c("mfrow", "mar", "mfg"))
    plot(fixture$energy, highlight = fixture$planted)
    plot(fixture$distances)
    plot(fixture$coefficients)
    plot(fixture$uncertainty)
    expect_identical(graphics::par(c("mfrow", "mar", "mfg")), before)
  })
})
