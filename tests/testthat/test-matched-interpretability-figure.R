if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

canonical_figure_environment <- function() {
  environment <- new.env(parent = globalenv())
  for (file in c(
    "00-mixture-generator.R", "01-multiscale-scenarios.R",
    "02-paired-observations.R", "03-conventional-baselines.R",
    "04-canonical-figure.R"
  )) {
    sys.source(testthat::test_path(
      "..", "..", "benchmarks", "matched-interpretability", file
    ), envir = environment)
  }
  environment
}

canonical_figure_artifact <- function(name) {
  testthat::test_path("..", "..", "inst", "extdata", "certification", name)
}

test_that("canonical data contains every panel scenario scale and interval", {
  path <- canonical_figure_artifact("matched-interpretability-figure-data.csv")
  skip_if_not(file.exists(path), "canonical figure data has not been generated")
  data <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_setequal(
    data$metric,
    c("ground_truth_pattern", "regional_activation",
      "aggregate_crossvalidated_magnitude", "total", "coherent",
      "configuration", "coherence_fraction")
  )
  expect_setequal(data$scenario,
                  c("broad_coherent", "mixed_broad_fine",
                    "fine_configuration"))
  spectrum <- data[data$metric == "coherence_fraction", ]
  expect_setequal(spectrum$scale, c(0.01, 1.01, 2.01, 4.01))
  estimates <- data[is.finite(data$estimate_mean), ]
  expect_true(all(estimates$n_replications == 24L))
  expect_true(all(estimates$n_valid > 0L))
  expect_true(all(estimates$interval_lower <= estimates$estimate_mean))
  expect_true(all(estimates$estimate_mean <= estimates$interval_upper))
})

test_that("equal totals and organization ordering are numerical facts", {
  data <- utils::read.csv(canonical_figure_artifact(
    "matched-interpretability-figure-data.csv"
  ), stringsAsFactors = FALSE)
  aggregate <- data[data$metric == "aggregate_crossvalidated_magnitude", ]
  expect_lt(diff(range(aggregate$truth)), 1e-12)
  widest <- data[data$metric == "coherence_fraction" & data$scale == 4.01, ]
  widest <- widest[match(c("broad_coherent", "mixed_broad_fine",
                           "fine_configuration"), widest$scenario), ]
  expect_true(all(diff(widest$truth) < 0))
  expect_true(all(diff(widest$estimate_mean) < 0))
  expect_gt(widest$truth[[1L]] - widest$truth[[3L]], 0.9)

  coherent <- data[data$metric == "coherent" & data$scale == 4.01, ]
  configuration <- data[
    data$metric == "configuration" & data$scale == 4.01, ]
  coherent <- coherent[match(widest$scenario, coherent$scenario), ]
  configuration <- configuration[
    match(widest$scenario, configuration$scenario), ]
  expect_true(all(diff(coherent$truth) < 0))
  expect_true(all(diff(configuration$truth) > 0))
})

test_that("the PNG render has six governed panels and stable dimensions", {
  generator <- canonical_figure_environment()
  data <- utils::read.csv(canonical_figure_artifact(
    "matched-interpretability-figure-data.csv"
  ), stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".png")
  on.exit(unlink(path), add = TRUE)
  contract <- generator$render_canonical_figure(data, path)

  expect_identical(length(contract$panels), 6L)
  expect_identical(contract$dimensions,
                   c(width = 1800L, height = 1200L, resolution = 150L))
  expect_identical(contract$component_ylim,
                   range(contract$component_ylim))
  expect_setequal(contract$scenarios,
                  c("broad_coherent", "mixed_broad_fine",
                    "fine_configuration"))
  expect_true(file.exists(path))
  expect_gt(file.info(path)$size, 30000)

  header <- readBin(path, what = "raw", n = 24L)
  expect_identical(as.integer(header[1:8]),
                   c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L))
  big_endian <- function(bytes) sum(as.integer(bytes) * 256^(3:0))
  expect_identical(big_endian(header[17:20]), 1800)
  expect_identical(big_endian(header[21:24]), 1200)
})

test_that("the article caption states the estimand and ambiguity criterion", {
  path <- testthat::test_path("..", "..", "vignettes",
                             "matched-interpretability.Rmd")
  article <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(article, "average squared-Euclidean", fixed = TRUE)
  expect_match(article, "differs by at most `1e-12`", fixed = TRUE)
  expect_match(article, "differ by at least `0.2`", fixed = TRUE)
  expect_match(article, "D-E: coherent and configuration magnitude", fixed = TRUE)
  expect_match(article, "matched-simulation evidence", fixed = TRUE)
  expect_match(article, "Monte Carlo variation rather than population coverage",
               fixed = TRUE)
})
