if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

population_interpretability_path <- function(name) testthat::test_path(
  "..", "..", "inst", "extdata", "certification", name
)

test_that("matched population interpretability gates pass with limits retained", {
  raw <- utils::read.csv(population_interpretability_path(
    "population-interpretability-replicates.csv"), stringsAsFactors = FALSE)
  results <- utils::read.csv(population_interpretability_path(
    "population-interpretability-results.csv"), stringsAsFactors = FALSE)
  verdicts <- utils::read.csv(population_interpretability_path(
    "population-interpretability-verdicts.csv"), stringsAsFactors = FALSE)

  expect_true(all(verdicts$passes))
  expect_setequal(unique(raw$coverage_regime),
    c("supported_complete", "failure_informative", "influence_supported"))
  expect_setequal(unique(results$component),
                  c("total", "coherent", "configuration"))
  expect_setequal(unique(results$term), c("(Intercept)", "covariate_x"))
  expect_true(all(results$marginal_claim_supported[
    results$coverage_regime == "supported_complete"
  ]))
  expect_true(all(!results$marginal_claim_supported[
    results$coverage_regime == "failure_informative"
  ]))
  expect_false(any(grepl("ratio", names(raw), ignore.case = TRUE)))
})

test_that("population interpretability artifacts remain source bound", {
  manifest <- utils::read.csv(population_interpretability_path(
    "population-interpretability-checksums.csv"), stringsAsFactors = FALSE)
  files <- file.path(testthat::test_path("..", ".."), manifest$path)
  expect_true(all(file.exists(files)))
  expect_identical(unname(tools::md5sum(files)), manifest$digest)
})
