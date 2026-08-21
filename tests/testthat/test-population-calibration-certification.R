if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

population_calibration_path <- function(name) testthat::test_path(
  "..", "..", "inst", "extdata", "certification", name
)

test_that("population calibration artifacts pass their prespecified gates", {
  results <- utils::read.csv(population_calibration_path(
    "population-calibration-results.csv"), stringsAsFactors = FALSE)
  verdicts <- utils::read.csv(population_calibration_path(
    "population-calibration-verdicts.csv"), stringsAsFactors = FALSE)
  replicates <- utils::read.csv(population_calibration_path(
    "population-calibration-replicates.csv"), stringsAsFactors = FALSE)

  expect_true(all(verdicts$passes))
  expect_setequal(unique(results$method),
                  c("classical", "HC3", "wild_bootstrap"))
  expect_setequal(unique(results$transport_regime), c("fixed", "cross_fitted"))
  expect_setequal(unique(results$coverage_regime),
                  c("complete", "independent", "informative"))
  expect_true(all(table(replicates$dataset_id) == 3L))
  expect_identical(nrow(replicates), 8L * 500L * 3L)
  informative <- results$coverage_regime == "informative"
  expect_true(all(!results$marginal_claim_supported[informative]))
})

test_that("population calibration checksums bind sources and artifacts", {
  manifest <- utils::read.csv(population_calibration_path(
    "population-calibration-checksums.csv"), stringsAsFactors = FALSE)
  files <- file.path(testthat::test_path("..", ".."), manifest$path)
  expect_true(all(file.exists(files)))
  expect_identical(unname(tools::md5sum(files)), manifest$digest)
})

test_that("a smoke court is paired and exercises informative refusal", {
  environment <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "benchmarks",
    "population-calibration", "00-simulate.R"), envir = environment)
  smoke <- environment$population_calibration_replicates(
    replications = 4L, bootstrap_replicates = 39L, seed = 9901L
  )
  expect_identical(nrow(smoke), 8L * 4L * 3L)
  expect_true(all(table(smoke$dataset_id) == 3L))
  expect_true(any(smoke$coverage_regime == "informative" &
                    !smoke$marginal_claim_supported))
})
