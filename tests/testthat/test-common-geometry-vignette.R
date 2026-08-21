if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

test_that("the theorem article exposes executable evidence and bounded claims", {
  path <- testthat::test_path("..", "..", "vignettes",
                             "common-geometry-equivalence.Rmd")
  article <- paste(readLines(path, warn = FALSE), collapse = "\n")

  for (heading in c(
    "## The theorem",
    "## A hand-sized exact case",
    "## Randomized production-to-oracle evidence",
    "## External differential evidence",
    "## Estimator claim table",
    "## Limits and evidence status"
  )) {
    expect_match(article, heading, fixed = TRUE)
  }
  for (executable in c(
    "oracle_energy <- function",
    "crossnobis(plan, contrast)",
    "rdm(plan)",
    "rsa(plan, models",
    "vapply(court_seeds, run_case",
    "stopifnot(all(randomized_summary$passes))",
    "knitr::kable(claim_table"
  )) {
    expect_match(article, executable, fixed = TRUE)
  }
  expect_match(article, "unification of the named fixed estimand", fixed = TRUE)
  expect_match(article, "does not mean that Crossform implements every",
               fixed = TRUE)
  expect_match(article, "design/oracles/common-geometry-equivalence.R",
               fixed = TRUE)
  expect_match(article, "common-geometry-external-parity.csv", fixed = TRUE)
})

test_that("the package parity receipt is generated from the source exemplar", {
  source_path <- testthat::test_path(
    "..", "..", "exemplars", "rsatoolbox-parity", "results", "agreement.csv"
  )
  receipt_path <- testthat::test_path(
    "..", "..", "inst", "extdata", "certification",
    "common-geometry-external-parity.csv"
  )
  skip_if_not(file.exists(source_path) && file.exists(receipt_path),
              "source parity evidence is unavailable")
  source <- utils::read.csv(source_path, stringsAsFactors = FALSE)
  receipt <- utils::read.csv(receipt_path, stringsAsFactors = FALSE)

  expect_identical(
    receipt[names(source)],
    source
  )
  expect_true(all(receipt$fixture_id == "rsatoolbox-fixed-linear-v1"))
  expect_true(all(receipt$rsatoolbox_version == "0.3.2"))
  expect_true(all(receipt$passes))
  expect_lte(max(receipt$max_abs_diff), 1e-12)
})
