if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

test_that("README and package help share the same bounded core claim", {
  repo <- testthat::test_path("..", "..")
  read <- function(path) paste(readLines(file.path(repo, path), warn = FALSE),
                               collapse = "\n")
  readme <- read("README.md")
  package <- read("R/crossform-package.R")
  description <- read("DESCRIPTION")
  phrases <- c(
    "univariate contrasts, multivariate distances and fixed linear",
    "one declared crossvalidated bilinear geometry",
    "coherent and",
    "configurational estimand components"
  )
  for (phrase in phrases) {
    expect_match(readme, phrase, fixed = TRUE)
    expect_match(package, phrase, fixed = TRUE)
  }
  expect_match(description, "crossvalidated bilinear geometry", fixed = TRUE)
  expect_match(readme, "not a complete fMRI analysis solution", fixed = TRUE)
  expect_match(package, "not separate biological mechanisms", fixed = TRUE)
  expect_match(readme, "every RSA statistic", fixed = TRUE)
  expect_match(package, "No simultaneous, transport-marginal", fixed = TRUE)
})
