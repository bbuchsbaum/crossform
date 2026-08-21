if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

test_that("the terminology contract fixes components instruments and units", {
  path <- testthat::test_path("..", "..", "design", "terminology.md")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  for (term in c("crossvalidated bilinear geometry", "coherent component",
    "configurational component", "coherent share", "point scale",
    "radius N units", "detection map", "attribution ledger",
    "population coefficient", "ontological categories",
    "signed response units", "squared response units", "dimensionless")) {
    expect_match(text, term, fixed = TRUE)
  }
})

test_that("headline prose avoids biological mechanism reification", {
  repo <- testthat::test_path("..", "..")
  files <- c("README.md", "R/crossform-package.R",
    "vignettes/matched-interpretability.Rmd",
    "vignettes/population-form.Rmd")
  text <- paste(unlist(lapply(file.path(repo, files), readLines, warn = FALSE)),
                collapse = "\n")
  forbidden <- c("coherent neural mechanism", "configurational neural mechanism",
                 "coherent brain process", "configuration network",
                 "average of participant ratios")
  for (phrase in forbidden) {
    expect_false(grepl(tolower(phrase), tolower(text), fixed = TRUE),
                 info = phrase)
  }
  expect_match(text, "not separate biological mechanisms", fixed = TRUE)
  expect_match(text, "participant_ratio_average", fixed = TRUE)
})

test_that("public component captions state roles and squared units", {
  plot_source <- paste(readLines(testthat::test_path(
    "..", "..", "R", "plot-methods.R"), warn = FALSE), collapse = "\n")
  expect_match(plot_source, "Coherent component (squared response units)",
               fixed = TRUE)
  expect_match(plot_source, "Configurational component (squared response units)",
               fixed = TRUE)
})
