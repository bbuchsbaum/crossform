if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

test_that("public docs no longer deny the implemented population surface", {
  repo <- testthat::test_path("..", "..")
  files <- c("README.md", "R/crossform-package.R",
    "vignettes/introduction.Rmd", "vignettes/interpreting-results.Rmd",
    "vignettes/population-form.Rmd")
  text <- paste(unlist(lapply(file.path(repo, files), readLines, warn = FALSE)),
                collapse = "\n")
  stale <- c("does no spatial and no group inference",
             "performs no spatial and no group inference",
             "Within-participant results only. There is no group-inference path")
  for (phrase in stale) expect_false(grepl(phrase, text, fixed = TRUE))
  for (term in c("population_uncertainty()", "classical OLS or HC3",
    "population_wild_bootstrap()", "pointwise",
    "conditional on that transport", "exact cellwise coverage",
    "Simultaneous", "cross-node covariance", "informative coverage")) {
    expect_match(text, term, fixed = TRUE)
  }
})

test_that("documented population methods match live signatures", {
  expect_true("estimator" %in% names(formals(population_uncertainty)))
  expect_identical(eval(formals(population_uncertainty)$estimator),
                   c("classical", "HC3"))
  expect_true(all(c("replicates", "seed", "weights", "level") %in%
                    names(formals(population_wild_bootstrap))))
  expect_true(all(c("minimum_coverage", "minimum_transport_quality",
                    "material_change") %in%
                    names(formals(population_diagnostics))))
  expect_true("interval" %in% names(formals(population_scale_profile)))
})
