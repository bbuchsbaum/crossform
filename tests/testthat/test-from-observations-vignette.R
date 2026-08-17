vignette_source_path <- function(name) {
  # `R CMD build` keeps the vignette source in `inst/doc`, so the installed
  # package carries it under check. A source checkout has neither, and reads
  # `vignettes/` directly.
  installed <- tryCatch(
    system.file("doc", name, package = "crossform"),
    error = function(condition) ""
  )
  candidates <- c(
    if (nzchar(installed)) installed,
    testthat::test_path("..", "..", "vignettes", name),
    testthat::test_path("..", "..", "inst", "doc", name)
  )
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(candidate)
    }
  }
  NA_character_
}

test_that("the first-moment guide preserves the public identity ladder", {
  path <- vignette_source_path("from-observations.Rmd")
  if (is.na(path)) {
    skip(paste(
      "The from-observations vignette source is not available here:",
      "neither system.file('doc') nor the checkout's vignettes/ has it."
    ))
  }
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(text, "plan_relation\\(", fixed = FALSE)
  expect_match(text, "relation_plan_receipts\\(", fixed = FALSE)
  expect_match(text, "fit\\$signature", fixed = FALSE)
  expect_match(text, "estimate_relation\\(relation_request\\)", fixed = FALSE)
  expect_match(text, "plan_geometry\\(", fixed = FALSE)
  expect_match(text, 'independence = "independent"', fixed = TRUE)
  expect_match(text, 'generalizes_over = "run"', fixed = TRUE)
  expect_match(text, "raw_design_model\\(\\)", fixed = FALSE)
  expect_match(text, "BIDS-shaped", fixed = TRUE)
})
