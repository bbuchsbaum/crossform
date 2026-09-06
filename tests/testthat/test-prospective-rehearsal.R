if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

test_that("the frozen prospective configuration rehearses end to end", {
  output <- withr::local_tempdir()
  script <- testthat::test_path("..", "..", "protocols", "prospective",
                               "rehearsal", "run.R")
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c(shQuote(script), shQuote(output)))
  expect_identical(status, 0L)
  config <- jsonlite::read_json(testthat::test_path(
    "..", "..", "protocols", "prospective", "discovery-v1.json"),
    simplifyVector = TRUE)
  expect_true(all(file.exists(file.path(output, config$outputs))))
  manifest <- utils::read.csv(file.path(output, "execution-manifest.csv"),
                              stringsAsFactors = FALSE)
  hashed <- manifest$hash_algorithm == "md5"
  expect_identical(unname(tools::md5sum(file.path(output,
    manifest$path[hashed]))), manifest$digest[hashed])
  expect_true(all(!manifest$real_data_claim))
  failures <- utils::read.csv(file.path(output, "failure-states.csv"),
                              stringsAsFactors = FALSE)
  expect_identical(failures$observed, failures$expected)
  expect_setequal(failures$observed, c("success", "missing_partition",
    "coverage_below_floor", "transport_quality_below_floor"))
  primary <- utils::read.csv(file.path(output, "primary-results.csv"),
                            stringsAsFactors = FALSE)
  expect_true(all(!primary$real_data_claim))
  expect_setequal(unique(primary$component),
                  c("total", "coherent", "configuration"))
})
