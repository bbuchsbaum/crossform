if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

review_bundle_paths <- function() {
  repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  list(
    repo = repo,
    builder = file.path(repo, "exemplars", "review-bundle", "build-bundle.R"),
    readme = file.path(repo, "exemplars", "review-bundle", "README.md")
  )
}

test_that("review navigation keeps evidence classes and refusals distinct", {
  paths <- review_bundle_paths()
  text <- paste(readLines(paths$readme, warn = FALSE), collapse = "\n")
  for (label in c(
    "THEOREM", "INTERNAL ORACLE", "EXTERNAL PARITY", "MATCHED SIMULATION",
    "RETROSPECTIVE ILLUSTRATION", "PROSPECTIVE PROTOCOL",
    "COMPLETED REAL DATA", "INDEPENDENT REPLICATION"
  )) {
    expect_match(text, label, fixed = TRUE)
  }
  expect_match(text, "readiness is `BLOCKED`", fixed = TRUE)
  expect_match(text, "No completed prospective real-data result", fixed = TRUE)
  expect_match(text, "informative coverage", fixed = TRUE)
  expect_match(text, "failed its retrospective", fixed = TRUE)
})

test_that("review builder declares rendered and source-to-artifact surfaces", {
  paths <- review_bundle_paths()
  text <- paste(readLines(paths$builder, warn = FALSE), collapse = "\n")
  for (term in c(
    "render_review_articles", "SOURCE-TO-ARTIFACT.csv", "ENVIRONMENT.txt",
    "common-geometry-equivalence.Rmd", "matched-interpretability.Rmd",
    "population-form.Rmd", "protocols/prospective"
  )) {
    expect_match(text, term, fixed = TRUE)
  }
})

test_that("review bundle rebuilds, verifies, and rejects a substitution", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("CROSSFORM_RUN_REVIEW_BUNDLE_TESTS"),
                        "true"),
              "set CROSSFORM_RUN_REVIEW_BUNDLE_TESTS=true for bundle court")
  skip_if_not(nzchar(Sys.which("tar")), "tar is required")
  paths <- review_bundle_paths()
  output <- tempfile("crossform-review-output-")
  dir.create(output)
  rscript <- file.path(R.home("bin"), "Rscript")
  date <- "2099-01-01"
  tarball <- file.path(output,
                       paste0("crossform-review-bundle-", date, ".tar.gz"))

  build <- system2(rscript, c(
    paths$builder, paste0("--output-dir=", output), paste0("--date=", date),
    "--force"
  ), stdout = TRUE, stderr = TRUE)
  expect_null(attr(build, "status"), info = paste(build, collapse = "\n"))
  expect_true(file.exists(tarball))

  verified <- system2(rscript, c(paths$builder, paste0("--verify=", tarball)),
                      stdout = TRUE, stderr = TRUE)
  expect_null(attr(verified, "status"),
              info = paste(verified, collapse = "\n"))
  expect_true(any(grepl("RESULT: OK", verified, fixed = TRUE)))

  scratch <- tempfile("crossform-review-tamper-")
  dir.create(scratch)
  utils::untar(tarball, exdir = scratch)
  root <- list.files(scratch, full.names = TRUE)
  expect_length(root, 1L)
  target <- file.path(root, "PUBLIC-SCOPE.md")
  writeLines(c(readLines(target, warn = FALSE), "tamper-probe"), target)
  tampered <- file.path(output, "tampered.tar.gz")
  tar_status <- system2("tar", c(
    "-czf", shQuote(tampered), "-C", shQuote(scratch), shQuote(basename(root))
  ))
  expect_identical(tar_status, 0L)

  rejected <- suppressWarnings(system2(
    rscript, c(paths$builder, paste0("--verify=", tampered)),
    stdout = TRUE, stderr = TRUE
  ))
  expect_identical(attr(rejected, "status"), 1L)
  expect_true(any(grepl("MISMATCH: PUBLIC-SCOPE.md", rejected, fixed = TRUE)))
  expect_true(any(grepl("RESULT: FAILED", rejected, fixed = TRUE)))
})
