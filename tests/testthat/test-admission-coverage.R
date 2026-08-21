# The admission coverage list and the promote table must stay aligned.
#
# A certified verb without a promote row, or a shipped gate with no listed
# role, is the failure mode this test replaces. The coverage file lives in
# `benchmarks/`, which is not installed, so the check runs only when a source
# checkout is reachable.

admission_repo_root <- function() {
  candidates <- c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    testthat::test_path("..", "..", "00_pkg_src", "crossform")
  )
  for (root in candidates) {
    if (file.exists(file.path(root, "benchmarks", "admission-coverage.R"))) {
      return(normalizePath(root, mustWork = TRUE))
    }
  }
  NULL
}

admission_tokens <- function(verbs) {
  unique(unlist(strsplit(verbs, ",\\s*"), use.names = FALSE))
}

test_that("the promote table is exactly the certified and refused coverage rows", {
  root <- admission_repo_root()
  skip_if(is.null(root), "package sources are not available under this runner")

  environment <- new.env(parent = globalenv())
  sys.source(
    file.path(root, "benchmarks", "admission-coverage.R"),
    envir = environment
  )
  coverage <- environment$.crossform_admission_coverage()
  promotable <- environment$.crossform_promotable_artifacts(coverage)

  known_roles <- c("certified", "covered", "refused", "local", "gap")
  expect_length(setdiff(coverage$role, known_roles), 0L)
  expect_true(all(c("certified", "covered", "refused", "local") %in%
    coverage$role))

  shipped_roles <- coverage$role %in% c("certified", "refused")
  shipped <- coverage[shipped_roles & !is.na(coverage$artifact) &
    nzchar(coverage$artifact), ]
  expect_identical(anyDuplicated(shipped$artifact), 0L)
  expect_setequal(names(promotable), shipped$artifact)
  expect_identical(unname(promotable[shipped$artifact]), shipped$runner)

  covered <- coverage[coverage$role == "covered", ]
  expect_true(all(covered$artifact %in% shipped$artifact))

  gaps <- coverage[coverage$role == "gap", ]
  expect_true(all(is.na(gaps$artifact) | !nzchar(gaps$artifact)))
  expect_false(any(gaps$verbs %in% names(promotable)))
  expect_false(any(gaps$verbs == "measurement_form"))

  measurement <- coverage[coverage$verbs == "measurement_form", ]
  expect_identical(nrow(measurement), 1L)
  expect_identical(measurement$role, "certified")
  expect_identical(measurement$artifact, "measurement-profile.rds")
  expect_identical(measurement$runner, "run-measurement-profile.R")
  expect_identical(measurement$reader, "test-certification-artifacts.R")

  locals <- coverage[coverage$role == "local", ]
  expect_false(any(locals$artifact %in% names(promotable)))

  out_of_scope <- environment$.crossform_admission_out_of_scope()
  certified_verbs <- admission_tokens(
    coverage$verbs[coverage$role %in% c("certified", "covered", "gap")]
  )
  expect_length(intersect(out_of_scope, certified_verbs), 0L)
})

test_that("every coverage runner, reader, and shipped artifact is present", {
  root <- admission_repo_root()
  skip_if(is.null(root), "package sources are not available under this runner")

  environment <- new.env(parent = globalenv())
  sys.source(
    file.path(root, "benchmarks", "admission-coverage.R"),
    envir = environment
  )
  coverage <- environment$.crossform_admission_coverage()
  promotable <- environment$.crossform_promotable_artifacts(coverage)

  runners <- unique(coverage$runner[!is.na(coverage$runner)])
  expect_true(all(file.exists(file.path(root, "benchmarks", runners))))

  readers <- unique(coverage$reader[!is.na(coverage$reader)])
  expect_true(all(file.exists(file.path(root, "tests", "testthat", readers))))

  promote_text <- paste(
    readLines(file.path(root, "benchmarks", "promote-artifacts.R")),
    collapse = "\n"
  )
  expect_match(promote_text, "admission-coverage\\.R")
  expect_match(promote_text, "\\.crossform_promotable_artifacts")

  shipped <- list.files(
    file.path(root, "inst", "extdata", "certification"),
    pattern = "\\.rds$"
  )
  expect_setequal(shipped, names(promotable))

  readme <- paste(
    readLines(file.path(root, "benchmarks", "README.md")),
    collapse = "\n"
  )
  recertify <- paste(
    readLines(file.path(root, "benchmarks", "RECERTIFY.md")),
    collapse = "\n"
  )
  expect_match(recertify, "admission-coverage\\.R")
  expect_match(recertify, "measurement_form")

  expect_match(readme, "Map-scale admission coverage")
  expect_match(readme, "measurement_form")
  certified <- coverage$artifact[coverage$role == "certified"]
  missing <- certified[!vapply(certified, function(name) {
    grepl(name, readme, fixed = TRUE)
  }, logical(1))]
  expect_identical(missing, character())
})
