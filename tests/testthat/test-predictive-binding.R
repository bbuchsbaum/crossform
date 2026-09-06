test_that("[T61] the bare-R binding gate rejects a changed oracle at unchanged production source", {
  root <- .certification_source_root()
  skip_if(is.na(root), "Binding tooling is available only in a source checkout")
  x <- certified_artifact("predictive-geometry-validation.rds", "predictive-geometry/certify.R")
  temp <- tempfile("harness-binding-"); dir.create(temp)
  withr::defer(unlink(temp, recursive = TRUE))
  manifests <- c(x$harness_digests, x$sampling$harness, x$selection$harness,
    x$performance$harness_digests, x$mutations$harness_digests)
  files <- unique(c(names(manifests), "benchmarks/check-certification-binding.R",
    "benchmarks/provenance.R", paste0("R/", list.files(file.path(root, "R"),
      pattern = "[.]R$", full.names = FALSE))))
  # Use explicit relative paths so the gate cannot accidentally read the
  # original checkout's oracle after we alter this isolated copy.
  files <- unique(c(files, paste0("inst/extdata/certification/",
    list.files(file.path(root, "inst/extdata/certification"), full.names = FALSE))))
  for (file in files) {
    target <- file.path(temp, file); dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(file.path(root, file), target))
  }
  run <- function() processx::run(file.path(R.home("bin"), "Rscript"),
    c(file.path(temp, "benchmarks/check-certification-binding.R"), temp), error_on_status = FALSE)
  expect_identical(run()$status, 0L)
  cat("\n# deliberately changed oracle\n", file = file.path(temp, x$oracle$file), append = TRUE)
  changed <- run()
  expect_identical(changed$status, 1L)
  expect_match(changed$stdout, "stale or missing harness: design/oracles/predictive-geometry.R", fixed = TRUE)
})
