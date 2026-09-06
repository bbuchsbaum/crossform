test_that("[T45 T46 T47 T52 T54 T61] recorded predictive risk calibration binds and meets fixed precision", {
  x <- certified_artifact("predictive-geometry-validation.rds", "predictive-geometry/certify.R")
  expect_identical(x$schema, "predictive-certification-v1")
  expect_equal(unname(x$sampling$counts), c(20000, 27000, 20000, 20000))
  expect_equal(unname(x$selection$counts), c(2000, 5000))
  expect_identical(x$sampling$config$production_seed, 2026090403L)
  expect_identical(x$selection$config$production_seed, 2026090404L)
  for (a in list(x$sampling, x$selection)) {
    expect_identical(a$config_hash, digest::digest(a$config, algo = "sha256"))
    expect_equal(unname(a$seeds), a$config$production_seed + 100L * seq_along(a$counts))
    s <- a$summary; c <- a$config
    expect_identical(c$confidence_multiplier, 5)
    expect_identical(c$relative_equivalence_margin, .02)
    expect_true(all(is.finite(s$error) & is.finite(s$mcse) & s$mcse >= 0))
    expect_true(all(abs(s$error) <= 5 * s$mcse + c$numerical_tolerance))
    expect_true(all(abs(s$error) + 5 * s$mcse <= .02 * s$reference_scale))
    expect_equal(s$normalized_error, s$error / s$reference_scale, tolerance = 1e-12)
    expect_equal(s$n, unname(a$counts[s$arm]))
  }
  g <- x$sampling$generator_moments
  expect_equal(g$error, g$observed - g$exact, tolerance = 1e-12)
  expect_true(all(abs(g$error) <= 5 * g$mcse + 1e-10))
  null <- subset(x$sampling$summary, arm == "null")
  expect_true(all(null$mean_leakage_error - 5 * null$leakage_mcse >
    .005 * null$reference_scale))
  null <- subset(x$selection$summary, arm == "null")
  expect_true(all(null$leakage_bias - 5 * null$leakage_mcse > .001 * null$reference_scale))
  paired <- x$sampling$paired_comparisons
  expect_identical(nrow(paired), 8L)
  expect_true(all(is.finite(paired$mean) & is.finite(paired$mcse) & paired$mcse > 0))
})

test_that("[T25 T43 T53 T56 T57 T61] recorded scale and reuse evidence retain numerical parity and measured limits", {
  x <- certified_artifact("predictive-geometry-validation.rds", "predictive-geometry/certify.R")
  p <- x$performance; s <- p$summary
  expect_setequal(paste(s$case, s$route), c(paste(rep(c("tiny", "medium"), each = 3),
    rep(c("memory", "block", "dense"), 2)), "large block"))
  expect_true(all(s$form_error < 1e-8 & s$gain_error < 1e-8))
  expect_true(all(s$planned_bytes <= 512 * 1024^2))
  expect_true(all(s$seconds > 0 & s$retained_bytes > 0 & s$read_calls > 0 & s$kernels > 0))
  expect_true(all(is.finite(s$incremental_rss_bytes) & s$incremental_rss_bytes >= 0))
  rank <- c(tiny = 4, medium = 12, large = 20)[s$case]
  q <- c(tiny = 12, medium = 64, large = 200)[s$case]
  width <- ifelse(s$route == "dense", q * (q + 1)/2, rank * (rank + 1)/2)
  expect_equal(s$packed_width, unname(width))
  expect_true(all(s$maximum_read_width <= 128))
  for (name in c("tiny", "medium")) {
    rows <- s[s$case == name, ]
    expect_length(unique(rows$read_calls), 1L)
    expect_length(unique(rows$read_columns), 1L)
  }
  for (a in p$cases) {
    expect_true(all(a$reference_errors < 1e-8))
    expect_length(a$reference_rows, if (a$case$id == "large") 64L else a$case$nodes)
    if (!is.null(a$reuse)) {
      expect_identical(a$reuse$cached$calls$materializations, 2L)
      expect_identical(a$reuse$uncached$calls$materializations, 12L)
      expect_lt(a$reuse$cached$calls$paths, a$reuse$uncached$calls$paths)
      expect_true(all(a$reuse$errors < 1e-8))
    }
  }
  expect_true(any(grepl("not process RSS", x$limits, fixed = TRUE)))
})

test_that("[T61 T62] real production mutations have exercised positive controls and semantic kills", {
  x <- certified_artifact("predictive-geometry-validation.rds", "predictive-geometry/certify.R")
  m <- x$mutations
  expect_identical(nrow(m$matrix), 20L)
  expect_true(all(m$matrix$status == "KILLED"))
  expect_setequal(vapply(m$witnesses, `[[`, "", "id"), m$matrix$id)
  for (w in m$witnesses) {
    expect_identical(w$baseline_failures, 0L)
    expect_identical(w$control_failures, 0L)
    expect_gt(w$baseline_expectations, 0L)
    expect_gt(w$primary_failures, 0L)
    expect_true(all(w$control_hits > 0 & w$mutated_hits > 0))
    expect_setequal(names(w$mutated_hits), names(w$original_body_digests))
    expect_setequal(names(w$mutated_body_digests), names(w$original_body_digests))
  }
  expect_identical(x$oracle$exit_code, 0L)
  expect_true(any(grepl("PASS", x$oracle$output, fixed = TRUE)))
  root <- .certification_source_root()
  if (!is.na(root)) {
    hashes <- c(x$harness_digests, x$sampling$harness, x$selection$harness,
      x$performance$harness_digests, x$mutations$harness_digests)
    actual <- vapply(names(hashes), function(path)
      digest::digest(file = file.path(root, path), algo = "sha256"), "")
    expect_identical(actual, hashes)
    for (w in m$witnesses) expect_identical(w$test_digest,
      digest::digest(file = file.path(root, w$test_file), algo = "sha256"))
  }
})

test_that("[T61] promotion refuses stale source before copying an otherwise valid receipt", {
  root <- .certification_source_root()
  skip_if(is.na(root), "Promotion tooling is available only in a source checkout")
  temp <- tempfile("promotion-binding-"); dir.create(temp)
  withr::defer(unlink(temp, recursive = TRUE))
  repo <- file.path(temp, "repo"); results <- file.path(temp, "results")
  dir.create(repo); dir.create(results); dir.create(file.path(repo, "R"))
  dir.create(file.path(repo, "benchmarks"))
  writeLines("fixture <- 1", file.path(repo, "R", "fixture.R"))
  for (name in c("provenance.R", "admission-coverage.R"))
    expect_true(file.copy(file.path(root, "benchmarks", name), file.path(repo, "benchmarks", name)))
  env <- new.env(parent = globalenv())
  sys.source(file.path(repo, "benchmarks", "provenance.R"), env)
  name <- "predictive-geometry-validation.rds"
  artifact <- list(provenance = list(source_digest = "sha256:stale"))
  saveRDS(artifact, file.path(results, name))
  run <- function() processx::run(file.path(R.home("bin"), "Rscript"),
    c(file.path(root, "benchmarks", "promote-artifacts.R"), results, repo),
    error_on_status = FALSE)
  stale <- run()
  expect_identical(stale$status, 1L)
  expect_match(stale$stdout, "stale source", fixed = TRUE)
  expect_false(file.exists(file.path(repo, "inst", "extdata", "certification", name)))
  artifact$provenance$source_digest <- env$.crossform_source_tree_digest(repo)
  saveRDS(artifact, file.path(results, name))
  current <- run()
  expect_identical(current$status, 0L)
  expect_identical(readRDS(file.path(repo, "inst", "extdata", "certification", name)), artifact)
})
