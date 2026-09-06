if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

matched_certification_environment <- function() {
  environment <- new.env(parent = globalenv())
  for (file in c(
    "00-mixture-generator.R", "01-multiscale-scenarios.R",
    "02-paired-observations.R", "03-conventional-baselines.R",
    "05-certify.R"
  )) {
    sys.source(testthat::test_path(
      "..", "..", "benchmarks", "matched-interpretability", file
    ), envir = environment)
  }
  environment
}

matched_certification_artifact <- function(name) {
  testthat::test_path("..", "..", "inst", "extdata", "certification", name)
}

test_that("the committed certification passes every prespecified metric", {
  metrics_path <- matched_certification_artifact(
    "matched-interpretability-metrics.csv"
  )
  parameters_path <- matched_certification_artifact(
    "matched-interpretability-parameters.csv"
  )
  skip_if_not(file.exists(metrics_path) && file.exists(parameters_path),
              "matched interpretability certification is not present")
  metrics <- utils::read.csv(metrics_path, stringsAsFactors = FALSE)
  parameters <- utils::read.csv(parameters_path, stringsAsFactors = FALSE)
  parameters <- stats::setNames(parameters$value, parameters$key)

  expect_true(nrow(metrics) >= 20L)
  expect_true(all(metrics$passes))
  expect_setequal(
    unique(metrics$metric),
    c("conservation_max_abs", "recomposition_max_abs",
      "null_false_separation_max_abs", "recoverable_total_relative_bias",
      "recoverable_share_mae", "recoverable_ordering_rate",
      "aggregate_ambiguity_95_bound", "spectrum_separation",
      "activation_negative_control")
  )
  expect_identical(parameters[["schema_version"]],
                   "matched-interpretability-certification-v1")
  expect_identical(parameters[["evidence_claim"]], "CF-H10")
  expect_identical(as.integer(parameters[["replications"]]), 48L)
  expect_identical(parameters[["scenarios"]],
                   "broad_coherent,mixed_broad_fine,fine_configuration")
})

test_that("numeric evidence remains bound to source contract and ledger", {
  path <- matched_certification_artifact(
    "matched-interpretability-checksums.csv"
  )
  skip_if_not(file.exists(path), "certification checksums are not present")
  manifest <- utils::read.csv(path, stringsAsFactors = FALSE)
  repo <- testthat::test_path("..", "..")
  files <- file.path(repo, manifest$path)
  hashed <- manifest$hash_algorithm == "md5"

  expect_true(all(file.exists(files)))
  expect_identical(
    unname(tools::md5sum(files[hashed])),
    manifest$digest[hashed]
  )
  expect_identical(sum(manifest$role == "semantic_plot"), 1L)
  expect_identical(manifest$hash_algorithm[manifest$role == "semantic_plot"],
                   "semantic_render_check")
  expect_true(all(c("contract", "evidence_ledger", "source",
                    "numeric_artifact") %in% manifest$role))
})

test_that("the evidence ledger promotes CF-H10 only to matched simulation", {
  ledger <- utils::read.csv(matched_certification_artifact(
    "evidence-status-ledger.csv"
  ), stringsAsFactors = FALSE)
  all_claim <- ledger[ledger$claim_id == "CF-H10", ]
  expect_true(nrow(all_claim) >= 1L)
  expect_true(all(all_claim$evidence_class == "matched_simulation"))
  claim <- all_claim[
    all_claim$evidence_artifact ==
      "inst/extdata/certification/matched-interpretability-metrics.csv",
    , drop = FALSE
  ]
  expect_identical(nrow(claim), 1L)
  expect_identical(claim$evidence_class, "matched_simulation")
  expect_identical(
    claim$evidence_artifact,
    "inst/extdata/certification/matched-interpretability-metrics.csv"
  )
  expect_true(claim$completed)
  expect_false(claim$protocol_frozen_before_analysis)
  expect_false(claim$external_route)
})

test_that("the fast full-grid smoke reconstruction preserves core laws", {
  generator <- matched_certification_environment()
  smoke <- generator$matched_certification_smoke()
  expect_true(smoke$passes)
  expect_lte(smoke$conservation, 1e-10)
  expect_lte(smoke$recomposition, 1e-10)
  expect_lte(smoke$null_separation, 1e-12)
  expect_identical(smoke$cells, 54L)
  expect_identical(smoke$seeds, c(6201L, 6202L))
})

test_that("certification shares retain signed crossvalidated components", {
  generator <- matched_certification_environment()
  raw <- generator$matched_certification_raw(6201L)
  defined <- is.finite(raw$total) & abs(raw$total) > .Machine$double.eps

  expect_equal(raw$coherent_share[defined],
               raw$coherent[defined] / raw$total[defined])
  expect_true(any(raw$coherent[defined] < 0 |
                    raw$configuration[defined] < 0))
  expect_true(all(is.finite(raw$coherent_share[defined])))
})

test_that("the optional full court rebuilds the committed numeric evidence", {
  skip_if_not(identical(
    Sys.getenv("CROSSFORM_FULL_MATCHED_SIMULATION"), "true"
  ), "set CROSSFORM_FULL_MATCHED_SIMULATION=true for the 48-seed court")
  generator <- matched_certification_environment()
  raw <- generator$matched_certification_raw(5001:5048)
  rebuilt_metrics <- generator$matched_certification_metrics(raw)
  attr(rebuilt_metrics, "thresholds") <- NULL
  rebuilt_summary <- generator$.certification_long_summary(raw)
  recorded_metrics <- utils::read.csv(matched_certification_artifact(
    "matched-interpretability-metrics.csv"
  ), stringsAsFactors = FALSE)
  recorded_summary <- utils::read.csv(matched_certification_artifact(
    "matched-interpretability-summary.csv"
  ), stringsAsFactors = FALSE)
  expect_equal(rebuilt_metrics, recorded_metrics, tolerance = 1e-13)
  expect_equal(rebuilt_summary, recorded_summary, tolerance = 1e-13)
})
