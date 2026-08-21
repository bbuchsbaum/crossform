if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

prospective_readiness_environment <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "protocols", "prospective",
                                "readiness.R"), envir = environment)
  environment
}

readiness_checks <- function(value = TRUE) stats::setNames(
  as.list(rep(value, 6L)), c("eligible_data", "frozen_configuration",
    "environment_locked", "provenance_bound", "analyst_identified",
    "artifact_storage_reserved"))

test_that("readiness reports only blocked ready or executed with reasons", {
  gate <- prospective_readiness_environment()$prospective_readiness
  blocked_checks <- readiness_checks()
  blocked_checks$eligible_data <- FALSE
  blocked <- gate(blocked_checks)
  ready <- gate(readiness_checks())
  execution <- list(manifest_verified = TRUE,
    protocol_hash_predates_outcomes = TRUE, deviations_logged = TRUE,
    artifact_root = "artifacts/discovery-001")
  executed <- gate(readiness_checks(), execution = execution)

  expect_identical(blocked$state, "BLOCKED")
  expect_identical(blocked$reasons, "eligible_data")
  expect_identical(ready$state, "READY")
  expect_false(ready$executed)
  expect_identical(ready$evidence_state, "prospective_protocol")
  expect_identical(executed$state, "EXECUTED")
  expect_identical(executed$evidence_state, "completed_real_data_result")
})

test_that("replication requires separate complete execution evidence", {
  gate <- prospective_readiness_environment()$prospective_readiness
  execution <- list(manifest_verified = TRUE,
    protocol_hash_predates_outcomes = TRUE, deviations_logged = TRUE,
    artifact_root = "artifacts/discovery-001")
  same <- gate(readiness_checks(), phase = "replication",
    execution = execution, discovery_artifact_root = "artifacts/discovery-001")
  execution$artifact_root <- "artifacts/replication-001"
  replicated <- gate(readiness_checks(), phase = "replication",
    execution = execution, discovery_artifact_root = "artifacts/discovery-001")

  expect_identical(same$state, "BLOCKED")
  expect_identical(same$reasons, "replication_artifact_must_be_distinct")
  expect_identical(replicated$state, "EXECUTED")
  expect_identical(replicated$evidence_state, "independent_replication")
})

test_that("the current state is blocked and makes no result claim", {
  current <- jsonlite::read_json(testthat::test_path(
    "..", "..", "protocols", "prospective", "readiness-current.json"),
    simplifyVector = TRUE)
  expect_identical(current$state, "BLOCKED")
  expect_false(current$ready_implies_result)
  expect_false(current$executed)
  expect_identical(current$evidence_state, "prospective_protocol")
  expect_true(length(current$reasons) > 0)
})
