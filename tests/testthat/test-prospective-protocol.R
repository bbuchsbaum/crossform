if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

prospective_protocol_text <- function(name) paste(readLines(testthat::test_path(
  "..", "..", "protocols", "prospective", name
), warn = FALSE), collapse = "\n")

test_that("dataset eligibility is outcome blind and complete", {
  text <- prospective_protocol_text("eligibility-v1.md")
  for (term in c("At least 24", "four independent evaluation partitions",
    "three prespecified conditions", "80 percent", "70 percent",
    "Participant metadata", "Documented preprocessing", "license",
    "circular or outcome-selected", "unavailable independent partitions",
    "incompatible condition coding", "immutable dataset/version")) {
    expect_match(text, term, fixed = TRUE)
  }
  expect_match(text, "before any desired contrast", fixed = TRUE)
  expect_match(text, "new dated file", fixed = TRUE)
})

test_that("existing real-data examples are not prospectively promoted", {
  text <- prospective_protocol_text("eligibility-v1.md")
  expect_match(text, "six participants", fixed = TRUE)
  expect_match(text, "only 12 participants", fixed = TRUE)
  expect_match(text, "remains existing illustration", fixed = TRUE)
  expect_match(text, "cannot retroactively", fixed = TRUE)
})

test_that("comparator ambiguity is quantitative matched and falsifiable", {
  text <- prospective_protocol_text("comparators-v1.md")
  for (term in c("identical admitted participants", "Activation:",
    "Aggregate MVPA:", "Fixed linear RSA:", "Crossform:",
    "all-unordered-partition-pairs crossnobis", "95 percent pointwise HC3",
    "15 percent", "at least 0.20", "two adjacent nonpoint scales",
    "participant ratios are not averaged", "Statistical significance alone",
    "Crossform fails", "Condition-label", "organization-label",
    "one-subject influence")) {
    expect_match(text, term, fixed = TRUE)
  }
})

test_that("the discovery specification is machine readable and frozen", {
  path <- testthat::test_path("..", "..", "protocols", "prospective",
                             "discovery-v1.json")
  config <- jsonlite::read_json(path, simplifyVector = TRUE)
  expect_identical(config$schema_version, "crossform-prospective-discovery-v1")
  expect_identical(config$status, "BLOCKED_DATASET_NOT_SELECTED")
  expect_null(config$dataset$dataset_id)
  expect_null(config$dataset$manifest_sha256)
  expect_identical(config$condition_roles,
                   c("condition_a", "condition_b", "negative_control"))
  expect_identical(config$partitions$minimum_partitions, 4L)
  expect_identical(config$spatial$scales_mm, c(0L, 4L, 8L, 12L))
  expect_identical(config$population$formula, "~ 1")
  expect_identical(config$population$coverage_policy, "all_planned")
  expect_identical(config$uncertainty$primary, "HC3_pointwise_95_percent")
  expect_identical(config$uncertainty$bootstrap_replicates, 1999L)
  expect_identical(config$uncertainty$simultaneous_coverage, "not_available")
  expect_true(all(c("primary-results.csv", "primary-figure.png",
                    "failure-states.csv", "execution-manifest.csv",
                    "deviations.md") %in% config$outputs))
  expect_true(all(c("missing_partition", "coverage_below_floor",
                    "transport_evaluation_overlap", "conservation_failed",
                    "interpretive_gain_not_established") %in%
                    config$failure_states))
  expect_false(config$exploratory$headline_eligible)
})

test_that("frozen discovery prose forbids overwrite and fabricated results", {
  text <- prospective_protocol_text("discovery-v1.md")
  expect_match(text, "not fabricated result values", fixed = TRUE)
  expect_match(text, "cannot overwrite", fixed = TRUE)
  expect_match(text, "never headline eligible", fixed = TRUE)
  amendment <- prospective_protocol_text("amendments/README.md")
  expect_match(amendment, "append-only", fixed = TRUE)
  expect_match(amendment, "never", fixed = TRUE)
})

test_that("replication protocol has two genuinely distinct ready routes", {
  text <- prospective_protocol_text("replication-v1.md")
  for (term in c("no replication has occurred", "Route A:", "Route B:",
    "no participant, acquisition session", "outside the discovery analysis team",
    "Optional blinding", "named custodian", "Shared code is disclosed",
    "Two named reviewers", "without overwriting", "exact analyst wording",
    "different artifact roots", "does not license the word replicated")) {
    expect_match(text, term, fixed = TRUE)
  }
  expect_match(text, "same direction and scale ordering", fixed = TRUE)
  log <- prospective_protocol_text("external-feedback-log.md")
  expect_match(log, "Append-only", fixed = TRUE)
  expect_match(log, "No entries", fixed = TRUE)
})
