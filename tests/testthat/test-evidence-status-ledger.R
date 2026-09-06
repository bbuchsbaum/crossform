if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

evidence_ledger_path <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "inst", "extdata", "certification",
                       "evidence-status-ledger.csv"),
    system.file("extdata", "certification", "evidence-status-ledger.csv",
                package = "crossform")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) NA_character_ else candidates[[1L]]
}

evidence_registry_path <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "inst", "extdata", "certification",
                       "evidence-claim-registry.csv"),
    system.file("extdata", "certification", "evidence-claim-registry.csv",
                package = "crossform")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) NA_character_ else candidates[[1L]]
}

test_that("the evidence ledger uses only the eight governed classes", {
  path <- evidence_ledger_path()
  skip_if(is.na(path), "evidence status ledger is not present")
  ledger <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_identical(
    names(ledger),
    c("claim_id", "headline_claim", "evidence_class", "evidence_artifact",
      "boundary", "headline_documents", "protocol_frozen_before_analysis",
      "external_route", "completed")
  )
  allowed <- c(
    "algebraic_theorem", "internal_oracle", "external_parity",
    "matched_simulation", "existing_illustration", "prospective_protocol",
    "completed_real_data_result", "independent_replication"
  )
  expect_true(all(ledger$evidence_class %in% allowed))
  expect_true(all(grepl("^CF-H[0-9]{2}$", ledger$claim_id)))
  expect_false(any(!nzchar(ledger$headline_claim)))
  expect_false(any(!nzchar(ledger$boundary)))
  expect_identical(anyDuplicated(paste(
    ledger$claim_id, ledger$evidence_class, ledger$evidence_artifact
  )), 0L)
  repo <- testthat::test_path("..", "..")
  expect_true(all(file.exists(file.path(repo, ledger$evidence_artifact))))
})

test_that("prospective evidence cannot present itself as completed", {
  ledger <- utils::read.csv(evidence_ledger_path(), stringsAsFactors = FALSE)
  prospective <- ledger$evidence_class == "prospective_protocol"
  expect_true(any(prospective))
  expect_true(all(!ledger$completed[prospective]))
  expect_true(all(!ledger$protocol_frozen_before_analysis[prospective]))
  prohibited <- "demonstrated|validated|confirmed|replicated"
  prospective_text <- paste(ledger$headline_claim[prospective],
                            ledger$boundary[prospective])
  expect_false(any(grepl(prohibited, prospective_text, ignore.case = TRUE)))

  # No real-data or replication class may appear without its prerequisites.
  completed <- ledger$evidence_class == "completed_real_data_result"
  replicated <- ledger$evidence_class == "independent_replication"
  expect_true(all(ledger$completed[completed]))
  expect_true(all(ledger$protocol_frozen_before_analysis[completed]))
  expect_true(all(ledger$completed[replicated]))
  expect_true(all(ledger$protocol_frozen_before_analysis[replicated]))
  expect_true(all(ledger$external_route[replicated]))
})

test_that("Haxby and ds003745 stay illustrative unless prospectively frozen", {
  ledger <- utils::read.csv(evidence_ledger_path(), stringsAsFactors = FALSE)
  haxby <- grepl("haxby", ledger$evidence_artifact, ignore.case = TRUE) |
    grepl("haxby", ledger$boundary, ignore.case = TRUE)
  ds003745 <- grepl("population-slice2", ledger$evidence_artifact,
                   ignore.case = TRUE) |
    grepl("ds003745", ledger$boundary, ignore.case = TRUE)
  expect_true(any(haxby))
  expect_true(any(ds003745))
  expect_true(all(ledger$evidence_class[haxby] %in%
                    c("external_parity", "existing_illustration")))
  expect_true(all(ledger$evidence_class[ds003745] == "existing_illustration"))
  expect_true(all(!ledger$protocol_frozen_before_analysis[haxby | ds003745]))
  expect_false(any(ledger$evidence_class %in%
                     c("completed_real_data_result", "independent_replication")))
})

test_that("every headline document is mapped to a governed claim", {
  ledger <- utils::read.csv(evidence_ledger_path(), stringsAsFactors = FALSE)
  covered <- unique(unlist(strsplit(ledger$headline_documents, ";",
                                    fixed = TRUE)))
  required <- c(
    "README.md",
    "vignettes/introduction.Rmd",
    "vignettes/from-observations.Rmd",
    "vignettes/evidence-pairing.Rmd",
    "vignettes/conservative-frames.Rmd",
    "vignettes/population-form.Rmd",
    "vignettes/common-geometry-equivalence.Rmd",
    "vignettes/novelty.Rmd",
    "design/unification-contract.md",
    "design/population-estimand-contract.md"
  )
  expect_setequal(intersect(covered, required), required)
  repo <- testthat::test_path("..", "..")
  expect_true(all(file.exists(file.path(repo, required))))
})

test_that("the normative prose includes every required non-goal and rule", {
  path <- testthat::test_path("..", "..", "design",
                             "evidence-status-ledger.md")
  contract <- paste(readLines(path, warn = FALSE), collapse = "\n")
  for (term in c(
    "Preprocessing", "Universal HRF or GLM modeling",
    "Classification and prediction", "Nonlinear or adaptive distances",
    "Every RSA method", "Transport learning",
    "Universal downstream inference"
  )) {
    expect_match(contract, term, fixed = TRUE)
  }
  expect_match(contract, "A `prospective_protocol` becomes a `completed_real_data_result`",
               fixed = TRUE)
  expect_match(contract, "completed result becomes an `independent_replication`",
               fixed = TRUE)
  expect_match(contract, "support no claim that functional", fixed = TRUE)
})

test_that("every governed claim has one owner and one current status", {
  ledger <- utils::read.csv(evidence_ledger_path(), stringsAsFactors = FALSE)
  registry <- utils::read.csv(evidence_registry_path(), stringsAsFactors = FALSE)
  expect_identical(
    names(registry),
    c("claim_id", "owner", "current_status", "strongest_current_class",
      "strongest_current_artifact", "limitation", "promotion_history")
  )
  expect_identical(anyDuplicated(registry$claim_id), 0L)
  expect_setequal(registry$claim_id, unique(ledger$claim_id))
  expect_equal(nrow(registry), 14L)
  expect_true(all(nzchar(registry$owner)))
  expect_true(all(nzchar(registry$limitation)))
  expect_true(all(registry$current_status %in% c(
    "supported_bounded", "supported_synthetic_regimes",
    "unsupported_retrospective", "blocked_not_executed",
    "blocked_not_implemented"
  )))
  repo <- testthat::test_path("..", "..")
  expect_true(all(file.exists(file.path(
    repo, registry$strongest_current_artifact
  ))))
  expect_true(all(registry$strongest_current_class %in%
                    ledger$evidence_class))
})

test_that("population and protocol current states preserve their boundaries", {
  registry <- utils::read.csv(evidence_registry_path(), stringsAsFactors = FALSE)
  status <- stats::setNames(registry$current_status, registry$claim_id)
  expect_identical(unname(status[c("CF-H09", "CF-H10")]),
                   rep("supported_synthetic_regimes", 2L))
  expect_identical(unname(status[["CF-H11"]]), "unsupported_retrospective")
  expect_identical(unname(status[["CF-H12"]]), "blocked_not_executed")
  expect_identical(unname(status[["CF-H14"]]), "blocked_not_implemented")

  h06 <- registry[registry$claim_id == "CF-H06", , drop = FALSE]
  expect_identical(h06$strongest_current_class, "external_parity")
  expect_match(h06$limitation, "version-pinned", fixed = TRUE)

  history_path <- testthat::test_path("..", "..", "design",
                                     "evidence-promotion-history.md")
  history <- paste(readLines(history_path, warn = FALSE), collapse = "\n")
  for (term in c("CF-H09", "CF-H10", "CF-H11", "CF-H12", "CF-H14",
                 "blocked_not_executed", "independent_replication")) {
    expect_match(history, term, fixed = TRUE)
  }
})
