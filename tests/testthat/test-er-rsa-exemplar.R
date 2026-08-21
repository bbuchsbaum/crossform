# The recorded rectangular encoding-retrieval RSA result, held as a ratchet.
#
# `exemplars/er-rsa` simulates a designed encoding-retrieval experiment with
# unequal item sets, missing matches, a trial covariate, run nuisance, and
# planted regional structure, then recovers all of it through a rectangular
# crossform plan. This file does not recompute the exemplar (it takes about
# twenty seconds and writes to the exemplar directory); it asserts that the
# committed record still says what the exemplar README and the novelty
# ledger's gate 3 claim it says. Loosening a verdict, dropping a readout, or
# pasting a worse number into the CSVs fails here.
#
# The file is skipped when the exemplar results are absent (an installed
# package, or a checkout where the exemplar has not been run).

er_rsa_result <- function(file) {
  candidates <- c(
    testthat::test_path("..", "..", "exemplars", "er-rsa", "results", file),
    system.file("exemplars", "er-rsa", "results", file, package = "crossform")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) NA_character_ else normalizePath(candidates[[1L]])
}

er_rsa_read <- function(file) {
  path <- er_rsa_result(file)
  skip_if(is.na(path), "er-rsa exemplar results are not present")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

test_that("every recorded er-rsa recovery verdict still passes", {
  verdicts <- er_rsa_read("recovery-verdicts.csv")

  expect_setequal(
    names(verdicts),
    c("claim", "statistic", "criterion", "passed", "detail")
  )
  expect_gte(nrow(verdicts), 8L)
  expect_true(all(verdicts$passed))
  expect_true(all(nzchar(verdicts$detail)))

  # The verdicts that carry the scientific content of gate 3 must be present
  # by name, so a future edit cannot pass by deleting the hard ones.
  expect_true(any(grepl("item-specific reinstatement", verdicts$claim)))
  expect_true(any(grepl("no item reinstatement", verdicts$claim)))
  expect_true(any(grepl("carries nothing", verdicts$claim)))
  expect_true(any(grepl("study-duration covariate", verdicts$claim)))
  expect_true(any(grepl("same-run pairing", verdicts$claim)))
})

test_that("the recorded er-rsa estimates recover their planted values", {
  recovery <- er_rsa_read("planted-vs-estimated.csv")
  cross <- recovery[recovery$plan == "cross", ]

  expect_gte(nrow(cross), 30L)
  expect_true(all(is.finite(cross$planted)))
  expect_true(all(cross$se > 0))

  # Every cross-run readout's planted value lies inside its subject CI, and
  # no readout is biased by as much as two standard errors.
  covered <- cross$planted >= cross$ci_low & cross$planted <= cross$ci_high
  expect_true(all(covered))
  expect_lt(max(abs(cross$bias_t)), 2)
  expect_lt(max(abs(cross$bias)), 0.02)
})

test_that("the recorded er-rsa regional dissociation is intact", {
  levels_table <- er_rsa_read("coupling-levels.csv")
  rownames(levels_table) <- levels_table$region
  expect_setequal(levels_table$region, c("regionA", "regionB", "regionC"))

  # regionA: match coupling exceeds both control sets, and the
  # category-matched contrast is a large positive item-specific effect.
  expect_gt(levels_table["regionA", "match_coupling"],
    levels_table["regionA", "control_coupling_category"])
  expect_gt(levels_table["regionA", "item_specific_contrast"], 0.45)
  expect_gt(levels_table["regionA", "item_specific_t"], 20)

  # regionB: category structure only. The naive contrast is large, the
  # category-matched contrast is not.
  expect_gt(levels_table["regionB", "naive_contrast"], 0.8)
  expect_lt(abs(levels_table["regionB", "item_specific_contrast"]), 0.05)

  # regionC: nothing, on either control set.
  expect_lt(abs(levels_table["regionC", "item_specific_contrast"]), 0.05)
  expect_lt(abs(levels_table["regionC", "naive_contrast"]), 0.05)

  # The recorded estimates still match the planted values they claim.
  expect_lt(max(abs(levels_table$item_specific_contrast -
    levels_table$item_specific_planted)), 0.02)
  expect_lt(max(abs(levels_table$naive_contrast -
    levels_table$naive_planted)), 0.02)
})

test_that("all five pair-query constructors are exercised on record", {
  diagnostics <- er_rsa_read("query-diagnostics.csv")

  expect_true(any(grepl("match_coupling", diagnostics$constructor)))
  expect_true(any(grepl("control_coupling", diagnostics$constructor)))
  expect_true(any(diagnostics$constructor == "coupling_contrast"))
  expect_true(any(diagnostics$constructor == "match_control"))
  expect_true(any(diagnostics$constructor == "pair_lm_query"))

  # Rectangular, unequal, with missing matches: the plain normalized
  # difference is not additive-baseline invariant here, and the
  # nuisance-adjusted regression forms are. That asymmetry is the reason
  # both constructors exist, so it is asserted rather than described.
  contrast <- diagnostics[diagnostics$constructor == "coupling_contrast", ]
  expect_true(all(!contrast$additive_baseline_invariant))
  adjusted <- diagnostics[diagnostics$constructor == "pair_lm_query", ]
  expect_true(all(adjusted$additive_baseline_invariant))
  expect_true(all(adjusted$rank == adjusted$design_columns))
})

test_that("the recorded er-rsa routes and refusals are unchanged", {
  route <- er_rsa_read("route-check.csv")
  expect_lt(max(route$max_abs_difference), 1e-12)
  expect_true(any(route$query == "total - (coherent + configuration)"))

  refusals <- er_rsa_read("refusals.csv")
  expect_true("rectangular_fixed_metric" %in% refusals$capability)
  expect_true("symmetric_self_form" %in% refusals$capability)
  expect_true(any(grepl("rank deficient", refusals$message)))
})
