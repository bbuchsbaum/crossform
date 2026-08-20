# Population slice 2 on twelve ds003745 subjects, held as a ratchet.
#
# `exemplars/population-slice2/` carries twelve subjects' native searchlight
# geometries onto 1241 shared consensus grey-matter group nodes under two
# declared transports -- an anatomical nearest-centre P^A and a functionally
# fingerprinted P^F cross-fitted on a different task -- fits one group model
# under each, and commits its receipts. This file does not recompute any of it:
# that needs 7.21 GB of downloaded fMRIPrep derivatives which are deliberately
# gitignored. It asserts that the committed record still says what the exemplar
# README claims.
#
# What slice 1's ratchet could not check, and this one can, is section 7 of
# population-form-v1: `eta_transport` with its null band, and the six
# transport diagnostics. At region level with an identity transport there was
# nothing to diagnose. Here the transport displaces mass by 5-7 mm, spreads it
# over ~2.5 effective nodes, and sinks a fifth of every subject's territory.
#
# THE THING THIS FILE MOST EXISTS TO PROTECT is not the headline eta. It is the
# pair of numbers underneath it. Slice 2's eta came out at +0.93 past all 200
# null draws, and the components say that happened because V^W collapsed by
# 75 % while V^C sat at the null median -- i.e. eta was extreme for a reason
# its own name gets wrong. A future edit that kept `eta_transport` and dropped
# `V_C_PF_percentile_in_null` / `V_W_PF_percentile_in_null` would restore
# exactly the uninterpretable headline contract 7.4 forbids, so those two are
# recomputed here from the committed null band rather than trusted.
#
# What else fails here: dropping a subject, dropping an identity row, loosening
# the tolerance column, pasting a worse identity number in, letting the two
# transports' sinks drift apart (which would let eta be won by discarding
# territory, contract 7.4), letting P^F collapse the group grid onto a few
# nodes (the concentration analogue, which none of the six catches), reporting
# eta without its null band or without the 7.5 diagnostics, or clamping a
# negative eta at zero.
#
# The file skips cleanly when the exemplar results are absent (an installed
# package, or a checkout where the exemplar has not been run).

SLICE2_SUBJECT_IDS <- c("sub-104", "sub-105", "sub-107", "sub-108",
                        "sub-111", "sub-112", "sub-113", "sub-127",
                        "sub-128", "sub-129", "sub-130", "sub-131")
SLICE2_QUERIES <- c("choice-vs-outcome", "social-vs-computer-choice",
                    "recip-vs-defect", "friend-vs-stranger",
                    "friend-vs-computer-choice")
SLICE2_TOLERANCE <- 1e-12

slice2_results <- function(file) {
  candidates <- c(
    testthat::test_path("..", "..", "exemplars", "population-slice2",
                        "results", file),
    system.file("exemplars", "population-slice2", "results", file,
                package = "crossform")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  skip_if(!length(candidates),
          paste("the population slice 2 exemplar result", file,
                "is not present"))
  utils::read.csv(normalizePath(candidates[[1L]]), stringsAsFactors = FALSE)
}

slice2_receipts <- function() {
  slice2_results("population-slice2-receipts.csv")
}

slice2_value <- function(receipts, quantity, subject = "group") {
  rows <- receipts[receipts$quantity == quantity &
                     receipts$subject == subject, ]
  expect_equal(nrow(rows), 1L,
               info = paste("receipt", quantity, "for", subject))
  rows$value[[1L]]
}

test_that("the slice 2 receipts cover twelve subjects and the group", {
  receipts <- slice2_receipts()

  expect_setequal(
    names(receipts),
    c("subject", "quantity", "value", "tolerance", "passes", "note")
  )
  expect_setequal(unique(receipts$subject), c(SLICE2_SUBJECT_IDS, "group"))
  expect_equal(slice2_value(receipts, "n_subjects"),
               length(SLICE2_SUBJECT_IDS))
  # N - rank(X) for the group model ~ 1, under both transports.
  for (tag in c("_PA", "_PF")) {
    expect_equal(slice2_value(receipts, paste0("residual_df", tag)),
                 length(SLICE2_SUBJECT_IDS) - 1L)
  }

  # Every participant carries the same battery, so a subject cannot be kept in
  # the table by contributing only the rows that happen to be easy.
  per_subject <- receipts[receipts$subject %in% SLICE2_SUBJECT_IDS, ]
  split_quantities <- split(per_subject$quantity, per_subject$subject)
  reference <- sort(split_quantities[[SLICE2_SUBJECT_IDS[[1L]]]])
  for (subj in SLICE2_SUBJECT_IDS) {
    expect_identical(sort(split_quantities[[subj]]), reference,
                     info = paste("quantities recorded for", subj))
  }
})

test_that("every asserted slice 2 identity is inside its acceptance", {
  receipts <- slice2_receipts()

  asserted <- receipts[!is.na(receipts$tolerance), ]
  expect_gt(nrow(asserted), 0L)
  # `passes` is exactly `abs(value) <= tolerance`; recomputing it here is what
  # stops a hand-edited TRUE from surviving.
  expect_true(all(asserted$passes))
  expect_equal(abs(asserted$value) <= asserted$tolerance, asserted$passes)

  # The population-layer tolerance is the contract's 1e-12 and may not be
  # loosened. The frame's conservative column-mass check accumulates over
  # ~75 voxels per searchlight and keeps the exemplar's own 1e-10.
  population_rows <- asserted[!grepl("frame_column_mass", asserted$quantity), ]
  expect_true(all(population_rows$tolerance <= SLICE2_TOLERANCE))
})

test_that("the commutation acceptance holds under both transports", {
  receipts <- slice2_receipts()

  # THE section 3 acceptance: query-then-transport equals
  # transport-then-query, at 1241 group nodes under a transport that is not
  # the identity. This is what slice 1 could only demonstrate on a 3x3
  # permutation.
  for (tag in c("_PA", "_PF")) {
    commutation <- slice2_value(receipts, paste0("commutation_claim3", tag))
    expect_lt(abs(commutation), SLICE2_TOLERANCE)
    theta <- slice2_value(receipts, paste0("theta_sum_over_nodes", tag))
    expect_lt(abs(theta), SLICE2_TOLERANCE)
    worst <- slice2_value(receipts, paste0("budget_preservation_worst", tag))
    expect_lt(abs(worst), SLICE2_TOLERANCE)
  }
})

test_that("the sink is real, and identical under the two transports", {
  receipts <- slice2_receipts()

  # Contract 7.4: a transport can raise the consensus share by discarding the
  # nodes that disagree. Slice 2 removes that degree of freedom by pinning
  # P^F's sink to P^A's row by row. If these two ever drift apart, eta stops
  # being a controlled comparison and this test is the thing that says so.
  for (subj in SLICE2_SUBJECT_IDS) {
    gap <- slice2_value(receipts, "sink_identical_max_gap", subj)
    expect_lt(abs(gap), SLICE2_TOLERANCE)
    expect_equal(slice2_value(receipts, "sink_territory_PA", subj),
                 slice2_value(receipts, "sink_territory_PF", subj),
                 tolerance = 1e-12)
  }

  # ...and the sink is carrying real mass, not rounding error: E12 chose this
  # dataset because a third of the union coverage is subject-variable, and
  # that has to survive into the operator or the choice was pointless.
  sink <- vapply(SLICE2_SUBJECT_IDS,
                 function(s) slice2_value(receipts, "sink_territory_PA", s),
                 numeric(1))
  expect_true(all(sink > 0.10))
  expect_true(all(sink < 0.50))

  disagreement <- slice2_value(receipts, "coverage_disagreement_share")
  expect_gt(disagreement, 0.25)

  # Deleting the sink row breaks the two conservation identities by roughly
  # the sink's own size. Slice 1 recorded these as still-zero because its
  # region map had full coverage; here they must NOT be zero, or the sink is
  # not doing the work the exemplar says it does.
  for (tag in c("_PA", "_PF")) {
    expect_gt(abs(slice2_value(receipts,
                               paste0("theta_sum_without_sink", tag))),
              1e-6)
  }
})

test_that("the section 7.5 diagnostics are all present and coherent", {
  receipts <- slice2_receipts()

  # Contract 7.4 is normative: eta may not be reported without these six in
  # the same object. The test enforces the conjunction, because the failure
  # mode it guards against is exactly a future edit that keeps the headline
  # and drops the caveats.
  expect_true("eta_transport" %in% receipts$quantity)
  required_group <- c("V_C_PA", "V_W_PA", "V_W_PA_se", "R_PA",
                      "V_C_PF", "V_W_PF", "V_W_PF_se", "R_PF",
                      "eta_null_sd", "eta_null_q95", "eta_null_draws",
                      "eta_rank_in_null", "coverage_min_subjects_PA",
                      "coverage_nodes_below_floor_PA")
  for (q in required_group) {
    expect_true(q %in% receipts$quantity, info = paste("missing receipt", q))
  }

  for (subj in SLICE2_SUBJECT_IDS) {
    for (q in c("displacement_median_PA: mm", "displacement_p90_PA: mm",
                "displacement_max_PA: mm", "displacement_masswt_mean_PA: mm",
                "displacement_median_PF: mm", "entropy_mean_PA: nats",
                "entropy_mean_PF: nats", "perplexity_mean_PF",
                "exp_mean_entropy_PF", "all_sink_rows",
                "sink_territory_PA")) {
      expect_true(any(receipts$quantity == q & receipts$subject == subj),
                  info = paste("missing", q, "for", subj))
    }

    # P^A is a hard assignment, so its row entropy is identically zero. P^F is
    # not, and if its entropy ever collapsed to zero the "functional" transport
    # would have silently become another nearest-centre map.
    expect_equal(slice2_value(receipts, "entropy_mean_PA: nats", subj), 0)
    expect_gt(slice2_value(receipts, "entropy_mean_PF: nats", subj), 0.1)

    # The two entropy summaries the contract insists are labelled: the row mean
    # of exp(H) and exp(row mean H) are not equal, and by Jensen the first is
    # the larger.
    expect_gt(slice2_value(receipts, "perplexity_mean_PF", subj),
              slice2_value(receipts, "exp_mean_entropy_PF", subj))

    # The transport moves mass a real distance, on the scale of the group grid.
    expect_gt(slice2_value(receipts, "displacement_median_PA: mm", subj), 1)
    expect_lt(slice2_value(receipts, "displacement_max_PA: mm", subj),
              12)
  }
})

test_that("eta_transport is reported signed, with its null band", {
  receipts <- slice2_receipts()
  eta_table <- slice2_results("population-slice2-eta.csv")
  null <- slice2_results("population-slice2-eta-null.csv")

  eta <- slice2_value(receipts, "eta_transport")
  # Contract 7.3: eta is reported unclamped and may be negative. The assertion
  # is that it is FINITE and that it equals R(P^F) - R(P^A) as recorded --
  # never that it is positive. A test that required eta > 0 would be the
  # clamping the contract forbids, moved into the test suite.
  expect_true(is.finite(eta))
  expect_equal(eta,
               slice2_value(receipts, "R_PF") - slice2_value(receipts, "R_PA"),
               tolerance = 1e-10)

  # Both components are reported separately, not only as their ratio (7.3).
  for (tag in c("PA", "PF")) {
    vw <- slice2_value(receipts, paste0("V_W_", tag))
    vc <- slice2_value(receipts, paste0("V_C_", tag))
    r  <- slice2_value(receipts, paste0("R_", tag))
    # V^W must clear the declared floor before the share is formed, and its
    # standard error must be reported so a reader can see by how much: a V^W
    # that clears the floor by less than its own noise is not evidence that
    # anything reproduced.
    expect_gt(vw, 0)
    expect_equal(r, vc / vw, tolerance = 1e-10)
    expect_gt(slice2_value(receipts, paste0("V_W_", tag, "_se")), 0)
  }

  # The null band exists, has the declared number of draws, and eta's rank
  # inside it is recorded. Without this, eta is the uninterpretable headline
  # 7.4 warns about.
  expect_equal(nrow(null), slice2_value(receipts, "eta_null_draws"))
  expect_setequal(names(null), c("draw", "V_C", "V_W", "R", "eta"))
  expect_gt(nrow(null), 100L)
  expect_equal(slice2_value(receipts, "eta_null_sd"), stats::sd(null$eta),
               tolerance = 1e-8)
  rank <- slice2_value(receipts, "eta_rank_in_null")
  expect_equal(rank, sum(null$eta < eta))

  # THE two numbers that say what eta is made of. eta's rank in the band says
  # the transport is unusual; these say which way. They are recomputed from the
  # committed band rather than trusted, because the whole point of this slice's
  # finding is that a headline eta can be extreme for a reason its own name
  # gets wrong -- and a future edit that dropped these would put the headline
  # back on its own.
  expect_equal(slice2_value(receipts, "V_C_PF_percentile_in_null"),
               mean(null$V_C < slice2_value(receipts, "V_C_PF")),
               tolerance = 1e-10)
  expect_equal(slice2_value(receipts, "V_W_PF_percentile_in_null"),
               mean(null$V_W < slice2_value(receipts, "V_W_PF")),
               tolerance = 1e-10)

  expect_setequal(eta_table$quantity,
                  c("V_C_PA", "V_W_PA", "R_PA", "V_W_PA_se",
                    "V_C_PF", "V_W_PF", "R_PF", "V_W_PF_se",
                    "eta_transport", "eta_null_mean", "eta_null_sd",
                    "eta_null_q95", "eta_null_max",
                    "eta_null_negative_share", "eta_null_draws",
                    "eta_rank_in_null", "eta_one_sided_p",
                    "V_C_PF_percentile_in_null",
                    "V_W_PF_percentile_in_null"))
})

test_that("the group ledger is recorded for every query in the bank", {
  receipts <- slice2_receipts()

  # The bank is read in one pass, so a query cannot be dropped from the report
  # by being expensive. `friend-vs-computer-choice` in particular must survive:
  # it is the bare pairwise difference the commutation acceptance compares
  # against, and deleting it would silently disarm that check.
  for (q in SLICE2_QUERIES) {
    for (tag in c("_PA", "_PF")) {
      expect_true(
        any(receipts$quantity == paste0("group_ledger_over_all_nodes", tag,
                                        ": ", q)),
        info = paste("missing group ledger for", q, tag))
      expect_true(
        any(receipts$quantity == paste0("group_ledger_sink", tag, ": ", q)),
        info = paste("missing sink ledger for", q, tag))
    }
  }
})

test_that("the second cross-fit axis is reported beside the first", {
  receipts <- slice2_receipts()
  across <- slice2_results("population-slice2-eta-across-run.csv")
  null <- slice2_results("population-slice2-eta-across-run-null.csv")

  # DECISION.md section 5.3 promised two independence levels and risk 8 says
  # what their disagreement would mean. Reporting only the axis that came out
  # better is the failure mode this guards, so the axis must be present
  # whatever it says -- including when it says nothing.
  expect_true("eta_transport_across_run" %in% across$quantity)
  expect_true("eta_transport_across_run" %in% receipts$quantity)

  # Both energy components are reported for both transports regardless, with
  # the standard error that says whether V^W is distinguishable from zero.
  # Contract 7.3 wants the components, not only their ratio, and here the
  # components are the entire content of the result.
  for (q in c("V_C_PA_acrossrun", "V_W_PA_acrossrun", "V_W_PA_acrossrun_se",
              "V_C_PF2", "V_W_PF2", "V_W_PF2_se")) {
    expect_true(q %in% receipts$quantity, info = paste("missing", q))
    expect_true(is.finite(slice2_value(receipts, q)), info = q)
  }

  # Two structural assertions that hold whether or not a share was formable:
  # the sink still matches P^A, and the permuted rebuild still reproduces the
  # sealed operator under the identity permutation.
  expect_lt(abs(slice2_value(receipts, "eta_across_run_sink_gap")),
            SLICE2_TOLERANCE)
  expect_lt(abs(slice2_value(receipts,
                             "eta_across_run_identity_permutation_gap")),
            SLICE2_TOLERANCE)

  eta2 <- slice2_value(receipts, "eta_transport_across_run")
  formable <- slice2_value(receipts, "eta_across_run_shares_formable") == 1

  if (formable) {
    expect_true(is.finite(eta2))
    expect_equal(eta2,
                 slice2_value(receipts, "R_PF2") -
                   slice2_value(receipts, "R_PA_acrossrun"),
                 tolerance = 1e-10)
    expect_gt(slice2_value(receipts, "V_W_PF2"), 0)
    expect_gt(slice2_value(receipts, "V_W_PA_acrossrun"), 0)
    expect_gt(nrow(null), 100L)
    expect_equal(slice2_value(receipts, "eta_across_run_rank_in_null"),
                 sum(null$eta < eta2))
    expect_equal(slice2_value(receipts, "eta_across_run_null_sd"),
                 stats::sd(null$eta), tolerance = 1e-8)
    # R(P^A) is recomputed on the across-run partitions rather than reused
    # from the across-task ones -- different held-out data, so a shared value
    # would mean one eta was assembled from mismatched halves.
    expect_false(isTRUE(all.equal(slice2_value(receipts, "R_PA"),
                                  slice2_value(receipts, "R_PA_acrossrun"))))
  } else {
    # Contract 7.1: when V^W does not clear the declared floor the share is NA,
    # NOT a large number, and no null band is drawn around it. This branch
    # asserts that the NA is a *declared* NA rather than a crash artefact --
    # eta is NA, the band is empty, and no p-value was manufactured for it.
    expect_true(is.na(eta2))
    expect_true(is.na(slice2_value(receipts, "R_PA_acrossrun")) ||
                  is.na(slice2_value(receipts, "R_PF2")))
    expect_equal(nrow(null), 0L)
    expect_true(is.na(slice2_value(receipts, "eta_across_run_one_sided_p")))
    expect_true(slice2_value(receipts, "V_W_PA_acrossrun") <= 0 ||
                  slice2_value(receipts, "V_W_PF2") <= 0)
  }
})

test_that("the slice 2 transport diagnostics table matches the receipts", {
  receipts <- slice2_receipts()
  diagnostics <- slice2_results(
    "population-slice2-transport-diagnostics.csv")

  expect_setequal(diagnostics$subject, SLICE2_SUBJECT_IDS)
  expect_true(all(diagnostics$n_group_nodes ==
                    slice2_value(receipts, "n_group_nodes")))
  for (subj in SLICE2_SUBJECT_IDS) {
    row <- diagnostics[diagnostics$subject == subj, ]
    expect_equal(row$sink_territory,
                 slice2_value(receipts, "sink_territory_PA", subj),
                 tolerance = 1e-10)
    expect_equal(row$entropy_mean_F,
                 slice2_value(receipts, "entropy_mean_PF: nats", subj),
                 tolerance = 1e-10)
  }

  # Group nodes are consensus grey matter, so every one of them is reachable
  # by every subject; the coverage floor is met with nothing below it.
  expect_equal(slice2_value(receipts, "group_node_min_subject_reach"),
               length(SLICE2_SUBJECT_IDS))
  expect_equal(slice2_value(receipts, "coverage_nodes_below_floor_PA"), 0)
})

test_that("P^F does not win eta by concentrating the territory", {
  diagnostics <- slice2_results(
    "population-slice2-transport-diagnostics.csv")
  n_nodes <- unique(diagnostics$n_group_nodes)
  expect_equal(length(n_nodes), 1L)

  # A SEVENTH diagnostic, past the six section 7.5 requires. The six catch a
  # transport that wins by DISCARDING territory; none catches one that wins by
  # CONCENTRATING it. A P^F that piled the brain onto a handful of nodes would
  # have full subject coverage, bounded displacement, respectable entropy and
  # P^A's exact sink -- and would pass every one of the six while making the
  # group grid a fiction.
  for (subj in SLICE2_SUBJECT_IDS) {
    row <- diagnostics[diagnostics$subject == subj, ]

    # Total arriving mass is identical under the two transports -- the same
    # statement as the sink control, read from the group-node side.
    expect_gt(row$group_node_mass_total, 0)

    # No group node is starved of mass entirely under either transport.
    expect_equal(row$group_nodes_with_zero_mass_A, 0L)
    expect_equal(row$group_nodes_with_zero_mass_F, 0L)

    # P^F is ALLOWED to concentrate -- choosing a destination is what it is
    # for -- but the effective number of group nodes carrying the territory
    # must stay a large fraction of the grid. P^A sits near 0.94 of the grid
    # and P^F near 0.82; a floor of 0.5 is far below both and would catch a
    # transport that had started collapsing onto attractor nodes.
    expect_gt(row$group_node_effective_A / n_nodes, 0.5)
    expect_gt(row$group_node_effective_F / n_nodes, 0.5)
    expect_lte(row$group_node_effective_F, row$group_node_effective_A)
  }
})
