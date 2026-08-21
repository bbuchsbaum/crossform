# Population slice 1 on six Haxby subjects, held as a ratchet.
#
# `exemplars/haxby2001/10-population-slice1.R` carries six subjects' native VT
# territories onto three shared group ROI nodes with a declared transport,
# fits one group model over them, and commits
# `exemplars/haxby2001/results/population-slice1-receipts.csv`. This file does
# not recompute any of it -- that needs ~1.8 GB of downloaded raw data which is
# deliberately gitignored. It asserts that the committed record still says what
# the exemplar README claims: six subjects, three group nodes, and every
# population-form-v1 identity inside its acceptance.
#
# What fails here: dropping a subject, dropping an identity row, loosening the
# tolerance column, pasting a worse number in, quietly clipping the indefinite
# heterogeneity spectrum, or emitting a coherence fraction on the signed
# estimation layer.
#
# The file skips cleanly when the exemplar results are absent (an installed
# package, or a checkout where the exemplar has not been run).

POPULATION_SUBJECT_IDS <- paste0("subj", 1:6)
POPULATION_GROUP_NODES <- c("face-territory", "house-territory", "other-VT")
POPULATION_QUERIES <- c("face-house", "animate-inanimate")

population_slice1_receipts <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "exemplars", "haxby2001", "results",
                        "population-slice1-receipts.csv"),
    system.file("exemplars", "haxby2001", "results",
                "population-slice1-receipts.csv", package = "crossform")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  skip_if(!length(candidates),
          "the Haxby population slice 1 exemplar results are not present")
  utils::read.csv(normalizePath(candidates[[1L]]), stringsAsFactors = FALSE)
}

population_value <- function(receipts, quantity, subject = "group") {
  rows <- receipts[receipts$quantity == quantity &
                     receipts$subject == subject, ]
  expect_equal(nrow(rows), 1L,
               info = paste("receipt", quantity, "for", subject))
  rows$value[[1L]]
}

test_that("the population slice 1 receipts cover six subjects and the group", {
  receipts <- population_slice1_receipts()

  expect_setequal(
    names(receipts),
    c("subject", "quantity", "value", "tolerance", "passes", "note")
  )
  expect_setequal(unique(receipts$subject),
                  c(POPULATION_SUBJECT_IDS, "group"))
  expect_equal(population_value(receipts, "n_subjects"),
               length(POPULATION_SUBJECT_IDS))
  expect_equal(population_value(receipts, "n_group_nodes"),
               length(POPULATION_GROUP_NODES))
  # N - rank(X) for the group model ~ 1.
  expect_equal(population_value(receipts, "residual_df"),
               length(POPULATION_SUBJECT_IDS) - 1L)

  # Every participant carries the same battery, so a subject cannot be kept in
  # the table by contributing only the rows that happen to be easy.
  per_subject <- receipts[receipts$subject %in% POPULATION_SUBJECT_IDS, ]
  split_quantities <- split(per_subject$quantity, per_subject$subject)
  reference <- sort(split_quantities[[POPULATION_SUBJECT_IDS[[1L]]]])
  for (subj in POPULATION_SUBJECT_IDS) {
    expect_identical(sort(split_quantities[[subj]]), reference,
                     info = paste("quantities recorded for", subj))
  }
})

test_that("every recorded population identity passes", {
  receipts <- population_slice1_receipts()
  asserted <- receipts[!is.na(receipts$tolerance), ]

  expect_gt(nrow(asserted), 0L)
  expect_true(all(asserted$passes))
  # `passes` must be exactly the comparison it claims to be, not a column of
  # TRUEs someone typed.
  expect_identical(asserted$passes, abs(asserted$value) <= asserted$tolerance)
  expect_true(all(asserted$tolerance <= 1e-10))
  expect_true(all(is.finite(asserted$value)))
})

test_that("the transport identities hold for every participant", {
  receipts <- population_slice1_receipts()

  # population-form-v1 claim 2: the transported total, sink included, is the
  # native total. The tolerance is the contract's, relative to the ledger's L1
  # norm rather than to a signed total that may sit near zero.
  budget <- receipts[receipts$quantity == "transport_budget_preservation", ]
  expect_setequal(budget$subject, POPULATION_SUBJECT_IDS)
  expect_true(all(budget$tolerance == 1e-12))
  expect_true(all(abs(budget$value) <= 1e-12))

  # The executor-plumbing check: the group fit's response really is
  # `transport_values()` of the subject's own `contrast_energy()` ledger. Both
  # sides query before transporting, so this is exact by construction -- and
  # the test pins it AT exactly zero rather than merely inside 1e-12, because
  # anything else would mean the executor started doing something extra.
  plumbing <- receipts[receipts$quantity == "transport_executor_agreement", ]
  expect_setequal(plumbing$subject, POPULATION_SUBJECT_IDS)
  expect_true(all(plumbing$tolerance == 1e-12))
  expect_true(all(plumbing$value == 0))

  # Density conserves nothing. This is asserted as a LARGE inequality and
  # never as an equality -- but the threshold is only meaningful against the
  # territory sizes, because the gap is about 1 - 1/(mean territory size). The
  # honest pin is that the arm is genuinely a density arm: its column sum is
  # far from the native total AND the per-voxel readings differ from the
  # budget ones.
  gap <- receipts[receipts$quantity == "density_nonconservation_relative", ]
  expect_setequal(gap$subject, POPULATION_SUBJECT_IDS)
  expect_true(all(is.na(gap$tolerance)))
  expect_true(all(gap$value > 0.9))

  node_sum <- receipts[
    receipts$quantity == "transport_node_sum_plus_sink", ]
  expect_setequal(node_sum$subject, POPULATION_SUBJECT_IDS)
  expect_true(all(abs(node_sum$value) <= 1e-12))

  # The region transport has full coverage, so the sink is materialized and
  # exactly empty. It is still a column: the contract refuses a transport that
  # drops it (section 1.1).
  for (quantity in c("transport_sink_budget", "transport_sink_territory")) {
    rows <- receipts[receipts$quantity == quantity, ]
    expect_setequal(rows$subject, POPULATION_SUBJECT_IDS)
    expect_true(all(abs(rows$value) <= 1e-12), info = quantity)
  }

  # The region frame must actually partition each subject's VT, and its three
  # territory totals must sum to the whole-VT budget 09 reports. Without this,
  # everything above would be exact arithmetic on the wrong ledger.
  for (quantity in c("frame_column_mass_max_deviation",
                     "region_ledger_vs_whole_vt")) {
    rows <- receipts[receipts$quantity == quantity, ]
    expect_setequal(rows$subject, POPULATION_SUBJECT_IDS)
    expect_true(all(abs(rows$value) <= 1e-10), info = quantity)
  }

  # An identity on an all-zero ledger is free. The budgets are large.
  budgets <- receipts[
    receipts$quantity == "whole_vt_total: face-house", ]
  expect_setequal(budgets$subject, POPULATION_SUBJECT_IDS)
  expect_true(all(abs(budgets$value) > 1))
})

test_that("the group-level identities hold", {
  receipts <- population_slice1_receipts()

  for (quantity in c("budget_preservation_worst", "commutation_claim3",
                     "commutation_claim3_density", "executor_agreement_worst",
                     "sum_over_nodes_identity",
                     "ledger_identity_coherent_plus_configuration",
                     "rdm_edge_matches_query_bank",
                     "budget_agrees_with_script_09")) {
    rows <- receipts[receipts$quantity == quantity &
                       receipts$subject == "group", ]
    expect_equal(nrow(rows), 1L, info = quantity)
    expect_true(is.finite(rows$tolerance[[1L]]), info = quantity)
    expect_lte(abs(rows$value[[1L]]), rows$tolerance[[1L]])
  }

  # THE acceptance of ticket E11. `commutation_claim3` compares two genuinely
  # different orders -- estimate_population() contracts the contrast into one
  # packed operator and transports one number per node, while
  # materialize_population() transports all 36 packed coordinates and rdm()
  # contracts the (face, house) edge out of them afterwards.
  #
  # It must NOT be exactly zero. An exact zero would mean the two sides had
  # collapsed onto the same arithmetic, which is precisely the failure mode
  # that made an earlier draft of this receipt vacuous: comparing
  # `transport_values(P, contrast_energy(...))` against the query-bank fit
  # compares a call with itself, because contrast_energy() lowers to exactly
  # the operator the bank carries. A real order reversal leaves rounding.
  claim3 <- population_value(receipts, "commutation_claim3")
  expect_gt(abs(claim3), 0)
  expect_lt(abs(claim3), 1e-12)
  # The readout-named row is the same number.
  expect_equal(population_value(receipts, "rdm_edge_matches_query_bank"),
               claim3)
  # Density is a declared row-mass ratio, so it is still a fixed linear map
  # and the same order reversal closes for it too.
  expect_lt(abs(population_value(receipts, "commutation_claim3_density")),
            1e-12)

  # Every subject reaches every group node, so no group node is a group
  # estimate resting on a handful of participants (contract 7.5), and none
  # falls below the declared floor.
  expect_equal(population_value(receipts, "group_node_min_subject_coverage"),
               length(POPULATION_SUBJECT_IDS))
  expect_equal(population_value(receipts,
                                "group_nodes_below_coverage_floor"), 0)
})

test_that("the group readouts are complete and well formed", {
  receipts <- population_slice1_receipts()

  for (query in POPULATION_QUERIES) {
    for (node in POPULATION_GROUP_NODES) {
      estimate <- population_value(
        receipts, paste0("group_contrast: ", node, ": ", query))
      se <- population_value(
        receipts, paste0("group_contrast_se: ", node, ": ", query))
      expect_true(is.finite(estimate))
      # A zero SE would make the t and the interval NA; the group nodes all
      # carry real between-subject spread.
      expect_gt(se, 0)

      # The interval the README prints must bracket the estimate and be the
      # symmetric t interval the SE and df imply -- not a wider band someone
      # pasted in.
      lower <- population_value(
        receipts, paste0("group_contrast_lower: ", node, ": ", query))
      upper <- population_value(
        receipts, paste0("group_contrast_upper: ", node, ": ", query))
      expect_lt(lower, estimate)
      expect_gt(upper, estimate)
      half <- stats::qt(0.975, df = length(POPULATION_SUBJECT_IDS) - 1L) * se
      expect_equal(upper - lower, 2 * half, tolerance = 1e-8)
      expect_equal((lower + upper) / 2, estimate, tolerance = 1e-8)

      # The two transported component ledgers add to the transported total.
      # This is contract 8.1's positive statement, held row by row rather than
      # only as the aggregate identity above.
      coherent <- population_value(
        receipts, paste0("native_coherent_ledger: ", node, ": ", query))
      configuration <- population_value(
        receipts, paste0("native_configuration_ledger: ", node, ": ", query))
      expect_equal(coherent + configuration, estimate, tolerance = 1e-10)

      prevalence <- population_value(
        receipts, paste0("prevalence_positive: ", node, ": ", query))
      expect_gte(prevalence, 0)
      expect_lte(prevalence, 1)
      # A prevalence is a count over six subjects, so it lands on a sixth.
      expect_equal(prevalence * length(POPULATION_SUBJECT_IDS),
                   round(prevalence * length(POPULATION_SUBJECT_IDS)),
                   tolerance = 1e-9)
    }
  }

  # Prevalence is read on the latent descriptive layer and carries the
  # fraction a pure-noise cell reports. Without that reference a count of
  # participants reads as if 0 were its null, which it is not: thresholding a
  # signed crossvalidated estimate keeps the sign and discards the magnitude.
  expect_equal(population_value(receipts, "prevalence_noise_reference"), 0.5)

  # The alignment reading, one per group node: agreement in direction with the
  # LEAVE-ONE-OUT mean of the other participants.
  for (node in POPULATION_GROUP_NODES) {
    alignment <- population_value(receipts,
      paste0("prevalence_alignment: ", node))
    expect_gte(alignment, 0)
    expect_lte(alignment, 1)
  }
  # The alignment inner product is Euclidean in the query readout, not
  # Frobenius in the forms, and how far the bank's Gram sits from the identity
  # is the number that says how unequally the queries are weighted. It must be
  # recorded; it is not required to be small.
  expect_true(is.finite(
    population_value(receipts, "prevalence_readout_gram_deviation")))

  # The group RDM edge at each node is the face - house contrast, and the
  # complete form was read over all 28 condition pairs.
  expect_equal(population_value(receipts, "group_rdm_pairs"),
               choose(8L, 2L))
  for (node in POPULATION_GROUP_NODES) {
    edge <- population_value(receipts,
      paste0("group_rdm_edge: ", node, ": face-house"))
    contrast <- population_value(receipts,
      paste0("group_contrast: ", node, ": face-house"))
    expect_equal(edge, contrast, tolerance = 1e-10)
  }

  # No coherence fraction is emitted anywhere. A fraction of signed estimates
  # is a latent-layer object (contract 8.1); a receipt carrying one on the
  # estimation layer is the exact error the transported names exist to
  # prevent, and it would not announce itself. The pattern is deliberately
  # loose because the offending name is not known in advance.
  # ("ratio" is anchored to a separator: "configuRATIOn" is not a fraction,
  # and an unanchored pattern would flag the ledger identity row.)
  coherence_words <- grepl("coheren", receipts$quantity, ignore.case = TRUE)
  fraction_words <- grepl(
    paste0("fraction|share|percent|pct|proportion|over_total|_of_total",
           "|(^|_)ratio($|_)"),
    receipts$quantity, ignore.case = TRUE)
  expect_false(any(coherence_words & fraction_words))
  # The pattern must be capable of firing, or it is decoration.
  expect_true(any(grepl("(^|_)ratio($|_)", c("coherence_ratio", "ratio_x"))))
  expect_false(any(grepl("(^|_)ratio($|_)", "configuration")))

  # And the transported components are named for what they are. Contract 8.1
  # forbids the bare words on a transported result: the coherent ledger is
  # native-node coherence carried to a group location, not a group-node
  # common mode, and a `group_coherent_*` receipt would assert the reading
  # the name exists to prevent.
  expect_true(all(grepl("^native_coherent_ledger: ",
    grep("coherent_ledger", receipts$quantity, value = TRUE))))
  expect_false(any(grepl("^group_coherent|^group_configuration|^coherent:|^configuration:",
                         receipts$quantity)))
})

test_that("the committed group numbers are the ones the README reports", {
  # The point of a ratchet is that a degraded rerun fails it. Structure and
  # tolerances alone do not do that: every number below is load-bearing in the
  # README's tables, so each is pinned to the value that run produced. A
  # genuine change of estimand should fail here and be re-recorded knowingly.
  receipts <- population_slice1_receipts()

  expected <- list(
    "group_contrast: face-territory: face-house" = 10.9006,
    "group_contrast: house-territory: face-house" = 117.0872,
    "group_contrast: other-VT: face-house" = 49.7202,
    "group_contrast: face-territory: animate-inanimate" = 7.9417,
    "group_contrast: house-territory: animate-inanimate" = 12.6870,
    "group_contrast: other-VT: animate-inanimate" = 18.7523,
    "group_contrast_se: face-territory: face-house" = 4.0372,
    "group_contrast_se: house-territory: face-house" = 27.7427,
    "group_contrast_se: other-VT: face-house" = 10.7514,
    # The density arm reorders face-territory above other-VT per voxel. That
    # reordering is the README's point, so both numbers are pinned.
    "group_density: face-territory: face-house" = 0.4310,
    "group_density: house-territory: face-house" = 1.0475,
    "group_density: other-VT: face-house" = 0.1467,
    "native_coherent_ledger: other-VT: face-house" = 12.2843,
    "native_configuration_ledger: other-VT: face-house" = 37.4359,
    "heterogeneity_eigenvalue: mode1" = 7373.75,
    "heterogeneity_eigenvalue: mode2" = 2390.13,
    "heterogeneity_n_eff" = 2.5876
  )
  for (quantity in names(expected)) {
    expect_equal(population_value(receipts, quantity), expected[[quantity]],
                 tolerance = 1e-3, info = quantity)
  }

  # Per-voxel, face-territory carries more face - house evidence than generic
  # VT while carrying far less budget. That inversion is the README's whole
  # argument for reporting both semantics, so it is asserted as an ordering
  # and not only as two numbers.
  budget <- vapply(POPULATION_GROUP_NODES, function(node)
    population_value(receipts, paste0("group_contrast: ", node,
                                      ": face-house")), numeric(1))
  density <- vapply(POPULATION_GROUP_NODES, function(node)
    population_value(receipts, paste0("group_density: ", node,
                                      ": face-house")), numeric(1))
  expect_lt(budget[["face-territory"]], budget[["other-VT"]])
  expect_gt(density[["face-territory"]], density[["other-VT"]])
  # house-territory leads under both, which is what the argmax receipt says.
  expect_equal(population_value(receipts,
                                "budget_argmax_matches_density_argmax"), 1)
  expect_identical(names(which.max(budget)), "house-territory")
  expect_identical(names(which.max(density)), "house-territory")

  # The face-wins overlap rule never fires on this distribution, which is why
  # the README can say the tie-break is immaterial here.
  overlap <- receipts[receipts$quantity == "roi_overlap_voxels", ]
  expect_setequal(overlap$subject, POPULATION_SUBJECT_IDS)
  expect_true(all(overlap$value == 0))

  # The subject loadings the README names.
  expect_equal(population_value(receipts, "heterogeneity_loading: mode1",
                                "subj4"), 0.8184, tolerance = 1e-3)
  expect_equal(population_value(receipts, "heterogeneity_loading: mode2",
                                "subj6"), 0.7906, tolerance = 1e-3)

  # The identities are pinned at their own scale, not merely inside the loose
  # 1e-10 the generic test allows: a rerun that degraded them by three orders
  # of magnitude would still be "passing" without this.
  expect_lt(abs(population_value(receipts, "sum_over_nodes_identity")), 1e-12)
  expect_lt(max(abs(receipts$value[
    receipts$quantity == "region_ledger_vs_whole_vt"])), 1e-11)
  expect_lt(abs(population_value(receipts,
                                 "budget_agrees_with_script_09")), 1e-11)
})

test_that("the heterogeneity spectrum is reported signed and unclipped", {
  receipts <- population_slice1_receipts()

  # The Gram is N x N, so it has N eigenvalues; its RANK is N - rank(X).
  n_modes <- population_value(receipts, "heterogeneity_gram_size")
  expect_equal(n_modes, length(POPULATION_SUBJECT_IDS))

  eigenvalues <- vapply(seq_len(n_modes), function(k)
    population_value(receipts, paste0("heterogeneity_eigenvalue: mode", k)),
    numeric(1))
  expect_true(all(is.finite(eigenvalues)))
  expect_true(all(diff(eigenvalues) <= 1e-6))
  expect_gt(eigenvalues[[1L]], 0)

  # The Gram is centered across subjects, so rank is at most N - rank(X) = 5
  # and the last eigenvalue sits at numerical zero. The cross-fitted Gram is
  # indefinite by construction, so the count of negative modes is recorded and
  # not asserted to be zero.
  expect_lt(abs(eigenvalues[[n_modes]]), 1e-6 * eigenvalues[[1L]])
  negative <- population_value(receipts, "heterogeneity_negative_modes")
  expect_gte(negative, 0)
  expect_lte(negative, n_modes)

  # n_eff and the moved mass are nonnegative functionals and therefore live on
  # the projected latent layer; both must be present, and the moved mass is
  # what says how much the projection cost.
  n_eff <- population_value(receipts, "heterogeneity_n_eff")
  expect_gte(n_eff, 1)
  expect_lte(n_eff, n_modes)
  expect_gte(population_value(receipts, "heterogeneity_moved_mass"), 0)

  # Subject loadings on the two leading modes, one row per participant.
  for (k in 1:2) {
    rows <- receipts[
      receipts$quantity == paste0("heterogeneity_loading: mode", k), ]
    expect_setequal(rows$subject, POPULATION_SUBJECT_IDS)
    expect_true(all(abs(rows$value) <= 1 + 1e-9))
    # A mode is a unit vector over subjects, so its loadings cannot all be
    # negligible.
    expect_gt(sum(rows$value^2), 0.5)
  }
})

test_that("the population descriptive readouts stay well formed", {
  receipts <- population_slice1_receipts()
  value_of <- function(quantity) {
    rows <- receipts[receipts$quantity == quantity &
                       receipts$subject %in% POPULATION_SUBJECT_IDS, ]
    stats::setNames(rows$value, rows$subject)[POPULATION_SUBJECT_IDS]
  }

  # 09's domain, subject for subject: the two scripts must be reading the same
  # voxels or their budgets are not comparable.
  voxels <- value_of("vt_voxels_analyzed")
  expect_true(all(voxels > 100))

  # Three territories partition VT, so their voxel counts add to the domain.
  counts <- vapply(POPULATION_GROUP_NODES, function(node)
    value_of(paste0("native_region_voxels: ", node)), numeric(6))
  expect_equal(unname(rowSums(counts)), unname(voxels))
  expect_true(all(counts > 0))

  # One native node per nonempty territory, and therefore one transport row.
  expect_equal(unname(value_of("native_nodes")),
               unname(rowSums(counts > 0)))

  runs <- value_of("n_runs")
  expect_equal(unname(value_of("n_run_pairs")), unname(choose(runs, 2)))
  # 08's preparation drops subj5's ninth run, which carries no labelled
  # volumes in this distribution. A rerun that silently restores it is a
  # change of preparation, not a fix.
  expect_equal(unname(runs[["subj5"]]), 11)
  expect_true(all(runs[setdiff(POPULATION_SUBJECT_IDS, "subj5")] == 12))
  # The cross-fitted Gram needs four partitions per participant.
  expect_true(all(runs >= 4))

  # Every quirk row carries a sentence, whether or not it counts anything.
  quirks <- receipts[receipts$quantity == "n_quirks", ]
  expect_equal(nrow(quirks), 1L)
  expect_true(nzchar(quirks$note[[1L]]))
})
