# The recorded conservative-geometry result on Haxby subject 1, held as a
# ratchet.
#
# `exemplars/haxby2001/07-conservative-geometry.R` runs the detection map, the
# conservative attribution map, the latent layer and the coherence spectrum on
# 300 MB of downloaded fMRI data and writes
# `results/conservative-geometry-receipts.csv`. Rerunning it needs that data,
# so this file recomputes nothing: it asserts that the committed record still
# says what `design/conservative-geometry-contract.md` and the exemplar README
# claim it says. Loosening a tolerance, dropping an identity, deleting a
# refusal receipt, or pasting a worse number into the CSV fails here.
#
# The file is skipped when the exemplar results are absent (an installed
# package, or a checkout where the exemplar has not been run) -- `exemplars/`
# is Rbuildignored, so that is the normal case off a source tree.

conservative_receipts_path <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "exemplars", "haxby2001", "results",
      "conservative-geometry-receipts.csv"),
    system.file("exemplars", "haxby2001", "results",
      "conservative-geometry-receipts.csv", package = "crossform")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) NA_character_ else normalizePath(candidates[[1L]])
}

conservative_receipts <- function() {
  path <- conservative_receipts_path()
  testthat::skip_if(is.na(path),
    "haxby2001 conservative-geometry exemplar results are not present")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

# One value, addressed the way the exemplar recorded it. A missing row is a
# failure rather than an NA, because every assertion below is about a claim
# the README makes.
receipt_value <- function(receipts, quantity, contrast = NULL, group = NULL) {
  rows <- receipts[receipts$quantity == quantity, , drop = FALSE]
  if (!is.null(contrast)) {
    rows <- rows[!is.na(rows$contrast) & rows$contrast == contrast, ,
      drop = FALSE]
  }
  if (!is.null(group)) {
    rows <- rows[!is.na(rows$group) & rows$group == group, , drop = FALSE]
  }
  testthat::expect_equal(nrow(rows), 1L,
    info = paste("receipt", quantity, contrast, group))
  rows$value[[1L]]
}

FAMILY_SCALES <- c("radius-4", "radius-8", "radius-12")
CONTRAST_LABELS <- c("face - house", "animate - inanimate")

test_that("the recorded conservative receipts are internally consistent", {
  receipts <- conservative_receipts()

  expect_setequal(
    names(receipts),
    c("panel", "quantity", "contrast", "group", "value", "tolerance",
      "passes", "note")
  )
  expect_true(all(receipts$panel %in% c("A", "B", "C", "D")))
  expect_true(all(is.finite(receipts$value)))
  expect_true(all(nzchar(receipts$note)))

  # `passes` is not an opinion: it is exactly the tolerance test, and a row
  # without a tolerance does not get to claim one.
  asserted <- !is.na(receipts$tolerance)
  expect_true(any(asserted))
  expect_identical(
    receipts$passes[asserted],
    abs(receipts$value[asserted]) <= receipts$tolerance[asserted]
  )
  expect_true(all(receipts$passes[asserted]))
  expect_true(all(is.na(receipts$passes[!asserted])))

  # Widening a tolerance to make a worse number pass is what this catches.
  expect_lte(max(receipts$tolerance[asserted]), 1e-10)
})

test_that("the conservative frame still conserves the whole-VT budget", {
  receipts <- conservative_receipts()

  # Contract section 2, claim 2, measured on real data: the per-node totals of
  # a column-normalized frame sum exactly to the total under
  # `whole_brain("none")`.
  for (label in CONTRAST_LABELS) {
    identity <- receipt_value(receipts, "conservation_identity", label)
    expect_lte(abs(identity), 1e-10)

    node_sum <- receipt_value(receipts, "conservative_total_sum", label)
    whole <- receipt_value(receipts, "whole_vt_total", label, "whole VT")
    expect_equal(node_sum - whole, identity, tolerance = 1e-12)
    expect_gt(whole, 0)

    # The ledger reading only adds: territories partition the same budget.
    expect_lte(abs(receipt_value(receipts, "ledger_closes", label)), 1e-10)
    territories <- receipts[receipts$quantity == "ledger_total" &
      receipts$contrast == label, ]
    expect_gt(nrow(territories), 1L)
    expect_equal(sum(territories$value), node_sum, tolerance = 1e-9)
  }

  # The frame's own certificate, not just the contracted numbers.
  expect_lte(
    abs(receipt_value(receipts, "frame_column_mass_max_deviation")), 1e-10
  )
})

test_that("the recorded detection map is still not a ledger", {
  receipts <- conservative_receipts()

  # A local frame has no budget (contract section 1.1). The exemplar records
  # both halves of that: the certificate says FALSE, and the sum of the
  # detection map is nowhere near the total it is not an attribution of.
  expect_equal(receipt_value(receipts, "local_frame_conserved"), 0)
  expect_gt(
    receipt_value(receipts, "local_frame_feature_mass_max_deviation"), 1e-6
  )
  for (label in CONTRAST_LABELS) {
    ratio <- receipt_value(receipts, "local_map_sum_over_whole_vt_total",
      label)
    expect_false(isTRUE(all.equal(ratio, 1, tolerance = 1e-3)))
  }

  # And the package refuses to add it up rather than returning the number.
  expect_equal(
    receipt_value(receipts, "contribution_on_local_refused",
      group = "conservative_frame"),
    1
  )
})

test_that("the recorded per-scale energy is fixed by alpha, and refused", {
  receipts <- conservative_receipts()

  for (label in CONTRAST_LABELS) {
    whole <- receipt_value(receipts, "whole_vt_total", label, "whole VT")
    totals <- vapply(FAMILY_SCALES, function(scale) {
      receipt_value(receipts, "per_scale_total", label, scale)
    }, numeric(1))
    alphas <- vapply(FAMILY_SCALES, function(scale) {
      receipt_value(receipts, "family_alpha", label, scale)
    }, numeric(1))

    # Contract section 3.1: the rows of scale s sum to alpha_s * G_Omega
    # whatever the data say. Equal weights therefore give three identical
    # numbers, which is precisely why the panel is not a finding.
    expect_equal(sum(alphas), 1, tolerance = 1e-12)
    expect_equal(unname(totals), unname(alphas * whole), tolerance = 1e-10)
    for (scale in FAMILY_SCALES) {
      expect_lte(
        abs(receipt_value(receipts, "per_scale_total_minus_alpha_times_whole",
          label, scale)),
        1e-10
      )
    }
  }

  # The panel that would draw that column is refused by name.
  expect_equal(
    receipt_value(receipts, "scale_energy_panel_refused",
      group = "scale_energy_panel"),
    1
  )
})

test_that("the recorded coherent share is alpha-invariant, and the energy is not", {
  receipts <- conservative_receipts()

  # Contract section 3.2, measured on real data: rerunning the same radii
  # under a lopsided alpha must leave the share alone and move the energy.
  # Both halves matter -- an invariant share with an invariant energy would
  # mean the second family was never actually reweighted.
  for (label in CONTRAST_LABELS) {
    expect_lte(
      abs(receipt_value(receipts, "coherent_share_alpha_invariance", label)),
      1e-10
    )
    expect_gt(
      receipt_value(receipts, "per_scale_total_alpha_dependence", label), 1
    )
  }
})

test_that("the recorded coherence spectrum is the alpha-invariant finding", {
  receipts <- conservative_receipts()

  for (label in CONTRAST_LABELS) {
    shares <- vapply(FAMILY_SCALES, function(scale) {
      receipt_value(receipts, "coherent_share", label, scale)
    }, numeric(1))
    expect_true(all(shares > 0 & shares < 1))

    # The recorded fact about this subject's VT: the coherent share of the
    # fixed budget falls as the neighbourhood grows. Unlike the energy column
    # above, nothing forces this, which is what makes it worth recording.
    expect_true(all(diff(shares) < 0))

    node_medians <- vapply(FAMILY_SCALES, function(scale) {
      receipt_value(receipts, "node_coherent_share_median", label, scale)
    }, numeric(1))
    expect_true(all(node_medians > 0 & node_medians < 1))
    expect_true(all(diff(node_medians) < 0))

    # Shares are masked where the components are not a nonnegative partition,
    # never clamped, so the valid count is positive and below the node count.
    valid <- vapply(FAMILY_SCALES, function(scale) {
      receipt_value(receipts, "node_coherent_share_n_valid", label, scale)
    }, numeric(1))
    expect_true(all(valid > 0))
  }

  # Every member of the family is column normalized on its own, which the
  # stacked family conserving does not imply (contract 3.1, precondition 1).
  for (scale in FAMILY_SCALES) {
    expect_lte(
      abs(receipt_value(receipts,
        "family_member_column_mass_max_deviation", group = scale)),
      1e-10
    )
  }
})

test_that("the recorded latent layer reports what its projection cost", {
  receipts <- conservative_receipts()

  # Contract section 6: C(k) and n_eff exist only on the latent layer, and the
  # projection is never silent about the mass it moved.
  n_eff <- receipt_value(receipts, "latent_n_eff", group = "whole VT")
  expect_gt(n_eff, 1)

  cumulative <- receipts[receipts$quantity == "latent_cumulative_contribution",
    , drop = FALSE]
  expect_gt(nrow(cumulative), 1L)
  order_k <- order(as.integer(sub("^k=", "", cumulative$group)))
  curve <- cumulative$value[order_k]
  expect_true(all(diff(curve) >= -1e-12))
  expect_equal(curve[[length(curve)]], 1, tolerance = 1e-10)
  expect_lt(n_eff, length(curve))

  # The empirical form really does clip: this is not a formality on real data.
  expect_gte(receipt_value(receipts, "latent_clipped_measurements",
    group = "whole VT"), 1)
  expect_gt(receipt_value(receipts, "latent_moved_mass", group = "whole VT"), 0)
  moved_share <- receipt_value(receipts, "latent_moved_share",
    group = "whole VT")
  expect_gt(moved_share, 0)
  expect_lt(moved_share, 0.5)
  expect_gte(
    receipt_value(receipts, "signed_spectrum_negative_roots",
      group = "whole VT"),
    1
  )
})
