# The six-subject Haxby conservative ledger, held as a ratchet.
#
# `exemplars/haxby2001/09-six-subject-conservative.R` runs the conservative
# attribution arm of `07-conservative-geometry.R` over all six Haxby 2001
# subjects and commits
# `exemplars/haxby2001/results/six-subject-conservative-receipts.csv`. This
# file does not recompute any of it -- that needs ~1.8 GB of downloaded raw
# data which is deliberately gitignored. It asserts that the committed record
# still says what the exemplar README claims: six subjects, and every
# conservation identity inside its acceptance.
#
# Deleting the hard subjects, loosening the tolerance column, or pasting a
# worse identity into the CSV fails here.
#
# The file skips cleanly when the exemplar results are absent (an installed
# package, or a checkout where the exemplar has not been run).

SIX_SUBJECT_IDS <- paste0("subj", 1:6)

six_subject_receipts <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "exemplars", "haxby2001", "results",
                        "six-subject-conservative-receipts.csv"),
    system.file("exemplars", "haxby2001", "results",
                "six-subject-conservative-receipts.csv", package = "crossform")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  skip_if(!length(candidates),
          "the six-subject Haxby exemplar results are not present")
  utils::read.csv(normalizePath(candidates[[1L]]), stringsAsFactors = FALSE)
}

test_that("the six-subject Haxby receipts cover all six subjects", {
  receipts <- six_subject_receipts()

  expect_setequal(
    names(receipts),
    c("subject", "quantity", "value", "tolerance", "passes", "note")
  )
  expect_setequal(unique(receipts$subject), SIX_SUBJECT_IDS)

  # Every subject carries the same battery, so a subject cannot be kept in the
  # table by contributing only the rows that happen to be easy.
  per_subject <- split(receipts$quantity, receipts$subject)
  reference <- sort(per_subject[[SIX_SUBJECT_IDS[[1L]]]])
  for (subj in SIX_SUBJECT_IDS) {
    expect_identical(sort(per_subject[[subj]]), reference,
                     info = paste("quantities recorded for", subj))
  }
})

test_that("every recorded six-subject conservation identity passes", {
  receipts <- six_subject_receipts()
  asserted <- receipts[!is.na(receipts$tolerance), ]

  expect_gt(nrow(asserted), 0L)
  expect_true(all(asserted$passes))
  # `passes` must be exactly the comparison it claims to be, not a column of
  # TRUEs someone typed.
  expect_identical(asserted$passes, abs(asserted$value) <= asserted$tolerance)
  expect_true(all(asserted$tolerance <= 1e-10))
  expect_true(all(is.finite(asserted$value)))
})

test_that("the conservation identity is recorded for every subject", {
  receipts <- six_subject_receipts()

  # THE acceptance of ticket E10: sum_x total_x - whole_VT_total, per subject.
  identity <- receipts[receipts$quantity == "conservation_identity", ]
  expect_setequal(identity$subject, SIX_SUBJECT_IDS)
  expect_true(all(identity$tolerance == 1e-10))
  expect_true(all(abs(identity$value) <= 1e-10))

  # The identity is only meaningful against a nonzero budget: an all-zero
  # ledger would satisfy it trivially.
  budget <- receipts[receipts$quantity == "whole_vt_total", ]
  expect_setequal(budget$subject, SIX_SUBJECT_IDS)
  expect_true(all(abs(budget$value) > 1))

  totals <- receipts[receipts$quantity == "conservative_total_sum", ]
  merged <- merge(budget[, c("subject", "value")],
                  totals[, c("subject", "value")], by = "subject",
                  suffixes = c("_budget", "_sum"))
  expect_equal(nrow(merged), 6L)
  expect_equal(merged$value_sum, merged$value_budget, tolerance = 1e-10)

  # The frame each subject's ledger rests on must actually partition its
  # voxels; conservation of the total follows from that column mass.
  column_mass <- receipts[
    receipts$quantity == "frame_column_mass_max_deviation", ]
  expect_setequal(column_mass$subject, SIX_SUBJECT_IDS)
  expect_true(all(abs(column_mass$value) <= 1e-10))
})

test_that("the six-subject descriptive readouts stay well formed", {
  receipts <- six_subject_receipts()
  value_of <- function(quantity) {
    rows <- receipts[receipts$quantity == quantity, ]
    stats::setNames(rows$value, rows$subject)[SIX_SUBJECT_IDS]
  }

  voxels <- value_of("vt_voxels_analyzed")
  expect_true(all(voxels > 100))
  # One conservative searchlight per analyzed voxel: centres are the domain.
  expect_equal(unname(value_of("n_searchlights")), unname(voxels))

  runs <- value_of("n_runs")
  expect_true(all(runs >= 2))
  expect_true(all(runs <= 12))
  # Every off-diagonal ordered run pair, so the pair count is choose(runs, 2).
  expect_equal(unname(value_of("n_run_pairs")), unname(choose(runs, 2)))

  # subj5's ninth run carries no labelled volumes in this distribution and is
  # dropped, so it is the one subject with fewer than twelve. If a future
  # rerun silently restores it, that is a change of preparation, not a fix.
  expect_equal(unname(runs[["subj5"]]), 11)
  expect_true(all(runs[setdiff(SIX_SUBJECT_IDS, "subj5")] == 12))

  # The coherent share is a fraction of a frame-relative mass, so it is only
  # interpretable inside [0, 1]; a value outside it means the decomposition
  # stopped being a decomposition.
  share <- value_of("coherent_share_base_radius")
  expect_true(all(share > 0 & share < 1))

  # Every quirk row carries a sentence, whether or not it counts anything.
  quirks <- receipts[receipts$quantity == "n_quirks", ]
  expect_setequal(quirks$subject, SIX_SUBJECT_IDS)
  expect_true(all(nzchar(quirks$note)))
  expect_true(all(quirks$value >= 0))
})
