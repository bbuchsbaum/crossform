# The recorded rsatoolbox parity result, held as a ratchet.
#
# `exemplars/rsatoolbox-parity` runs R -> Python -> R and writes
# `results/agreement.csv`. Rerunning it needs a version-pinned Python
# environment, so this file does not recompute anything: it asserts that the
# committed record still says what the novelty ledger and the exemplar README
# claim it says. Loosening a tolerance, dropping a comparison, or pasting a
# worse number into the CSV fails here.
#
# The file is skipped when the exemplar results are absent (an installed
# package, or a checkout where the exemplar has not been run).

agreement_path <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "exemplars", "rsatoolbox-parity",
      "results", "agreement.csv"),
    system.file("exemplars", "rsatoolbox-parity", "results", "agreement.csv",
      package = "crossform")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) NA_character_ else normalizePath(candidates[[1L]])
}

test_that("the recorded rsatoolbox agreement is still within tolerance", {
  path <- agreement_path()
  skip_if(is.na(path), "rsatoolbox parity exemplar results are not present")
  agreement <- utils::read.csv(path, stringsAsFactors = FALSE)

  expect_setequal(
    names(agreement),
    c("quantity", "comparator", "n_values", "max_abs_diff", "max_rel_diff",
      "tolerance", "passes", "note")
  )
  expect_true(all(is.finite(agreement$max_abs_diff)))
  expect_true(all(is.finite(agreement$tolerance)))
  expect_true(all(agreement$n_values > 0))

  # The recorded run must be self-consistent: `passes` is not an opinion.
  expect_identical(
    agreement$passes,
    agreement$max_abs_diff <= agreement$tolerance
  )
  expect_true(all(agreement$passes))

  # Tolerances are the ones the ledger and README quote. Widening one to make
  # a worse number pass is exactly what this guard exists to catch.
  expect_lte(max(agreement$tolerance), 1e-8)
  expect_lte(
    max(agreement$tolerance[agreement$comparator !=
      "rsatoolbox ModelWeighted + fit_regress(cosine), rescaled"]),
    1e-10
  )
})

test_that("the recorded agreement still covers both parity claims", {
  path <- agreement_path()
  skip_if(is.na(path), "rsatoolbox parity exemplar results are not present")
  agreement <- utils::read.csv(path, stringsAsFactors = FALSE)

  rdm_row <- agreement[
    agreement$quantity == "crossnobis_rdm" &
      agreement$comparator == "rsatoolbox::calc_rdm_crossnobis", ]
  expect_identical(nrow(rdm_row), 1L)
  expect_gte(rdm_row$n_values, 60)
  # The ledger records 3.8e-15; keep a little headroom for a rerun on other
  # hardware without letting an estimator-scale regression through.
  expect_lt(rdm_row$max_abs_diff, 1e-12)

  oracle_row <- agreement[
    agreement$quantity == "crossnobis_rdm" &
      agreement$comparator == "explicit all-pairs numpy oracle", ]
  expect_identical(nrow(oracle_row), 1L)
  expect_lt(oracle_row$max_abs_diff, 1e-12)

  rsa_row <- agreement[
    agreement$quantity == "linear_rsa_coefficients", ]
  expect_identical(nrow(rsa_row), 1L)
  expect_gte(rsa_row$n_values, 20)
  expect_lt(rsa_row$max_abs_diff, 1e-12)
})
