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

parity_exemplar_path <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "exemplars", "rsatoolbox-parity"),
    system.file("exemplars", "rsatoolbox-parity", package = "crossform")
  )
  candidates <- candidates[nzchar(candidates) & dir.exists(candidates)]
  if (!length(candidates)) NA_character_ else normalizePath(candidates[[1L]])
}

parity_manifest_mismatches <- function(manifest, repo) {
  paths <- file.path(repo, manifest$path)
  exists <- file.exists(paths)
  actual_size <- rep(NA_real_, length(paths))
  actual_digest <- rep(NA_character_, length(paths))
  actual_size[exists] <- unname(file.info(paths[exists])$size)
  actual_digest[exists] <- unname(tools::md5sum(paths[exists]))
  manifest$path[
    !exists |
      actual_size != manifest$size_bytes |
      is.na(actual_digest) |
      actual_digest != manifest$digest
  ]
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

test_that("the parity manifest binds every recorded result to its sources", {
  exemplar <- parity_exemplar_path()
  skip_if(is.na(exemplar), "rsatoolbox parity exemplar is not present")
  path <- file.path(exemplar, "results", "parity-manifest.csv")
  skip_if_not(file.exists(path), "rsatoolbox parity manifest is not present")
  manifest <- utils::read.csv(path, stringsAsFactors = FALSE)
  repo <- normalizePath(file.path(exemplar, "..", ".."))

  expect_identical(
    names(manifest),
    c("schema_version", "fixture_id", "hash_algorithm", "role", "path",
      "size_bytes", "digest")
  )
  expect_true(nrow(manifest) > 10L)
  expect_true(all(manifest$schema_version == 1L))
  expect_true(all(manifest$fixture_id == "rsatoolbox-fixed-linear-v1"))
  expect_true(all(manifest$hash_algorithm == "md5"))
  expect_identical(anyDuplicated(manifest$path), 0L)
  expect_setequal(
    unique(manifest$role),
    c("fixture_source", "external_implementation", "comparison_source",
      "environment_lock", "algebraic_claim", "certification_copy",
      "fixture_contract", "recorded_output")
  )
  expect_identical(parity_manifest_mismatches(manifest, repo), character())

  # Exercise the diagnostic itself: a forged expected digest must identify
  # exactly the path whose source/artifact binding is stale.
  forged <- manifest
  forged$digest[[1L]] <- paste0("0", substring(forged$digest[[1L]], 2L))
  if (identical(forged$digest[[1L]], manifest$digest[[1L]])) {
    forged$digest[[1L]] <- paste0("1", substring(forged$digest[[1L]], 2L))
  }
  expect_identical(parity_manifest_mismatches(forged, repo), forged$path[1L])
})

test_that("the recorded external case names its complete convention mapping", {
  exemplar <- parity_exemplar_path()
  skip_if(is.na(exemplar), "rsatoolbox parity exemplar is not present")
  fixture_path <- file.path(exemplar, "results", "fixture-meta.csv")
  external_path <- file.path(exemplar, "results", "rsatoolbox-meta.csv")
  skip_if_not(file.exists(fixture_path) && file.exists(external_path),
              "rsatoolbox parity metadata is not present")

  fixture <- utils::read.csv(fixture_path, stringsAsFactors = FALSE)
  fixture <- stats::setNames(fixture$value, fixture$key)
  external <- utils::read.csv(external_path, stringsAsFactors = FALSE)
  external <- stats::setNames(external$value, external$key)

  expect_identical(external[["python"]], "3.12.11")
  expect_identical(external[["rsatoolbox"]], "0.3.2")
  expect_identical(external[["numpy"]], "2.5.2")
  expect_identical(external[["scipy"]], "1.18.0")
  expect_identical(fixture[["metric_role"]], "fixed_noise_precision")
  expect_identical(
    fixture[["metric_normalization"]],
    "frame_local_divide_by_support_size"
  )
  expect_identical(
    fixture[["effect_centering"]],
    "none_pair_differences_are_zero_sum"
  )
  expect_identical(
    fixture[["partition_scheme"]],
    "uniform_unordered_cross_run_pairs"
  )
  expect_lt(abs(as.numeric(fixture[["partition_weight"]]) - 1 / 6), 1e-15)
  expect_identical(fixture[["pair_order"]], "row_major_upper_triangle")
  expect_identical(fixture[["rsa_objective"]],
                   "fixed_ols_on_vectorized_rdm")
  expect_identical(fixture[["claim_scope"]],
                   "crossnobis_and_fixed_linear_rsa_only")

  agreement <- utils::read.csv(file.path(exemplar, "results", "agreement.csv"),
                               stringsAsFactors = FALSE)
  expect_lte(max(agreement$tolerance), 1e-8)
  expect_true(all(agreement$passes))
})
