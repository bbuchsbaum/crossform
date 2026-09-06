if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

test_that("the common statistic reproduces the exact hand geometry", {
  environment <- new.env(parent = globalenv())
  sys.source(test_path("..", "..", "design", "oracles",
    "common-geometry-equivalence.R"), envir = environment)
  court <- environment$common_geometry_equivalence_oracle()

  expect_identical(court$geometry, matrix(c(
    4, 0, 3, 0, 2, 2.5, 3, 2.5, 5
  ), 3L, 3L, byrow = TRUE))
  expect_identical(unname(court$rdm), c(6, 3, 2))
  expect_identical(court$rdm, court$direct_pairs)
  expect_identical(court$contrast_ab, 6)
})

test_that("fixed linear RSA is exactly the RDM-adjoint geometry query", {
  environment <- new.env(parent = globalenv())
  sys.source(test_path("..", "..", "design", "oracles",
    "common-geometry-equivalence.R"), envir = environment)
  court <- environment$common_geometry_equivalence_oracle()

  expect_equal(unname(court$rsa), c(3, 1), tolerance = 1e-14)
  expect_equal(court$rsa_adjoint, court$rsa, tolerance = 1e-14)
  expect_equal(court$centred_rdm_adjoint, court$centred_operator,
    tolerance = 1e-14)
  expect_identical(qr(court$rsa_design)$rank, 2L)
})

test_that("metric and partition assumptions are load-bearing", {
  environment <- new.env(parent = globalenv())
  sys.source(test_path("..", "..", "design", "oracles",
    "common-geometry-equivalence.R"), envir = environment)
  court <- environment$common_geometry_equivalence_oracle()

  euclidean <- environment$oracle_energy(tcrossprod(c(1, -1, 0)),
    court$blocks, court$gamma, diag(2))
  self_paired <- environment$oracle_energy(tcrossprod(c(1, -1, 0)),
    court$blocks, diag(2) / 2, court$metric)
  expect_false(isTRUE(all.equal(euclidean, court$contrast_ab)))
  expect_false(isTRUE(all.equal(self_paired, court$contrast_ab)))
  expect_identical(sum(court$gamma), 1)
  expect_identical(diag(court$gamma), c(0, 0))
})

test_that("the equivalence boundary is explicit", {
  contract <- paste(readLines(test_path("..", "..", "design",
    "common-geometry-equivalence.md"), warn = FALSE), collapse = "\n")
  for (phrase in c(
    "Spearman or Kendall RSA", "rank transforms",
    "data-adaptive model selection", "full column rank",
    "declared independent partitions", "positive definite",
    "fixed *linear* RSA"
  )) {
    expect_match(contract, phrase, fixed = TRUE)
  }
})
