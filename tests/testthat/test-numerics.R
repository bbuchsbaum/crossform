test_that("the contract names three different numerical guarantees", {
  contract <- numerical_contract()

  expect_identical(contract$scheduling$guarantee, "bitwise")
  expect_identical(contract$block_partition$guarantee, "tolerance")
  expect_identical(contract$cross_platform$guarantee, "tolerance")
  expect_false(contract$bitwise_across_blocking)
  expect_false(contract$bitwise_across_platforms)
})
test_that("canonical reduction is bitwise invariant to completion order", {
  values <- list(
    a = matrix(c(1, 2, 3, 4), 2),
    b = matrix(c(0.25, 0.5, 0.75, 1), 2),
    c = matrix(c(-2, 1, -1, 2), 2)
  )
  task_id <- c(30L, 10L, 20L)
  oracle <- effectagram:::.canonical_reduce(values, task_id)

  set.seed(20260812)
  for (iteration in seq_len(25)) {
    order <- sample(seq_along(values))
    got <- effectagram:::.canonical_reduce(values[order], task_id[order])
    agreement <- numerical_agreement(got, oracle, guarantee = "scheduling")
    expect_true(agreement$passed)
    expect_identical(got, oracle)
  }
})

test_that("block repartitioning is tolerance-qualified, not bitwise-promised", {
  x <- seq(-1, 1, length.out = 1001)^3
  one_block <- sum(x)
  ten_blocks <- sum(vapply(split(x, ceiling(seq_along(x) / 101)), sum, numeric(1)))

  agreement <- numerical_agreement(
    one_block,
    ten_blocks,
    guarantee = "block_partition",
    contract = numerical_contract(atol = 1e-13, rtol = 1e-12)
  )

  expect_true(agreement$passed)
  expect_identical(agreement$comparison, "tolerance")
})

test_that("combined tolerance handles scale and rejects real disagreement", {
  contract <- numerical_contract(atol = 1e-12, rtol = 1e-8)
  near <- numerical_agreement(
    c(0, 1e8),
    c(5e-13, 1e8 + 0.5),
    guarantee = "cross_platform",
    contract = contract
  )
  far <- numerical_agreement(
    c(0, 1),
    c(1e-6, 1.1),
    guarantee = "cross_platform",
    contract = contract
  )

  expect_true(near$passed)
  expect_false(far$passed)
})

test_that("non-finite comparisons and invalid tolerances are rejected", {
  expect_error(numerical_contract(atol = -1), "nonnegative")
  expect_error(numerical_contract(rtol = Inf), "finite")
  expect_error(numerical_agreement(NA_real_, 0), "finite numeric")
  expect_error(
    effectagram:::.canonical_reduce(list(1, 2), c(1, 1)),
    "unique"
  )
})

test_that("numerical agreement revalidates mutated contracts", {
  contract <- numerical_contract()
  contract$block_partition$guarantee <- "bitwise"

  expect_error(numerical_agreement(1, 1, contract = contract),
    "inconsistent or noncanonical")
})
