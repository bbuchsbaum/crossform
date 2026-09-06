if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

common_oracle_environment <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(test_path("..", "..", "design", "oracles",
    "common-geometry-equivalence.R"), envir = environment)
  environment
}

test_that("explicit matrices and direct indices agree on rectangular forms", {
  oracle <- common_oracle_environment()
  left <- list(
    l1 = matrix(c(1, 2, -1, 0), 2L, 2L, byrow = TRUE),
    l2 = matrix(c(0, 1, 2, -2), 2L, 2L, byrow = TRUE),
    l3 = matrix(c(3, -1, 1, 1), 2L, 2L, byrow = TRUE)
  )
  right <- list(
    r1 = matrix(c(1, 0, 0, 1, 2, 1), 3L, 2L, byrow = TRUE),
    r2 = matrix(c(-1, 1, 1, 2, 0, 3), 3L, 2L, byrow = TRUE)
  )
  gamma <- matrix(c(0.2, -0.1, 0.3, 0.4, -0.2, 0.4), 3L, 2L,
    dimnames = list(names(left), names(right)))
  metric <- matrix(c(2, 0.25, 0.5, 1), 2L, 2L, byrow = TRUE)
  operator <- matrix(c(1, -1, 0.5, 0, 2, -0.25), 2L, 3L,
    byrow = TRUE)

  explicit <- oracle$oracle_bilinear(
    operator, left, right, gamma, metric
  )
  direct <- oracle$oracle_bilinear_direct(
    operator, left, right, gamma, metric
  )
  comparison <- oracle$oracle_numeric_comparison(
    explicit, direct, matrices = list(metric), operations = 72L
  )
  expect_true(comparison$pass)
  expect_lte(comparison$error, comparison$bound)
  expect_gt(comparison$condition, 1)
})

test_that("random square fixtures agree within conditioning-aware bounds", {
  oracle <- common_oracle_environment()
  for (seed in 201:205) {
    set.seed(seed)
    partitions <- paste0("run", 1:4)
    blocks <- stats::setNames(lapply(partitions, function(run) {
      matrix(stats::rnorm(12), 4L, 3L)
    }), partitions)
    factor <- matrix(stats::rnorm(9), 3L, 3L)
    metric <- crossprod(factor) + diag(0.25, 3L)
    gamma <- matrix(1, 4L, 4L,
      dimnames = list(partitions, partitions))
    diag(gamma) <- 0
    gamma <- gamma / sum(gamma)
    contrast <- c(1, -1, 0.5, -0.5)
    operator <- tcrossprod(contrast)

    explicit <- oracle$oracle_bilinear(
      operator, blocks, blocks, gamma, metric
    )
    direct <- oracle$oracle_bilinear_direct(
      operator, blocks, blocks, gamma, metric
    )
    court <- oracle$oracle_numeric_comparison(
      explicit, direct, matrices = list(metric), operations = 576L
    )
    expect_true(court$pass, info = paste("seed", seed))
    expect_true(is.finite(court$condition))
  }
})

test_that("centering is exactly the common-shift boundary", {
  oracle <- common_oracle_environment()
  set.seed(310)
  blocks <- lapply(1:3, function(run) matrix(stats::rnorm(12), 4L, 3L))
  gamma <- matrix(1 / 6, 3L, 3L)
  diag(gamma) <- 0
  metric <- diag(c(1, 2, 3))
  shift <- c(0.7, -1.2, 0.4)
  shifted <- lapply(blocks, function(value) value +
    matrix(shift, nrow(value), length(shift), byrow = TRUE))

  centred <- c(1, -1, 0, 0)
  uncentred <- c(1, 0, 0, 0)
  expect_equal(
    oracle$oracle_energy(tcrossprod(centred), blocks, gamma, metric),
    oracle$oracle_energy(tcrossprod(centred), shifted, gamma, metric),
    tolerance = 1e-13
  )
  expect_false(isTRUE(all.equal(
    oracle$oracle_energy(tcrossprod(uncentred), blocks, gamma, metric),
    oracle$oracle_energy(tcrossprod(uncentred), shifted, gamma, metric),
    tolerance = 1e-8
  )))
})

test_that("H K Gamma and partition-label perturbations are detected", {
  oracle <- common_oracle_environment()
  court <- oracle$common_geometry_equivalence_oracle()
  operator <- tcrossprod(c(1, -1, 0))
  reference <- oracle$oracle_energy(
    operator, court$blocks, court$gamma, court$metric
  )
  detect <- function(operator, gamma, metric) {
    candidate <- oracle$oracle_energy(
      operator, court$blocks, gamma, metric
    )
    !oracle$oracle_numeric_comparison(reference, candidate)$pass
  }

  changed_h <- operator
  changed_h[1, 1] <- changed_h[1, 1] + 0.01
  changed_k <- court$metric
  changed_k[1, 1] <- changed_k[1, 1] + 0.01
  changed_gamma <- court$gamma
  changed_gamma[1, 2] <- changed_gamma[1, 2] + 0.01
  expect_true(detect(changed_h, court$gamma, court$metric))
  expect_true(detect(operator, court$gamma, changed_k))
  expect_true(detect(operator, changed_gamma, court$metric))

  blocks <- stats::setNames(court$blocks, c("run1", "run2"))
  gamma <- court$gamma
  dimnames(gamma) <- list(c("run1", "wrong"), c("run1", "run2"))
  expect_error(oracle$oracle_bilinear(
    operator, blocks, blocks, gamma, court$metric
  ), "partition labels do not bind")
})

test_that("the oracle remains independent of production execution helpers", {
  source <- paste(readLines(test_path("..", "..", "design", "oracles",
    "common-geometry-equivalence.R"), warn = FALSE), collapse = "\n")
  for (forbidden in c(
    ".run_geometry_compiler", "evaluate_geometry", "bilinear_query",
    ".pair_difference_query", ".thin_qr_coefficient_map"
  )) {
    expect_false(grepl(forbidden, source, fixed = TRUE), info = forbidden)
  }
})
