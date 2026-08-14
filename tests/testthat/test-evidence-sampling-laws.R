test_that("Diedrichsen Eq. 10 is exactly the all-unordered-pairs estimator", {
  set.seed(81311)
  observed <- vapply(2:12, function(partitions) {
    values <- matrix(rnorm(partitions * 17), partitions, 17)
    abs(
      sampling_oracle_leave_one_out(values) -
        sampling_oracle_all_pairs(values)
    )
  }, numeric(1))

  expect_lt(max(observed), 5e-16)
})

test_that("Eq. 13 equals an independent endpoint-overlap enumeration", {
  set.seed(81312)
  patterns <- matrix(rnorm(5 * 13), 5, 13)
  sigma_k_raw <- matrix(rnorm(5 * 5), 5, 5)
  sigma_k <- tcrossprod(sigma_k_raw) / 5 + diag(0.3, 5)
  sigma_r <- toeplitz(0.42^(0:12))
  components <- sampling_oracle_components(patterns, sigma_k, sigma_r)

  for (partitions in c(2L, 3L, 4L, 7L, 12L)) {
    equation <- sampling_oracle_eq13(
      components$delta, components$xi, sigma_r, partitions
    )
    enumerated <- sampling_oracle_endpoint_enumeration(
      components$delta, components$xi, sigma_r, partitions
    )

    expect_equal(
      enumerated$signal_coefficient, 4 / partitions,
      tolerance = 2e-15
    )
    expect_equal(
      enumerated$noise_coefficient,
      2 / (partitions * (partitions - 1L)),
      tolerance = 2e-15
    )
    expect_equal(
      enumerated$covariance, equation$covariance,
      tolerance = 3e-14
    )
  }
})

test_that("Eq. 13 is one complete distance-by-distance covariance", {
  set.seed(81313)
  patterns <- matrix(rnorm(4 * 9), 4, 9)
  sigma_k <- matrix(c(
    1.3, 0.2, 0.1, 0.0,
    0.2, 1.1, 0.3, 0.1,
    0.1, 0.3, 1.4, 0.2,
    0.0, 0.1, 0.2, 1.2
  ), 4, 4, byrow = TRUE)
  sigma_r <- toeplitz(0.3^(0:8))
  components <- sampling_oracle_components(patterns, sigma_k, sigma_r)
  result <- sampling_oracle_eq13(
    components$delta, components$xi, sigma_r, partitions = 6L
  )

  expect_identical(dim(result$covariance), c(6L, 6L))
  expect_equal(result$covariance, t(result$covariance), tolerance = 2e-15)
  expect_gt(max(abs(result$covariance[row(result$covariance) !=
    col(result$covariance)])), 0)
  expect_equal(result$covariance, result$signal + result$noise,
    tolerance = 0)
  expect_gt(sum(abs(result$signal)), 0)
  expect_gt(sum(abs(result$noise)), 0)
  expect_gte(min(eigen(result$covariance, symmetric = TRUE,
    only.values = TRUE)$values), -1e-12)
})

test_that("signal and noise terms have their distinct partition scaling", {
  delta <- matrix(c(2.0, 0.7, 0.7, 1.5), 2, 2)
  xi <- matrix(c(1.2, 0.4, 0.4, 0.9), 2, 2)
  sigma_r <- diag(11)

  small <- sampling_oracle_eq13(delta, xi, sigma_r, partitions = 3L)
  large <- sampling_oracle_eq13(delta, xi, sigma_r, partitions = 12L)

  expect_equal(small$signal / large$signal, matrix(4, 2, 2),
    tolerance = 2e-15)
  expect_equal(
    small$noise / large$noise,
    matrix((12 * 11) / (3 * 2), 2, 2), tolerance = 2e-15
  )
})

test_that("the two naive signal-dominated SE factors are not conflated", {
  for (partitions in c(3L, 4L, 6L, 8L, 12L)) {
    factors <- sampling_oracle_naive_signal_factors(partitions)
    expect_equal(
      unname(factors[["sample_sd"]]), sqrt(partitions + 1),
      tolerance = 2e-15
    )
    expect_equal(
      unname(factors[["marginal_edge"]]), sqrt(partitions - 1),
      tolerance = 2e-15
    )
  }
})
