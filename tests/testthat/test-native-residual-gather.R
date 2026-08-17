# Prototype CSC gather from the canonical upper pair graph.
# Isolated timing on 800 supports of mean size 61 was 1.69x vs the R
# gather, below the 2x admission bar, so production still uses the R path.

residual_gather_fixture <- function() {
  domain <- volume_domain(array(TRUE, c(4, 3, 2)), id = "residual-gather")
  index <- crossform:::.support_index_materialize_pair_pattern(
    crossform:::.euclidean_support_index(domain, 1.25)
  )
  pattern <- index$pair_pattern
  set.seed(19)
  covariance <- rnorm(length(pattern@i))
  list(index = index, covariance = covariance, domain = domain)
}

test_that("native gather matches the retained R gather", {
  fixture <- residual_gather_fixture()
  supports <- list(
    crossform:::.support_index_support(fixture$index, 1L)[[1L]],
    crossform:::.support_index_support(fixture$index, 7L)[[1L]],
    seq_len(min(3L, fixture$domain$n_features))
  )
  pattern <- fixture$index$pair_pattern
  for (support in supports) {
    native <- crossform:::.local_residual_gather_cpp(
      as.integer(pattern@p), as.integer(pattern@i),
      as.numeric(fixture$covariance), as.integer(support)
    )
    oracle <- crossform:::.local_residual_covariance_oracle(
      fixture$index, support, fixture$covariance
    )
    expect_equal(native, oracle, tolerance = 1e-12)
    expect_equal(native, t(native), tolerance = 0)
  }
})

test_that("disconnected and singleton supports still gather exactly", {
  fixture <- residual_gather_fixture()
  singleton <- 5L
  pattern <- fixture$index$pair_pattern
  native <- crossform:::.local_residual_gather_cpp(
    as.integer(pattern@p), as.integer(pattern@i),
    as.numeric(fixture$covariance), as.integer(singleton)
  )
  oracle <- crossform:::.local_residual_covariance_oracle(
    fixture$index, singleton, fixture$covariance
  )
  expect_equal(native, oracle, tolerance = 1e-12)
  expect_identical(dim(native), c(1L, 1L))
})
