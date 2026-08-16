local_fixture <- function() {
  array(
    c(
      2, 4, 10, 14, 18, 22,
      6, 8, 12, 16, 20, 24
    ),
    dim = c(2, 2, 3),
    dimnames = list(c("m1", "m2"), c("a", "b"), c("r1", "r2", "r3"))
  )
}

test_that("cross-partition independence is never silently granted", {
  undeclared <- cross_partitions(c("run-1", "run-2", "run-3"))
  declared <- cross_partitions(
    c("run-1", "run-2", "run-3"), independence = "independent"
  )

  expect_identical(attr(undeclared, "independence"), "undeclared")
  expect_identical(attr(undeclared, "estimate"),
    "independence_undeclared")
  expect_identical(attr(declared, "independence"), "independent")
  expect_identical(attr(declared, "estimate"), "cross_generalized")
  expect_false(identical(
    crossform:::.metric_pairing_identity(undeclared),
    crossform:::.metric_pairing_identity(declared)
  ))
})

test_that("undirected marginals are invariant to every stored edge orientation", {
  local <- local_fixture()
  p1 <- pairing(
    left = c("r1", "r1", "r2"),
    right = c("r2", "r3", "r3"),
    weight = c(1, 2, 3),
    directed = FALSE
  )
  p2 <- pairing(
    left = p1$right,
    right = p1$left,
    weight = p1$weight,
    directed = FALSE
  )

  a <- crossform:::pairing_marginals(local, p1, mass = c(2, 4))
  b <- crossform:::pairing_marginals(local, p2, mass = c(2, 4))

  expect_named(a, "endpoint")
  expect_equal(a$endpoint, b$endpoint, tolerance = 1e-14)
  expect_identical(attr(a, "semantics"), "undirected_endpoint")
})

test_that("undirected marginals survive arbitrary endpoint flips and edge order", {
  set.seed(20260812)
  local <- array(
    rnorm(4 * 3 * 5),
    dim = c(4, 3, 5),
    dimnames = list(paste0("m", 1:4), paste0("e", 1:3), paste0("r", 1:5))
  )
  base <- cross_partitions(dimnames(local)[[3]])
  oracle <- crossform:::pairing_marginals(local, base, mass = c(1, 2, 3, 4))$endpoint

  for (iteration in seq_len(25)) {
    flip <- sample(c(FALSE, TRUE), nrow(base), replace = TRUE)
    order <- sample(seq_len(nrow(base)))
    left <- ifelse(flip, base$right, base$left)
    right <- ifelse(flip, base$left, base$right)
    transformed <- pairing(
      left[order], right[order], base$weight[order], directed = FALSE
    )

    expect_equal(
      crossform:::pairing_marginals(local, transformed, mass = c(1, 2, 3, 4))$endpoint,
      oracle,
      tolerance = 1e-13
    )
  }
})

test_that("undirected endpoint marginal agrees with an explicit oracle", {
  local <- local_fixture()
  over <- pairing("r1", "r3", directed = FALSE)

  got <- crossform:::pairing_marginals(local, over, mass = c(2, 4))$endpoint
  r1 <- local[, , "r1"] / c(2, 4)
  r3 <- local[, , "r3"] / c(2, 4)
  oracle <- 0.5 * (r1 + r3)

  expect_equal(got, oracle, tolerance = 1e-14)
})

test_that("directed pairings preserve left and right roles", {
  local <- local_fixture()
  over <- pairing("r1", "r3", directed = TRUE)

  got <- crossform:::pairing_marginals(local, over, mass = 2)

  expect_named(got, c("left", "right"))
  expect_equal(got$left, local[, , "r1"] / 2, tolerance = 1e-14)
  expect_equal(got$right, local[, , "r3"] / 2, tolerance = 1e-14)
  expect_identical(attr(got, "semantics"), "directed_roles")
})

test_that("undirected reverse duplicates are rejected", {
  expect_error(
    pairing(c("r1", "r2"), c("r2", "r1"), directed = FALSE),
    "duplicate"
  )
})

test_that("self-products require explicit biased non-independent semantics", {
  expect_error(pairing("r1", "r1"), "allow_biased")
  expect_error(
    pairing("r1", "r1", self_pairs = "allow_biased"),
    "not_independent"
  )
  self <- pairing(
    "r1", "r1",
    self_pairs = "allow_biased",
    independence = "not_independent"
  )
  expect_identical(attr(self, "estimate"), "self_product_biased")
  expect_error(
    pairing(
      c("r1", "r1"), c("r1", "r2"),
      self_pairs = "allow_biased",
      independence = "not_independent"
    ),
    "cannot mix"
  )
})

test_that("weight normalization is overflow safe", {
  over <- pairing(
    c("r1", "r2"), c("r2", "r3"),
    weight = rep(.Machine$double.xmax, 2),
    directed = TRUE
  )
  expect_equal(over$weight, c(0.5, 0.5), tolerance = 0)
  expect_true(all(is.finite(over$weight)))
  expect_equal(sum(over$weight), 1, tolerance = 1e-15)
})

test_that("cross_partitions accepts a relation directly", {
  rel <- relation(list(
    run1 = matrix(1, 2, 3),
    run2 = matrix(2, 2, 3),
    run3 = matrix(3, 2, 3)
  ), effects = c("a", "b"))
  expect_identical(cross_partitions(rel), cross_partitions(rel$partitions))
})

test_that("pairing use sites reject forged and mutated edge tables", {
  local <- local_fixture()
  valid <- pairing("r1", "r2")

  negative <- valid
  negative$weight <- -1
  expect_error(crossform:::pairing_marginals(local, negative), "nonnegative")

  unnormalized <- valid
  unnormalized$weight <- 0.5
  expect_error(crossform:::pairing_marginals(local, unnormalized), "unit mass")

  forged <- structure(
    data.frame(left = "r1", right = "r1", weight = 1),
    directed = FALSE,
    self_pairs = "forbid",
    independence = "independent",
    estimate = "cross_generalized",
    class = c("effect_pairing", "data.frame")
  )
  expect_error(crossform:::pairing_marginals(local, forged), "explicit biased")

  duplicate <- pairing(c("r1", "r2"), c("r2", "r3"))
  duplicate$left[[2]] <- "r1"
  duplicate$right[[2]] <- "r2"
  expect_error(crossform:::pairing_marginals(local, duplicate), "duplicate")
})
