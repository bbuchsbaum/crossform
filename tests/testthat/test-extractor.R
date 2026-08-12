test_that("explicit extractors are named finite linear maps", {
  extractor <- effect_extractor(matrix(c(1, -1, 0, 1), 2), c("a", "b"))

  expect_s3_class(extractor, "effect_extractor")
  expect_identical(extractor$n_observations, 2L)
  expect_identical(rownames(extractor$map), c("a", "b"))
  expect_error(effect_extractor(matrix(Inf, 1, 1)), "finite")
  expect_error(effect_extractor(diag(2), c("same", "same")), "unique")
})

test_that("full-rank lm extraction agrees with explicit OLS", {
  set.seed(51)
  design <- cbind(1, x = seq(-1, 1, length.out = 12))
  target <- rbind(intercept = c(1, 0), slope = c(0, 1))
  response <- matrix(rnorm(12 * 5), 12, 5)
  extractor <- lm_extractor(design, target)

  expect_identical(extractor$diagnostics$solver, "pivoted_qr")
  expect_equal(extractor$map %*% response,
    target %*% qr.coef(qr(design), response), tolerance = 1e-12)
})

test_that("whitened lm extraction agrees with transformed least squares", {
  set.seed(52)
  design <- cbind(1, x = seq_len(10))
  target <- matrix(c(0, 1), 1, dimnames = list("slope", NULL))
  whiten <- diag(seq(0.5, 1.5, length.out = 10))
  response <- matrix(rnorm(30), 10, 3)
  extractor <- lm_extractor(design, target, whiten = whiten)

  oracle <- target %*% qr.coef(qr(whiten %*% design), whiten %*% response)
  expect_equal(extractor$map %*% response, oracle, tolerance = 1e-12)
})

test_that("rank-deficient extraction admits only estimable targets", {
  x <- seq_len(8)
  design <- cbind(x, duplicate = x)
  estimable <- matrix(c(1, 1), 1, dimnames = list("sum", NULL))
  extractor <- lm_extractor(design, estimable)

  expect_true(extractor$diagnostics$rank_deficient)
  expect_identical(extractor$diagnostics$solver, "svd_estimable_fallback")
  expect_error(
    lm_extractor(design, matrix(c(1, 0), 1,
      dimnames = list("first_only", NULL))),
    "not estimable"
  )
})

test_that("mutated extractor contracts fail closed", {
  extractor <- effect_extractor(diag(2), c("a", "b"))
  extractor$n_observations <- 3L
  expect_error(effectagram:::.validate_effect_extractor(extractor),
    "inconsistent")
})
