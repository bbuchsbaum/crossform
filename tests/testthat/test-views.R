view_geometry_fixture <- function() {
  effects <- c("a", "b", "c", "d")
  total_matrices <- list(
    matrix(c(
      4, 1, 0, 0,
      1, 3, 0.5, 0,
      0, 0.5, 2, 0.25,
      0, 0, 0.25, 1
    ), 4, 4, byrow = TRUE),
    diag(c(2, 1, -0.5, -2))
  )
  coherent_matrices <- list(
    tcrossprod(c(1, 0.5, 0, 0)),
    diag(c(0.5, 0.25, -0.1, -0.25))
  )
  pack <- function(values) do.call(rbind,
    lapply(values, effectagram:::.svec_symmetric))
  endpoint <- matrix(c(
    1, 2, 3, 4,
    -1, 0, 1, 2
  ), 2, 4, byrow = TRUE, dimnames = list(NULL, effects))
  marginals <- structure(list(endpoint = endpoint),
    semantics = "undirected_endpoint", class = c("effect_marginals", "list"))
  geometry <- effect_geometry(pack(total_matrices), pack(coherent_matrices),
    marginals, effects, receipt_fixture(), index = c("left", "right"))
  list(geometry = geometry, total = total_matrices,
    coherent = coherent_matrices, endpoint = endpoint)
}

test_that("contrast returns one exact decomposition and signed marginal", {
  fixture <- view_geometry_fixture()
  weights <- c(a = 1, b = -1, c = 0, d = 0)
  got <- contrast(fixture$geometry, weights)
  expected_total <- vapply(fixture$total,
    function(value) drop(weights %*% value %*% weights), numeric(1))
  expected_coherent <- vapply(fixture$coherent,
    function(value) drop(weights %*% value %*% weights), numeric(1))

  expect_equal(got$signed, drop(fixture$endpoint %*% weights), tolerance = 0)
  expect_equal(got$total, expected_total, tolerance = 1e-12)
  expect_equal(got$coherent, expected_coherent, tolerance = 1e-12)
  expect_equal(got$configuration, expected_total - expected_coherent,
    tolerance = 1e-12)
  expect_equal(got$total, got$coherent + got$configuration,
    tolerance = 0)
  expect_true(all(is.na(got$coherence_fraction[!got$coherence_fraction_valid])))
})

test_that("contrast names are semantic rather than positional when supplied", {
  fixture <- view_geometry_fixture()
  ordered <- contrast(fixture$geometry, c(a = 1, b = -1, c = 0, d = 0))
  permuted <- contrast(fixture$geometry, c(d = 0, b = -1, a = 1, c = 0))
  expect_equal(permuted$total, ordered$total, tolerance = 0)
  expect_error(contrast(fixture$geometry, c(a = 1, b = -1, c = 0, x = 0)),
    "identify every effect")
})

test_that("RDM view equals explicit squared geometry distances", {
  fixture <- view_geometry_fixture()
  got <- rdm(fixture$geometry)
  pairs <- utils::combn(seq_len(4), 2)
  expected <- t(vapply(fixture$total, function(value) {
    vapply(seq_len(ncol(pairs)), function(pair) {
      delta <- numeric(4)
      delta[pairs[, pair]] <- c(1, -1)
      drop(delta %*% value %*% delta)
    }, numeric(1))
  }, numeric(ncol(pairs))))

  expect_equal(got$values, expected, tolerance = 1e-12,
    ignore_attr = TRUE)
  expect_identical(got$pairs$left, c("a", "a", "a", "b", "b", "c"))
  expect_identical(got$pairs$right, c("b", "c", "d", "c", "d", "d"))
})

test_that("compiled multiple-regression RSA equals explicit RDM regression", {
  fixture <- view_geometry_fixture()
  model <- matrix(c(
    0, 1, 2, 3,
    1, 0, 1, 2,
    2, 1, 0, 1,
    3, 2, 1, 0
  ), 4, 4, byrow = TRUE)
  nuisance <- matrix(c(
    0, 0, 1, 1,
    0, 0, 1, 1,
    1, 1, 0, 0,
    1, 1, 0, 0
  ), 4, 4, byrow = TRUE)
  got <- rsa(fixture$geometry, models = list(ordinal = model),
    nuisance = list(group = nuisance), intercept = TRUE)
  distances <- rdm(fixture$geometry)$values
  design <- cbind(`(Intercept)` = 1,
    ordinal = effectagram:::.rdm_vector(model),
    group = effectagram:::.rdm_vector(nuisance))
  expected <- t(apply(distances, 1, function(value) {
    qr.coef(qr(design), value)
  }))

  expect_equal(got$coefficients, expected, tolerance = 1e-12,
    ignore_attr = TRUE)
  expect_identical(got$terms$role, c("intercept", "model", "nuisance"))
})

test_that("thin QR coefficient maps equal the dense-identity oracle", {
  set.seed(2026081502L)
  design <- cbind(
    intercept = 1,
    weak = stats::rnorm(300L, sd = 0.01),
    strong = stats::rnorm(300L, sd = 10)
  )
  decomposition <- qr(design, LAPACK = FALSE)
  got <- effectagram:::.thin_qr_coefficient_map(decomposition)
  oracle <- qr.coef(decomposition, diag(nrow(design)))

  expect_equal(got, oracle, tolerance = 1e-12)
  expect_identical(dim(got), c(ncol(design), nrow(design)))
})

test_that("named RSA axes align by experimental identity", {
  fixture <- view_geometry_fixture()
  model <- as.matrix(dist(seq_len(4)))
  dimnames(model) <- list(fixture$geometry$effects, fixture$geometry$effects)
  permutation <- c("d", "b", "a", "c")
  ordered <- rsa(fixture$geometry, list(model = model), intercept = TRUE)
  permuted <- rsa(fixture$geometry,
    list(model = model[permutation, permutation]), intercept = TRUE)
  expect_equal(permuted$coefficients, ordered$coefficients, tolerance = 0)
})

test_that("RSA rejects rank-deficient and malformed RDM designs", {
  fixture <- view_geometry_fixture()
  model <- matrix(1, 4, 4)
  diag(model) <- 0
  expect_error(rsa(fixture$geometry,
    models = list(first = model, duplicate = model)), "rank deficient")
  asymmetric <- model
  asymmetric[1, 2] <- 2
  expect_error(rsa(fixture$geometry, list(bad = asymmetric)), "symmetric")
})

test_that("geometry spectrum preserves negative cross-generalized roots", {
  fixture <- view_geometry_fixture()
  got <- geometry_spectrum(fixture$geometry, row_block = 1)
  expected <- t(vapply(fixture$total, function(value) {
    eigen(value, symmetric = TRUE, only.values = TRUE)$values
  }, numeric(4)))

  expect_equal(got$values, expected, tolerance = 1e-12,
    ignore_attr = TRUE)
  expect_true(any(got$values < 0))
  expect_true(got$indefinite_estimates_preserved)
})
