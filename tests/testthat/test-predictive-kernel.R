# Independent input/kernel checks for predictive-geometry-v1. No package
# Gram helper constructs the expected kernel or distances.
pk_features <- function() {
  rbind(a = c(2, 0), b = c(0, 1), c = c(-1, 0), d = c(0, -2))
}
pk_gram <- function(F) {
  H <- diag(nrow(F)) - matrix(1/nrow(F), nrow(F), nrow(F))
  K <- H %*% F %*% t(F) %*% H
  dimnames(K) <- list(rownames(F), rownames(F))
  K
}
pk_rdm <- function(F) {
  D <- matrix(0, nrow(F), nrow(F), dimnames = list(rownames(F), rownames(F)))
  for (i in seq_len(nrow(F))) for (j in seq_len(nrow(F))) {
    D[i, j] <- sum((F[i, ] - F[j, ])^2)
  }
  D
}
pk_value <- function(basis) tcrossprod(basis$Q %*% basis$R)

test_that("[T01] explicit kernels, features and squared distances retain the same geometry", {
  F <- pk_features()
  expected <- pk_gram(F)
  values <- list(
    model_basis(features = list(m = F)),
    model_basis(models = list(m = pk_rdm(F)), distance = "squared_euclidean"),
    model_basis(kernels = list(m = expected))
  )
  for (value in values) {
    expect_equal(pk_value(value), expected, tolerance = 1e-12)
    expect_identical(value$models$rank_effective, 2L)
    expect_lt(max(abs(colSums(value$Q))), 1e-12)
  }
  expect_identical(values[[3L]]$models$kind, "kernel")
  expect_length(unique(vapply(values, `[[`, "", "signature")), 3L)
  expect_error(model_basis(models = pk_rdm(F)), "distance")
})

test_that("[T02] kernel axes align independently and ambiguous labels refuse", {
  K <- pk_gram(pk_features())
  nm <- rownames(K)
  scrambled <- K[c(3, 1, 4, 2), c(2, 4, 1, 3)]
  expect_equal(pk_value(model_basis(kernels = scrambled, conditions = nm)), K,
    tolerance = 1e-12)
  expect_equal(pk_value(model_basis(kernels = unname(K), conditions = nm)), K,
    tolerance = 1e-12)
  expect_error(model_basis(kernels = unname(K)), "conditions",
    class = "effect_input_error")
  one_axis <- K; colnames(one_axis) <- NULL
  expect_error(model_basis(kernels = one_axis), "only its rows",
    class = "effect_input_error")
  bad <- K; rownames(bad)[2] <- rownames(bad)[1]
  expect_error(model_basis(kernels = bad, conditions = nm), "axis names",
    class = "effect_input_error")
  expect_error(model_basis(kernels = K[-1, -1], conditions = nm), "missing",
    class = "effect_input_error")
  expect_error(model_basis(kernels = K[, -1]), "square",
    class = "effect_input_error")
})

test_that("[T03] kernel PSD admission happens before centering", {
  n <- 4L; H <- diag(n) - 1/n
  bad <- H - 2 * matrix(1/n, n, n)
  nm <- letters[1:n]; dimnames(bad) <- list(nm, nm)
  refused <- catch_refusal(model_basis(kernels = bad))
  expect_s3_class(refused, "effect_capability_refusal")
  expect_identical(refused$reasons[[1]], "model_gram_indefinite")
  expect_match(conditionMessage(refused), "before centering")
  repaired <- model_basis(kernels = bad, negative_share = 1)
  expect_equal(unname(pk_value(repaired)), H, tolerance = 1e-12)
  expect_equal(repaired$models$input_negative_mass, 2, tolerance = 1e-12)
  expect_lt(repaired$models$dropped_negative_mass, 1e-10)
  expect_error(model_basis(kernels = matrix(1, n, n), conditions = nm),
    "no geometry", class = "effect_input_error")
  expect_error(model_basis(kernels = matrix(0, n, n), conditions = nm),
    "no geometry", class = "effect_input_error")
})

test_that("[T03] invalid kernels fail at the input boundary and numerical roots are recorded", {
  K <- pk_gram(pk_features())
  for (bad_number in c(NA_real_, NaN, Inf, -Inf)) {
    bad <- K; bad[1, 1] <- bad_number
    expect_error(model_basis(kernels = bad), "non-finite", class = "effect_input_error")
  }
  bad <- K; bad[1, 2] <- bad[1, 2] + 0.01
  expect_error(model_basis(kernels = bad), "symmetric", class = "effect_input_error")
  V <- qr.Q(qr(cbind(1, diag(4))))[, -1]
  K <- V %*% diag(c(3, 1, -1e-12)) %*% t(V)
  value <- model_basis(kernels = K, conditions = letters[1:4])
  expect_identical(value$models$input_negative_mass, 0)
  expect_gt(value$models$input_numerical_negative_mass, 0)
  expect_lt(value$models$input_numerical_negative_mass, 1e-9)
  expect_identical(value$models$rank_effective, 2L)
  expect_error(model_basis(features = pk_features() * 1e200), "overflows",
    class = "effect_input_error")
})

test_that("[T04] model coordinates require shared declared effect units and scales", {
  K <- pk_gram(pk_features())
  expect_error(model_basis(kernels = K,
    conditions = effect_space(rownames(K), units = c("m", "m", "m", "s"))),
    "different units", class = "effect_input_error")
  expect_error(model_basis(kernels = K,
    conditions = effect_space(rownames(K), scale = c(1, 1, 1, 2))),
    "different scales", class = "effect_input_error")
  value <- model_basis(kernels = K,
    conditions = effect_space(rownames(K), units = "percent", scale = 2))
  expect_true(all(value$effect_space$units == "percent"))
  expect_true(all(value$effect_space$scale == 2))
})

# An explicit orthonormal frame keeps analytic eigenvalues independent of
# the constructor's eigenvectors, signs and centered-frame implementation.
pk_frame <- function() {
  cbind(c(1, -1, 0, 0)/sqrt(2), c(1, 1, -2, 0)/sqrt(6),
    c(1, 1, 1, -3)/sqrt(12))
}
pk_kernel <- function(roots) {
  V <- pk_frame()
  K <- V %*% diag(roots, 3) %*% t(V)
  dimnames(K) <- list(letters[1:4], letters[1:4])
  K
}

test_that("[T05] trace normalization follows truncation and records every scale", {
  K <- pk_kernel(c(8, 2, 0))
  b <- model_basis(kernels = K, rank = 1, normalize = "trace")
  expected <- tcrossprod(pk_frame()[, 1])
  expect_equal(unname(pk_value(b)), expected, tolerance = 1e-12)
  expect_equal(b$models$raw_trace, 10, tolerance = 1e-12)
  expect_equal(b$models$retained_trace, 8, tolerance = 1e-12)
  expect_equal(b$models$normalization_scale, 1/8, tolerance = 1e-12)
  expect_equal(b$models$truncated_positive_mass, 2, tolerance = 1e-12)
  expect_equal(sum(b$R^2), 1, tolerance = 1e-12)
  expect_identical(b$normalize, "trace")
  expect_identical(crossform:::.validate_model_basis(b), b)
  expect_gt(max(abs(unname(pk_value(b)) - expected * 0.8)), 0.05)
})

test_that("[T06 T10] normalized model geometry obeys scale and baseline laws", {
  F <- pk_features()
  base <- model_basis(features = F, rank = 1, normalize = "trace")
  for (scale in c(1e-8, 1e-4, 1, 1e4, 1e8)) {
    b <- model_basis(kernels = pk_gram(F) * scale, rank = 1, normalize = "trace")
    expect_equal(pk_value(b), pk_value(base), tolerance = 1e-11)
    expect_identical(b$models$rank_effective, base$models$rank_effective)
    expect_identical(crossform:::.validate_model_basis(b), b)
  }
  shifted <- sweep(F, 2, c(300, -100), "+")
  expect_equal(pk_value(model_basis(features = shifted, normalize = "trace")),
    pk_value(model_basis(features = F, normalize = "trace")), tolerance = 1e-11)
  neural <- matrix(seq_len(20), 4, 5)
  expect_equal(crossprod(base$Q, neural),
    crossprod(base$Q, sweep(neural, 2, c(-100, 20, 4, 300, 2), "+")),
    tolerance = 1e-11)
  poor <- pk_frame() %*% diag(c(1e-3, 1e-4, 1e-5))
  rownames(poor) <- letters[1:4]
  b <- model_basis(features = poor, tolerance = 1e-12, rank = 2, normalize = "trace")
  expect_identical(b$dimension, 2L)
  expect_lt(max(abs(colSums(b$Q))), 1e-12)
})

test_that("[T07] model ranks clamp explicitly and reject invalid or overflowing requests", {
  K <- pk_kernel(c(8, 2, 0))
  b <- model_basis(kernels = K, rank = 99, normalize = "trace")
  expect_identical(b$models$rank_requested, 99L)
  expect_identical(b$models$rank_effective, 2L)
  expect_true(b$models$rank_clamped)
  for (bad in list(0, -1, 1.5, NA_real_, Inf, 2^31, "2", c(1, 2))) {
    expect_error(model_basis(kernels = K, rank = bad), class = "effect_input_error")
  }
  expect_error(model_basis(kernels = K * 0, normalize = "trace"), "no geometry")
})

test_that("[T11] positive model ties require complete retained groups", {
  K <- pk_kernel(c(2, 2, 0))
  refused <- catch_refusal(model_basis(kernels = K, rank = 1))
  expect_identical(refused$reasons[[1]], "ambiguous_rank_cut")
  b <- model_basis(kernels = K, rank = 2, normalize = "trace")
  V <- pk_frame()[, 1:2]
  O <- matrix(c(0.6, -0.8, 0.8, 0.6), 2)
  F <- sqrt(2) * V %*% O; rownames(F) <- letters[1:4]
  rotated <- model_basis(features = F, rank = 2, normalize = "trace")
  expect_equal(pk_value(b), pk_value(rotated), tolerance = 1e-12)
  expect_equal(unname(pk_value(b)), tcrossprod(V)/2, tolerance = 1e-12)
  # A resolved gap is not mistaken for a tie.
  expect_identical(model_basis(kernels = pk_kernel(c(2, 2 - 1e-6, 0)), rank = 1)$dimension, 1L)
})

test_that("[T03 T23] numerical roots and normalization metadata stay distinct and sealed", {
  b <- model_basis(kernels = pk_kernel(c(3, 1e-12, -1e-12)), normalize = "trace")
  expect_identical(b$models$rank_effective, 1L)
  expect_lt(abs(b$models$numerical_positive_mass - 1e-12), 1e-14)
  expect_lt(abs(b$models$numerical_negative_mass - 1e-12), 1e-14)
  expect_identical(b$models$dropped_negative_mass, 0)
  other <- model_basis(kernels = pk_kernel(c(3, 1e-12, -1e-12)), normalize = "none")
  expect_false(identical(b$signature, other$signature))
  for (field in c("normalize", "centering", "spectra")) {
    forged <- b; forged[[field]] <- other[[field]]
    if (identical(forged, b)) next
    expect_error(crossform:::.validate_model_basis(forged), class = "effect_contract_error")
  }
  forged <- b; forged$models$retained_trace <- 2
  expect_error(crossform:::.validate_model_basis(forged), class = "effect_contract_error")
  # Even a self-consistent digest cannot turn a false mass record into a fact.
  forged$signature <- crossform:::.sha256_signature(crossform:::.model_basis_semantic(forged),
    "model-basis-sha256:")
  expect_error(crossform:::.validate_model_basis(forged), "normalization or spectral mass")
  units <- model_basis(kernels = pk_kernel(c(3, 1, 0)),
    conditions = effect_space(letters[1:4], units = "percent", scale = 2))
  arbitrary <- model_basis(kernels = pk_kernel(c(3, 1, 0)))
  expect_false(identical(units$signature, arbitrary$signature))
  expect_identical(units$centering$method, "condition_mean")
})

pk_pool_value <- function(b, weights = NULL) {
  p <- crossform:::.model_basis_pool(b, weights)
  b$Q %*% p$effective_kernel %*% t(b$Q)
}

test_that("[T08] duplicates with split weights preserve the pooled kernel", {
  F <- pk_features()
  one <- model_basis(features = list(m = F), normalize = "trace")
  duplicate <- model_basis(features = list(a = F, b = F), normalize = "trace")
  expect_identical(dim(duplicate$R), c(2L, 4L))
  expect_identical(duplicate$span_overlap_rank, 2L)
  expect_equal(pk_pool_value(duplicate, c(a = .3, b = .7)),
    pk_pool_value(one), tolerance = 1e-12)
  expect_equal(pk_pool_value(duplicate, c(a = 1, b = 0)),
    pk_pool_value(one), tolerance = 1e-12)
  expect_identical(crossform:::.validate_model_basis(duplicate), duplicate)
  V <- pk_frame()
  A <- V[, 1:2]; B <- V[, 2:3]
  rownames(A) <- rownames(B) <- letters[1:4]
  overlap <- model_basis(features = list(a = A, b = B), normalize = "trace")
  expect_identical(overlap$span_overlap_rank, 1L)
  expected <- .25 * tcrossprod(V[, 1]) + .5 * tcrossprod(V[, 2]) + .25 * tcrossprod(V[, 3])
  expect_equal(unname(pk_pool_value(overlap, c(a = .5, b = .5))), expected, tolerance = 1e-12)
})

test_that("[T09] simplex vertices reduce positive support and malformed weights refuse", {
  V <- pk_frame()
  A <- V[, 1, drop = FALSE]; B <- V[, 2, drop = FALSE]
  rownames(A) <- rownames(B) <- letters[1:4]
  b <- model_basis(features = list(a = A, b = B), normalize = "trace")
  for (weights in list(c(a = 1, b = 0), c(a = 0, b = 1))) {
    p <- crossform:::.model_basis_pool(b, weights)
    expect_identical(p$dimension, 1L)
    expected <- tcrossprod(if (weights[[1]] == 1) A else B)
    expect_equal(pk_pool_value(b, weights), expected, tolerance = 1e-12)
    expect_equal(p$effective_kernel %*% p$vectors, p$vectors %*% diag(p$values, 1),
      tolerance = 1e-12)
  }
  invalid <- list(NULL, c(.5, .5), c(a = 0, b = 0), c(a = 1, b = 1),
    c(a = -1, b = 2), c(a = NA_real_, b = 1), c(a = Inf, b = 0),
    c(a = .5, other = .5), c(a = .5, a = .5), c(a = 1),
    c(a = .5, b = .5 + 1e-10), matrix(c(.5, .5), 1))
  for (weights in invalid) {
    expect_error(crossform:::.model_basis_pool(b, weights), class = "effect_input_error")
  }
  canonical <- crossform:::.model_basis_pool(b, c(a = .2, b = .8))
  expect_identical(canonical, crossform:::.model_basis_pool(b, c(b = .8, a = .2)))
  # Admission tolerance is not an instruction to renormalize user weights.
  almost <- c(a = .5, b = .5 + 1e-13)
  expect_identical(crossform:::.model_basis_pool(b, almost)$weights, almost)
})

test_that("[T06 T18 T23] pooling obeys rescaling and reports union and support conditioning", {
  V <- pk_frame()
  A <- V[, 1:2] %*% diag(c(2, .3)); B <- V[, 2:3] %*% diag(c(.8, .1))
  rownames(A) <- rownames(B) <- letters[1:4]
  base <- model_basis(features = list(a = A, b = B), normalize = "trace")
  expected <- pk_pool_value(base, c(a = .4, b = .6))
  for (scale in c(1e-8, 1e-4, 1e4, 1e8)) {
    b <- model_basis(features = list(a = A * sqrt(scale), b = B / sqrt(scale)),
      normalize = "trace")
    expect_equal(pk_pool_value(b, c(a = .4, b = .6)), expected, tolerance = 1e-11)
  }
  near <- cbind(V[, 1] + 1e-4 * V[, 2])
  rownames(near) <- letters[1:4]
  a <- V[, 1, drop = FALSE]; rownames(a) <- letters[1:4]
  b <- model_basis(features = list(a = a, near = near), normalize = "trace")
  p <- crossform:::.model_basis_pool(b, c(a = .5, near = .5))
  expect_identical(b$dimension, 2L)
  expect_identical(p$dimension, 2L)
  expect_gt(b$basis_condition, 1e4)
  expect_gt(p$condition, 1e8)
  expect_equal(p$condition, b$basis_condition^2, tolerance = 1e-7)
  # Eigenvalue support and factor singular-value support use their declared
  # thresholds in their own spaces, not a hidden square-root substitution.
  p <- crossform:::.model_basis_pool(b, c(a = 1, near = 0))
  expect_identical(p$dimension, 1L)
  expect_identical(p$union_dimension, 2L)
  expect_lt(max(abs(p$effective_kernel - t(p$effective_kernel))), 1e-12)
})
