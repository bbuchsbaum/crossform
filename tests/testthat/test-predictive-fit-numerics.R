# Exact fixtures mirror the independent base-R court under design/oracles,
# but these tests also ship in the source tarball (design/ does not).
# The 2x2 expected solution uses the quadratic formula, never package eigen
# helpers. Larger spectral and numerical-optimization checks are secondary.
pf_norm <- function(A) sum(A * A)
pf_positive2 <- function(B) {
  a <- B[1, 1]; b <- B[1, 2]; d <- B[2, 2]
  gap <- sqrt((a - d)^2 + 4*b^2)
  high <- (a + d + gap)/2; low <- (a + d - gap)/2
  high * (B - low * diag(2))/gap
}
pf_case <- function(S, J, penalty, rank, tolerance = 1e-10) {
  r <- nrow(J); n <- r + 1L
  V <- qr.Q(qr(cbind(1, diag(n))))[, -1, drop = FALSE]
  K <- V %*% J %*% t(V); union <- tcrossprod(V)
  nm <- paste0("c", seq_len(n))
  dimnames(K) <- dimnames(union) <- list(nm, nm)
  # The zero-weight union model makes the null-space trap testable even when
  # the active kernel is singular. It is not the predictor's support.
  basis <- model_basis(kernels = list(active = K, union = union), tolerance = tolerance)
  pool <- crossform:::.model_basis_pool(basis, c(active = 1, union = 0))
  map <- crossprod(basis$Q, V)
  reduced <- map %*% S %*% t(map)
  path <- crossform:::.geometry_fit_path(reduced, pool, penalty)
  fit <- crossform:::.geometry_fit_rank(path, rank)
  list(fit = fit, path = path, pool = pool, basis = basis, map = map,
    form = t(map) %*% fit$form %*% map)
}
pf_dense <- function(S, J, rank, penalty, tolerance = 1e-10) {
  k <- eigen(J, symmetric = TRUE)
  keep <- k$values > tolerance * max(abs(k$values))
  Q <- k$vectors[, keep, drop = FALSE]
  shifted <- crossprod(Q, S %*% Q) - diag(penalty/k$values[keep], sum(keep))
  e <- eigen(shifted, symmetric = TRUE)
  mu <- pmax(e$values, 0)
  mu[seq_along(mu) > rank] <- 0
  Z <- Q %*% e$vectors
  Z %*% diag(mu, length(mu)) %*% t(Z)
}

test_that("[T12 T16] inverse-kernel shift has the exact F01 optimum and complete rank path", {
  S <- diag(c(6, 3, -1)); J <- diag(c(.5, .25, .25))
  c <- pf_case(S, J, .5, 2)
  expect_equal(c$path$signed_spectrum, c(6, 3, -1), tolerance = 1e-12)
  expect_equal(c$path$shifted_spectrum, c(5, 1, -3), tolerance = 1e-12)
  expect_equal(c$form, diag(c(5, 1, 0)), tolerance = 1e-12)
  expect_equal(c$fit$objective, 10, tolerance = 1e-12)
  expect_equal(c$fit$source_negative_mass, 1, tolerance = 1e-12)
  expect_equal(c$fit$shifted_negative_mass, 3, tolerance = 1e-12)
  expect_identical(c$fit$rank_effective, 2L)
  expected <- list(matrix(0, 3, 3), diag(c(5, 0, 0)), diag(c(5, 1, 0)), diag(c(5, 1, 0)))
  frozen <- c$path
  for (rank in 0:3) {
    fit <- crossform:::.geometry_fit_rank(c$path, rank)
    expect_equal(t(c$map) %*% fit$form %*% c$map, expected[[rank + 1L]], tolerance = 1e-12)
    expect_identical(c$path, frozen)
  }
  oversized <- crossform:::.geometry_fit_rank(c$path, 99)
  expect_identical(oversized$rank_requested, 99L)
  expect_identical(oversized$rank_budget, 3L)
  expect_true(oversized$rank_clamped)
  expect_equal(oversized$form, c$fit$form)
})

test_that("[T09 T13 T15] noncommuting penalties and singular pooled support use the exact objective", {
  S <- matrix(c(1, 2, 2, 1), 2); J <- diag(c(.8, .2))
  expected <- pf_positive2(matrix(c(.75, 2, 2, 0), 2))
  c <- pf_case(S, J, .2, 1)
  expect_equal(c$form, expected, tolerance = 1e-12)
  wrong <- pf_positive2(matrix(1.5, 2, 2) - diag(c(.25, 1)))
  expect_gt(max(abs(c$form - wrong)), .01)
  O <- matrix(c(.6, -.8, .8, .6), 2)
  rotated <- pf_case(O %*% S %*% t(O), O %*% J %*% t(O), .2, 1)
  expect_equal(t(O) %*% rotated$form %*% O, expected, tolerance = 1e-12)
  singular <- pf_case(diag(c(2, 7)), diag(c(1, 0)), 1, 2)
  expect_equal(singular$form, diag(c(1, 0)), tolerance = 1e-12)
  expect_identical(singular$pool$dimension, 1L)
  expect_identical(singular$pool$union_dimension, 2L)
  expect_gt(max(abs(singular$form - diag(c(1, 7)))), 6)
})

test_that("[T14 T16] signed compression and zero budgets never invent a positive mode", {
  G <- matrix(c(1, 2, 2, 1), 2)
  Q <- matrix(c(1, 0), 2)
  signed <- crossprod(Q, G %*% Q)
  c <- pf_case(signed, matrix(1), 0, 1)
  expect_equal(c$form, matrix(1), tolerance = 1e-12)
  expect_gt(abs(c$form[[1]] - 1.5), .4)
  for (S in list(matrix(0, 2, 2), -diag(c(1, 2)))) {
    c <- pf_case(S, diag(c(.8, .2)), .1, 2)
    expect_identical(c$fit$rank_effective, 0L)
    expect_identical(c$fit$prediction_norm_sq, 0)
    expect_equal(c$fit$form, matrix(0, 2, 2), ignore_attr = TRUE)
    expect_identical(dim(c$fit$factor), c(2L, 0L))
  }
  c <- pf_case(diag(c(2, 1)), diag(c(.8, .2)), 0, 0)
  expect_identical(c$fit$prediction_norm_sq, 0)
  expect_identical(c$fit$rank_budget, 0L)
  expect_identical(c$fit$penalty_trace, 0)
  expect_identical(dim(c$fit$modes), c(2L, 0L))
  expect_equal(crossform:::.latent_rank_psd_form(diag(2), 0)$form, matrix(0, 2, 2))
  for (rank in list(-1, 1.2, Inf, NA_real_, 2^31, c(1, 2))) {
    expect_error(crossform:::.geometry_fit_rank(c$path, rank), class = "effect_input_error")
  }
})

test_that("[T17 T18] rank and penalty paths obey optimization and scale laws", {
  S <- matrix(c(3, 1, .2, 1, 2, .5, .2, .5, -.5), 3)
  J <- diag(c(.6, .3, .1))
  c <- pf_case(S, J, .15, 2)
  objective <- vapply(0:3, function(r) crossform:::.geometry_fit_rank(c$path, r)$objective, 0.)
  expect_true(all(diff(objective) <= 1e-12))
  penalties <- c(0, .01, .1, .3, 1)
  fits <- lapply(penalties, function(p) pf_case(S, J, p, 2)$fit)
  expect_true(all(diff(vapply(fits, `[[`, 0., "objective")) >= -1e-12))
  expect_true(all(diff(vapply(fits, `[[`, 0., "penalty_trace")) <= 1e-12))
  for (scale in c(1e-4, .3, 20, 1e4)) {
    response <- pf_case(S * scale, J, .15 * scale, 2)
    expect_equal(response$form / scale, c$form, tolerance = 1e-10)
    model <- pf_case(S, J * scale, .15 * scale, 2)
    expect_equal(model$form, c$form, tolerance = 1e-10)
    test <- diag(c(2, -1, .4))
    gain <- function(A, B) 2*sum(A*B) - pf_norm(A)
    expect_equal(gain(response$form, test * scale)/scale^2, gain(c$form, test), tolerance = 1e-10)
  }
  perm <- c(3, 1, 2)
  permuted <- pf_case(S[perm, perm], J[perm, perm], .15, 2)
  expect_equal(permuted$form[order(perm), order(perm)], c$form, tolerance = 1e-11)
})

test_that("[T19] full and compressed objectives differ by an identified constant only", {
  G <- matrix(c(3, 1, 2, 1, -1, .4, 2, .4, 5), 3)
  Q <- cbind(c(1, 0, 0), c(0, 1, 0)); J <- diag(c(.8, .2))
  S <- crossprod(Q, G %*% Q); penalty <- .2
  B <- S - penalty * solve(J)
  candidates <- list(matrix(0, 2, 2), diag(c(.4, 1)), tcrossprod(matrix(c(.3, -.6), 2)))
  compressed <- full <- numeric(length(candidates))
  for (i in seq_along(candidates)) {
    A <- candidates[[i]]
    compressed[i] <- .5*pf_norm(S - A) + penalty * sum(solve(J) * A)
    full[i] <- .5*pf_norm(G - Q %*% A %*% t(Q)) + penalty * sum(solve(J) * A)
    expect_equal(compressed[i], .5*pf_norm(A - B) + .5*pf_norm(S) - .5*pf_norm(B), tolerance = 1e-12)
    expect_equal(full[i], .5*pf_norm(A - B) + .5*pf_norm(G) - .5*pf_norm(B), tolerance = 1e-12)
  }
  expect_equal(diff(full), diff(compressed), tolerance = 1e-12)
  expect_gt(full[1] - compressed[1], 1)
  fitted <- pf_case(S, J, penalty, 1)
  expect_equal(fitted$fit$objective, .5*pf_norm(S - fitted$form) +
    penalty*sum(solve(J) * fitted$form), tolerance = 1e-12)
})

test_that("[T20] redundant factors give the minimum-norm penalty and the global optimum", {
  T <- cbind(diag(c(sqrt(.8), sqrt(.2)))/sqrt(2),
    diag(c(sqrt(.8), sqrt(.2)))/sqrt(2))
  J <- tcrossprod(T); S <- matrix(c(1, 2, 2, 1), 2); penalty <- .2
  result <- pf_case(S, J, penalty, 1)
  e <- eigen(result$form, symmetric = TRUE)
  Z <- sqrt(e$values[1]) * e$vectors[, 1, drop = FALSE]
  W <- t(T) %*% solve(tcrossprod(T), Z)
  expect_equal(T %*% W %*% t(W) %*% t(T), result$form, tolerance = 1e-12)
  expect_equal(sum(W^2), sum(solve(J) * result$form), tolerance = 1e-12)
  objective <- function(w) {
    F <- tcrossprod(T %*% matrix(w, ncol = 1))
    .5*pf_norm(S - F) + penalty*sum(w^2)
  }
  starts <- list(c(.1, .2, .3, .4), c(1, -.4, .2, .8), c(-.5, 1, .1, -.3))
  optimized <- vapply(starts, function(w) optim(w, objective, method = "BFGS",
    control = list(maxit = 3000, reltol = 1e-12))$value, 0.)
  expect_true(all(optimized >= result$fit$objective - 1e-10))
  expect_lt(max(abs(optimized - result$fit$objective)), 1e-7)
})

test_that("[T21 T22] model preference survives full span and ties have an explicit policy", {
  S <- diag(c(3, 3)); J1 <- diag(c(.8, .2)); J2 <- diag(c(.2, .8))
  unregularized <- lapply(list(J1, J2), function(J) pf_case(S, J, 0, 2))
  expect_equal(unregularized[[1]]$form, unregularized[[2]]$form, tolerance = 1e-12)
  expect_true(unregularized[[1]]$basis$saturated)
  expect_equal(pf_case(S, J1, .4, 2)$form, diag(c(2.5, 1)), tolerance = 1e-12)
  expect_equal(pf_case(S, J2, .4, 2)$form, diag(c(1, 2.5)), tolerance = 1e-12)
  for (gap in c(0, 1e-12)) {
    refusal <- catch_refusal(pf_case(diag(c(2, 2 - gap)), diag(2), 0, 1))
    expect_identical(refusal$reasons[[1]], "ambiguous_rank_cut")
  }
  retained <- pf_case(diag(c(2, 2)), diag(2), 0, 2)
  expect_equal(retained$form, diag(c(2, 2)), tolerance = 1e-12)
  expect_identical(retained$fit$groups, c(1L, 1L))
  expect_identical(pf_case(diag(c(2, 2 - 1e-6)), diag(2), 0, 1)$fit$rank_effective, 1L)
  expect_identical(pf_case(diag(c(2, 2)), diag(2), 0, 0)$fit$rank_effective, 0L)
})

test_that("[T23] controlled conditioning and signed-source stress agree with the dense reference", {
  set.seed(42901)
  for (r in c(1L, 3L, 7L)) for (condition in c(1, 1e4, 1e9)) {
    O <- qr.Q(qr(matrix(rnorm(r*r), r)))
    J <- O %*% diag(exp(seq(0, -log(condition), length.out = r)), r) %*% t(O)
    S <- matrix(rnorm(r*r), r); S <- (S + t(S))/2
    for (scale in c(1e-8, 1, 1e8)) {
      rank <- min(2L, r); penalty <- .03 * scale
      c <- pf_case(S * scale, J, penalty, rank)
      expected <- pf_dense(S * scale, J, rank, penalty)
      tol <- 1e-12 + 1e-6 * sqrt(pf_norm(S * scale))
      expect_lt(sqrt(pf_norm(c$form - expected)), tol)
      expect_true(all(is.finite(c$fit$factor)))
      expect_lt(max(abs(c$form - t(c$form))), tol)
      expect_gte(min(eigen(c$form, symmetric = TRUE)$values), -tol)
      expect_lte(c$fit$rank_effective, rank)
      P <- tcrossprod(c$pool$vectors)
      expect_lt(sqrt(pf_norm((diag(r) - P) %*% c$fit$form)), tol)
      expect_lt(abs(pf_norm(c$fit$form) - c$fit$prediction_norm_sq),
        1e-6 * max(1, pf_norm(c$fit$form)))
    }
  }
})

test_that("[T16 T23] numeric boundaries refuse before invalid spectral results can escape", {
  c <- pf_case(diag(c(2, 1)), diag(c(.8, .2)), .1, 1)
  for (penalty in list(-1, Inf, NA_real_, c(.1, .2))) {
    expect_error(crossform:::.geometry_fit_path(diag(2), c$pool, penalty), class = "effect_input_error")
  }
  for (S in list(matrix(1, 2, 3), matrix(NA_real_, 2, 2), matrix(c(1, .2, 0, 1), 2))) {
    expect_error(crossform:::.geometry_fit_path(S, c$pool, .1), class = "effect_input_error")
  }
  expect_error(crossform:::.geometry_fit_path(diag(2), c$pool, 1e308), "overflows")
})
