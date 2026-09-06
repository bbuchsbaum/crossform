# Independent algebra checks for the predictive-geometry design review.
# Base R only. This is evidence for a proposed estimator, not production code.
# Run: Rscript .planning/2026-09-04-predictive-geometry-checks.R
set.seed(20260904)
cat("Predictive geometry design checks; seed = 20260904\n")
sym <- function(x) (x + t(x)) / 2
norm2 <- function(x) sum(x * x)
inner <- function(x, y) sum(x * y)
near <- function(label, x, y, tol = 1e-10) {
  error <- max(abs(x - y))
  stopifnot(is.finite(error), error <= tol * max(1, abs(x), abs(y)))
  cat(sprintf("%-49s max error %.6g\n", label, error))
}
psd <- function(S, rank) {
  e <- eigen(sym(S), symmetric = TRUE)
  mu <- pmax(e$values, 0)
  if (rank < length(mu)) mu[seq.int(rank + 1L, length(mu))] <- 0
  list(form = e$vectors %*% (mu * t(e$vectors)), mu = mu, V = e$vectors)
}
support <- function(K) {
  e <- eigen(sym(K), symmetric = TRUE)
  keep <- e$values > 1e-10 * max(abs(e$values))
  list(Q = e$vectors[, keep, drop = FALSE], d = e$values[keep])
}
fit <- function(G, K, rank, penalty) {
  k <- support(K)
  r <- length(k$d)
  if (!r) return(matrix(0, nrow(G), ncol(G)))
  shifted <- crossprod(k$Q, G %*% k$Q) - diag(penalty / k$d, r)
  k$Q %*% psd(shifted, min(rank, r))$form %*% t(k$Q)
}
pinv <- function(T) {
  e <- svd(T)
  keep <- e$d > 1e-10 * max(e$d)
  e$v[, keep, drop = FALSE] %*%
    ((1 / e$d[keep]) * t(e$u[, keep, drop = FALSE]))
}
n <- 6L
V <- qr.Q(qr(cbind(1, diag(n))))[, -1L]
# Two model feature sets sharing exactly one direction.
T1 <- V[, 1:2] %*% diag(c(2, 0.6))
T2 <- V[, 2:3] %*% matrix(c(1, 0.2, 0.3, 1.4), 2L)
T1 <- T1 / sqrt(norm2(T1))
T2 <- T2 / sqrt(norm2(T2))
T <- cbind(T1, T2)
Q <- support(tcrossprod(T))$Q
R1 <- crossprod(Q, T1)
R2 <- crossprod(Q, T2)
G <- sym(matrix(rnorm(n * n), n)) + 2 * tcrossprod(V[, 1:3])
S <- crossprod(Q, G %*% Q)
O <- qr.Q(qr(matrix(rnorm(9), 3)))
cat(sprintf("Overlapping factors: %d columns, union rank %d\n", ncol(T), ncol(Q)))
for (w in list(c(0.3, 0.7), c(1, 0), c(0, 1))) {
  K <- w[1] * tcrossprod(T1) + w[2] * tcrossprod(T2)
  J <- w[1] * tcrossprod(R1) + w[2] * tcrossprod(R2)
  A <- fit(S, J, 2L, 0.15)
  F <- fit(G, K, 2L, 0.15)
  near(paste("pooled small/full, weights", paste(w, collapse = "/")),
    Q %*% A %*% t(Q), F)
  rotated <- fit(crossprod(O, S %*% O), crossprod(O, J %*% O), 2L, 0.15)
  near("orthogonal coordinate change", O %*% rotated %*% t(O), A)
  near("zero-rank baseline", fit(G, K, 0L, 0.15), matrix(0, n, n))
  near("kernel/penalty joint scaling", fit(G, 7 * K, 2L, 7 * 0.15), F)
  near("response/penalty joint scaling", fit(4 * G, K, 2L, 4 * 0.15), 4 * F)
}
K <- 0.3 * tcrossprod(T1) + 0.7 * tcrossprod(T2)
F <- fit(G, K, 2L, 0.15)
Ts <- cbind(sqrt(0.3) * T1, sqrt(0.7) * T2)
ef <- support(F)
Z <- ef$Q %*% diag(sqrt(ef$d), length(ef$d))
W <- pinv(Ts) %*% Z
near("redundant-feature minimum-norm extension", Ts %*% W, Z)
near("feature penalty / kernel penalty", norm2(W), inner(pinv(K), F))

# The loss identity is exact for a noncommuting kernel and response.
k <- support(K)
S <- crossprod(k$Q, G %*% k$Q)
L <- diag(0.15 / k$d, length(k$d))
X <- matrix(rnorm(length(k$d) * 2), length(k$d))
A <- tcrossprod(X)
constant <- 0.5 * (norm2(G) - norm2(S)) +
  0.5 * (norm2(S) - norm2(S - L))
near("completed square, noncommuting response/kernel",
  0.5 * norm2(G - k$Q %*% A %*% t(k$Q)) + inner(L, A),
  0.5 * norm2(A - (S - L)) + constant)
cat(sprintf("Commutator norm (nonzero): %.6g\n", sqrt(norm2(S %*% L - L %*% S))))

# Independent factor optimization is a numerical cross-check, not a proof.
objective <- function(z) {
  Z <- matrix(z, length(k$d), 2L)
  A <- tcrossprod(Z)
  0.5 * norm2(S - A) + inner(L, A)
}
exact <- psd(S - L, 2L)$form
exact_value <- 0.5 * norm2(S - exact) + inner(L, exact)
values <- replicate(12, optim(rnorm(length(k$d) * 2), objective,
  method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))$value)
stopifnot(min(values) >= exact_value - 1e-8)
near("best of 12 factor optimizations / spectral fit", min(values), exact_value, 1e-7)

# Full-span models coincide only at zero penalty, away from rank-cut ties.
K1 <- V %*% diag(c(0.5, 0.2, 0.15, 0.1, 0.05)) %*% t(V)
K2 <- V %*% diag(c(0.05, 0.1, 0.15, 0.2, 0.5)) %*% t(V)
near("full-span invariance at penalty zero", fit(G, K1, 2L, 0), fit(G, K2, 2L, 0))
gap <- sqrt(norm2(fit(G, K1, 2L, 0.15) - fit(G, K2, 2L, 0.15)))
stopifnot(gap > 0.01)
cat(sprintf("Full-span model difference with penalty: %.6g\n", gap))

# Pullback-adjoint and packed-codec identities without package internals.
Ge <- sym(matrix(rnorm(n * n), n))
small <- crossprod(Q, F %*% Q)
near("pullback-adjoint", inner(F, Ge), inner(small, crossprod(Q, Ge %*% Q)))
pack <- function(M) {
  multiplier <- matrix(sqrt(2), nrow(M), ncol(M)); diag(multiplier) <- 1
  (M * multiplier)[lower.tri(M, diag = TRUE)]
}
near("isometric packed dot product", sum(pack(F) * pack(Ge)), inner(F, Ge))
gain <- function(F, Gtest) 2 * inner(F, Gtest) - norm2(F)
mu <- ef$d
evidence <- diag(crossprod(ef$Q, Ge %*% ef$Q))
near("mode-wise gains sum to total", sum(2 * mu * evidence - mu^2), gain(F, Ge))
Gstar <- tcrossprod(V[, 1:2])
N <- sym(matrix(rnorm(n * n), n))
near("conditional risk identity (exact paired noise)",
  mean(c(gain(F, Gstar + N), gain(F, Gstar - N))),
  norm2(Gstar) - norm2(Gstar - F))
near("null gain is minus squared prediction norm",
  mean(c(gain(F, N), gain(F, -N))), -norm2(F))

# Both training and test geometries are averages of independent cross-products.
# No package executor or fitted-spectrum routine is used in this experiment.
reps <- 5000L
null_readings <- replicate(reps, {
  B <- replicate(4, matrix(rnorm(9), 3), simplify = FALSE)
  S_train <- sym(B[[1]] %*% t(B[[2]])) / 3
  S_test <- sym(B[[3]] %*% t(B[[4]])) / 3
  trained <- psd(S_train - diag(0.2, 3), 1L)$form
  c(gain = gain(trained, S_test), squared_norm = norm2(trained))
})
calibration <- null_readings["gain", ] + null_readings["squared_norm", ]
mcse <- sd(calibration) / sqrt(reps)
stopifnot(abs(mean(calibration)) < 5 * mcse)
cat(sprintf("Independent null (%d reps): gain %.6g; -norm %.6g; identity error %.6g; MCSE %.6g\n",
  reps, mean(null_readings["gain", ]), -mean(null_readings["squared_norm", ]),
  mean(calibration), mcse))
cat("All design checks passed. No production API was implemented.\n")
