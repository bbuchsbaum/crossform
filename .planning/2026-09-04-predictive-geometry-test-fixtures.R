# Exact fixtures for the proposed predictive-geometry test specification.
# Base R only; no production estimator and no package internals.
# Run: Rscript .planning/2026-09-04-predictive-geometry-test-fixtures.R
# This checks the fixture arithmetic and sensitivity of selected assertions.
# It does NOT establish that a future implementation passes these tests.
cat("Predictive geometry exact fixtures\n", R.version.string, "\n", sep = "")
norm2 <- function(A) sum(A * A)
inner <- function(A, B) sum(A * B)
gain <- function(A, S) 2 * inner(A, S) - norm2(A)
fixtures <- list()
mutants <- 0L
near <- function(label, actual, expected, tol = 1e-12) {
  error <- max(abs(actual - expected))
  stopifnot(is.finite(error), error <= tol * max(1, abs(expected)))
  cat(sprintf("PASS %-47s error %.3g\n", label, error))
}
separates <- function(label, correct, wrong) {
  gap <- max(abs(correct - wrong))
  stopifnot(is.finite(gap), gap > 1e-5)
  mutants <<- mutants + 1L
  cat(sprintf("KILL %-47s gap %.6g\n", label, gap))
}

# F01: inverse-eigenvalue penalty, rank path, and objective arithmetic.
S <- diag(c(6, 3, -1)); J <- diag(c(0.5, 0.25, 0.25)); penalty <- 0.5
shifted <- diag(c(5, 1, -3)); A <- diag(c(5, 1, 0))
fixtures$F01 <- list(S = S, J = J, penalty = penalty, shifted = shifted,
  rank_path = list(matrix(0, 3, 3), diag(c(5, 0, 0)), A, A))
near("F01 inverse penalty", S - penalty * solve(J), shifted)
near("F01 objective", 0.5 * norm2(S - A) + penalty * inner(solve(J), A), 10)
separates("penalize by D instead of inverse D", shifted, S - penalty * J)
separates("wrong regularization sign", shifted, S + penalty * solve(J))

# F02: pseudoinverse alone is not a range constraint.
S <- diag(c(2, 7)); J <- diag(c(1, 0)); A <- diag(c(1, 0))
fixtures$F02 <- list(S = S, J = J, penalty = 1, rank = 2L, A = A)
separates("fit unsupported null-space directions", A, diag(c(1, 7)))

# Closed form for the rank-one positive part of an indefinite symmetric 2x2.
# No eigendecomposition in the oracle.
positive2 <- function(B) {
  a <- B[1, 1]; b <- B[1, 2]; d <- B[2, 2]
  gap <- sqrt((a - d)^2 + 4 * b^2)
  high <- (a + d + gap) / 2; low <- (a + d - gap) / 2
  stopifnot(high > 0, low < 0)
  high * (B - low * diag(2)) / gap
}
# F03: noncommuting source and penalty; preclipping the source is wrong.
S <- matrix(c(1, 2, 2, 1), 2)
J <- diag(c(0.8, 0.2)); penalty <- 0.2
shifted <- S - diag(c(0.25, 1))
A <- positive2(shifted)
fixtures$F03 <- list(S = S, J = J, penalty = penalty, rank = 1L, A = A)
near("F03 characteristic equation", shifted %*% A, sum(diag(A)) * A)
near("F03 rank one", det(A), 0)
wrong <- positive2(matrix(1.5, 2, 2) - diag(c(0.25, 1)))
separates("clip compressed source before penalty", A, wrong)
# Numerical eigensolver is only the candidate, not the oracle.
ev <- eigen(shifted, symmetric = TRUE)
candidate <- ev$values[1] * tcrossprod(ev$vectors[, 1, drop = FALSE])
near("F03 spectral candidate / analytic oracle", candidate, A)

# F04: clipping full neural geometry before model compression changes it.
fixtures$F04 <- list(G = S, Q = matrix(c(1, 0), 2), expected_S = matrix(1),
  prematurely_clipped_S = matrix(1.5))
separates("clip full source before compression", 1, 1.5)

# F05: amplitude-aware gain, with trace and Frobenius norm distinguished.
A <- diag(c(2, 1)); S <- diag(c(1, -1))
fixtures$F05 <- list(A = A, S = S, inner = 1, squared_norm = 5, gain = -3)
near("F05 predictive gain", gain(A, S), -3)
separates("omit prediction norm", gain(A, S), 2 * inner(A, S))
separates("use trace in place of squared norm", gain(A, S),
  2 * inner(A, S) - sum(diag(A)))
separates("use trace squared in place of squared norm", gain(A, S),
  2 * inner(A, S) - sum(diag(A))^2)
separates("omit leading factor two", gain(A, S), inner(A, S) - norm2(A))

# F06: off-diagonal codec scaling. Independent loop-free packing oracle.
A <- matrix(c(1, 0.5, 0.5, 1), 2); S <- matrix(c(0, 1, 1, 0), 2)
pack2 <- function(B) c(B[1, 1], sqrt(2) * B[2, 1], B[2, 2])
fixtures$F06 <- list(A = A, S = S, inner = 1, squared_norm = 2.5, gain = -0.5)
near("F06 packed/Frobenius inner product", sum(pack2(A) * pack2(S)), 1)
near("F06 gain", gain(A, S), -0.5)
separates("omit sqrt(2) packed multiplier", inner(A, S),
  sum(A[lower.tri(A, diag = TRUE)] * S[lower.tri(S, diag = TRUE)]))

# F07: held-out gain need not improve as rank increases.
A <- diag(c(2, 1)); S <- diag(c(2, 0))
fixtures$F07 <- list(A = A, S = S, mode_gain = c(4, -1), gain = 3)
near("F07 mode sum", sum(2 * diag(A) * diag(S) - diag(A)^2), 3)
separates("discard negative modes using test outcomes", gain(A, S), 4)
fixtures$F07$overestimated_mode <- c(amplitude = 3, evidence = 1, gain = -3)
near("F07 reproducible but overestimated mode", 2 * 3 * 1 - 3^2, -3)

# F08: coherent/configuration evidence adds, prediction cost is counted once.
A <- diag(c(2, 1)); coherent <- diag(c(1, 0)); configuration <- diag(c(0, -1))
fixtures$F08 <- list(A = A, coherent = coherent, configuration = configuration,
  gain = -3)
near("F08 additive linear evidence", inner(A, coherent + configuration),
  inner(A, coherent) + inner(A, configuration))
separates("subtract prediction cost for each component", gain(A, coherent + configuration),
  gain(A, coherent) + gain(A, configuration))

# F09: exact conditional risk expectation for a finite mean-zero noise law.
truth <- matrix(c(2, 0.5, 0.5, 1), 2)
N <- matrix(c(1, 2, 2, -3), 2)
fixtures$F09 <- list(A = A, truth = truth, noise_support = list(N, -N))
near("F09 conditional risk identity", mean(c(gain(A, truth + N), gain(A, truth - N))),
  norm2(truth) - norm2(truth - A))
near("F09 nonzero predictor under null", mean(c(gain(A, N), gain(A, -N))), -5)

# F10: unit trace of the retained kernel, not the original kernel.
fixtures$F10 <- list(original = diag(c(8, 2)), rank = 1L,
  retained = diag(c(8, 0)), normalized = diag(c(1, 0)))
separates("normalize by original pre-truncation trace", diag(c(1, 0)), diag(c(0.8, 0)))

# F11: compression cannot identify whole-space residual magnitude.
G1 <- diag(c(1, 0)); G2 <- diag(c(1, 4)); F <- diag(c(1, 0))
fixtures$F11 <- list(G1 = G1, G2 = G2, F = F, residuals = c(0, 16))
near("F11 equal compressed input", G1[1, 1], G2[1, 1])
near("F11 unequal whole-space residuals", c(norm2(G1 - F), norm2(G2 - F)), c(0, 16))

# F12: nonlinear fitting does not preserve component additivity.
S1 <- diag(c(1, -2)); S2 <- diag(c(-2, 1))
fixtures$F12 <- list(S1 = S1, S2 = S2, separate_sum = diag(2), joint = matrix(0, 2, 2))
separates("assume PSD fits add", fixtures$F12$joint, fixtures$F12$separate_sum)

# F13: same full span, different preferences at positive penalty.
fixtures$F13 <- list(S = diag(c(3, 3)), J1 = diag(c(0.8, 0.2)),
  J2 = diag(c(0.2, 0.8)), penalty = 0.4, A1 = diag(c(2.5, 1)), A2 = diag(c(1, 2.5)))
separates("discard model eigenvalues and retain only span", fixtures$F13$A1, fixtures$F13$A2)

# F14: tied rank boundary, equal training loss but unequal test evidence.
B <- diag(c(2, 2)); A1 <- diag(c(2, 0)); A2 <- diag(c(0, 2)); S <- diag(c(1, 0))
fixtures$F14 <- list(shifted = B, rank = 1L, A1 = A1, A2 = A2, test = S)
near("F14 equal optimum objectives", norm2(B - A1), norm2(B - A2))
separates("assume tied rank-cut prediction is unique", gain(A1, S), gain(A2, S))
cat(sprintf("Validated %d fixtures and %d deliberately wrong alternatives.\n", length(fixtures), mutants))
cat("These are design witnesses; production tests remain to be implemented.\n")
