# Independent executable court for `design/model-coordinate-geometry-contract.md`.
#
# This script deliberately does not load crossform.  It constructs the
# cross-partition geometry, a declared model basis, the lowered (model
# coordinate) relations, and the laws the contract states, from base matrix
# algebra only.  A child ticket of the epic compares the public package route
# against these values in `tests/testthat/test-model-coordinate-contract.R`.
#
# Sections mirror the contract:
#   O1  lowering commutes with the cross-partition estimator (total, coherent,
#       configuration)
#   O2  baseline invariance holds iff the basis is centered, and a thin QR of a
#       rank-deficient factor is not guaranteed to be centered
#   O3  the exact Frobenius solution (Theorem 4) and its Pythagorean split
#   O4  saturation (Theorem 3): a full-span basis reproduces any rank-s target
#   O5  the diagonal fit is a quadratic program in the model-only Gram, and it
#       needs only the compressed form
#   O6  the trace decomposition: addressable + orthogonal = centered total
#   O7  the compressed form is unbiased; the rank-s latent energy is not; the
#       edge-disjoint cross-fit restores unbiasedness under pure noise
#   O8  exact block coordinate descent on the block structure is monotone,
#       reaches a blockwise minimizer, and nests between the shared and the
#       diagonal fits

set.seed(20260904)
tol <- 1e-12
n <- 6L; p <- 12L; n_part <- 4L
conditions <- paste0("c", seq_len(n))
H <- diag(n) - 1 / n
frob <- function(x) sqrt(sum(x^2))

# --- model side --------------------------------------------------------------

gram <- function(D) -0.5 * H %*% D %*% H

# Rank-revealing model factor: keep only roots above tolerance, clamp the
# requested rank to the effective rank, and record the negative mass dropped.
model_factor <- function(D, rank, rel_tol = 1e-10) {
  e <- eigen(gram(D), symmetric = TRUE)
  keep <- which(e$values > rel_tol * max(abs(e$values)))
  keep <- keep[seq_len(min(rank, length(keep)))]
  list(
    T = e$vectors[, keep, drop = FALSE] %*%
      diag(sqrt(e$values[keep]), length(keep)),
    requested = rank, effective = length(keep),
    dropped_negative = sum(pmin(e$values, 0))
  )
}

# Rank-revealing orthonormal basis of col(T): T = Q R with Q'Q = I.
basis <- function(T, rel_tol = 1e-10) {
  s <- svd(T)
  keep <- s$d > rel_tol * s$d[[1L]]
  Q <- s$u[, keep, drop = FALSE]
  R <- diag(s$d[keep], sum(keep)) %*% t(s$v[, keep, drop = FALSE])
  list(Q = Q, R = R)
}

category <- rep(c(0, 1), each = n / 2)
D_category <- outer(category, category, function(a, b) (a - b)^2)
features <- matrix(rnorm(n * 3L), n, 3L)
D_feature <- as.matrix(dist(features))^2

m_category <- model_factor(D_category, rank = 2L)   # effective rank 1
m_feature <- model_factor(D_feature, rank = 2L)     # effective rank 2
stopifnot(m_category$effective == 1L, m_feature$effective == 2L)
T <- cbind(m_category$T, m_feature$T)
b <- basis(T); Q <- b$Q; R <- b$R
q <- ncol(Q)
stopifnot(
  max(abs(T - Q %*% R)) < tol,
  max(abs(crossprod(Q) - diag(q))) < tol,
  max(abs(colSums(Q))) < 1e-10      # centered: Q'1 = 0
)
P <- Q %*% t(Q)

# --- neural side -------------------------------------------------------------

# Planted signal: rank-2 inside the model span, plus one centered direction
# orthogonal to the model span, plus a condition-independent baseline.
W_star <- matrix(rnorm(ncol(T) * 2L), ncol(T), 2L)
Z_star <- T %*% W_star
u <- (H - P) %*% rnorm(n); u <- u / frob(u)
B_star <- Z_star %*% t(matrix(rnorm(p * 2L), p, 2L)) +
  u %*% t(rnorm(p)) + rep(1, n) %*% t(rnorm(p))
sigma <- 0.5
w <- runif(p, 0.5, 1.5)                # one frame node, positive weights
a <- sum(w)
draw_partitions <- function(B) {
  lapply(seq_len(n_part), function(r) B + matrix(rnorm(n * p, sd = sigma), n, p))
}
B <- draw_partitions(B_star)

# Undirected all-pairs pairing: six unordered edges of weight 1/6, each the two
# ordered half-edges of weight 1/12.
Gamma <- matrix(1 / 12, n_part, n_part); diag(Gamma) <- 0
stopifnot(abs(sum(Gamma) - 1) < tol)

geometry <- function(B, Gamma) {
  k <- nrow(B[[1L]])
  total <- coherent <- matrix(0, k, k)
  for (r in seq_along(B)) for (s in seq_along(B)) if (Gamma[r, s] != 0) {
    total <- total + Gamma[r, s] * B[[r]] %*% (w * t(B[[s]]))
    coherent <- coherent + Gamma[r, s] *
      tcrossprod(B[[r]] %*% w, B[[s]] %*% w) / a
  }
  list(total = total, coherent = coherent, configuration = total - coherent)
}

G <- geometry(B, Gamma)
lower <- function(B, Q) lapply(B, function(b) t(Q) %*% b)
G_m <- geometry(lower(B, Q), Gamma)

# --- O1: lowering commutes ---------------------------------------------------

o1 <- vapply(names(G), function(component) {
  max(abs(t(Q) %*% G[[component]] %*% Q - G_m[[component]]))
}, numeric(1))
stopifnot(all(o1 < tol))
cat("O1 lowering commutes (max abs diff):",
  paste(names(o1), formatC(o1, format = "e", digits = 2)), "\n")

# --- O2: baseline invariance -------------------------------------------------

shifted <- lapply(B, function(b) b + rep(1, n) %*% t(rnorm(p)))
G_shift <- geometry(shifted, Gamma)$total
o2_raw <- max(abs(G_shift - G$total))
o2_lowered <- max(abs(t(Q) %*% G_shift %*% Q - G_m$total))
stopifnot(o2_raw > 1e-3, o2_lowered < 1e-10)
cat("O2 baseline shift moves the raw form by", formatC(o2_raw, digits = 3),
  "and the compressed form by", formatC(o2_lowered, format = "e", digits = 2), "\n")

# The trap: the literal factor U_r Lambda_r^{1/2} with r above the effective
# rank has a zero column, and a thin QR of it is not a basis of col(T).
e_cat <- eigen(gram(D_category), symmetric = TRUE)
T_bad <- e_cat$vectors[, 1:2] %*% diag(sqrt(pmax(e_cat$values[1:2], 0)), 2L)
Q_bad <- qr.Q(qr(T_bad))
o2_trap_centering <- max(abs(colSums(Q_bad)))
o2_trap_shift <- max(abs(t(Q_bad) %*% G_shift %*% Q_bad -
  t(Q_bad) %*% G$total %*% Q_bad))
stopifnot(o2_trap_centering > 1e-8, o2_trap_shift > 1e-8)
cat("O2 trap: thin QR of the rank-deficient factor gives |Q'1| =",
  formatC(o2_trap_centering, digits = 3), "and baseline leakage",
  formatC(o2_trap_shift, digits = 3), "\n")

# --- O3: Theorem 4 -----------------------------------------------------------

S <- t(Q) %*% G$total %*% Q            # signed, indefinite in general
s_rank <- 1L
psd_rank <- function(S, s) {
  e <- eigen(S, symmetric = TRUE)
  keep <- seq_len(s)
  lambda <- pmax(e$values[keep], 0)
  e$vectors[, keep, drop = FALSE] %*% (lambda * t(e$vectors[, keep, drop = FALSE]))
}
A_star <- psd_rank(S, s_rank)
R_inv <- solve(R)
C_star <- R_inv %*% A_star %*% t(R_inv)
objective <- function(C) frob(G$total - T %*% C %*% t(T))^2
f_star <- objective(C_star)
constant <- frob(G$total - Q %*% S %*% t(Q))^2
o3_split <- abs(f_star - (constant + frob(S - A_star)^2))
stopifnot(o3_split < 1e-10)
random_feasible <- replicate(2000, {
  Wc <- matrix(rnorm(ncol(T) * s_rank), ncol(T), s_rank)
  objective(tcrossprod(Wc))
})
perturbed <- replicate(2000, {
  A_pert <- psd_rank(A_star + 0.05 * frob(A_star) * {
    E <- matrix(rnorm(q * q), q, q); (E + t(E)) / 2
  }, s_rank)
  objective(R_inv %*% A_pert %*% t(R_inv))
})
stopifnot(f_star <= min(random_feasible) + 1e-10,
  f_star <= min(perturbed) + 1e-10)
cat("O3 Pythagorean split error", formatC(o3_split, format = "e", digits = 2), "\n")
cat("O3 closed form objective", formatC(f_star, digits = 6),
  "vs best of 2000 random feasible", formatC(min(random_feasible), digits = 6),
  "and best of 2000 local perturbations", formatC(min(perturbed), digits = 6), "\n")
cat("O3 signed compressed spectrum:", formatC(eigen(S, symmetric = TRUE)$values, digits = 4), "\n")

# --- O4: saturation ----------------------------------------------------------

Q_full <- eigen(H, symmetric = TRUE)$vectors[, seq_len(n - 1L)]  # basis of 1-perp
stopifnot(max(abs(colSums(Q_full))) < 1e-10)
L <- H %*% matrix(rnorm(n * 2L), n, 2L)
G_target <- tcrossprod(L)                # any centered PSD rank-2 target
fit_full <- Q_full %*% psd_rank(t(Q_full) %*% G_target %*% Q_full, 2L) %*% t(Q_full)
fit_model <- Q %*% psd_rank(t(Q) %*% G_target %*% Q, 2L) %*% t(Q)
o4_full <- frob(G_target - fit_full)
o4_model <- frob(G_target - fit_model)
stopifnot(o4_full < 1e-10, o4_model > 1e-3)
df <- function(q, s) q * s - s * (s - 1) / 2
cat("O4 saturation: full-span residual", formatC(o4_full, format = "e", digits = 2),
  "; model-span residual", formatC(o4_model, digits = 3),
  "; df(q=", q, ", s=2) =", df(q, 2L), "vs unconstrained", df(n - 1L, 2L), "\n")

# --- O5: the diagonal fit ----------------------------------------------------

g <- vapply(seq_len(ncol(T)), function(j) drop(t(T[, j]) %*% G$total %*% T[, j]), numeric(1))
g_from_S <- vapply(seq_len(ncol(T)), function(j) drop(t(R[, j]) %*% S %*% R[, j]), numeric(1))
M <- crossprod(T)^2
c_try <- runif(ncol(T))
lhs <- 0.5 * frob(G$total - T %*% (c_try * t(T)))^2
rhs <- 0.5 * frob(G$total)^2 - sum(g * c_try) + 0.5 * drop(t(c_try) %*% M %*% c_try)
stopifnot(abs(lhs - rhs) < 1e-10, max(abs(g - g_from_S)) < 1e-10)
cat("O5 diagonal objective identity error", formatC(abs(lhs - rhs), format = "e", digits = 2),
  "; g from compressed form error", formatC(max(abs(g - g_from_S)), format = "e", digits = 2), "\n")

# --- O6: trace decomposition -------------------------------------------------

centered_total <- sum(diag(H %*% G$total %*% H))
addressable <- sum(diag(P %*% G$total))
orthogonal <- sum(diag((H - P) %*% G$total))
stopifnot(abs(centered_total - addressable - orthogonal) < 1e-10,
  abs(addressable - sum(diag(S))) < 1e-10)
cat("O6 centered total", formatC(centered_total, digits = 5), "= addressable",
  formatC(addressable, digits = 5), "+ orthogonal", formatC(orthogonal, digits = 5),
  "; addressable = tr(S) to", formatC(abs(addressable - sum(diag(S))), format = "e", digits = 1), "\n")

# --- O7: unbiasedness, latent bias, and the edge-disjoint cross-fit ----------

reps <- 600L
edge <- function(r, s) { E <- matrix(0, n_part, n_part); E[r, s] <- E[s, r] <- 0.5; E }
disjoint <- list(c(1, 2, 3, 4), c(1, 3, 2, 4), c(1, 4, 2, 3))
mc <- function(B_true) {
  out <- matrix(NA_real_, reps, 3L,
    dimnames = list(NULL, c("addressable", "latent_rank1", "cross_fit_rank1")))
  for (i in seq_len(reps)) {
    Bi <- lower(draw_partitions(B_true), Q)
    Si <- geometry(Bi, Gamma)$total
    cf <- mean(vapply(disjoint, function(d) {
      S_train <- geometry(Bi, edge(d[[3]], d[[4]]))$total
      v <- eigen(S_train, symmetric = TRUE)$vectors[, 1L]
      S_eval <- geometry(Bi, edge(d[[1]], d[[2]]))$total
      drop(t(v) %*% S_eval %*% v)
    }, numeric(1)))
    out[i, ] <- c(sum(diag(Si)), sum(diag(psd_rank(Si, 1L))), cf)
  }
  out
}
null <- mc(matrix(0, n, p))
null_mean <- colMeans(null); null_se <- apply(null, 2, sd) / sqrt(reps)
stopifnot(
  abs(null_mean[["addressable"]]) < 4 * null_se[["addressable"]],
  null_mean[["latent_rank1"]] > 4 * null_se[["latent_rank1"]],
  abs(null_mean[["cross_fit_rank1"]]) < 4 * null_se[["cross_fit_rank1"]]
)
cat("O7 pure noise, mean (se):",
  sprintf("%s %.4f (%.4f)", names(null_mean), null_mean, null_se), "\n")
signal <- mc(B_star)
truth <- t(Q) %*% B_star %*% (w * t(B_star)) %*% Q
o7_bias <- mean(signal[, "addressable"]) - sum(diag(truth))
stopifnot(abs(o7_bias) < 4 * sd(signal[, "addressable"]) / sqrt(reps))
cat("O7 with signal: addressable energy bias", formatC(o7_bias, digits = 3),
  "(se", formatC(sd(signal[, "addressable"]) / sqrt(reps), digits = 3), ") against truth",
  formatC(sum(diag(truth)), digits = 5), "\n")

# --- O8: block coordinate descent on the block structure ----------------------
#
# Contract section 4.2. With Q_i an orthonormal basis of col(R[, block i]) in
# Q-coordinates, the exact block update is
#   A_i <- [Q_i' (S - sum_{j != i} Q_j A_j Q_j') Q_i]_{+, s_i},
# the minimizer of the objective over PSD rank-s_i A_i with the other blocks
# held fixed. Three things are measured: the objective never increases across
# updates; the sweep converges to a point every block update leaves fixed (a
# blockwise minimizer); and the nesting of the feasible sets shows in the
# objectives, shared <= block <= diagonal when each block's rank budget is
# its own dimension and the shared budget is q.

blocks <- list(category = 1L, feature = 2:3)
block_basis <- lapply(blocks, function(cols) {
  sv <- svd(R[, cols, drop = FALSE])
  sv$u[, sv$d > 1e-10 * sv$d[[1L]], drop = FALSE]
})
block_objective <- function(A) {
  fit <- Reduce(`+`, Map(function(Qi, Ai) Qi %*% Ai %*% t(Qi), block_basis, A))
  frob(S - fit)^2
}
block_descent <- function(S, s, sweeps = 20000L, tol = 1e-11) {
  A <- lapply(block_basis, function(Qi) matrix(0, ncol(Qi), ncol(Qi)))
  history <- block_objective(A)
  for (sweep in seq_len(sweeps)) {
    moved <- 0
    for (i in seq_along(A)) {
      others <- Reduce(`+`, Map(function(Qj, Aj) Qj %*% Aj %*% t(Qj),
        block_basis[-i], A[-i]), matrix(0, q, q))
      updated <- psd_rank(t(block_basis[[i]]) %*% (S - others) %*% block_basis[[i]],
        s[[i]])
      moved <- max(moved, max(abs(updated - A[[i]])))
      A[[i]] <- updated
      history <- c(history, block_objective(A))
    }
    if (moved < tol * max(1, max(abs(S)))) break
  }
  list(A = A, history = history, sweeps = sweep, moved = moved)
}
o8 <- block_descent(S, s = c(1L, 2L))
o8_increase <- max(diff(o8$history))
stopifnot(o8_increase <= 1e-10 * max(1, o8$history[[1L]]))
# Fixed point: one more exact update of every block changes nothing.
o8_refit <- block_descent_refit <- {
  A <- o8$A
  moved <- 0
  for (i in seq_along(A)) {
    others <- Reduce(`+`, Map(function(Qj, Aj) Qj %*% Aj %*% t(Qj),
      block_basis[-i], A[-i]), matrix(0, q, q))
    updated <- psd_rank(t(block_basis[[i]]) %*% (S - others) %*% block_basis[[i]],
      c(1L, 2L)[[i]])
    moved <- max(moved, max(abs(updated - A[[i]])))
  }
  moved
}
stopifnot(o8_refit < 1e-8 * max(1, max(abs(S))))
# Nesting: the shared fit (rank q) is at least as good as the block fit, and
# the block fit at least as good as the nonnegative diagonal fit, because
# diagonal C is block-diagonal with diagonal blocks and every block-diagonal
# PSD C is a PSD C.
f_shared_full <- frob(S - psd_rank(S, q))^2
# Diagonal fit by projected gradient on the O5 quadratic program; the
# objective is 0.5 * ||G - T diag(c) T'||^2, reported here in S-coordinates
# as the same residual minus the constant outside the model span.
constant_out <- frob(G$total - Q %*% S %*% t(Q))^2
diag_c <- rep(0, ncol(T))
step <- 1 / max(eigen(M, symmetric = TRUE)$values)
for (it in seq_len(20000L)) diag_c <- pmax(diag_c - step * (M %*% diag_c - g), 0)
f_diag <- 0.5 * frob(G$total - T %*% (as.numeric(diag_c) * t(T)))^2 * 2 - constant_out
f_block <- o8$history[[length(o8$history)]]
stopifnot(f_shared_full <= f_block + 1e-8, f_block <= f_diag + 1e-8)
cat("O8 block descent: objective", formatC(o8$history[[1L]], digits = 6), "->",
  formatC(f_block, digits = 6), "in", o8$sweeps, "sweeps; max increase",
  formatC(o8_increase, format = "e", digits = 1), "; fixed-point move",
  formatC(o8_refit, format = "e", digits = 1), "\n")
cat("O8 nesting in S-coordinates: shared(q)", formatC(f_shared_full, digits = 6),
  "<= block", formatC(f_block, digits = 6), "<= diagonal", formatC(f_diag, digits = 6), "\n")

cat("model-coordinate-geometry oracle: all laws hold\n")
