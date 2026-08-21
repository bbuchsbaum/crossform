#!/usr/bin/env Rscript
# Oracle: the subject-Gram trick and the consensus / heterogeneity split
# ======================================================================
#
# Supports claims 5 and 6 of `design/population-form-contract.md`.
#
#   P5.a  The packed codec must carry sqrt(2) on off-diagonals to be
#         Frobenius consistent; a naive `upper.tri` packing is NOT, and it
#         changes the population spectrum.
#   P5.b  Geometry-space covariance over N subjects has rank <= N-1, and the
#         N x N subject Gram reproduces its nonzero spectrum EXACTLY.
#   P5.c  The eigenvector (the geometry mode) is recovered from the Gram
#         eigenvector by C^T u, so the D x D covariance is never formed.
#   P5.d  The Gram is additive over group nodes: Gram = sum_j Gram_j, so it
#         streams node by node and the N x D stack is never materialized.
#   P6.a  V = Q^C + Q^H is an EXACT algebraic identity (mean square = squared
#         mean + variance), verified to machine precision.
#   P6.b  The plug-in estimator of that identity is BIASED: within-subject
#         sampling noise lands in Q^H.  Monte Carlo bias measured.
#   P6.c  The cross-fitted split (cross-partition for Q^W, cross-participant
#         for Q^C) is unbiased for Sigma_B and can be indefinite.  Measured.
#   P6.d  Subject loadings are read off the Gram modes; recovery measured.
#   P6.e  A PSD projection of an indefinite Q^H moves mass; the amount is
#         measured and must be recorded (contract G7 discipline).
#   P6.f  The CROSS-FITTED subject Gram is the object P5's covariance trick
#         must be built on: unbiased for the noiseless subject Gram, still
#         N x N and node-additive, indefinite by construction, and its
#         diagonal / off-diagonal contrast IS the V^W / V^C split of P7.
#
# Pure matrix algebra; the crossform package is NOT loaded.
#
# Run:  Rscript design/oracles/population-geometry-split.R

options(digits = 17)

say <- function(...) cat(..., "\n", sep = "")
rule <- function(title) say("\n", title, "\n", strrep("-", nchar(title)))
sci <- function(x) format(x, scientific = TRUE, digits = 4)

set.seed(20260820)

svec_index <- function(q) {
  A <- matrix(0, q, q)
  idx <- which(upper.tri(A, diag = TRUE), arr.ind = TRUE)
  list(idx = idx, w = ifelse(idx[, 1L] == idx[, 2L], 1, sqrt(2)))
}
svec <- function(A, ix = svec_index(nrow(A))) A[ix$idx] * ix$w
svec_naive <- function(A, ix = svec_index(nrow(A))) A[ix$idx]
unsvec <- function(v, q, ix = svec_index(q)) {
  A <- matrix(0, q, q)
  A[ix$idx] <- v / ix$w
  A <- A + t(A)
  diag(A) <- diag(A) / 2
  A
}
sym_draw <- function(q) {
  m <- matrix(rnorm(q * q), q, q)
  0.5 * (m + t(m))
}

# ---------------------------------------------------------------------------
# P5.a  the sqrt(2) is load-bearing
# ---------------------------------------------------------------------------
rule("P5.a  sqrt(2) off-diagonals are required for Frobenius consistency")

q <- 6L
ix <- svec_index(q)
Pdim <- q * (q + 1L) / 2L
worst_ok <- 0; worst_naive <- 0
for (i in seq_len(200L)) {
  A <- sym_draw(q); B <- sym_draw(q)
  rhs <- sum(A * B)
  worst_ok <- max(worst_ok, abs(sum(svec(A, ix) * svec(B, ix)) - rhs) /
    max(1, abs(rhs)))
  worst_naive <- max(worst_naive, abs(sum(svec_naive(A, ix) * svec_naive(B, ix)) -
    rhs) / max(1, abs(rhs)))
}
say("  q = ", q, ", packed width = ", Pdim, " (vs ", q * q, " rectangular)")
say("  worst relative |<svec A, svec B> - <A,B>_F|, sqrt(2) codec = ",
  sci(worst_ok))
say("  worst relative |<naive A, naive B> - <A,B>_F|              = ",
  sci(worst_naive))
say("  round trip: max |unsvec(svec(A)) - A| = ",
  sci(max(abs(unsvec(svec(A, ix), q, ix) - A))))

# ---------------------------------------------------------------------------
# P5.b  rank <= N-1 and the subject-Gram spectrum
# ---------------------------------------------------------------------------
rule("P5.b  geometry covariance has rank <= N-1; the subject Gram recovers it")

N <- 5L
mg <- 4L                                   # group nodes carrying transported forms
D <- mg * Pdim                             # full population coordinate width
pack_subject <- function(codec) {
  # one subject's transported forms at mg group nodes, concatenated
  as.numeric(vapply(seq_len(mg), function(j) codec(sym_draw(q), ix),
    numeric(Pdim)))
}
Y <- t(vapply(seq_len(N), function(i) pack_subject(svec), numeric(D)))
C <- sweep(Y, 2L, colMeans(Y))
Sigma <- crossprod(C) / (N - 1L)           # D x D
Gram <- tcrossprod(C) / (N - 1L)           # N x N
ev_s <- sort(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values,
  decreasing = TRUE)
ev_g <- sort(eigen(Gram, symmetric = TRUE, only.values = TRUE)$values,
  decreasing = TRUE)
say("  N = ", N, "   group nodes = ", mg, "   D = ", D)
say("  rank(Sigma)                          = ", qr(Sigma)$rank,
  "   (bound N-1 = ", N - 1L, ")")
say("  max |lambda_k(Sigma) - lambda_k(Gram)|, k = 1..", N - 1L, " : ",
  sci(max(abs(ev_s[seq_len(N - 1L)] - ev_g[seq_len(N - 1L)]))))
say("  max |lambda_k(Sigma)| for k > N-1                     : ",
  sci(max(abs(ev_s[N:D]))))
say("  eigen cost: ", N, "^3 instead of ", D, "^3  (ratio ",
  sprintf("%.3g", (D / N)^3), ")")

# The naive packing changes the answer.
Yn <- t(vapply(seq_len(N), function(i) pack_subject(svec_naive), numeric(D)))
Cn <- sweep(Yn, 2L, colMeans(Yn))
ev_n <- sort(eigen(tcrossprod(Cn) / (N - 1L), symmetric = TRUE,
  only.values = TRUE)$values, decreasing = TRUE)
say("  top eigenvalue, sqrt(2) codec = ", sprintf("%.9f", ev_g[1L]),
  " ; naive codec = ", sprintf("%.9f", ev_n[1L]))
say("  -> the naive packing under-weights every off-diagonal by sqrt(2), so",
  " its `covariance` is not a covariance of forms at all.")

# ---------------------------------------------------------------------------
# P5.c  recovering the geometry mode without forming Sigma
# ---------------------------------------------------------------------------
rule("P5.c  the mode is C^T u, so the D x D covariance is never formed")

eg <- eigen(Gram, symmetric = TRUE)
u1 <- eg$vectors[, 1L]; lam1 <- eg$values[1L]
v1 <- as.numeric(crossprod(C, u1)); v1 <- v1 / sqrt(sum(v1^2))
es <- eigen(Sigma, symmetric = TRUE)
say("  |lambda_1(Gram) - lambda_1(Sigma)|          = ",
  sci(abs(lam1 - es$values[1L])))
say("  max |Sigma v1 - lambda_1 v1|                = ",
  sci(max(abs(as.numeric(Sigma %*% v1) - lam1 * v1))))
say("  |<v1, first eigenvector of Sigma>|          = ",
  sprintf("%.15f", abs(sum(v1 * es$vectors[, 1L]))))
mode1 <- unsvec(v1[seq_len(Pdim)], q, ix)
say("  the mode unpacks to a symmetric q x q form: max |M - M^T| = ",
  sci(max(abs(mode1 - t(mode1)))))
say("  subject loadings on mode 1 = ",
  paste(sprintf("%+.4f", u1), collapse = " "))

# ---------------------------------------------------------------------------
# P5.d  the Gram is additive over group nodes
# ---------------------------------------------------------------------------
rule("P5.d  Gram = sum_j Gram_j  (streaming, node by node)")

Gsum <- matrix(0, N, N)
for (j in seq_len(mg)) {
  cols <- ((j - 1L) * Pdim + 1L):(j * Pdim)
  Gsum <- Gsum + tcrossprod(C[, cols, drop = FALSE]) / (N - 1L)
}
say("  max |Gram - sum_j Gram_j| = ", sci(max(abs(Gram - Gsum))))
say("  -> the population layer accumulates an N x N matrix while streaming")
say("     nodes; it never materializes the N x ", D, " stack.")

# ---------------------------------------------------------------------------
# P6.a  V = Q^C + Q^H is an exact identity
# ---------------------------------------------------------------------------
rule("P6.a  V = Q^C + Q^H is exact (mean square = squared mean + variance)")

V <- crossprod(Y) / N                              # (1/N) sum_i Y_i Y_i^T
Ybar <- colMeans(Y)
QC <- tcrossprod(Ybar)                             # consensus
QH <- crossprod(C) / N                             # heterogeneity, 1/N divisor
say("  max |V - (Q^C + Q^H)|              = ", sci(max(abs(V - QC - QH))))
say("  trace V = ", sprintf("%.12f", sum(diag(V))),
  " = ", sprintf("%.12f", sum(diag(QC))), " + ",
  sprintf("%.12f", sum(diag(QH))))
say("  consensus share tr(Q^C)/tr(V)      = ",
  sprintf("%.9f", sum(diag(QC)) / sum(diag(V))))

# ---------------------------------------------------------------------------
# P6.b / P6.c  the identity is exact, the plug-in estimator is not
# ---------------------------------------------------------------------------
rule("P6.b/c  plug-in split is biased; the cross-fitted split is not")

# Planted population: Y_ir = mu + U_i + eps_ir, two independent runs r.
Nsub <- 12L
Dp <- Pdim
mu_true <- svec(sym_draw(q), ix)
Lb <- matrix(rnorm(Dp * 2L), Dp, 2L) * 0.6         # rank-2 between-subject
SigmaB <- tcrossprod(Lb)
sd_eps <- 0.8
trB <- sum(diag(SigmaB))
trE <- Dp * sd_eps^2

one_rep <- function() {
  Uload <- matrix(rnorm(Nsub * 2L), Nsub, 2L)
  U <- Uload %*% t(Lb)                             # Nsub x Dp
  Y1 <- sweep(U, 2L, mu_true, "+") + matrix(rnorm(Nsub * Dp, sd = sd_eps), Nsub, Dp)
  Y2 <- sweep(U, 2L, mu_true, "+") + matrix(rnorm(Nsub * Dp, sd = sd_eps), Nsub, Dp)
  Ybar_i <- (Y1 + Y2) / 2
  Cc <- sweep(Ybar_i, 2L, colMeans(Ybar_i))
  QH_plug <- crossprod(Cc) / (Nsub - 1L)           # unbiased-divisor plug-in
  # cross-fitted: within = cross-partition, consensus = cross-participant
  QW <- (crossprod(Y1, Y2) + crossprod(Y2, Y1)) / (2 * Nsub)
  S <- colSums(Ybar_i)
  cross <- crossprod(Ybar_i, matrix(S, Nsub, Dp, byrow = TRUE)) -
    crossprod(Ybar_i)                              # sum_{i != i'} Y_i Y_i'^T
  QC_cf <- (cross + t(cross)) / (2 * Nsub * (Nsub - 1L))
  QH_cf <- QW - QC_cf
  list(plug = sum(diag(QH_plug)), cf = sum(diag(QH_cf)),
       cfC = sum(diag(QC_cf)), cfW = sum(diag(QW)),
       min_ev = min(eigen(QH_cf, symmetric = TRUE, only.values = TRUE)$values),
       Uload = Uload, QH_cf = QH_cf, Cc = Cc)
}
reps <- 2000L
res <- replicate(reps, one_rep()[c("plug", "cf", "cfC", "cfW", "min_ev")],
  simplify = FALSE)
plug <- vapply(res, `[[`, numeric(1), "plug")
cf <- vapply(res, `[[`, numeric(1), "cf")
mev <- vapply(res, `[[`, numeric(1), "min_ev")
say("  planted: tr(Sigma_B) = ", sprintf("%.6f", trB),
  " ; tr(Sigma_eps) = ", sprintf("%.6f", trE),
  " ; within-subject mean of 2 runs halves it -> ",
  sprintf("%.6f", trE / 2))
say("  Monte Carlo over ", reps, " replications, N = ", Nsub, ", D = ", Dp, ":")
say(sprintf("    mean tr(Q^H) plug-in     = %+.6f   bias = %+.6f  (predicted %+.6f)",
  mean(plug), mean(plug) - trB, trE / 2))
say(sprintf("    mean tr(Q^H) cross-fitted= %+.6f   bias = %+.6f  (MC se %.4f)",
  mean(cf), mean(cf) - trB, sd(cf) / sqrt(reps)))
say("    plug-in inflation of the heterogeneity trace = ",
  sprintf("%+.1f%%", 100 * (mean(plug) - trB) / trB))
say("  -> within-subject sampling noise is entirely booked as heterogeneity by")
say("     the plug-in split.  Only the cross-fitted split estimates Sigma_B.")
say("  cross-fitted Q^H is indefinite in ",
  sprintf("%.1f%%", 100 * mean(mev < 0)), " of replications; ",
  "worst min eigenvalue = ", sprintf("%+.6f", min(mev)))

# ---------------------------------------------------------------------------
# P6.d  subject loadings from the Gram modes
# ---------------------------------------------------------------------------
rule("P6.d  subject loadings are read off the Gram modes")

set.seed(4242)
r1 <- one_rep()
Gr <- tcrossprod(r1$Cc) / (Nsub - 1L)
eg2 <- eigen(Gr, symmetric = TRUE)
plane <- r1$Uload                                   # planted 2-d loading plane
fit_r2 <- function(u) summary(lm(u ~ plane))$r.squared
say("  planted heterogeneity is rank 2; recovered Gram spectrum (top 4) = ",
  paste(sprintf("%.4f", eg2$values[1:4]), collapse = " "))
say("  R^2 of Gram mode 1 on the planted loading plane = ",
  sprintf("%.6f", fit_r2(eg2$vectors[, 1L])))
say("  R^2 of Gram mode 2 on the planted loading plane = ",
  sprintf("%.6f", fit_r2(eg2$vectors[, 2L])))
say("  R^2 of Gram mode 3 (noise) on the same plane    = ",
  sprintf("%.6f", fit_r2(eg2$vectors[, 3L])))

# ---------------------------------------------------------------------------
# P6.e  what a PSD projection of Q^H costs
# ---------------------------------------------------------------------------
rule("P6.e  a PSD projection moves mass, and the amount must be recorded")

QH_cf <- r1$QH_cf
e <- eigen(QH_cf, symmetric = TRUE)
neg <- e$values[e$values < 0]
QH_psd <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors)
say("  eigenvalues of the cross-fitted Q^H: ", sum(e$values < 0), " of ",
  length(e$values), " negative; min = ", sprintf("%+.6f", min(e$values)))
say("  tr(Q^H) signed   = ", sprintf("%+.9f", sum(diag(QH_cf))))
say("  tr(Q^H) PSD-clip = ", sprintf("%+.9f", sum(diag(QH_psd))))
say("  moved mass = sum|negative eigenvalues| = ", sprintf("%.9f", sum(abs(neg))),
  "   (", sprintf("%+.2f%%", 100 * sum(abs(neg)) / abs(sum(diag(QH_cf)))),
  " of the signed trace)")
say("  -> eigenvalue truncation and a per-node total clamp are DIFFERENT")
say("     operators moving different mass; the latent layer names which one")
say("     and the receipt records the moved mass (contract sec.6 / G7).")

# ---------------------------------------------------------------------------
# P6.f  the CROSS-FITTED subject Gram: one N x N matrix behind P5, P6 and P7
# ---------------------------------------------------------------------------
rule("P6.f  the cross-fitted subject Gram unifies the trick and the split")

# Gamma_hat_{ii'} = sym cross-PARTITION inner product.  Its expectation is the
# Gram of the NOISELESS subject vectors -- on the diagonal too, because the two
# partitions are independent.  This is the object P5's covariance trick must be
# built on; P5's plug-in Gram carries the +62.7% inflation P6.b measures.
gram_cf <- function(Y1, Y2) (tcrossprod(Y1, Y2) + tcrossprod(Y2, Y1)) / 2

set.seed(31415)
Uload <- matrix(rnorm(Nsub * 2L), Nsub, 2L)
U <- Uload %*% t(Lb)
Ytrue <- sweep(U, 2L, mu_true, "+")
Y1 <- Ytrue + matrix(rnorm(Nsub * Dp, sd = sd_eps), Nsub, Dp)
Y2 <- Ytrue + matrix(rnorm(Nsub * Dp, sd = sd_eps), Nsub, Dp)
G_cf <- gram_cf(Y1, Y2)
G_true <- tcrossprod(Ytrue)
Ybar_i <- (Y1 + Y2) / 2
G_plug <- tcrossprod(Ybar_i)

say("  unbiasedness (2000 reps, Frobenius error against the true Gram):")
mc <- replicate(2000L, {
  Ul <- matrix(rnorm(Nsub * 2L), Nsub, 2L)
  Yt <- sweep(Ul %*% t(Lb), 2L, mu_true, "+")
  A <- Yt + matrix(rnorm(Nsub * Dp, sd = sd_eps), Nsub, Dp)
  B <- Yt + matrix(rnorm(Nsub * Dp, sd = sd_eps), Nsub, Dp)
  Gt <- tcrossprod(Yt)
  c(cf = mean(diag(gram_cf(A, B) - Gt)),
    plug = mean(diag(tcrossprod((A + B) / 2) - Gt)))
})
say("    mean diagonal bias, cross-fitted = ", sprintf("%+.6f", mean(mc["cf", ])),
  "  (MC se ", sprintf("%.4f", sd(mc["cf", ]) / sqrt(2000)), ")")
say("    mean diagonal bias, plug-in      = ", sprintf("%+.6f", mean(mc["plug", ])),
  "  (predicted +", sprintf("%.4f", Dp * sd_eps^2 / 2), ")")

# The scalar split of P6/P7 IS the diagonal-vs-off-diagonal contrast of Gamma.
off <- function(G) (sum(G) - sum(diag(G))) / (Nsub * (Nsub - 1L))
VW <- mean(rowSums(Y1 * Y2))
say("  mean(diag Gamma_cf) equals V^W of the population contract sec.7.1: ",
  sci(abs(mean(diag(G_cf)) - VW)))
say("  mean(offdiag Gamma_cf) = ", sprintf("%+.6f", off(G_cf)),
  " ; ||mu||^2 = ", sprintf("%+.6f", sum(mu_true^2)))
say("  consensus share from ONE N x N matrix = ",
  sprintf("%+.6f", off(G_cf) / mean(diag(G_cf))))

# Node-additive, so it still streams.
say("  node-additivity: max |Gamma_cf - sum_j Gamma_cf,j| = ", sci({
  acc <- matrix(0, Nsub, Nsub)
  for (j in seq_len(3L)) {
    cols <- ((j - 1L) * 7L + 1L):(j * 7L)
    acc <- acc + gram_cf(Y1[, cols, drop = FALSE], Y2[, cols, drop = FALSE])
  }
  max(abs(G_cf - acc))
}))

# Centered, it estimates the heterogeneity spectrum -- and is indefinite.
cen <- diag(Nsub) - matrix(1 / Nsub, Nsub, Nsub)
Hc <- cen %*% G_cf %*% cen / (Nsub - 1L)
Hp <- cen %*% G_plug %*% cen / (Nsub - 1L)
ev_c <- eigen(Hc, symmetric = TRUE)
ev_p <- eigen(Hp, symmetric = TRUE, only.values = TRUE)$values
say("  centered cross-fitted Gram spectrum (top 4) = ",
  paste(sprintf("%+.4f", ev_c$values[1:4]), collapse = " "))
say("  centered plug-in      Gram spectrum (top 4) = ",
  paste(sprintf("%+.4f", ev_p[1:4]), collapse = " "))
say("    trace, THIS DRAW: cross-fitted ", sprintf("%+.4f", sum(diag(Hc))),
  " vs plug-in ", sprintf("%+.4f", sum(diag(Hp))),
  " ; difference ", sprintf("%+.4f", sum(diag(Hp)) - sum(diag(Hc))),
  " (predicted noise term ", sprintf("%.4f", Dp * sd_eps^2 / 2), ")")
say("    a single draw of either trace is NOT an estimate of tr(Sigma_B):")
mc2 <- replicate(2000L, {
  Ul <- matrix(rnorm(Nsub * 2L), Nsub, 2L)
  Yt <- sweep(Ul %*% t(Lb), 2L, mu_true, "+")
  A <- Yt + matrix(rnorm(Nsub * Dp, sd = sd_eps), Nsub, Dp)
  B <- Yt + matrix(rnorm(Nsub * Dp, sd = sd_eps), Nsub, Dp)
  c(cf = sum(diag(cen %*% gram_cf(A, B) %*% cen)) / (Nsub - 1L),
    plug = sum(diag(cen %*% tcrossprod((A + B) / 2) %*% cen)) / (Nsub - 1L))
})
say("      over 2000 reps: cross-fitted mean ", sprintf("%+.4f", mean(mc2["cf", ])),
  " (sd ", sprintf("%.4f", sd(mc2["cf", ])), ") vs planted ",
  sprintf("%.4f", trB))
say("      over 2000 reps: plug-in      mean ", sprintf("%+.4f", mean(mc2["plug", ])),
  " (sd ", sprintf("%.4f", sd(mc2["plug", ])), ")")
say("    negative eigenvalues: cross-fitted ", sum(ev_c$values < 0),
  " of ", Nsub, " ; plug-in ", sum(ev_p < 0),
  "   -> the cross-fitted Gram is indefinite BY CONSTRUCTION")
plane2 <- Uload
say("  R^2 of cross-fitted Gram mode 1 on the planted loading plane = ",
  sprintf("%.6f", summary(lm(ev_c$vectors[, 1L] ~ plane2))$r.squared))
say("  R^2 of cross-fitted Gram mode 2 on the planted loading plane = ",
  sprintf("%.6f", summary(lm(ev_c$vectors[, 2L] ~ plane2))$r.squared))

say("\nDONE")
