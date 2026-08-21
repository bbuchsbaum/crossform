#!/usr/bin/env Rscript
# Oracle: conservation under a fixed neural metric
# ================================================
#
# Supports claim 5 of `design/conservative-geometry-contract.md`.
#
#   O2.a  Identity metric: sum_x G_x = G_Omega (machine precision).
#   O2.b  Diagonal metric Q = D(q): conservation survives, because the
#         native composition D(sqrt(w)) Q D(sqrt(w)) = D(w q) is still
#         feature additive and sum_x w_xv = 1.
#   O2.c  Dense metric Q: conservation FAILS under the native composition.
#         The exact algebraic law for the failure is
#            sum_x G_x = B_L (S o Q) B_R^T,   S_uv = sum_x sqrt(w_xu w_xv),
#         with S_uu = 1 but S_uv != 1 off the diagonal.  The oracle checks the
#         prediction against the executed package result.
#   O2.d  The alternative *whitened* composition Q^(1/2) D(w) Q^(1/2)
#         conserves exactly, because sum_x Q^(1/2) D(w_x) Q^(1/2)
#         = Q^(1/2) D(1) Q^(1/2) = Q.  It is a different estimand: the
#         per-node values differ from the native composition.  It is computed
#         from first principles here and, since D6, also checked against the
#         executed package path `plan_geometry(..., composition = "whitened")`.
#   O2.e  `.metric_frame_conservation()` refuses a non-diagonal schedule.
#   O2.f  The diagonal-metric route rebuilds the frame with
#         `normalization = "none"`, discarding the declared normalization from
#         provenance even though the numbers survive.
#
# Run:  Rscript design/oracles/conservative-metric-composition.R

suppressMessages(pkgload::load_all(quiet = TRUE))
options(digits = 17)

say <- function(...) cat(..., "\n", sep = "")
rule <- function(title) say("\n", title, "\n", strrep("-", nchar(title)))
sci <- function(x) format(x, scientific = TRUE, digits = 4)

# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------
set.seed(20260817)
n_features <- 9L
q <- 3L
domain <- abstract_domain(
  n_features,
  coordinates = cbind(seq_len(n_features) - 1, 0),
  id = "oracle:metric-composition:v1"
)
mk <- function() {
  m <- matrix(rnorm(q * n_features), q, n_features)
  rownames(m) <- c("a", "b", "c")
  m
}
effects <- list(run1 = mk(), run2 = mk())
rel <- relation(effects, domain = domain)
over <- cross_partitions(rel)
B1 <- effects$run1
B2 <- effects$run2

conservative <- compile_frame(searchlights(1.01, "conservative"), domain)
W <- as.matrix(conservative$weights)
global_frame <- compile_frame(whole_brain("none"), domain)

svec <- function(M) {
  M <- 0.5 * (M + t(M))
  out <- numeric(q * (q + 1) / 2)
  k <- 1L
  for (j in seq_len(q)) for (i in j:q) {
    out[k] <- if (i == j) M[i, j] else sqrt(2) * M[i, j]
    k <- k + 1L
  }
  out
}
# The cross-partition estimator with two runs: sym(B1 K B2^T).
geom_K <- function(K) svec(0.5 * (B1 %*% K %*% t(B2) + B2 %*% K %*% t(B1)))

run_total <- function(frame, metric = NULL, composition = "native") {
  plan <- plan_geometry(rel, frame, over, metric = metric,
    composition = composition)
  geometry_component(materialize_geometry(plan), "total")
}

report <- function(label, local_total, global_total) {
  observed <- colSums(local_total)
  abs_err <- max(abs(observed - global_total))
  rel_err <- max(abs(observed - global_total)) / max(abs(global_total))
  # A signed scalar readout of the same statement, so the sign of the failure
  # is visible: trace of the summed local form vs the global form.
  tr_local <- sum(observed[c(1L, 4L, 6L)])
  tr_global <- sum(global_total[c(1L, 4L, 6L)])
  say(sprintf("  %-34s max abs = %-11s  max rel = %-11s  signed trace err = %+.4f%%",
    label, sci(abs_err), sci(rel_err), 100 * (tr_local - tr_global) / tr_global))
  invisible(list(abs = abs_err, rel = rel_err,
    signed = (tr_local - tr_global) / tr_global))
}

# ---------------------------------------------------------------------------
# O2.a  Identity metric
# ---------------------------------------------------------------------------
rule("O2.a  identity metric (implicit)")
G_omega_I <- drop(run_total(global_frame))
res_identity <- report("implicit identity", run_total(conservative), G_omega_I)
say("  algebraic check |G_Omega - sym(B1 B2^T)| = ",
  sci(max(abs(G_omega_I - geom_K(diag(n_features))))))

# ---------------------------------------------------------------------------
# O2.b  Diagonal metric
# ---------------------------------------------------------------------------
rule("O2.b  diagonal metric Q = D(q), q ~ U(0.5, 2.5)")
qdiag <- runif(n_features, 0.5, 2.5)
Qd <- diag(qdiag)
metric_d <- neural_metric(Qd, domain)
say("  native_diagonal capability  = ", metric_d$capabilities$native_diagonal)
say("  feature_additive capability = ", metric_d$capabilities$feature_additive)
G_omega_d <- drop(run_total(global_frame, metric_d))
res_diag <- report("diagonal metric", run_total(conservative, metric_d), G_omega_d)
say("  algebraic check |G_Omega(Qd) - sym(B1 Qd B2^T)| = ",
  sci(max(abs(G_omega_d - geom_K(Qd)))))

# ---------------------------------------------------------------------------
# O2.c  Dense metric -- native composition D(sqrt(w)) Q D(sqrt(w))
# ---------------------------------------------------------------------------
rule("O2.c  dense metric Q = crossprod(A)/9 + I  (native composition)")
A <- matrix(rnorm(n_features * n_features), n_features, n_features)
Q <- crossprod(A) / n_features + diag(n_features)
metric_Q <- neural_metric(Q, domain)
say("  native_diagonal capability  = ", metric_Q$capabilities$native_diagonal)
say("  feature_additive capability = ", metric_Q$capabilities$feature_additive)
say("  condition number            = ", sci(kappa(Q, exact = TRUE)))
G_omega_Q <- drop(run_total(global_frame, metric_Q))
local_Q <- run_total(conservative, metric_Q)
res_dense <- report("dense metric, native", local_Q, G_omega_Q)

# The exact algebraic prediction for the failure.
S <- matrix(0, n_features, n_features)
for (x in seq_len(nrow(W))) {
  rw <- sqrt(W[x, ])
  S <- S + tcrossprod(rw)
}
predicted <- geom_K(S * Q)
say("  |sum_x G_x - B(S o Q)B^T| (prediction)   = ",
  sci(max(abs(colSums(local_Q) - predicted))))
say("  max |S_uu - 1| (diagonal of S)           = ", sci(max(abs(diag(S) - 1))))
say("  range of off-diagonal S_uv               = [",
  sprintf("%.4f", min(S[upper.tri(S)])), ", ",
  sprintf("%.4f", max(S[upper.tri(S)])), "]")

# ---------------------------------------------------------------------------
# O2.d  Whitened composition Q^(1/2) D(w) Q^(1/2)  (no package path today)
# ---------------------------------------------------------------------------
rule("O2.d  whitened composition Q^(1/2) D(w) Q^(1/2)  (first principles)")
eig <- eigen(Q, symmetric = TRUE)
Qhalf <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
say("  |Qhalf %*% Qhalf - Q| = ", sci(max(abs(Qhalf %*% Qhalf - Q))))

whitened <- t(vapply(seq_len(nrow(W)), function(x) {
  geom_K(Qhalf %*% diag(W[x, ]) %*% Qhalf)
}, numeric(q * (q + 1) / 2)))
native <- t(vapply(seq_len(nrow(W)), function(x) {
  rw <- sqrt(W[x, ])
  geom_K(Q * tcrossprod(rw))
}, numeric(q * (q + 1) / 2)))

say("  first-principles native reproduces package: ",
  sci(max(abs(native - local_Q))))
# Since D6 the package has a whitened code path, so the first-principles values
# above are no longer the only evidence: they are now a check on it.
local_whitened <- run_total(conservative, metric_Q, composition = "whitened")
say("  first-principles whitened reproduces package: ",
  sci(max(abs(whitened - local_whitened))))
res_whitened <- report("dense metric, whitened", whitened, G_omega_Q)
res_whitened_pkg <- report("dense metric, whitened (package)",
  local_whitened, G_omega_Q)
say("  per-node disagreement max|whitened - native| = ",
  sci(max(abs(whitened - native))),
  "   (relative to max|native| = ", sci(max(abs(native))), ")")
say("  -> the two compositions are different estimands at the node level,")
say("     and only the whitened one closes the global ledger.")

# ---------------------------------------------------------------------------
# O2.d'  The ROOT is part of the whitened estimand.
#
# The conservation argument never uses symmetry of the root:
#     sum_x R D(w_x) R^T = R (sum_x D(w_x)) R^T = R I R^T = Q
# holds for ANY R with R R^T = Q. So conservation is root invariant. The
# per-node values are NOT. Writing the law as Q^(1/2) D(w) Q^(1/2) silently
# selects the symmetric PSD root; "whitened" alone does not name an estimand,
# and a conservation certificate cannot tell two roots apart.
# Supports claim 5g of the contract (added by the 2026-08-20 review).
# ---------------------------------------------------------------------------
rule("O2.d'  the root is part of the estimand: symmetric vs Cholesky")
Qchol <- t(chol(Q))                        # lower factor, Qchol %*% t(Qchol) = Q
say("  |Qsym  Qsym^T  - Q| = ", sci(max(abs(Qhalf %*% t(Qhalf) - Q))),
  "   (symmetric PSD root)")
say("  |Qchol Qchol^T - Q| = ", sci(max(abs(Qchol %*% t(Qchol) - Q))),
  "   (lower Cholesky factor)")

whitened_chol <- t(vapply(seq_len(nrow(W)), function(x) {
  geom_K(Qchol %*% diag(W[x, ]) %*% t(Qchol))
}, numeric(q * (q + 1) / 2)))

rel_sym <- max(abs(colSums(whitened) - G_omega_Q)) / max(abs(G_omega_Q))
rel_chol <- max(abs(colSums(whitened_chol) - G_omega_Q)) / max(abs(G_omega_Q))
say("  conservation, symmetric root : max rel |sum_x G_x - G_Omega| = ",
  sci(rel_sym))
say("  conservation, Cholesky root  : max rel |sum_x G_x - G_Omega| = ",
  sci(rel_chol))
node_gap <- max(abs(whitened - whitened_chol))
say("  per-node disagreement between roots = ", sci(node_gap),
  "   (relative to max|symmetric| = ", sci(max(abs(whitened))), ")")
say(sprintf("  -> both roots conserve to machine precision, and their node"))
say(sprintf("     values differ by %.1f%% of the largest node value.",
  100 * node_gap / max(abs(whitened))))
say("     `composition = \"whitened\"` must therefore name the ROOT.")
say("     The contract pins the symmetric PSD root (section 5.2.1); the root")
say("     identity has to enter plan identity, because the conservation")
say("     certificate is blind to this choice.")

# ---------------------------------------------------------------------------
# O2.c'  The SIZE and SIGN of the dense-metric failure are fixture dependent.
#
# `.planning/2026-08-17-feedback-assessment.md` reports -6.6% for one draw of
# `crossprod(A)/9 + I`; this oracle's default draw gives a different, positive
# number.  Both are instances of the same law.  The sweep below shows the
# signed error crossing zero across draws, so no single percentage is the
# contract claim -- the claim is "not zero, and not bounded".
# ---------------------------------------------------------------------------
rule("O2.c'  sweep over 12 draws of Q = crossprod(A)/p + I")
signed_errors <- vapply(seq_len(12L), function(draw) {
  set.seed(1000L + draw)
  Ad <- matrix(rnorm(n_features * n_features), n_features, n_features)
  Qd2 <- crossprod(Ad) / n_features + diag(n_features)
  Sd <- matrix(0, n_features, n_features)
  for (x in seq_len(nrow(W))) Sd <- Sd + tcrossprod(sqrt(W[x, ]))
  observed <- geom_K(Sd * Qd2)
  expected <- geom_K(Qd2)
  tr <- c(1L, 4L, 6L)
  (sum(observed[tr]) - sum(expected[tr])) / sum(expected[tr])
}, numeric(1))
say("  signed trace errors (%): ",
  paste(sprintf("%+.2f", 100 * signed_errors), collapse = "  "))
say(sprintf("  range = [%+.2f%%, %+.2f%%]   any negative: %s   any |.| < 1%%: %s",
  100 * min(signed_errors), 100 * max(signed_errors),
  any(signed_errors < 0), any(abs(signed_errors) < 0.01)))
say("  -> the assessment's -6.6% and this oracle's +21.2% are both draws from")
say("     this distribution.  The invariant claim is the algebraic law above,")
say("     not a percentage.")

# ---------------------------------------------------------------------------
# O2.e  The package refuses a non-diagonal schedule
# ---------------------------------------------------------------------------
rule("O2.e  .metric_frame_conservation() certificates")
node_metric <- function(frame, node, builder) {
  pos <- which(as.matrix(frame$weights)[node, ] > 0)
  support <- frame$domain$feature_ids[pos]
  neural_metric(builder(length(pos), pos), frame$domain, support = support)
}
cert_identity <- crossform:::.metric_frame_conservation(conservative)
say("  identity  : feature_additive=", cert_identity$feature_additive,
  "  identity_conservation=", cert_identity$identity_conservation)
say("              global_metric_kind=", cert_identity$global_metric_kind)
say("              reason: ", cert_identity$reason)

diag_metrics <- lapply(seq_len(nrow(W)), function(node) {
  node_metric(conservative, node, function(k, pos) diag(qdiag[pos], nrow = k))
})
cert_diag <- crossform:::.metric_frame_conservation(conservative, diag_metrics)
say("  diagonal  : feature_additive=", cert_diag$feature_additive,
  "  identity_conservation=", cert_diag$identity_conservation)
say("              reason: ", cert_diag$reason)

dense_metrics <- lapply(seq_len(nrow(W)), function(node) {
  node_metric(conservative, node, function(k, pos) Q[pos, pos, drop = FALSE])
})
cert_dense <- crossform:::.metric_frame_conservation(conservative, dense_metrics)
say("  dense     : feature_additive=", cert_dense$feature_additive,
  "  identity_conservation=", cert_dense$identity_conservation)
say("              global_metric_kind=", cert_dense$global_metric_kind)
say("              reason: ", cert_dense$reason)
say("  .require_metric_conservation(dense, 'feature_additive') -> ",
  tryCatch({
    crossform:::.require_metric_conservation(cert_dense, "feature_additive")
    "ACCEPTED (unexpected)"
  }, error = function(e) paste0("refused: ", conditionMessage(e))))

# ---------------------------------------------------------------------------
# O2.f  The diagonal-metric frame fold and its provenance
# ---------------------------------------------------------------------------
rule("O2.f  diagonal-metric frame fold: weights change, provenance is kept")
schedule <- crossform:::.geometry_metric_schedule(conservative, metric_d)
rebuilt <- crossform:::.metric_additive_frame(conservative, schedule)
say("  declared frame normalization = ", conservative$normalization)
say("  folded   frame normalization = ", rebuilt$normalization,
  "   (the truth about the weights it now carries)")
say("  folded has $index ?          ", !is.null(rebuilt$index),
  "   folded has $specification ? ", !is.null(rebuilt$specification))
say("  folded has $metric_folded ?  ", !is.null(rebuilt$metric_folded))
if (!is.null(rebuilt$metric_folded)) {
  say("    declared_normalization = ", rebuilt$metric_folded$declared_normalization)
  say("    composition            = ", rebuilt$metric_folded$composition)
  say("    metric_kind            = ", rebuilt$metric_folded$metric_kind)
  say("    |reference_mass - q|   = ",
    sci(max(abs(rebuilt$metric_folded$reference_mass - qdiag))))
}
say("  folded colSums range         = [",
  sprintf("%.6f", min(Matrix::colSums(rebuilt$weights))), ", ",
  sprintf("%.6f", max(Matrix::colSums(rebuilt$weights))), "]  (== q_v)")
say("  numbers survive: |folded colSums - q| = ",
  sci(max(abs(Matrix::colSums(rebuilt$weights) - qdiag))))
say("  -> the declared normalization is recoverable from provenance, so the")
say("     fold is legible; what a *family* row still cannot recover is which")
say("     scale it came from (see conservative-transport-readiness.R O3.f).")

# ---------------------------------------------------------------------------
rule("Summary")
say(sprintf("  identity            rel err = %-11s   signed trace err = %+.4f%%",
  sci(res_identity$rel), 100 * res_identity$signed))
say(sprintf("  diagonal            rel err = %-11s   signed trace err = %+.4f%%",
  sci(res_diag$rel), 100 * res_diag$signed))
say(sprintf("  dense (native)      rel err = %-11s   signed trace err = %+.4f%%",
  sci(res_dense$rel), 100 * res_dense$signed))
say(sprintf("  dense (whitened)    rel err = %-11s   signed trace err = %+.4f%%",
  sci(res_whitened$rel), 100 * res_whitened$signed))
say("\nDONE")
