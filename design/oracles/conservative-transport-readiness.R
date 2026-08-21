#!/usr/bin/env Rscript
# Oracle: what a population layer needs from a conservative result
# ================================================================
#
# Supports claims 6 and 7 of `design/conservative-geometry-contract.md`.
# Nothing here proposes an API; it establishes which properties WS-E may
# rely on and which it must supply for itself.
#
#   O3.a  The `symmetric_packed` codec is Frobenius consistent:
#           <svec(A), svec(B)> = <A, B>_F.
#         This is what licenses the subject-Gram eigen trick -- Euclidean
#         geometry in packed coordinates IS Frobenius geometry on forms.
#   O3.b  Geometry-space covariance over N subjects has rank <= N - 1, and
#         the N x N subject Gram reproduces its nonzero spectrum exactly.
#   O3.c  Row-stochastic transport P preserves each subject's budget exactly;
#         WITHOUT a sink node, a non-surjective P silently loses mass.
#   O3.d  Budget (sum) and density (mean) semantics give different group
#         numbers from the same P; neither is a default.
#   O3.e  Signed layer: crossvalidated contributions go negative, so shares,
#         cumulative curves and n_eff are undefined on the raw layer.
#         PSD clipping changes the total; it must be declared, never silent.
#   O3.f  What a conservative result carries today, and what it does not.
#
# Run:  Rscript design/oracles/conservative-transport-readiness.R

suppressMessages(pkgload::load_all(quiet = TRUE))
options(digits = 17)

say <- function(...) cat(..., "\n", sep = "")
rule <- function(title) say("\n", title, "\n", strrep("-", nchar(title)))
sci <- function(x) format(x, scientific = TRUE, digits = 4)

set.seed(20260817)

# ---------------------------------------------------------------------------
# O3.a  Frobenius consistency of the packed codec
# ---------------------------------------------------------------------------
rule("O3.a  symmetric_packed is Frobenius consistent")
q <- 5L
sym <- function() {
  m <- matrix(rnorm(q * q), q, q)
  0.5 * (m + t(m))
}
worst <- 0
for (i in seq_len(200L)) {
  A <- sym(); B <- sym()
  lhs <- sum(crossform:::.svec_symmetric(A) * crossform:::.svec_symmetric(B))
  rhs <- sum(A * B)
  worst <- max(worst, abs(lhs - rhs) / max(1, abs(rhs)))
}
say("  worst relative |<svec A, svec B> - <A,B>_F| over 200 draws = ", sci(worst))
say("  packed width for q=", q, " is ", q * (q + 1L) / 2L,
  " (vs ", q * q, " rectangular)")
say("  -> a packed row is an isometric embedding, so Euclidean distance,")
say("     inner product, covariance and PCA in packed coordinates are the")
say("     Frobenius versions of the same operations on the forms.")

# ---------------------------------------------------------------------------
# O3.b  Subject-Gram trick
# ---------------------------------------------------------------------------
rule("O3.b  geometry-space covariance has rank <= N - 1; subject Gram recovers it")
N <- 6L
P <- q * (q + 1L) / 2L
Y <- t(vapply(seq_len(N), function(s) crossform:::.svec_symmetric(sym()),
  numeric(P)))
centered <- sweep(Y, 2L, colMeans(Y))
Sigma <- crossprod(centered) / (N - 1L)              # P x P
Gram <- tcrossprod(centered) / (N - 1L)              # N x N
ev_sigma <- sort(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values,
  decreasing = TRUE)
ev_gram <- sort(eigen(Gram, symmetric = TRUE, only.values = TRUE)$values,
  decreasing = TRUE)
say("  N = ", N, "  P = ", P)
say("  rank(Sigma)                     = ", qr(Sigma)$rank, "  (bound: N-1 = ",
  N - 1L, ")")
say("  top ", N - 1L, " eigenvalues agree, max abs diff = ",
  sci(max(abs(ev_sigma[seq_len(N - 1L)] - ev_gram[seq_len(N - 1L)]))))
say("  Sigma eigenvalues ", N, "..", P, " max abs = ",
  sci(max(abs(ev_sigma[N:P]))))
say("  cost: N x N eigen instead of P x P; for a real study P = q(q+1)/2 with")
say("  q in the tens or hundreds, so this is the difference between feasible")
say("  and not.")

# ---------------------------------------------------------------------------
# O3.c  Row-stochastic transport, and why the sink is required
# ---------------------------------------------------------------------------
rule("O3.c  row-stochastic transport preserves the subject budget")
n_native <- 10L
n_group <- 4L
assign_to <- c(1, 1, 2, 2, 2, 3, 3, 4, 4, 4)
P_full <- matrix(0, n_native, n_group)
P_full[cbind(seq_len(n_native), assign_to)] <- 1
say("  P row sums: all one? ", all(abs(rowSums(P_full) - 1) < 1e-15))

# Per-node contributions of one subject: a conservative frame's total field.
contrib <- rnorm(n_native)
budget <- sum(contrib)
group_budget <- as.numeric(crossprod(P_full, contrib))
say("  subject budget            = ", sprintf("%.15f", budget))
say("  sum of group node budgets = ", sprintf("%.15f", sum(group_budget)))
say("  |difference|              = ", sci(abs(sum(group_budget) - budget)))

# Now drop the last group node's mass -- a partial atlas with no sink.
P_partial <- P_full
P_partial[assign_to == 4L, ] <- 0                    # rows no longer stochastic
say("\n  partial atlas (nodes 8-10 unmapped, NO sink):")
say("    P row sums in {0,1}: ", all(rowSums(P_partial) %in% c(0, 1)))
lost <- budget - sum(crossprod(P_partial, contrib))
say("    mass silently lost        = ", sprintf("%+.15f", lost),
  "   (", sprintf("%.1f%%", 100 * abs(lost / budget)), " of the budget)")

P_sink <- cbind(P_partial, sink = 0)
P_sink[rowSums(P_partial) == 0, n_group + 1L] <- 1
say("  with an explicit sink column:")
say("    P row sums: all one? ", all(abs(rowSums(P_sink) - 1) < 1e-15))
say("    |sum over group+sink - budget| = ",
  sci(abs(sum(crossprod(P_sink, contrib)) - budget)))
say("    sink holds ", sprintf("%+.6f", crossprod(P_sink, contrib)[n_group + 1L]),
  " -- visible, not lost.")

# ---------------------------------------------------------------------------
# O3.d  Budget vs density
# ---------------------------------------------------------------------------
rule("O3.d  budget (sum) and density (mean) are different group quantities")
subjects <- list(
  fine   = list(contrib = rep(1, 10) / 10, map = assign_to),
  coarse = list(contrib = rep(1, 4) / 4, map = c(1, 2, 3, 4))
)
say(sprintf("  %-8s %-10s %s", "subject", "budget", "per-group-node budget"))
for (nm in names(subjects)) {
  s <- subjects[[nm]]
  Pm <- matrix(0, length(s$contrib), n_group)
  Pm[cbind(seq_along(s$contrib), s$map)] <- 1
  b <- as.numeric(crossprod(Pm, s$contrib))
  d <- b / pmax(1, colSums(Pm))
  say(sprintf("  %-8s %-10.6f %s", nm, sum(s$contrib),
    paste(sprintf("%.4f", b), collapse = " ")))
  say(sprintf("  %-8s %-10s %s", "", "(density)",
    paste(sprintf("%.4f", d), collapse = " ")))
}
say("  Both subjects carry budget 1.  Under BUDGET semantics the fine subject")
say("  puts 0.3 into group node 2 and the coarse one 0.25, purely because its")
say("  native frame is finer.  Under DENSITY semantics both give 0.1 and 0.25.")
say("  Neither is a default; the choice is a declared field on the transport.")

# ---------------------------------------------------------------------------
# O3.e  Signed layer, and what PSD clipping costs
# ---------------------------------------------------------------------------
rule("O3.e  crossvalidated contributions are signed")
n_features <- 8L
qq <- 3L
domain <- abstract_domain(n_features,
  coordinates = cbind(seq_len(n_features) - 1, 0), id = "oracle:transport:v1")
# Pure noise: the crossvalidated estimator is unbiased for zero, so roughly
# half the per-node values must be negative.
mk <- function() {
  m <- matrix(rnorm(qq * n_features), qq, n_features)
  rownames(m) <- c("a", "b", "c")
  m
}
rel <- relation(list(run1 = mk(), run2 = mk()), domain = domain)
over <- cross_partitions(rel)
frame <- compile_frame(searchlights(1.01, "conservative"), domain)
plan <- plan_geometry(rel, frame, over)
view <- contrast_energy(plan, c(a = 1, b = -1, c = 0))
say("  per-node `total` under a pure-noise fixture:")
say("    ", paste(sprintf("%+.4f", view$total), collapse = "  "))
say("  negative nodes: ", sum(view$total < 0), " of ", length(view$total))
say("  sum (the conserved budget) = ", sprintf("%+.9f", sum(view$total)))
say("  coherence_fraction valid at ", sum(view$coherence_fraction_valid),
  " of ", length(view$coherence_fraction), " nodes; NA elsewhere:")
say("    ", paste(ifelse(is.na(view$coherence_fraction), "NA",
  sprintf("%.3f", view$coherence_fraction)), collapse = "  "))
say("  -> a contribution SHARE (x / sum) on this layer is meaningless: the")
say("     denominator is a signed budget and shares need not lie in [0,1].")
shares <- view$total / sum(view$total)
say("    share range = [", sprintf("%.3f", min(shares)), ", ",
  sprintf("%.3f", max(shares)), "]   any outside [0,1]: ",
  any(shares < 0 | shares > 1))

say("\n  cost of clipping to a latent PSD layer:")
clipped <- pmax(view$total, 0)
say("    signed total  = ", sprintf("%+.9f", sum(view$total)))
say("    clipped total = ", sprintf("%+.9f", sum(clipped)))
say("    inflation     = ", sprintf("%+.2f%%",
  100 * (sum(clipped) - sum(view$total)) / abs(sum(view$total))))
say("    -> clipping reintroduces exactly the noise bias the cross-partition")
say("       pairing removes, and it breaks conservation.  A latent layer is a")
say("       separate, declared object, not a display option.")

# ---------------------------------------------------------------------------
# O3.f  What the result actually carries today
# ---------------------------------------------------------------------------
rule("O3.f  metadata a conservative result exposes today")
labelled <- abstract_domain(n_features,
  coordinates = cbind(seq_len(n_features) - 1, 0),
  feature_ids = paste0("vox", seq_len(n_features)),
  id = "oracle:transport:labelled:v1")
rel_l <- relation(list(run1 = mk(), run2 = mk()), domain = labelled)
frame_l <- compile_frame(searchlights(1.01, "conservative"), labelled)
view_l <- contrast_energy(
  plan_geometry(rel_l, frame_l, cross_partitions(rel_l)),
  c(a = 1, b = -1, c = 0)
)
say("  frame$index columns        : ",
  paste(names(frame_l$index), collapse = ", "))
say("  frame$specification fields : ",
  paste(names(frame_l$specification), collapse = ", "))
say("  frame$specification$radius : ", frame_l$specification$radius,
  "   (frame-wide, NOT per row)")
say("  frame$normalization        : ", frame_l$normalization)
say("  as.data.frame(view) columns: ",
  paste(names(as.data.frame(view_l)), collapse = ", "))
say("  first measurement labels   : ",
  paste(utils::head(as.data.frame(view_l)$measurement, 3), collapse = ", "))
say("  -> node LABELS do reach the result when the frame was compiled.")
say("     What is absent is any per-row family / scale / center column: the")
say("     radius lives once on $specification, so a stacked family cannot say")
say("     which scale a row belongs to.")

family <- additive_frame(
  rbind(0.5 * as.matrix(compile_frame(voxelwise("conservative"), labelled)$weights),
        0.5 * as.matrix(frame_l$weights)),
  normalization = "conservative", domain = labelled
)
say("\n  a two-scale family built by rbind(...) then additive_frame():")
say("    rows                 = ", nrow(family$weights),
  "   conserved: ", frame_conservation(family)$conserved)
say("    has $index ?         ", !is.null(family$index))
say("    has $specification ? ", !is.null(family$specification))
family_view <- contrast_energy(
  plan_geometry(rel_l, family, cross_partitions(rel_l)),
  c(a = 1, b = -1, c = 0)
)
say("    measurement column   : ",
  paste(utils::head(as.data.frame(family_view)$measurement, 4), collapse = ", "),
  " ... (positions, not labels)")
say("  -> the construction works numerically and is provenance blind: the")
say("     family route loses even the node label it would have had, and no")
say("     row can say which scale, center or family it came from.  That is")
say("     the single largest gap between today's frame and a transportable")
say("     one.")

say("\nDONE")
