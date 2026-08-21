#!/usr/bin/env Rscript
# Oracle: the conservative multiscale ledger
# ==========================================
#
# Supports claims 2 and 3 of `design/conservative-geometry-contract.md`.
#
#   O1.a  A conservative (column-normalized) frame conserves the `total`
#         component exactly:  sum_x G_x = G_Omega.
#   O1.b  The `total` component of a conservative frame is a *linear
#         reallocation of one fixed per-voxel ledger*:
#            G_x^total = sum_v w_xv b_v^(L) b_v^(R)T
#         reproduced here from first principles without calling any crossform
#         contraction.
#   O1.c  Consequently, for an alpha-weighted frame family stacking S scales,
#            sum_{x in scale s} G_{s,x} = alpha_s * G_Omega
#         for EVERY scale, to machine precision. Per-scale total energy is
#         therefore fixed by the analyst's choice of alpha and carries no
#         multivoxel information.
#   O1.d  The coherent component does NOT conserve, and the coherent share
#         *does* vary with scale. That is the informative object.
#
# Run:  Rscript design/oracles/conservative-multiscale-ledger.R

suppressMessages(pkgload::load_all(quiet = TRUE))
options(digits = 17)

say <- function(...) cat(..., "\n", sep = "")
rule <- function(title) say("\n", title, "\n", strrep("-", nchar(title)))

# ---------------------------------------------------------------------------
# Fixture: a 1-D line domain so that Euclidean searchlights are intervals.
# ---------------------------------------------------------------------------
set.seed(20260817)
n_features <- 12L
q <- 3L
domain <- abstract_domain(
  n_features,
  coordinates = cbind(seq_len(n_features) - 1, 0),
  id = "oracle:multiscale:v1"
)

effect_names <- c("a", "b", "c")
mk <- function() {
  m <- matrix(rnorm(q * n_features), q, n_features)
  rownames(m) <- effect_names
  m
}
effects <- list(run1 = mk(), run2 = mk())
rel <- relation(effects, domain = domain)
over <- cross_partitions(rel)

# The global comparator must be the *unnormalized* whole-brain operator:
# `whole_brain()` defaults to "local", which divides by the feature count.
global <- materialize_geometry(
  rel, compile_frame(whole_brain("none"), domain), over
)
G_omega <- drop(geometry_component(global, "total"))
coh_omega <- drop(geometry_component(global, "coherent"))

# ---------------------------------------------------------------------------
# O1.b  The per-voxel ledger, from first principles.
#
# `cross_partitions()` on two partitions gives the symmetrized estimator
#   G = 1/2 ( B_1 K B_2^T + B_2 K B_1^T ),
# so the per-voxel atom is the symmetrized outer product of the two runs'
# columns.  We build the ledger in the packed codec crossform uses
# (lower triangle, column major, sqrt(2) off-diagonals) so it can be compared
# to `geometry_component()` output directly.
# ---------------------------------------------------------------------------
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

ledger <- t(vapply(seq_len(n_features), function(v) {
  b1 <- effects$run1[, v]
  b2 <- effects$run2[, v]
  svec(0.5 * (tcrossprod(b1, b2) + tcrossprod(b2, b1)))
}, numeric(q * (q + 1) / 2)))
stopifnot(dim(ledger) == c(n_features, q * (q + 1) / 2))

rule("O1.b  total is the frame contraction of one fixed per-voxel ledger")
ledger_global <- colSums(ledger)
say("max |ledger colSums - G_Omega|            = ",
  format(max(abs(ledger_global - G_omega)), scientific = TRUE))

# ---------------------------------------------------------------------------
# O1.a / O1.c  The alpha-weighted frame family.
# ---------------------------------------------------------------------------
scales <- list(
  point       = voxelwise("conservative"),
  searchlight = searchlights(1.01, "conservative"),
  wide        = searchlights(2.01, "conservative")
)
alpha <- c(point = 0.2, searchlight = 0.5, wide = 0.3)
stopifnot(abs(sum(alpha) - 1) < 1e-15)

frames <- lapply(scales, compile_frame, domain = domain)
blocks <- vapply(frames, function(f) nrow(f$weights), integer(1))

# Each individual scale conserves on its own (O1.a).
rule("O1.a  each conservative scale conserves `total` on its own")
for (s in names(frames)) {
  g <- materialize_geometry(rel, frames[[s]], over)
  W <- as.matrix(frames[[s]]$weights)
  tot <- geometry_component(g, "total")
  say(sprintf("  %-12s nodes=%2d  max|colSums W - 1| = %s", s,
    nrow(W), format(max(abs(colSums(W) - 1)), scientific = TRUE)))
  say(sprintf("               max|sum_x G_x - G_Omega|      = %s",
    format(max(abs(colSums(tot) - G_omega)), scientific = TRUE)))
  # O1.b again, per scale: the frame contraction of the ledger IS the total.
  say(sprintf("               max|G_x - (W %%*%% ledger)_x|   = %s",
    format(max(abs(tot - W %*% ledger)), scientific = TRUE)))
}

# Stack the scales with alpha weights into a single conservative family frame.
family_weights <- do.call(rbind, lapply(names(frames), function(s) {
  alpha[[s]] * as.matrix(frames[[s]]$weights)
}))
family_frame <- additive_frame(
  family_weights, normalization = "conservative", domain = domain
)
say("\nfamily frame: ", nrow(family_weights), " rows, ",
  "max|colSums - 1| = ",
  format(max(abs(colSums(family_weights) - 1)), scientific = TRUE))

family_geom <- materialize_geometry(rel, family_frame, over)
tot_fam <- geometry_component(family_geom, "total")
coh_fam <- geometry_component(family_geom, "coherent")

row_of <- function(s) {
  offset <- if (s == names(blocks)[1]) 0L else sum(blocks[seq_len(match(s, names(blocks)) - 1L)])
  offset + seq_len(blocks[[s]])
}

rule("O1.c  per-scale totals equal alpha_s * G_Omega (by construction)")
worst_abs <- 0
worst_rel <- 0
for (s in names(blocks)) {
  rows <- row_of(s)
  observed <- colSums(tot_fam[rows, , drop = FALSE])
  expected <- alpha[[s]] * G_omega
  abs_err <- max(abs(observed - expected))
  rel_err <- abs_err / max(abs(expected))
  worst_abs <- max(worst_abs, abs_err)
  worst_rel <- max(worst_rel, rel_err)
  say(sprintf("  scale %-12s alpha=%.2f  max abs err = %s   rel = %s",
    s, alpha[[s]], format(abs_err, scientific = TRUE),
    format(rel_err, scientific = TRUE)))
}
say(sprintf("  worst absolute = %s  |  worst relative = %s",
  format(worst_abs, scientific = TRUE), format(worst_rel, scientific = TRUE)))
say("  PASS (<= 1e-12): ", worst_abs <= 1e-12)

# The whole family conserves too: sum over ALL rows == G_Omega, since
# sum_s alpha_s = 1.
say(sprintf("  whole family: max|sum_all G_x - G_Omega| = %s",
  format(max(abs(colSums(tot_fam) - G_omega)), scientific = TRUE)))

# ---------------------------------------------------------------------------
# The scalar readout a Panel-D style figure would plot.
# ---------------------------------------------------------------------------
H <- diag(q)                                   # any fixed query
Hp <- svec(H)
rule("O1.c  a fixed query <H, .> reproduces the same statement")
say(sprintf("  <H, G_Omega> = %.12f", sum(Hp * G_omega)))
for (s in names(blocks)) {
  rows <- row_of(s)
  energy <- sum(tot_fam[rows, , drop = FALSE] %*% Hp)
  say(sprintf("  scale %-12s  E_s = %.12f   alpha_s * <H,G_Omega> = %.12f",
    s, energy, alpha[[s]] * sum(Hp * G_omega)))
}

# ---------------------------------------------------------------------------
# O1.d  Coherent does not conserve; the coherent SHARE varies with scale.
# ---------------------------------------------------------------------------
rule("O1.d  coherent does not conserve; the coherent share is scale dependent")
say(sprintf("  <H, coherent_Omega>       = %.12f", sum(Hp * coh_omega)))
for (s in names(blocks)) {
  rows <- row_of(s)
  coh <- sum(coh_fam[rows, , drop = FALSE] %*% Hp)
  tot <- sum(tot_fam[rows, , drop = FALSE] %*% Hp)
  say(sprintf(
    "  scale %-12s  E_coh = %+.9f  (alpha_s*<H,coh_Omega> = %+.9f)  share = %+.6f",
    s, coh, alpha[[s]] * sum(Hp * coh_omega), coh / tot))
}
say("\n  Interpretation: the E_s column above is alpha_s times a constant.")
say("  The `share` column is not, and it is the only scale-resolved quantity")
say("  in this table that depends on the data rather than on alpha.")

# Single-scale coherent totals versus the global coherent component.
rule("O1.d  per-scale coherent sums vs the global coherent component")
for (s in names(frames)) {
  g <- materialize_geometry(rel, frames[[s]], over)
  coh <- colSums(geometry_component(g, "coherent"))
  say(sprintf("  %-12s max|sum_x G^coh_x - G^coh_Omega| = %s",
    s, format(max(abs(coh - coh_omega)), scientific = TRUE)))
}
say("\n(voxelwise is the degenerate case: every node is a singleton, so its")
say(" coherent component equals its total and the deviation is the whole")
say(" global configuration mass.)")

say("\nDONE")
