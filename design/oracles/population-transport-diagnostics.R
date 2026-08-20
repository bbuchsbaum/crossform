#!/usr/bin/env Rscript
# Oracle: eta_transport, the transport diagnostics slice, and labelling
# =====================================================================
#
# Supports claims 7, 8 and 9 of `design/population-form-contract.md`.
#
#   P7.a  Cross-fitted consensus share: V^C from cross-participant products,
#         V^W from cross-partition products.  eta_transport is the change in
#         consensus share on HELD-OUT data.
#   P7.b  A functional transport estimated on an independent run gives a
#         positive eta; the SAME estimate re-used on the evaluation run gives
#         a much larger, circular eta.  Both measured.
#   P7.c  A null transport gives a NEGATIVE eta.  It is reported as-is.
#   P7.d  eta is not interpretable alone: a transport that sinks its noisy
#         periphery raises the consensus share without aligning anything.
#         This is why the diagnostics slice is required alongside eta.
#   P7.e  The diagnostics: displacement, entropy/perplexity, sink mass
#         (data-free coverage and data-dependent budget), group-node subject
#         coverage.  Each computed from its stated definition.
#   P8    Transported coherent is a LEDGER of native-node coherence, not the
#         coherence of the group node's own geometry.  On a fixture where both
#         are computable they differ by a large margin, while the transported
#         TOTAL matches the group node's own total exactly.
#   P9    Executable readiness predicates: what a conservative result must
#         supply, and what crossform refuses to supply.
#
# Pure matrix algebra; the crossform package is NOT loaded.
#
# Run:  Rscript design/oracles/population-transport-diagnostics.R

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

# ---------------------------------------------------------------------------
# Fixture: subjects whose functional bump is displaced from anatomy
# ---------------------------------------------------------------------------
rule("fixture: N subjects, a group node line, per-subject functional shift")

N <- 16L
m <- 12L                                   # group nodes at coordinates 1..m
qq <- 3L
ixq <- svec_index(qq)
Pd <- qq * (qq + 1L) / 2L                  # packed width 6
g_dir <- svec(crossprod(matrix(rnorm(qq * qq), qq, qq)) / qq, ixq)
g_dir <- g_dir / sqrt(sum(g_dir^2))
amp <- 1.5                                 # bump amplitude
tau <- 0.35                                # between-subject form heterogeneity
sd_noise <- 0.45                           # per-run measurement noise
bump <- function(x, c0) exp(-(x - c0)^2 / (2 * 1.5^2))

delta_true <- sample(-2:2, N, replace = TRUE)
u_subj <- lapply(seq_len(N), function(i) rnorm(Pd, sd = tau))

# Three independent runs: run 1 trains the transport, runs 2-3 are held out
# and supply the two partitions the cross-fitted estimators need.
draw_run <- function() {
  lapply(seq_len(N), function(i) {
    prof <- bump(seq_len(m), 6.5 + delta_true[i])
    signal <- outer(amp * prof, g_dir + u_subj[[i]])
    signal + matrix(rnorm(m * Pd, sd = sd_noise), m, Pd)
  })
}
R1 <- draw_run(); R2 <- draw_run(); R3 <- draw_run()
say("  N = ", N, " subjects, ", m, " native nodes each, packed width ", Pd)
say("  true shifts delta_i = ", paste(sprintf("%+d", delta_true), collapse = " "))
say("  between-subject form sd = ", tau, " ; per-run noise sd = ", sd_noise)

# ---------------------------------------------------------------------------
# P7.a  the transports, and the cross-fitted consensus share
# ---------------------------------------------------------------------------
rule("P7.a  cross-fitted consensus share and eta_transport")

shift_transport <- function(delta) {
  M <- matrix(0, m, m + 1L)
  for (x in seq_len(m)) {
    tgt <- x - delta
    if (tgt >= 1L && tgt <= m) M[x, tgt] <- 1 else M[x, m + 1L] <- 1
  }
  M
}
shift_profile <- function(v, d) {
  out <- rep(NA_real_, m)
  for (x in seq_len(m)) {
    tgt <- x - d
    if (tgt >= 1L && tgt <= m) out[tgt] <- v[x]
  }
  out
}
# Leave-one-out iterative cross-correlation alignment.  The reference is built
# from the OTHER subjects only, so the estimate of subject i never sees its own
# contribution -- the cross-fit discipline applies within the training run too.
est_shift <- function(runs, cand = -3:3, iters = 4L) {
  S <- t(vapply(runs, function(Z) rowSums(Z^2), numeric(m)))
  S <- sweep(S, 1L, rowMeans(S))
  d <- rep(0L, N)
  score <- function(v, ref, dd) {
    sv <- shift_profile(v, dd)
    ok <- !is.na(sv)
    if (sum(ok) < 5L) return(-Inf)
    den <- sqrt(sum(sv[ok]^2) * sum(ref[ok]^2))
    if (den <= 0) return(-Inf)
    sum(sv[ok] * ref[ok]) / den
  }
  for (it in seq_len(iters)) {
    A <- t(vapply(seq_len(N), function(i) {
      v <- shift_profile(S[i, ], d[i]); v[is.na(v)] <- 0; v
    }, numeric(m)))
    for (i in seq_len(N)) {
      ref <- colMeans(A[-i, , drop = FALSE])
      d[i] <- cand[which.max(vapply(cand, function(dd) score(S[i, ], ref, dd),
        numeric(1)))]
    }
  }
  as.integer(d)
}
transport_stack <- function(deltas) lapply(deltas, shift_transport)

# A realistic anatomical warp does not land a native node on exactly one group
# node: it spreads partial volume over neighbours.  Mass falling outside the
# group node set goes to the sink, so the row is still stochastic.
soft_transport <- function(delta, kappa = 0.7) {
  M <- matrix(0, m, m + 1L)
  for (x in seq_len(m)) {
    c0 <- x - delta
    w <- exp(-(seq_len(m) - c0)^2 / (2 * kappa^2))
    outside <- sum(exp(-(setdiff(-6:(m + 6L), seq_len(m)) - c0)^2 / (2 * kappa^2)))
    tot <- sum(w) + outside
    M[x, seq_len(m)] <- w / tot
    M[x, m + 1L] <- outside / tot
  }
  M
}
soft_stack <- function(deltas, kappa = 0.7) lapply(deltas, soft_transport, kappa = kappa)

flatten <- function(Pl, runs) {
  t(vapply(seq_len(N), function(i) {
    as.numeric(crossprod(Pl[[i]], runs[[i]])[seq_len(m), , drop = FALSE])
  }, numeric(m * Pd)))
}
consensus_share <- function(Pl, runA, runB) {
  TA <- flatten(Pl, runA); TB <- flatten(Pl, runB)
  VW <- mean(rowSums(TA * TB))                            # cross-partition
  Tb <- (TA + TB) / 2
  VC <- (sum(colSums(Tb)^2) - sum(rowSums(Tb^2))) / (N * (N - 1L))
  c(VC = VC, VW = VW, share = VC / VW)
}

P_anat <- transport_stack(rep(0L, N))                     # anatomy = identity
d_hat1 <- est_shift(R1)                                   # trained on run 1
P_func <- transport_stack(d_hat1)
d_hat2 <- est_shift(R2)                                   # trained on run 2 (circular)
P_circ <- transport_stack(d_hat2)
P_null <- transport_stack(sample(-2:2, N, replace = TRUE))
P_oracle <- transport_stack(delta_true)
P_soft <- soft_stack(d_hat1)                              # fractional P^F

# A global shift of every subject is unidentifiable (it only moves which nodes
# fall in the sink), so alignment quality is measured on CENTERED shifts.
centered_err <- function(d) (d - mean(d)) - (delta_true - mean(delta_true))
say("  shift recovery from run 1: mean |centered error| = ",
  sprintf("%.3f", mean(abs(centered_err(d_hat1)))), " nodes ; ",
  sum(abs(centered_err(d_hat1)) < 0.5), " of ", N, " within half a node")
say("  (raw agreement is not the right metric: a common offset is unidentifiable)")
sa <- consensus_share(P_anat, R2, R3)
sf <- consensus_share(P_func, R2, R3)
so <- consensus_share(P_oracle, R2, R3)
sn <- consensus_share(P_null, R2, R3)
report <- function(nm, s) say(sprintf("  %-22s V^C = %+10.4f   V^W = %+10.4f   share = %+.6f",
  nm, s[["VC"]], s[["VW"]], s[["share"]]))
report("anatomical P^A", sa)
report("functional P^F (run 1)", sf)
report("oracle P (true shifts)", so)
report("null P (random shifts)", sn)
ssoft <- consensus_share(P_soft, R2, R3)
report("fractional P^F (soft)", ssoft)
say("  eta_transport = share(P^F) - share(P^A) = ",
  sprintf("%+.6f", sf[["share"]] - sa[["share"]]))
say("  ceiling: share(oracle) - share(P^A)     = ",
  sprintf("%+.6f", so[["share"]] - sa[["share"]]))

# ---------------------------------------------------------------------------
# P7.b  cross-fit provenance is not optional
# ---------------------------------------------------------------------------
rule("P7.b  re-using the evaluation run to fit P inflates eta")

sc <- consensus_share(P_circ, R2, R3)
report("circular P^F (run 2)", sc)
say("  eta cross-fitted (trained on run 1, evaluated on runs 2-3) = ",
  sprintf("%+.6f", sf[["share"]] - sa[["share"]]))
say("  eta circular     (trained on run 2, evaluated on runs 2-3) = ",
  sprintf("%+.6f", sc[["share"]] - sa[["share"]]))
say("  inflation factor = ", sprintf("%.2fx",
  (sc[["share"]] - sa[["share"]]) / (sf[["share"]] - sa[["share"]])))
say("  -> without cross-fit provenance on the record, the two numbers are")
say("     indistinguishable in a result object.  Provenance is REQUIRED.")

# ---------------------------------------------------------------------------
# P7.c  eta may be negative
# ---------------------------------------------------------------------------
rule("P7.c  a null transport gives a negative eta, reported as-is")

say("  eta(null P) = ", sprintf("%+.6f", sn[["share"]] - sa[["share"]]))
etas_null <- replicate(200L, {
  consensus_share(transport_stack(sample(-2:2, N, replace = TRUE)), R2, R3)[["share"]]
}) - sa[["share"]]
say("  over 200 random null transports: mean eta = ",
  sprintf("%+.6f", mean(etas_null)),
  " ; ", sprintf("%.1f%%", 100 * mean(etas_null < 0)), " negative")
say("  null band: sd = ", sprintf("%.6f", sd(etas_null)),
  " , q95 = ", sprintf("%+.6f", unname(quantile(etas_null, 0.95))),
  " , max = ", sprintf("%+.6f", max(etas_null)))
say("  the honest eta = ", sprintf("%+.6f", eta_honest <- sf[["share"]] - sa[["share"]]),
  " exceeds ", sum(eta_honest > etas_null), " of 200 null draws")
say("  -> THIS is what makes +0.05 mean something.  A bare eta with no null")
say("     band is the uninterpretable headline P7.d warns about.")
say("  -> clipping eta at zero would turn a failed alignment into a null")
say("     result and repeat the sec.6 error at the population level.")

# ---------------------------------------------------------------------------
# P7.d  eta alone is not interpretable
# ---------------------------------------------------------------------------
rule("P7.d  a sink-heavy transport raises the share without aligning anything")

# Keep only the four central group nodes; everything else goes to the sink.
P_sinkheavy <- lapply(seq_len(N), function(i) {
  M <- matrix(0, m, m + 1L)
  for (x in seq_len(m)) {
    if (x >= 6L && x <= 7L) M[x, x] <- 1 else M[x, m + 1L] <- 1
  }
  M
})
ss <- consensus_share(P_sinkheavy, R2, R3)
report("sink-heavy P (no shift)", ss)
eta_sink <- ss[["share"]] - sa[["share"]]
eta_func <- sf[["share"]] - sa[["share"]]
say("  eta(sink-heavy) = ", sprintf("%+.6f", eta_sink),
  " vs eta(honest functional) = ", sprintf("%+.6f", eta_func),
  "   -> sink-heavy is ", sprintf("%.2fx", eta_sink / eta_func), " the honest gain")
sink_frac <- function(Pl) mean(vapply(Pl, function(M) mean(M[, m + 1L]), numeric(1)))
say("  but its unmapped native territory is ",
  sprintf("%.1f%%", 100 * sink_frac(P_sinkheavy)), " vs ",
  sprintf("%.1f%%", 100 * sink_frac(P_func)), " for P^F.")
say("  -> eta must NEVER be reported without the coverage diagnostics; a")
say("     transport can buy consensus by discarding the disagreeing nodes.")

# ---------------------------------------------------------------------------
# P7.e  the diagnostics slice
# ---------------------------------------------------------------------------
rule("P7.e  transport diagnostics: displacement, entropy, sink, coverage")

native_coord <- seq_len(m)                 # 1-d template coordinates
group_coord <- seq_len(m)

diagnostics <- function(Pl, ledger) {
  # ledger: N x m matrix of native per-node budgets (for the data-dependent parts)
  disp <- entro <- perp <- numeric(0)
  wts <- numeric(0)
  sink_terr <- sink_budget <- numeric(N)
  reach <- matrix(0, N, m)
  for (i in seq_len(N)) {
    M <- Pl[[i]]
    Gm <- M[, seq_len(m), drop = FALSE]
    rs <- rowSums(Gm)
    sink_terr[i] <- mean(M[, m + 1L])                       # data-free
    sink_budget[i] <- sum(M[, m + 1L] * ledger[i, ])        # data-dependent
    keep <- rs > 0
    Pt <- Gm[keep, , drop = FALSE] / rs[keep]               # renormalized rows
    cen <- as.numeric(Pt %*% group_coord)
    disp <- c(disp, abs(cen - native_coord[keep]))
    wts <- c(wts, rs[keep])
    h <- apply(Pt, 1L, function(p) {
      p <- p[p > 0]
      -sum(p * log(p))
    })
    entro <- c(entro, h); perp <- c(perp, exp(h))
    reach[i, ] <- colSums(Gm) > 0
  }
  list(
    displacement = c(median = median(disp), p90 = unname(quantile(disp, 0.9)),
      max = max(disp), mass_weighted_mean = sum(wts * disp) / sum(wts)),
    entropy_nats = mean(entro), perplexity = mean(perp),
    rows_all_sink = sum(vapply(Pl, function(M) sum(M[, m + 1L] == 1), numeric(1))),
    sink_territory = mean(sink_terr), sink_budget = mean(sink_budget),
    group_node_subject_coverage = colSums(reach)
  )
}
ledger2 <- t(vapply(R2, function(Z) rowSums(Z^2), numeric(m)))  # a nonneg ledger
for (nm in c("P_anat", "P_func", "P_soft", "P_sinkheavy")) {
  d <- diagnostics(get(nm), ledger2)
  say("\n  ", nm, ":")
  say("    displacement (template units): median ",
    sprintf("%.3f", d$displacement[["median"]]), " , p90 ",
    sprintf("%.3f", d$displacement[["p90"]]), " , max ",
    sprintf("%.3f", d$displacement[["max"]]), " , mass-weighted mean ",
    sprintf("%.3f", d$displacement[["mass_weighted_mean"]]))
  say("    row entropy = ", sprintf("%.4f", d$entropy_nats),
    " nats ; perplexity = ", sprintf("%.4f", d$perplexity),
    "  (effective group nodes per native node)")
  say("    unmapped native territory (data-free) = ",
    sprintf("%.2f%%", 100 * d$sink_territory),
    " ; rows entirely in sink = ", d$rows_all_sink, " of ", N * m)
  say("    sink budget (data-dependent, mean over subjects) = ",
    sprintf("%.4f", d$sink_budget))
  say("    subjects reaching each group node: min ",
    min(d$group_node_subject_coverage), " , max ",
    max(d$group_node_subject_coverage), " of ", N,
    "   (nodes below 2 subjects: ",
    sum(d$group_node_subject_coverage < 2L), ")")
}

# ---------------------------------------------------------------------------
# P8  transported coherent is a ledger, not group-node coherence
# ---------------------------------------------------------------------------
rule("P8  the transported coherent ledger is not group-node coherence")

nv <- 20L                                  # template voxels
qc <- 3L
B <- matrix(rnorm(qc * nv), qc, nv)        # one relation estimate over the template
# native frame: radius-1 searchlights, column-normalized (conservative)
Wn <- outer(seq_len(nv), seq_len(nv), function(a, b) as.numeric(abs(a - b) <= 1))
Wn <- Wn / matrix(colSums(Wn), nv, nv, byrow = TRUE)
say("  native frame: ", nv, " radius-1 searchlights, conservative")
say("    max |colSums(W) - 1| = ", sci(max(abs(colSums(Wn) - 1))))

# group nodes: four contiguous regions of five template voxels
region <- rep(1:4, each = 5L)
Pt <- matrix(0, nv, 5L)                                   # 4 group nodes + sink
Pt[cbind(seq_len(nv), region)] <- 1                       # native node -> region
say("  transport is row-stochastic with the required sink column (P1.a): ",
  "max |rowSums - 1| = ", sci(max(abs(rowSums(Pt) - 1))),
  " ; sink mass = ", sum(Pt[, 5L]))

tot_x <- lapply(seq_len(nv), function(x) B %*% diag(Wn[x, ]) %*% t(B))
coh_x <- lapply(seq_len(nv), function(x) {
  a <- as.numeric(B %*% Wn[x, ])
  tcrossprod(a) / sum(Wn[x, ])
})
cfg_x <- Map(function(t, c) t - c, tot_x, coh_x)

acc <- function(L, j) Reduce(`+`, Map(function(M, w) w * M, L, Pt[, j]))
Wg <- crossprod(Pt[, 1:4, drop = FALSE], Wn)              # group frame: 4 x nv
say("  group frame w^G_j = sum_x P_xj w_x is itself conservative: ",
  "max |colSums - 1| = ", sci(max(abs(colSums(Wg) - 1))))

rel <- function(A, Bm) max(abs(A - Bm)) / max(abs(Bm))
for (j in 1:4) {
  tot_led <- acc(tot_x, j)
  coh_led <- acc(coh_x, j)
  cfg_led <- acc(cfg_x, j)
  tot_grp <- B %*% diag(Wg[j, ]) %*% t(B)
  a_g <- as.numeric(B %*% Wg[j, ])
  coh_grp <- tcrossprod(a_g) / sum(Wg[j, ])
  say(sprintf("  group node %d:", j))
  say("    max |transported TOTAL - group-node own total|   = ",
    sci(max(abs(tot_led - tot_grp))), "   (exact: total is linear in w)")
  say("    ledger identity |coh + cfg - total|              = ",
    sci(max(abs(coh_led + cfg_led - tot_led))))
  say("    transported coherent vs group-node coherent: rel = ",
    sprintf("%.2f%%", 100 * rel(coh_led, coh_grp)),
    "   (tr ledger ", sprintf("%.4f", sum(diag(coh_led))),
    " vs tr group ", sprintf("%.4f", sum(diag(coh_grp))), ")")
}
say("  -> the transported TOTAL is the group node's own total, exactly.")
say("     The transported COHERENT is a ledger of native-node common modes")
say("     carried to a group location; it is a different object from the")
say("     group node's own common mode even where the latter is defined.")
say("     In crossform the latter is not defined at all: transport maps NODES,")
say("     so no group frame and no group-node `a` exists.  Hence the names")
say("     must differ (`native_coherent_ledger`, not `coherent`).")

# ---------------------------------------------------------------------------
# P9  executable readiness predicates
# ---------------------------------------------------------------------------
rule("P9  what a conservative result must supply for WS-E to consume it")

required <- c(
  "per_row_family", "per_row_scale", "per_row_center", "per_row_label",
  "per_row_alpha", "packed_codec_frobenius", "normalization_declared",
  "conservation_certificate", "composition_and_root", "signed_layer_declared",
  "latent_projection_named"
)
ready <- function(x) vapply(required, function(f) isTRUE(x[[f]]), logical(1))

# NOTE: this block STIPULATES the interface; it does not measure the package.
# The two "today" columns are transcribed from `conservative-geometry-v1` sec.8
# (as amended by D2) and are asserted here so the required set is executable.
compliant <- as.list(setNames(rep(TRUE, length(required)), required))
# `frame_family(..., alpha, normalization = "conservative")` -- delivered by D2.
# Carries measurement/family/node/scale/center/alpha per row (R/frame.R,
# `.frame_family_member_index()`); still owes composition+root (G5/D6) and the
# named latent projection (G7/D7).
frame_family_route <- modifyList(compliant, list(
  composition_and_root = FALSE, latent_projection_named = FALSE
))
# The bare rbind() + additive_frame(..., "conservative") route: still works,
# still provenance blind (`test-frame-family.R`).
rbind_route <- modifyList(frame_family_route, list(
  per_row_family = FALSE, per_row_scale = FALSE, per_row_center = FALSE,
  per_row_label = FALSE, per_row_alpha = FALSE
))
tab <- rbind(required = ready(compliant),
  frame_family_today = ready(frame_family_route),
  bare_rbind_today = ready(rbind_route))
print(t(tab))
say("  unmet items: frame_family() ", sum(!ready(frame_family_route)),
  " of ", length(required), " ; bare rbind ", sum(!ready(rbind_route)),
  " of ", length(required))
say("  a WS-E plan must REFUSE a frame that fails any row of this table.")
say("  crossform refuses to supply, and requires as typed input: image")
say("  registration, functional-transport learning, resampling of subject")
say("  images, and any group frame over group features.")

say("\nDONE")
