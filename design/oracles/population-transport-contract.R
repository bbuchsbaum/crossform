#!/usr/bin/env Rscript
# Oracle: the transport object, budget preservation, commutation, normalization
# ============================================================================
#
# Supports claims 1-4 of `design/population-form-contract.md`.  Nothing here
# proposes an API surface; it fixes what the transport object must carry and
# which algebraic laws the population layer may rely on.
#
#   P1.a  A well-formed transport is nonnegative, row-stochastic INCLUDING the
#         sink column, and has a group-node index plus a sink.  Three
#         well-formedness failures are exhibited and detected.
#   P1.b  `semantics = "budget"` and `semantics = "density"` are two different
#         linear maps built from the same P; density is defined by a declared
#         row-mass vector mu.  Both are linear in the data.
#   P2    Budget preservation: for budget semantics the transported group total
#         equals the native total minus the sink mass, exactly, for SIGNED c.
#   P3.a  Query, transport and the group fit act on three different tensor
#         axes, so with a subject-constant weight operator (OLS, or GLS with a
#         common Omega) every evaluation order agrees.
#   P3.b  Counterexample 1: weights that vary along the NODE axis break
#         (node-map o fit) = (fit o node-map).
#   P3.c  Counterexample 2: weights that vary along the QUERY-COORDINATE axis
#         break (query o fit) = (fit o query).
#   P4    Per-subject budget normalization: `none`, `unit_budget`,
#         `precision_weighted` are three different group estimands built from
#         the same subject ledgers; all three conserve against their own total.
#
# Pure matrix algebra; the crossform package is NOT loaded and no production
# contraction path is exercised.
#
# Run:  Rscript design/oracles/population-transport-contract.R

options(digits = 17)

say <- function(...) cat(..., "\n", sep = "")
rule <- function(title) say("\n", title, "\n", strrep("-", nchar(title)))
sci <- function(x) format(x, scientific = TRUE, digits = 4)

set.seed(20260820)

# ---------------------------------------------------------------------------
# Shared helpers: the packed symmetric codec (sqrt(2) off-diagonals)
# ---------------------------------------------------------------------------
# Frobenius-consistent: <svec A, svec B> = <A, B>_F.  Equivalence with the
# package codec `crossform:::.svec_symmetric` is established by
# `conservative-transport-readiness.R` O3.a (1.78e-15); it is re-derived
# from first principles here so this oracle stands alone.
svec_index <- function(q) {
  A <- matrix(0, q, q)
  idx <- which(upper.tri(A, diag = TRUE), arr.ind = TRUE)
  list(idx = idx, w = ifelse(idx[, 1L] == idx[, 2L], 1, sqrt(2)))
}
svec <- function(A, ix = svec_index(nrow(A))) {
  A[ix$idx] * ix$w
}
sym_draw <- function(q) {
  m <- matrix(rnorm(q * q), q, q)
  0.5 * (m + t(m))
}

# ---------------------------------------------------------------------------
# P1.a  A well-formed transport object
# ---------------------------------------------------------------------------
rule("P1.a  well-formedness of the transport object")

# `group` labels the m group nodes; column m+1 is the REQUIRED sink.
make_transport <- function(rows, m, semantics = "budget", mu = NULL,
                           provenance = list(kind = "anatomical")) {
  structure(list(
    matrix = rows,
    n_native = nrow(rows), n_group = m,
    sink_column = m + 1L,
    semantics = match.arg(semantics, c("budget", "density")),
    row_mass = if (is.null(mu)) rep(1, nrow(rows)) else mu,
    provenance = provenance
  ), class = "location_transport")
}

check_transport <- function(P, tol = 1e-12) {
  M <- P$matrix
  c(
    nonnegative = all(M >= 0),
    has_sink_column = ncol(M) == P$n_group + 1L,
    row_stochastic = max(abs(rowSums(M) - 1)) <= tol,
    semantics_declared = P$semantics %in% c("budget", "density"),
    row_mass_positive = all(P$row_mass > 0),
    provenance_kind = !is.null(P$provenance$kind),
    cross_fit_when_functional =
      !identical(P$provenance$kind, "functional") ||
        !is.null(P$provenance$cross_fit)
  )
}

n_native <- 12L
m <- 4L
# A fractional (non-indicator) transport: native node 5 splits across two group
# nodes, nodes 11-12 are only partially covered and put the remainder in sink.
M <- matrix(0, n_native, m + 1L)
assign_to <- c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4)
M[cbind(seq_len(n_native), assign_to)] <- 1
M[5L, ] <- c(0, 0.6, 0.4, 0, 0)
M[11L, ] <- c(0, 0, 0, 0.7, 0.3)   # 30% of this node's territory unmapped
M[12L, ] <- c(0, 0, 0, 0, 1)       # entirely unmapped
P_ok <- make_transport(M, m, provenance = list(
  kind = "anatomical", method = "nearest group node under a fixed warp",
  warp_id = "oracle:warp:v1"
))
say("  a well-formed transport (n_native = ", n_native, ", m = ", m,
  " + sink):")
print(check_transport(P_ok))
say("  max |rowSums - 1| = ", sci(max(abs(rowSums(M) - 1))))
say("  nonzeros = ", sum(M != 0), " of ", length(M),
  "  (sparsity ", sprintf("%.1f%%", 100 * mean(M == 0)), ")")

say("\n  four well-formedness failures, each detected:")
bad_no_sink <- make_transport(M[, seq_len(m), drop = FALSE], m)
say("    (a) sink column dropped   -> has_sink_column = ",
  check_transport(bad_no_sink)[["has_sink_column"]],
  " , row_stochastic = ", check_transport(bad_no_sink)[["row_stochastic"]],
  " ; worst row deficit ",
  sci(max(abs(rowSums(M[, seq_len(m), drop = FALSE]) - 1))))
M_short <- M; M_short[3L, ] <- M_short[3L, ] * 0.8
say("    (b) rows summing to 0.8   -> row_stochastic = ",
  check_transport(make_transport(M_short, m))[["row_stochastic"]])
M_neg <- M; M_neg[7L, 3L] <- -0.2; M_neg[7L, 4L] <- 1.2
say("    (c) a negative entry      -> nonnegative    = ",
  check_transport(make_transport(M_neg, m))[["nonnegative"]])
P_fun_bad <- make_transport(M, m, provenance = list(kind = "functional"))
say("    (d) functional P with no cross-fit provenance -> ",
  "cross_fit_when_functional = ",
  check_transport(P_fun_bad)[["cross_fit_when_functional"]])

# ---------------------------------------------------------------------------
# P1.b  budget and density are two linear maps built from the same P
# ---------------------------------------------------------------------------
rule("P1.b  budget vs density: two linear maps, one P")

# Declared row masses: native node territory.  mu == 1 recovers the
# `conservative-transport-readiness.R` O3.d "count" density exactly.
mu_fine <- rep(1, n_native)
transport_budget <- function(P, c_native) as.numeric(crossprod(P$matrix, c_native))
transport_density <- function(P, c_native) {
  num <- as.numeric(crossprod(P$matrix, c_native))
  den <- as.numeric(crossprod(P$matrix, P$row_mass))
  out <- ifelse(den > 0, num / den, NA_real_)
  # the sink is an accounting column and is ALWAYS reported in budget units
  out[P$sink_column] <- num[P$sink_column]
  out
}
c_native <- rnorm(n_native)
b <- transport_budget(P_ok, c_native)
d <- transport_density(P_ok, c_native)
say("  native ledger (signed)  : ", paste(sprintf("%+.3f", c_native), collapse = " "))
say("  budget  into g1..g4|sink: ", paste(sprintf("%+.4f", b), collapse = " "))
say("  density into g1..g4|sink: ", paste(sprintf("%+.4f", d), collapse = " "),
  "   (sink in budget units)")
# Density is still a LINEAR map of the data: D(1/(P^T mu)) P^T.
Lin <- diag(1 / as.numeric(crossprod(P_ok$matrix, mu_fine))) %*% t(P_ok$matrix)
lin_gap <- max(abs(as.numeric(Lin %*% c_native)[seq_len(m)] - d[seq_len(m)]))
say("  density equals the fixed linear map D(1/(P^T mu)) P^T applied to c: ",
  "max diff = ", sci(lin_gap))
say("  -> BOTH semantics are linear and data-independent, so both commute with")
say("     the query and the fit (P3).  What density gives up is conservation.")
say("  budget totals over group nodes = ", sprintf("%+.15f", sum(b[seq_len(m)])))
say("  density totals over group nodes= ", sprintf("%+.15f", sum(d[seq_len(m)])),
  "   (no conservation law)")

# ---------------------------------------------------------------------------
# P2  budget preservation, exactly, for signed ledgers
# ---------------------------------------------------------------------------
rule("P2  transported total = native total - sink mass")

native_total <- sum(c_native)
sink_mass <- b[P_ok$sink_column]
group_total <- sum(b[seq_len(m)])
say("  native total            = ", sprintf("%+.17f", native_total))
say("  sink mass               = ", sprintf("%+.17f", sink_mass))
say("  group total             = ", sprintf("%+.17f", group_total))
say("  |group + sink - native| = ",
  sci(abs(group_total + sink_mass - native_total)))
say("  |group - (native-sink)| = ",
  sci(abs(group_total - (native_total - sink_mass))))
say("  negative entries in the native ledger: ", sum(c_native < 0), " of ",
  n_native, "  -> the law is a SIGNED-sum law, not a mass law.")

# Over 500 random (P, c) draws, including fractional rows and sink mass.
worst <- 0
for (rep in seq_len(500L)) {
  nn <- sample(6:20, 1L)
  mm <- sample(2:6, 1L)
  R <- matrix(rexp(nn * (mm + 1L)), nn, mm + 1L)
  R[runif(length(R)) < 0.6] <- 0
  bad <- rowSums(R) == 0
  R[bad, mm + 1L] <- 1
  R <- R / rowSums(R)
  cc <- rnorm(nn, sd = 3)
  bb <- as.numeric(crossprod(R, cc))
  worst <- max(worst, abs(sum(bb[seq_len(mm)]) + bb[mm + 1L] - sum(cc)))
}
say("  worst |group + sink - native| over 500 random transports = ", sci(worst))

# Without the sink the loss is silent and proportional to the missing row mass.
R_nosink <- P_ok$matrix[, seq_len(m), drop = FALSE]
lost <- sum(c_native) - sum(crossprod(R_nosink, c_native))
say("\n  the same P with the sink column DELETED (rows no longer stochastic):")
say("    mass silently lost = ", sprintf("%+.9f", lost),
  "   (", sprintf("%.1f%%", 100 * abs(lost / native_total)),
  " of the native total)")
say("    unmapped native row mass = ",
  sprintf("%.4f", sum(P_ok$matrix[, P_ok$sink_column])), " of ", n_native,
  " rows -> visible in the sink, invisible without it.")

# ---------------------------------------------------------------------------
# P3.a  query / transport / fit act on three axes and commute
# ---------------------------------------------------------------------------
rule("P3.a  commutation under a subject-constant weight operator")

N <- 8L                                  # subjects
q <- 4L                                  # query coordinates (conditions)
Pdim <- q * (q + 1L) / 2L                # packed form width
mg <- 5L                                 # group nodes
ix <- svec_index(q)
H <- sym_draw(q); h <- svec(H, ix)       # an arbitrary linear query

# Each subject has its own native node count and its own transport.
subjects <- lapply(seq_len(N), function(i) {
  ni <- sample(9:16, 1L)
  R <- matrix(rexp(ni * (mg + 1L)), ni, mg + 1L)
  R[runif(length(R)) < 0.7] <- 0
  R[rowSums(R) == 0, mg + 1L] <- 1
  R <- R / rowSums(R)
  # Signed, indefinite native forms: the estimation layer of the contract.
  Y <- t(vapply(seq_len(ni), function(x) svec(sym_draw(q), ix), numeric(Pdim)))
  list(P = R, Y = Y, n = ni)
})
X <- cbind(intercept = 1, covariate = scale(rnorm(N))[, 1L])   # group design
ols <- function(Y) solve(crossprod(X), crossprod(X, Y))

# Order A: query -> transport -> fit
TA <- t(vapply(subjects, function(s) as.numeric(crossprod(s$P, s$Y %*% h)),
  numeric(mg + 1L)))
beta_A <- ols(TA)                                   # k x (mg+1)

# Order B: transport -> fit -> query
Z <- array(0, c(N, mg + 1L, Pdim))
for (i in seq_len(N)) Z[i, , ] <- crossprod(subjects[[i]]$P, subjects[[i]]$Y)
beta_B <- matrix(0, ncol(X), mg + 1L)
for (j in seq_len(mg + 1L)) beta_B[, j] <- ols(Z[, j, ]) %*% h

# Order C: transport -> query -> fit  (the third route)
TC <- t(vapply(seq_len(N), function(i) as.numeric(Z[i, , ] %*% h),
  numeric(mg + 1L)))
beta_C <- ols(TC)

say("  max |beta(query->transport->fit) - beta(transport->fit->query)| = ",
  sci(max(abs(beta_A - beta_B))))
say("  max |beta(query->transport->fit) - beta(transport->query->fit)| = ",
  sci(max(abs(beta_A - beta_C))))

# GLS with a common (subject-constant) Omega commutes too.
Om <- diag(N) + 0.3 * tcrossprod(rnorm(N)) / N
Om <- Om %*% t(Om)
Oi <- solve(Om)
gls <- function(Y) solve(crossprod(X, Oi %*% X), crossprod(X, Oi %*% Y))
gA <- gls(TA)
gB <- matrix(0, ncol(X), mg + 1L)
for (j in seq_len(mg + 1L)) gB[, j] <- gls(Z[, j, ]) %*% h
say("  GLS with a common Omega: max |orderA - orderB| = ",
  sci(max(abs(gA - gB))))
say("  -> the commuting class is `the weight operator is constant across BOTH")
say("     the node axis and the coordinate axis`; OLS is its default member.")

# ---------------------------------------------------------------------------
# P3.b  counterexample 1: weights varying along the NODE axis
# ---------------------------------------------------------------------------
rule("P3.b  per-node weights break (node-map o fit) = (fit o node-map)")

# Slice-1 geometry: all subjects share a native region set, and the node-axis
# map A aggregates regions into networks.  This is the only setting in which
# both orders are even defined, which is why the counterexample lives here.
nr <- 6L; nk <- 3L
A <- matrix(0, nr, nk)
A[cbind(seq_len(nr), c(1, 1, 2, 2, 3, 3))] <- 1                 # regions->networks
Yreg <- matrix(rnorm(N * nr), N, nr)                            # subject x region
w_common <- runif(N, 0.5, 2)                                    # per-SUBJECT only
w_node <- matrix(runif(N * nr, 0.2, 5), N, nr)                  # per-SUBJECT x NODE

wls <- function(y, w) {
  W <- diag(w)
  as.numeric(solve(crossprod(X, W %*% X), crossprod(X, W %*% y)))
}
fit_then_map <- function(wmat) {
  B <- vapply(seq_len(nr), function(x) wls(Yreg[, x], wmat[, x]),
    numeric(ncol(X)))                                            # k x nr
  B %*% A
}
map_then_fit <- function(wmat) {
  Yk <- Yreg %*% A
  wk <- (wmat %*% A) / matrix(colSums(A), N, nk, byrow = TRUE)   # mass-mean weight
  vapply(seq_len(nk), function(j) wls(Yk[, j], wk[, j]), numeric(ncol(X)))
}
Wc <- matrix(w_common, N, nr)
gap_common <- max(abs(fit_then_map(Wc) - map_then_fit(Wc)))
gap_node <- max(abs(fit_then_map(w_node) - map_then_fit(w_node)))
scale_node <- max(abs(fit_then_map(w_node)))
say("  weights constant across nodes: max |order gap| = ", sci(gap_common))
say("  weights varying across nodes : max |order gap| = ", sci(gap_node),
  "   (", sprintf("%.1f%%", 100 * gap_node / scale_node),
  " of the largest coefficient)")
say("  -> with per-node weights the ORDER is part of the estimand and must be")
say("     pinned in plan identity.  No choice of aggregated weight repairs it.")

# ---------------------------------------------------------------------------
# P3.c  counterexample 2: weights varying along the COORDINATE axis
# ---------------------------------------------------------------------------
rule("P3.c  per-coordinate weights break (query o fit) = (fit o query)")

w_coord <- matrix(runif(N * Pdim, 0.2, 5), N, Pdim)
j0 <- 2L                                   # one group node suffices
Zj <- Z[, j0, ]                            # N x Pdim
# fit then query, with a weight per coordinate
B_pc <- vapply(seq_len(Pdim), function(p) wls(Zj[, p], w_coord[, p]),
  numeric(ncol(X)))
lhs <- as.numeric(B_pc %*% h)
# query then fit, with the mass-weighted collapse of those same weights
w_q <- as.numeric((w_coord %*% (h^2)) / sum(h^2))
rhs <- wls(as.numeric(Zj %*% h), w_q)
say("  fit-then-query coefficients : ", paste(sprintf("%+.6f", lhs), collapse = " "))
say("  query-then-fit coefficients : ", paste(sprintf("%+.6f", rhs), collapse = " "))
say("  max |gap| = ", sci(max(abs(lhs - rhs))), "   (",
  sprintf("%.1f%%", 100 * max(abs(lhs - rhs)) / max(abs(lhs))),
  " of the largest coefficient)")
# and the same computation with a coordinate-constant weight agrees exactly
w_flat <- matrix(runif(N, 0.2, 5), N, Pdim)
B_flat <- vapply(seq_len(Pdim), function(p) wls(Zj[, p], w_flat[, p]),
  numeric(ncol(X)))
say("  same code with coordinate-CONSTANT weights: max |gap| = ",
  sci(max(abs(as.numeric(B_flat %*% h) -
    wls(as.numeric(Zj %*% h), w_flat[, 1L])))))

# ---------------------------------------------------------------------------
# P4  per-subject budget normalization (gap G9)
# ---------------------------------------------------------------------------
rule("P4  none / unit_budget / precision_weighted are three estimands")

Ns <- 6L
mg4 <- 4L
# Subjects with very different native budgets and different spatial profiles.
profile <- rbind(
  c(0.50, 0.25, 0.15, 0.10),
  c(0.45, 0.30, 0.15, 0.10),
  c(0.40, 0.30, 0.20, 0.10),
  c(0.05, 0.10, 0.25, 0.60),
  c(0.08, 0.12, 0.20, 0.60),
  c(0.05, 0.15, 0.25, 0.55)
)
budget <- c(4.0, 3.5, 3.0, 0.6, 0.5, 0.4)       # incommensurable totals
Ledger <- profile * budget                       # subject x group node
se_total <- c(0.9, 1.0, 0.8, 0.10, 0.09, 0.11)   # per-subject SE of the total
prec <- 1 / se_total^2

norm_none <- colMeans(Ledger)
norm_unit <- colMeans(Ledger / rowSums(Ledger))
norm_prec <- colSums(prec * Ledger) / sum(prec)

show_row <- function(nm, v) {
  say(sprintf("  %-19s %s   total = %+.9f", nm,
    paste(sprintf("%+.6f", v), collapse = " "), sum(v)))
}
say("  subject native totals: ", paste(sprintf("%.2f", rowSums(Ledger)),
  collapse = " "))
say(sprintf("  %-19s %s", "", "g1         g2         g3         g4"))
show_row("none", norm_none)
show_row("unit_budget", norm_unit)
show_row("precision_weighted", norm_prec)
say("  argmax group node: none = ", which.max(norm_none),
  " ; unit_budget = ", which.max(norm_unit),
  " ; precision_weighted = ", which.max(norm_prec))
rel <- function(a, b) max(abs(a - b)) / max(abs(a))
say("  max relative difference none vs unit_budget        = ",
  sprintf("%.2f%%", 100 * rel(norm_none, norm_unit)))
say("  max relative difference none vs precision_weighted = ",
  sprintf("%.2f%%", 100 * rel(norm_none, norm_prec)))

# All three conserve -- against their OWN total.
say("\n  each normalization conserves against its own total:")
say("    none:       |sum(ledger) - mean_i T_i|            = ",
  sci(abs(sum(norm_none) - mean(rowSums(Ledger)))))
say("    unit:       |sum(ledger) - 1|                     = ",
  sci(abs(sum(norm_unit) - 1)))
say("    precision:  |sum(ledger) - sum(pi_i T_i)/sum(pi)| = ",
  sci(abs(sum(norm_prec) - sum(prec * rowSums(Ledger)) / sum(prec))))

# precision_weighted degenerates to none when precisions are equal
# (Definitional, not a measurement: with pi_i constant the weighted mean IS the
# plain mean.  Recorded so the degenerate case is pinned, not as evidence.)
say("  precision_weighted with equal precisions == none (definitional): ",
  sci(max(abs(colSums(rep(1, Ns) * Ledger) / Ns - norm_none))))

# unit_budget is a ratio of correlated estimates and is unstable near T_i = 0
Led2 <- Ledger
Led2[6L, ] <- Led2[6L, ] * (0.02 / sum(Led2[6L, ]))   # a near-null subject
say("\n  unit_budget failure mode 1 -- a subject with almost no evidence")
say("  (estimated total 0.02 instead of 0.40) still contributes budget 1:")
say("    none        -> ", paste(sprintf("%+.5f", colMeans(Led2)), collapse = " "))
say("    unit_budget -> ",
  paste(sprintf("%+.5f", colMeans(Led2 / rowSums(Led2))), collapse = " "),
  "   (bit-identical to the un-weakened fixture)")
say("    max |unit_budget(weakened) - unit_budget(original)| = ",
  sci(max(abs(colMeans(Led2 / rowSums(Led2)) - norm_unit))))

say("\n  unit_budget failure mode 2 -- a signed subject total near zero.")
Led3 <- Ledger
Led3[6L, ] <- c(+0.30, +0.25, -0.50, -0.30)           # signed ledger, total < 0
say("    subject 6 native ledger = ",
  paste(sprintf("%+.4f", Led3[6L, ]), collapse = " "),
  "   total = ", sprintf("%+.4f", sum(Led3[6L, ])))
say("    its unit_budget shares  = ",
  paste(sprintf("%+.4f", Led3[6L, ] / sum(Led3[6L, ])), collapse = " "),
  "   (sum ", sprintf("%+.1f", sum(Led3[6L, ] / sum(Led3[6L, ]))), ")")
say("    -> nodes g1,g2 hold POSITIVE evidence and receive NEGATIVE shares,")
say("       and |share| reaches ", sprintf("%.1f",
    max(abs(Led3[6L, ] / sum(Led3[6L, ])))),
  ".  The divisor is a signed estimate (contract sec.6), so unit_budget is a")
say("       ratio of correlated estimates with an unbounded denominator.")

# The per-subject scalar commutes with transport.  This needs a REAL transport
# on real native ledgers, not a rescaling of the already-transported table --
# otherwise the check is vacuous.
set.seed(99)
Pn <- lapply(seq_len(Ns), function(i) {
  ni <- 6L + i
  R <- matrix(rexp(ni * (mg4 + 1L)), ni, mg4 + 1L)
  R[runif(length(R)) < 0.5] <- 0
  R[rowSums(R) == 0, mg4 + 1L] <- 1
  R / rowSums(R)
})
cn <- lapply(Pn, function(R) rnorm(nrow(R)))
sc <- vapply(cn, function(v) 1 / sum(v), numeric(1))     # unit_budget scalars
native_side <- t(vapply(seq_len(Ns), function(i) {
  as.numeric(crossprod(Pn[[i]], sc[i] * cn[[i]]))        # scale THEN transport
}, numeric(mg4 + 1L)))
group_side <- t(vapply(seq_len(Ns), function(i) {
  sc[i] * as.numeric(crossprod(Pn[[i]], cn[[i]]))        # transport THEN scale
}, numeric(mg4 + 1L)))
say("\n  a per-subject scalar commutes with transport and with the query")
say("  (", Ns, " subjects, native node counts ",
  paste(vapply(Pn, nrow, integer(1)), collapse = "/"), ", real sparse P):")
say("    max |scale-then-transport - transport-then-scale| = ",
  sci(max(abs(native_side - group_side))))
say("    it does NOT commute with the fit -- that is exactly why it is part of")
say("    the estimand and must be declared in plan identity.")

say("\nDONE")
