# Package-route witness for `design/model-coordinate-geometry-contract.md`.
#
# Unlike the other scripts in this directory, this one LOADS crossform: it is
# not an independent court but the record of how the contract's package-route
# numbers were produced, on the README fixture, through public entry points
# only.  The base-R court is `model-coordinate-geometry.R`.  Ticket M7 promotes
# this script into `tests/testthat/test-model-coordinate-contract.R`.
#
# What it shows:
#   R1  lowering the relation through an effect_extractor with map t(Q)
#       reproduces t(Q) %*% G_x %*% Q from the materialized full form, for the
#       total, coherent and configuration components;
#   R2  one operator of the proposal's query bank reproduces one entry;
#   R3  geometry_spectrum() and latent_geometry() run unchanged on the
#       lowered form;
#   R4  the trace split at the peak searchlight is exact and the addressable
#       term equals tr(S_x).
#
# The basis is rank-revealing (contract section 1.1): the README category RDM
# has a rank-1 Gram, so a requested rank of 2 is clamped to 1 and no
# uncentered direction can enter Q.

# Run from the repository root, or point CROSSFORM_ROOT at it.
suppressPackageStartupMessages(pkgload::load_all(
  Sys.getenv("CROSSFORM_ROOT", "."), quiet = TRUE
))

example <- example_fmri_effects()
rel <- example$fit$relation
over <- cross_partitions(rel, independence = "independent", generalizes_over = "run")
plan <- plan_geometry(rel, at = example$frame, over = over)
G <- materialize_geometry(plan)
q_full <- length(rel$effects)
unsvec <- function(v, q) crossform:::.unsvec_symmetric(v, q)
upper <- function(S) S[upper.tri(S, diag = TRUE)]

# --- model basis: centered Gram -> rank-revealing factor -> rank-revealing SVD

D <- example$model_rdm
n <- nrow(D); H <- diag(n) - 1 / n
K <- -0.5 * H %*% D %*% H
e <- eigen(K, symmetric = TRUE)
keep <- which(e$values > 1e-10 * max(abs(e$values)))
keep <- keep[seq_len(min(2L, length(keep)))]          # requested 2, clamped
Tm <- e$vectors[, keep, drop = FALSE] %*% diag(sqrt(e$values[keep]), length(keep))
s <- svd(Tm); rr <- s$d > 1e-10 * s$d[[1L]]
Q <- s$u[, rr, drop = FALSE]
q <- ncol(Q)
cat("model Gram eigenvalues:", round(e$values, 6), "\n")
cat("requested rank 2, effective rank", length(keep), ", basis dimension q =", q,
  ", max |Q'1| =", formatC(max(abs(colSums(Q))), format = "e", digits = 1), "\n")
stopifnot(max(abs(colSums(Q))) < 1e-10)

# --- R1: lowering commutes -----------------------------------------------------

betas <- lapply(rel$partitions, function(p) relation_block(rel, p, seq_len(rel$n_features)))
names(betas) <- rel$partitions
basis_space <- effect_space(paste0("m", seq_len(q)), basis_id = "model-coordinates:readme-category")
extractor <- effect_extractor(t(Q), effects = basis_space, estimator = "model-basis")
rel_m <- relation(betas, extract = extractor, domain = example$domain)
over_m <- cross_partitions(rel_m, independence = "independent", generalizes_over = "run")
G_m <- materialize_geometry(plan_geometry(rel_m, at = example$frame, over = over_m))

r1 <- vapply(c("total", "coherent", "configuration"), function(component) {
  full <- geometry_component(G, component)
  lowered <- geometry_component(G_m, component)
  a <- t(apply(full, 1, function(v) upper(t(Q) %*% unsvec(v, q_full) %*% Q)))
  b <- t(apply(lowered, 1, function(v) upper(unsvec(v, q))))
  max(abs(a - b))
}, numeric(1))
stopifnot(all(r1 < 1e-12))
cat("R1 lowering vs projected full form, max abs diff:",
  paste(names(r1), formatC(r1, format = "e", digits = 2)), "\n")

# --- R2: one operator of the bank ---------------------------------------------

M11 <- Q[, 1] %*% t(Q[, 1])
v11 <- evaluate_geometry(plan, query = bilinear_query(M11))
direct <- apply(geometry_component(G, "total"), 1, function(v) drop(t(Q[, 1]) %*% unsvec(v, q_full) %*% Q[, 1]))
r2 <- max(abs(as.numeric(v11$values) - direct))
stopifnot(r2 < 1e-12)
cat("R2 bank operator (1,1) vs projected entry:", formatC(r2, format = "e", digits = 2), "\n")

# --- R3: existing readers on the lowered form ---------------------------------

spectrum <- geometry_spectrum(G_m)
latent <- latent_geometry(G_m)
stopifnot(inherits(spectrum, "effect_spectrum_view"), inherits(latent, "effect_latent_geometry"))
cat("R3 geometry_spectrum() and latent_geometry() ran on the lowered form;",
  "latent measurements clipped:", latent$projection$clipped, "of",
  latent$projection$measurements, "\n")

# --- R4: trace split at the peak ---------------------------------------------

packed <- geometry_component(G, "total")
packed_m <- geometry_component(G_m, "total")
peak <- which.max(rowSums(abs(packed_m)))
Gx <- unsvec(packed[peak, ], q_full)
S <- unsvec(packed_m[peak, ], q)
P <- Q %*% t(Q)
centered_total <- sum(diag(H %*% Gx %*% H))
addressable <- sum(diag(P %*% Gx))
orthogonal <- sum(diag((H - P) %*% Gx))
stopifnot(abs(centered_total - addressable - orthogonal) < 1e-10,
  abs(addressable - sum(diag(S))) < 1e-10)
cat(sprintf("R4 peak searchlight %d: centered total %.5f = addressable %.5f + orthogonal %.5f; tr(S_x) = %.5f\n",
  peak, centered_total, addressable, orthogonal, sum(diag(S))))
cat("R4 compressed spectrum at the peak:", formatC(eigen(S, symmetric = TRUE)$values, digits = 5), "\n")
cat("model-coordinate package route: all checks hold\n")
