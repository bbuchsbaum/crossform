#!/usr/bin/env Rscript
# 04-eta-transport.R -- V^C, V^W, eta_transport, and the null band.
#
# This is the acceptance for population-form-v1 section 7, and it is the whole
# reason slice 2 exists. Slice 1 could not run it: a transport that is the
# identity on named regions has no alternative to be compared against.
#
# THE DEFINITION (contract section 7.1), implemented literally. For a transport
# P and two independent partitions A, B of the HELD-OUT data, with Z_i^(r) the
# concatenated svec'd transported forms over the m group nodes and the sink
# EXCLUDED:
#
#   V^W(P) = (1/N) sum_i   <Z_i^A, Z_i^B>            cross-partition
#   V^C(P) = (1/(N(N-1))) sum_{i != i'} <Zbar_i, Zbar_i'>   cross-participant
#   R(P)   = V^C(P) / V^W(P)                          consensus share
#   eta    = R(P^F) - R(P^A)                          signed, unclamped
#
# THE PARTITIONS AND WHY THEY ARE LEGAL. P^F was fitted on task `sharedreward`
# (02), so every one of the five `trust` runs is held out. A = runs {1,2},
# B = runs {3,4,5}. Both sides have at least two runs, which is what makes a
# cross-validated conservative form computable on each; that is precisely the
# property section 4 of DECISION.md used to eliminate the AOMIC datasets, and
# it is why the partition is 2/3 rather than something tidier. The plan's
# `cross_fit` provenance names `task-sharedreward`, and this script asserts
# that no partition it evaluates on appears there.
#
# THE NULL BAND (contract section 7.3; the randomization is open maintainer
# decision 14.5, so this slice picks one and records it). 200 draws, each
# permuting the group-node fingerprint atlas across group nodes: node j is
# given node pi(j)'s fingerprint. The permuted transport keeps the same 12 mm
# support, the same softmax, the same temperature, the same sink and the same
# amount of spreading -- everything except the correspondence between a
# voxel's function and its destination's function. So the band prices the
# smoothing that P^F's wider support buys it, and eta's rank within the band
# prices the fingerprint itself.
#
# eta IS REPORTED SIGNED AND UNCLAMPED. Contract section 7.3 measures that
# 89.5 % of null transports give a negative eta; a negative honest eta here
# would be a finding about the fingerprint, not a bug, and clamping it at zero
# would convert a demonstrated alignment failure into a null result.
#
# eta MAY NOT BE READ ALONE. Section 7.4 measures an adversarial transport
# reporting eta = +0.167 -- 3.3x the honest functional gain in the contract's
# own fixture -- while sending 83 % of native territory to the sink. The six
# section 7.5 diagnostics are computed in 02 and written beside these numbers
# in the same receipts file, and this slice additionally pins P^F's sink to
# P^A's so that particular cheat is unavailable by construction.
#
# Environment: SLICE2_DIR, ETA_DRAWS (override the draw count for a smoke run).
# Output: results/population-slice2-eta.csv        (committed)
#         results/population-slice2-eta-null.csv   (committed, the 200 draws)
#         appends to results/population-slice2-receipts.csv

SLICE2 <- if (nzchar(Sys.getenv("SLICE2_DIR"))) Sys.getenv("SLICE2_DIR") else
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])))
source(file.path(SLICE2, "00-common.R"))
crossform_version <- load_crossform()
# `held_out_forms()`, `consensus_share()`, `softmax_transport_matrix()` and
# `summarize_null()` are shared with 05, which measures the same eta on the
# across-run axis. They live in one file so that the two eta values are the
# same estimator on different data rather than two estimators.
source(file.path(SLICE2, "eta-common.R"))
script_t0 <- Sys.time()

PARTITION_A <- c(1L, 2L)
PARTITION_B <- c(3L, 4L, 5L)
# Declared criterion for "V^W bounded away from zero" (contract 7.1). V^W is
# itself a cross-partition product and can be zero or negative when nothing
# reproduces; the share is then NA, not a large number.
VW_FLOOR <- 0

n_draws <- if (nzchar(Sys.getenv("ETA_DRAWS"))) {
  as.integer(Sys.getenv("ETA_DRAWS"))
} else {
  ETA_NULL_DRAWS
}

TRANSPORTS <- readRDS(file.path(DERIVED_DIR, "transports.rds"))
present <- TRANSPORTS$subjects
group_nodes <- TRANSPORTS$group$nodes
n_group <- length(group_nodes)
N <- length(present)
say("Subjects: %d | group nodes: %d | held-out partitions: {%s} vs {%s}",
    N, n_group, paste(PARTITION_A, collapse = ","),
    paste(PARTITION_B, collapse = ","))

## ---- Cross-fit provenance is checked, not assumed -------------------------
# The plan refuses to evaluate eta on any partition named in `cross_fit`
# (contract 7.2). Here the evaluation partitions are `trust` runs and the
# fitting partition is a different task, so the check is a string comparison --
# but it is the string comparison that would have caught the circular fixture
# section 7.2 measures at 3.15x the honest gain, so it is made explicitly.
evaluation_partitions <- paste0("task-trust-run-", c(PARTITION_A, PARTITION_B))
for (s in present) {
  cf <- TRANSPORTS$transports_F[[s]]$provenance$cross_fit
  stopifnot(is.character(cf), length(cf) >= 1L, nzchar(cf))
  if (any(evaluation_partitions %in% cf)) {
    stop("P^F for ", s, " declares cross_fit on a partition eta evaluates on: ",
         paste(intersect(evaluation_partitions, cf), collapse = ", "))
  }
}
say("Cross-fit provenance: P^F fitted on %s; eta evaluated on %s. Disjoint.",
    paste(unique(unlist(lapply(TRANSPORTS$transports_F,
                               function(x) x$provenance$cross_fit))),
          collapse = ", "),
    paste(evaluation_partitions, collapse = ", "))

## ---- The held-out native forms, one per subject per partition -------------
# `materialize_geometry()` gives the complete packed form at every native
# searchlight: a V x 45 matrix of svec'd 9 x 9 condition geometries. Transport
# is linear in these, so once they exist ANY transport's Z is one sparse
# product away -- which is what makes a 200-draw null band affordable at all.
say("\n== held-out native forms ==")
ref_space <- reference_space()
forms <- held_out_forms(present, TRANSPORTS$native_nodes,
                        list(A = PARTITION_A, B = PARTITION_B), ref_space,
                        reporter = message)
packed_width <- ncol(forms[[present[1L]]][[1L]])

## ---- V^C, V^W, R ----------------------------------------------------------
# Z_i^(r) = P_i' F_i^(r) over the m group nodes, sink excluded.
group_block <- function(P) P$matrix[, seq_len(n_group), drop = FALSE]
PA_mats <- lapply(present, function(s) group_block(TRANSPORTS$transports_A[[s]]))
PF_mats <- lapply(present, function(s) group_block(TRANSPORTS$transports_F[[s]]))
names(PA_mats) <- names(PF_mats) <- present

say("\n== consensus shares ==")
t0 <- Sys.time()
sA <- consensus_share(PA_mats, forms, present, VW_FLOOR)
sF <- consensus_share(PF_mats, forms, present, VW_FLOOR)
eta <- sF$R - sA$R
say("  P^A   V_C %+12.4f   V_W %+12.4f +/- %.4f   R %s",
    sA$V_C, sA$V_W, sA$V_W_se, format(sA$R, digits = 6))
say("  P^F   V_C %+12.4f   V_W %+12.4f +/- %.4f   R %s",
    sF$V_C, sF$V_W, sF$V_W_se, format(sF$R, digits = 6))
say("  eta_transport = R(P^F) - R(P^A) = %+.6f   (%.0fs per pair)",
    eta, elapsed(t0))

## ---- The null band --------------------------------------------------------
# One draw = one permutation of the group-node fingerprint atlas. Everything
# else about P^F is held fixed, including the support, the sink and the
# temperature, so the band is a band around a transport that is P^F in every
# respect except knowing which node is which.
say("\n== null band (%d draws, permuted group-node fingerprints) ==", n_draws)
set.seed(ETA_NULL_SEED)

permuted_PF <- function(s, perm) {
  softmax_transport_matrix(
    support = TRANSPORTS$support[[s]],
    fingerprint = TRANSPORTS$fingerprint[[s]],
    atlas = TRANSPORTS$atlas[[s]]$atlas,
    defined = TRANSPORTS$atlas[[s]]$defined,
    degenerate = TRANSPORTS$degenerate_fp[[s]],
    n_native = length(TRANSPORTS$native_nodes[[s]]),
    n_group = n_group, temperature = PF_TEMPERATURE, perm = perm)
}

## Before drawing anything: the null machinery is checked against the real
## thing. Rebuilding P^F under the IDENTITY permutation must reproduce, to the
## bit, the operator 02 built and `location_transport()` sealed. Without this
## the null band could be a band around a subtly different estimator and eta's
## rank inside it would mean nothing.
identity_gap <- max(vapply(present, function(s)
  max(abs(permuted_PF(s, seq_len(n_group)) - PF_mats[[s]])), numeric(1)))
say("  identity-permutation check: %.2e  (must be <= %g)",
    identity_gap, POPULATION_TOLERANCE)
stopifnot(identity_gap <= POPULATION_TOLERANCE)

null_rows <- vector("list", n_draws)
t0 <- Sys.time()
for (d in seq_len(n_draws)) {
  perm <- sample.int(n_group)
  mats <- lapply(stats::setNames(present, present), permuted_PF, perm = perm)
  sN <- consensus_share(mats, forms, present, VW_FLOOR)
  null_rows[[d]] <- data.frame(draw = d, V_C = sN$V_C, V_W = sN$V_W,
                               R = sN$R, eta = sN$R - sA$R)
  if (d %% 25L == 0L) {
    say("  %3d/%d draws  %.1f min elapsed", d, n_draws, elapsed(t0) / 60)
  }
}
null <- do.call(rbind, null_rows)
nb <- summarize_null(null, eta)

# WHERE THE HONEST TRANSPORT SITS IN THE NULL, COMPONENT BY COMPONENT.
# eta's rank in the band says the transport is unusual. These two say WHICH
# WAY, and they are the difference between "P^F carries the population better"
# and "P^F suppressed the denominator". A transport whose V^C is at the null
# median has done nothing for cross-participant consensus that a scrambled
# fingerprint would not have done; if its eta is nonetheless extreme, the eta
# is coming from V^W.
vc_pctl <- mean(null$V_C < sF$V_C)
vw_pctl <- mean(null$V_W < sF$V_W)
say("\n  P^F V_C %+.2f sits at the %.1f%% percentile of the null V_C",
    sF$V_C, 100 * vc_pctl)
say("  P^F V_W %+.2f sits at the %.1f%% percentile of the null V_W",
    sF$V_W, 100 * vw_pctl)
null_mean <- nb$mean; null_sd <- nb$sd; null_q95 <- nb$q95
null_max <- nb$max; null_neg_share <- nb$negative_share
eta_rank <- nb$rank; eta_p <- nb$p

say("\n  null eta: mean %+.6f  sd %.6f  q95 %+.6f  max %+.6f  (%.1f%% negative)",
    null_mean, null_sd, null_q95, null_max, 100 * null_neg_share)
say("  honest eta %+.6f exceeds %d of %d null draws  (one-sided p = %.4f)",
    eta, eta_rank, nb$draws, eta_p)
say("  null band computed in %.1f min", elapsed(t0) / 60)

## ---- Write ----------------------------------------------------------------
eta_table <- data.frame(
  quantity = c("V_C_PA", "V_W_PA", "R_PA",
               "V_W_PA_se",
               "V_C_PF", "V_W_PF", "R_PF", "V_W_PF_se",
               "eta_transport",
               "eta_null_mean", "eta_null_sd", "eta_null_q95", "eta_null_max",
               "eta_null_negative_share", "eta_null_draws",
               "eta_rank_in_null", "eta_one_sided_p",
               "V_C_PF_percentile_in_null", "V_W_PF_percentile_in_null"),
  value = c(sA$V_C, sA$V_W, sA$R, sA$V_W_se,
            sF$V_C, sF$V_W, sF$R, sF$V_W_se, eta,
            null_mean, null_sd, null_q95, null_max, null_neg_share,
            nb$draws, eta_rank, eta_p, vc_pctl, vw_pctl),
  stringsAsFactors = FALSE
)
utils::write.csv(eta_table,
                 file.path(RESULTS_DIR, "population-slice2-eta.csv"),
                 row.names = FALSE)
utils::write.csv(null,
                 file.path(RESULTS_DIR, "population-slice2-eta-null.csv"),
                 row.names = FALSE)

## Append to the receipts written by 03.
receipts <- new_receipts()
note_common <- sprintf(
  paste0("held-out task-trust, partitions A={%s} B={%s}; P^F cross-fitted on ",
         "task-sharedreward (contract 7.1)"),
  paste(PARTITION_A, collapse = ","), paste(PARTITION_B, collapse = ","))
record(receipts, "group", "V_C_PA", sA$V_C,
       note = paste("cross-participant energy under P^A;", note_common))
record(receipts, "group", "V_W_PA", sA$V_W,
       note = paste("cross-partition energy under P^A, reported separately from the ratio (contract 7.3);",
                    note_common))
record(receipts, "group", "V_W_PA_se", sA$V_W_se,
       note = "standard error of V_W_PA over the 12 per-participant cross-partition products. V^W is a MEAN, so whether it clears the contract 7.1 floor is a claim about a noisy average, and the share R inherits whatever that average is worth")
record(receipts, "group", "R_PA", sA$R,
       note = "consensus share of reproducible energy under P^A = V_C/V_W")
record(receipts, "group", "V_C_PF", sF$V_C,
       note = paste("cross-participant energy under P^F;", note_common))
record(receipts, "group", "V_W_PF", sF$V_W,
       note = paste("cross-partition energy under P^F;", note_common))
record(receipts, "group", "V_W_PF_se", sF$V_W_se,
       note = "standard error of V_W_PF over the 12 per-participant cross-partition products")
record(receipts, "group", "R_PF", sF$R,
       note = "consensus share of reproducible energy under P^F = V_C/V_W")
record(receipts, "group", "eta_transport", eta,
       note = paste0("THE SECTION 7 ACCEPTANCE. R(P^F) - R(P^A), signed and unclamped. ",
                     "Uninterpretable without the six section 7.5 diagnostics recorded ",
                     "beside it in this same file (displacement, entropy/perplexity, ",
                     "sink territory and sink budget, group-node subject coverage, and ",
                     "V_C/V_W separately per transport)."))
record(receipts, "group", "eta_null_mean", null_mean,
       note = sprintf("mean eta over %d draws that permute the group-node fingerprint atlas across group nodes, holding support, softmax, temperature and sink fixed (contract 7.3; randomization method is open maintainer decision 14.5 and this is the choice made here)",
                      nb$draws))
record(receipts, "group", "eta_null_draws", nb$draws,
       note = "number of null draws with a finite eta; the band is only a band if this is the declared draw count")
record(receipts, "group", "eta_null_sd", null_sd,
       note = "sd of the null eta distribution")
record(receipts, "group", "eta_null_q95", null_q95,
       note = "95th percentile of the null eta distribution")
record(receipts, "group", "eta_null_max", null_max,
       note = "largest null eta drawn")
record(receipts, "group", "eta_null_negative_share", null_neg_share,
       note = "share of null draws with eta < 0; contract 7.3 measures 89.5 % on its own fixture")
record(receipts, "group", "eta_rank_in_null", eta_rank,
       note = "number of null draws the honest eta exceeds; this is what makes the honest eta mean something")
record(receipts, "group", "eta_one_sided_p", eta_p,
       note = "(1 + #{null >= eta}) / (1 + #null): a permutation p-value against the fingerprint-correspondence null, descriptive at n = 12")
record(receipts, "group", "V_C_PF_percentile_in_null", vc_pctl,
       note = paste0("share of null draws whose cross-participant energy is below the honest P^F's. ",
                     "THIS IS THE NUMBER THAT SAYS WHAT ETA IS MADE OF: a V^C at the null median means ",
                     "the true fingerprint built no more cross-participant consensus than a scrambled ",
                     "one, and any extreme eta beside it is coming from the denominator"))
record(receipts, "group", "V_W_PF_percentile_in_null", vw_pctl,
       note = paste0("share of null draws whose cross-partition energy is below the honest P^F's. ",
                     "A value at 0 means the true fingerprint destroyed more within-subject ",
                     "reproducible energy than any of the randomized transports did"))
record(receipts, "group", "eta_packed_width", packed_width,
       note = "svec width of the transported form at each group node; Z is this many coordinates times the group nodes, sink excluded")
record(receipts, "group", "eta_null_identity_permutation_gap", identity_gap,
       tolerance = POPULATION_TOLERANCE,
       note = "max over subjects and entries of |P^F rebuilt under the identity permutation - the P^F location_transport() sealed in 02|. Pins the null band to the same estimator the honest eta uses, so eta's rank inside the band is a rank in its own distribution")

path <- file.path(RESULTS_DIR, "population-slice2-receipts.csv")
new_rows <- receipts_frame(receipts)
if (file.exists(path)) {
  old <- utils::read.csv(path, stringsAsFactors = FALSE)
  old <- old[!(old$quantity %in% new_rows$quantity & old$subject == "group"), ]
  new_rows <- rbind(old, new_rows)
}
utils::write.csv(new_rows, path, row.names = FALSE)

say("\n04-eta-transport.R done in %.1f min", elapsed(script_t0) / 60)
