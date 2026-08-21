#!/usr/bin/env Rscript
# 05-eta-across-run.R -- the second independence axis.
#
# 04 measures eta_transport across TASK: the fingerprint comes from
# `sharedreward` and eta is evaluated on `trust`. That is the stronger
# independence -- different task, different block, no shared trials -- and it
# is the headline. This script measures the same quantity across RUN inside
# `trust`: the fingerprint comes from `trust` run 1, and eta is evaluated on
# runs {2,3} against {4,5}, with run 1 excluded from the evaluation entirely.
#
# WHY BOTH. DECISION.md section 5.3 promised the two levels, and risk 8 says
# what their disagreement would mean: "A negative eta on the across-task fit
# with a positive eta on the across-run fit would indicate the fingerprint is
# capturing task-specific rather than anatomical idiosyncrasy -- a finding, not
# a bug." The across-run fingerprint is allowed to know things about `trust`
# that the across-task fingerprint cannot: it was measured on the same nine
# conditions, in the same session, in the same subject. If a functional
# transport is going to help this dataset at all, this is the axis where it
# should show, and if it does not help here either then the honest reading is
# that a seven- or nine-condition fingerprint at this grid spacing does not
# carry cross-subject correspondence, rather than that the two tasks disagree.
#
# WHY THE PARTITION IS {2,3} vs {4,5} AND NOT SOMETHING TIDIER. Both sides of
# the held-out split must support a CROSS-VALIDATED conservative form, which
# needs at least two runs each. Five runs, one spent on the fingerprint, leaves
# four: 2 + 2 is the only split that works. Fitting on runs {1,2} as
# DECISION.md section 5.3 sketched would leave {3,4,5}, which cannot be halved
# into two two-run partitions -- so this script spends one run on the
# fingerprint rather than two. That is a weaker fingerprint and it is the
# honest price of keeping both evaluation sides cross-validatable.
#
# Everything else -- the support, the sink, the softmax, the temperature, the
# leave-one-subject-out atlas, the permutation null -- is 04's, from
# `eta-common.R`, so the two eta values are the same estimator on different
# data rather than two estimators.
#
# Environment: SLICE2_DIR, ETA_DRAWS.
# Output: results/population-slice2-eta-across-run.csv       (committed)
#         results/population-slice2-eta-across-run-null.csv  (committed)
#         appends to results/population-slice2-receipts.csv

SLICE2 <- if (nzchar(Sys.getenv("SLICE2_DIR"))) Sys.getenv("SLICE2_DIR") else
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])))
source(file.path(SLICE2, "00-common.R"))
crossform_version <- load_crossform()
source(file.path(SLICE2, "eta-common.R"))
script_t0 <- Sys.time()

FIT_RUN     <- 1L                 # the fingerprint source; never evaluated on
PARTITION_A <- c(2L, 3L)
PARTITION_B <- c(4L, 5L)
GROUP_FINGERPRINT_RADIUS_MM <- 6  # 02's value, restated so this script is readable alone
VW_FLOOR <- 0

n_draws <- if (nzchar(Sys.getenv("ETA_DRAWS"))) {
  as.integer(Sys.getenv("ETA_DRAWS"))
} else {
  ETA_NULL_DRAWS
}

TRANSPORTS <- readRDS(file.path(DERIVED_DIR, "transports.rds"))
present <- TRANSPORTS$subjects
n_group <- length(TRANSPORTS$group$nodes)
N <- length(present)
say("Subjects: %d | group nodes: %d", N, n_group)
say("P^F2 fitted on task-trust run %d; eta evaluated on {%s} vs {%s}",
    FIT_RUN, paste(PARTITION_A, collapse = ","),
    paste(PARTITION_B, collapse = ","))

## ---- The run-1 fingerprint ------------------------------------------------
# The nine `trust` condition betas from run 1 alone, centred and unit-normed.
# Run 1 appears nowhere in the evaluation, so this is a genuine cross-fit --
# but a within-task one, which is the point of running it beside 04.
say("\n== run-%d fingerprints ==", FIT_RUN)
fingerprints <- list(); degenerate <- list()
for (s in present) {
  f <- readRDS(path_betas(s))
  fp <- t(f$trust[[FIT_RUN]])
  fp <- fp - rowMeans(fp)
  nrm <- sqrt(rowSums(fp^2))
  bad <- !is.finite(nrm) | nrm <= .Machine$double.eps
  nrm[bad] <- 1
  fingerprints[[s]] <- fp / nrm
  degenerate[[s]] <- bad
  stopifnot(nrow(fp) == length(TRANSPORTS$native_nodes[[s]]))
}
say("  %d subjects x %d conditions", N, ncol(fingerprints[[present[1L]]]))

## The fingerprint-neighbourhood support (6 mm) is not the transport support
## (12 mm) and 02 did not save it, so it is rebuilt here from the group grid.
ref_space <- reference_space()
DIM3 <- as.integer(dim(ref_space))
SPACING <- as.numeric(neuroim2::spacing(ref_space))
gid <- integer(prod(DIM3))
gid[TRANSPORTS$group$lin] <- seq_len(n_group)
fp_offsets <- local({
  lim <- floor(GROUP_FINGERPRINT_RADIUS_MM / SPACING)
  g <- expand.grid(di = seq(-lim[1L], lim[1L]), dj = seq(-lim[2L], lim[2L]),
                   dk = seq(-lim[3L], lim[3L]))
  d <- sqrt((g$di * SPACING[1L])^2 + (g$dj * SPACING[2L])^2 +
            (g$dk * SPACING[3L])^2)
  cbind(as.matrix(g[d <= GROUP_FINGERPRINT_RADIUS_MM, , drop = FALSE]),
        dist = d[d <= GROUP_FINGERPRINT_RADIUS_MM])
})
fp_supports <- list()
for (s in present) {
  vox <- readRDS(path_betas(s))$voxel_index
  pieces <- vector("list", nrow(fp_offsets))
  for (o in seq_len(nrow(fp_offsets))) {
    ii <- vox[, 1L] + fp_offsets[o, "di"]
    jj <- vox[, 2L] + fp_offsets[o, "dj"]
    kk <- vox[, 3L] + fp_offsets[o, "dk"]
    ok <- ii >= 1L & ii <= DIM3[1L] & jj >= 1L & jj <= DIM3[2L] &
          kk >= 1L & kk <= DIM3[3L]
    if (!any(ok)) next
    g <- gid[ii[ok] + (jj[ok] - 1L) * DIM3[1L] +
             (kk[ok] - 1L) * DIM3[1L] * DIM3[2L]]
    hit <- g > 0L
    if (!any(hit)) next
    pieces[[o]] <- cbind(row = which(ok)[hit], col = g[hit])
  }
  fp_supports[[s]] <- do.call(rbind, pieces)
}
atlases <- build_loo_atlas(present, fingerprints, fp_supports, n_group)

## ---- Build P^F2 -----------------------------------------------------------
# Same support as P^F, same rows P^A placed, same temperature, same sink.
say("\n== P^F2 ==")
PA_mats <- lapply(present, function(s)
  TRANSPORTS$transports_A[[s]]$matrix[, seq_len(n_group), drop = FALSE])
names(PA_mats) <- present

make_PF2 <- function(s, perm = NULL) {
  softmax_transport_matrix(
    support = TRANSPORTS$support[[s]],
    fingerprint = fingerprints[[s]],
    atlas = atlases[[s]]$atlas, defined = atlases[[s]]$defined,
    degenerate = degenerate[[s]],
    n_native = length(TRANSPORTS$native_nodes[[s]]),
    n_group = n_group, temperature = PF_TEMPERATURE, perm = perm)
}
PF2_mats <- lapply(stats::setNames(present, present), make_PF2)

# The controlled-pair assertion, restated for this transport: P^F2's sink must
# equal P^A's row by row, or eta stops being a comparison and becomes a
# comparison-plus-a-discard (contract 7.4).
sink_gap <- max(vapply(present, function(s)
  max(abs(Matrix::rowSums(PF2_mats[[s]]) - Matrix::rowSums(PA_mats[[s]]))),
  numeric(1)))
say("  sink identical to P^A: max gap %.2e", sink_gap)
stopifnot(sink_gap <= POPULATION_TOLERANCE)

# Every P^F2 is declared to crossform with its own cross-fit provenance, so
# that the operators carry the record even though this script's arithmetic goes
# through the raw matrices. A transport fitted on data is `functional`
# whichever program built it (R/transport.R), and `location_transport()`
# REFUSES a functional method that does not name the partitions it was fitted
# on -- capability `cross_fit_provenance`. Declaring all twelve rather than one
# means all twelve go through that refusal path, which is the point of
# declaring them at all.
evaluation_partitions <- sprintf("task-trust-run-%d",
                                 c(PARTITION_A, PARTITION_B))
PF2_declared <- lapply(stats::setNames(present, present), function(s)
  location_transport(
    matrix = PF2_mats[[s]],
    native_index = data.frame(node = TRANSPORTS$native_nodes[[s]],
                              stringsAsFactors = FALSE),
    group_index = TRANSPORTS$group$index,
    semantics = "budget",
    provenance = list(
      method = "functional",
      details = sprintf(
        paste0("softmax(%g * r) over the group nodes within %g mm, where r is ",
               "the correlation between this voxel's nine-condition task-trust ",
               "run-%d fingerprint and the leave-one-subject-out group node ",
               "fingerprint within %g mm. Restricted to the rows P^A placed ",
               "within %g mm, so the sink is P^A's."),
        PF_TEMPERATURE, PF_SUPPORT_RADIUS_MM, FIT_RUN,
        GROUP_FINGERPRINT_RADIUS_MM, TRANSPORT_RADIUS_MM),
      cross_fit = sprintf("task-trust-run-%d", FIT_RUN),
      fitted_on = sprintf("ds003745 task-trust run %d", FIT_RUN),
      evaluated_on = sprintf("ds003745 task-trust runs %s",
                             paste(c(PARTITION_A, PARTITION_B),
                                   collapse = ", ")))))
for (s in present) {
  cf <- PF2_declared[[s]]$provenance$cross_fit
  stopifnot(is.character(cf), nzchar(cf),
            !any(evaluation_partitions %in% cf))
}
# The group nodes must be one shared object across subjects, exactly as
# `plan_population()` requires of the transports it admits.
stopifnot(all(vapply(present, function(s)
  identical(PF2_declared[[s]]$group_index, TRANSPORTS$group$index),
  logical(1))))
say("  %d operators declared; cross_fit = %s; evaluated on %s. Disjoint.",
    length(PF2_declared), PF2_declared[[present[1L]]]$provenance$cross_fit,
    paste(evaluation_partitions, collapse = ", "))

## ---- Held-out forms and the shares ----------------------------------------
say("\n== held-out native forms (runs {%s} vs {%s}) ==",
    paste(PARTITION_A, collapse = ","), paste(PARTITION_B, collapse = ","))
forms <- held_out_forms(present, TRANSPORTS$native_nodes,
                        list(A = PARTITION_A, B = PARTITION_B), ref_space,
                        reporter = message)

say("\n== consensus shares ==")
sA <- consensus_share(PA_mats, forms, present, VW_FLOOR)
s2 <- consensus_share(PF2_mats, forms, present, VW_FLOOR)
eta2 <- s2$R - sA$R
say("  P^A    V_C %+12.4f   V_W %+12.4f +/- %.4f   R %s",
    sA$V_C, sA$V_W, sA$V_W_se, format(sA$R, digits = 6))
say("  P^F2   V_C %+12.4f   V_W %+12.4f +/- %.4f   R %s",
    s2$V_C, s2$V_W, s2$V_W_se, format(s2$R, digits = 6))
say("  eta_transport (across run) = %s", format(eta2, digits = 6))

## ---- Null band ------------------------------------------------------------
# The identity-permutation check is structural and runs regardless of whether
# the shares turn out formable: it is a statement about the operator, not about
# the data.
identity_gap <- max(vapply(present, function(s)
  max(abs(make_PF2(s, seq_len(n_group)) - PF2_mats[[s]])), numeric(1)))
say("\n  identity-permutation check: %.2e", identity_gap)
stopifnot(identity_gap <= POPULATION_TOLERANCE)

formable <- is.finite(sA$R) && is.finite(s2$R)
null <- data.frame(draw = integer(0), V_C = numeric(0), V_W = numeric(0),
                   R = numeric(0), eta = numeric(0))
if (!formable) {
  # Contract 7.1, literally: V^W is a cross-partition product and can be zero
  # or negative when nothing reproduces, and the share is then NA rather than a
  # large number. A null band around an NA is a band of NaNs with a p-value
  # attached, which is worse than no band, so none is drawn -- and V^C and V^W
  # are reported anyway, because they are the evidence for the NA.
  say("\n== null band NOT FORMED ==")
  say("  V^W does not clear the declared floor of %g for at least one", VW_FLOOR)
  say("  transport, so R and eta are NA by contract 7.1 and no band is drawn.")
  say("  V^W per subject under P^A ranges %+.1f to %+.1f with se %.1f:",
      min(sA$V_W_per_subject), max(sA$V_W_per_subject), sA$V_W_se)
  say("  the mean is %.2f standard errors from zero, i.e. this axis has no",
      sA$V_W / sA$V_W_se)
  say("  reproducible energy to take a share OF, not a share that came out low.")
} else {
  say("\n== null band (%d draws) ==", n_draws)
  set.seed(ETA_NULL_SEED + 1L)
  null_rows <- vector("list", n_draws)
  t0 <- Sys.time()
  for (d in seq_len(n_draws)) {
    perm <- sample.int(n_group)
    mats <- lapply(stats::setNames(present, present), make_PF2, perm = perm)
    sN <- consensus_share(mats, forms, present, VW_FLOOR)
    null_rows[[d]] <- data.frame(draw = d, V_C = sN$V_C, V_W = sN$V_W,
                                 R = sN$R, eta = sN$R - sA$R)
    if (d %% 25L == 0L) {
      say("  %3d/%d draws  %.1f min", d, n_draws, elapsed(t0) / 60)
    }
  }
  null <- do.call(rbind, null_rows)
}
nb <- summarize_null(null, eta2)
if (formable) {
  say("\n  null eta: mean %+.6f  sd %.6f  q95 %+.6f  max %+.6f  (%.1f%% negative)",
      nb$mean, nb$sd, nb$q95, nb$max, 100 * nb$negative_share)
  say("  eta2 %+.6f exceeds %d of %d null draws  (one-sided p = %.4f)",
      eta2, nb$rank, nb$draws, nb$p)
}

## ---- Write ----------------------------------------------------------------
utils::write.csv(
  data.frame(
    quantity = c("V_C_PA_acrossrun", "V_W_PA_acrossrun", "V_W_PA_acrossrun_se",
                 "R_PA_acrossrun",
                 "V_C_PF2", "V_W_PF2", "V_W_PF2_se", "R_PF2",
                 "eta_transport_across_run", "shares_formable",
                 "eta_null_mean", "eta_null_sd", "eta_null_q95",
                 "eta_null_max", "eta_null_negative_share", "eta_null_draws",
                 "eta_rank_in_null", "eta_one_sided_p"),
    value = c(sA$V_C, sA$V_W, sA$V_W_se, sA$R,
              s2$V_C, s2$V_W, s2$V_W_se, s2$R, eta2, as.numeric(formable),
              nb$mean, nb$sd, nb$q95, nb$max, nb$negative_share, nb$draws,
              nb$rank, nb$p),
    stringsAsFactors = FALSE),
  file.path(RESULTS_DIR, "population-slice2-eta-across-run.csv"),
  row.names = FALSE)
utils::write.csv(null,
                 file.path(RESULTS_DIR,
                           "population-slice2-eta-across-run-null.csv"),
                 row.names = FALSE)

receipts <- new_receipts()
axis <- sprintf("across-run axis: P^F2 fitted on task-trust run %d, evaluated on runs {%s} vs {%s}",
                FIT_RUN, paste(PARTITION_A, collapse = ","),
                paste(PARTITION_B, collapse = ","))
record(receipts, "group", "V_C_PA_acrossrun", sA$V_C,
       note = paste("cross-participant energy under P^A on the across-run partitions;", axis))
record(receipts, "group", "V_W_PA_acrossrun", sA$V_W,
       note = paste("cross-partition energy under P^A on the across-run partitions;", axis))
record(receipts, "group", "V_W_PA_acrossrun_se", sA$V_W_se,
       note = "standard error of V_W_PA_acrossrun over the 12 per-participant cross-partition products. V^W is a MEAN, and whether it clears the contract 7.1 floor is a claim about a noisy average")
record(receipts, "group", "R_PA_acrossrun", sA$R,
       note = "consensus share under P^A on the across-run partitions; differs from R_PA because it is a different held-out split")
record(receipts, "group", "V_C_PF2", s2$V_C,
       note = paste("cross-participant energy under P^F2;", axis))
record(receipts, "group", "V_W_PF2", s2$V_W,
       note = paste("cross-partition energy under P^F2;", axis))
record(receipts, "group", "V_W_PF2_se", s2$V_W_se,
       note = "standard error of V_W_PF2 over the 12 per-participant cross-partition products")
record(receipts, "group", "eta_across_run_shares_formable",
       as.numeric(formable),
       note = "1 when V^W cleared the declared floor for both transports and a share could be formed, 0 when contract 7.1 made R and eta NA. When 0 there is no null band, and its absence is the finding")
record(receipts, "group", "R_PF2", s2$R,
       note = "consensus share under P^F2")
record(receipts, "group", "eta_transport_across_run", eta2,
       note = paste0("the SECOND independence axis (DECISION.md 5.3), signed and unclamped. ",
                     axis, ". Read beside eta_transport: DECISION.md risk 8 says a ",
                     "disagreement in sign between the two axes would say the fingerprint ",
                     "captures task-specific rather than anatomical idiosyncrasy"))
record(receipts, "group", "eta_across_run_null_sd", nb$sd,
       note = "sd of the across-run null band, same permutation randomization as 04")
record(receipts, "group", "eta_across_run_null_q95", nb$q95,
       note = "95th percentile of the across-run null band")
record(receipts, "group", "eta_across_run_rank_in_null", nb$rank,
       note = "number of across-run null draws eta2 exceeds")
record(receipts, "group", "eta_across_run_one_sided_p", nb$p,
       note = "(1 + #{null >= eta2}) / (1 + #null), descriptive at n = 12")
record(receipts, "group", "eta_across_run_sink_gap", sink_gap,
       tolerance = POPULATION_TOLERANCE,
       note = "max over native rows of |P^F2 group mass - P^A group mass|: the controlled-pair assertion for the second axis")
record(receipts, "group", "eta_across_run_identity_permutation_gap",
       identity_gap, tolerance = POPULATION_TOLERANCE,
       note = "P^F2 rebuilt under the identity permutation reproduces P^F2")

path <- file.path(RESULTS_DIR, "population-slice2-receipts.csv")
new_rows <- receipts_frame(receipts)
if (file.exists(path)) {
  old <- utils::read.csv(path, stringsAsFactors = FALSE)
  old <- old[!(old$quantity %in% new_rows$quantity & old$subject == "group"), ]
  new_rows <- rbind(old, new_rows)
}
utils::write.csv(new_rows, path, row.names = FALSE)

say("\n05-eta-across-run.R done in %.1f min", elapsed(script_t0) / 60)
