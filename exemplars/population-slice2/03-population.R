#!/usr/bin/env Rscript
# 03-population.R -- two population fits over the same twelve subjects, one
# under P^A and one under P^F, with the identity checks that make the carrying
# auditable.
#
# This is slice 1's script at searchlight scale: 1241 group nodes instead of 3,
# ~65,000 native searchlight nodes per subject instead of 3 regions, and a
# transport that actually does something -- 5.3 mm of median displacement, a
# fifth of every subject's territory in the sink -- instead of a 3x3 identity.
# The algebra is the same algebra, which is the point: the identities below are
# the same identities slice 1 asserts, at the same 1e-12, on an operator that
# is no longer trivial.
#
# WHAT IS CLAIMED. The group form at 1241 consensus grey-matter nodes for five
# contrasts over the nine `trust` conditions, with between-subject standard
# errors; the two transported component ledgers under the names contract
# section 8.1 requires; a descriptive cross-fitted heterogeneity spectrum; and
# the same fit repeated under a functionally-informed transport so that 04 can
# ask whether the second one carries the population better.
#
# WHAT IS NOT. Twelve subjects is twelve subjects. The intervals are
# UNCALIBRATED and `population_uncertainty()` says so in every row it emits.
# The 6 + 6 age split is a design choice and no younger-versus-older contrast
# is computed anywhere. The older group carries substantially more head motion
# (01 censors up to 112 of 217 volumes for `sub-129`), which lowers their
# residual degrees of freedom and inflates their beta variance; this is
# recorded per subject and is a reason to read the group fit, not the subject
# spread, as the result.
#
# Environment: SLICE2_DIR.
# Output: data/derived/population.rds  (git-ignored; the two fits, for 04)
#         results/population-slice2-receipts.csv  (committed evidence)

source(file.path(
  if (nzchar(Sys.getenv("SLICE2_DIR"))) Sys.getenv("SLICE2_DIR") else
    normalizePath(dirname(sub("^--file=", "",
      grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))),
  "00-common.R"))
crossform_version <- load_crossform()
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
script_t0 <- Sys.time()

tr_path <- file.path(DERIVED_DIR, "transports.rds")
if (!file.exists(tr_path)) stop("run 02-transports.R first")
TR <- readRDS(tr_path)
present <- TR$subjects
group_nodes <- TR$group$nodes
n_group <- length(group_nodes)
say("Subjects: %s", paste(present, collapse = ", "))
say("Group nodes: %d", n_group)

receipts <- new_receipts()

## ---- Per-subject conservative geometry plans ------------------------------
# Rebuilt rather than carried from 02: a geometry plan holds the compiled frame
# and the relation, which together are ~100 MB per subject, and rebuilding is
# two seconds. What IS carried from 02 is the transport, so the operator the
# population layer sees is byte-identical to the one 02 certified.
say("\n== per-subject geometry plans ==")
ref_space <- reference_space()
DIM3 <- as.integer(dim(ref_space))
plans <- list(); fitinfo <- list()
for (s in present) {
  t0 <- Sys.time()
  f <- readRDS(path_betas(s))
  mask_arr <- array(FALSE, DIM3)
  mask_arr[f$voxel_linear] <- TRUE
  dom <- neuroim2_volume_domain(
    neuroim2::LogicalNeuroVol(mask_arr, ref_space),
    id = paste0("ds003745-", s))
  stopifnot(identical(as.character(dom$feature_ids), TR$native_nodes[[s]]))

  runs <- lapply(f$trust, function(B) {
    dimnames(B) <- list(TRUST_CONDITIONS, as.character(dom$feature_ids))
    B
  })
  eff <- effect_space(TRUST_CONDITIONS,
                      basis_id = "ds003745-trust-condition-means:v1",
                      units = "percent-signal-change")
  rel <- relation(runs, effects = eff, domain = dom)
  over <- cross_partitions(rel, independence = "independent",
                           generalizes_over = "run")
  frame <- compile_frame(
    searchlights(radius = SEARCHLIGHT_RADIUS_MM,
                 normalization = "conservative"), dom)
  cons <- frame_conservation(frame)
  stopifnot(cons$conserved)
  plans[[s]] <- plan_geometry(rel, at = frame, over = over)
  fitinfo[[s]] <- f
  record(receipts, s, "native_searchlight_nodes", length(dom$feature_ids),
         note = "conservative searchlights, one per covered voxel; rows of this subject's transport")
  record(receipts, s, "coverage_voxels_seven_run_intersection",
         f$n_coverage_voxels,
         note = "voxels inside all five trust and both sharedreward run brain masks")
  record(receipts, s, "degenerate_voxels_dropped", f$n_degenerate,
         note = "voxels with no signal or no variance in some run; removed from native territory")
  record(receipts, s, "frame_column_mass_max_deviation", cons$max_deviation,
         tolerance = FRAME_TOLERANCE,
         note = "max |sum_x w_xv - 1|: each voxel's mass is split exactly once among the searchlights containing it")
  record(receipts, s, "glm_residual_df_min", min(f$resid_df),
         note = "smallest per-run residual df across the five trust runs")
  record(receipts, s, "glm_censored_volumes_total", sum(f$n_spike),
         note = "volumes with framewise displacement > 0.5 mm, modelled as one indicator each, summed over five runs")
  record(receipts, s, "glm_mean_framewise_displacement", mean(f$mean_fd),
         note = "mean FD in mm averaged over the five trust runs")
  record(receipts, s, "glm_trust_nuisance_levels", max(f$trust_nuisance_levels),
         note = "trial types modelled but not read off as conditions in the trust runs; only missed_trial exists, so this is 0 or 1")
  record(receipts, s, "glm_shared_nuisance_levels",
         max(f$shared_nuisance_levels),
         note = "trial types modelled but not read off as conditions in the sharedreward runs: the six block levels, the two too-sparse neutral event types, and missed_trial. Fitting the seven conditions alone instead changes the resulting fingerprint beyond recognition (median voxelwise correlation 0.13 on sub-104), so these are modelled rather than left to contaminate the regressors of interest")
  say("  %-9s V=%6d  %.0fs", s, length(dom$feature_ids), elapsed(t0))
}

## ---- The two population plans ---------------------------------------------
say("\n== population plans ==")
plan_A <- plan_population(plans, TR$transports_A[present], model = ~ 1,
                          normalization = "none")
plan_F <- plan_population(plans, TR$transports_F[present], model = ~ 1,
                          normalization = "none")
# `plan_population()` re-sorts subjects by label; every downstream array uses
# that order, so it is read back rather than assumed.
subjects <- plan_A$subject_index$subject
stopifnot(identical(subjects, plan_F$subject_index$subject))
node_ids <- group_nodes
query_names <- rownames(QUERY_BANK)

say("  P^A: %d subjects, %d group nodes, semantics %s",
    nrow(plan_A$subject_index), n_group, plan_A$semantics)
say("  P^F: cross_fit provenance = %s",
    paste(plan_F$subject_index$cross_fit, collapse = "/"))
stopifnot(all(plan_F$subject_index$cross_fit == "task-sharedreward"))
stopifnot(all(plan_F$subject_index$provenance == "functional"))
stopifnot(all(plan_A$subject_index$provenance == "anatomical"))

## ---- Estimate --------------------------------------------------------------
say("\n== estimate_population ==")
t0 <- Sys.time()
fit_A <- estimate_population(plan_A, QUERY_BANK)
say("  P^A total component  %.0fs", elapsed(t0))
t0 <- Sys.time()
fit_F <- estimate_population(plan_F, QUERY_BANK)
say("  P^F total component  %.0fs", elapsed(t0))

coherent_A <- estimate_population(plan_A, QUERY_BANK, component = "coherent")
configuration_A <- estimate_population(plan_A, QUERY_BANK,
                                       component = "configuration")
say("  ledger names: %s / %s / %s", fit_A$ledger, coherent_A$ledger,
    configuration_A$ledger)

unc_A <- population_uncertainty(fit_A)
unc_F <- population_uncertainty(fit_F)

## ---- Whole-form readout and the commutation acceptance --------------------
say("\n== materialize_population ==")
t0 <- Sys.time()
form_A <- materialize_population(plan_A)
say("  P^A complete form  %.0fs", elapsed(t0))
t0 <- Sys.time()
form_F <- materialize_population(plan_F)
say("  P^F complete form  %.0fs", elapsed(t0))

rdm_edge_estimate <- function(form) {
  tab <- as.data.frame(rdm(form))
  # The readout emits one row per node PLUS the sink, and orients the pair
  # however the packed coordinate ran. A squared distance is symmetric in the
  # pair, so both orientations are the same edge and both are accepted; the
  # sink is dropped because it is not a group node.
  hit <- (tab$left == COMMUTATION_EDGE[1L] & tab$right == COMMUTATION_EDGE[2L]) |
         (tab$left == COMMUTATION_EDGE[2L] & tab$right == COMMUTATION_EDGE[1L])
  edge <- tab[hit & !tab$sink, , drop = FALSE]
  # With a richer group model `rdm()` emits one row per node per term and
  # `match()` would silently take the first; ~ 1 has one term, and this is what
  # says so.
  stopifnot(nrow(edge) == length(node_ids), !anyDuplicated(edge$node))
  edge$estimate[match(node_ids, edge$node)]
}

## ---- Identity checks ------------------------------------------------------
# Every number is read out of the objects above, never re-derived, so a change
# in the package moves these numbers.
say("\n== identity checks ==")

# Each subject's own native ledger, one column per query, computed once and
# reused for both transports: `contrast_energy()` on a subject plan takes one
# contrast at a time, and running the bank twice over would double the most
# expensive check here for no additional evidence.
t0 <- Sys.time()
native_ledger <- lapply(subjects, function(s) {
  cols <- lapply(query_names, function(q) {
    v <- contrast_energy(plans[[s]], QUERY_BANK[q, ])
    # `$index` on a contrast view is the measurement id vector itself, not a
    # data frame with a `measurement` column -- unlike the frame index and
    # unlike the population result's `$index`. Three shapes, one word.
    ids <- if (is.data.frame(v$index)) v$index$measurement else v$index
    stats::setNames(v$total, as.character(ids))
  })
  M <- do.call(cbind, cols)
  colnames(M) <- query_names
  M
})
names(native_ledger) <- subjects
say("  native ledgers (%d subjects x %d queries)  %.0fs",
    length(subjects), length(query_names), elapsed(t0))

identities <- list()
for (nm in c("A", "F")) {
  fit  <- if (nm == "A") fit_A  else fit_F
  plan <- if (nm == "A") plan_A else plan_F
  form <- if (nm == "A") form_A else form_F
  P    <- if (nm == "A") TR$transports_A else TR$transports_F

  ## (1) Budget preservation through the transport, per subject: the fit's own
  ##     certificate, which `estimate_population()` enforces at fit time.
  budget_deviation <- stats::setNames(fit$receipt$subjects$budget_deviation,
                                      fit$receipt$subjects$subject)

  ## (2) Executor plumbing. The group fit's per-subject response really is
  ##     `transport_values()` of that subject's own native ledger. Exact by
  ##     construction -- both sides query before transporting -- so this pins
  ##     the executor's row handling and is NOT the commutation claim.
  executor <- vapply(subjects, function(s) {
    nl <- native_ledger[[s]]
    carried <- transport_values(P[[s]],
      nl[as.character(P[[s]]$native_index$node), , drop = FALSE])
    max(abs(carried[node_ids, query_names, drop = FALSE] -
            fit$values[node_ids, query_names, s]))
  }, numeric(1))

  ## (3) Sum over group nodes AND the sink equals the transported native total.
  ##     Unlike slice 1 the sink is not empty here, so dropping it breaks this
  ##     identity by a fifth of the budget -- which is reported as (3b).
  ##
  ##     Asserted RELATIVELY, on the L1 norm of the carried ledger, which is
  ##     the scale `estimate_population()`'s own budget certificate uses
  ##     (`receipt$budget$scale == "relative_to_ledger_l1_norm"`). At slice 1
  ##     the region totals were O(1) and an absolute 1e-12 was the same
  ##     statement; here they run to O(100) over 1242 rows, and an absolute
  ##     1e-12 would be asserting ~1e-14 relative -- below the accumulation
  ##     floor of the sum itself, so it would be a tolerance on the summation
  ##     order rather than on the transport.
  relative_gap <- function(s, rows) {
    carried <- fit$values[rows, query_names, s, drop = FALSE]
    dev <- abs(colSums(carried[, , 1L, drop = FALSE]) -
               fit$receipt$native_total[s, query_names])
    scale <- colSums(abs(carried[, , 1L, drop = FALSE]))
    max(dev / pmax(scale, .Machine$double.eps))
  }
  all_rows <- seq_len(dim(fit$values)[1L])
  node_sum <- vapply(subjects, relative_gap, numeric(1), rows = all_rows)
  node_sum_no_sink <- vapply(subjects, relative_gap, numeric(1),
                             rows = node_ids)

  ## (4) The E6 sum-over-nodes identity at the group.
  ##     Also asserted relatively, for the reason given at (3).
  whole <- contribution(fit, by = rep("brain", n_group))
  group_response <- apply(fit$values, c(3L, 2L), sum)
  group_global <- qr.coef(qr(plan$model$matrix),
                          group_response)["(Intercept)", query_names]
  wv <- whole$values[, query_names, drop = FALSE]
  theta <- max(abs(colSums(wv) - group_global) /
               pmax(colSums(abs(wv)), .Machine$double.eps))
  wv_no_sink <- wv[!whole$index$sink, , drop = FALSE]
  theta_no_sink <- max(abs(colSums(wv_no_sink) - group_global) /
                       pmax(colSums(abs(wv_no_sink)),
                            .Machine$double.eps))

  ## (5) THE COMMUTATION ACCEPTANCE (contract claim 3).
  ##     `estimate_population()` contracts the friend - computer choice
  ##     contrast into ONE packed operator and transports one number per node:
  ##     query, then transport. `materialize_population()` transports all 45
  ##     packed coordinates and `rdm()` contracts that same edge out of them
  ##     afterwards: transport, then query. Different arithmetic, same answer.
  commutation <- max(abs(
    rdm_edge_estimate(form) -
      fit$coefficients[node_ids, COMMUTATION_QUERY, "(Intercept)"]))

  identities[[nm]] <- list(
    budget = budget_deviation, executor = executor,
    node_sum = node_sum, node_sum_no_sink = node_sum_no_sink,
    theta = theta, theta_no_sink = theta_no_sink,
    commutation = commutation, whole = whole)

  say("  [P^%s] budget %.2e | executor %.2e | node sum + sink %.2e | Theta %.2e | COMMUTATION %.2e",
      nm, max(abs(budget_deviation)), max(executor), max(node_sum), theta,
      commutation)
  say("        (sink deleted: node sum %.3f, Theta %.3f (relative) -- both broken, because",
      max(node_sum_no_sink), theta_no_sink)
  say("         the sink carries %.1f%% of the territory and is not decoration)",
      100 * mean(plan$subject_index$sink_territory))
}

## (6) coherent + configuration = total, on the transported ledgers. Cannot
##     fail: configuration is computed as total - coherent inside each
##     subject's own execution, before any transport. Contract 8.1 requires the
##     sentence to be stated, so it is stated as documentation of the
##     decomposition's shape, not as evidence about its content.
ledger_identity <- max(abs(coherent_A$values + configuration_A$values -
                           fit_A$values))
say("  coherent + configuration = total  %.2e  (rounding; cannot fail)",
    ledger_identity)

## ---- Prevalence and coverage ----------------------------------------------
say("\n== prevalence and group-node coverage ==")
prev_A <- population_prevalence(fit_A, coverage_floor = length(subjects))
prev_F <- population_prevalence(fit_F, coverage_floor = length(subjects))
# `$below_floor` is the VECTOR OF NODE NAMES under the floor, not a logical
# mask, so it is counted with length() and not summed.
say("  P^A coverage: min %d subjects, %d node(s) below the floor of %d",
    min(prev_A$coverage$minimum), length(prev_A$coverage$below_floor),
    length(subjects))
say("  P^F coverage: min %d subjects, %d node(s) below the floor of %d",
    min(prev_F$coverage$minimum), length(prev_F$coverage$below_floor),
    length(subjects))

## ---- Heterogeneity --------------------------------------------------------
# Cross-fitted, the only estimator the contract admits for a spectrum
# (section 6.4). Five runs per subject clears the four-partition floor.
say("\n== heterogeneity (cross-fitted) ==")
t0 <- Sys.time()
het_A <- heterogeneity(plan_A, modes = 2L)
say("  P^A  %.0fs; n_eff %.2f", elapsed(t0), het_A$latent$n_eff)
t0 <- Sys.time()
het_F <- heterogeneity(plan_F, modes = 2L)
say("  P^F  %.0fs; n_eff %.2f", elapsed(t0), het_F$latent$n_eff)

## ---- Record ---------------------------------------------------------------
diag <- TR$diagnostics
for (s in subjects) {
  d <- diag[diag$subject == s, ]
  i <- match(s, plan_A$subject_index$subject)
  record(receipts, s, "sink_territory_PA",
         plan_A$subject_index$sink_territory[i],
         note = "unmapped native territory under P^A: a property of the operator alone (contract 7.5), computable before any data")
  record(receipts, s, "sink_territory_PF",
         plan_F$subject_index$sink_territory[
           match(s, plan_F$subject_index$subject)],
         note = "the same under P^F; equal to P^A's by construction because P^F carries exactly the rows P^A placed")
  record(receipts, s, "sink_identical_max_gap", d$sink_identical_max_gap,
         tolerance = POPULATION_TOLERANCE,
         note = "max over native rows of |P^A sink - P^F sink|: the controlled-pair assertion that eta cannot be won by discarding territory")
  record(receipts, s, "all_sink_rows", d$all_sink_rows,
         note = "native rows with no group mass at all; excluded from displacement and entropy and counted here (contract 7.5)")
  record(receipts, s, "transport_support_mean_nodes", d$support_mean_nodes,
         note = "mean number of group nodes P^F may spread a placed row over")
  record(receipts, s, "fingerprint_similarity_sd", d$fingerprint_similarity_sd,
         note = "sd of the voxel-to-group-node fingerprint correlation over the support; with temperature 6 this sets how peaked P^F is")
  for (nm in c("A", "F")) {
    record(receipts, s, sprintf("displacement_median_P%s: mm", nm),
           d[[paste0("displacement_median_", nm)]],
           note = "median over native rows of the distance from the row's centre to the mass-weighted centroid of its group nodes (contract 7.5)")
    record(receipts, s, sprintf("displacement_p90_P%s: mm", nm),
           d[[paste0("displacement_p90_", nm)]], note = "90th percentile of the same")
    record(receipts, s, sprintf("displacement_max_P%s: mm", nm),
           d[[paste0("displacement_max_", nm)]], note = "maximum of the same")
    record(receipts, s, sprintf("displacement_masswt_mean_P%s: mm", nm),
           d[[paste0("displacement_masswt_mean_", nm)]],
           note = "row-mass-weighted mean of the same")
    record(receipts, s, sprintf("entropy_mean_P%s: nats", nm),
           d[[paste0("entropy_mean_", nm)]],
           note = "mean over native rows of -sum_j ptilde log ptilde over the renormalized group columns; identically 0 for the hard P^A")
  }
  record(receipts, s, "perplexity_mean_PF", d$perplexity_mean_F,
         note = "MEAN of exp(H) over rows -- the effective number of group nodes a row spreads into. Not exp(mean H), which is reported separately; contract 7.5 requires the label")
  record(receipts, s, "exp_mean_entropy_PF", d$exp_mean_entropy_F,
         note = "exp(MEAN H). Differs from the row mean of exp(H) by Jensen; both are reported because the contract says an implementation must say which it means")

  for (nm in c("A", "F")) {
    id <- identities[[nm]]
    record(receipts, s, sprintf("transport_budget_preservation_P%s", nm),
           id$budget[[s]], tolerance = POPULATION_TOLERANCE,
           note = "|transported total incl. sink - native total| / L1 norm of the native ledger (contract claim 2)")
    record(receipts, s, sprintf("transport_executor_agreement_P%s", nm),
           id$executor[[s]], tolerance = POPULATION_TOLERANCE,
           note = "|transport_values() of this subject's own contrast_energy() ledger - the group fit's response|. Exact by construction; pins the executor's plumbing, NOT the commutation claim")
    record(receipts, s, sprintf("transport_node_sum_plus_sink_P%s", nm),
           id$node_sum[[s]], tolerance = POPULATION_TOLERANCE,
           note = "max over queries of |sum over group nodes AND the sink - the transported native total|, relative to the L1 norm of the carried ledger (the scale the fit's own budget certificate uses)")
    record(receipts, s, sprintf("transport_node_sum_without_sink_P%s", nm),
           id$node_sum_no_sink[[s]],
           note = "the same relative gap with the sink row deleted. NOT an identity and not asserted: it is the size of the error a reader makes by ignoring the sink, and here that is about a fifth of the budget")
    fit <- if (nm == "A") fit_A else fit_F
    record(receipts, s, sprintf("transport_sink_budget_P%s: %s", nm,
                                HEADLINE_QUERY),
           fit$receipt$sink_budget[s, HEADLINE_QUERY],
           note = "sink budget in ledger units for the headline query: the part of this subject's contrast total that reached no group node")
  }
}

for (nm in c("A", "F")) {
  id <- identities[[nm]]
  fit <- if (nm == "A") fit_A else fit_F
  unc <- if (nm == "A") unc_A else unc_F
  prev <- if (nm == "A") prev_A else prev_F
  het <- if (nm == "A") het_A else het_F
  tag <- paste0("_P", nm)
  record(receipts, "group", paste0("budget_preservation_worst", tag),
         fit$receipt$budget$max_relative_deviation,
         tolerance = fit$receipt$budget$tolerance,
         note = sprintf("the fit's own certificate over every participant, scale \"%s\"",
                        fit$receipt$budget$scale))
  record(receipts, "group", paste0("commutation_claim3", tag),
         id$commutation, tolerance = POPULATION_TOLERANCE,
         note = paste0("THE acceptance. |query-then-transport - transport-then-query| at the group nodes: ",
                       "estimate_population() contracts ", COMMUTATION_QUERY,
                       " into one packed operator and transports one number per node; ",
                       "materialize_population() transports all 45 packed coordinates and rdm() ",
                       "contracts the (", paste(COMMUTATION_EDGE, collapse = ", "),
                       ") edge afterwards (contract claim 3)"))
  record(receipts, "group", paste0("theta_sum_over_nodes", tag),
         id$theta, tolerance = POPULATION_TOLERANCE,
         note = "aggregated ledger over group nodes and sink reproduces the group coefficient of the summed response (E6), relative to the L1 norm of the aggregated ledger")
  record(receipts, "group", paste0("theta_sum_without_sink", tag),
         id$theta_no_sink,
         note = "the same relative gap with the sink row deleted. NOT an identity: at slice 1 this was still zero because the region map had full coverage; here it is the share of the group ledger the sink carries")
  record(receipts, "group", paste0("residual_df", tag),
         unc$between$residual_df,
         note = "N - rank(X) for the group model ~ 1")
  record(receipts, "group", paste0("coverage_min_subjects", tag),
         min(prev$coverage$minimum),
         note = paste0("group-node subject coverage from population_prevalence(): ",
                       prev$coverage$definition))
  record(receipts, "group", paste0("coverage_nodes_below_floor", tag),
         length(prev$coverage$below_floor),
         note = sprintf("group nodes below the declared floor of %d subjects (contract 7.5; the floor itself is open maintainer decision 14.3)",
                        length(subjects)))
  record(receipts, "group", paste0("heterogeneity_n_eff", tag),
         het$latent$n_eff,
         note = "effective number of between-subject modes, cross-fitted; one draw, descriptive, not an estimate of a between-subject trace")
  for (q in query_names) {
    record(receipts, "group", paste0("group_ledger_over_all_nodes", tag, ": ", q),
           sum(id$whole$values[!id$whole$index$sink, q]),
           note = "the group contrast total aggregated over all group nodes, sink excluded")
    record(receipts, "group", paste0("group_ledger_sink", tag, ": ", q),
           id$whole$values[id$whole$index$sink, q],
           note = "the same aggregation's sink row: group-level territory that reached no node")
  }
}
record(receipts, "group", "n_subjects", length(subjects),
       note = "participants in the population plan")
record(receipts, "group", "n_group_nodes", n_group,
       note = "consensus grey-matter searchlight group nodes, sink excluded")
record(receipts, "group", "coverage_union_voxels", TR$group$union_voxels,
       note = "voxels covered by at least one subject's seven-run intersection mask")
record(receipts, "group", "coverage_consensus_voxels", TR$group$consensus_voxels,
       note = "voxels covered by all subjects; the group grid is drawn from these")
record(receipts, "group", "coverage_disagreement_share",
       1 - TR$group$consensus_voxels / TR$group$union_voxels,
       note = "share of the union covered for some subjects and not others: the coverage heterogeneity the sink exists for")
record(receipts, "group", "group_node_min_subject_reach",
       min(TR$group$subject_count),
       note = "subjects putting nonzero mass on the least-reached group node, computed from the operators alone (the data-free half of contract 7.5)")
record(receipts, "group", "coherent_plus_configuration_PA", ledger_identity,
       tolerance = POPULATION_TOLERANCE,
       note = "coherent + configuration - total on the transported ledgers; a rounding identity that cannot fail, recorded because contract 8.1 requires the shape to be stated")

utils::write.csv(receipts_frame(receipts),
                 file.path(RESULTS_DIR, "population-slice2-receipts.csv"),
                 row.names = FALSE)

saveRDS(list(subjects = subjects, group_nodes = group_nodes,
             query_names = query_names,
             fit_A = fit_A, fit_F = fit_F,
             unc_A = unc_A, unc_F = unc_F,
             het_A = het_A, het_F = het_F,
             identities = identities,
             crossform_version = crossform_version),
        file.path(DERIVED_DIR, "population.rds"))

say("\n03-population.R done in %.1f min", elapsed(script_t0) / 60)
