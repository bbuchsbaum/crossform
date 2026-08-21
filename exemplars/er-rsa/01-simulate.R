#!/usr/bin/env Rscript
# 01-simulate.R -- build the deterministic encoding-retrieval simulation and
# write both the data and the closed-form ground truth.
#
# Writes:
#   data/simulation.rds            items, probes, design, truth, responses
#                                  (git-ignored; regenerate by re-running)
#   results/design-items.csv       the 36-item study list and its covariate
#   results/design-probes.csv      the 30 retrieval probes, old and lure
#   results/design-summary.csv     the design in one table
#   results/planted-cell-summary.csv
#                                  planted pair values by region and cell
#                                  class, for cross-run and same-run edges
#
# Nothing here touches crossform. The simulation is the exemplar's own, so
# the recovery in 03-recover.R compares crossform against an independent
# construction rather than against itself.

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
source(file.path(exemplar_dir, "00-common.R"))
results <- file.path(exemplar_dir, "results")
data_dir <- file.path(exemplar_dir, "data")
dir.create(results, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Design -------------------------------------------------------------
items <- er_items()
probes <- er_probes(items)
design <- er_design(items, probes)
truth <- er_truth(items, probes)

message("items encoded: ", nrow(items),
        "  |  retrieval probes: ", nrow(probes),
        " (", sum(probes$old), " old, ", sum(!probes$old), " lures)")
message("encoded but never retrieved: ", sum(!items$retrieved),
        "  |  probes never encoded: ", sum(!probes$old))
message("runs: ", N_RUNS, " study-test cycles, ", N_TR, " TRs each; ",
        "residual df per run = ", design$residual_df)

## ---- Data ---------------------------------------------------------------
factor <- er_noise_factor()
responses <- lapply(seq_len(N_SUBJECTS), function(subject) {
  er_simulate_subject(design, truth, subject, factor)
})
names(responses) <- sprintf("sub-%02d", seq_len(N_SUBJECTS))
message("simulated subjects: ", N_SUBJECTS,
        "  |  voxels: ", N_VOXELS,
        " (", paste(sprintf("%s=%d", names(REGION_SIZES), REGION_SIZES),
                    collapse = ", "), ")")

saveRDS(
  list(
    items = items, probes = probes, design = design, truth = truth,
    responses = responses, seed = SEED, n_subjects = N_SUBJECTS
  ),
  file.path(data_dir, "simulation.rds")
)

## ---- Ground-truth tables ------------------------------------------------
items_out <- items
items_out$reinstatement_gain <- ifelse(
  items$retrieved, truth$gain[items$item], NA_real_)
items_out$study_duration_centered <-
  items$study_duration - truth$duration_center
write.csv(items_out, file.path(results, "design-items.csv"), row.names = FALSE)
write.csv(probes, file.path(results, "design-probes.csv"), row.names = FALSE)

summary_rows <- data.frame(
  quantity = c(
    "encoding effects (left axis)", "retrieval effects (right axis)",
    "encoded items never retrieved", "retrieval probes never encoded",
    "study-test cycles (runs)", "TRs per run", "design columns per run",
    "residual df per run", "voxels", "regionA voxels", "regionB voxels",
    "regionC voxels", "simulated subjects",
    "cross-run directed edges", "same-run directed edges",
    "reinstatement intercept (rho0)", "reinstatement slope (rho1, per s)",
    "study duration levels (s)", "study duration centering constant (s)",
    "same-run item state sd"
  ),
  value = c(
    nrow(items), nrow(probes), sum(!items$retrieved), sum(!probes$old),
    N_RUNS, N_TR, ncol(design$runs[[1L]]$design), design$residual_df,
    N_VOXELS, REGION_SIZES[["regionA"]], REGION_SIZES[["regionB"]],
    REGION_SIZES[["regionC"]], N_SUBJECTS,
    N_RUNS * (N_RUNS - 1L), N_RUNS,
    REINSTATEMENT_INTERCEPT, REINSTATEMENT_SLOPE,
    paste(range(STUDY_DURATIONS), collapse = "-"),
    round(truth$duration_center, 4), STATE_SD
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_rows, file.path(results, "design-summary.csv"),
  row.names = FALSE)

## ---- Planted structure, by region and cell class ------------------------
# The three cell classes are the ones the analysis has to tell apart:
#   match                   a probe and the encoding trial of the same item
#   control_same_category   a different item of the same category
#   control_other_category  a different item of a different category
cell_class <- outer(items$item, probes$source_item, function(a, b) {
  ifelse(!is.na(b) & a == b, "match", NA_character_)
})
same_category <- outer(items$category, probes$category, `==`)
cell_class[is.na(cell_class) & same_category] <- "control_same_category"
cell_class[is.na(cell_class)] <- "control_other_category"

planted_rows <- list()
for (region in names(REGION_SIZES)) {
  weights <- as.numeric(REGION_LABELS == region) / REGION_SIZES[[region]]
  for (edges in c("cross", "same")) {
    G <- er_planted_geometry(truth, weights, edges)
    for (class in c("match", "control_same_category",
                    "control_other_category")) {
      planted_rows[[length(planted_rows) + 1L]] <- data.frame(
        region = region, edges = edges, cell_class = class,
        cells = sum(cell_class == class),
        planted_mean = mean(G[cell_class == class]),
        stringsAsFactors = FALSE
      )
    }
  }
}
planted <- do.call(rbind, planted_rows)
planted$planted_mean <- round(planted$planted_mean, 6)
write.csv(planted, file.path(results, "planted-cell-summary.csv"),
  row.names = FALSE)

message("\nPlanted pair values (cross-run edges), mean over cells:")
print(planted[planted$edges == "cross", c("region", "cell_class", "cells",
  "planted_mean")], row.names = FALSE)
message("\nSame-run edges, the comparator the estimand excludes:")
print(planted[planted$edges == "same", c("region", "cell_class",
  "planted_mean")], row.names = FALSE)
message("\nwrote ", file.path(data_dir, "simulation.rds"))
