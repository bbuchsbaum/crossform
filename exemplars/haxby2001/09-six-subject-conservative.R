#!/usr/bin/env Rscript
# 09-six-subject-conservative.R -- the conservative ledger, six times.
#
# `07-conservative-geometry.R` runs the whole conservative apparatus on Haxby
# subject 1: detection vs attribution, the latent layer, the multiscale
# coherence spectrum, and two refusals. It is a four-panel argument about one
# subject. This script asks the narrower question that argument leaves open:
#
#   does the conservation identity hold on every subject, or did subject 1
#   happen to be well behaved?
#
# So it runs 07's panel B and nothing else, at one radius, for the face - house
# contrast, on each subject `08-prepare-subjects.R` prepared:
#
#   * the global budget: the contrast total over whole VT under
#     `whole_brain("none")`, the UNNORMALIZED comparator (contract section 2,
#     precondition 1 -- the default "local" would manufacture a spurious
#     failure);
#   * the attribution map: conservative searchlights at PAIR_RADIUS_MM, whose
#     per-node totals must sum to that budget;
#   * the receipt `sum_x total_x - whole_VT_total`, accepted at 1e-10;
#   * the aggregated coherent share at that radius, via `contribution()` over
#     the whole mask -- the alpha-invariant quantity of contract 3.2, read here
#     at a single scale rather than across a family;
#   * Haxby's own functional ROIs as territories, when the subject's
#     distribution ships them.
#
# What this script is NOT: it is not a group analysis. Six subjects analyzed
# separately with no group model, no permutation test and no p-value are six
# separate point estimates. Conservation is a per-subject algebraic identity
# and holding six times is six checks of the algebra, not evidence about a
# population (contract 7.6a: a conserved ledger buys no error bars). The
# scientific spread across subjects is reported as description only.
#
# Idempotent and resumable: subjects whose prepared object is absent are
# skipped and named, never invented. Runtime is a few seconds per subject.
#
# Environment: SUBJECTS, EXEMPLAR_DIR.
#
# Output: results/six-subject-conservative-receipts.csv  (committed evidence),
#         results/six-subject-conservative.rds,
#         results/six-subject-conservative.png.

suppressMessages(library(neuroim2))

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
source(file.path(exemplar_dir, "models.R"))
source(file.path(exemplar_dir, "subjects.R"))
crossform_source <- load_crossform(exemplar_dir)
paths <- exemplar_paths(exemplar_dir)
dir.create(paths$results, showWarnings = FALSE, recursive = TRUE)

## ---- Settings, matched to 07 --------------------------------------------
PAIR_RADIUS_MM <- 8                # 07's panel A/B radius
CONSERVATION_TOLERANCE <- 1e-10    # 07's acceptance
CONTRAST_LABEL <- "face - house"

named_contrast <- function(positive, negative) {
  w <- stats::setNames(numeric(length(CONDITIONS)), CONDITIONS)
  w[positive] <- 1 / length(positive)
  w[negative] <- -1 / length(negative)
  w
}
CONTRAST <- named_contrast("face", "house")
stopifnot(abs(sum(CONTRAST)) < 1e-12)

requested <- Sys.getenv("SUBJECTS")
requested <- if (nzchar(requested)) {
  trimws(strsplit(requested, "[,[:space:]]+")[[1L]])
} else {
  SUBJECT_IDS
}
requested <- requested[nzchar(requested)]

subjects <- Filter(function(s)
  file.exists(subject_paths(exemplar_dir, s)$prepared), requested)
missing_prep <- setdiff(requested, subjects)
if (length(missing_prep)) {
  message("No prepared object for: ", paste(missing_prep, collapse = ", "),
          "\n  (run 00-download.R then 08-prepare-subjects.R)")
}
if (!length(subjects)) {
  stop("Nothing prepared. Run 00-download.R and 08-prepare-subjects.R first.")
}
message("Subjects in scope: ", paste(subjects, collapse = ", "))

script_t0 <- Sys.time()

## ---- Receipts ------------------------------------------------------------
# One flat table, `subject` as the leading key. `tolerance` is finite only on
# rows that assert an identity, and on those rows `passes` is exactly
# `abs(value) <= tolerance`. Descriptive rows carry NA for both.
receipts <- list()
record <- function(subject, quantity, value, tolerance = NA_real_, note = "") {
  passes <- if (is.na(tolerance)) NA else abs(value) <= tolerance
  receipts[[length(receipts) + 1L]] <<- data.frame(
    subject = subject, quantity = quantity, value = as.numeric(value),
    tolerance = as.numeric(tolerance), passes = passes, note = note,
    stringsAsFactors = FALSE
  )
  invisible(value)
}
node_centres <- function(view) {
  index <- view$index
  if (is.data.frame(index)) as.integer(index$measurement) else as.integer(index)
}

## ---- Per subject ---------------------------------------------------------
per_subject <- list()
failures <- list()

for (subj in subjects) {
  message("\n== ", subj, " ==")
  sp <- subject_paths(exemplar_dir, subj)
  subject_t0 <- Sys.time()

  out <- tryCatch({
    prep <- readRDS(sp$prepared)
    stopifnot(identical(prep$conditions, CONDITIONS))

    ## Domain: the VT mask exactly as the preparation left it, degenerate
    ## voxels removed. Building it from the mask file plus the recorded
    ## exclusion keeps the frame's feature order and the pattern columns in
    ## lockstep; the stopifnot below is what proves it rather than assumes it.
    vt_mask <- vt_mask_volume(sp$mask, drop = prep$dropped_index)
    domain <- neuroim2_volume_domain(vt_mask,
                                     id = paste0("haxby-", subj, "-vt"))
    stopifnot(identical(as.integer(domain$feature_ids),
                        as.integer(prep$vt_index)))
    n_voxels <- domain$n_features

    effects <- effect_space(CONDITIONS, basis_id = "haxby-condition-means:v1",
                            units = "within-run-z")
    rel <- relation(prep$runs, effects = effects, domain = domain)
    over <- cross_partitions(rel, independence = "independent",
                            generalizes_over = "run")
    stopifnot(nrow(over) == choose(length(prep$runs), 2L))

    ## ---- The global budget: whole VT, UNNORMALIZED ------------------------
    whole_frame <- compile_frame(whole_brain("none"), domain)
    stopifnot(as.numeric(Matrix::rowSums(whole_frame$weights)) == n_voxels)
    whole_plan <- plan_geometry(rel, at = whole_frame, over = over)
    whole_view <- contrast_energy(whole_plan, CONTRAST)
    whole_total <- whole_view$total[[1L]]

    ## ---- The attribution map ----------------------------------------------
    cons_frame <- neuroim2_searchlights(vt_mask, radius = PAIR_RADIUS_MM,
                                        domain = domain,
                                        normalization = "conservative")
    cons_conservation <- frame_conservation(cons_frame)
    stopifnot(cons_conservation$conserved)
    cons_plan <- plan_geometry(rel, at = cons_frame, over = over)
    cons_view <- contrast_energy(cons_plan, CONTRAST)
    n_nodes <- length(cons_view$total)
    node_voxels <- as.numeric(Matrix::rowSums(cons_frame$weights > 0))

    identity_value <- sum(cons_view$total) - whole_total

    ## ---- The aggregated ledger over the whole mask ------------------------
    # `contribution()` is the package's own aggregation, so the coherent share
    # is read the way the contract says to read it rather than recomputed here.
    whole_ledger <- as.data.frame(contribution(cons_view,
                                               by = rep("VT", n_nodes)))
    coherent_share <- whole_ledger$coherence_fraction[[1L]]
    ledger_closes <- sum(whole_ledger$total) - sum(cons_view$total)

    ## ---- Haxby's functional ROIs, when this subject has them --------------
    roi_files <- c(`face-selective` = "mask8_face_vt.nii.gz",
                   `house-selective` = "mask8_house_vt.nii.gz")
    have_rois <- all(file.exists(file.path(sp$subj, roi_files)))
    territory_ledger <- NULL
    territory_counts <- NULL
    if (have_rois) {
      roi_index <- function(f) which(as.array(
        neuroim2::read_vol(file.path(sp$subj, f))) > 0)
      face_roi <- roi_index(roi_files[["face-selective"]])
      house_roi <- roi_index(roi_files[["house-selective"]])
      overlap <- length(intersect(face_roi, house_roi))
      centres <- node_centres(cons_view)
      territory <- ifelse(centres %in% face_roi, "face-selective",
                   ifelse(centres %in% house_roi, "house-selective", "other VT"))
      territory <- factor(territory, levels = c("face-selective",
                                                "house-selective", "other VT"))
      territory_counts <- stats::setNames(as.integer(tabulate(territory)),
                                          levels(territory))
      territory_ledger <- as.data.frame(contribution(cons_view,
                                                     by = territory))
      territory_ledger$overlap_voxels <- overlap
    }

    ## ---- Quirks ------------------------------------------------------------
    prov_path <- file.path(sp$data, sprintf("download-provenance-%s.rds", subj))
    download_quirk <- if (file.exists(prov_path)) {
      readRDS(prov_path)$quirk
    } else {
      ""
    }
    quirks <- c(prep$quirks,
                if (nzchar(download_quirk)) download_quirk else NULL,
                if (!have_rois)
                  "no mask8_face_vt / mask8_house_vt: territory ledger not read"
                else NULL)

    list(prep = prep, n_voxels = n_voxels, n_nodes = n_nodes,
         n_run_pairs = nrow(over), node_voxels = node_voxels,
         whole_total = whole_total, cons_total = sum(cons_view$total),
         identity = identity_value,
         column_mass_max_deviation = cons_conservation$max_deviation,
         coherent_share = coherent_share, ledger_closes = ledger_closes,
         node_totals = cons_view$total, node_centres = node_centres(cons_view),
         node_coherence_fraction = cons_view$coherence_fraction,
         whole_ledger = whole_ledger, territory_ledger = territory_ledger,
         territory_counts = territory_counts, have_rois = have_rois,
         quirks = quirks,
         scientific_plan_ids = list(whole = whole_plan$scientific_plan_id,
                                    conservative = cons_plan$scientific_plan_id),
         seconds = as.numeric(difftime(Sys.time(), subject_t0, units = "secs")))
  }, error = function(e) e)

  if (inherits(out, "error")) {
    message("  FAILED: ", conditionMessage(out))
    failures[[subj]] <- conditionMessage(out)
    next
  }
  per_subject[[subj]] <- out

  ## ---- Record ------------------------------------------------------------
  record(subj, "vt_voxels_in_mask", out$prep$n_vt_full,
         note = "voxels in mask4_vt.nii.gz as shipped")
  record(subj, "vt_voxels_analyzed", out$n_voxels,
         note = "domain features: mask voxels that vary in every run")
  record(subj, "vt_voxels_dropped", out$prep$n_vt_full - out$n_voxels,
         note = "constant in at least one run; z-scoring undefined, so excluded")
  record(subj, "n_runs", out$prep$n_runs,
         note = "runs contributing a full eight-condition mean matrix")
  record(subj, "n_run_pairs", out$n_run_pairs,
         note = "off-diagonal ordered run pairs the crossvalidated estimate averages")
  record(subj, "n_searchlights", out$n_nodes,
         note = sprintf("conservative searchlights, r = %g mm", PAIR_RADIUS_MM))
  record(subj, "median_searchlight_voxels", stats::median(out$node_voxels),
         note = "median node support at that radius")

  record(subj, "frame_column_mass_max_deviation",
         out$column_mass_max_deviation, tolerance = CONSERVATION_TOLERANCE,
         note = "max |sum_x w_xv - 1| over voxels: each voxel's evidence partitioned exactly once")
  record(subj, "whole_vt_total", out$whole_total,
         note = "the budget: face - house contrast total under whole_brain(\"none\")")
  record(subj, "conservative_total_sum", out$cons_total,
         note = "sum over the conservative searchlights")
  record(subj, "conservation_identity", out$identity,
         tolerance = CONSERVATION_TOLERANCE,
         note = "sum_x total_x - whole_VT_total (contract section 2, claim 2) -- THE acceptance")
  record(subj, "conservation_identity_relative",
         out$identity / out$whole_total, tolerance = CONSERVATION_TOLERANCE,
         note = "the same identity as a fraction of the budget")
  record(subj, "ledger_closes", out$ledger_closes,
         tolerance = CONSERVATION_TOLERANCE,
         note = "sum over contribution() groups minus sum over nodes; contribution() only adds")

  record(subj, "coherent_share_base_radius", out$coherent_share,
         note = sprintf("coherent share of the budget at r = %g mm, aggregated over VT by contribution(); frame-relative (contract 4)",
                        PAIR_RADIUS_MM))
  record(subj, "attribution_map_negative_nodes", sum(out$node_totals < 0),
         note = "signed estimation layer: crossvalidated node totals may be negative (contract 6)")
  record(subj, "half_budget_nodes",
         which(cumsum(sort(out$node_totals, decreasing = TRUE)) /
                 out$whole_total >= 0.5)[1L],
         note = "top-ranked nodes holding half the budget; a reading only an attribution map admits")

  if (out$have_rois) {
    for (i in seq_len(nrow(out$territory_ledger))) {
      row <- out$territory_ledger[i, ]
      record(subj, paste0("territory_budget_share: ", row$measurement),
             row$total / out$whole_total,
             note = sprintf("%d searchlights centred in Haxby's %s ROI", row$n_rows,
                            as.character(row$measurement)))
    }
    record(subj, "territory_ledger_closes",
           sum(out$territory_ledger$total) - sum(out$node_totals),
           tolerance = CONSERVATION_TOLERANCE,
           note = "sum over face/house/other territories minus sum over nodes")
  }

  record(subj, "n_quirks", length(out$quirks),
         note = if (length(out$quirks)) paste(out$quirks, collapse = "; ")
                else "none: 12 runs, all 8 conditions, no degenerate voxels, all masks present")

  message(sprintf("  VT %d voxels (%d dropped) | %d nodes | budget %.6f | sum %.6f",
                  out$n_voxels, out$prep$n_vt_full - out$n_voxels,
                  out$n_nodes, out$whole_total, out$cons_total))
  message(sprintf("  identity %s  [%s 1e-10] | coherent share %.4f | %.1f s",
                  format(out$identity, digits = 3),
                  if (abs(out$identity) <= CONSERVATION_TOLERANCE) "<=" else "FAILS >",
                  out$coherent_share, out$seconds))
  if (length(out$quirks)) message("  quirk: ", paste(out$quirks, collapse = "; "))
}

if (!length(per_subject)) stop("No subject completed; nothing to write.")

## ---- Receipts CSV --------------------------------------------------------
receipts <- do.call(rbind, receipts)
rownames(receipts) <- NULL
tolerance_rows <- !is.na(receipts$tolerance)
stopifnot(identical(receipts$passes[tolerance_rows],
                    abs(receipts$value[tolerance_rows]) <=
                      receipts$tolerance[tolerance_rows]))
if (!all(receipts$passes[tolerance_rows])) {
  print(receipts[tolerance_rows & !receipts$passes, ])
  stop("a recorded identity is outside its tolerance")
}
# A run deliberately narrowed with SUBJECTS is a debugging run, and must not
# quietly replace the committed six-subject evidence with a one-subject table.
# A run that asked for all six and got fewer -- a download that failed -- DOES
# write the canonical file, with only the subjects that completed in it: that
# is the honest record, and the ratchet test is supposed to notice.
narrowed <- !setequal(requested, SUBJECT_IDS)
receipts_path <- file.path(
  paths$results,
  if (narrowed) "six-subject-conservative-receipts-partial.csv"
  else "six-subject-conservative-receipts.csv")
utils::write.csv(receipts, receipts_path, row.names = FALSE)
message("\nWrote ", receipts_path, " (", nrow(receipts), " receipts over ",
        length(per_subject), " subjects, ", sum(tolerance_rows),
        " of them identities, all within tolerance)")
if (narrowed) {
  message("  SUBJECTS was narrowed to ", paste(requested, collapse = ", "),
          ", so the committed six-subject receipts were NOT overwritten.\n",
          "  Rerun without SUBJECTS to refresh them.")
}

## ---- Figure --------------------------------------------------------------
suffix <- if (narrowed) "-partial" else ""
png_path <- file.path(paths$results,
                      sprintf("six-subject-conservative%s.png", suffix))
grDevices::png(png_path, width = 1500, height = 620, res = 130)
op <- graphics::par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.4, 1.2),
                    mgp = c(2.7, 0.8, 0), cex.main = 1.0, font.main = 1)
palette6 <- c("#2166ac", "#4393c3", "#92c5de", "#d6604d", "#b2182b", "#762a83")
cols <- stats::setNames(palette6[seq_along(per_subject)], names(per_subject))

## Left: every subject's ledger closes at exactly 1.
plot(NA, xlim = c(0, max(vapply(per_subject, function(p) p$n_nodes, numeric(1)))),
     ylim = c(0, 1.05), xlab = "VT searchlights, ranked by their share",
     ylab = "cumulative share of that subject's whole-VT budget",
     main = sprintf("Attribution ledger, conservative r = %g mm", PAIR_RADIUS_MM))
graphics::abline(h = 1, lty = 2, col = "#9a9a9a")
for (subj in names(per_subject)) {
  p <- per_subject[[subj]]
  cum <- cumsum(sort(p$node_totals, decreasing = TRUE)) / p$whole_total
  graphics::lines(seq_along(cum), cum, col = cols[[subj]], lwd = 2)
}
graphics::legend("bottomright", bty = "n", cex = 0.72,
  legend = sprintf("%s  (%d nodes)", names(per_subject),
                   vapply(per_subject, function(p) p$n_nodes, numeric(1))),
  col = cols, lwd = 2)
graphics::mtext(
  sprintf("every curve closes at 1: max |sum_x total_x - whole VT| = %.1e over %d subjects",
          max(abs(vapply(per_subject, function(p) p$identity, numeric(1)))),
          length(per_subject)),
  side = 3, line = 0.25, cex = 0.62, col = "#b2182b")

## Right: the alpha-invariant share, per subject, with VT size alongside.
shares <- vapply(per_subject, function(p) p$coherent_share, numeric(1))
sizes <- vapply(per_subject, function(p) p$n_voxels, numeric(1))
bp <- graphics::barplot(shares, col = cols, border = NA, ylim = c(0, 1),
  ylab = "coherent share of the budget", names.arg = names(per_subject),
  main = sprintf("Coherent share at r = %g mm", PAIR_RADIUS_MM))
graphics::text(bp, shares, sprintf("%.3f", shares), pos = 3, cex = 0.72)
graphics::text(bp, 0.04, sprintf("%d vox", sizes), cex = 0.66, col = "white")
graphics::mtext(
  "frame-relative and descriptive: six separate point estimates, no group model, no inference",
  side = 3, line = 0.25, cex = 0.62, col = "#b2182b")
graphics::par(op)
invisible(grDevices::dev.off())
message("Wrote ", png_path)

## ---- Save ----------------------------------------------------------------
script_seconds <- as.numeric(difftime(Sys.time(), script_t0, units = "secs"))
identities <- vapply(per_subject, function(p) p$identity, numeric(1))
out <- list(
  arm = "crossform conservative geometry, six subjects (attribution only)",
  contract = "design/conservative-geometry-contract.md (conservative-geometry-v1)",
  estimand = paste(
    "crossvalidated face - house contrast energy over the off-diagonal run",
    "pairs, identity metric, read under conservative normalization at one",
    "radius, per subject, with the whole-VT total under whole_brain(\"none\")",
    "as the budget it must sum to"),
  caveat = paste(
    "Six subjects analyzed separately is not a group analysis. Conservation",
    "is a per-subject algebraic identity; holding six times checks the",
    "algebra six times and says nothing about a population. The spread of",
    "coherent shares across subjects is description, not a random effect",
    "(contract 7.6a: a conserved ledger buys no error bars)."),
  contrast = CONTRAST, contrast_label = CONTRAST_LABEL,
  radius_mm = PAIR_RADIUS_MM, tolerance = CONSERVATION_TOLERANCE,
  subjects = names(per_subject),
  requested = requested, missing_prepared = missing_prep, failures = failures,
  max_abs_conservation_identity = max(abs(identities)),
  receipts = receipts,
  per_subject = lapply(per_subject, function(p) {
    p$prep <- NULL
    p
  }),
  caveats = c("six subjects, analyzed one at a time; no group model",
              "block design, per-run condition means rather than GLM betas",
              "within-run z-scored time series, 2 TR haemodynamic shift",
              "no inference: no permutation test, no p-values, no error bars",
              "identity metric, native composition",
              "one radius; the multiscale reading and its alpha caveat are in 07"),
  session = list(crossform = as.character(utils::packageVersion("crossform")),
                 crossform_source = crossform_source,
                 neuroim2 = as.character(utils::packageVersion("neuroim2")),
                 when = Sys.time()),
  seconds = script_seconds
)
rds_path <- file.path(paths$results,
                      sprintf("six-subject-conservative%s.rds", suffix))
saveRDS(out, rds_path)
message("Wrote ", rds_path, " (", round(file.size(rds_path) / 1024), " KB)")

message("\n== six-subject summary ==")
print(data.frame(
  subject = names(per_subject),
  vt = vapply(per_subject, function(p) p$n_voxels, numeric(1)),
  nodes = vapply(per_subject, function(p) p$n_nodes, numeric(1)),
  budget = round(vapply(per_subject, function(p) p$whole_total, numeric(1)), 4),
  identity = signif(identities, 3),
  coherent_share = round(vapply(per_subject, function(p) p$coherent_share,
                                numeric(1)), 4),
  row.names = NULL), row.names = FALSE)
message("max |conservation identity| = ", format(max(abs(identities)), digits = 3),
        "  (acceptance ", CONSERVATION_TOLERANCE, ")")
if (length(failures)) {
  message("FAILED subjects: ", paste(names(failures), collapse = ", "))
}
message("Total: ", round(script_seconds, 1), " s")
