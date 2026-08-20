#!/usr/bin/env Rscript
# 07-conservative-geometry.R -- attribution, not detection.
#
# Scripts 02-04 established parity. 05 added the error channel and 06 ran the
# coherent/configuration decomposition on real data. Every one of them read a
# *detection* map: `normalization = "local"` searchlights, whose node value is
# a weighted average over the node's own support and therefore an intensive
# quantity, comparable across nodes and meaningless to add
# (`design/conservative-geometry-contract.md` section 1.1).
#
# This script runs the other instrument. A conservative frame is column
# normalized: every voxel's unit of evidence is partitioned among the nodes
# that see it, so a node value is its *share of one fixed global budget* --
# extensive, addable over territories, and not comparable node-for-node with a
# detection map (contract section 1.2). It draws four panels on Haxby subject
# 1's VT mask, for face - house and animate - inanimate:
#
#   A  the LOCAL detection map at one radius. What 02 and 06 already show,
#      relabelled as what it is: peaks, no budget. Its node values do not sum
#      to anything, and `contribution()` refuses it.
#
#   B  the CONSERVATIVE attribution map at the same radius and the same
#      neighbourhoods -- only the normalization differs. Its per-node totals
#      sum EXACTLY to the whole-VT total under `whole_brain("none")`, the
#      unnormalized global comparator (contract section 2, claim 2). The
#      identity `sum_x total_x - whole_VT_total` is recorded, with an
#      acceptance of 1e-10. `contribution()` then adds the ledger up over
#      Haxby's own functional ROIs (mask8_face_vt, mask8_house_vt).
#
#   C  the LATENT layer at whole VT. Crossvalidated estimates are signed, so
#      cumulative-contribution curves, effective counts and any other
#      functional that needs a nonnegative partition live only on a declared
#      PSD projection of them, labelled as such and reporting the mass it
#      moved (contract section 6). `latent_geometry()` is that projection.
#
#   D  the COHERENCE SPECTRUM over radii {4, 8, 12} mm, via a multiscale
#      conservative `frame_family()`. This is the panel the whole apparatus
#      exists for, and it is *not* a plot of energy against scale: see the
#      caveat below.
#
# THE SMOOTHED-LEDGER CAVEAT (contract section 3, and 3.1).
#   Under an identity metric the conservative total is
#   `G_x = sum_v w_xv b_v^L b_v^R'`: the per-voxel outer products form a fixed
#   ledger that the frame merely REALLOCATES. The total field at any scale is
#   a spatially smoothed version of one univariate per-voxel map and carries
#   no cross-voxel information at all. Worse, under a family whose members are
#   each column normalized and whose weights sum to one, the rows of scale s
#   sum to exactly `alpha_s * G_Omega` -- the analyst's own weight times a
#   constant -- WHATEVER THE DATA SAY. A multiscale panel of total energy per
#   scale is therefore a plot of `alpha`, not a finding, and the contract
#   makes it normative that no such panel may be presented as evidence about
#   spatial scale. crossform now refuses to draw it: this script exercises
#   that refusal and records it as a receipt.
#
#   What is NOT fixed by alpha is the split of each scale's fixed budget into
#   coherent and configuration parts. The coherent SHARE is exactly invariant
#   to alpha (contract section 3.2), so it is the one scale-resolved quantity
#   a conservative family can report. Panel D plots that and nothing else.
#
# Runtime: about 4 s of analysis on the prepared 12-run x 8-condition x 577-
# voxel object, plus R startup and mask reading. Nothing is subsampled.
#
# Output: results/conservative-geometry.png,
#         results/conservative-geometry-receipts.csv  (the committed receipts),
#         results/conservative-geometry.rds.

suppressMessages(library(neuroim2))

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
source(file.path(exemplar_dir, "models.R"))
crossform_source <- load_crossform(exemplar_dir)
paths <- exemplar_paths(exemplar_dir)
dir.create(paths$results, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(paths$prepared)) {
  stop("Prepared data not found at ", paths$prepared,
       "\nRun 01-prepare-data.R first (00-download.R before it if the raw ",
       "subject directory is absent).")
}
prep <- readRDS(paths$prepared)
stopifnot(identical(prep$conditions, CONDITIONS))

script_t0 <- Sys.time()

## ---- Scales -------------------------------------------------------------
# One radius for the A/B pair, so the only thing that differs between the
# detection map and the attribution map is the normalization; and three for
# the family. 8 mm sits inside the family, which keeps B a member of D rather
# than a fourth unrelated analysis. 06 reads its detection map at
# RADIUS_MM = 11.25 mm; the instrument, not the radius, is what A demonstrates.
PAIR_RADIUS_MM <- 8
FAMILY_RADII_MM <- c(4, 8, 12)
CONSERVATION_TOLERANCE <- 1e-10

## ---- Contrasts (identical coding to 06) ---------------------------------
named_contrast <- function(positive, negative) {
  w <- stats::setNames(numeric(length(CONDITIONS)), CONDITIONS)
  w[positive] <- 1 / length(positive)
  w[negative] <- -1 / length(negative)
  w
}
animate <- ANIMATE
inanimate <- setdiff(CONDITIONS, c(ANIMATE, "scrambledpix"))
CONTRASTS <- list(
  `face - house` = named_contrast("face", "house"),
  `animate - inanimate` = named_contrast(animate, inanimate)
)
stopifnot(all(vapply(CONTRASTS, function(w) abs(sum(w)) < 1e-12, logical(1))))

## ---- Domain, relation, pairing (identical to 02 and 06) -----------------
# `prep$mask_file` is an absolute path recorded when the object was built. The
# repository has since moved, so resolve against this checkout's own data
# directory and fall back only if that is absent.
mask_file <- file.path(paths$subj, basename(prep$mask_file))
if (!file.exists(mask_file)) mask_file <- prep$mask_file
if (!file.exists(mask_file)) {
  stop("VT mask not found at ", mask_file,
       "\nRun 00-download.R to fetch the subject directory.")
}
vt_vol <- neuroim2::read_vol(mask_file)
vt_mask <- neuroim2::LogicalNeuroVol(as.array(vt_vol) > 0,
                                     neuroim2::space(vt_vol))
domain <- neuroim2_volume_domain(vt_mask, id = "haxby-subj1-vt")
stopifnot(identical(domain$feature_ids, prep$vt_index))
n_voxels <- domain$n_features

effects <- effect_space(CONDITIONS, basis_id = "haxby-condition-means:v1",
                        units = "within-run-z")
rel <- relation(prep$runs, effects = effects, domain = domain)
over <- cross_partitions(rel, independence = "independent",
                         generalizes_over = "run")
stopifnot(nrow(over) == choose(length(prep$runs), 2L))
message("VT voxels: ", n_voxels, " | run pairs: ", nrow(over))

## ---- Receipts ------------------------------------------------------------
# Every number this script claims is appended here and written out as a CSV.
# `tolerance` is finite only for rows that assert an identity; those rows'
# `passes` is exactly `abs(value) <= tolerance` and nothing else. Descriptive
# rows carry NA for both.
receipts <- list()
# A view's `$index` is a bare vector of measurement identifiers for a
# single-frame result and a table for a family; take the centres either way.
node_centres <- function(view) {
  index <- view$index
  if (is.data.frame(index)) as.integer(index$measurement) else as.integer(index)
}
record <- function(panel, quantity, value, contrast = NA_character_,
                   group = NA_character_, tolerance = NA_real_,
                   note = "") {
  passes <- if (is.na(tolerance)) NA else abs(value) <= tolerance
  receipts[[length(receipts) + 1L]] <<- data.frame(
    panel = panel, quantity = quantity, contrast = contrast, group = group,
    value = as.numeric(value), tolerance = as.numeric(tolerance),
    passes = passes, note = note, stringsAsFactors = FALSE
  )
  invisible(value)
}

## ---- The global comparator: whole VT, unnormalized -----------------------
# Contract section 2, precondition 1: the comparator against which a
# conservative frame is checked is `whole_brain("none")`, NOT the default
# "local", which divides by the feature count and would manufacture a
# spurious conservation failure.
whole_frame <- compile_frame(whole_brain("none"), domain)
stopifnot(as.numeric(Matrix::rowSums(whole_frame$weights)) == n_voxels)
whole_plan <- plan_geometry(rel, at = whole_frame, over = over)
whole_view <- lapply(CONTRASTS, function(w) contrast_energy(whole_plan, w))
whole_total <- vapply(whole_view, function(v) v$total[[1L]], numeric(1))
for (label in names(CONTRASTS)) {
  record("B", "whole_vt_total", whole_total[[label]], contrast = label,
    group = "whole VT",
    note = "the global budget: contrast total under whole_brain(\"none\")")
}
message("whole-VT total (whole_brain(\"none\")): ",
        paste(sprintf("%s = %.6f", names(whole_total), whole_total),
              collapse = " | "))

## =========================================================================
## Panel A -- the LOCAL detection map
## =========================================================================
local_frame <- neuroim2_searchlights(vt_mask, radius = PAIR_RADIUS_MM,
                                     domain = domain, normalization = "local")
local_plan <- plan_geometry(rel, at = local_frame, over = over)
t0 <- Sys.time()
local_view <- lapply(CONTRASTS, function(w) contrast_energy(local_plan, w))
local_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
n_nodes <- length(local_view[[1L]]$total)

# A detection map has no budget. `frame_conservation()` says so rather than
# being asked politely, and the deviation it reports is the per-feature mass
# the local normalization does not control.
local_conservation <- frame_conservation(local_frame)
stopifnot(!local_conservation$conserved)
record("A", "local_frame_conserved", 0,
  note = "frame_conservation(local)$conserved is FALSE: a detection map has no budget (contract 1.1)")
record("A", "local_frame_feature_mass_max_deviation",
  local_conservation$max_deviation,
  note = "max |sum_x w_xv - 1| for the local frame; overlapping nodes double-count")

for (label in names(CONTRASTS)) {
  v <- local_view[[label]]
  record("A", "local_map_sum", sum(v$total), contrast = label,
    note = "NOT an estimate of anything: summing a detection map double-counts shared voxels (contract 1.1)")
  record("A", "local_map_sum_over_whole_vt_total",
    sum(v$total) / whole_total[[label]], contrast = label,
    note = "how far the meaningless sum lands from the budget it is not")
  record("A", "local_map_peak_total", max(v$total), contrast = label,
    group = as.character(node_centres(v)[which.max(v$total)]),
    note = "peak of the detection map; the reading a detection map admits")
  record("A", "local_map_median_total", stats::median(v$total),
    contrast = label, note = "median node of the detection map (intensive)")
}

# The refusal is the point, not an accident: adding up a detection map is the
# error a conservative frame exists to prevent.
territory_probe <- rep("VT", n_nodes)
local_refusal <- catch_refusal(contribution(local_view[[1L]],
                                            by = territory_probe))
stopifnot(inherits(local_refusal, "effect_capability_refusal"),
          identical(local_refusal$capability, "conservative_frame"))
record("A", "contribution_on_local_refused", 1,
  group = local_refusal$capability,
  note = paste0("catch_refusal(contribution(local_view, by = ...))$capability == \"",
                local_refusal$capability, "\""))
message("A: detection map, ", n_nodes, " local searchlights at r = ",
        PAIR_RADIUS_MM, " mm, ", round(local_seconds, 2),
        " s; contribution() refused with capability \"",
        local_refusal$capability, "\"")

## =========================================================================
## Panel B -- the CONSERVATIVE attribution map, and its ledger
## =========================================================================
cons_frame <- neuroim2_searchlights(vt_mask, radius = PAIR_RADIUS_MM,
                                    domain = domain,
                                    normalization = "conservative")
cons_conservation <- frame_conservation(cons_frame)
stopifnot(cons_conservation$conserved)
record("B", "frame_column_mass_max_deviation",
  cons_conservation$max_deviation, tolerance = CONSERVATION_TOLERANCE,
  note = "max |sum_x w_xv - 1| over voxels; the frame partitions each voxel's evidence exactly once")

cons_plan <- plan_geometry(rel, at = cons_frame, over = over)
t0 <- Sys.time()
cons_view <- lapply(CONTRASTS, function(w) contrast_energy(cons_plan, w))
cons_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# THE identity. Claim 2 of the contract, on real data.
conservation_identity <- numeric(0)
for (label in names(CONTRASTS)) {
  v <- cons_view[[label]]
  identity <- sum(v$total) - whole_total[[label]]
  conservation_identity[[label]] <- identity
  record("B", "conservative_total_sum", sum(v$total), contrast = label,
    note = "sum over the 577 conservative searchlights")
  record("B", "conservation_identity", identity, contrast = label,
    tolerance = CONSERVATION_TOLERANCE,
    note = "sum_x total_x - whole_VT_total under whole_brain(\"none\") (contract 2, claim 2)")
  record("B", "conservation_identity_relative",
    identity / whole_total[[label]], contrast = label,
    tolerance = CONSERVATION_TOLERANCE,
    note = "the same identity as a fraction of the budget")
  record("B", "attribution_map_negative_nodes",
    sum(v$total < 0), contrast = label,
    note = "signed estimation layer: crossvalidated node totals may be negative (contract 6)")
  message("B: ", label, " -- sum_x total_x = ",
          format(sum(v$total), digits = 15), ", whole VT = ",
          format(whole_total[[label]], digits = 15), ", identity = ",
          format(identity, digits = 3),
          if (abs(identity) <= CONSERVATION_TOLERANCE) "  [<= 1e-10]" else
            "  [FAILS]")
}
stopifnot(all(abs(conservation_identity) <= CONSERVATION_TOLERANCE))

# The ledger reading. Haxby's own functionally defined VT ROIs are the
# territories; a searchlight belongs to the territory its CENTRE falls in,
# because `contribution()` groups by row and a node is not divisible.
roi_index <- function(file) {
  vol <- neuroim2::read_vol(file.path(paths$subj, file))
  which(as.array(vol) > 0)
}
face_roi <- roi_index("mask8_face_vt.nii.gz")
house_roi <- roi_index("mask8_house_vt.nii.gz")
stopifnot(!length(intersect(face_roi, house_roi)))
centres <- node_centres(cons_view[[1L]])
territory <- ifelse(centres %in% face_roi, "face-selective",
             ifelse(centres %in% house_roi, "house-selective", "other VT"))
territory <- factor(territory,
                    levels = c("face-selective", "house-selective", "other VT"))
message("territories (by searchlight centre): ",
        paste(sprintf("%s %d", levels(territory), tabulate(territory)),
              collapse = " | "))

ledgers <- list()
for (label in names(CONTRASTS)) {
  ledger <- contribution(cons_view[[label]], by = territory)
  ledgers[[label]] <- as.data.frame(ledger)
  record("B", "ledger_closes", sum(ledger$total) - sum(cons_view[[label]]$total),
    contrast = label, tolerance = CONSERVATION_TOLERANCE,
    note = "sum over territories minus sum over nodes; contribution() only adds")
  for (i in seq_len(nrow(ledgers[[label]]))) {
    row <- ledgers[[label]][i, ]
    record("B", "ledger_total", row$total, contrast = label,
      group = as.character(row$measurement),
      note = sprintf("%d searchlights centred in this territory", row$n_rows))
    record("B", "ledger_budget_share", row$total / whole_total[[label]],
      contrast = label, group = as.character(row$measurement),
      note = "this territory's share of the whole-VT budget")
    record("B", "ledger_coherence_fraction", row$coherence_fraction,
      contrast = label, group = as.character(row$measurement),
      note = "frame-relative: a share of THIS frame's coherent mass (contract 4)")
  }
}
message("B: attribution map, ", n_nodes, " conservative searchlights, ",
        round(cons_seconds, 2), " s")

## =========================================================================
## Panel C -- the latent PSD descriptive layer at whole VT
## =========================================================================
# Crossvalidated estimates are signed (contract section 6). Cumulative
# contribution curves, effective counts and "top-k explains X%" all need a
# nonnegative partition, so they live only on a declared projection, which
# records how much mass it moved. Read at whole VT, so the object projected is
# exactly the budget panel B partitions.
t0 <- Sys.time()
whole_geometry <- materialize_geometry(whole_plan)
signed_spectrum <- as.numeric(geometry_spectrum(whole_geometry)$values)
latent <- latent_geometry(whole_geometry)
latent_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cumulative <- as.numeric(latent$cumulative[1L, ])
n_roots <- length(signed_spectrum)

record("C", "signed_spectrum_negative_roots", sum(signed_spectrum < 0),
  group = "whole VT",
  note = "the signed layer keeps negative roots; it never truncates them")
record("C", "signed_spectrum_min_root", min(signed_spectrum),
  group = "whole VT", note = "most negative crossvalidated root at whole VT")
record("C", "latent_n_eff", latent$n_eff[[1L]], group = "whole VT",
  note = paste0("participation ratio of the PSD-projected spectrum over ",
                n_roots, " roots; latent layer only"))
record("C", "latent_moved_mass", latent$moved_mass[[1L]], group = "whole VT",
  note = "absolute mass the psd_projection removed (sum of |negative roots|)")
record("C", "latent_moved_share", latent$moved_share[[1L]], group = "whole VT",
  note = "moved mass as a fraction of the source spectrum's total absolute mass")
record("C", "latent_clipped_measurements", latent$projection$clipped,
  group = "whole VT",
  note = "the empirical form DOES clip: the projection is not a formality here")
for (k in seq_along(cumulative)) {
  record("C", "latent_cumulative_contribution", cumulative[[k]],
    group = sprintf("k=%d", k),
    note = "C(k) on the latent PSD layer; forbidden on the signed layer (contract 6)")
}
message("C: latent layer at whole VT -- n_eff = ",
        round(latent$n_eff[[1L]], 4), ", moved mass = ",
        round(latent$moved_mass[[1L]], 4), " (",
        round(100 * latent$moved_share[[1L]], 3), "% of absolute mass), ",
        latent$projection$clipped, " of ", latent$projection$measurements,
        " measurements clipped, ", round(latent_seconds, 2), " s")

## =========================================================================
## Panel D -- the coherence spectrum over scale
## =========================================================================
family <- neuroim2_searchlights(vt_mask, FAMILY_RADII_MM, domain = domain,
                                normalization = "conservative")
family_conservation <- frame_conservation(family)
stopifnot(family_conservation$conserved)
record("D", "family_column_mass_max_deviation",
  family_conservation$max_deviation, tolerance = CONSERVATION_TOLERANCE,
  note = "the stacked family conserves")
for (i in seq_len(nrow(family_conservation$members))) {
  member <- family_conservation$members[i, ]
  record("D", "family_member_column_mass_max_deviation", member$max_deviation,
    group = as.character(member$family), tolerance = CONSERVATION_TOLERANCE,
    note = "each member is column normalized ON ITS OWN; a conserving stack does not imply it (contract 3.1, precondition 1)")
}
node_voxels <- as.numeric(Matrix::rowSums(family$weights > 0))
for (scale_name in unique(family$index$family)) {
  keep <- family$index$family == scale_name
  record("D", "median_node_voxels", stats::median(node_voxels[keep]),
    group = scale_name, note = "median searchlight size at this radius")
}

family_plan <- plan_geometry(rel, at = family, over = over)
t0 <- Sys.time()
spectra <- lapply(CONTRASTS, function(w) as.data.frame(
  coherence_spectrum(family_plan, w)))
located <- lapply(CONTRASTS, function(w) as.data.frame(
  coherence_spectrum(family_plan, w, by_location = TRUE)))
spectrum_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

for (label in names(CONTRASTS)) {
  sp <- spectra[[label]]
  lo <- located[[label]]
  for (i in seq_len(nrow(sp))) {
    scale_name <- as.character(sp$family[i])
    # The alpha-fixed column. Recorded so the caveat is checkable rather than
    # asserted -- and explicitly NOT a finding.
    record("D", "per_scale_total", sp$total[i], contrast = label,
      group = scale_name,
      note = "FIXED BY ALPHA, not a finding: alpha_s times the whole-VT total whatever the data say (contract 3.1)")
    record("D", "per_scale_total_minus_alpha_times_whole",
      sp$total[i] - sp$alpha[i] * whole_total[[label]], contrast = label,
      group = scale_name, tolerance = CONSERVATION_TOLERANCE,
      note = "the block law of contract 3.1, measured")
    record("D", "family_alpha", sp$alpha[i], contrast = label,
      group = scale_name, note = "the analyst's own weight; equal by default")
    # The alpha-invariant column. This is the finding.
    record("D", "coherent_share", sp$coherence_fraction[i], contrast = label,
      group = scale_name,
      note = "THE scale-resolved quantity: exactly invariant to alpha (contract 3.2)")
    at_scale <- lo[as.character(lo$family) == scale_name, ]
    valid <- !is.na(at_scale$coherence_fraction)
    record("D", "node_coherent_share_median",
      stats::median(at_scale$coherence_fraction[valid]), contrast = label,
      group = scale_name,
      note = "median over nodes with a valid fraction; the location-wise spectrum")
    record("D", "node_coherent_share_n_valid", sum(valid), contrast = label,
      group = scale_name,
      note = "nodes whose components form a nonnegative partition; the rest are masked, never clamped")
  }
  message("D: ", label, " -- coherent share ",
          paste(sprintf("%s: %.4f", sp$family, sp$coherence_fraction),
                collapse = ", "),
          " | per-scale total ",
          paste(sprintf("%.4f", sp$total), collapse = ", "),
          " (identical by construction)")
}

## ---- alpha-invariance, demonstrated rather than cited --------------------
# Contract section 3.2 claims the coherent share is EXACTLY invariant to the
# family weights while the energy is exactly proportional to them. On a
# fixture that is a theorem; here it is measured. The same three radii under a
# lopsided alpha must give identical shares and different energies, which is
# the whole reason the spectrum may be reported without disclosing alpha.
SKEWED_ALPHA <- c(0.6, 0.3, 0.1)
skewed_family <- neuroim2_searchlights(vt_mask, FAMILY_RADII_MM,
                                       domain = domain,
                                       normalization = "conservative",
                                       weights = SKEWED_ALPHA)
skewed_plan <- plan_geometry(rel, at = skewed_family, over = over)
skewed_spectra <- list()
for (label in names(CONTRASTS)) {
  skewed <- as.data.frame(coherence_spectrum(skewed_plan, CONTRASTS[[label]]))
  skewed_spectra[[label]] <- skewed
  equal <- spectra[[label]]
  stopifnot(identical(as.character(skewed$family), as.character(equal$family)))
  record("D", "coherent_share_alpha_invariance",
    max(abs(skewed$coherence_fraction - equal$coherence_fraction)),
    contrast = label, tolerance = CONSERVATION_TOLERANCE,
    note = "max |share under alpha = (0.6, 0.3, 0.1) - share under equal alpha|; contract 3.2, measured not cited")
  record("D", "per_scale_total_alpha_dependence",
    max(abs(skewed$total - equal$total)), contrast = label,
    note = "the same comparison for the energy column, which moves with alpha exactly as 3.1 says")
  message("D: ", label, " -- alpha-invariance of the share: ",
          format(max(abs(skewed$coherence_fraction -
                           equal$coherence_fraction)), digits = 3),
          "; the energy column moves by ",
          format(max(abs(skewed$total - equal$total)), digits = 6))
}

# The refusal, exercised rather than described. `plot(which = "profile")`
# walks `total` along the group index; on a coherence spectrum that column is
# alpha times a constant, so the package refuses to draw the panel and names
# the reason. Run inside a null device so no file is written for a plot that
# never happens.
spectrum_view <- coherence_spectrum(family_plan, CONTRASTS[[1L]])
grDevices::pdf(NULL)
panel_refusal <- catch_refusal(plot(spectrum_view, which = "profile"))
invisible(grDevices::dev.off())
stopifnot(inherits(panel_refusal, "effect_capability_refusal"),
          identical(panel_refusal$capability, "scale_energy_panel"))
record("D", "scale_energy_panel_refused", 1,
  group = panel_refusal$capability,
  note = paste0("catch_refusal(plot(spectrum, which = \"profile\"))$capability == \"",
                panel_refusal$capability, "\"; reason: ",
                paste(panel_refusal$reasons, collapse = ", ")))
message("D: plot(spectrum, which = \"profile\") refused -- capability \"",
        panel_refusal$capability, "\"; reason: ",
        paste(panel_refusal$reasons, collapse = ", "))

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
receipts_path <- file.path(paths$results, "conservative-geometry-receipts.csv")
utils::write.csv(receipts, receipts_path, row.names = FALSE)
message("\nWrote ", receipts_path, " (", nrow(receipts), " receipts, ",
        sum(tolerance_rows), " of them identities, all within tolerance)")

## ---- Figure --------------------------------------------------------------
png_path <- file.path(paths$results, "conservative-geometry.png")
grDevices::png(png_path, width = 1500, height = 1150, res = 130)
op <- graphics::par(mfrow = c(2, 2), mar = c(4.4, 4.6, 3.6, 1.4),
                    mgp = c(2.7, 0.8, 0), cex.main = 1.0, font.main = 1)
blue <- "#2166ac"
rust <- "#b2182b"
grey <- "#9a9a9a"
label_contrast <- "face - house"

## A -- detection map: peaks, deliberately drawn as a ranked profile so the
## picture itself says "read the top of this, do not add it up".
a_values <- sort(local_view[[label_contrast]]$total, decreasing = TRUE)
plot(seq_along(a_values), a_values, type = "h", col = grey, lwd = 1,
     xlab = "VT searchlights, ranked by their own value",
     ylab = "mean evidence density (intensive)",
     main = sprintf("A  DETECTION  local searchlights, r = %g mm",
                    PAIR_RADIUS_MM))
graphics::abline(h = 0, col = "black", lwd = 1)
top20 <- seq_len(20L)
graphics::points(top20, a_values[top20], pch = 16, col = blue, cex = 0.7)
graphics::legend("topright", bty = "n", cex = 0.8,
  legend = c("top 20 peaks",
             sprintf("peak %.3f, median %.3f",
                     a_values[1L], stats::median(a_values))),
  pch = c(16, NA), col = c(blue, NA))
graphics::mtext(
  "a weighted average over each node's own support: comparable across nodes, NOT addable",
  side = 3, line = 0.25, cex = 0.62, col = rust)

## B -- attribution map: shares of one budget, drawn as the cumulative ledger
## the instrument is for, so the picture closes at exactly 1.
b_values <- sort(cons_view[[label_contrast]]$total, decreasing = TRUE)
b_cumulative <- cumsum(b_values) / whole_total[[label_contrast]]
plot(seq_along(b_cumulative), b_cumulative, type = "l", col = blue, lwd = 2,
     ylim = range(0, b_cumulative, 1.02),
     xlab = "VT searchlights, ranked by their share",
     ylab = "cumulative share of the whole-VT budget",
     main = sprintf("B  ATTRIBUTION  conservative, r = %g mm",
                    PAIR_RADIUS_MM))
graphics::abline(h = 1, lty = 2, col = grey)
half <- which(b_cumulative >= 0.5)[1L]
graphics::segments(half, 0, half, 0.5, lty = 3, col = rust)
graphics::segments(0, 0.5, half, 0.5, lty = 3, col = rust)
graphics::points(half, 0.5, pch = 16, col = rust, cex = 0.8)
graphics::text(half, 0.42, sprintf("half the budget\nin %d of %d nodes",
                                   half, n_nodes),
               col = rust, cex = 0.72, adj = c(-0.08, 1))
graphics::legend("bottomright", bty = "n", cex = 0.72, text.col = "black",
  legend = c(
    sprintf("sum over nodes  %.9f", sum(cons_view[[label_contrast]]$total)),
    sprintf("whole VT        %.9f", whole_total[[label_contrast]]),
    sprintf("identity        %.2e  (<= 1e-10)",
            conservation_identity[[label_contrast]])))
graphics::mtext(
  "each voxel's evidence partitioned among the nodes that see it: shares of one fixed budget",
  side = 3, line = 0.25, cex = 0.62, col = rust)

## C -- the latent layer.
plot(seq_along(cumulative), cumulative, type = "b", pch = 16, col = blue,
     lwd = 2, ylim = c(0, 1.04), xaxt = "n",
     xlab = "roots retained, k", ylab = "cumulative contribution  C(k)",
     main = "C  LATENT LAYER  whole VT, psd_projection")
graphics::axis(1, at = seq_along(cumulative))
graphics::abline(h = 1, lty = 2, col = grey)
graphics::abline(v = latent$n_eff[[1L]], lty = 3, col = rust, lwd = 2)
graphics::text(latent$n_eff[[1L]], 0.22,
               sprintf("n_eff = %.3f", latent$n_eff[[1L]]),
               col = rust, cex = 0.78, adj = c(-0.1, 0.5))
graphics::legend("bottomright", bty = "n", cex = 0.72,
  legend = c(
    sprintf("%d of %d signed roots negative", sum(signed_spectrum < 0),
            n_roots),
    sprintf("moved mass %.3f (%.2f%% of |mass|)",
            latent$moved_mass[[1L]], 100 * latent$moved_share[[1L]])))
graphics::mtext(
  "C(k) and n_eff exist only here: the signed layer admits no nonnegative partition (contract 6)",
  side = 3, line = 0.25, cex = 0.62, col = rust)

## D -- the coherence spectrum. Share only. The energy panel is refused.
d_x <- FAMILY_RADII_MM
plot(d_x, spectra[[1L]]$coherence_fraction, type = "b", pch = 16, col = blue,
     lwd = 2, ylim = c(0, 1), xlim = range(d_x) + c(-1, 1),
     xlab = "searchlight radius (mm)", ylab = "coherent share of the budget",
     main = "D  COHERENCE SPECTRUM  conservative family {4, 8, 12} mm")
graphics::lines(d_x, spectra[[2L]]$coherence_fraction, type = "b", pch = 17,
                col = rust, lwd = 2, lty = 2)
node_median <- vapply(seq_along(d_x), function(i) {
  at_scale <- located[[1L]]
  at_scale <- at_scale[as.character(at_scale$family) ==
                         as.character(spectra[[1L]]$family[i]), ]
  stats::median(at_scale$coherence_fraction, na.rm = TRUE)
}, numeric(1))
graphics::lines(d_x, node_median, type = "b", pch = 1, col = blue, lty = 3)
# The same two spectra recomputed under alpha = (0.6, 0.3, 0.1) land exactly
# on the equal-alpha lines: the share is invariant to the weighting, while the
# energy the refused panel would have drawn is not.
graphics::points(d_x, skewed_spectra[[1L]]$coherence_fraction, pch = 3,
                 col = "black", cex = 1.1, lwd = 1.4)
graphics::points(d_x, skewed_spectra[[2L]]$coherence_fraction, pch = 3,
                 col = "black", cex = 1.1, lwd = 1.4)
graphics::legend("bottomleft", bty = "n", cex = 0.7,
  legend = c("face - house (aggregated)", "animate - inanimate (aggregated)",
             "face - house (median node)",
             sprintf("both, under alpha = (%s): identical to %.0e",
                     paste(SKEWED_ALPHA, collapse = ", "),
                     max(abs(skewed_spectra[[1L]]$coherence_fraction -
                               spectra[[1L]]$coherence_fraction)))),
  pch = c(16, 17, 1, 3), lty = c(1, 2, 3, NA),
  col = c(blue, rust, blue, "black"))
graphics::mtext(
  sprintf("per-scale TOTAL is %.4f at every radius = alpha_s x %.4f: plot(which=\"profile\") refused",
          spectra[[1L]]$total[1L], whole_total[[label_contrast]]),
  side = 3, line = 0.25, cex = 0.62, col = rust)
graphics::par(op)
invisible(grDevices::dev.off())
message("Wrote ", png_path)

## ---- Save ----------------------------------------------------------------
script_seconds <- as.numeric(difftime(Sys.time(), script_t0, units = "secs"))
out <- list(
  arm = "crossform conservative geometry (attribution)",
  contract = "design/conservative-geometry-contract.md (conservative-geometry-v1)",
  estimand = paste(
    "crossvalidated contrast energy over 66 off-diagonal run pairs, identity",
    "metric; read once under local normalization (detection: an intensive",
    "weighted average per node) and once under conservative normalization",
    "(attribution: each node's share of the whole-VT budget), plus a",
    "multiscale conservative family at 4/8/12 mm read as a coherence",
    "spectrum"),
  caveat = paste(
    "SMOOTHED LEDGER (contract section 3): under an identity metric the",
    "conservative total at any scale is a spatially smoothed version of one",
    "fixed per-voxel ledger, and under a family it is alpha_s times the",
    "whole-domain total by construction. Per-scale energy is the analyst's",
    "alpha vector, not a finding; only the coherent share, which is exactly",
    "alpha-invariant, is scale-resolved evidence (contract section 3.2)."),
  contrasts = CONTRASTS,
  animate = animate, inanimate = inanimate,
  n_voxels = n_voxels, n_nodes = n_nodes, n_run_pairs = nrow(over),
  pair_radius_mm = PAIR_RADIUS_MM, family_radii_mm = FAMILY_RADII_MM,
  whole_total = whole_total,
  receipts = receipts,
  panel_a = list(
    frame = "neuroim2_searchlights(normalization = \"local\")",
    conserved = local_conservation$conserved,
    feature_mass_max_deviation = local_conservation$max_deviation,
    map_sum = vapply(local_view, function(v) sum(v$total), numeric(1)),
    values = lapply(local_view, function(v) {
      data.frame(center = node_centres(v), total = v$total,
                 coherent = v$coherent, configuration = v$configuration)
    }),
    contribution_refusal = list(capability = local_refusal$capability,
                                namespace = local_refusal$namespace,
                                message = conditionMessage(local_refusal))
  ),
  panel_b = list(
    frame = "neuroim2_searchlights(normalization = \"conservative\")",
    column_mass_max_deviation = cons_conservation$max_deviation,
    conservation_identity = conservation_identity,
    tolerance = CONSERVATION_TOLERANCE,
    territory = list(levels = levels(territory),
                     counts = as.integer(tabulate(territory)),
                     source = "Haxby mask8_face_vt / mask8_house_vt, assigned by searchlight centre"),
    ledgers = ledgers,
    values = lapply(cons_view, function(v) {
      data.frame(center = node_centres(v), total = v$total,
                 coherent = v$coherent, configuration = v$configuration,
                 coherence_fraction = v$coherence_fraction)
    })
  ),
  panel_c = list(
    at = "whole VT, whole_brain(\"none\")",
    signed_spectrum = signed_spectrum,
    method = latent$method, component = latent$component,
    n_eff = latent$n_eff[[1L]], moved_mass = latent$moved_mass[[1L]],
    moved_share = latent$moved_share[[1L]],
    cumulative = cumulative,
    projection = latent$projection
  ),
  panel_d = list(
    family = "neuroim2_searchlights(c(4, 8, 12), normalization = \"conservative\")",
    skewed_alpha = SKEWED_ALPHA,
    skewed_spectra = skewed_spectra,
    members = family_conservation$members,
    median_node_voxels = tapply(node_voxels, family$index$family,
                                stats::median),
    spectra = spectra,
    node_shares = lapply(located, function(df) {
      df[, c("scale", "center", "coherent", "configuration", "total",
             "coherence_fraction")]
    }),
    scale_energy_panel_refusal = list(capability = panel_refusal$capability,
                                      namespace = panel_refusal$namespace,
                                      reasons = panel_refusal$reasons,
                                      message = conditionMessage(panel_refusal))
  ),
  scientific_plan_ids = list(
    whole = whole_plan$scientific_plan_id,
    local = local_plan$scientific_plan_id,
    conservative = cons_plan$scientific_plan_id,
    family = family_plan$scientific_plan_id
  ),
  timing = list(local_seconds = local_seconds,
                conservative_seconds = cons_seconds,
                latent_seconds = latent_seconds,
                spectrum_seconds = spectrum_seconds,
                script_seconds = script_seconds),
  caveats = c("one subject (Haxby 2001 subject 1)",
              "block design, per-run condition means rather than GLM betas",
              "within-run z-scored time series, 2 TR haemodynamic shift",
              "no inference: no permutation test, no group model, no p-values",
              "identity metric, native composition",
              "conservation is a point-estimate law and buys no error bars (contract 7.6a)"),
  session = list(crossform = as.character(utils::packageVersion("crossform")),
                 crossform_source = crossform_source,
                 neuroim2 = as.character(utils::packageVersion("neuroim2")),
                 when = Sys.time())
)
rds_path <- file.path(paths$results, "conservative-geometry.rds")
saveRDS(out, rds_path)
message("Wrote ", rds_path, " (",
        round(file.size(rds_path) / 1024), " KB)")
message("Total script time: ", round(script_seconds, 1), " s")
