#!/usr/bin/env Rscript
# 06-coherent-configuration.R -- the coherent/configuration decomposition on
# real data.
#
# Scripts 02-04 established that crossform's cross-validated distances agree
# with an independent implementation. Parity is the floor. This script runs
# the claim that has, until now, only been demonstrated on a generated
# fixture: that one plan reports a contrast's crossvalidated energy split
# exactly into a coherent part (the searchlight's own weighted common spatial
# mode) and a configuration part (the orthogonal pattern remainder), with
#
#     total = coherent + configuration
#
# holding exactly at every measurement, and with the signed marginal retained
# alongside so the direction of the effect is still readable.
#
# It consumes the same prepared object as 02 (data/prepared-smoke.rds: 12 runs
# x 8 conditions x VT voxels), builds the same VT searchlight frame, and reads
# two contrasts:
#
#   face - house              the canonical Haxby category contrast
#   animate - inanimate       (face, cat) vs (house, chair, bottle, scissors,
#                             shoe); scrambledpix is a visual control, not an
#                             object category, and carries weight 0
#
# The same plan is then read once more at a single whole-VT region, so the
# three-number summary can be read at one scale without refitting anything.
#
# Caveats, restated in the report this writes: one subject, a block design,
# per-run condition means rather than GLM betas, and no inference of any kind.
# The numbers here describe this subject's VT; they are not a group result.
#
# Output: results/coherent-configuration.rds, results/coherent-configuration.png,
#         and a dated section appended to README.md by hand from the console
#         summary this prints.

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

## ---- Contrasts ----------------------------------------------------------
# Weights are named, so contrast_energy() aligns them to the relation's own
# effect order rather than to the order written here.
named_contrast <- function(positive, negative) {
  w <- stats::setNames(numeric(length(CONDITIONS)), CONDITIONS)
  w[positive] <- 1 / length(positive)
  w[negative] <- -1 / length(negative)
  w
}

contrast_face_house <- named_contrast("face", "house")

animate <- ANIMATE                                   # cat, face
inanimate <- setdiff(CONDITIONS, c(ANIMATE, "scrambledpix"))
contrast_animacy <- named_contrast(animate, inanimate)
stopifnot(abs(sum(contrast_face_house)) < 1e-12,
          abs(sum(contrast_animacy)) < 1e-12)
message("animate: ", paste(animate, collapse = ", "),
        " | inanimate: ", paste(inanimate, collapse = ", "),
        " | zero-weighted: scrambledpix")

## ---- Domain and searchlight frame (identical to 02) ---------------------
# `prep$mask_file` is an absolute path recorded when the object was built and
# goes stale when the checkout moves; resolve against this checkout first.
mask_file <- file.path(paths$subj, basename(prep$mask_file))
if (!file.exists(mask_file)) mask_file <- prep$mask_file
if (!file.exists(mask_file)) {
  stop("VT mask not found at ", mask_file,
       "; re-run 00-download.R and 01-prepare-data.R.")
}
vt_vol <- neuroim2::read_vol(mask_file)
vt_mask <- neuroim2::LogicalNeuroVol(as.array(vt_vol) > 0,
                                     neuroim2::space(vt_vol))
domain <- neuroim2_volume_domain(vt_mask, id = "haxby-subj1-vt")
stopifnot(identical(domain$feature_ids, prep$vt_index))

frame <- neuroim2_searchlights(vt_mask, radius = RADIUS_MM, domain = domain,
                               normalization = "local")
n_centers <- nrow(frame$weights)
message("frame: ", n_centers, " VT searchlights at r = ", RADIUS_MM, " mm")

## ---- Relation, pairing, plan -------------------------------------------
effects <- effect_space(CONDITIONS, basis_id = "haxby-condition-means:v1",
                        units = "within-run-z")
rel <- relation(prep$runs, effects = effects, domain = domain)
over <- cross_partitions(rel, independence = "independent",
                         generalizes_over = "run")
stopifnot(nrow(over) == choose(length(prep$runs), 2L),
          !any(over$left == over$right))
message("pairing: ", nrow(over), " unordered run pairs, cross-run only")

plan <- plan_geometry(rel, at = frame, over = over)

## ---- The decomposition at every VT searchlight --------------------------
read_contrast <- function(weights, label) {
  t0 <- Sys.time()
  view <- contrast_energy(plan, weights)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # The partition is exact by construction; check it rather than assert it.
  recompose <- max(abs(view$total - (view$coherent + view$configuration)))
  message(label, ": ", n_centers, " searchlights in ", round(secs, 2),
          " s, max |total - (coherent + configuration)| = ",
          format(recompose, digits = 3))
  if (recompose > 1e-12) {
    stop("coherent + configuration does not recompose the total for ", label)
  }
  list(view = view, seconds = secs, recompose_max_abs = recompose)
}

fh <- read_contrast(contrast_face_house, "face - house")
an <- read_contrast(contrast_animacy, "animate - inanimate")

summarize <- function(res, label) {
  v <- res$view
  peak <- which.max(v$total)
  valid <- v$coherence_fraction_valid
  share <- v$coherence_fraction[valid]
  out <- list(
    label = label,
    n_measurements = length(v$total),
    peak_position = peak,
    peak_center = as.integer(v$index[peak]),
    peak = c(signed = v$signed[peak], coherent = v$coherent[peak],
             configuration = v$configuration[peak], total = v$total[peak],
             coherence_fraction = v$coherence_fraction[peak]),
    n_valid_fraction = sum(valid),
    share_median = if (any(valid)) stats::median(share) else NA_real_,
    share_quartiles = if (any(valid)) {
      stats::quantile(share, c(0.25, 0.75), names = FALSE)
    } else c(NA_real_, NA_real_),
    frac_total_positive = mean(v$total > 0),
    frac_configuration_positive = mean(v$configuration > 0),
    frac_coherent_positive = mean(v$coherent > 0),
    median_total = stats::median(v$total),
    median_coherent = stats::median(v$coherent),
    median_configuration = stats::median(v$configuration),
    # The coherent part is the crossvalidated energy of the same regional-mean
    # difference that $signed reports, so it should track signed^2 while
    # sitting below it: signed^2 squares one average and keeps that average's
    # noise, whereas the coherent term multiplies estimates from different
    # runs and is bias-free around zero. Recording the ratio makes that
    # relationship checkable rather than asserted.
    coherent_vs_signed_squared_median = stats::median(v$coherent / v$signed^2),
    coherent_vs_signed_squared_cor = stats::cor(v$coherent, v$signed^2),
    seconds = res$seconds,
    recompose_max_abs = res$recompose_max_abs
  )
  message("\n== ", label, " ==")
  message("  peak searchlight (position ", out$peak_position,
          ", centre voxel ", out$peak_center, "): signed ",
          round(out$peak[["signed"]], 4), ", coherent ",
          round(out$peak[["coherent"]], 4), ", configuration ",
          round(out$peak[["configuration"]], 4), ", total ",
          round(out$peak[["total"]], 4), ", coherent share ",
          if (is.na(out$peak[["coherence_fraction"]])) "NA (not a nonnegative partition)"
          else round(out$peak[["coherence_fraction"]], 3))
  message("  coherent share over VT: median ", round(out$share_median, 3),
          " IQR [", round(out$share_quartiles[1], 3), ", ",
          round(out$share_quartiles[2], 3), "] over ", out$n_valid_fraction,
          " of ", out$n_measurements, " searchlights with a valid fraction")
  message("  total > 0 at ", round(100 * out$frac_total_positive, 1),
          "% of searchlights; configuration > 0 at ",
          round(100 * out$frac_configuration_positive, 1),
          "%; coherent > 0 at ", round(100 * out$frac_coherent_positive, 1), "%")
  message("  coherent vs signed^2: median ratio ",
          round(out$coherent_vs_signed_squared_median, 3), ", r = ",
          round(out$coherent_vs_signed_squared_cor, 3))
  out
}

sum_fh <- summarize(fh, "face - house")
sum_an <- summarize(an, "animate - inanimate")

## ---- The same plan read at one whole-VT region --------------------------
# Not a second analysis: the same relation and pairing, measured at a frame
# with a single row instead of 577.
region_labels <- rep("VT", length(domain$feature_ids))
region_frame <- compile_frame(regions(region_labels), domain)
region_plan <- plan_geometry(rel, at = region_frame, over = over)
t0 <- Sys.time()
region_fh <- contrast_energy(region_plan, contrast_face_house)
region_an <- contrast_energy(region_plan, contrast_animacy)
region_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

region_row <- function(v) {
  c(signed = v$signed[1], coherent = v$coherent[1],
    configuration = v$configuration[1], total = v$total[1],
    coherence_fraction = v$coherence_fraction[1])
}
region_table <- rbind(`face - house` = region_row(region_fh),
                      `animate - inanimate` = region_row(region_an))
message("\n== whole VT as one region (", length(domain$feature_ids),
        " voxels), ", round(region_seconds, 2), " s ==")
print(round(region_table, 5))

## ---- Figure -------------------------------------------------------------
top_k <- 20L
top_fh <- order(fh$view$total, decreasing = TRUE)[seq_len(top_k)]
png_path <- file.path(paths$results, "coherent-configuration.png")
grDevices::png(png_path, width = 1400, height = 640, res = 130)
plot(fh$view, highlight = top_fh,
     highlight_label = paste0("top ", top_k, " searchlights by total"),
     main = c("Haxby subj-1 VT: face - house, coherent/configuration plane",
              "Haxby subj-1 VT: face - house, total energy by searchlight"))
invisible(grDevices::dev.off())
message("\nWrote ", png_path)

## ---- Save ---------------------------------------------------------------
script_seconds <- as.numeric(difftime(Sys.time(), script_t0, units = "secs"))
out <- list(
  arm = "crossform coherent/configuration",
  estimand = paste("crossvalidated contrast energy over 66 off-diagonal run",
                   "pairs, identity metric, local frame normalisation,",
                   "split into the frame's weighted common spatial mode",
                   "(coherent) and its orthogonal remainder (configuration)"),
  contrasts = list(`face - house` = contrast_face_house,
                   `animate - inanimate` = contrast_animacy),
  animate = animate, inanimate = inanimate,
  centers = as.integer(fh$view$index),
  n_voxels = as.numeric(Matrix::rowSums(frame$weights > 0)),
  values = list(
    `face - house` = data.frame(
      center = as.integer(fh$view$index), signed = fh$view$signed,
      coherent = fh$view$coherent, configuration = fh$view$configuration,
      total = fh$view$total, coherence_fraction = fh$view$coherence_fraction,
      coherence_fraction_valid = fh$view$coherence_fraction_valid),
    `animate - inanimate` = data.frame(
      center = as.integer(an$view$index), signed = an$view$signed,
      coherent = an$view$coherent, configuration = an$view$configuration,
      total = an$view$total, coherence_fraction = an$view$coherence_fraction,
      coherence_fraction_valid = an$view$coherence_fraction_valid)
  ),
  summaries = list(`face - house` = sum_fh, `animate - inanimate` = sum_an),
  region = list(n_voxels = length(domain$feature_ids), table = region_table,
                seconds = region_seconds),
  radius_mm = RADIUS_MM,
  n_runs = length(prep$runs),
  scientific_plan_id = plan$scientific_plan_id,
  region_scientific_plan_id = region_plan$scientific_plan_id,
  timing = list(face_house_seconds = fh$seconds,
                animacy_seconds = an$seconds,
                region_seconds = region_seconds,
                script_seconds = script_seconds),
  caveats = c("one subject (Haxby 2001 subject 1)",
              "block design, per-run condition means rather than GLM betas",
              "within-run z-scored time series, 2 TR haemodynamic shift",
              "no inference: no permutation test, no group model, no p-values",
              "identity metric, not a whitened or learned metric"),
  session = list(crossform = as.character(utils::packageVersion("crossform")),
                 crossform_source = crossform_source,
                 neuroim2 = as.character(utils::packageVersion("neuroim2")),
                 when = Sys.time())
)
saveRDS(out, file.path(paths$results, "coherent-configuration.rds"))
message("Wrote ", file.path(paths$results, "coherent-configuration.rds"))
message("Total script time: ", round(script_seconds, 1), " s")
