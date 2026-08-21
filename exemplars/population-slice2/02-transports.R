#!/usr/bin/env Rscript
# 02-transports.R -- the group searchlight grid, P^A, P^F, and the six
# data-free diagnostics of population-form-v1 section 7.5.
#
# THE GEOMETRY. Every subject's functional data already lives on one shared
# 66 x 78 x 61 MNI lattice (DECISION.md check 13), so unlike slice 1 there IS
# voxel correspondence between subjects and a searchlight-level transport is
# meaningful. What differs between subjects is coverage: each subject's native
# territory is the intersection of their own seven run brain masks, and those
# masks disagree substantially.
#
#   native nodes   one conservative searchlight per covered voxel, radius
#                  8 mm, centred on that voxel. The frame is conservative in
#                  the sense conservative-geometry-v1 section 2 means: each
#                  voxel's mass is split among the searchlights containing it,
#                  so the column masses sum to one and the searchlight totals
#                  sum to the subject's whole-coverage budget.
#   group nodes    every third lattice voxel in each direction (8.92 x 8.92 x
#                  9.66 mm), restricted to the 12-subject CONSENSUS coverage
#                  and to grey matter (across-subject mean fMRIPrep GM
#                  probseg, box-averaged onto the functional grid, > 0.25).
#
# THE TWO TRANSPORTS ARE A CONTROLLED PAIR IN THE SINK.
#
#   P^A  hard nearest-centre assignment (crossform's `anatomical_transport()`),
#        mass 1 on the nearest group node within 9 mm, otherwise all sink.
#        Row entropy is identically zero and displacement is the distance to
#        the nearest node: that is what "anatomical" means here.
#   P^F  softmax over a WIDER 12 mm neighbourhood, weighted by the similarity
#        between the voxel's own functional fingerprint and each candidate
#        group node's fingerprint -- but restricted to the rows P^A placed, so
#        that its sink mass is P^A's sink mass, row by row and bit for bit.
#
# Section 7.4 measures an adversarial transport that "wins" on eta purely by
# sinking the territory that disagrees, at eta = +0.167 while discarding 83 %
# of the brain. Pinning P^F's sink to P^A's removes that degree of freedom
# entirely: the two transports carry exactly the same territory, and eta can
# only be responding to where the mass goes among the nodes they both reach.
#
# The 12 mm support is a real asymmetry -- P^F may redistribute where P^A may
# not, and spreading mass smooths, which can raise a consensus share on its
# own. At 9 mm the mean support is 1.8 nodes, which is too few destinations for
# a transport to express any preference at all, so the asymmetry buys the
# question its subject matter. It is also exactly what the section 7.3 null
# band controls: a permuted-fingerprint P^F has the same 12 mm support, the
# same softmax, the same sink and the same amount of smoothing, so the null
# band prices the smoothing and eta's rank within it prices the fingerprint.
# The six diagnostics are still all reported -- the control is an argument for
# believing eta, not a licence to stop measuring.
#
# THE FINGERPRINT AND ITS CROSS-FIT. Fingerprints come from task
# `sharedreward` -- a different task, different runs, no shared trials with
# `trust` -- using only the seven event-level conditions present in all 24 runs
# (DECISION.md risk 7). A voxel's fingerprint is its 7-vector of condition
# betas, centred and unit-normed. A group node's fingerprint is the mean of the
# OTHER ELEVEN subjects' fingerprints over the voxels within 6 mm of it:
# leave-one-subject-out, so that a subject's transport is never built from that
# subject's own data even within the fitting task. `provenance$cross_fit`
# records `"task-sharedreward"`, and 04 refuses to evaluate eta on it.
#
# The softmax temperature is FIXED at 6, not tuned. Tuning it against the
# held-out `trust` data would be exactly the circularity contract section 7.2
# measures at 3.15x the honest gain.
#
# Environment: SLICE2_DIR.
# Output: data/derived/transports.rds   (git-ignored; P^A, P^F, the support
#           and the fingerprints 04 needs to rebuild permuted nulls)
#         results/population-slice2-transport-diagnostics.csv  (committed)

source(file.path(
  if (nzchar(Sys.getenv("SLICE2_DIR"))) Sys.getenv("SLICE2_DIR") else
    normalizePath(dirname(sub("^--file=", "",
      grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))),
  "00-common.R"))
crossform_version <- load_crossform()
# `build_loo_atlas()`; also used by 05 to build the across-run atlas.
source(file.path(SLICE_DIR, "eta-common.R"))
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

GROUP_FINGERPRINT_RADIUS_MM <- 6
script_t0 <- Sys.time()

## ---- Load what 01 fitted --------------------------------------------------
present <- SUBJECTS[file.exists(vapply(SUBJECTS, path_betas, ""))]
if (length(present) < 2L) {
  stop("Need at least two fitted subjects; found ", length(present),
       ". Run ./fetch.sh then 01-glm.R.")
}
if (length(present) < length(SUBJECTS)) {
  say("NOTE: %d of %d subjects fitted; running on those and saying so.",
      length(present), length(SUBJECTS))
}
say("Subjects: %s", paste(present, collapse = ", "))

fits <- lapply(present, function(s) readRDS(path_betas(s)))
names(fits) <- present

ref_space <- reference_space()
DIM3    <- as.integer(dim(ref_space))
SPACING <- as.numeric(neuroim2::spacing(ref_space))
ORIGIN  <- as.numeric(neuroim2::origin(ref_space))
NVOX3   <- prod(DIM3)

# Coordinates in the domain's own frame: mm from the grid origin, which is what
# `neuroim2_volume_domain()` reports and therefore what `compile_frame()` and
# `anatomical_transport()` measure their radii in. MNI world coordinates are
# these plus ORIGIN; the affine is diagonal, so the two frames differ by a
# translation and nothing else.
lattice_mm <- function(vox) {
  cbind((vox[, 1L] - 1) * SPACING[1L],
        (vox[, 2L] - 1) * SPACING[2L],
        (vox[, 3L] - 1) * SPACING[3L])
}

## ---- The group grid -------------------------------------------------------
say("\n== group grid ==")

coverage_count <- integer(NVOX3)
for (s in present) {
  coverage_count[fits[[s]]$voxel_linear] <-
    coverage_count[fits[[s]]$voxel_linear] + 1L
}
consensus <- coverage_count == length(present)
any_cover <- coverage_count > 0L
say("  coverage: union %d, consensus (all %d) %d, %.1f%% of the union is",
    sum(any_cover), length(present), sum(consensus),
    100 * (1 - sum(consensus) / sum(any_cover)))
say("            covered for some subjects and not others")

# Grey matter: the across-subject mean of the 1 mm fMRIPrep GM probseg, box
# averaged onto the functional grid. Box averaging rather than interpolation
# because the target voxel is ~28x the volume of the source voxel; every anat
# voxel is assigned to the functional voxel whose centre is nearest, and the
# functional value is the mean over its members (risk 2: resample once, and
# record it).
gm_cache <- file.path(DERIVED_DIR, "gm-mean-func-grid.rds")
if (file.exists(gm_cache)) {
  gm_func <- readRDS(gm_cache)
} else {
  t0 <- Sys.time()
  anat_sum <- NULL
  anat_space <- NULL
  for (s in present) {
    v <- neuroim2::read_vol(path_gm_probseg(s))
    if (is.null(anat_sum)) {
      anat_sum <- as.array(v) * 0
      anat_space <- neuroim2::space(v)
    }
    anat_sum <- anat_sum + as.array(v)
  }
  anat_mean <- anat_sum / length(present)
  adim <- dim(anat_mean)
  aspc <- as.numeric(neuroim2::spacing(anat_space))
  aorg <- as.numeric(neuroim2::origin(anat_space))
  ai <- arrayInd(seq_len(prod(adim)), adim)
  world <- cbind((ai[, 1L] - 1) * aspc[1L] + aorg[1L],
                 (ai[, 2L] - 1) * aspc[2L] + aorg[2L],
                 (ai[, 3L] - 1) * aspc[3L] + aorg[3L])
  fi <- round((world[, 1L] - ORIGIN[1L]) / SPACING[1L]) + 1L
  fj <- round((world[, 2L] - ORIGIN[2L]) / SPACING[2L]) + 1L
  fk <- round((world[, 3L] - ORIGIN[3L]) / SPACING[3L]) + 1L
  ok <- fi >= 1L & fi <= DIM3[1L] & fj >= 1L & fj <= DIM3[2L] &
        fk >= 1L & fk <= DIM3[3L]
  flin <- fi[ok] + (fj[ok] - 1L) * DIM3[1L] +
          (fk[ok] - 1L) * DIM3[1L] * DIM3[2L]
  vals <- as.numeric(anat_mean)[ok]
  o <- order(flin); fl <- flin[o]; vv <- vals[o]
  ends   <- c(which(diff(fl) != 0L), length(fl))
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  cs <- cumsum(vv)
  sums <- cs[ends] - c(0, cs[utils::head(ends, -1L)])
  gm_func <- numeric(NVOX3)
  gm_func[fl[ends]] <- sums / (ends - starts + 1L)
  saveRDS(gm_func, gm_cache)
  say("  GM probseg resampled 1 mm -> functional grid in %.0fs", elapsed(t0))
}

grid_vox <- arrayInd(seq_len(NVOX3), DIM3)
on_grid <- (grid_vox[, 1L] %% GRID_STEP == 0L) &
           (grid_vox[, 2L] %% GRID_STEP == 0L) &
           (grid_vox[, 3L] %% GRID_STEP == 0L)
group_lin <- which(on_grid & consensus & gm_func > GM_THRESHOLD)
if (!length(group_lin)) stop("group grid is empty; loosen GM_THRESHOLD")
group_vox <- arrayInd(group_lin, DIM3)
group_mm  <- lattice_mm(group_vox)
n_group   <- length(group_lin)
group_nodes <- sprintf("g%05d", seq_len(n_group))
say("  group nodes: %d  (every %dth voxel, consensus & mean GM > %.2f)",
    n_group, GRID_STEP, GM_THRESHOLD)
say("  spacing: %.2f x %.2f x %.2f mm; support radius %g mm",
    GRID_STEP * SPACING[1L], GRID_STEP * SPACING[2L], GRID_STEP * SPACING[3L],
    TRANSPORT_RADIUS_MM)

# The group-node id at each lattice position, 0 where there is none. Lets the
# support be looked up by lattice arithmetic instead of an n x m distance
# matrix (64k x 2.5k would be 1.3 GB of distances to find ~4 neighbours each).
gid <- integer(NVOX3)
gid[group_lin] <- seq_len(n_group)

## ---- Neighbourhood offsets ------------------------------------------------
# Offsets are computed in MILLIMETRES from the anisotropic spacing, never in
# voxels (risk 3): 3 steps in x is 8.92 mm but 3 steps in z is 9.66 mm, so a
# voxel-radius ball would be squashed along z.
support_offsets <- function(radius_mm) {
  lim <- floor(radius_mm / SPACING)
  g <- expand.grid(di = seq(-lim[1L], lim[1L]),
                   dj = seq(-lim[2L], lim[2L]),
                   dk = seq(-lim[3L], lim[3L]))
  d <- sqrt((g$di * SPACING[1L])^2 + (g$dj * SPACING[2L])^2 +
            (g$dk * SPACING[3L])^2)
  keep <- d <= radius_mm
  cbind(as.matrix(g[keep, , drop = FALSE]), dist = d[keep])
}

# For each native voxel, every group node within `offsets` of it.
gather_support <- function(vox, offsets) {
  n <- nrow(vox)
  pieces <- vector("list", nrow(offsets))
  for (o in seq_len(nrow(offsets))) {
    ii <- vox[, 1L] + offsets[o, "di"]
    jj <- vox[, 2L] + offsets[o, "dj"]
    kk <- vox[, 3L] + offsets[o, "dk"]
    ok <- ii >= 1L & ii <= DIM3[1L] & jj >= 1L & jj <= DIM3[2L] &
          kk >= 1L & kk <= DIM3[3L]
    if (!any(ok)) next
    lin <- ii[ok] + (jj[ok] - 1L) * DIM3[1L] +
           (kk[ok] - 1L) * DIM3[1L] * DIM3[2L]
    g <- gid[lin]
    hit <- g > 0L
    if (!any(hit)) next
    pieces[[o]] <- cbind(row = which(ok)[hit], col = g[hit],
                         dist = offsets[o, "dist"])
  }
  s <- do.call(rbind, pieces)
  s[order(s[, "row"], s[, "col"]), , drop = FALSE]
}

SUPPORT_OFFSETS <- support_offsets(PF_SUPPORT_RADIUS_MM)
FP_OFFSETS      <- support_offsets(GROUP_FINGERPRINT_RADIUS_MM)
say("  P^F support offsets: %d lattice positions within %g mm (P^A assigns within %g mm)",
    nrow(SUPPORT_OFFSETS), PF_SUPPORT_RADIUS_MM, TRANSPORT_RADIUS_MM)

## ---- Per-subject domains, frames, fingerprints ----------------------------
say("\n== per-subject frames and fingerprints ==")
subj <- list()
for (s in present) {
  t0 <- Sys.time()
  f <- fits[[s]]
  mask_arr <- array(FALSE, DIM3)
  mask_arr[f$voxel_linear] <- TRUE
  lmask <- neuroim2::LogicalNeuroVol(mask_arr, ref_space)
  dom <- neuroim2_volume_domain(lmask, id = paste0("ds003745-", s))

  # The domain's feature order and this script's voxel order have to be the
  # same object or every row of the transport is mislabelled. `which()` is
  # column-major and so is the domain, but that is asserted against the
  # domain's own coordinates rather than believed.
  stopifnot(length(dom$feature_ids) == length(f$voxel_linear))
  expect_mm <- lattice_mm(f$voxel_index)
  stopifnot(max(abs(dom$coordinates - expect_mm)) < 1e-8)

  frame <- compile_frame(
    searchlights(radius = SEARCHLIGHT_RADIUS_MM, normalization = "conservative"),
    dom)
  cons <- frame_conservation(frame)
  stopifnot(cons$conserved, cons$max_deviation <= FRAME_TOLERANCE)
  native_nodes <- as.character(frame$index$measurement)
  stopifnot(identical(native_nodes, as.character(dom$feature_ids)))

  # Fingerprint: the 7 sharedreward condition betas, centred then unit-normed,
  # so that a dot product between two fingerprints IS their correlation.
  fp <- t(f$shared)                                  # V x 7
  fp <- fp - rowMeans(fp)
  nrm <- sqrt(rowSums(fp^2))
  degenerate_fp <- nrm <= .Machine$double.eps
  nrm[degenerate_fp] <- 1
  fp <- fp / nrm

  support <- gather_support(f$voxel_index, SUPPORT_OFFSETS)
  fp_support <- gather_support(f$voxel_index, FP_OFFSETS)

  subj[[s]] <- list(
    subject = s, group = f$group, domain = dom, frame = frame,
    voxel_index = f$voxel_index, coords = dom$coordinates,
    native_nodes = native_nodes, n_native = length(native_nodes),
    fingerprint = fp, degenerate_fp = degenerate_fp,
    support = support, fp_support = fp_support,
    frame_max_deviation = cons$max_deviation,
    mean_searchlight_voxels = length(frame$weights@x) / nrow(frame$index)
  )
  say("  %-9s V=%6d  searchlight mean %.1f voxels  support mean %.2f nodes  %.0fs",
      s, length(native_nodes), subj[[s]]$mean_searchlight_voxels,
      nrow(support) / length(native_nodes), elapsed(t0))
}

## ---- Group-node fingerprint atlas, leave one subject out ------------------
# Per subject and group node, the mean fingerprint over that subject's voxels
# within 6 mm; then the atlas a subject sees is the mean over the OTHER
# subjects only.
#
# `build_loo_atlas()` lives in eta-common.R because 05 builds a second atlas
# the same way from `trust` run-1 fingerprints. Sharing the code is not tidiness
# here: the README claims the two eta axes are the same estimator on different
# data, and two hand-written copies of the atlas would be the first place that
# claim quietly stopped being true.
say("\n== group fingerprint atlas (leave-one-subject-out) ==")
atlases <- build_loo_atlas(present,
                           lapply(subj, `[[`, "fingerprint"),
                           lapply(subj, `[[`, "fp_support"),
                           n_group)
say("  group nodes reached by every subject: %d of %d",
    sum(attr(atlases, "reach") == length(present)), n_group)

## ---- Build the two transports ---------------------------------------------
say("\n== transports ==")
group_index_shared <- NULL
transports_A <- list(); transports_F <- list()
diag_rows <- list(); pf_support <- list()

# Row entropy and the mass-weighted displacement, computed on the RENORMALIZED
# group part of each row exactly as contract section 7.5 defines them. Rows
# with no group mass are excluded from both and counted separately.
row_diagnostics <- function(P, native_mm, group_mm) {
  M <- P$matrix
  m <- nrow(P$group_index)
  G <- M[, seq_len(m), drop = FALSE]
  gmass <- Matrix::rowSums(G)
  live <- which(gmass > 0)
  Tg <- as(G[live, , drop = FALSE], "TsparseMatrix")
  r <- Tg@i + 1L; cc <- Tg@j + 1L; p <- Tg@x / gmass[live][Tg@i + 1L]
  # displacement: || center(x) - sum_j ptilde_xj center(j) ||
  bary <- matrix(0, length(live), 3L)
  for (d in 1:3) bary[, d] <- as.vector(rowsum(p * group_mm[cc, d], r,
                                               reorder = TRUE))
  disp <- sqrt(rowSums((native_mm[live, , drop = FALSE] - bary)^2))
  # entropy in nats over the renormalized group columns
  plp <- ifelse(p > 0, p * log(p), 0)
  H <- -as.vector(rowsum(plp, r, reorder = TRUE))
  list(displacement = disp, entropy = H, perplexity = exp(H),
       live = live, gmass = gmass,
       n_all_sink = length(gmass) - length(live))
}

for (s in present) {
  ss <- subj[[s]]
  n <- ss$n_native
  native_index <- data.frame(node = ss$native_nodes, stringsAsFactors = FALSE)

  ## P^A -- crossform's own anatomical constructor. Hard nearest centre with a
  ## radius; everything further than 9 mm from every group node is all sink.
  PA <- anatomical_transport(
    native_coords = ss$coords,
    group_coords  = group_mm,
    semantics     = "budget",
    radius        = TRANSPORT_RADIUS_MM,
    native_index  = native_index,
    group_index   = data.frame(node = group_nodes, stringsAsFactors = FALSE),
    provenance    = list(
      details = sprintf(
        paste0("nearest of %d consensus grey-matter group nodes within %g mm, ",
               "on the shared ds003745 MNI functional lattice; native nodes are ",
               "conservative %g mm searchlights centred on this subject's ",
               "seven-run coverage"),
        n_group, TRANSPORT_RADIUS_MM, SEARCHLIGHT_RADIUS_MM),
      dataset = "ds003745 snapshot 2.1.1",
      preprocessing = "fMRIPrep 21.0.2 MNI152NLin2009cAsym",
      grid_step_voxels = GRID_STEP, gm_threshold = GM_THRESHOLD))

  if (is.null(group_index_shared)) group_index_shared <- PA$group_index
  # `plan_population()` requires every subject's `$group_index` to be
  # `identical()`; reusing the first one verbatim is what guarantees it, and
  # also makes P^F's index identical to P^A's so the two runs are comparable
  # node for node.
  stopifnot(identical(PA$group_index, group_index_shared))

  ## P^F's support is the wider PF_SUPPORT_RADIUS_MM neighbourhood, restricted
  ## to the rows P^A actually placed. A row P^A sent entirely to the sink stays
  ## entirely in the sink under P^F -- that restriction is what makes the two
  ## sinks identical while letting P^F redistribute more widely among the rows
  ## it does carry.
  placed <- Matrix::rowSums(PA$matrix[, seq_len(n_group), drop = FALSE]) > 0
  sup <- ss$support[placed[ss$support[, "row"]], , drop = FALSE]
  la <- atlases[[s]]
  sim <- rowSums(ss$fingerprint[sup[, "row"], , drop = FALSE] *
                 la$atlas[sup[, "col"], , drop = FALSE])
  sim[!la$defined[sup[, "col"]]] <- 0
  sim[ss$degenerate_fp[sup[, "row"]]] <- 0
  sim_sd <- stats::sd(sim)

  ## P^F -- weights exp(lambda * similarity), normalized within the row. Each
  ## placed row therefore carries group mass exactly 1, which is what P^A gives
  ## it too, so the sink is bit-for-bit the same; the assertion below is what
  ## makes that a fact rather than an intention.
  w <- softmax_rows(sim, sup[, "row"], PF_TEMPERATURE)

  PFmat <- Matrix::sparseMatrix(i = sup[, "row"], j = sup[, "col"], x = w,
                                dims = c(n, n_group),
                                dimnames = list(ss$native_nodes, group_nodes))
  PF <- location_transport(
    matrix       = PFmat,
    native_index = native_index,
    group_index  = group_index_shared,
    semantics    = "budget",
    provenance   = list(
      method = "functional",
      details = sprintf(
        paste0("softmax(%g * r) over the group nodes within %g mm, where r is ",
               "the correlation between this voxel's %d-condition task-",
               "sharedreward fingerprint and the group node's fingerprint. The ",
               "fingerprint is read off a GLM in which every trial type in the ",
               "events file is modelled and only these %d are read; the group ",
               "fingerprint is the mean over the other %d subjects of their ",
               "fingerprints within %g mm of the node (leave-one-subject-out). ",
               "Restricted to the rows P^A placed within %g mm, so the sink is ",
               "P^A's sink row for row and the two transports differ only in ",
               "where mass goes among the nodes both reach. GLM model id: %s."),
        PF_TEMPERATURE, PF_SUPPORT_RADIUS_MM, length(SHARED_CONDITIONS),
        length(SHARED_CONDITIONS), length(present) - 1L,
        GROUP_FINGERPRINT_RADIUS_MM, TRANSPORT_RADIUS_MM, GLM_MODEL_ID),
      cross_fit = "task-sharedreward",
      fitted_on = "ds003745 task-sharedreward runs 1-2",
      evaluated_on = "ds003745 task-trust runs 1-5",
      temperature = PF_TEMPERATURE,
      leave_one_subject_out = TRUE))

  ## `plan_population()` refuses a set of transports whose group indices are
  ## not `identical()`. Both runs share one index object, which also means the
  ## two population fits are comparable node for node rather than merely
  ## node-count for node-count.
  stopifnot(identical(PF$group_index, group_index_shared))

  ## The controlled-pair assertion: identical sink, node for node.
  sinkA <- PA$matrix[, n_group + 1L]
  sinkF <- PF$matrix[, n_group + 1L]
  sink_gap <- max(abs(sinkA - sinkF))
  stopifnot(sink_gap <= POPULATION_TOLERANCE)

  dA <- row_diagnostics(PA, ss$coords, group_mm)
  dF <- row_diagnostics(PF, ss$coords, group_mm)
  terrA <- crossform:::.transport_sink_territory(PA)
  terrF <- crossform:::.transport_sink_territory(PF)

  ## A SEVENTH DIAGNOSTIC, beyond the six section 7.5 requires.
  ##
  ## The six catch a transport that wins eta by DISCARDING territory. They do
  ## not catch one that wins by CONCENTRATING it: a P^F that piles most of the
  ## brain onto a handful of group nodes has full subject coverage, bounded
  ## displacement, respectable entropy and P^A's exact sink, and would pass all
  ## six while making the group grid a fiction. So the arriving mass per group
  ## node is summarized by its inverse participation ratio -- (sum x)^2 /
  ## sum(x^2), the effective number of group nodes actually carrying the
  ## territory. Reported for both transports, because the comparison is the
  ## point: P^F must be allowed to concentrate somewhat (that is what choosing
  ## a destination means) and must not be allowed to concentrate onto nothing.
  mass_A <- Matrix::colSums(PA$matrix[, seq_len(n_group), drop = FALSE])
  mass_F <- Matrix::colSums(PF$matrix[, seq_len(n_group), drop = FALSE])
  stopifnot(abs(sum(mass_A) - sum(mass_F)) <= POPULATION_TOLERANCE * sum(mass_A))
  ipr <- function(x) {
    x <- x[x > 0]
    if (!length(x)) return(0)
    sum(x)^2 / sum(x^2)
  }

  transports_A[[s]] <- PA
  transports_F[[s]] <- PF
  pf_support[[s]] <- sup
  diag_rows[[s]] <- data.frame(
    subject = s, group = ss$group, n_native = n,
    n_group_nodes = n_group,
    n_placed_rows = sum(placed),
    support_mean_nodes = nrow(sup) / sum(placed),
    fingerprint_similarity_sd = sim_sd,
    sink_territory = terrA$share,
    sink_territory_PF = terrF$share,
    sink_identical_max_gap = sink_gap,
    all_sink_rows = dA$n_all_sink,
    displacement_median_A = stats::median(dA$displacement),
    displacement_p90_A = unname(stats::quantile(dA$displacement, 0.90)),
    displacement_max_A = max(dA$displacement),
    displacement_masswt_mean_A =
      sum(dA$displacement * dA$gmass[dA$live]) / sum(dA$gmass[dA$live]),
    displacement_median_F = stats::median(dF$displacement),
    displacement_p90_F = unname(stats::quantile(dF$displacement, 0.90)),
    displacement_max_F = max(dF$displacement),
    displacement_masswt_mean_F =
      sum(dF$displacement * dF$gmass[dF$live]) / sum(dF$gmass[dF$live]),
    group_node_mass_total = sum(mass_A),
    group_node_effective_A = ipr(mass_A),
    group_node_effective_F = ipr(mass_F),
    group_node_mass_max_A = max(mass_A),
    group_node_mass_max_F = max(mass_F),
    group_node_mass_min_A = min(mass_A),
    group_node_mass_min_F = min(mass_F),
    group_nodes_with_zero_mass_A = sum(mass_A <= 0),
    group_nodes_with_zero_mass_F = sum(mass_F <= 0),
    entropy_mean_A = mean(dA$entropy),
    entropy_mean_F = mean(dF$entropy),
    perplexity_mean_F = mean(dF$perplexity),
    exp_mean_entropy_F = exp(mean(dF$entropy)),
    stringsAsFactors = FALSE
  )
  say("  %-9s sink %.1f%%  sd(sim) %.3f  disp A/F %.2f/%.2f mm  H(F) %.3f  pplx %.2f  eff nodes A/F %.0f/%.0f of %d",
      s, 100 * terrA$share, sim_sd,
      stats::median(dA$displacement), stats::median(dF$displacement),
      mean(dF$entropy), mean(dF$perplexity),
      ipr(mass_A), ipr(mass_F), n_group)
}

diagnostics <- do.call(rbind, diag_rows)

## ---- Group-node subject coverage (data-free half of section 7.5) ----------
# The number of subjects putting nonzero mass on each group node, computed from
# the operators alone. `population_prevalence()` reports a data-dependent proxy
# for the same quantity in 03; both are recorded, and they answer different
# questions -- this one asks whether a subject can reach the node at all.
node_subject_count <- integer(n_group)
for (s in present) {
  reached <- Matrix::colSums(transports_A[[s]]$matrix[, seq_len(n_group),
                                                      drop = FALSE]) > 0
  node_subject_count <- node_subject_count + as.integer(reached)
}
say("\n== group-node subject coverage (from the operators) ==")
say("  min %d, median %d, max %d subjects per node; %d node(s) below the declared floor of %d",
    min(node_subject_count), stats::median(node_subject_count),
    max(node_subject_count), sum(node_subject_count < COVERAGE_FLOOR),
    COVERAGE_FLOOR)

## ---- Save -----------------------------------------------------------------
saveRDS(list(
  subjects = present,
  crossform_version = crossform_version,
  group = list(nodes = group_nodes, lin = group_lin, vox = group_vox,
               mm = group_mm, index = group_index_shared,
               subject_count = node_subject_count,
               consensus_voxels = sum(consensus), union_voxels = sum(any_cover)),
  settings = list(searchlight_radius_mm = SEARCHLIGHT_RADIUS_MM,
                  transport_radius_mm = TRANSPORT_RADIUS_MM,
                  grid_step = GRID_STEP, gm_threshold = GM_THRESHOLD,
                  pf_temperature = PF_TEMPERATURE,
                  group_fingerprint_radius_mm = GROUP_FINGERPRINT_RADIUS_MM),
  transports_A = transports_A,
  transports_F = transports_F,
  # Everything 04 needs to rebuild a permuted P^F without redoing the geometry.
  support = pf_support,
  fingerprint = lapply(subj, function(x) x$fingerprint),
  degenerate_fp = lapply(subj, function(x) x$degenerate_fp),
  atlas = unclass(atlases),
  native_nodes = lapply(subj, function(x) x$native_nodes),
  diagnostics = diagnostics
), file.path(DERIVED_DIR, "transports.rds"))

utils::write.csv(diagnostics,
                 file.path(RESULTS_DIR,
                           "population-slice2-transport-diagnostics.csv"),
                 row.names = FALSE)

say("\n02-transports.R done in %.1f min", elapsed(script_t0) / 60)
