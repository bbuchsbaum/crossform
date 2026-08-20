#!/usr/bin/env Rscript
# 10-population-slice1.R -- the first real population analysis: six subjects,
# region-level nodes, group form.
#
# `09-six-subject-conservative.R` runs six subjects one at a time and says so
# loudly: six point estimates are not a group analysis. This script is the
# thing 09 refused to be. It carries each subject's conservative geometry onto
# a shared set of group nodes with a declared transport, fits one group model
# over the six carried subjects, and reports the group form with the identity
# checks that make the carrying auditable.
#
# It is slice 1 of `design/population-form-contract.md` (population-form-v1):
# region-level nodes, a trivial identity transport, and only sections 1-6 of
# that contract. Section 7 -- `eta_transport`, the null band and the six
# transport diagnostics -- is slice 2's job and is deliberately absent here,
# because a transport that is the identity on named regions has nothing to
# diagnose: no displacement, no spreading, no sink. The two data-free
# diagnostics that ARE meaningful at region level (unmapped territory, group
# node subject coverage) are recorded, and they are recorded as the trivial
# numbers they are rather than dressed up.
#
# THE GEOMETRY. Each subject's VT mask is its own native space -- there is no
# voxel correspondence between Haxby subjects and this exemplar performs no
# registration (contract section 9.2 refuses it). So the nodes are REGIONS:
# each subject's analyzed VT is partitioned into three territories using the
# subject's own Haxby-distributed functional ROIs,
#
#     face-territory   mask8_face_vt   (intersected with the analyzed VT)
#     house-territory  mask8_house_vt  (minus any face-territory voxel)
#     other-VT         everything else in mask4_vt
#
# and the transport maps each native territory wholly onto the same-named
# group node. Three native nodes per subject, three group nodes, an identity
# 3x3 operator, an empty sink. The territories are a PARTITION of VT, so the
# region frame is conservative in exactly the sense
# `conservative-geometry-v1` section 2 means: the three region totals sum to
# the subject's whole-VT budget. That is the number 09 already reports, and
# this script reads 09's committed CSV and asserts the two agree rather than
# asserting it of itself.
#
# THE TERRITORIES ARE FUNCTIONALLY DEFINED, AND THAT IS A LIMITATION.
# `mask8_face_vt` and `mask8_house_vt` were derived by Haxby's authors from
# these same subjects' responses, over runs the distribution does not name.
# The transport therefore declares `provenance$method = "external"` -- the
# contract's name for an operator whose construction crossform did not observe
# -- and can declare no `cross_fit`. Reading a face - house contrast inside a
# face-selective territory is consequently NOT protected against circularity,
# and no number here should be read as evidence that face-selective cortex
# discriminates faces from houses. What the slice is evidence about is the
# population ALGEBRA.
#
# WHAT IS CLAIMED. The group face - house and animacy contrasts per group ROI,
# with between-subject standard errors; a group RDM over the eight conditions
# at each node; the two transported component ledgers under the names contract
# section 8.1 requires (`native_coherent_ledger`, `native_configuration_ledger`
# -- never the bare words, which would invite reading them as group-node
# quantities they are not); a descriptive cross-fitted heterogeneity spectrum
# with subject loadings. WHAT IS NOT: six subjects is six subjects, the
# intervals are UNCALIBRATED (`population_uncertainty()` says so in every row
# it emits), the heterogeneity spectrum is one draw and not an estimate of a
# between-subject trace, the territory definition is circular in the sense
# above, and the transport is a research decision that is part of the estimand
# (contract section 1.5), not a neutral preprocessing step.
#
# Idempotent and resumable: subjects whose prepared object is absent, or who
# are missing one of the three territories, are skipped and named -- never
# carried with a fabricated zero at the node they cannot reach. Runtime is a
# few seconds -- region-level nodes are three rows per subject, not six
# hundred.
#
# Environment: SUBJECTS, EXEMPLAR_DIR.
#
# Output: results/population-slice1-receipts.csv  (committed evidence),
#         results/population-slice1.rds,
#         results/population-slice1.png.

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

## ---- Settings -------------------------------------------------------------
# The population layer's own tolerance (contract section 11): budget
# preservation and commutation are asserted at 1e-12, scaled where the
# contract scales them.
POPULATION_TOLERANCE <- 1e-12
# The conservative-geometry tie back to 09 is an accumulation over a few
# hundred voxels and keeps 09's own acceptance.
CONSERVATION_TOLERANCE <- 1e-10

REGION_LEVELS <- c("face-territory", "house-territory", "other-VT")
ROI_FILES <- c(`face-territory` = "mask8_face_vt.nii.gz",
               `house-territory` = "mask8_house_vt.nii.gz")

named_contrast <- function(positive, negative) {
  w <- stats::setNames(numeric(length(CONDITIONS)), CONDITIONS)
  w[positive] <- 1 / length(positive)
  w[negative] <- -1 / length(negative)
  w
}
# A two-query bank, read in one pass at every group node. Both rows are
# zero-sum, which is what the distance basis requires of a contrast.
QUERY_BANK <- rbind(
  `face-house` = named_contrast("face", "house"),
  `animate-inanimate` = named_contrast(ANIMATE, setdiff(CONDITIONS, ANIMATE))
)
stopifnot(max(abs(rowSums(QUERY_BANK))) < 1e-12)
HEADLINE_QUERY <- "face-house"

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
if (length(subjects) < 2L) {
  stop("A population plan needs at least two subjects; ", length(subjects),
       " prepared. Run 00-download.R and 08-prepare-subjects.R first.")
}
message("Subjects in scope: ", paste(subjects, collapse = ", "))

script_t0 <- Sys.time()

## ---- Receipts -------------------------------------------------------------
# 09's shape, with `subject` as the leading key and "group" for the rows that
# belong to the population fit rather than to one participant. `tolerance` is
# finite only on rows that assert an identity, and on those rows `passes` is
# exactly `abs(value) <= tolerance`.
receipts <- list()
record <- function(subject, quantity, value, tolerance = NA_real_, note = "") {
  # A receipt row is one number. A field that silently arrived empty would
  # otherwise drop the row and leave the CSV looking complete.
  if (length(value) != 1L || length(tolerance) != 1L) {
    stop("receipt \"", quantity, "\" for ", subject, ": value has length ",
         length(value), " and tolerance length ", length(tolerance),
         "; both must be scalars.")
  }
  passes <- if (is.na(tolerance)) NA else abs(value) <= tolerance
  receipts[[length(receipts) + 1L]] <<- data.frame(
    subject = subject, quantity = quantity, value = as.numeric(value),
    tolerance = as.numeric(tolerance), passes = passes, note = note,
    stringsAsFactors = FALSE
  )
  invisible(value)
}

## ---- Per subject: region frame, geometry plan, transport ------------------
message("\n== per-subject region frames ==")
built <- list()
failures <- list()

for (subj in subjects) {
  out <- tryCatch({
    sp <- subject_paths(exemplar_dir, subj)
    prep <- readRDS(sp$prepared)
    stopifnot(identical(prep$conditions, CONDITIONS))

    ## The domain is 09's domain, voxel for voxel: the VT mask with the
    ## preparation's degenerate voxels removed. Building it the same way is
    ## what lets 10's whole-VT budget be compared with 09's.
    vt_mask <- vt_mask_volume(sp$mask, drop = prep$dropped_index)
    domain <- neuroim2_volume_domain(vt_mask,
                                     id = paste0("haxby-", subj, "-vt"))
    stopifnot(identical(as.integer(domain$feature_ids),
                        as.integer(prep$vt_index)))
    features <- as.integer(domain$feature_ids)

    ## The territories. A subject missing an ROI file, or whose territory is
    ## empty after intersecting with the analyzed VT, is EXCLUDED by name.
    ## It is tempting to let such a subject carry two nodes instead of three
    ## and call the third "not fed", but that is not what happens: under
    ## budget semantics `transport_values()` returns a hard 0 at the unreached
    ## group node, and a 0 is an observation -- it enters the group OLS stack
    ## and pulls the mean toward zero. Under density it returns NA, which
    ## propagates silently into the fit. Neither is "simply not fed", so the
    ## subject is dropped and said out loud instead.
    roi_present <- file.exists(file.path(sp$subj, ROI_FILES))
    names(roi_present) <- names(ROI_FILES)
    if (!all(roi_present)) {
      stop("missing ROI mask(s): ",
           paste(ROI_FILES[!roi_present], collapse = ", "),
           " -- the three-territory partition is not available for this subject")
    }
    roi_index <- function(name) {
      intersect(which(as.array(neuroim2::read_vol(
        file.path(sp$subj, ROI_FILES[[name]]))) > 0), features)
    }
    face <- roi_index("face-territory")
    house <- roi_index("house-territory")
    ## `regions()` takes one label per feature, so an overlapping voxel needs
    ## a rule rather than a split. Face wins, which is 09's convention, and
    ## the overlap is counted so the rule is visible rather than implied.
    overlap <- length(intersect(face, house))
    labels <- ifelse(features %in% face, "face-territory",
              ifelse(features %in% house, "house-territory", "other-VT"))
    labels <- factor(labels, levels = REGION_LEVELS)
    region_counts <- stats::setNames(as.integer(tabulate(labels,
      length(REGION_LEVELS))), REGION_LEVELS)
    if (any(region_counts == 0L)) {
      stop("empty territory after intersecting with the analyzed VT: ",
           paste(REGION_LEVELS[region_counts == 0L], collapse = ", "),
           " -- see the note on excluded subjects in the script header")
    }

    effects <- effect_space(CONDITIONS, basis_id = "haxby-condition-means:v1",
                            units = "within-run-z")
    rel <- relation(prep$runs, effects = effects, domain = domain)
    over <- cross_partitions(rel, independence = "independent",
                             generalizes_over = "run")

    ## The conservative region frame: three rows, and every VT voxel's
    ## evidence assigned to exactly one of them.
    region_frame <- compile_frame(regions(labels, normalization = "conservative"),
                                  domain)
    region_conservation <- frame_conservation(region_frame)
    stopifnot(region_conservation$conserved)
    region_plan <- plan_geometry(rel, at = region_frame, over = over)

    ## The global comparator, exactly as 09 reads it: whole VT, UNNORMALIZED.
    whole_plan <- plan_geometry(rel,
      at = compile_frame(whole_brain("none"), domain), over = over)

    ## Native ledgers, per query, read BEFORE any transport. These are the
    ## left-hand side of the commutation check.
    query_names <- rownames(QUERY_BANK)
    native_order <- as.character(region_frame$index$measurement)
    native_ledger <- vapply(query_names, function(q)
      as.numeric(contrast_energy(region_plan, QUERY_BANK[q, ])$total),
      numeric(length(native_order)))
    rownames(native_ledger) <- native_order
    whole_total <- vapply(query_names, function(q)
      contrast_energy(whole_plan, QUERY_BANK[q, ])$total[[1L]], numeric(1))

    list(subject = subj, prep = prep, domain = domain,
         plan = region_plan, frame = region_frame,
         conservation = region_conservation,
         labels = labels, region_counts = region_counts,
         roi_present = roi_present,
         overlap = overlap, native_order = native_order,
         native_ledger = native_ledger, whole_total = whole_total,
         n_runs = prep$n_runs, n_run_pairs = nrow(over),
         n_voxels = domain$n_features,
         partitions = rel$partitions)
  }, error = function(e) e)

  if (inherits(out, "error")) {
    message("  ", subj, ": FAILED -- ", conditionMessage(out))
    failures[[subj]] <- conditionMessage(out)
    next
  }
  built[[subj]] <- out
  message(sprintf("  %s: VT %d vox -> %s | %d runs, %d pairs",
    subj, out$n_voxels,
    paste(sprintf("%s %d", names(out$region_counts), out$region_counts),
          collapse = ", "),
    out$n_runs, out$n_run_pairs))
}

if (length(built) < 2L) {
  stop("Fewer than two subjects completed; a population plan needs two.")
}
subjects <- names(built)

## ---- The transport --------------------------------------------------------
# One group index, built once and shared by every subject: `plan_population()`
# compares the group indices with `identical()`, which is the right strictness
# -- two subjects carried onto "the same" node set that disagree on a column
# are not on the same node set.
group_index <- data.frame(node = REGION_LEVELS,
                          territory = c("Haxby mask8_face_vt",
                                        "Haxby mask8_house_vt",
                                        "mask4_vt minus both"),
                          stringsAsFactors = FALSE)

region_transport <- function(entry, semantics = "budget") {
  native <- entry$native_order
  operator <- matrix(0, length(native), length(REGION_LEVELS),
                     dimnames = list(native, REGION_LEVELS))
  # `regions()` does not order its rows by `REGION_LEVELS`, and the order it
  # does use differs between subjects. Matching by NAME rather than by
  # position is what makes the map the identity it claims to be; a positional
  # assignment here would silently relabel two subjects' territories.
  operator[cbind(seq_along(native), match(native, REGION_LEVELS))] <- 1
  stopifnot(!anyNA(operator), all(rowSums(operator) == 1))
  location_transport(
    operator,
    native_index = data.frame(node = native,
      voxels = as.integer(entry$region_counts[native]),
      stringsAsFactors = FALSE),
    group_index = group_index,
    semantics = semantics,
    # The declared row mass, and the reason density is worth reading here: a
    # territory's own size. Under "budget" this vector is carried but unused;
    # under "density" it is the denominator (contract section 1.3).
    row_mass = as.numeric(entry$region_counts[native]),
    # `method = "external"`, NOT "anatomical", and the distinction is the
    # honest one. `mask8_face_vt` and `mask8_house_vt` are FUNCTIONAL ROIs:
    # Haxby's authors defined them from these same subjects' responses. A
    # transport built from data is `"functional"` under contract section 1.4
    # and then owes a `cross_fit` naming the runs it was fitted on -- and the
    # distribution does not say which runs those were, so that field cannot
    # be filled honestly. `"external"` is the contract's name for exactly this
    # case: an operator whose construction crossform did not observe. Section
    # 7.2 measures what an undeclared circular transport buys (3.15x the
    # honest gain), which is why this is stated here and again in the
    # caveats rather than quietly rounded to "anatomical".
    provenance = list(
      method = "external",
      details = paste(
        "identity region map: each subject's native territory maps wholly to",
        "the same-named group ROI. Territories are the subject's own",
        "Haxby-distributed mask8_face_vt / mask8_house_vt intersected with",
        "mask4_vt, and their complement inside mask4_vt. crossform performs",
        "no registration and no resampling here, and the map itself uses no",
        "response data -- but the ROIs that DEFINE the territories were",
        "derived by the original authors from these subjects' own responses,",
        "over runs this distribution does not name, so no cross-fit can be",
        "declared and the face/house territory readings are not protected",
        "against circularity (population-form-v1 sections 1.4 and 7.2)."),
      atlas = "Haxby 2001 subject-native FUNCTIONAL ROIs (mask8_*_vt)",
      circularity = paste(
        "unquantified: the territory definition and the contrast evaluated",
        "in it come from the same responses")
    )
  )
}
transports <- lapply(built, region_transport)
# The SAME operator read under the other semantics. Contract section 1.3:
# budget and density are two linear maps built from one P, and they answer
# different questions -- here, "how much of this subject's VT evidence sits in
# this territory" versus "how much per voxel of it". At region level the
# territories differ in size by more than an order of magnitude, so the two
# readings are not a presentation choice and the budget ranking on its own
# would be read as a claim about location when it is partly a claim about
# size.
density_transports <- lapply(built, region_transport, semantics = "density")

## ---- The population plan --------------------------------------------------
message("\n== population plan ==")
population_plan <- plan_population(
  subjects = lapply(built, `[[`, "plan"),
  transport = transports,
  model = ~ 1,
  normalization = "none"
)
print(population_plan)

## `plan_population()` sorts its participants by name so that plan identity is
## order-invariant, and every array it returns carries that order. Adopting it
## here keeps the figure colours, the legend and the receipts indexed the same
## way as `fit$values`; using the insertion order instead would mislabel the
## per-subject points whenever SUBJECTS was given out of order.
stopifnot(setequal(population_plan$subject_index$subject, names(built)))
subjects <- population_plan$subject_index$subject

## ---- The group fit --------------------------------------------------------
message("\n== group fit ==")
fit <- estimate_population(population_plan, QUERY_BANK)
coherent <- estimate_population(population_plan, QUERY_BANK,
                                component = "coherent")
configuration <- estimate_population(population_plan, QUERY_BANK,
                                     component = "configuration")
print(fit)

node_ids <- dimnames(fit$values)$node          # group nodes, then "<sink>"
query_names <- rownames(QUERY_BANK)
group_nodes <- node_ids[!fit$index$sink]
sink_id <- node_ids[fit$index$sink]

## Uncertainty: between-subject only, and labelled uncalibrated by the verb
## itself. There is no within-subject layer here because no per-subject
## sampling covariance was supplied, and the two are never pooled anyway.
uncertainty <- population_uncertainty(fit)

## The same six subjects, the same operator, read under density semantics.
## `plan_population()` requires one semantics per plan, so this is a second
## plan with a different scientific identity -- which is the point: two
## transports are two estimands (contract section 1.5).
density_plan <- plan_population(
  subjects = lapply(built, `[[`, "plan"),
  transport = density_transports,
  model = ~ 1,
  normalization = "none"
)
density_fit <- estimate_population(density_plan, QUERY_BANK)
density_uncertainty <- population_uncertainty(density_fit)

## The whole-form readout, and the group RDM off it. `materialize_population()`
## reads the complete q x q form, so every symmetric query -- including all 28
## condition pairs -- is in span rather than being approximated.
form <- materialize_population(population_plan)
group_rdm <- rdm(form)
rdm_table <- as.data.frame(group_rdm)

## The aggregated ledger over all three group nodes: the package's own
## aggregation, used here as the E6 sum-over-nodes identity.
whole_ledger <- contribution(fit, by = rep("VT", length(group_nodes)))

## Heterogeneity: cross-fitted, which is the only estimator the contract
## admits for a spectrum (section 6.4). Each subject's own runs are
## interleaved into two halves; 11-12 runs per subject is far above the
## four-partition floor.
message("\n== heterogeneity (cross-fitted) ==")
het <- heterogeneity(population_plan, modes = min(2L, length(built) - 1L))
print(het)

## Prevalence. Two readings, both on the latent descriptive layer: the sign
## count (how many participants carry a positive ledger at this node and
## query) and the alignment count (how many agree in direction with the
## leave-one-out mean of the others). The second is the more informative one
## precisely because the reference excludes the participant being scored, so
## every inner product is cross-participant. A pure-noise cell reports 0.5 in
## both, which is the reference the verb carries rather than a null model.
##
## The verb also supplies the two section 7.5 coverage diagnostics. An earlier
## draft of this script carried a fallback that recomputed the sign fraction by
## hand when the verb was absent -- but that fallback also had to invent the
## coverage numbers, and recording "every node is fully covered" without
## computing it is worse than not recording it at all. The verb is exported,
## so it is simply required.
if (!is.function(get0("population_prevalence"))) {
  stop("population_prevalence() is not available in this crossform build; ",
       "slice 1 reports prevalence and the section 7.5 coverage diagnostics ",
       "through it and does not approximate them.")
}
prevalence_record <- population_prevalence(fit,
                                           coverage_floor = length(subjects))
prevalence <- prevalence_record$sign$fraction

## ---- Identity checks ------------------------------------------------------
# Every one is computed from the objects above rather than re-derived, so a
# change in the package moves these numbers.
#
# A WORD ON WHAT AN EXACT ZERO WOULD MEAN. Both population executors query
# first and transport afterwards: `estimate_population()` runs the geometry
# compiler on the packed contrast and hands the resulting per-node numbers to
# `transport_values()`, and `materialize_population()` does the same with the
# packed-coordinate probe. So a "commutation" check written as
# `transport_values(P, contrast_energy(plan, c))` against
# `fit$values` is not a check at all -- `contrast_energy()` lowers to exactly
# the packed operator the bank carries, so both sides are the same call on the
# same numbers and the answer is 0 by construction. That comparison is worth
# recording, but as what it is: an executor-plumbing check (2). The genuine
# order reversal is (5b).
message("\n== identity checks ==")

## (1) Budget preservation through the transport, per subject.
budget_deviation <- stats::setNames(
  fit$receipt$subjects$budget_deviation, fit$receipt$subjects$subject)

## (2) Executor plumbing: the group fit's per-subject response really is
##     `transport_values()` applied to the subject's own `contrast_energy()`
##     ledger, and not something else that happens to be close. Exact by
##     construction, and recorded so that a future executor that started
##     rescaling or reordering rows would show up here. NOT a commutation
##     check -- see the note above.
executor_agreement <- vapply(subjects, function(subj) {
  entry <- built[[subj]]
  carrier <- transports[[subj]]
  carried <- transport_values(carrier,
    entry$native_ledger[as.character(carrier$native_index$node), ,
                        drop = FALSE])
  max(abs(carried[node_ids, query_names, drop = FALSE] -
          fit$values[node_ids, query_names, subj]))
}, numeric(1))

## (2c) Density does NOT conserve. This is arithmetic, not a discovery: each
##      group column is divided by its own territory size, so the column sum
##      misses the native total by roughly `1 - 1/(mean territory size)`. It is
##      recorded as an inequality rather than an equality precisely because
##      there is no law here to assert -- contract 1.3 says a density column
##      sum satisfies no conservation law, and this is what that looks like.
##      (The sum is also unit-inhomogeneous: the sink row stays in budget units
##      under either semantics. Harmless here only because the sink is zero.)
density_nonconservation <- vapply(subjects, function(subj)
  max(abs(colSums(density_fit$values[, , subj]) -
          fit$receipt$native_total[subj, query_names]) /
      pmax(abs(fit$receipt$native_total[subj, query_names]),
           .Machine$double.eps)),
  numeric(1))

## (3) Sum over group nodes + sink equals the transported total, per subject.
node_sum_deviation <- vapply(subjects, function(subj)
  max(abs(colSums(fit$values[, query_names, subj]) -
          fit$receipt$native_total[subj, query_names])),
  numeric(1))

## (4) The E6 sum-over-nodes identity at the group: the aggregated ledger over
##     group nodes AND the sink reproduces the group coefficient of the summed
##     response, and dropping the sink breaks it.
group_response <- apply(fit$values, c(3L, 2L), sum)
group_global <- qr.coef(qr(population_plan$model$matrix),
                        group_response)["(Intercept)", ]
theta_identity <- max(abs(colSums(whole_ledger$values) - group_global))
theta_without_sink <- max(abs(
  colSums(whole_ledger$values[!whole_ledger$index$sink, , drop = FALSE]) -
    group_global))

## (5) The transported ledger identity: coherent + configuration = total.
##     This one CANNOT fail, and saying so is part of recording it: the
##     configuration component is computed as `total - coherent` inside each
##     subject's own execution, before any transport, so this is a rounding
##     identity and an error in `coherent` would be absorbed by
##     `configuration` and still pass. Contract 8.1 requires the sentence to
##     be stated, so it is stated -- as documentation of the decomposition's
##     shape, not as evidence about its content.
ledger_identity <- max(abs(
  coherent$values + configuration$values - fit$values))

## (5b) THE COMMUTATION ACCEPTANCE (contract claim 3, query <-> transport).
##      This is the one place the two operations genuinely change order.
##      `estimate_population()` contracts the face - house contrast into a
##      single packed operator FIRST and transports one number per node:
##      query, then transport. `materialize_population()` transports all 36
##      packed coordinates of the form and the `rdm()` view contracts the
##      (face, house) edge out of them AFTERWARDS: transport, then query. The
##      contrast IS that RDM edge, so the two orders must land on the same
##      number -- and unlike (2) they are not the same arithmetic, which is
##      why this one lands at a rounding scale rather than at exactly zero.
face_house_edge <- rdm_table[rdm_table$left == "face" &
                             rdm_table$right == "house", ]
# One row per node: with a richer group model `rdm()` emits one row per node
# per term, and `match()` would silently take the first.
stopifnot(nrow(face_house_edge) == length(node_ids),
          !anyDuplicated(face_house_edge$node))
commutation <- max(abs(
  face_house_edge$estimate[match(node_ids, face_house_edge$node)] -
    fit$coefficients[node_ids, "face-house", "(Intercept)"]))

## (5c) The same order reversal under DENSITY semantics. Contract 1.3's first
##      consequence is that density is still a fixed linear map -- `P` and the
##      row mass are declared, not estimated -- so everything section 3 proves
##      about commutation holds for it too. What density gives up is
##      conservation (2c), not commutation.
density_form <- materialize_population(density_plan)
density_rdm_table <- as.data.frame(rdm(density_form))
density_edge <- density_rdm_table[density_rdm_table$left == "face" &
                                  density_rdm_table$right == "house", ]
stopifnot(nrow(density_edge) == length(node_ids))
density_commutation <- max(abs(
  density_edge$estimate[match(node_ids, density_edge$node)] -
    density_fit$coefficients[node_ids, "face-house", "(Intercept)"]))

## (6) The tie back to 09, asserted against 09's OWN committed number rather
##     than against a number this script also computed. Two separate scripts,
##     two separate frames (three regions here, one whole-VT comparator
##     there), and the budget has to be the same or they are not reading the
##     same subjects.
conservation_deviation <- vapply(subjects, function(subj)
  max(abs(colSums(built[[subj]]$native_ledger) - built[[subj]]$whole_total)),
  numeric(1))
nine_receipts_path <- file.path(paths$results,
                                "six-subject-conservative-receipts.csv")
nine_agreement <- if (file.exists(nine_receipts_path)) {
  nine <- utils::read.csv(nine_receipts_path, stringsAsFactors = FALSE)
  nine <- nine[nine$quantity == "whole_vt_total", ]
  shared <- intersect(subjects, nine$subject)
  if (length(shared)) {
    max(abs(vapply(shared, function(subj)
      built[[subj]]$whole_total[[HEADLINE_QUERY]], numeric(1)) -
      nine$value[match(shared, nine$subject)]))
  } else {
    NA_real_
  }
} else {
  NA_real_
}

message(sprintf("  budget preservation      %.2e  (<= %g)",
                max(abs(budget_deviation)), POPULATION_TOLERANCE))
message(sprintf("  executor plumbing        %.2e  (exact by construction)",
                max(executor_agreement)))
message(sprintf("  COMMUTATION (claim 3)    %.2e  (<= %g)",
                commutation, POPULATION_TOLERANCE))
message(sprintf("  commutation (density)    %.2e  (<= %g)",
                density_commutation, POPULATION_TOLERANCE))
message(sprintf("  node sum + sink          %.2e  (<= %g)",
                max(node_sum_deviation), POPULATION_TOLERANCE))
message(sprintf("  sum_u Theta identity     %.2e  (<= %g)",
                theta_identity, POPULATION_TOLERANCE))
message(sprintf("    (with the sink row deleted: %.2e -- still an identity, because",
                theta_without_sink))
message("     the region map has full coverage and the sink is exactly empty)")
message(sprintf("  coherent+configuration   %.2e  (rounding; cannot fail)",
                ledger_identity))
message(sprintf("  region ledger vs whole   %.2e  (<= %g)",
                max(conservation_deviation), CONSERVATION_TOLERANCE))
message(sprintf("  budget agrees with 09    %s",
                if (is.na(nine_agreement)) "09 receipts absent, not checked"
                else sprintf("%.2e  (<= %g)", nine_agreement,
                             CONSERVATION_TOLERANCE)))

## ---- Record ---------------------------------------------------------------
for (subj in subjects) {
  entry <- built[[subj]]
  record(subj, "vt_voxels_analyzed", entry$n_voxels,
         note = "domain features: mask4_vt voxels that vary in every run (09's domain)")
  for (region in REGION_LEVELS) {
    record(subj, paste0("native_region_voxels: ", region),
           entry$region_counts[[region]],
           note = sprintf("voxels in this subject's %s; a native node exists only when this is > 0",
                          region))
  }
  record(subj, "native_nodes", length(entry$native_order),
         note = "region rows in the conservative frame = rows of this subject's transport")
  record(subj, "roi_overlap_voxels", entry$overlap,
         note = "voxels in both mask8_face_vt and mask8_house_vt; assigned to face-territory (09's convention)")
  record(subj, "n_runs", entry$n_runs,
         note = "runs contributing a full eight-condition mean matrix")
  record(subj, "n_run_pairs", entry$n_run_pairs,
         note = "unordered run pairs the crossvalidated estimate averages over, choose(runs, 2)")

  record(subj, "frame_column_mass_max_deviation",
         entry$conservation$max_deviation, tolerance = CONSERVATION_TOLERANCE,
         note = "max |sum_x w_xv - 1| over voxels: the three territories partition VT exactly once")
  for (q in query_names) {
    record(subj, paste0("whole_vt_total: ", q), entry$whole_total[[q]],
           note = "the budget: contrast total over whole VT under whole_brain(\"none\")")
    record(subj, paste0("native_region_ledger_sum: ", q),
           sum(entry$native_ledger[, q]),
           note = "sum over the subject's three territory nodes, before transport")
  }
  record(subj, "region_ledger_vs_whole_vt", conservation_deviation[[subj]],
         tolerance = CONSERVATION_TOLERANCE,
         note = "max over queries of sum_x region_x - whole_VT_total: the tie back to 09")

  record(subj, "transport_budget_preservation", budget_deviation[[subj]],
         tolerance = POPULATION_TOLERANCE,
         note = "|transported total incl. sink - native total| / L1 norm of the native ledger (contract claim 2)")
  record(subj, "transport_executor_agreement", executor_agreement[[subj]],
         tolerance = POPULATION_TOLERANCE,
         note = paste0("|transport_values() of this subject's own contrast_energy() ledger - the group fit's ",
                       "response for it|. Exact by construction (both sides query before transporting), so it ",
                       "pins the executor's plumbing and is NOT the commutation claim; that is commutation_claim3"))
  record(subj, "density_nonconservation_relative",
         density_nonconservation[[subj]],
         note = paste0("max relative gap between the density column sum and the native total. Arithmetic, not a ",
                       "law: dividing each column by its territory size misses the total by about ",
                       "1 - 1/(mean territory size). Recorded because contract 1.3 says a density conserves nothing"))
  record(subj, "transport_node_sum_plus_sink", node_sum_deviation[[subj]],
         tolerance = POPULATION_TOLERANCE,
         note = "max over queries of sum over group nodes and sink minus the transported total")
  record(subj, "transport_sink_budget",
         max(abs(fit$receipt$sink_budget[subj, ])),
         tolerance = POPULATION_TOLERANCE,
         note = "sink budget in budget units; exactly zero because the region map has full coverage")
  record(subj, "transport_sink_territory",
         population_plan$subject_index$sink_territory[
           match(subj, population_plan$subject_index$subject)],
         tolerance = POPULATION_TOLERANCE,
         note = "unmapped native territory, a property of P alone (contract 7.5); zero at region level")

  for (q in query_names) {
    for (node in group_nodes) {
      record(subj, paste0("subject_transported_value: ", node, ": ", q),
             fit$values[node, q, subj],
             note = "this subject's carried territory ledger at that group node; the group fit's response")
    }
  }
}

## Group rows.
record("group", "n_subjects", length(subjects),
       note = "participants in the population plan")
record("group", "n_group_nodes", length(group_nodes),
       note = "group ROI nodes, sink excluded")
record("group", "residual_df", uncertainty$between$residual_df,
       note = "N - rank(X) for the group model ~ 1")
record("group", "budget_preservation_worst",
       fit$receipt$budget$max_relative_deviation,
       tolerance = fit$receipt$budget$tolerance,
       note = paste0("the fit's own certificate over every participant, scale \"",
                     fit$receipt$budget$scale, "\""))
record("group", "commutation_claim3", commutation,
       tolerance = POPULATION_TOLERANCE,
       note = paste0("THE acceptance. |query-then-transport - transport-then-query| at the group nodes: ",
                     "estimate_population() contracts face-house into one packed operator and transports one ",
                     "number per node; materialize_population() transports all 36 packed coordinates and rdm() ",
                     "contracts the (face,house) edge afterwards (contract claim 3)"))
record("group", "commutation_claim3_density", density_commutation,
       tolerance = POPULATION_TOLERANCE,
       note = "the same order reversal under density semantics: a declared row-mass ratio is still a fixed linear map, so commutation survives it (contract 1.3, consequence 1)")
record("group", "executor_agreement_worst", max(executor_agreement),
       tolerance = POPULATION_TOLERANCE,
       note = "worst executor-plumbing deviation over participants; exact by construction, not a commutation result")
record("group", "budget_argmax_matches_density_argmax",
       as.numeric(identical(
         group_nodes[which.max(fit$coefficients[group_nodes, HEADLINE_QUERY,
                                                "(Intercept)"])],
         group_nodes[which.max(density_fit$coefficients[group_nodes,
                                                        HEADLINE_QUERY,
                                                        "(Intercept)"])])),
       note = paste0("1 if budget and density rank the same group node first for ",
                     HEADLINE_QUERY, "; the two semantics are two estimands and need not agree (contract 1.3, 1.5)"))
record("group", "sum_over_nodes_identity", theta_identity,
       tolerance = POPULATION_TOLERANCE,
       note = "sum over group nodes and sink of contribution() minus the group coefficient of the summed response (E6)")
record("group", "sum_over_nodes_without_sink", theta_without_sink,
       tolerance = POPULATION_TOLERANCE,
       note = paste0("the same sum with the sink row deleted. It passes HERE only because ",
                     "the region map has full coverage and the sink is exactly empty; on a ",
                     "transport that loses territory this is where the loss shows up, which ",
                     "is why the sink is a row and not an option (contract 1.1)"))
record("group", "ledger_identity_coherent_plus_configuration", ledger_identity,
       tolerance = POPULATION_TOLERANCE,
       note = paste0("native_coherent_ledger + native_configuration_ledger - transported_total (contract 8.1). ",
                     "A ROUNDING identity that cannot fail: configuration is computed as total - coherent before ",
                     "transport, so an error in coherent is absorbed by configuration. Recorded because 8.1 ",
                     "requires the statement, not as evidence about the decomposition"))
record("group", "budget_agrees_with_script_09",
       if (is.na(nine_agreement)) 0 else nine_agreement,
       tolerance = if (is.na(nine_agreement)) NA_real_ else CONSERVATION_TOLERANCE,
       note = if (is.na(nine_agreement))
         "09's committed receipts were not on disk, so the cross-script tie was NOT checked; value is a placeholder"
       else paste0("max |this script's whole-VT ", HEADLINE_QUERY,
                   " budget - the value committed in six-subject-conservative-receipts.csv|, over the shared ",
                   "subjects: two scripts, two frames, one budget"))

for (q in query_names) {
  for (node in group_nodes) {
    record("group", paste0("group_contrast: ", node, ": ", q),
           fit$coefficients[node, q, "(Intercept)"],
           note = "OLS group mean of the carried territory ledger, native evidence units")
    record("group", paste0("group_contrast_se: ", node, ": ", q),
           uncertainty$between$se[node, q, "(Intercept)"],
           note = "between-subject SE, UNCALIBRATED (see benchmarks/POPULATION-NULL-COVERAGE.md)")
    record("group", paste0("group_contrast_t: ", node, ": ", q),
           uncertainty$between$t[node, q, "(Intercept)"],
           note = "estimate / SE against t on the residual df; uncalibrated, and not a p-value")
    record("group", paste0("group_contrast_lower: ", node, ": ", q),
           uncertainty$between$lower[node, q, "(Intercept)"],
           note = sprintf("lower bound of the %g%% nominal between-subject interval; UNCALIBRATED",
                          uncertainty$between$level * 100))
    record("group", paste0("group_contrast_upper: ", node, ": ", q),
           uncertainty$between$upper[node, q, "(Intercept)"],
           note = sprintf("upper bound of the %g%% nominal between-subject interval; UNCALIBRATED",
                          uncertainty$between$level * 100))
    # The two transported component ledgers, under the names contract 8.1
    # requires. Their SUM is the transported total exactly (asserted above);
    # what is deliberately NOT recorded is their ratio. A coherence fraction
    # is a nonnegative functional of signed estimates and lives only on the
    # latent PSD layer (contract 8.1, and conservative-geometry-v1 section 6),
    # so emitting one here as a receipt would be the exact error the naming
    # rule exists to prevent.
    record("group", paste0("native_coherent_ledger: ", node, ": ", q),
           coherent$coefficients[node, q, "(Intercept)"],
           note = "native_coherent_ledger: native-node coherent evidence carried to this group node. NOT a group-node common mode (contract 8)")
    record("group", paste0("native_configuration_ledger: ", node, ": ", q),
           configuration$coefficients[node, q, "(Intercept)"],
           note = "native_configuration_ledger; adds with the coherent ledger to the transported total exactly")
    record("group", paste0("group_density: ", node, ": ", q),
           density_fit$coefficients[node, q, "(Intercept)"],
           note = "the SAME operator read under density semantics: group mean of transported budget per native voxel (contract 1.3)")
    record("group", paste0("group_density_se: ", node, ": ", q),
           density_uncertainty$between$se[node, q, "(Intercept)"],
           note = "between-subject SE of the density reading, UNCALIBRATED")
    record("group", paste0("prevalence_positive: ", node, ": ", q),
           prevalence[node, q],
           note = paste0("fraction of participants whose carried ledger exceeds 0 here; ",
                         "latent descriptive layer (population_prevalence()$sign), ",
                         "and a pure-noise cell reads 0.5, not 0"))
  }
  record("group", paste0("group_ledger_over_vt: ", q),
         whole_ledger$values[which(!whole_ledger$index$sink)[[1L]], q],
         note = "contribution() over all three group nodes: the group ledger for whole VT")
}
record("group", "group_node_min_subject_coverage",
       min(prevalence_record$coverage$minimum),
       note = "fewest participants contributing a finite nonzero transported value to any group node (contract 7.5)")
record("group", "group_nodes_below_coverage_floor",
       length(prevalence_record$coverage$below_floor),
       note = paste0("group nodes whose contributing-participant count falls below the declared floor of ",
                     length(subjects), "; a node below it is not a group estimate and the record marks it"))
record("group", "prevalence_noise_reference", prevalence_record$reference,
       note = "the fraction a pure-noise cell reports; this is what makes a prevalence readable without a permutation test")
record("group", "prevalence_readout_gram_deviation",
       prevalence_record$alignment$readout_gram_deviation,
       note = paste0("how far the two-query bank's Gram sits from the identity: the alignment ",
                     "inner product is Euclidean in the READOUT, not Frobenius in the forms, ",
                     "and a non-orthonormal bank weights the queries unequally"))
for (node in group_nodes) {
  record("group", paste0("prevalence_alignment: ", node),
         prevalence_record$alignment$fraction[[node]],
         note = paste0("fraction of participants whose carried form agrees in direction with the ",
                       "LEAVE-ONE-OUT mean of the others (",
                       prevalence_record$alignment$inner_product,
                       "); the reference excludes the participant being scored, so every product is cross-participant"))
}

## Heterogeneity rows.
record("group", "heterogeneity_gram_size", length(het$spectrum),
       note = paste0("the cross-fitted subject Gram is ", length(subjects),
                     " x ", length(subjects), ", so it has that many ",
                     "eigenvalues; its RANK is at most N - rank(X) = ",
                     length(subjects) - 1L))
for (k in seq_along(het$spectrum)) {
  record("group", paste0("heterogeneity_eigenvalue: mode", k),
         het$spectrum[[k]],
         note = "signed eigenvalue of the CROSS-FITTED subject Gram; one draw, not an estimate of a between-subject trace")
}
record("group", "heterogeneity_negative_modes", sum(het$spectrum < 0),
       note = paste0("negative eigenvalues, reported as-is and never clipped. The cross-fitted Gram MAY be ",
                     "indefinite (contract 6.4 measures it so in 100% of its own simulations); on this data ",
                     "it is not -- read this count together with heterogeneity_moved_share, which says how ",
                     "much mass the projection actually had to move"))
record("group", "heterogeneity_n_eff", het$latent$n_eff,
       note = paste0("effective mode count after ", het$latent$method,
                     "; a nonnegative functional, so latent layer only"))
record("group", "heterogeneity_moved_mass", het$latent$moved_mass,
       note = "mass the PSD projection moved (contract 6.5); recorded because the projection is never free -- even when, as here, it is nearly free")
record("group", "heterogeneity_moved_share", het$latent$moved_share,
       note = "the same, as a share of the signed trace; at machine epsilon this says the projection cost nothing on this data, NOT that projections are free in general")
for (subj in subjects) {
  for (k in seq_len(min(2L, ncol(het$loadings)))) {
    record(subj, paste0("heterogeneity_loading: mode", k),
           het$loadings[subj, k],
           note = "subject loading on the cross-fitted heterogeneity mode; sign is oriented, not meaningful in itself")
  }
}

## RDM rows: the face/house edge at every group node, which is the minimum the
## slice claims, plus the whole 28-pair matrix in the .rds.
record("group", "rdm_edge_matches_query_bank", commutation,
       tolerance = POPULATION_TOLERANCE,
       note = "the same number as commutation_claim3, recorded under the readout's own name: the (face,house) RDM edge from the materialized form IS the face-house coefficient from the query bank")
for (i in seq_len(nrow(face_house_edge))) {
  row <- face_house_edge[i, ]
  record("group", paste0("group_rdm_edge: ", row$node, ": face-house"),
         row$estimate,
         note = "group RDM edge from materialize_population(); the complete form, so the pair is in span exactly")
}
record("group", "group_rdm_pairs", nrow(group_rdm$columns),
       note = "condition pairs read from the materialized group form")

## The estimand identities, so a rerun that changes the plan changes the CSV.
record("group", "n_quirks",
       sum(vapply(built, function(e) length(e$prep$quirks), integer(1))),
       note = paste(c(unlist(lapply(built, function(e)
         if (length(e$prep$quirks))
           paste0(e$subject, ": ", paste(e$prep$quirks, collapse = "; "))
         else NULL)),
         "region transport is the identity on named territories: no displacement, no spreading, no sink"),
         collapse = " | "))

## ---- Receipts CSV ---------------------------------------------------------
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
narrowed <- !setequal(requested, SUBJECT_IDS) || !setequal(subjects, SUBJECT_IDS)
receipts_path <- file.path(
  paths$results,
  if (narrowed) "population-slice1-receipts-partial.csv"
  else "population-slice1-receipts.csv")
utils::write.csv(receipts, receipts_path, row.names = FALSE)
message("\nWrote ", receipts_path, " (", nrow(receipts), " receipts over ",
        length(subjects), " subjects, ", sum(tolerance_rows),
        " of them identities, all within tolerance)")
if (narrowed) {
  message("  Not all six subjects were fitted, so the committed six-subject\n",
          "  receipts were NOT overwritten. Rerun without SUBJECTS to refresh them.")
}

## ---- Figure ---------------------------------------------------------------
suffix <- if (narrowed) "-partial" else ""
png_path <- file.path(paths$results,
                      sprintf("population-slice1%s.png", suffix))
grDevices::png(png_path, width = 1500, height = 640, res = 130)
op <- graphics::par(mfrow = c(1, 2), mar = c(4.6, 4.8, 3.6, 1.2),
                    mgp = c(2.9, 0.8, 0), cex.main = 1.0, font.main = 1)
palette6 <- c("#2166ac", "#4393c3", "#92c5de", "#d6604d", "#b2182b", "#762a83")
cols <- stats::setNames(palette6[seq_along(subjects)], subjects)

## Left: the group contrast at each group node, with the between-subject
## interval, and every subject's own carried value drawn beside it.
est <- fit$coefficients[group_nodes, HEADLINE_QUERY, "(Intercept)"]
lo <- uncertainty$between$lower[group_nodes, HEADLINE_QUERY, "(Intercept)"]
hi <- uncertainty$between$upper[group_nodes, HEADLINE_QUERY, "(Intercept)"]
per_subject <- fit$values[group_nodes, HEADLINE_QUERY, , drop = TRUE]
ylim <- range(0, per_subject, lo, hi, finite = TRUE)
at <- seq_along(group_nodes)
plot(NA, xlim = c(0.5, length(group_nodes) + 0.5), ylim = ylim, xaxt = "n",
     xlab = "", ylab = sprintf("%s contrast energy (native evidence units)",
                               HEADLINE_QUERY),
     main = sprintf("Group %s per ROI, %d subjects", HEADLINE_QUERY,
                    length(subjects)))
## The territory sizes belong on the axis, not in a footnote: under budget
## semantics a bigger territory holds more budget, so the ranking below is
## partly a ranking of size. The density reading that removes it is in the
## receipts (`group_density: ...`).
voxel_span <- vapply(group_nodes, function(node) {
  counts <- vapply(built, function(e) e$region_counts[[node]], integer(1))
  sprintf("%d-%d vox", min(counts), max(counts))
}, character(1))
graphics::axis(1, at = at, labels = FALSE)
graphics::mtext(group_nodes, side = 1, at = at, line = 0.7, cex = 0.8)
graphics::mtext(voxel_span, side = 1, at = at, line = 1.7, cex = 0.62,
                col = "#666666")
graphics::abline(h = 0, lty = 2, col = "#9a9a9a")
for (i in at) {
  graphics::points(rep(i - 0.22, length(subjects)), per_subject[i, ],
                   pch = 16, cex = 0.9, col = cols[subjects])
  graphics::arrows(i + 0.12, lo[[i]], i + 0.12, hi[[i]], code = 3, angle = 90,
                   length = 0.05, col = "#222222", lwd = 1.6)
  graphics::points(i + 0.12, est[[i]], pch = 18, cex = 1.5, col = "#222222")
}
graphics::legend("topleft", bty = "n", cex = 0.7, ncol = 2,
  legend = c(subjects, "group mean"), pch = c(rep(16, length(subjects)), 18),
  col = c(cols[subjects], "#222222"))
graphics::mtext(
  sprintf(paste0("bars are 95%% between-subject intervals on %d df and are ",
                 "UNCALIBRATED; N = %d; budget semantics, so size counts"),
          uncertainty$between$residual_df, length(subjects)),
  side = 3, line = 0.25, cex = 0.62, col = "#b2182b")

## Right: the subject loadings on the two leading cross-fitted heterogeneity
## modes. This is the descriptive half of the slice.
n_modes <- min(2L, ncol(het$loadings))
load_mat <- t(het$loadings[subjects, seq_len(n_modes), drop = FALSE])
bp <- graphics::barplot(load_mat, beside = TRUE, col = c("#4d4d4d", "#bdbdbd"),
  border = NA, ylim = range(-1, 1, load_mat) * 1.15,
  ylab = "loading on the cross-fitted heterogeneity mode",
  names.arg = subjects, cex.names = 0.85,
  main = sprintf("Heterogeneity modes (cross-fitted Gram, %d nodes)",
                 length(group_nodes)))
graphics::abline(h = 0, col = "#9a9a9a")
graphics::legend("topright", bty = "n", cex = 0.72, fill = c("#4d4d4d", "#bdbdbd"),
  border = NA, legend = sprintf("mode %d  (eigenvalue %s)", seq_len(n_modes),
    format(signif(het$spectrum[seq_len(n_modes)], 3), big.mark = "")))
graphics::mtext(
  sprintf(paste0("cross-fitted, reported signed and unclipped (%d negative ",
                 "eigenvalue%s here); one draw, not an estimate of a ",
                 "between-subject trace"),
          sum(het$spectrum < 0), if (sum(het$spectrum < 0) == 1L) "" else "s"),
  side = 3, line = 0.25, cex = 0.62, col = "#b2182b")
graphics::par(op)
invisible(grDevices::dev.off())
message("Wrote ", png_path)

## ---- Save -----------------------------------------------------------------
script_seconds <- as.numeric(difftime(Sys.time(), script_t0, units = "secs"))
out <- list(
  arm = "crossform population form, slice 1: six subjects, region-level nodes",
  contract = "design/population-form-contract.md (population-form-v1), sections 1-6",
  estimand = paste(
    "the OLS group mean, over participants, of each participant's",
    "crossvalidated contrast energy read at three native VT territories under",
    "conservative normalization and carried onto three same-named group ROI",
    "nodes by an identity transport of declared external provenance under",
    "budget semantics, with",
    "budget normalization \"none\" (mean subject ledger, native evidence units)"),
  caveat = paste(
    "Six subjects. The between-subject intervals are UNCALIBRATED -- the",
    "package's own null-coverage benchmark measures 0.918 coverage for a",
    "nominal 0.95 interval at this residual df when participant noise is",
    "linked to group covariates, and it degrades further as N grows. The heterogeneity spectrum is ONE DRAW of an indefinite",
    "cross-fitted Gram and is not an estimate of a between-subject trace",
    "(population-form-v1 section 6.4). The transport is part of the estimand",
    "(section 1.5), not neutral preprocessing. No permutation test, no",
    "p-value, no multiple-comparison correction anywhere in this script."),
  query_bank = QUERY_BANK, headline_query = HEADLINE_QUERY,
  region_levels = REGION_LEVELS, group_index = group_index,
  population_tolerance = POPULATION_TOLERANCE,
  conservation_tolerance = CONSERVATION_TOLERANCE,
  subjects = subjects, requested = requested,
  missing_prepared = missing_prep, failures = failures,
  receipts = receipts,
  scientific_plan_ids = list(
    population = population_plan$scientific_plan_id,
    query_bank = fit$scientific_plan_id,
    complete_form = form$scientific_plan_id,
    heterogeneity = het$scientific_plan_id,
    density = density_plan$scientific_plan_id
  ),
  transport = list(
    semantics = population_plan$semantics,
    provenance = transports[[1L]]$provenance,
    signatures = vapply(transports, function(p) p$signature, character(1)),
    subject_index = population_plan$subject_index
  ),
  group = list(
    coefficients = fit$coefficients,
    values = fit$values,
    index = fit$index,
    ledger = fit$ledger,
    native_coherent_ledger_name = coherent$ledger,
    native_configuration_ledger_name = configuration$ledger,
    native_coherent_ledger_values = coherent$values,
    native_configuration_ledger_values = configuration$values,
    between = uncertainty$between[c("estimate", "se", "t", "lower", "upper",
                                    "residual_sd", "residual_df", "level",
                                    "calibration")],
    prevalence = prevalence,
    prevalence_record = prevalence_record[
      c("layer", "reference", "sign", "alignment", "coverage", "threshold")],
    whole_vt_ledger = as.data.frame(whole_ledger),
    rdm = rdm_table,
    heterogeneity = list(gram = het$gram, spectrum = het$spectrum,
                         loadings = het$loadings, latent = het$latent,
                         cross_fit = het$receipt$cross_fit)
  ),
  density = list(
    semantics = density_plan$semantics,
    row_mass = lapply(density_transports, `[[`, "row_mass"),
    coefficients = density_fit$coefficients,
    values = density_fit$values,
    between = density_uncertainty$between[c("estimate", "se", "t", "lower",
                                            "upper", "residual_df", "level",
                                            "calibration")]
  ),
  identities = list(
    budget_preservation = budget_deviation,
    commutation_density = density_commutation,
    density_nonconservation = density_nonconservation,
    node_sum_plus_sink = node_sum_deviation,
    sum_over_nodes = theta_identity,
    sum_over_nodes_without_sink = theta_without_sink,
    ledger_identity = ledger_identity,
    commutation_claim3 = commutation,
    executor_agreement = executor_agreement,
    budget_agrees_with_script_09 = nine_agreement,
    region_ledger_vs_whole_vt = conservation_deviation
  ),
  per_subject = lapply(built, function(e) e[c("subject", "region_counts",
    "roi_present", "overlap", "native_order",
    "native_ledger", "whole_total", "n_runs", "n_run_pairs", "n_voxels",
    "partitions")]),
  caveats = c(
    "six subjects, region-level nodes; not a searchlight population analysis",
    "the transport is an identity map on named territories -- it is the",
    "  simplest possible P, and slice 2 is where a transport has to be earned",
    "uncalibrated between-subject intervals; no permutation test, no p-values",
    "descriptive heterogeneity: one draw of an indefinite cross-fitted Gram",
    "the face and house territories are FUNCTIONALLY defined from these same",
    "  subjects; the transport declares external provenance and no cross-fit,",
    "  so those two nodes are not protected against circularity",
    "the headline is BUDGET semantics, so a bigger territory holds more of",
    "  it; the density reading is reported beside it and is a different",
    "  estimand, not a correction",
    "prevalence is a same-data fraction over six subjects, not an estimate",
    "eta_transport and the section 7 diagnostics belong to slice 2 and are",
    "  absent here: an identity transport has nothing to diagnose",
    "block design, per-run condition means rather than GLM betas",
    "within-run z-scored time series, 2 TR haemodynamic shift",
    "identity metric, native composition"),
  session = list(crossform = as.character(utils::packageVersion("crossform")),
                 crossform_source = crossform_source,
                 neuroim2 = as.character(utils::packageVersion("neuroim2")),
                 when = Sys.time()),
  seconds = script_seconds
)
rds_path <- file.path(paths$results,
                      sprintf("population-slice1%s.rds", suffix))
saveRDS(out, rds_path)
message("Wrote ", rds_path, " (", round(file.size(rds_path) / 1024), " KB)")

## ---- Summary --------------------------------------------------------------
message("\n== group ", HEADLINE_QUERY, " per group node ==")
print(data.frame(
  node = group_nodes,
  voxels = voxel_span,
  estimate = round(est, 3),
  se = round(uncertainty$between$se[group_nodes, HEADLINE_QUERY, "(Intercept)"], 3),
  t = round(uncertainty$between$t[group_nodes, HEADLINE_QUERY, "(Intercept)"], 3),
  lower = round(lo, 3), upper = round(hi, 3),
  positive_subjects = sprintf("%d/%d",
    round(prevalence[group_nodes, HEADLINE_QUERY] * length(subjects)),
    length(subjects)),
  per_voxel = round(density_fit$coefficients[group_nodes, HEADLINE_QUERY,
                                             "(Intercept)"], 4),
  per_voxel_se = round(density_uncertainty$between$se[group_nodes,
    HEADLINE_QUERY, "(Intercept)"], 4),
  row.names = NULL), row.names = FALSE)
message("df ", uncertainty$between$residual_df, ", ",
        uncertainty$between$level * 100, "% nominal, ",
        uncertainty$between$calibration)

message("\n== identities ==")
print(data.frame(
  quantity = c("budget preservation", "COMMUTATION (claim 3)",
               "commutation (density)", "executor plumbing (exact)",
               "node sum + sink", "sum over nodes (E6)",
               "coherent + configuration (rounding)",
               "region ledger vs whole VT", "budget agrees with 09"),
  worst = signif(c(max(abs(budget_deviation)), commutation,
                   density_commutation, max(executor_agreement),
                   max(node_sum_deviation), theta_identity, ledger_identity,
                   max(conservation_deviation),
                   if (is.na(nine_agreement)) NA_real_ else nine_agreement), 3),
  tolerance = c(rep(POPULATION_TOLERANCE, 7L),
                CONSERVATION_TOLERANCE, CONSERVATION_TOLERANCE),
  row.names = NULL), row.names = FALSE)
message("density conserves nothing, by construction: worst relative gap ",
        format(max(density_nonconservation), digits = 3),
        " between its column sum and the native total")

message("\n== prevalence (latent descriptive; noise reference ",
        prevalence_record$reference, ") ==")
print(prevalence_record)

message("\nheterogeneity spectrum: ",
        paste(signif(het$spectrum, 4), collapse = ", "),
        "  (", sum(het$spectrum < 0), " negative, n_eff ",
        signif(het$latent$n_eff, 4), ")")
if (length(failures)) {
  message("FAILED subjects: ", paste(names(failures), collapse = ", "))
}
message("Total: ", round(script_seconds, 1), " s")
