#!/usr/bin/env Rscript
# 02-analyze.R -- the crossform analysis: two relations over one neural
# domain, a rectangular encoding-by-retrieval plan, and the five pair-space
# readouts evaluated query-first for every simulated subject.
#
# Writes:
#   results/estimates-by-subject.csv  one row per subject x plan x query x region
#   results/query-diagnostics.csv     what each compiled operator is
#   results/route-check.csv           query-first vs materialize-then-project
#   results/refusals.csv              the boundaries this analysis ran into
#   data/analysis.rds                 queries + planted values for 03 (ignored)

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
source(file.path(exemplar_dir, "00-common.R"))
mode <- load_crossform(exemplar_dir)
results <- file.path(exemplar_dir, "results")
data_dir <- file.path(exemplar_dir, "data")
dir.create(results, showWarnings = FALSE, recursive = TRUE)
message("crossform loaded from ", mode)

simulation <- readRDS(file.path(data_dir, "simulation.rds"))
items <- simulation$items
probes <- simulation$probes
design <- simulation$design
truth <- simulation$truth
subjects <- names(simulation$responses)

## ---- Neural domain and spatial frame ------------------------------------
domain <- abstract_domain(
  N_VOXELS,
  coordinates = cbind(x = seq_len(N_VOXELS), y = 0, z = 0),
  id = "er-rsa-domain", coordinate_units = "mm"
)
frame <- compile_frame(regions(REGION_LABELS, normalization = "local"), domain)
measurements <- as.character(frame$index[[1L]])

## ---- Two relations over the same domain ---------------------------------
# Each run's responses are fitted once per axis. The encoding relation asks
# for the 36 encoding-trial coefficients, the retrieval relation for the 30
# retrieval-probe coefficients, out of the same design. Run baseline, two
# drift polynomials, and a motion trace are design columns in both fits and
# belong to neither effect space: they are nuisance, not effects.
fit_subject <- function(subject) {
  sources <- simulation$responses[[subject]]
  designs <- lapply(design$runs, `[[`, "design")
  encoding <- lm_relation_fit(
    sources, designs, design$encoding_target,
    sampling_unit = "trial", domain = domain,
    provenance = list(axis = "encoding", subject = subject)
  )
  retrieval <- lm_relation_fit(
    sources, designs, design$retrieval_target,
    sampling_unit = "trial", domain = domain,
    provenance = list(axis = "retrieval", subject = subject)
  )
  list(encoding = encoding$relation, retrieval = retrieval$relation,
       encoding_fit = encoding, retrieval_fit = retrieval)
}

first <- fit_subject(subjects[[1L]])
left_space <- first$encoding$effect_space
right_space <- first$retrieval$effect_space
stopifnot(
  identical(left_space$coordinates, paste0("enc_", items$item)),
  identical(right_space$coordinates, probes$probe)
)
message("left axis: ", length(left_space$coordinates), " encoding effects; ",
        "right axis: ", length(right_space$coordinates),
        " retrieval effects; rectangular, unequal, directed")

## ---- Generalization: the estimand excludes same-run products ------------
# Encoding and retrieval coefficients from the *same* run come out of one
# GLM, so their errors are correlated, and this simulation additionally
# plants a same-run item state shared by an item's two trials. The scientific
# estimand therefore pairs runs only across cycles. Both relations use the
# same partition names, so the ordered edge set says exactly that.
cross_edges <- expand.grid(left = names(design$runs), right = names(design$runs),
  stringsAsFactors = FALSE)
cross_edges <- cross_edges[cross_edges$left != cross_edges$right, ]
over_cross <- pairing(
  cross_edges$left, cross_edges$right,
  directed = TRUE, independence = "independent", generalizes_over = "run"
)
# The comparator, kept only to show what the estimand refuses to include.
over_same <- pairing(
  names(design$runs), names(design$runs),
  directed = TRUE, self_pairs = "allow_biased",
  independence = "not_independent"
)
message("cross-run edges: ", nrow(over_cross),
        " (", attr(over_cross, "estimate"), ")")
message("same-run comparator edges: ", nrow(over_same),
        " (", attr(over_same, "estimate"), ")")

## ---- Couplings: matches, controls, eligibility --------------------------
matched <- probes[probes$old, ]
match_left <- paste0("enc_", matched$source_item)
match_right <- matched$probe

# (1) match_coupling(): the declared encoding-retrieval correspondence.
matches_all <- match_coupling(match_left, match_right, left_space, right_space)

# The same matches, restricted to a category-matched eligible set. Matched
# pairs are always same-category, so an unrestricted control set compares
# within-category matches against mostly across-category controls. Any
# category signal then reads as reinstatement.
eligible_category <- expand.grid(
  left = seq_len(nrow(items)), right = seq_len(nrow(probes))
)
eligible_category <- eligible_category[
  items$category[eligible_category$left] ==
    probes$category[eligible_category$right], ]
eligible_category <- data.frame(
  left = left_space$coordinates[eligible_category$left],
  right = right_space$coordinates[eligible_category$right],
  stringsAsFactors = FALSE
)
matches_category <- match_coupling(match_left, match_right, left_space,
  right_space, eligible = eligible_category)

# (2) control_coupling(): the comparison cells, and the marginal baseline.
controls_all <- control_coupling(matches_all)
controls_category <- control_coupling(matches_category)
baseline_all <- control_coupling(matches_all, include_matches = TRUE)

## ---- Queries ------------------------------------------------------------
normalized_query <- function(coupling, constructor) {
  pair_query(coupling$value / sum(coupling$value),
    coupling$left_space, coupling$right_space,
    metadata = list(constructor = constructor))
}

queries <- list()
# Levels: the mean coupling over match cells and over control cells. These
# are the "match coupling > control coupling" numbers themselves, read as
# plain pair_query() operators over normalized coupling masses.
queries$level_match <- normalized_query(matches_all, "match level")
queries$level_control_all <- normalized_query(controls_all, "control level")
queries$level_control_category <-
  normalized_query(controls_category, "category-matched control level")
queries$level_baseline <- normalized_query(baseline_all, "eligible baseline")

# (3) coupling_contrast(): the difference of the normalized couplings.
queries$contrast_all <- coupling_contrast(matches_all, controls_all)
queries$contrast_category <-
  coupling_contrast(matches_category, controls_category)

# (4) match_control(): the same comparison as a nuisance-adjusted weighted
# least-squares coefficient over the eligible pair set.
queries$match_control_all <- match_control(matches_all)

refusals <- list()
record_refusal <- function(label, expression) {
  caught <- tryCatch(expression, error = function(e) e)
  refusals[[length(refusals) + 1L]] <<- data.frame(
    situation = label,
    class = paste(setdiff(class(caught), c("condition", "error")),
      collapse = "/"),
    capability = if (!is.null(caught$capability)) caught$capability else NA,
    message = conditionMessage(caught),
    stringsAsFactors = FALSE
  )
  invisible(caught)
}

# Both item nuisance families over a category-blocked eligible set are not
# jointly identified: inside each category block the encoding and retrieval
# level shifts trade off exactly. crossform refuses rather than silently
# dropping a column, so the exemplar drops one family explicitly.
category_both <- record_refusal(
  "match_control() with both nuisance families on category-blocked cells",
  match_control(matches_category)
)
queries$match_control_category <- if (inherits(category_both, "error")) {
  match_control(matches_category, retrieval_nuisance = FALSE)
} else {
  category_both
}

# (5) pair_lm_query(): the pair-space regression with the trial covariate.
# One row per eligible pair; `match` is the item-specific reinstatement
# indicator, `match_depth` carries the encoding study-duration covariate on
# matched cells only, `same_category` adjusts for the category confound the
# unrestricted control set otherwise absorbs, and both item nuisance
# families take out per-item level shifts.
pair_design <- expand.grid(
  left_index = seq_len(nrow(items)), right_index = seq_len(nrow(probes))
)
pair_design$left <- left_space$coordinates[pair_design$left_index]
pair_design$right <- right_space$coordinates[pair_design$right_index]
pair_design$match <- as.numeric(
  !is.na(probes$source_item[pair_design$right_index]) &
    items$item[pair_design$left_index] ==
      probes$source_item[pair_design$right_index]
)
centered_duration <- items$study_duration - truth$duration_center
pair_design$match_depth <-
  pair_design$match * centered_duration[pair_design$left_index]
pair_design$same_category <- as.numeric(
  items$category[pair_design$left_index] ==
    probes$category[pair_design$right_index])
pair_design <- pair_design[, c("left", "right", "match", "match_depth",
  "same_category")]

queries$lm_match <- pair_lm_query(
  pair_design, "match", left_space, right_space,
  encoding_nuisance = TRUE, retrieval_nuisance = TRUE
)
queries$lm_match_depth <- pair_lm_query(
  pair_design, "match_depth", left_space, right_space,
  encoding_nuisance = TRUE, retrieval_nuisance = TRUE
)

query_labels <- c(
  level_match = "mean match coupling",
  level_control_all = "mean control coupling (all eligible)",
  level_control_category = "mean control coupling (category-matched)",
  level_baseline = "mean coupling over every eligible pair",
  contrast_all = "match - control, all eligible pairs",
  contrast_category = "match - control, category-matched pairs",
  match_control_all = "match coefficient, item nuisance adjusted",
  match_control_category =
    "match coefficient, category-matched, encoding nuisance",
  lm_match = "match coefficient adjusting for same_category",
  lm_match_depth = "study-duration slope on matched pairs"
)
query_constructors <- c(
  level_match = "pair_query(match_coupling)",
  level_control_all = "pair_query(control_coupling)",
  level_control_category = "pair_query(control_coupling, eligible)",
  level_baseline = "pair_query(control_coupling(include_matches))",
  contrast_all = "coupling_contrast",
  contrast_category = "coupling_contrast",
  match_control_all = "match_control",
  match_control_category = "match_control",
  lm_match = "pair_lm_query",
  lm_match_depth = "pair_lm_query"
)

## ---- Compiled-operator diagnostics --------------------------------------
diagnostic_rows <- lapply(names(queries), function(name) {
  query <- queries[[name]]
  H <- as.matrix(query$operator)
  balance <- query$metadata$balance
  if (is.null(balance)) balance <- query$metadata$diagnostics$balance
  data.frame(
    query = name, readout = query_labels[[name]],
    constructor = query_constructors[[name]],
    eligible_cells = sum(H != 0),
    rank = if (is.null(query$metadata$diagnostics)) NA_integer_ else
      query$metadata$diagnostics$rank,
    design_columns = if (is.null(query$metadata$diagnostics)) NA_integer_ else
      length(query$metadata$diagnostics$columns),
    zero_row_marginals = if (is.null(balance)) NA else
      balance$zero_row_marginals,
    zero_column_marginals = if (is.null(balance)) NA else
      balance$zero_column_marginals,
    additive_baseline_invariant = if (is.null(balance)) NA else
      balance$additive_baseline_invariant,
    claim = if (is.null(query$metadata$claim)) NA_character_ else
      query$metadata$claim,
    stringsAsFactors = FALSE
  )
})
diagnostics <- do.call(rbind, diagnostic_rows)
write.csv(diagnostics, file.path(results, "query-diagnostics.csv"),
  row.names = FALSE)
message("\ncompiled pair operators:")
print(diagnostics[, c("query", "constructor", "eligible_cells", "rank",
  "additive_baseline_invariant")], row.names = FALSE)

## ---- Planted value of every query, per region and pairing ---------------
planted_rows <- list()
for (measurement in measurements) {
  weights <- as.numeric(frame$weights[
    match(measurement, measurements), , drop = TRUE])
  for (edges in c("cross", "same")) {
    G <- er_planted_geometry(truth, weights, edges)
    for (name in names(queries)) {
      planted_rows[[length(planted_rows) + 1L]] <- data.frame(
        plan = edges, query = name, region = measurement,
        planted = er_planted_value(queries[[name]], G),
        stringsAsFactors = FALSE
      )
    }
  }
}
planted <- do.call(rbind, planted_rows)

## ---- Execute, query-first, for every subject ----------------------------
plans_for <- function(relations) {
  list(
    cross = plan_geometry(relations$encoding, frame, over_cross,
      right = relations$retrieval),
    same = plan_geometry(relations$encoding, frame, over_same,
      right = relations$retrieval)
  )
}

estimate_rows <- list()
for (subject in subjects) {
  relations <- if (identical(subject, subjects[[1L]])) first else
    fit_subject(subject)
  plans <- plans_for(relations)
  for (plan_name in names(plans)) {
    for (name in names(queries)) {
      view <- evaluate_geometry(plans[[plan_name]], query = queries[[name]])
      estimate_rows[[length(estimate_rows) + 1L]] <- data.frame(
        subject = subject, plan = plan_name, query = name,
        region = as.character(view$index),
        estimate = as.numeric(view$values),
        stringsAsFactors = FALSE
      )
    }
  }
  message("subject ", subject, ": ", length(queries) * length(plans),
          " query-first evaluations")
}
estimates <- do.call(rbind, estimate_rows)
estimates <- merge(estimates, planted, by = c("plan", "query", "region"),
  all.x = TRUE)
estimates$readout <- query_labels[estimates$query]
estimates <- estimates[order(estimates$plan, estimates$query,
  estimates$region, estimates$subject), ]
write.csv(estimates, file.path(results, "estimates-by-subject.csv"),
  row.names = FALSE)

## ---- Route check: query-first equals materialize-then-project ----------
form <- materialize_geometry(plans_for(first)$cross)
total <- geometry_component(form, "total")
coherent <- geometry_component(form, "coherent")
configuration <- geometry_component(form, "configuration")
recomposition <- max(abs(total - (coherent + configuration)))
route_rows <- lapply(names(queries), function(name) {
  direct <- evaluate_geometry(plans_for(first)$cross, query = queries[[name]])
  projected <- query_geometry(form, queries[[name]])
  data.frame(
    query = name,
    max_abs_difference = max(abs(as.numeric(direct$values) -
      as.numeric(projected$values))),
    stringsAsFactors = FALSE
  )
})
route <- do.call(rbind, route_rows)
route <- rbind(route, data.frame(
  query = "total - (coherent + configuration)",
  max_abs_difference = recomposition, stringsAsFactors = FALSE))
write.csv(route, file.path(results, "route-check.csv"), row.names = FALSE)
message("\nquery-first vs materialized route: max abs difference = ",
        format(max(route$max_abs_difference[
          route$query != "total - (coherent + configuration)"]), digits = 3),
        "; rectangular recomposition = ", format(recomposition, digits = 3))

## ---- Recorded boundaries ------------------------------------------------
metric_refusal <- catch_refusal(plan_geometry(
  first$encoding, frame, over_cross, right = first$retrieval,
  metric = neural_metric(diag(N_VOXELS), domain)
))
refusals[[length(refusals) + 1L]] <- data.frame(
  situation = "a fixed noise metric on the rectangular plan",
  class = paste(setdiff(class(metric_refusal), c("condition", "error")),
    collapse = "/"),
  capability = metric_refusal$capability,
  message = conditionMessage(metric_refusal), stringsAsFactors = FALSE
)
rdm_refusal <- catch_refusal(rdm(plans_for(first)$cross))
refusals[[length(refusals) + 1L]] <- data.frame(
  situation = "rdm() on a rectangular cross-axis plan",
  class = paste(setdiff(class(rdm_refusal), c("condition", "error")),
    collapse = "/"),
  capability = rdm_refusal$capability,
  message = conditionMessage(rdm_refusal), stringsAsFactors = FALSE
)
refusal_table <- do.call(rbind, refusals)
write.csv(refusal_table, file.path(results, "refusals.csv"), row.names = FALSE)
message("\nrecorded boundaries:")
print(refusal_table[, c("situation", "capability")], row.names = FALSE)

saveRDS(
  list(queries = queries, planted = planted, diagnostics = diagnostics,
       labels = query_labels, constructors = query_constructors,
       measurements = measurements),
  file.path(data_dir, "analysis.rds")
)
message("\nwrote ", file.path(results, "estimates-by-subject.csv"))
