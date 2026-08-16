#!/usr/bin/env Rscript

# First-moment relation gate: executes every new public path from BIDS-shaped
# facts through a study-bound semantic model, relation plan, independent
# fmrireg point route, and second-moment point and uncertainty views.

if (!requireNamespace("devtools", quietly = TRUE) ||
    !requireNamespace("ps", quietly = TRUE) ||
    !requireNamespace("fmridesign", quietly = TRUE) ||
    !requireNamespace("fmrireg", quietly = TRUE)) {
  stop("The first-moment gate requires devtools, ps, fmridesign, and fmrireg.")
}

arguments <- commandArgs(trailingOnly = TRUE)
repo <- if (length(arguments)) {
  normalizePath(arguments[[1L]], mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
output_dir <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(repo, "benchmark-results")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(devtools::load_all(repo, quiet = TRUE))
source(file.path(repo, "tests", "testthat",
  "helper-relation-plan-fixtures.R"), local = TRUE)
source(file.path(repo, "tests", "testthat",
  "helper-first-moment-vertical-slice.R"), local = TRUE)

measure <- function(path, runner, repetitions = 3L) {
  elapsed <- resident <- result_size <- numeric(repetitions)
  last <- NULL
  for (iteration in seq_len(repetitions)) {
    gc()
    before <- unname(ps::ps_memory_info(ps::ps_handle(Sys.getpid()))[["rss"]])
    timing <- system.time(last <- runner())[["elapsed"]]
    after <- unname(ps::ps_memory_info(ps::ps_handle(Sys.getpid()))[["rss"]])
    elapsed[[iteration]] <- unname(timing)
    resident[[iteration]] <- max(0, after - before)
    result_size[[iteration]] <- as.numeric(utils::object.size(last))
  }
  list(
    value = last,
    record = data.frame(
      path = path,
      repetitions = repetitions,
      median_elapsed_seconds = stats::median(elapsed),
      max_resident_delta_bytes = max(resident),
      result_bytes = max(result_size),
      stringsAsFactors = FALSE
    )
  )
}

cell_build <- measure("plan_relation_cell_fixed", function() {
  first_moment_vertical_fixture("cell", "qr", "fixed_gls")
})
cell <- cell_build$value
treatment_build <- measure("plan_relation_treatment_fixed", function() {
  first_moment_vertical_fixture("treatment", "svd", "fixed_gls")
})
treatment <- treatment_build$value
native_run <- measure(
  "estimate_native_fixed", function() estimate_relation(cell$plan)
)
native <- native_run$value
conformance_run <- measure(
  "compiler_conformance", function() compiler_conformance(cell$plan)
)

ols <- first_moment_vertical_fixture("cell", "qr", "ols")
fmrireg_run <- measure(
  "fmrireg_relation_ols", function() fmrireg_relation(ols$plan)
)
fmrireg_fit <- fmrireg_run$value

geometry_plan_run <- measure(
  "plan_geometry",
  function() plan_geometry(native$relation, cell$frame, cell$over)
)
geometry_plan <- geometry_plan_run$value
contrast_weights <- c(face = 0.5, body = 0.5, house = -0.5, tool = -0.5)
contrast_run <- measure("contrast_query_first", function() {
  contrast_energy(geometry_plan, contrast_weights)
})
rdm_run <- measure("rdm_query_first", function() rdm(geometry_plan))
labels <- cell$conditions$coordinates
category <- outer(c(0, 0, 1, 1), c(0, 0, 1, 1), function(x, y) x != y)
storage.mode(category) <- "double"
dimnames(category) <- list(labels, labels)
rsa_run <- measure("rsa_query_first", function() {
  rsa(geometry_plan, models = list(category = category))
})
covariance_run <- measure("rdm_sampling_covariance", function() {
  rdm_sampling_covariance(geometry_plan, native, target = "plugin", at = 3L)
})

temporary <- tempfile("first-moment-bids-")
dir.create(temporary)
event_files <- confound_files <- stats::setNames(
  character(length(cell$partitions)), cell$partitions
)
for (partition in cell$partitions) {
  event_selected <- cell$study$events$data$partition == partition
  event_table <- cell$study$events$data[event_selected,
    c("onset", "duration", "condition")]
  event_files[[partition]] <- file.path(temporary,
    paste0(partition, "_events.tsv"))
  utils::write.table(event_table, event_files[[partition]], sep = "\t",
    row.names = FALSE, quote = FALSE)
  confound_selected <- cell$study$confounds$data$partition == partition
  confound_table <- cell$study$confounds$data[confound_selected,
    c("motion", "retained")]
  confound_files[[partition]] <- file.path(temporary,
    paste0(partition, "_desc-confounds_timeseries.tsv"))
  utils::write.table(confound_table, confound_files[[partition]], sep = "\t",
    row.names = FALSE, quote = FALSE)
}
bids_events_run <- measure(
  "bids_events", function() bids_events(event_files, cell$partitions)
)
bids_confounds_run <- measure("bids_confounds", function() {
  bids_confounds(
    confound_files,
    cell$partitions,
    lapply(cell$study$observations$indexes, `[[`, "observation_id"),
    censor = "retained"
  )
})
bids_study_run <- measure("bids_study", function() {
  bids_study(
    cell$study$observations,
    event_files,
    confound_files,
    partitions = cell$partitions,
    censor = "retained",
    hierarchy = cell$study$hierarchy
  )
})
imported_study <- bids_study_run$value

event_data <- cell$study$events$data
block <- match(event_data$partition, cell$partitions)
external_model <- fmridesign::event_model(
  onset ~ fmridesign::hrf(condition),
  data = event_data,
  block = block,
  sampling_frame = fmridesign::sampling_frame(
    blocklens = cell$counts, TR = 2, start_time = 0
  ),
  durations = event_data$duration
)
fmridesign_run <- measure("fmridesign_design_model", function() {
  fmridesign_design_model(
    external_model,
    imported_study,
    basis_id = "canonical-hrf-amplitude",
    units = "arbitrary-BOLD"
  )
})

direct_blocks <- first_moment_direct_blocks(cell)
direct_forms <- first_moment_direct_forms(direct_blocks, cell$frame, cell$over)
direct_rdm <- first_moment_direct_rdm(direct_forms)
native_block_error <- max(vapply(cell$partitions, function(partition) {
  max(abs(relation_block(native, partition, 1:60) -
    direct_blocks[[partition]]))
}, numeric(1)))
coding_fit <- estimate_relation(treatment$plan)
coding_error <- max(vapply(cell$partitions, function(partition) {
  max(abs(relation_block(native, partition, 1:60) -
    relation_block(coding_fit, partition, 1:60)))
}, numeric(1)))
fmrireg_error <- max(vapply(cell$partitions, function(partition) {
  max(abs(relation_block(estimate_relation(ols$plan), partition, 1:60) -
    relation_block(fmrireg_fit, partition, 1:60)))
}, numeric(1)))
rdm_error <- max(abs(as.matrix(rdm_run$value$values) - direct_rdm))

timings <- do.call(rbind, lapply(list(
  cell_build, treatment_build, native_run, conformance_run, fmrireg_run,
  geometry_plan_run, contrast_run, rdm_run, rsa_run, covariance_run,
  bids_events_run, bids_confounds_run, bids_study_run, fmridesign_run
), `[[`, "record"))
checks <- list(
  relation_plan_identity_stable = identical(
    cell$plan$relation_plan_id, treatment$plan$relation_plan_id
  ),
  design_receipts_route_specific = !identical(
    cell$plan$design_receipts$`run-1`$design_receipt_id,
    treatment$plan$design_receipts$`run-1`$design_receipt_id
  ),
  conformance_pass = all(as.matrix(conformance_run$value[-1L])),
  native_oracle_error = native_block_error,
  coding_error = coding_error,
  fmrireg_error = fmrireg_error,
  rdm_oracle_error = rdm_error,
  covariance_finite = all(is.finite(sampling_covariance(covariance_run$value))),
  timing_complete = nrow(timings) == 14L &&
    all(is.finite(timings$median_elapsed_seconds)),
  memory_complete = all(is.finite(timings$max_resident_delta_bytes)) &&
    all(is.finite(timings$result_bytes))
)
gate <- isTRUE(checks$relation_plan_identity_stable) &&
  isTRUE(checks$design_receipts_route_specific) &&
  isTRUE(checks$conformance_pass) &&
  checks$native_oracle_error < 1e-10 &&
  checks$coding_error < 1e-10 &&
  checks$fmrireg_error < 1e-10 &&
  checks$rdm_oracle_error < 1e-10 &&
  isTRUE(checks$covariance_finite) &&
  isTRUE(checks$timing_complete) &&
  isTRUE(checks$memory_complete)

result <- list(
  schema_version = 1L,
  fixture_version = cell$version,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  package_version = as.character(utils::packageVersion("crossform")),
  adapter_versions = c(
    fmridesign = as.character(utils::packageVersion("fmridesign")),
    fmrireg = as.character(utils::packageVersion("fmrireg"))
  ),
  dimensions = list(
    input_rows = unname(cell$counts),
    retained_rows = rep(cell$retained_count, length(cell$partitions)),
    partitions = length(cell$partitions),
    effects = length(cell$conditions$coordinates),
    neural_features = cell$study$observations$n_features,
    measurements = nrow(cell$frame$weights)
  ),
  identities = list(
    relation_plan_id = cell$plan$relation_plan_id,
    study_id = cell$study$study_id,
    design_model_id = cell$model$design_model_id,
    effect_map_id = cell$effects$effect_map_id,
    fixed_observation_model_id = cell$observation$observation_model_id
  ),
  timings = timings,
  checks = checks,
  gate = list(passed = gate)
)
if (!gate) stop("The first-moment vertical-slice gate failed.")
saveRDS(result, file.path(output_dir, "first-moment-vertical-slice.rds"),
  version = 3)
utils::write.csv(timings,
  file.path(output_dir, "first-moment-vertical-slice-summary.csv"),
  row.names = FALSE)
print(timings, row.names = FALSE)
cat("gate: PASS\n")
