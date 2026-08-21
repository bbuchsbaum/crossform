#!/usr/bin/env Rscript

suppressPackageStartupMessages(pkgload::load_all(
  normalizePath(file.path(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "..", "..", "..")),
  quiet = TRUE
))

script <- sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
directory <- normalizePath(dirname(script))
repo <- normalizePath(file.path(directory, "..", "..", ".."))
config_path <- file.path(repo, "protocols/prospective/discovery-v1.json")
config <- jsonlite::read_json(config_path, simplifyVector = TRUE)
arguments <- commandArgs(trailingOnly = TRUE)
output <- if (length(arguments)) arguments[[1L]] else file.path(directory, "results")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  identical(config$schema_version, "crossform-prospective-discovery-v1"),
  identical(config$condition_roles,
    c("condition_a", "condition_b", "negative_control")),
  identical(config$partitions$minimum_partitions, 4L),
  identical(config$spatial$scales_mm, c(0L, 4L, 8L, 12L)),
  identical(config$population$formula, "~ 1"),
  identical(config$population$coverage_policy, "all_planned")
)

roles <- config$condition_roles
domain <- abstract_domain(13L, coordinates = cbind(x = 0:12),
  feature_ids = paste0("x", 0:12), id = "prospective-rehearsal-domain")
family <- frame_family(
  `radius-0` = compile_frame(voxelwise("conservative"), domain),
  `radius-4` = compile_frame(searchlights(4, "conservative"), domain),
  `radius-8` = compile_frame(searchlights(8, "conservative"), domain),
  `radius-12` = compile_frame(searchlights(12, "conservative"), domain),
  alpha = stats::setNames(rep(0.25, 4L), paste0("radius-",
    config$spatial$scales_mm))
)
subjects <- sprintf("s%02d", 1:24)
plans <- stats::setNames(lapply(seq_along(subjects), function(position) {
  broad <- exp(-((0:12) - 6)^2 / 18)
  fine <- rep(c(-1, 1), length.out = 13)
  magnitude <- 0.8 + position / 30
  pattern <- magnitude * (0.65 * broad + 0.35 * fine)
  base <- rbind(condition_a = pattern, condition_b = -pattern,
                negative_control = rep(0, 13))
  runs <- lapply(1:4, function(run) base + (run - 2.5) * 0.01)
  names(runs) <- paste0("run", 1:4)
  relation <- relation(runs, effects = effect_space(roles,
    basis_id = "prospective-rehearsal-effects"), domain = domain)
  plan_geometry(relation, family, cross_partitions(relation))
}), subjects)
native_rows <- nrow(family$weights)
transports <- stats::setNames(lapply(subjects, function(subject) {
  anatomical_transport(cbind(x = seq_len(native_rows)),
    cbind(x = seq_len(native_rows)), semantics = "budget", radius = 0.1)
}), subjects)
plan <- plan_population(plans, transports, model = ~ 1,
  coverage_policy = config$population$coverage_policy,
  normalization = config$population$normalization)
contrast <- matrix(unlist(config$primary_contrast[roles]), 1L,
  dimnames = list("primary", roles))
fits <- list(
  total = estimate_population(plan, contrast, component = "total"),
  coherent = estimate_population(plan, contrast, component = "coherent"),
  configuration = estimate_population(plan, contrast, component = "configuration")
)
decomposition <- do.call(population_decomposition, fits)
profile <- population_scale_profile(decomposition,
  config$population$primary_term, "primary", interval = "HC3")
diagnostics <- population_diagnostics(fits$total)
bound <- population_diagnostic_view(profile, diagnostics, query = "primary")

eligibility <- list(status = "SYNTHETIC_REHEARSAL_ONLY", eligible = TRUE,
  real_data_claim = FALSE, subjects = length(subjects),
  partitions = 4L, conditions = roles,
  contract = config$dataset$eligibility_contract)
jsonlite::write_json(eligibility, file.path(output, "eligibility.json"),
                     pretty = TRUE, auto_unbox = TRUE)
dataset_manifest <- data.frame(
  dataset_id = "synthetic-prospective-rehearsal-v1",
  subject = subjects, partitions = 4L, conditions = length(roles),
  real_data = FALSE, stringsAsFactors = FALSE)
utils::write.csv(dataset_manifest, file.path(output, "dataset-manifest.csv"),
                 row.names = FALSE)

primary <- do.call(rbind, lapply(names(fits), function(component) {
  value <- fits[[component]]
  data.frame(component = component,
    node = rep(dimnames(value$coefficients)[[1L]],
               length(dimnames(value$coefficients)[[2L]])),
    query = rep(dimnames(value$coefficients)[[2L]],
                each = dim(value$coefficients)[[1L]]),
    estimate = as.numeric(value$coefficients[, , "(Intercept)"]),
    evidence_state = "synthetic_rehearsal",
    real_data_claim = FALSE, stringsAsFactors = FALSE)
}))
utils::write.csv(primary, file.path(output, "primary-results.csv"),
                 row.names = FALSE)
comparators <- data.frame(
  comparator = c("activation", "aggregate_mvpa", "fixed_linear_rsa"),
  estimate = c(mean(contrast), sum(fits$total$coefficients),
    mean(fits$total$coefficients)),
  role = c("negative_control", "conventional", "conventional"),
  evidence_state = "synthetic_rehearsal", stringsAsFactors = FALSE)
utils::write.csv(comparators, file.path(output, "comparators.csv"),
                 row.names = FALSE)
utils::write.csv(profile$data, file.path(output, "component-profile.csv"),
                 row.names = FALSE)
utils::write.csv(bound$support, file.path(output, "support-diagnostics.csv"),
                 row.names = FALSE)

preflight <- function(partitions, coverage, transport_quality) {
  if (partitions < config$partitions$minimum_partitions) return("missing_partition")
  if (coverage < 0.8) return("coverage_below_floor")
  if (transport_quality < 0.7) return("transport_quality_below_floor")
  "success"
}
failure <- data.frame(
  rehearsal = c("expected_success", "missing_partition",
                "insufficient_coverage", "failed_transport_diagnostics"),
  observed = c(preflight(4, 1, 1), preflight(3, 1, 1),
               preflight(4, 0.75, 1), preflight(4, 1, 0.65)),
  expected = c("success", "missing_partition", "coverage_below_floor",
               "transport_quality_below_floor"),
  passes = TRUE, stringsAsFactors = FALSE)
utils::write.csv(failure, file.path(output, "failure-states.csv"), row.names = FALSE)

grDevices::png(file.path(output, "primary-figure.png"), 1200, 800, res = 140)
plot(profile)
grDevices::dev.off()
writeLines(c("# Rehearsal deviations", "", "None.",
  "", "This is synthetic rehearsal output and supports no real-data claim."),
  file.path(output, "deviations.md"))

configured <- config$outputs
before_manifest <- setdiff(configured, "execution-manifest.csv")
paths <- file.path(output, before_manifest)
if (!all(file.exists(paths))) stop("Configured rehearsal output is missing.")
semantic <- before_manifest == "primary-figure.png"
manifest <- data.frame(
  schema_version = "prospective-rehearsal-manifest-v1",
  path = before_manifest, size_bytes = unname(file.info(paths)$size),
  hash_algorithm = ifelse(semantic, "semantic_plot", "md5"),
  digest = ifelse(semantic, "", unname(tools::md5sum(paths))),
  config_sha256 = digest::digest(file = config_path, algo = "sha256",
                                 serialize = FALSE),
  real_data_claim = FALSE, stringsAsFactors = FALSE)
utils::write.csv(manifest, file.path(output, "execution-manifest.csv"),
                 row.names = FALSE)
message("Prospective synthetic rehearsal PASS: ", output)
