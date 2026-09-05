# Compact, source-bound predictive evidence. Run from the repository root:
# Rscript benchmarks/predictive-geometry/certify.R
# This verifies existing full records; it never re-stamps stale measurements.
source("benchmarks/provenance.R")
root <- normalizePath(".")
current <- .crossform_source_tree_digest(root)
file_hash <- function(path) digest::digest(file = path, algo = "sha256")
inputs <- character()
read_bound <- function(path) {
  inputs <<- c(inputs, path)
  x <- readRDS(path)
  bound <- if (is.null(x$source_digest)) x$provenance$source_digest else x$source_digest
  if (!identical(bound, current)) stop("Stale input: ", path)
  hashes <- if (is.null(x$harness)) x$harness_digests else x$harness
  stopifnot(length(hashes) > 0L,
    identical(hashes, vapply(names(hashes), file_hash, "")))
  x
}
near <- function(x, y) stopifnot(isTRUE(all.equal(x, y, tolerance = 1e-12,
  check.attributes = FALSE)))
statistical <- function(kind) {
  path <- paste0("benchmark-results/predictive-geometry-", kind, ".rds")
  x <- read_bound(path)
  pilot_path <- sub(".rds", "-pilot.rds", path, fixed = TRUE)
  pilot <- readRDS(pilot_path); inputs <<- c(inputs, pilot_path)
  stopifnot(identical(x$mode, "run"), identical(x$config_hash, pilot$config_hash),
    identical(x$config_hash, digest::digest(x$config, algo = "sha256")),
    identical(x$harness, pilot$harness))
  for (i in seq_along(x$results)) {
    result <- x$results[[i]]; arm <- names(x$results)[[i]]
    stopifnot(result$count == pilot$counts[[arm]],
      result$seed == x$config$production_seed + 100L * i)
    rows <- x$summary[x$summary$arm == arm, , drop = FALSE]
    for (j in seq_len(nrow(rows))) {
      errors <- if (kind == "sampling") result$values[, "error", rows$predictor[j]] else
        result$values[, "error"]
      near(rows$error[j], mean(errors)); near(rows$mcse[j], sd(errors)/sqrt(length(errors)))
      near(rows$n[j], length(errors))
    }
  }
  s <- x$summary; c <- x$config
  stopifnot(all(abs(s$error) <= c$confidence_multiplier * s$mcse + c$numerical_tolerance),
    all(abs(s$error) + c$confidence_multiplier * s$mcse <=
      c$relative_equivalence_margin * s$reference_scale))
  null <- s[s$arm == "null", ]
  if (kind == "sampling") {
    g <- x$generator_moments
    near(g$error, g$observed - g$exact)
    stopifnot(all(abs(g$error) <= c$confidence_multiplier * g$mcse + c$numerical_tolerance),
      all(null$mean_leakage_error - c$confidence_multiplier * null$leakage_mcse >
        c$null_leakage_minimum * null$reference_scale))
  } else stopifnot(all(null$leakage_bias - c$confidence_multiplier * null$leakage_mcse >
    c$leakage_minimum * null$reference_scale))
  list(config = c, config_hash = x$config_hash, harness = x$harness,
    counts = pilot$counts, seeds = vapply(x$results, `[[`, 0L, "seed"),
    summary = s, generator_moments = x$generator_moments,
    paired_comparisons = x$paired_comparisons, runtime = x$runtime)
}
sampling <- statistical("sampling")
selection <- statistical("selection")

directory <- "benchmark-results/predictive-geometry-performance"
performance <- read_bound(file.path(directory, "predictive-performance-all.rds"))
raw_values <- list()
cases <- lapply(performance$case_artifacts, function(name) {
  x <- read_bound(file.path(directory, name))
  raw_values[[paste(x$case$id, x$route)]] <<- x$value
  stopifnot(identical(x$config, performance$config),
    all(x$reference_errors <= performance$config$relative_tolerance),
    x$admission$planned_workspace_bytes <= performance$config$workspace_bytes)
  row <- performance$summary[performance$summary$case == x$case$id &
    performance$summary$route == x$route, ]
  stopifnot(nrow(row) == 1L)
  near(row$seconds, x$elapsed_seconds)
  near(row$planned_bytes, x$admission$planned_workspace_bytes)
  near(row$incremental_rss_bytes, x$incremental_rss_bytes)
  near(row$largest_R_allocation, x$allocation$largest_R_allocation)
  near(row$retained_bytes, x$value$stored_bytes)
  near(row$read_calls, sum(x$source_reads$calls))
  near(row$read_columns, sum(x$source_reads$columns))
  near(row$form_error, x$reference_errors[["form"]])
  near(row$gain_error, x$reference_errors[["gain"]])
  reuse <- if (is.null(x$reuse)) NULL else {
    y <- x$reuse
    stopifnot(all(y$errors <= performance$config$relative_tolerance),
      y$cached$calls$materializations == 2L, y$uncached$calls$materializations == 12L)
    y
  }
  list(case = x$case, route = x$route, reference_errors = x$reference_errors,
    reference_rows = x$reference_rows, calls = x$calls, source_reads = x$source_reads,
    admission = x$admission, reuse = reuse, legacy_overlap = x$legacy_overlap,
    baseline_scope = x$baseline_scope, system = x$system)
})
names(cases) <- sub(".rds", "", performance$case_artifacts, fixed = TRUE)
stopifnot(nrow(performance$summary) == 7L, isTRUE(performance$parity))
matched_errors <- list()
for (id in c("tiny", "medium")) for (route in c("memory", "block")) {
  a <- raw_values[[paste(id, route)]]; b <- raw_values[[paste(id, "dense")]]
  errors <- vapply(c("forms", "spectra", "modes", "table", "diagnostics"), function(field) {
    actual <- a[[field]]; reference <- b[[field]]
    if (is.data.frame(actual)) {
      actual <- as.matrix(actual[, -1]); reference <- as.matrix(reference[, -1])
    }
    max(abs(actual - reference)) / max(1, max(abs(reference)))
  }, 0)
  stopifnot(all(errors <= performance$config$relative_tolerance))
  matched_errors[[paste(id, route)]] <- errors
}
rm(raw_values)
performance$cases <- cases
performance$matched_errors <- matched_errors
performance$case_artifacts <- NULL

directory <- "benchmark-results/predictive-geometry-mutations"
mutations <- read_bound(file.path(directory, "mutation-matrix-all.rds"))
stopifnot(nrow(mutations$matrix) == 20L, all(mutations$matrix$status == "KILLED"))
mutations$witnesses <- lapply(mutations$matrix$id, function(id) {
  x <- read_bound(file.path(directory, paste0(id, ".rds")))
  stopifnot(identical(x$test_digest, file_hash(x$test_file)),
    nrow(x$baseline$failures) == 0L, nrow(x$instrumented_control$failures) == 0L,
    x$baseline$expectations > 0L, all(x$instrumented_control$hits > 0L),
    all(x$mutated$hits > 0L), nrow(x$primary_failures) > 0L,
    all(grepl(x$primary, x$primary_failures$test, fixed = TRUE)))
  list(id = id, primary = x$primary, test_file = x$test_file, test_digest = x$test_digest,
    original_body_digests = x$original_body_digests,
    mutated_body_digests = vapply(x$mutated_bodies, digest::digest, "", algo = "sha256"),
    baseline_expectations = x$baseline$expectations, baseline_failures = nrow(x$baseline$failures),
    control_failures = nrow(x$instrumented_control$failures),
    control_hits = x$instrumented_control$hits, mutated_hits = x$mutated$hits,
    primary_failures = nrow(x$primary_failures))
})

oracle_file <- "design/oracles/predictive-geometry.R"
oracle_output <- system2(file.path(R.home("bin"), "Rscript"), oracle_file,
  stdout = TRUE, stderr = TRUE)
oracle_status <- attr(oracle_output, "status")
if (is.null(oracle_status)) oracle_status <- 0L
stopifnot(identical(oracle_status, 0L))
test_files <- c(list.files("tests/testthat", pattern = "^(test|helper)-predictive.*\\.R$",
  full.names = TRUE), "tests/testthat/test-architecture.R", "tests/testthat/test-api-surface.R")
harness_files <- c("benchmarks/predictive-geometry/certify.R", oracle_file,
  "benchmarks/predictive-geometry/restore-counts.R",
  "benchmarks/check-certification-binding.R", "benchmarks/promote-artifacts.R",
  "benchmarks/admission-coverage.R", test_files)
record <- list(schema = "predictive-certification-v1",
  provenance = crossform_benchmark_provenance(root, "predictive-geometry/certify.R"),
  input_hashes = vapply(unique(inputs), file_hash, ""),
  harness_digests = vapply(harness_files, file_hash, ""),
  sampling = sampling, selection = selection, performance = performance, mutations = mutations,
  oracle = list(file = oracle_file, digest = file_hash(oracle_file),
    exit_code = oracle_status, output = oracle_output),
  limits = c("same-condition independent declared origins; no automatic inference",
    "workspace_bytes bounds planned numerical buffers, not process RSS",
    "timings describe recorded workloads, not universal speedup"))
stopifnot(identical(current, .crossform_source_tree_digest(root)))
path <- "benchmark-results/predictive-geometry-validation.rds"
saveRDS(record, path, version = 2L, compress = "xz")
stopifnot(file.info(path)$size <= 64 * 1024)
cat("Verified and compacted predictive evidence for", current, "\n",
  path, file.info(path)$size, "bytes\n")
