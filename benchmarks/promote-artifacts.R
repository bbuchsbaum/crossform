#!/usr/bin/env Rscript

# Promote freshly recorded certification artifacts into the package.
#
# Runners write to `benchmark-results/`, which is not shipped. The small gate
# artifacts that the test suite reads as evidence live in
# `inst/extdata/certification/`, so that certification runs under `R CMD check`
# instead of skipping. This script is the only supported way to move one into
# the other.
#
#   Rscript benchmarks/promote-artifacts.R [results-dir] [repo]
#
# Promotion refuses an artifact that carries no provenance: an unbound
# artifact cannot be tied to the source that produced it, so it is not
# evidence. Large validation records stay out of the tarball by design and are
# never promoted; the tests that read them skip when they are absent.

arguments <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(arguments) >= 1L) arguments[[1L]] else
  "benchmark-results"
repo <- if (length(arguments) >= 2L) {
  normalizePath(arguments[[2L]], mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
destination <- file.path(repo, "inst", "extdata", "certification")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)

# Artifact name -> runner that re-certifies it. Derived from the admission
# coverage list; anything absent from that list is never promoted.
source(file.path(repo, "benchmarks", "admission-coverage.R"), local = TRUE)
promotable <- .crossform_promotable_artifacts()

# Compact tables promoted alongside their artifact. These are human-readable
# receipts, not test evidence.
promotable_summaries <- c(
  "public-map-scale-gate-summary.csv",
  "query-first-scale-gate-summary.csv",
  "sampling-covariance-scale-summary.csv",
  "first-moment-vertical-slice-summary.csv",
  "crossnobis-scale-gate-summary.csv",
  "shard-admission-summary.csv",
  "memory-benchmark-summary.csv",
  "learned-metric-policy-summary.csv",
  "sampling-covariance-validation-summary.csv"
)

maximum_shipped_bytes <- 64 * 1024

promote <- function(name, runner) {
  source_path <- file.path(results_dir, name)
  if (!file.exists(source_path)) {
    return(sprintf("skip   %-38s (not present in %s)", name, results_dir))
  }
  size <- file.info(source_path)$size
  if (size > maximum_shipped_bytes) {
    return(sprintf(
      "REFUSE %-38s (%.1f KiB exceeds the %.0f KiB shipped cap)",
      name, size / 1024, maximum_shipped_bytes / 1024
    ))
  }
  if (!is.na(runner)) {
    artifact <- tryCatch(readRDS(source_path), error = function(condition) NULL)
    provenance <- if (is.list(artifact)) artifact$provenance else NULL
    digest <- if (is.list(provenance)) provenance$source_digest else NULL
    # A recorded digest of NA is as unbound as a missing one.
    if (!is.character(digest) || length(digest) != 1L || is.na(digest) ||
        !nzchar(digest)) {
      return(sprintf(
        "REFUSE %-38s (no source digest: re-run benchmarks/%s)", name, runner
      ))
    }
  }
  target <- file.path(destination, name)
  if (grepl("\\.rds$", name)) {
    # Shipped artifacts are re-serialized at serialization version 2 so that
    # `R CMD build` does not add an implicit `R (>= 3.5.0)` dependency to a
    # package that has no other reason for one. The object is unchanged.
    saveRDS(readRDS(source_path), target, version = 2L, compress = "xz")
  } else if (!file.copy(source_path, target, overwrite = TRUE)) {
    return(sprintf("FAILED %-38s (copy failed)", name))
  }
  sprintf("promote %-37s (%.1f KiB)", name, file.info(target)$size / 1024)
}

report <- c(
  vapply(names(promotable), function(name) promote(name, promotable[[name]]),
    character(1), USE.NAMES = FALSE),
  vapply(promotable_summaries, function(name) promote(name, NA_character_),
    character(1), USE.NAMES = FALSE)
)
cat(report, sep = "\n")
cat("\n")
if (any(grepl("^(REFUSE|FAILED)", report))) {
  stop("Promotion refused at least one artifact; see the report above.",
    call. = FALSE)
}
