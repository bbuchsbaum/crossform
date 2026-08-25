#!/usr/bin/env Rscript
## Refuse a tree whose shipped certification artifacts do not bind to it.
##
##   Rscript benchmarks/check-certification-binding.R [repo-root]
##
## The test suite treats a stale artifact as a loud SKIP so that local work is
## never blocked by yesterday's benchmarks. CI wants the opposite polarity: a
## pull request that edits `R/` without re-running the runners and
## `benchmarks/promote-artifacts.R` must FAIL, or stale certification merges
## and the recorded evidence quietly stops being evidence of anything.
##
## This script is that gate. It recomputes the aggregate source digest with
## `benchmarks/provenance.R` — the same digest the artifacts record — and
## exits nonzero, naming every artifact and its re-certifying runner, when any
## shipped `.rds` -- or the shipped external-parity `.csv` receipt -- under
## `inst/extdata/certification/` records a different digest. Artifacts that record no source digest at all (the shard-admission
## record, whose unbound state is designed and documented) are reported and
## tolerated.
##
## Wired into CI by .github/workflows/certification-binding.yaml.

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else "."
root <- normalizePath(root, mustWork = TRUE)

environment <- new.env(parent = globalenv())
sys.source(file.path(root, "benchmarks", "provenance.R"), envir = environment)
current <- environment$.crossform_source_tree_digest(root)
cat("current R/ digest: ", current, "\n", sep = "")

shipped <- list.files(file.path(root, "inst", "extdata", "certification"),
  pattern = "\\.rds$", full.names = TRUE)
if (!length(shipped)) {
  stop("no shipped certification artifacts found; nothing to bind")
}

stale <- character()
unbound <- character()
for (path in shipped) {
  artifact <- readRDS(path)
  recorded <- artifact$provenance$source_digest
  name <- basename(path)
  if (is.null(recorded)) {
    unbound <- c(unbound, name)
  } else if (!identical(recorded, current)) {
    runner <- artifact$provenance$runner
    stale <- c(stale, sprintf("%s (recorded %s; re-run %s)", name,
      substr(sub("^sha256:", "", recorded), 1L, 12L),
      if (is.null(runner)) "its runner" else runner))
  }
}

## The external-parity receipt is a CSV, not an RDS, because it is produced by
## a two-language pipeline rather than an R runner. It was therefore invisible
## to the `\\.rds$` glob above and was the one shipped certification artifact
## with no binding at all. It records the digest in a column instead of a
## provenance list.
## Only standalone receipts are listed here. The other shipped CSVs are
## human-readable companions to an `.rds` that already binds, so checking them
## would double-count; this list is the set that has no `.rds` behind it. A
## receipt named here that loses its digest column is a hard failure, not a
## tolerated `unbound`, because silence is exactly how this gap arose.
csv_receipts <- file.path(root, "inst", "extdata", "certification",
  c("common-geometry-external-parity.csv"))
for (path in csv_receipts) {
  name <- basename(path)
  if (!file.exists(path)) {
    stale <- c(stale, sprintf("%s (missing)", name))
    next
  }
  receipt <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!"crossform_source_digest" %in% names(receipt)) {
    stale <- c(stale, sprintf(
      "%s (no crossform_source_digest column; regenerate with 05-manifest.R)",
      name))
    next
  }
  recorded <- unique(receipt$crossform_source_digest)
  if (length(recorded) != 1L) {
    stale <- c(stale, sprintf("%s (mixed source digests in one receipt)", name))
  } else if (!identical(recorded, current)) {
    stale <- c(stale, sprintf(
      "%s (recorded %s; re-run exemplars/rsatoolbox-parity/run-all.sh)", name,
      substr(sub("^sha256:", "", recorded), 1L, 12L)))
  }
}

if (length(unbound)) {
  cat("unbound by design (no recorded digest):\n",
    paste0("  - ", unbound, collapse = "\n"), "\n", sep = "")
}
if (length(stale)) {
  cat("STALE certification artifacts:\n",
    paste0("  - ", stale, collapse = "\n"), "\n", sep = "")
  cat("\nRe-certify on a frozen tree per benchmarks/RECERTIFY.md, then\n",
    "promote with benchmarks/promote-artifacts.R.\n", sep = "")
  quit(status = 1L)
}
cat("all ", length(shipped) - length(unbound) + length(csv_receipts),
  " digest-bound artifacts bind to the current tree\n", sep = "")
