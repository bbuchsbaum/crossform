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
## shipped `.rds` under `inst/extdata/certification/` records a different
## digest. Artifacts that record no source digest at all (the shard-admission
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
cat("all ", length(shipped) - length(unbound),
  " digest-bound artifacts bind to the current tree\n", sep = "")
