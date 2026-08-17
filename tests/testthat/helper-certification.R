# Binding recorded certification artifacts to the source that produced them.
#
# A recorded benchmark artifact is evidence only for the source tree that
# produced it. These helpers refuse to read a recorded gate as a boolean
# unless its provenance still binds it to the code under test.
#
# Two binding tiers, in decreasing strength:
#
# 1. `source_tree`  - the checkout is available, so the aggregate SHA-256 of
#    `R/*.R` recorded in the artifact must equal the current one.
# 2. `package_version` - under `R CMD check` the package sources are not
#    shipped, so only the recorded package version can be compared. The tier
#    is reported so that a version-level binding is never mistaken for a
#    source-level one.
#
# In every non-binding case the test skips with an instruction naming the
# runner that re-certifies it. A skip is loud and greppable; a silent pass on
# a stale artifact is the failure mode this replaces.

.certification_cache <- new.env(parent = emptyenv())

.certification_source_root <- function() {
  root <- testthat::test_path("..", "..")
  sources <- file.path(root, "R")
  helper <- file.path(root, "benchmarks", "provenance.R")
  if (!dir.exists(sources) || !file.exists(helper) ||
      !file.exists(file.path(root, "DESCRIPTION"))) {
    return(NA_character_)
  }
  if (!length(list.files(sources, pattern = "\\.[Rr]$"))) {
    return(NA_character_)
  }
  root
}

.certification_current_digest <- function() {
  if (!is.null(.certification_cache$digest)) {
    return(.certification_cache$digest)
  }
  root <- .certification_source_root()
  value <- if (is.na(root)) {
    NA_character_
  } else {
    environment <- new.env(parent = globalenv())
    sys.source(file.path(root, "benchmarks", "provenance.R"),
      envir = environment)
    environment$.crossform_source_tree_digest(root)
  }
  .certification_cache$digest <- value
  value
}

.certification_short <- function(digest) {
  if (!is.character(digest) || !length(digest) || is.na(digest[[1L]])) {
    return("<none>")
  }
  substr(sub("^sha256:", "", digest[[1L]]), 1L, 12L)
}

#' Locate a recorded certification artifact
#'
#' Shipped gate artifacts live in `inst/extdata/certification`, so they are
#' available under `R CMD check`. A fresh run written to `benchmark-results`
#' takes precedence when it exists, which is how a re-certification run is
#' picked up before it is promoted.
certification_artifact_path <- function(name) {
  fresh <- testthat::test_path("..", "..", "benchmark-results", name)
  if (file.exists(fresh)) {
    return(fresh)
  }
  installed <- tryCatch(
    system.file("extdata", "certification", name, package = "crossform"),
    error = function(condition) ""
  )
  if (nzchar(installed) && file.exists(installed)) {
    return(installed)
  }
  shipped <- testthat::test_path(
    "..", "..", "inst", "extdata", "certification", name
  )
  if (file.exists(shipped)) {
    return(shipped)
  }
  NA_character_
}

.certification_recorded_digest <- function(provenance) {
  digest <- provenance$source_digest
  if (!is.character(digest) || length(digest) != 1L || is.na(digest) ||
      !nzchar(digest)) {
    return(NA_character_)
  }
  digest
}

#' Refuse to read an artifact that is not bound to the current source
#'
#' @return Invisibly, the binding tier that was satisfied.
skip_unless_certified <- function(artifact, runner) {
  instruction <- paste0("re-run benchmarks/", runner,
    " and benchmarks/promote-artifacts.R to re-certify.")
  provenance <- if (is.list(artifact)) artifact$provenance else NULL
  recorded_digest <- if (is.list(provenance)) {
    .certification_recorded_digest(provenance)
  } else {
    NA_character_
  }
  # A recorded digest of NA is as unbound as a missing one: the runner could
  # not see an `R/` tree, so nothing ties the artifact to any source.
  if (!is.list(provenance) || is.na(recorded_digest)) {
    testthat::skip(paste(
      "CERTIFICATION UNBOUND: the recorded artifact carries no source digest,",
      "so it cannot be tied to any source tree;", instruction
    ))
  }
  current <- .certification_current_digest()
  if (!is.na(current)) {
    if (!identical(current, recorded_digest)) {
      testthat::skip(paste0(
        "CERTIFICATION STALE: the recorded artifact predates the current ",
        "source (recorded R/ digest ",
        .certification_short(recorded_digest), ", current ",
        .certification_short(current), "); ", instruction
      ))
    }
    return(invisible("source_tree"))
  }
  recorded <- as.character(provenance$package_version)
  installed <- tryCatch(
    as.character(utils::packageVersion("crossform")),
    error = function(condition) NA_character_
  )
  if (length(recorded) != 1L || is.na(recorded) || is.na(installed) ||
      !identical(recorded, installed)) {
    testthat::skip(paste0(
      "CERTIFICATION UNBOUND: package sources are unavailable here and the ",
      "recorded package version (", recorded, ") is not the installed one (",
      installed, "); ", instruction
    ))
  }
  # Say out loud that this is the weaker tier. A source edit that does not
  # bump the version leaves this comparison satisfied, so the pass below is
  # evidence about a version, not about the code that is running.
  message("CERTIFICATION TIER package_version (source unavailable): ",
    provenance$runner, " recorded against crossform ", recorded, ".")
  invisible("package_version")
}

#' Read a recorded certification artifact bound to the current source
#'
#' Skips, never silently passes, when the artifact is absent or unbound.
certified_artifact <- function(name, runner) {
  path <- certification_artifact_path(name)
  if (is.na(path)) {
    testthat::skip(paste0(
      "CERTIFICATION ABSENT: no recorded ", name, " is available here; ",
      "run benchmarks/", runner, " to produce one."
    ))
  }
  artifact <- readRDS(path)
  skip_unless_certified(artifact, runner)
  artifact
}
