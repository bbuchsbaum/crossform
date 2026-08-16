# Provenance binding for recorded certification artifacts.
#
# Every benchmarks/run-*.R runner records a `provenance` list inside the
# artifact it writes. Certification tests refuse to treat a recorded artifact
# as evidence unless that provenance still binds it to the source that
# produced it. See benchmarks/README.md for the two binding tiers.
#
# This file is sourced by the runners and, when a source checkout is
# available, by tests/testthat/helper-certification.R. It must stay free of
# package dependencies other than `digest`, which crossform already imports.

.crossform_source_tree_digest <- function(root = ".") {
  directory <- file.path(root, "R")
  if (!dir.exists(directory)) {
    return(NA_character_)
  }
  # Radix ordering keeps the digest independent of the collation locale.
  files <- sort(list.files(directory, pattern = "\\.[Rr]$"), method = "radix")
  if (!length(files)) {
    return(NA_character_)
  }
  entries <- vapply(files, function(file) {
    paste0(
      file, " ",
      digest::digest(file = file.path(directory, file), algo = "sha256")
    )
  }, character(1), USE.NAMES = FALSE)
  paste0("sha256:", digest::digest(
    paste(entries, collapse = "\n"),
    algo = "sha256", serialize = FALSE
  ))
}

.crossform_git_field <- function(root, arguments) {
  value <- tryCatch(
    suppressWarnings(system2(
      "git", c("-C", shQuote(root), arguments),
      stdout = TRUE, stderr = FALSE
    )),
    error = function(condition) character()
  )
  status <- attr(value, "status")
  if (!is.null(status) && !identical(as.integer(status), 0L)) {
    return(NA_character_)
  }
  if (!length(value)) {
    return(NA_character_)
  }
  paste(value, collapse = "\n")
}

.crossform_git_provenance <- function(root = ".") {
  commit <- .crossform_git_field(root, c("rev-parse", "HEAD"))
  if (is.na(commit)) {
    return(list(git_commit = NA_character_, git_dirty = NA))
  }
  status <- .crossform_git_field(root, c("status", "--porcelain", "--", "R"))
  list(
    git_commit = commit,
    git_dirty = if (is.na(status)) NA else nzchar(status)
  )
}

.crossform_description_version <- function(root = ".") {
  path <- file.path(root, "DESCRIPTION")
  if (!file.exists(path)) {
    return(NA_character_)
  }
  fields <- tryCatch(read.dcf(path, fields = "Version"),
    error = function(condition) NULL)
  if (is.null(fields) || is.na(fields[[1L, 1L]])) {
    return(NA_character_)
  }
  as.character(fields[[1L, 1L]])
}

#' Record the provenance of a certification artifact
#'
#' @param root Repository root of the source checkout being certified.
#' @param runner File name of the runner writing the artifact, used verbatim
#'   in the re-certification instruction a failing test prints.
#' @return A list stored as the `provenance` element of the artifact.
crossform_benchmark_provenance <- function(root = ".", runner = NA_character_) {
  root <- tryCatch(normalizePath(root, mustWork = TRUE),
    error = function(condition) root)
  git <- .crossform_git_provenance(root)
  version <- tryCatch(
    as.character(utils::packageVersion("crossform")),
    error = function(condition) .crossform_description_version(root)
  )
  blas <- tryCatch(unname(extSoftVersion()[["BLAS"]]),
    error = function(condition) NA_character_)
  lapack <- tryCatch(La_library(), error = function(condition) NA_character_)
  lapack_version <- tryCatch(La_version(),
    error = function(condition) NA_character_)
  list(
    schema = 1L,
    runner = runner,
    recorded_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package = "crossform",
    package_version = version,
    source_digest = .crossform_source_tree_digest(root),
    git_commit = git$git_commit,
    git_dirty = git$git_dirty,
    r_version = R.version.string,
    platform = R.version$platform,
    blas = blas,
    lapack = lapack,
    lapack_version = lapack_version
  )
}
