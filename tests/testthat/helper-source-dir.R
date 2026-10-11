# Locating the package's R/ sources for the static-analysis tests.
#
# test-architecture.R and test-adapter-protocol.R parse R/*.R. A candidate
# directory counts only when it holds .R files: under covr (and any runner
# that tests an installed package) "../../R" resolves to the *installed*
# package's R/, which carries crossform.rdb/.rdx and no sources. Accepting it
# yields an empty call graph and an error instead of a skip. This mirrors the
# `.R`-file check in `.certification_source_root()` (helper-certification.R).

.has_r_sources <- function(dir) {
  dir.exists(dir) && length(list.files(dir, pattern = "\\.[Rr]$")) > 0L
}

package_source_dir <- function() {
  candidates <- c(
    "../../R",                       # testthat::test_local / devtools::test
    "../../../R",
    "../../00_pkg_src/crossform/R"   # R CMD check on a built tarball
  )
  for (candidate in candidates) {
    if (.has_r_sources(candidate)) return(normalizePath(candidate))
  }
  NULL
}
