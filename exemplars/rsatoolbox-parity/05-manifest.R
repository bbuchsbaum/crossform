#!/usr/bin/env Rscript
# 05-manifest.R -- bind the recorded parity outputs to their producing sources.
#
# The manifest deliberately has no timestamp: the same bytes produce the same
# versioned record. It hashes the external implementation, the R fixture and
# comparison sources, the environment lock, the algebraic claim, and every
# tracked artifact needed to revalidate the mapped parity case.

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  normalizePath(Sys.getenv("EXEMPLAR_DIR"))
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
repo <- normalizePath(file.path(exemplar_dir, "..", ".."))

# Reuse benchmarks/provenance.R rather than reimplementing the digest: the
# canonical definition folds each file's *name* into the hash and sorts with
# method = "radix" so the result does not depend on collation locale. A
# private reimplementation that got either detail wrong would produce a digest
# that never matches the one CI checks, which is worse than no binding at all.
provenance_env <- new.env(parent = globalenv())
sys.source(file.path(repo, "benchmarks", "provenance.R"), envir = provenance_env)
source_digest <- provenance_env$.crossform_source_tree_digest(repo)
git_state <- provenance_env$.crossform_git_provenance(repo)
git_commit <- git_state$git_commit
git_dirty <- git_state$git_dirty

# The source exemplar is excluded from the package payload. Publish the small
# comparison receipt under inst/ so the executable article can read the same
# generated numbers from an installed package without copying them by hand.
agreement <- utils::read.csv(
  file.path(exemplar_dir, "results", "agreement.csv"),
  stringsAsFactors = FALSE
)
external_meta <- utils::read.csv(
  file.path(exemplar_dir, "results", "rsatoolbox-meta.csv"),
  stringsAsFactors = FALSE
)
external_meta <- stats::setNames(external_meta$value, external_meta$key)
certification <- data.frame(
  fixture_id = "rsatoolbox-standard-workflow-v2",
  crossform_source_digest = source_digest,
  crossform_git_commit = git_commit,
  python_version = external_meta[["python"]],
  rsatoolbox_version = external_meta[["rsatoolbox"]],
  numpy_version = external_meta[["numpy"]],
  scipy_version = external_meta[["scipy"]],
  agreement,
  stringsAsFactors = FALSE
)
certification_path <- file.path(
  repo, "inst", "extdata", "certification",
  "common-geometry-external-parity.csv"
)
dir.create(dirname(certification_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(certification, certification_path, row.names = FALSE)

entries <- data.frame(
  role = c(
    rep("fixture_source", 2L),
    "external_implementation",
    "comparison_source",
    "extension_source",
    "environment_lock",
    "algebraic_claim",
    "certification_copy",
    rep("fixture_contract", 4L),
    rep("recorded_output", 15L),
    rep("extension_output", 4L)
  ),
  path = c(
    "exemplars/rsatoolbox-parity/00-common.R",
    "exemplars/rsatoolbox-parity/01-fixture.R",
    "exemplars/rsatoolbox-parity/02-rsatoolbox.py",
    "exemplars/rsatoolbox-parity/03-compare.R",
    "exemplars/rsatoolbox-parity/04-extension.R",
    "exemplars/rsatoolbox-parity/requirements.txt",
    "design/common-geometry-equivalence.md",
    "inst/extdata/certification/common-geometry-external-parity.csv",
    "exemplars/rsatoolbox-parity/results/fixture-meta.csv",
    "exemplars/rsatoolbox-parity/results/model-rdms.csv",
    "exemplars/rsatoolbox-parity/results/model-rdms-extended.csv",
    "exemplars/rsatoolbox-parity/results/regions.csv",
    "exemplars/rsatoolbox-parity/results/crossform-rdm.csv",
    "exemplars/rsatoolbox-parity/results/crossform-rsa.csv",
    "exemplars/rsatoolbox-parity/results/crossform-distances.csv",
    "exemplars/rsatoolbox-parity/results/crossform-correlation.csv",
    "exemplars/rsatoolbox-parity/results/crossform-similarity.csv",
    "exemplars/rsatoolbox-parity/results/grand-means.csv",
    "exemplars/rsatoolbox-parity/results/precision-diagonal.csv",
    "exemplars/rsatoolbox-parity/results/rsatoolbox-rdm.csv",
    "exemplars/rsatoolbox-parity/results/rsatoolbox-rsa.csv",
    "exemplars/rsatoolbox-parity/results/rsatoolbox-distances.csv",
    "exemplars/rsatoolbox-parity/results/rsatoolbox-similarity.csv",
    "exemplars/rsatoolbox-parity/results/rsatoolbox-meta.csv",
    "exemplars/rsatoolbox-parity/results/agreement.csv",
    "exemplars/rsatoolbox-parity/results/coverage.csv",
    "exemplars/rsatoolbox-parity/results/comparison-meta.csv",
    "exemplars/rsatoolbox-parity/results/extension.csv",
    "exemplars/rsatoolbox-parity/results/extension-table.csv",
    "exemplars/rsatoolbox-parity/results/extension-uncertainty.csv",
    "exemplars/rsatoolbox-parity/results/extension-refusals.csv"
  ),
  stringsAsFactors = FALSE
)

absolute <- file.path(repo, entries$path)
missing <- entries$path[!file.exists(absolute)]
if (length(missing)) {
  stop("Cannot create parity manifest; missing: ",
       paste(missing, collapse = ", "))
}

# Bind the record to the crossform source that produced the R arm. Without
# this the manifest can detect a later edit to an output but not the fact that
# the outputs were generated by a different package source than the one in the
# tree. Same digest definition as benchmarks/provenance.R.
manifest <- data.frame(
  schema_version = 2L,
  crossform_source_digest = source_digest,
  crossform_git_commit = git_commit,
  crossform_git_dirty = git_dirty,
  fixture_id = "rsatoolbox-standard-workflow-v2",
  hash_algorithm = "md5",
  role = entries$role,
  path = entries$path,
  size_bytes = unname(file.info(absolute)$size),
  digest = unname(tools::md5sum(absolute)),
  stringsAsFactors = FALSE
)

output <- file.path(exemplar_dir, "results", "parity-manifest.csv")
utils::write.csv(manifest, output, row.names = FALSE)
message("Wrote ", nrow(manifest), " source/artifact bindings to ", output)
