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
  fixture_id = "rsatoolbox-fixed-linear-v1",
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
    "environment_lock",
    "algebraic_claim",
    "certification_copy",
    rep("fixture_contract", 3L),
    rep("recorded_output", 6L)
  ),
  path = c(
    "exemplars/rsatoolbox-parity/00-common.R",
    "exemplars/rsatoolbox-parity/01-fixture.R",
    "exemplars/rsatoolbox-parity/02-rsatoolbox.py",
    "exemplars/rsatoolbox-parity/03-compare.R",
    "exemplars/rsatoolbox-parity/requirements.txt",
    "design/common-geometry-equivalence.md",
    "inst/extdata/certification/common-geometry-external-parity.csv",
    "exemplars/rsatoolbox-parity/results/fixture-meta.csv",
    "exemplars/rsatoolbox-parity/results/model-rdms.csv",
    "exemplars/rsatoolbox-parity/results/regions.csv",
    "exemplars/rsatoolbox-parity/results/crossform-rdm.csv",
    "exemplars/rsatoolbox-parity/results/crossform-rsa.csv",
    "exemplars/rsatoolbox-parity/results/rsatoolbox-rdm.csv",
    "exemplars/rsatoolbox-parity/results/rsatoolbox-rsa.csv",
    "exemplars/rsatoolbox-parity/results/rsatoolbox-meta.csv",
    "exemplars/rsatoolbox-parity/results/agreement.csv"
  ),
  stringsAsFactors = FALSE
)

absolute <- file.path(repo, entries$path)
missing <- entries$path[!file.exists(absolute)]
if (length(missing)) {
  stop("Cannot create parity manifest; missing: ",
       paste(missing, collapse = ", "))
}

manifest <- data.frame(
  schema_version = 1L,
  fixture_id = "rsatoolbox-fixed-linear-v1",
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
