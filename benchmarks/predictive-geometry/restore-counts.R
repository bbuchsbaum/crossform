# Restore the frozen pilot-derived count contract in a fresh checkout.
# Rscript benchmarks/predictive-geometry/restore-counts.R [repo] [results-dir]
# This is metadata recovery, not a pilot run or new statistical evidence.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
output <- if (length(args) >= 2L) args[2] else file.path(root, "benchmark-results")
receipt_path <- file.path(root, "inst/extdata/certification/predictive-geometry-validation.rds")
receipt <- readRDS(receipt_path)
stopifnot(identical(receipt$schema, "predictive-certification-v1"))
contracts <- lapply(c("sampling", "selection"), function(kind) {
  x <- receipt[[kind]]
  stopifnot(identical(x$config_hash, digest::digest(x$config, algo = "sha256")),
    identical(x$harness, vapply(names(x$harness), function(path)
      digest::digest(file = file.path(root, path), algo = "sha256"), "")),
    identical(names(x$counts), x$config$arms),
    all(is.finite(x$counts) & x$counts >= x$config$minimum_replicates),
    all(x$counts %% x$config$count_rounding == 0))
  list(mode = "restored_pilot_count_contract", counts = x$counts,
    config = x$config, config_hash = x$config_hash, harness = x$harness,
    restored_from = list(receipt = receipt_path,
      receipt_digest = digest::digest(file = receipt_path, algo = "sha256"),
      source_digest = receipt$provenance$source_digest),
    interpretation = "Frozen counts recovered from prior certification; no pilot data were generated.")
})
names(contracts) <- c("sampling", "selection")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
for (kind in names(contracts)) {
  path <- file.path(output, paste0("predictive-geometry-", kind, "-pilot.rds"))
  if (file.exists(path)) {
    old <- readRDS(path)
    stopifnot(identical(old$counts, contracts[[kind]]$counts),
      identical(old$config_hash, contracts[[kind]]$config_hash),
      identical(old$harness, contracts[[kind]]$harness))
    cat("Retained existing matching count contract:", path, "\n")
  } else {
    saveRDS(contracts[[kind]], path, version = 2L)
    cat("Restored frozen count contract:", path, "\n")
  }
}
