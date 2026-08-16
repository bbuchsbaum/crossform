# Atomic integrity-checked checkpoints --------------------------------------

.checkpoint_schema_version <- 1L

.available_disk_bytes <- function(path) {
  lines <- tryCatch(
    system2("df", c("-Pk", shQuote(path)), stdout = TRUE, stderr = FALSE),
    error = function(error) character()
  )
  if (length(lines) < 2L) return(NA_real_)
  fields <- strsplit(trimws(lines[[length(lines)]]), "[[:space:]]+")[[1L]]
  value <- suppressWarnings(as.numeric(fields[[4L]]))
  if (is.finite(value)) value * 1024 else NA_real_
}

.checkpoint_identity <- function(receipt) {
  list(
    scientific_plan_id = receipt$scientific_plan_id,
    kernel_version = receipt$kernel_version,
    task_partition_id = receipt$task_partition_id,
    reduction_plan_id = receipt$reduction_plan_id,
    precision = receipt$precision,
    source_revisions = vapply(receipt$sources, `[[`, character(1), "stable_revision")
  )
}

.checkpoint_stop <- function(message, receipt, cleanup_success = TRUE) {
  stop(.execution_error(
    message,
    receipt = receipt,
    cleanup_status = list(attempted = TRUE, success = cleanup_success, message = NULL)
  ))
}

.write_checkpoint <- function(payload, path, receipt, available_bytes = NULL) {
  if (!inherits(receipt, "effect_execution_receipt")) {
    stop("`receipt` must be a crossform execution receipt.", call. = FALSE)
  }
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be one nonempty checkpoint directory.", call. = FALSE)
  }
  if (file.exists(path)) .checkpoint_stop("Refusing to overwrite a checkpoint.", receipt)
  parent <- dirname(path)
  if (!dir.exists(parent)) stop("Checkpoint parent directory does not exist.", call. = FALSE)
  required <- as.numeric(utils::object.size(payload)) * 3 + 65536
  if (is.null(available_bytes)) available_bytes <- .available_disk_bytes(parent)
  if (!is.numeric(available_bytes) || length(available_bytes) != 1L ||
      is.na(available_bytes) || !is.finite(available_bytes)) {
    .checkpoint_stop("Unable to establish available checkpoint disk space.", receipt)
  }
  if (available_bytes < required) {
    .checkpoint_stop("Insufficient disk space for checkpoint staging.", receipt)
  }

  temporary <- tempfile(pattern = paste0(".", basename(path), ".tmp-"), tmpdir = parent)
  if (!dir.create(temporary, mode = "0700")) {
    .checkpoint_stop("Could not create checkpoint staging directory.", receipt, FALSE)
  }
  committed <- FALSE
  on.exit(if (!committed && dir.exists(temporary)) unlink(temporary, recursive = TRUE), add = TRUE)
  payload_path <- file.path(temporary, "payload.rds")
  manifest_path <- file.path(temporary, "manifest.rds")
  saveRDS(payload, payload_path, version = 3)
  Sys.chmod(payload_path, mode = "0600", use_umask = FALSE)
  checksum <- unname(tools::md5sum(payload_path))
  manifest <- list(
    schema_version = .checkpoint_schema_version,
    checksum_algorithm = "md5",
    payload_checksum = checksum,
    payload_bytes = file.info(payload_path)$size,
    identity = .checkpoint_identity(receipt)
  )
  saveRDS(manifest, manifest_path, version = 3)
  Sys.chmod(manifest_path, mode = "0600", use_umask = FALSE)
  Sys.chmod(temporary, mode = "0700", use_umask = FALSE)
  if (!file.rename(temporary, path)) {
    .checkpoint_stop("Atomic checkpoint commit failed.", receipt, FALSE)
  }
  committed <- TRUE
  invisible(path)
}

.read_checkpoint <- function(path, expected_receipt) {
  if (!inherits(expected_receipt, "effect_execution_receipt")) {
    stop("`expected_receipt` must be a crossform execution receipt.", call. = FALSE)
  }
  manifest_path <- file.path(path, "manifest.rds")
  payload_path <- file.path(path, "payload.rds")
  if (!dir.exists(path) || !file.exists(manifest_path) || !file.exists(payload_path)) {
    .checkpoint_stop("Checkpoint is incomplete.", expected_receipt)
  }
  manifest <- tryCatch(readRDS(manifest_path), error = function(error) NULL)
  if (!is.list(manifest) || !identical(manifest$schema_version, .checkpoint_schema_version)) {
    .checkpoint_stop("Checkpoint schema version mismatch.", expected_receipt)
  }
  if (!identical(manifest$identity, .checkpoint_identity(expected_receipt))) {
    .checkpoint_stop("Checkpoint execution identity does not match resume plan.",
      expected_receipt)
  }
  checksum <- unname(tools::md5sum(payload_path))
  if (!identical(checksum, manifest$payload_checksum) ||
      !identical(file.info(payload_path)$size, manifest$payload_bytes)) {
    .checkpoint_stop("Checkpoint payload checksum or size mismatch.", expected_receipt)
  }
  if (.Platform$OS.type == "unix") {
    modes <- as.integer(file.info(c(path, manifest_path, payload_path))$mode)
    if (any(bitwAnd(modes, strtoi("077", base = 8L)) != 0L)) {
      .checkpoint_stop("Checkpoint permissions are not restrictive.", expected_receipt)
    }
  }
  tryCatch(
    readRDS(payload_path),
    error = function(error) .checkpoint_stop("Checkpoint payload cannot be decoded.",
      expected_receipt)
  )
}
