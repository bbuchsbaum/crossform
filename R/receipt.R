# Source and execution provenance -------------------------------------------

#' Declare source execution capabilities
#'
#' @param block_read Whether bounded feature-block reads are supported.
#' @param reopenable Whether a fresh read-only handle can be opened safely.
#' @param thread_safe Whether concurrent reads through one handle are supported.
#' @param stable_revision A strong immutable source revision or checksum.
#' @return A declarative source-capability value.
#' @export
source_capabilities <- function(block_read, reopenable = FALSE,
                                thread_safe = FALSE, stable_revision) {
  flags <- list(
    block_read = block_read,
    reopenable = reopenable,
    thread_safe = thread_safe
  )
  if (!all(vapply(flags, function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }, logical(1)))) {
    stop("Source capability flags must each be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(stable_revision) || length(stable_revision) != 1L ||
      is.na(stable_revision) ||
      !grepl("^sha256:[[:xdigit:]]{64}$", stable_revision)) {
    stop("`stable_revision` must be a sha256 identifier with 64 hexadecimal digits.",
      call. = FALSE)
  }
  value <- structure(
    c(flags, list(stable_revision = stable_revision)),
    class = "effect_source_capabilities"
  )
  .validate_source_capabilities(value)
}

.validate_source_capabilities <- function(value) {
  expected <- c("block_read", "reopenable", "thread_safe", "stable_revision")
  if (!inherits(value, "effect_source_capabilities") ||
      !identical(names(value), expected)) {
    stop("Source capabilities are missing required canonical fields.", call. = FALSE)
  }
  flags <- value[c("block_read", "reopenable", "thread_safe")]
  if (!all(vapply(flags, function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }, logical(1)))) {
    stop("Source capability flags must each be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(value$stable_revision) || length(value$stable_revision) != 1L ||
      is.na(value$stable_revision) ||
      !grepl("^sha256:[[:xdigit:]]{64}$", value$stable_revision)) {
    stop("Source capability revision must be one strong sha256 identifier.",
      call. = FALSE)
  }
  structure(as.list(value), class = "effect_source_capabilities")
}

#' Record the execution that produced a geometry or view
#'
#' A receipt separates scientific-plan identity, numerical execution policy,
#' source revisions, and reduction identity. Progress reporters and other
#' nonsemantic observers are intentionally absent.
#'
#' @param scientific_plan_id Stable identity of relation, frame, pairing, and
#'   query/materialization semantics.
#' @param domain_signature Optional exact neural-domain signature. Compiler
#'   receipts always provide it; standalone receipts may omit it.
#' @param compute An `effect_compute_policy`.
#' @param sources A nonempty list of `effect_source_capabilities` values.
#' @param memory An `effect_memory_plan`.
#' @param kernel_version,task_partition_id,reduction_plan_id,precision Nonempty
#'   execution identity strings.
#' @param numeric_contract An `effect_numerical_contract`.
#' @param completion_status One of `complete`, `planned`, `failed`, or
#'   `interrupted`.
#' @param task_count,completed_task_count Nonnegative exact whole task counts.
#' @param elapsed_seconds Nonnegative finite elapsed wall time.
#' @param blas A list naming the BLAS vendor and positive thread count.
#' @return An immutable-by-convention execution receipt.
#' @export
execution_receipt <- function(scientific_plan_id, compute, sources, memory,
                              kernel_version, task_partition_id,
                              reduction_plan_id, precision = "double",
                              numeric_contract = numerical_contract(),
                              completion_status = "complete",
                              task_count = 1L,
                              completed_task_count = NULL,
                              elapsed_seconds = 0,
                              blas = list(vendor = "unknown", threads = 1L),
                              domain_signature = NULL) {
  ids <- list(
    scientific_plan_id = scientific_plan_id,
    kernel_version = kernel_version,
    task_partition_id = task_partition_id,
    reduction_plan_id = reduction_plan_id,
    precision = precision
  )
  if (!all(vapply(ids, function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }, logical(1)))) {
    stop("Receipt identities must be nonempty character scalars.", call. = FALSE)
  }
  if (!identical(precision, "double")) {
    stop("effectagram 0.1 supports only `precision = \"double\"`.", call. = FALSE)
  }
  if (!is.null(domain_signature) &&
      (!is.character(domain_signature) || length(domain_signature) != 1L ||
       is.na(domain_signature) ||
       !grepl("^sha256:[[:xdigit:]]{64}$", domain_signature))) {
    stop("`domain_signature` must be NULL or one strong sha256 identifier.",
      call. = FALSE)
  }
  completion_status <- match.arg(completion_status,
    c("complete", "planned", "failed", "interrupted"))
  if (is.null(completed_task_count)) {
    completed_task_count <- if (completion_status == "complete") task_count else 0L
  }

  compute <- .validate_compute_policy(compute)
  if (!is.list(sources) || length(sources) < 1L ||
      !all(vapply(sources, inherits, logical(1), "effect_source_capabilities"))) {
    stop("`sources` must be a nonempty list of source capabilities.", call. = FALSE)
  }
  sources <- lapply(sources, .validate_source_capabilities)
  if (!all(vapply(sources, function(x) isTRUE(x$block_read), logical(1)))) {
    stop("Receipt sources must support bounded block reads.", call. = FALSE)
  }
  memory <- .validate_memory_plan_for_receipt(memory)
  numeric_contract <- .validate_numeric_contract_for_receipt(numeric_contract)
  blas <- .validate_blas_identity(blas)
  task_count <- .receipt_whole(task_count, "task_count")
  completed_task_count <- .receipt_whole(completed_task_count,
    "completed_task_count")
  if (completed_task_count > task_count) {
    stop("Completed task count cannot exceed total task count.", call. = FALSE)
  }
  if (completion_status == "complete" && completed_task_count != task_count) {
    stop("A complete receipt must report every task completed.", call. = FALSE)
  }
  if (completion_status == "planned" && completed_task_count != 0) {
    stop("A planned receipt cannot report completed tasks.", call. = FALSE)
  }
  if (!is.numeric(elapsed_seconds) || length(elapsed_seconds) != 1L ||
      is.na(elapsed_seconds) || !is.finite(elapsed_seconds) || elapsed_seconds < 0) {
    stop("`elapsed_seconds` must be one nonnegative finite number.", call. = FALSE)
  }
  if (!identical(memory$workers, compute$workers)) {
    stop("Receipt compute and memory plans must use the same worker count.",
      call. = FALSE)
  }
  if (!identical(is.null(compute$memory_bytes), is.null(memory$budget_bytes)) ||
      (!is.null(compute$memory_bytes) &&
       !isTRUE(all.equal(compute$memory_bytes, memory$budget_bytes,
         tolerance = 0)))) {
    stop("Receipt compute and memory plans must declare the same memory budget.",
      call. = FALSE)
  }

  receipt <- structure(
    c(ids["scientific_plan_id"], list(domain_signature = domain_signature),
      ids[-1L], list(
      compute = compute,
      sources = sources,
      memory = memory,
      numeric_contract = numeric_contract,
      completion_status = completion_status,
      task_count = task_count,
      completed_task_count = completed_task_count,
      elapsed_seconds = elapsed_seconds,
      blas = blas
    )),
    class = "effect_execution_receipt"
  )
  .validate_execution_receipt(receipt)
}

.receipt_whole <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0 || value %% 1 != 0 || value > 2^53) {
    stop(sprintf("`%s` must be one nonnegative exact whole number.", name),
      call. = FALSE)
  }
  value
}

.validate_memory_plan_for_receipt <- function(memory) {
  if (!inherits(memory, "effect_memory_plan") || !is.numeric(memory$categories)) {
    stop("`memory` must be an effectagram memory plan.", call. = FALSE)
  }
  expected_categories <- c(
    "shared_source", "runtime_reserve", "private_source_per_worker", "source_block",
    "relation_block", "atom_block", "output", "contraction",
    "replacement_copy", "serialization_overlap", "reorder_buffer",
    "checkpoint_buffer"
  )
  if (!identical(names(memory$categories), expected_categories)) {
    stop("Memory-plan categories are missing or noncanonical.", call. = FALSE)
  }
  c <- memory$categories
  rebuilt <- memory_plan(
    shared_source_bytes = c[["shared_source"]],
    runtime_reserve_bytes = c[["runtime_reserve"]],
    private_source_bytes = c[["private_source_per_worker"]],
    source_block_bytes = c[["source_block"]],
    relation_block_bytes = c[["relation_block"]],
    atom_block_bytes = c[["atom_block"]],
    output_bytes = c[["output"]],
    contraction_bytes = c[["contraction"]],
    replacement_copy_bytes = c[["replacement_copy"]],
    serialization_overlap_bytes = c[["serialization_overlap"]],
    reorder_buffer_bytes = c[["reorder_buffer"]],
    checkpoint_buffer_bytes = c[["checkpoint_buffer"]],
    workers = memory$workers,
    n_active = memory$n_active,
    safety_factor = memory$safety_factor,
    budget_bytes = memory$budget_bytes,
    measured_peak_bytes = memory$measured_peak_bytes
  )
  derived <- c(
    "worker_private_total_bytes", "task_private_per_active_bytes",
    "active_task_total_bytes", "modeled_peak_bytes", "conservative_peak_bytes",
    "fits_budget", "measurement_within_plan", "prediction_kind"
  )
  if (!all(vapply(derived, function(name) identical(memory[[name]], rebuilt[[name]]),
    logical(1)))) {
    stop("Memory-plan derived fields are inconsistent with its categories.",
      call. = FALSE)
  }
  rebuilt
}

.validate_numeric_contract_for_receipt <- function(contract) {
  if (!inherits(contract, "effect_numerical_contract")) {
    stop("`numeric_contract` must be an effectagram numerical contract.",
      call. = FALSE)
  }
  rebuilt <- numerical_contract(contract$atol, contract$rtol)
  if (!identical(contract, rebuilt)) {
    stop("Numerical contract fields are inconsistent or noncanonical.",
      call. = FALSE)
  }
  rebuilt
}

.validate_blas_identity <- function(blas) {
  if (!is.list(blas) || !identical(names(blas), c("vendor", "threads")) ||
      !is.character(blas$vendor) || length(blas$vendor) != 1L ||
      is.na(blas$vendor) || !nzchar(blas$vendor) ||
      !is.numeric(blas$threads) || length(blas$threads) != 1L ||
      is.na(blas$threads) || !is.finite(blas$threads) ||
      blas$threads < 1 || blas$threads %% 1 != 0) {
    stop("`blas` must name one vendor and one positive whole thread count.",
      call. = FALSE)
  }
  if (blas$threads != 1) {
    stop("effectagram 0.1 receipts require a one-thread BLAS identity.",
      call. = FALSE)
  }
  list(vendor = blas$vendor, threads = 1L)
}

.validate_execution_receipt <- function(receipt) {
  if (!inherits(receipt, "effect_execution_receipt")) {
    stop("`receipt` must be an effectagram execution receipt.", call. = FALSE)
  }
  expected <- c(
    "scientific_plan_id", "domain_signature", "kernel_version", "task_partition_id",
    "reduction_plan_id", "precision", "compute", "sources", "memory",
    "numeric_contract", "completion_status", "task_count",
    "completed_task_count", "elapsed_seconds", "blas"
  )
  if (!identical(names(receipt), expected)) {
    stop("Execution receipt fields are missing or noncanonical.", call. = FALSE)
  }
  ids <- receipt[c("scientific_plan_id", "kernel_version", "task_partition_id",
    "reduction_plan_id", "precision")]
  if (!all(vapply(ids, function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }, logical(1))) || !identical(receipt$precision, "double")) {
    stop("Execution receipt identities or precision are invalid.", call. = FALSE)
  }
  if (!is.null(receipt$domain_signature) &&
      (!is.character(receipt$domain_signature) ||
       length(receipt$domain_signature) != 1L || is.na(receipt$domain_signature) ||
       !grepl("^sha256:[[:xdigit:]]{64}$", receipt$domain_signature))) {
    stop("Execution receipt domain signature is invalid.", call. = FALSE)
  }
  canonical_compute <- .validate_compute_policy(receipt$compute)
  if (!is.list(receipt$sources) || length(receipt$sources) < 1L) {
    stop("Execution receipt sources must be a nonempty list.", call. = FALSE)
  }
  canonical_sources <- lapply(receipt$sources, .validate_source_capabilities)
  if (!all(vapply(canonical_sources, function(x) isTRUE(x$block_read), logical(1)))) {
    stop("Receipt sources must support bounded block reads.", call. = FALSE)
  }
  canonical_memory <- .validate_memory_plan_for_receipt(receipt$memory)
  .validate_numeric_contract_for_receipt(receipt$numeric_contract)
  .validate_blas_identity(receipt$blas)
  status <- receipt$completion_status
  if (!is.character(status) || length(status) != 1L ||
      !status %in% c("complete", "planned", "failed", "interrupted")) {
    stop("Execution receipt completion status is invalid.", call. = FALSE)
  }
  task_count <- .receipt_whole(receipt$task_count, "task_count")
  completed <- .receipt_whole(receipt$completed_task_count, "completed_task_count")
  if (completed > task_count || (status == "complete" && completed != task_count) ||
      (status == "planned" && completed != 0)) {
    stop("Execution receipt completion counts contradict its status.", call. = FALSE)
  }
  if (!is.numeric(receipt$elapsed_seconds) || length(receipt$elapsed_seconds) != 1L ||
      is.na(receipt$elapsed_seconds) || !is.finite(receipt$elapsed_seconds) ||
      receipt$elapsed_seconds < 0) {
    stop("Execution receipt elapsed time is invalid.", call. = FALSE)
  }
  if (!identical(canonical_memory$workers, canonical_compute$workers) ||
      !identical(is.null(canonical_compute$memory_bytes),
        is.null(canonical_memory$budget_bytes)) ||
      (!is.null(canonical_compute$memory_bytes) &&
       !isTRUE(all.equal(canonical_compute$memory_bytes,
         canonical_memory$budget_bytes, tolerance = 0)))) {
    stop("Execution receipt compute and memory policies are inconsistent.",
      call. = FALSE)
  }
  invisible(receipt)
}
