# Source and execution provenance -------------------------------------------

#' Declare source execution capabilities
#'
#' `source_capabilities()` states what an out-of-memory relation source can
#' actually do, so the compiler can choose a bounded execution route and
#' record the source revision in the receipt. Declare it when supplying a
#' custom source such as [file_matrix_source()].
#'
#' @param block_read Whether bounded feature-block reads are supported.
#' @param reopenable Whether a fresh read-only handle can be opened safely.
#' @param thread_safe Whether concurrent reads through one handle are supported.
#' @param stable_revision A strong immutable source revision or checksum.
#' @return An `effect_source_capabilities` value with the three logical flags
#'   `$block_read`, `$reopenable`, `$thread_safe`, and the
#'   `$stable_revision` identifier bound into the execution receipt.
#' @seealso [file_matrix_source()], which carries these capabilities, and
#'   [compute_policy()], which is checked against them before execution.
#' @family relation planning and fitting
#' @examples
#' # A block-readable, reopenable source identified by a content checksum.
#' revision <- paste0("sha256:", paste(rep("a", 64), collapse = ""))
#' capabilities <- source_capabilities(
#'   block_read = TRUE, reopenable = TRUE, stable_revision = revision
#' )
#' capabilities$block_read
#' capabilities$thread_safe
#'
#' # The revision must be a strong identifier: a mutable label such as a file
#' # modification time is refused, because receipts must stay verifiable.
#' weak <- try(
#'   source_capabilities(TRUE, stable_revision = "2026-08-15"), silent = TRUE
#' )
#' conditionMessage(attr(weak, "condition"))
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
#' @param observed Canonical observed execution facts. Planned receipts use an
#'   empty observation record that is populated by the coordinator.
#' @return An immutable-by-convention execution receipt.
#' @keywords internal
execution_receipt <- function(scientific_plan_id, compute, sources, memory,
                              kernel_version, task_partition_id,
                              reduction_plan_id, precision = "double",
                              numeric_contract = numerical_contract(),
                              completion_status = "complete",
                              task_count = 1L,
                              completed_task_count = NULL,
                              elapsed_seconds = 0,
                              blas = list(vendor = "unknown", requested_threads = 1L,
                                observed_threads = NA_integer_),
                              domain_signature = NULL,
                              observed = .empty_execution_observations()) {
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
    stop("crossform 0.1 supports only `precision = \"double\"`.", call. = FALSE)
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
  observed <- .validate_execution_observations(observed)
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
  if (!identical(is.null(compute$workspace_bytes), is.null(memory$budget_bytes)) ||
      (!is.null(compute$workspace_bytes) &&
       !isTRUE(all.equal(compute$workspace_bytes, memory$budget_bytes,
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
      blas = blas,
      observed = observed
    )),
    class = "effect_execution_receipt"
  )
  .validate_execution_receipt(receipt)
}

# Derive the receipt for a view projected from a materialized geometry. The
# projection inherits the parent execution's observed facts and marks itself
# in the task-partition identity; only the scientific id changes, to the
# route-stable view identity shared with fused query-first execution.
.projection_receipt <- function(receipt, scientific_plan_id) {
  .validate_execution_receipt(receipt)
  execution_receipt(
    scientific_plan_id = scientific_plan_id,
    compute = receipt$compute,
    sources = receipt$sources,
    memory = receipt$memory,
    kernel_version = receipt$kernel_version,
    task_partition_id = paste0(
      receipt$task_partition_id, "+projected_view"
    ),
    reduction_plan_id = receipt$reduction_plan_id,
    precision = receipt$precision,
    numeric_contract = receipt$numeric_contract,
    completion_status = receipt$completion_status,
    task_count = receipt$task_count,
    completed_task_count = receipt$completed_task_count,
    elapsed_seconds = receipt$elapsed_seconds,
    blas = receipt$blas,
    domain_signature = receipt$domain_signature,
    observed = receipt$observed
  )
}

.empty_execution_observations <- function() {
  list(
    task_counts = c(planned = 0, started = 0, completed = 0, failed = 0,
      retried = 0),
    features_completed = 0,
    bytes_read = 0,
    tiles = list(feature_block = NULL, row_tile = NULL, coordinate_tile = NULL),
    stage_seconds = numeric(),
    source_access = character(),
    memory = list(baseline_rss_bytes = NULL, incremental_peak_rss_bytes = NULL,
      absolute_peak_rss_bytes = NULL, measured_workspace_bytes = NULL),
    runtime = list(package_version = as.character(utils::packageVersion("crossform")),
      r_version = R.version.string, platform = R.version$platform),
    reporter_failures = character(),
    checkpoint = list(enabled = FALSE, resumed = FALSE, artifacts_written = 0L),
    cleanup = list(attempted = FALSE, success = NA, message = NULL)
  )
}

.validate_execution_observations <- function(x) {
  expected <- names(.empty_execution_observations())
  if (!is.list(x) || !identical(names(x), expected)) {
    stop("Execution observations are missing or noncanonical.", call. = FALSE)
  }
  expected_counts <- c("planned", "started", "completed", "failed", "retried")
  if (!is.numeric(x$task_counts) || !identical(names(x$task_counts), expected_counts) ||
      anyNA(x$task_counts) || any(!is.finite(x$task_counts)) ||
      any(x$task_counts < 0) || any(x$task_counts %% 1 != 0)) {
    stop("Observed task counts are invalid.", call. = FALSE)
  }
  .receipt_whole(x$features_completed, "features_completed")
  if (!is.numeric(x$bytes_read) || length(x$bytes_read) != 1L ||
      is.na(x$bytes_read) || !is.finite(x$bytes_read) || x$bytes_read < 0) {
    stop("Observed bytes read are invalid.", call. = FALSE)
  }
  tile_names <- c("feature_block", "row_tile", "coordinate_tile")
  if (!is.list(x$tiles) || !identical(names(x$tiles), tile_names) ||
      any(!vapply(x$tiles, function(value) is.null(value) ||
        (is.numeric(value) && length(value) == 1L && !is.na(value) &&
         is.finite(value) && value >= 1 && value %% 1 == 0), logical(1)))) {
    stop("Observed execution tiles are invalid.", call. = FALSE)
  }
  if (!is.numeric(x$stage_seconds) || anyNA(x$stage_seconds) ||
      any(!is.finite(x$stage_seconds)) || any(x$stage_seconds < 0) ||
      (length(x$stage_seconds) &&
       (is.null(names(x$stage_seconds)) || any(!nzchar(names(x$stage_seconds)))))) {
    stop("Observed stage timings are invalid.", call. = FALSE)
  }
  if (!is.character(x$source_access) || anyNA(x$source_access) ||
      (length(x$source_access) && (is.null(names(x$source_access)) ||
       any(!nzchar(names(x$source_access)))))) {
    stop("Observed source access modes are invalid.", call. = FALSE)
  }
  memory_names <- c("baseline_rss_bytes", "incremental_peak_rss_bytes",
    "absolute_peak_rss_bytes", "measured_workspace_bytes")
  if (!is.list(x$memory) || !identical(names(x$memory), memory_names) ||
      any(!vapply(x$memory, function(value) is.null(value) ||
        (is.numeric(value) && length(value) == 1L && !is.na(value) &&
         is.finite(value) && value >= 0), logical(1)))) {
    stop("Observed memory evidence is invalid.", call. = FALSE)
  }
  if (!is.list(x$runtime) ||
      !identical(names(x$runtime), c("package_version", "r_version", "platform")) ||
      any(!vapply(x$runtime, function(value) is.character(value) &&
        length(value) == 1L && !is.na(value) && nzchar(value), logical(1)))) {
    stop("Observed runtime identity is invalid.", call. = FALSE)
  }
  if (!is.character(x$reporter_failures) || anyNA(x$reporter_failures)) {
    stop("Observed reporter failures are invalid.", call. = FALSE)
  }
  if (!is.list(x$checkpoint) ||
      !identical(names(x$checkpoint), c("enabled", "resumed", "artifacts_written")) ||
      !all(vapply(x$checkpoint[c("enabled", "resumed")], function(value) {
        is.logical(value) && length(value) == 1L && !is.na(value)
      }, logical(1)))) {
    stop("Observed checkpoint state is invalid.", call. = FALSE)
  }
  .receipt_whole(x$checkpoint$artifacts_written, "artifacts_written")
  if (!is.list(x$cleanup) ||
      !identical(names(x$cleanup), c("attempted", "success", "message")) ||
      !is.logical(x$cleanup$attempted) || length(x$cleanup$attempted) != 1L ||
      is.na(x$cleanup$attempted) || !is.logical(x$cleanup$success) ||
      length(x$cleanup$success) != 1L ||
      !(is.na(x$cleanup$success) || !is.na(x$cleanup$success)) ||
      !(is.null(x$cleanup$message) || (is.character(x$cleanup$message) &&
        length(x$cleanup$message) == 1L && !is.na(x$cleanup$message)))) {
    stop("Observed cleanup outcome is invalid.", call. = FALSE)
  }
  x
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
    stop("`memory` must be a crossform memory plan.", call. = FALSE)
  }
  expected_categories <- c(
    "frame", "resident_source", "source_handles", "source_block",
    "relation_block", "atom_block", "local_state", "output", "contraction",
    "replacement_copy", "serialization_overlap", "reorder_buffer",
    "checkpoint_buffer"
  )
  if (!identical(names(memory$categories), expected_categories)) {
    stop("Memory-plan categories are missing or noncanonical.", call. = FALSE)
  }
  c <- memory$categories
  rebuilt <- memory_plan(
    frame_bytes = c[["frame"]],
    resident_source_bytes = c[["resident_source"]],
    source_handle_bytes = c[["source_handles"]],
    source_block_bytes = c[["source_block"]],
    relation_block_bytes = c[["relation_block"]],
    atom_block_bytes = c[["atom_block"]],
    local_state_bytes = c[["local_state"]],
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
    measured_workspace_bytes = memory$measured_workspace_bytes,
    baseline_rss_bytes = memory$baseline_rss_bytes,
    peak_rss_bytes = memory$absolute_peak_rss_bytes
  )
  derived <- c(
    "persistent_workspace_bytes", "task_workspace_per_active_bytes",
    "active_task_workspace_bytes", "modeled_workspace_bytes",
    "planned_workspace_bytes", "fits_budget", "measured_workspace_within_plan",
    "incremental_peak_rss_bytes", "prediction_kind"
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
    stop("`numeric_contract` must be a crossform numerical contract.",
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
  if (!is.list(blas) || !identical(names(blas),
      c("vendor", "requested_threads", "observed_threads")) ||
      !is.character(blas$vendor) || length(blas$vendor) != 1L ||
      is.na(blas$vendor) || !nzchar(blas$vendor) ||
      !is.numeric(blas$requested_threads) || length(blas$requested_threads) != 1L ||
      is.na(blas$requested_threads) || !is.finite(blas$requested_threads) ||
      blas$requested_threads != 1 ||
      !is.numeric(blas$observed_threads) || length(blas$observed_threads) != 1L ||
      (!is.na(blas$observed_threads) && (!is.finite(blas$observed_threads) ||
       blas$observed_threads < 1 || blas$observed_threads %% 1 != 0))) {
    stop("`blas` must separate requested and observed thread counts.",
      call. = FALSE)
  }
  list(vendor = blas$vendor, requested_threads = 1L,
    observed_threads = if (is.na(blas$observed_threads)) NA_integer_ else
      as.integer(blas$observed_threads))
}

.validate_execution_receipt <- function(receipt) {
  if (!inherits(receipt, "effect_execution_receipt")) {
    stop("`receipt` must be a crossform execution receipt.", call. = FALSE)
  }
  expected <- c(
    "scientific_plan_id", "domain_signature", "kernel_version", "task_partition_id",
    "reduction_plan_id", "precision", "compute", "sources", "memory",
    "numeric_contract", "completion_status", "task_count",
    "completed_task_count", "elapsed_seconds", "blas", "observed"
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
  .validate_execution_observations(receipt$observed)
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
      !identical(is.null(canonical_compute$workspace_bytes),
        is.null(canonical_memory$budget_bytes)) ||
      (!is.null(canonical_compute$workspace_bytes) &&
       !isTRUE(all.equal(canonical_compute$workspace_bytes,
         canonical_memory$budget_bytes, tolerance = 0)))) {
    stop("Execution receipt compute and memory policies are inconsistent.",
      call. = FALSE)
  }
  invisible(receipt)
}
