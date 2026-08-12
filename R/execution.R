# Sequential execution policy -----------------------------------------------

#' Construct the version 0.1 compute policy
#'
#' Version 0.1 deliberately owns no process pool and accepts exactly one
#' worker. Participant-level parallelism belongs outside the core geometry
#' call until a later executor passes its memory and determinism gates.
#'
#' @param workers Number of R workers. Must be exactly one in version 0.1.
#' @param block_features Optional positive feature-block size.
#' @param memory_bytes Optional positive conservative memory budget in bytes.
#' @return An immutable-by-convention declarative compute policy.
#' @export
compute_policy <- function(workers = 1L, block_features = NULL,
                           memory_bytes = NULL) {
  policy <- structure(
    list(
      workers = workers,
      block_features = block_features,
      memory_bytes = memory_bytes,
      process_backend = "sequential"
    ),
    class = "effect_compute_policy"
  )
  .validate_compute_policy(policy)
}

.validate_compute_policy <- function(policy) {
  if (!inherits(policy, "effect_compute_policy")) {
    stop("`compute` must be an effectagram compute policy.", call. = FALSE)
  }
  expected_names <- c("workers", "block_features", "memory_bytes", "process_backend")
  if (!identical(names(policy), expected_names)) {
    stop("Compute policy fields are missing, extra, or out of canonical order.",
      call. = FALSE)
  }
  workers <- policy$workers
  if (!is.numeric(workers) || length(workers) != 1L || is.na(workers) ||
      !is.finite(workers) || workers != 1) {
    stop("workers > 1 is not implemented in effectagram 0.1; `workers` must be 1.",
      call. = FALSE)
  }
  if (!is.null(policy$block_features) &&
      (!is.numeric(policy$block_features) || length(policy$block_features) != 1L ||
       is.na(policy$block_features) || !is.finite(policy$block_features) ||
       policy$block_features < 1 || policy$block_features %% 1 != 0)) {
    stop("`block_features` must be NULL or one positive integer.", call. = FALSE)
  }
  if (!is.null(policy$memory_bytes) &&
      (!is.numeric(policy$memory_bytes) || length(policy$memory_bytes) != 1L ||
       is.na(policy$memory_bytes) || !is.finite(policy$memory_bytes) ||
       policy$memory_bytes <= 0)) {
    stop("`memory_bytes` must be NULL or one positive finite number.", call. = FALSE)
  }
  if (!identical(policy$process_backend, "sequential")) {
    stop("effectagram 0.1 supports only the sequential process backend.",
      call. = FALSE)
  }
  structure(
    list(
      workers = 1L,
      block_features = if (is.null(policy$block_features)) NULL else
        as.integer(policy$block_features),
      memory_bytes = policy$memory_bytes,
      process_backend = "sequential"
    ),
    class = "effect_compute_policy"
  )
}

# This internal boundary encodes a critical ordering contract: compute-policy
# validation happens before even source-capability inspection. The eventual
# geometry compiler must call this function before opening or reading a source.
.execution_preflight <- function(compute, inspect_source) {
  compute <- .validate_compute_policy(compute)
  if (!is.function(inspect_source)) {
    stop("`inspect_source` must be a function.", call. = FALSE)
  }
  capabilities <- inspect_source()
  structure(
    list(
      compute = compute,
      workers = 1L,
      process_backend = "sequential",
      source_capabilities = capabilities
    ),
    class = "effect_execution_plan"
  )
}
