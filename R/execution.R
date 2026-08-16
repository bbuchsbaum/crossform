# Sequential execution policy -----------------------------------------------

#' Construct the version 0.1 compute policy
#'
#' Version 0.1 deliberately owns no process pool and accepts exactly one
#' worker. Participant-level parallelism belongs outside the core geometry
#' call until a later executor passes its memory and determinism gates.
#'
#' @param workers Number of R workers. Must be exactly one in version 0.1.
#' @param block_features Optional positive feature-block size.
#' @param workspace_bytes Optional positive budget for crossform-owned live
#'   workspace. Baseline and total process RSS are not charged to this budget.
#' @return An `effect_compute_policy` recording `$workers`,
#'   `$block_features`, `$workspace_bytes`, and the fixed
#'   `$process_backend`. It is a declaration consumed by [plan_geometry()]
#'   and friends; it never starts a worker itself.
#' @seealso [plan_geometry()] and [measurement_form()], which accept this
#'   policy; [catch_refusal()] for inspecting the refusal below.
#'
#' @section Refusal:
#' Requesting more than one worker signals an `effect_capability_refusal`
#' carrying capability `"parallel_execution"` in namespace `"compute_policy"`,
#' with reason `"worker_pool_not_implemented"`. Branch on it with
#' [catch_refusal()] rather than on the message text.
#' @family geometry plans and views
#' @examples
#' # The default: one worker, no declared block or workspace ceiling.
#' default <- compute_policy()
#' c(workers = default$workers, backend = default$process_backend)
#'
#' # Bound the feature block and the crossform-owned workspace. Neither
#' # choice changes the estimand, only the execution receipt.
#' policy <- compute_policy(block_features = 64L, workspace_bytes = 64 * 1024^2)
#' policy$block_features
#'
#' # Version 0.1 owns no process pool, so more than one worker is refused —
#' # as a classed capability refusal, not a message to match on.
#' refusal <- catch_refusal(compute_policy(workers = 4L))
#' refusal$capability
#' refusal$reasons
#' refusal$remedies
#' @export
compute_policy <- function(workers = 1L, block_features = NULL,
                           workspace_bytes = NULL) {
  policy <- structure(
    list(
      workers = workers,
      block_features = block_features,
      workspace_bytes = workspace_bytes,
      process_backend = "sequential"
    ),
    class = "effect_compute_policy"
  )
  .validate_compute_policy(policy)
}

.validate_compute_policy <- function(policy) {
  if (!inherits(policy, "effect_compute_policy")) {
    stop("`compute` must be a crossform compute policy.", call. = FALSE)
  }
  expected_names <- c("workers", "block_features", "workspace_bytes",
    "process_backend")
  if (!identical(names(policy), expected_names)) {
    stop("Compute policy fields are missing, extra, or out of canonical order.",
      call. = FALSE)
  }
  workers <- policy$workers
  if (!is.numeric(workers) || length(workers) != 1L || is.na(workers) ||
      !is.finite(workers) || workers != 1) {
    .capability_refusal(sprintf(paste0(
      "crossform 0.1 owns no process pool, so `workers` must be 1; received ",
      "%s. Sequential execution is a capability boundary, not a performance ",
      "default: an executor that spawns workers has to pass memory and ",
      "determinism gates first."
    ), if (is.numeric(workers) && length(workers) == 1L && !is.na(workers)) {
      paste0("`", format(workers), "`")
    } else {
      .msg_value(workers)
    }),
      capability = "parallel_execution",
      namespace = "compute_policy",
      reasons = "worker_pool_not_implemented",
      remedies = paste0(
        "Pass `workers = 1` and parallelize across participants outside the ",
        "geometry call, or bound the work with ",
        "`compute_policy(block_features = , workspace_bytes = )`."
      )
    )
  }
  if (!is.null(policy$block_features) &&
      (!is.numeric(policy$block_features) || length(policy$block_features) != 1L ||
       is.na(policy$block_features) || !is.finite(policy$block_features) ||
       policy$block_features < 1 || policy$block_features %% 1 != 0)) {
    stop("`block_features` must be NULL or one positive integer.", call. = FALSE)
  }
  if (!is.null(policy$workspace_bytes) &&
      (!is.numeric(policy$workspace_bytes) ||
       length(policy$workspace_bytes) != 1L ||
       is.na(policy$workspace_bytes) || !is.finite(policy$workspace_bytes) ||
       policy$workspace_bytes <= 0)) {
    stop("`workspace_bytes` must be NULL or one positive finite number.",
      call. = FALSE)
  }
  if (!identical(policy$process_backend, "sequential")) {
    stop("crossform 0.1 supports only the sequential process backend.",
      call. = FALSE)
  }
  structure(
    list(
      workers = 1L,
      block_features = if (is.null(policy$block_features)) NULL else
        as.integer(policy$block_features),
      workspace_bytes = policy$workspace_bytes,
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
