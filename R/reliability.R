# Failure receipts and observer isolation -----------------------------------

.safe_report <- function(reporter, event) {
  if (is.null(reporter)) return(NULL)
  if (!is.function(reporter)) {
    return("reporter is not a function")
  }
  tryCatch(
    {
      reporter(event)
      NULL
    },
    error = function(error) conditionMessage(error),
    interrupt = function(interrupt) conditionMessage(interrupt)
  )
}

.run_cleanup <- function(cleanup) {
  if (is.null(cleanup)) {
    return(list(attempted = FALSE, success = TRUE, message = NULL))
  }
  if (!is.function(cleanup)) {
    return(list(
      attempted = TRUE,
      success = FALSE,
      message = "cleanup is not a function"
    ))
  }
  tryCatch(
    {
      cleanup()
      list(attempted = TRUE, success = TRUE, message = NULL)
    },
    error = function(error) list(
      attempted = TRUE,
      success = FALSE,
      message = conditionMessage(error)
    ),
    interrupt = function(interrupt) list(
      attempted = TRUE,
      success = FALSE,
      message = paste0("cleanup interrupted: ", conditionMessage(interrupt))
    )
  )
}

.execution_error <- function(message, receipt, cleanup_status,
                             cause = NULL, observer_failures = character()) {
  structure(
    list(
      message = message,
      call = NULL,
      receipt = receipt,
      cleanup_status = cleanup_status,
      cause = cause,
      observer_failures = observer_failures
    ),
    class = c("effect_execution_error", "error", "condition")
  )
}

.execution_interrupt <- function(message, receipt, cleanup_status,
                                 cause = NULL, observer_failures = character()) {
  structure(
    list(
      message = message,
      call = NULL,
      receipt = receipt,
      cleanup_status = cleanup_status,
      cause = cause,
      observer_failures = observer_failures
    ),
    class = c("effect_execution_interrupt", "interrupt", "condition")
  )
}

.receipt_with_status <- function(receipt, status, elapsed_seconds) {
  receipt$completion_status <- status
  if (status == "complete") receipt$completed_task_count <- receipt$task_count
  receipt$elapsed_seconds <- elapsed_seconds
  .validate_execution_receipt(receipt)
  receipt
}

.execute_guarded <- function(compute, receipt, reporter = NULL, cleanup = NULL) {
  if (!is.function(compute)) {
    stop("`compute` must be a function.", call. = FALSE)
  }
  .validate_execution_receipt(receipt)
  state <- new.env(parent = emptyenv())
  state$finalized <- FALSE
  state$cleanup_status <- NULL
  state$receipt <- receipt
  state$observer_failures <- character()
  started <- proc.time()[["elapsed"]]
  report_failure <- .safe_report(reporter, list(type = "start", receipt = receipt))
  if (!is.null(report_failure)) {
    state$observer_failures <- c(state$observer_failures, report_failure)
  }

  finalize <- function(requested_status) {
    if (state$finalized) return(invisible(NULL))
    state$cleanup_status <- .run_cleanup(cleanup)
    final_status <- if (requested_status == "complete" &&
      !state$cleanup_status$success) "failed" else requested_status
    state$receipt <- .receipt_with_status(
      receipt, final_status, max(0, proc.time()[["elapsed"]] - started)
    )
    state$finalized <- TRUE
    report_failure <- .safe_report(reporter, list(
      type = final_status,
      receipt = state$receipt,
      cleanup_status = state$cleanup_status
    ))
    if (!is.null(report_failure)) {
      state$observer_failures <- c(state$observer_failures, report_failure)
    }
    invisible(NULL)
  }
  on.exit({
    if (!state$finalized) finalize("failed")
  }, add = TRUE)

  outcome <- tryCatch(
    list(kind = "success", value = compute()),
    interrupt = function(interrupt) list(kind = "interrupt", cause = interrupt),
    error = function(error) list(kind = "error", cause = error)
  )

  if (outcome$kind == "interrupt") {
    finalize("interrupted")
    stop(.execution_interrupt(
      conditionMessage(outcome$cause), state$receipt, state$cleanup_status,
      cause = outcome$cause, observer_failures = state$observer_failures
    ))
  }
  if (outcome$kind == "error") {
    finalize("failed")
    stop(.execution_error(
      conditionMessage(outcome$cause), state$receipt, state$cleanup_status,
      cause = outcome$cause, observer_failures = state$observer_failures
    ))
  }
  finalize("complete")
  if (!state$cleanup_status$success) {
    stop(.execution_error(
      paste0("Execution completed but cleanup failed: ",
        state$cleanup_status$message),
      state$receipt, state$cleanup_status,
      observer_failures = state$observer_failures
    ))
  }
  outcome$value
}
