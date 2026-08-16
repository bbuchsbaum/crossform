test_that("reporter failure cannot alter a successful scientific result", {
  receipt <- receipt_fixture()
  oracle <- matrix(c(1, 2, 3, 4), 2)
  got <- crossform:::.execute_guarded(
    compute = function() oracle,
    receipt = receipt,
    reporter = function(event) stop("display failed"),
    cleanup = function() invisible(TRUE)
  )

  expect_identical(got, oracle)
})

test_that("execution failure carries receipt and successful cleanup status", {
  receipt <- receipt_fixture()
  cleaned <- FALSE
  condition <- tryCatch(
    crossform:::.execute_guarded(
      compute = function() stop("kernel failed"),
      receipt = receipt,
      cleanup = function() cleaned <<- TRUE
    ),
    effect_execution_error = identity
  )

  expect_s3_class(condition, "effect_execution_error")
  expect_identical(condition$receipt$scientific_plan_id, receipt$scientific_plan_id)
  expect_true(condition$cleanup_status$success)
  expect_true(cleaned)
  expect_match(condition$message, "kernel failed")
  expect_identical(condition$receipt$completion_status, "failed")
  expect_true(condition$receipt$observed$cleanup$success)
})

test_that("cleanup failure is itself a receipt-bearing execution failure", {
  receipt <- receipt_fixture()
  condition <- tryCatch(
    crossform:::.execute_guarded(
      compute = function() 42,
      receipt = receipt,
      cleanup = function() stop("close failed")
    ),
    effect_execution_error = identity
  )

  expect_s3_class(condition, "effect_execution_error")
  expect_false(condition$cleanup_status$success)
  expect_match(condition$cleanup_status$message, "close failed")
  expect_identical(condition$receipt$scientific_plan_id, receipt$scientific_plan_id)
  expect_identical(condition$receipt$completion_status, "failed")
  expect_false(condition$receipt$observed$cleanup$success)
})

test_that("reporter diagnostics remain nonsemantic on computational failure", {
  receipt <- receipt_fixture()
  condition <- tryCatch(
    crossform:::.execute_guarded(
      compute = function() stop("science failed"),
      receipt = receipt,
      reporter = function(event) stop("observer failed")
    ),
    effect_execution_error = identity
  )

  expect_match(condition$message, "science failed")
  expect_true(length(condition$observer_failures) >= 1)
  expect_false("reporter" %in% names(condition$receipt))
  expect_true(length(condition$receipt$observed$reporter_failures) >= 1)
})

test_that("synthetic interrupts preserve interrupt semantics and clean exactly once", {
  cleaned <- 0L
  events <- list()
  interrupt <- structure(
    list(message = "user interrupt", call = NULL),
    class = c("interrupt", "condition")
  )
  condition <- tryCatch(
    crossform:::.execute_guarded(
      compute = function() signalCondition(interrupt),
      receipt = receipt_fixture(),
      reporter = function(event) events[[length(events) + 1L]] <<- event,
      cleanup = function() cleaned <<- cleaned + 1L
    ),
    interrupt = identity
  )

  expect_s3_class(condition, "effect_execution_interrupt")
  expect_s3_class(condition, "interrupt")
  expect_identical(cleaned, 1L)
  expect_true(condition$cleanup_status$success)
  expect_identical(condition$receipt$completion_status, "interrupted")
  expect_identical(events[[length(events)]]$type, "interrupted")
})

test_that("primary failure wins when cleanup and reporting also fail", {
  cleaned <- 0L
  condition <- tryCatch(
    crossform:::.execute_guarded(
      compute = function() stop("primary computation failed"),
      receipt = receipt_fixture(),
      reporter = function(event) stop("observer failed"),
      cleanup = function() {
        cleaned <<- cleaned + 1L
        stop("cleanup failed")
      }
    ),
    effect_execution_error = identity
  )

  expect_identical(cleaned, 1L)
  expect_match(condition$message, "primary computation failed")
  expect_match(condition$cleanup_status$message, "cleanup failed")
  expect_false(condition$cleanup_status$success)
  expect_identical(condition$receipt$completion_status, "failed")
  expect_gte(length(condition$observer_failures), 2L)
})

test_that("success cleanup and completion reporting occur exactly once", {
  cleaned <- 0L
  final_receipt <- NULL
  value <- crossform:::.execute_guarded(
    compute = function() 17,
    receipt = receipt_fixture(),
    reporter = function(event) {
      if (event$type == "complete") final_receipt <<- event$receipt
    },
    cleanup = function() cleaned <<- cleaned + 1L
  )

  expect_identical(value, 17)
  expect_identical(cleaned, 1L)
  expect_identical(final_receipt$completion_status, "complete")
  expect_identical(final_receipt$completed_task_count, final_receipt$task_count)
  expect_gte(final_receipt$elapsed_seconds, 0)
})

test_that("cleanup failure cannot replace interrupt semantics", {
  cleaned <- 0L
  interrupt <- structure(list(message = "stop now", call = NULL),
    class = c("interrupt", "condition"))
  condition <- tryCatch(
    crossform:::.execute_guarded(
      compute = function() signalCondition(interrupt),
      receipt = receipt_fixture(),
      reporter = function(event) stop("observer failed"),
      cleanup = function() {
        cleaned <<- cleaned + 1L
        stop("cleanup failed")
      }
    ),
    interrupt = identity
  )

  expect_s3_class(condition, "interrupt")
  expect_identical(cleaned, 1L)
  expect_match(condition$message, "stop now")
  expect_false(condition$cleanup_status$success)
  expect_identical(condition$receipt$completion_status, "interrupted")
})
