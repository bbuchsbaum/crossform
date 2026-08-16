test_that("version 0.1 accepts exactly one worker", {
  policy <- compute_policy(workers = 1)

  expect_identical(policy$workers, 1L)
  expect_identical(policy$process_backend, "sequential")
  expect_error(compute_policy(workers = 2), "not implemented in crossform 0.1")
  expect_error(compute_policy(workers = 0), "must be 1")
  expect_error(compute_policy(workers = Inf), "must be 1")
})

test_that("invalid worker policy fails before source inspection", {
  accessed <- 0L
  inspect_source <- function() {
    accessed <<- accessed + 1L
    list(block_read = TRUE)
  }
  invalid <- structure(
    list(
      workers = 4,
      block_features = NULL,
      workspace_bytes = NULL,
      process_backend = "sequential"
    ),
    class = "effect_compute_policy"
  )

  expect_error(
    crossform:::.execution_preflight(invalid, inspect_source),
    "not implemented in crossform 0.1"
  )
  expect_identical(accessed, 0L)
})

test_that("valid preflight is owned and sequential", {
  accessed <- 0L
  plan <- crossform:::.execution_preflight(
    compute_policy(block_features = 64, workspace_bytes = 1e7),
    function() {
      accessed <<- accessed + 1L
      list(reopenable = FALSE)
    }
  )

  expect_identical(accessed, 1L)
  expect_identical(plan$compute$block_features, 64L)
  expect_identical(plan$compute$workspace_bytes, 1e7)
  expect_identical(plan$workers, 1L)
  expect_identical(plan$process_backend, "sequential")
  expect_false(plan$source_capabilities$reopenable)
})

test_that("compute policies are canonicalized and revalidated at preflight", {
  policy <- compute_policy(block_features = 7, workspace_bytes = 1024)
  expect_identical(policy$block_features, 7L)

  mutated <- policy
  mutated$block_features <- 0
  expect_error(
    crossform:::.execution_preflight(mutated, function() stop("must not inspect")),
    "positive integer"
  )

  forged <- policy
  forged$extra <- TRUE
  expect_error(
    crossform:::.execution_preflight(forged, function() stop("must not inspect")),
    "canonical order"
  )
})
