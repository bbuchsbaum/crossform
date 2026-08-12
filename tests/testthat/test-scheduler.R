test_that("simulator bounds inflight work and reserved result memory", {
  values <- lapply(seq_len(8), function(i) matrix(i / 10, 2, 2))
  ids <- seq_len(8)
  got <- effectagram:::.simulate_bounded_schedule(
    values, ids, result_bytes = rep(100, 8),
    completion_order = rev(ids), max_inflight = 3, max_reorder_bytes = 300
  )
  oracle <- Reduce(`+`, values)

  expect_identical(got$value, oracle)
  expect_lte(got$diagnostics$max_inflight_seen, 3)
  expect_lte(got$diagnostics$max_reserved_bytes, 300)
  expect_gt(got$diagnostics$backpressure_count, 0)
  expect_identical(got$diagnostics$release_order, ids)
})
test_that("adversarial completion orders preserve canonical reduction", {
  values <- list(1e16, 1, -1e16, 3, 4)
  ids <- seq_along(values)
  oracle <- effectagram:::.simulate_bounded_schedule(
    values, ids, rep(8, 5), ids, 3, 40
  )$value

  set.seed(99)
  for (iteration in seq_len(25)) {
    order <- sample(ids)
    got <- effectagram:::.simulate_bounded_schedule(
      values, ids, rep(8, 5), order, 3, 40
    )$value
    expect_identical(got, oracle)
  }
})

test_that("simulated failure releases all scheduler-owned state", {
  condition <- tryCatch(
    effectagram:::.simulate_bounded_schedule(
      as.list(1:5), 1:5, rep(10, 5), rev(1:5), 3, 30,
      fail_task_id = 3L
    ),
    effect_scheduler_error = identity
  )

  expect_s3_class(condition, "effect_scheduler_error")
  expect_true(condition$cleanup_status$success)
  expect_gte(condition$cleanup_status$inflight_released, 0)
  expect_gte(condition$cleanup_status$buffered_released, 0)
})

test_that("tasks larger than the reorder reservation are rejected before dispatch", {
  expect_error(
    effectagram:::.simulate_bounded_schedule(list(1), 1, 101, 1, 1, 100),
    "cannot be dispatched"
  )
})
