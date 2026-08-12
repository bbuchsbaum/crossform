test_that("memory plan accounts for every category without double-counting shared source", {
  plan <- memory_plan(
    shared_source_bytes = 100,
    runtime_reserve_bytes = 0,
    private_source_bytes = 10,
    source_block_bytes = 20,
    relation_block_bytes = 30,
    atom_block_bytes = 40,
    output_bytes = 50,
    contraction_bytes = 60,
    replacement_copy_bytes = 70,
    serialization_overlap_bytes = 80,
    reorder_buffer_bytes = 90,
    checkpoint_buffer_bytes = 100,
    workers = 4,
    n_active = 4,
    safety_factor = 1
  )

  expect_named(plan$categories, c(
    "shared_source", "runtime_reserve", "private_source_per_worker", "source_block",
    "relation_block", "atom_block", "output", "contraction",
    "replacement_copy", "serialization_overlap", "reorder_buffer",
    "checkpoint_buffer"
  ))
  expect_equal(plan$worker_private_total_bytes, 40)
  expect_equal(plan$task_private_per_active_bytes, 490)
  expect_equal(plan$active_task_total_bytes, 1960)
  expect_equal(plan$modeled_peak_bytes, 2150)
  expect_equal(plan$conservative_peak_bytes, 2150)
})

test_that("every task-private category scales with active tasks", {
  arguments <- c(
    "source_block_bytes", "relation_block_bytes", "atom_block_bytes",
    "contraction_bytes", "replacement_copy_bytes",
    "serialization_overlap_bytes", "reorder_buffer_bytes",
    "checkpoint_buffer_bytes"
  )
  for (argument in arguments) {
    one <- do.call(memory_plan, c(setNames(list(100), argument),
      list(runtime_reserve_bytes = 0, workers = 4, n_active = 1,
        safety_factor = 1)))
    four <- do.call(memory_plan, c(setNames(list(100), argument),
      list(runtime_reserve_bytes = 0, workers = 4, n_active = 4,
        safety_factor = 1)))
    expect_equal(one$active_task_total_bytes, 100, info = argument)
    expect_equal(four$active_task_total_bytes, 400, info = argument)
  }
})

test_that("budget and measured peak are separate validations", {
  plan <- memory_plan(
    runtime_reserve_bytes = 0,
    output_bytes = 400,
    contraction_bytes = 200,
    replacement_copy_bytes = 200,
    safety_factor = 1.25,
    budget_bytes = 1100,
    measured_peak_bytes = 1050
  )

  expect_equal(plan$conservative_peak_bytes, 1000)
  expect_true(plan$fits_budget)
  expect_false(plan$measurement_within_plan)
  expect_identical(plan$prediction_kind, "conservative_upper_bound_to_validate")
})

test_that("memory plan rejects invalid and overconfident inputs", {
  expect_error(memory_plan(output_bytes = -1), "nonnegative")
  expect_error(memory_plan(workers = 1.5), "positive finite whole scalar")
  expect_error(memory_plan(workers = 2, n_active = 3), "cannot exceed")
  expect_error(memory_plan(safety_factor = 0.9), "greater than or equal")
  expect_error(memory_plan(budget_bytes = 0), "positive finite")
  expect_error(memory_plan(output_bytes = 0.5), "whole scalar")
  expect_error(memory_plan(private_source_bytes = 2^52, workers = 3),
    "overflows exact")
})

test_that("default plans include an explicit benchmarked runtime reserve", {
  plan <- memory_plan(safety_factor = 1)

  expect_equal(plan$categories[["runtime_reserve"]], 64 * 1024^2)
  expect_equal(plan$modeled_peak_bytes, 64 * 1024^2)
})
