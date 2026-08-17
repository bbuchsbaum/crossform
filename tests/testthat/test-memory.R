test_that("workspace plans account for persistent and active-task categories", {
  plan <- memory_plan(
    frame_bytes = 100,
    resident_source_bytes = 10,
    source_handle_bytes = 20,
    source_block_bytes = 30,
    relation_block_bytes = 40,
    atom_block_bytes = 50,
    local_state_bytes = 60,
    output_bytes = 70,
    contraction_bytes = 80,
    replacement_copy_bytes = 90,
    serialization_overlap_bytes = 100,
    reorder_buffer_bytes = 110,
    checkpoint_buffer_bytes = 120,
    workers = 4, n_active = 4, safety_factor = 1
  )

  expect_named(plan$categories, c(
    "frame", "resident_source", "source_handles", "source_block",
    "relation_block", "atom_block", "local_state", "output", "contraction",
    "replacement_copy", "serialization_overlap", "reorder_buffer",
    "checkpoint_buffer"
  ))
  expect_equal(plan$persistent_workspace_bytes, 260)
  expect_equal(plan$task_workspace_per_active_bytes, 620)
  expect_equal(plan$active_task_workspace_bytes, 2480)
  expect_equal(plan$modeled_workspace_bytes, 2740)
  expect_equal(plan$planned_workspace_bytes, 2740)
})

test_that("only per-task workspace scales with active tasks", {
  one <- memory_plan(frame_bytes = 100, source_block_bytes = 50,
    workers = 4, n_active = 1, safety_factor = 1)
  four <- memory_plan(frame_bytes = 100, source_block_bytes = 50,
    workers = 4, n_active = 4, safety_factor = 1)

  expect_equal(one$persistent_workspace_bytes, four$persistent_workspace_bytes)
  expect_equal(one$active_task_workspace_bytes, 50)
  expect_equal(four$active_task_workspace_bytes, 200)
})

test_that("workspace budgets and RSS observations remain distinct", {
  plan <- memory_plan(output_bytes = 400, contraction_bytes = 200,
    replacement_copy_bytes = 200, safety_factor = 1.25,
    budget_bytes = 1100, measured_workspace_bytes = 900,
    baseline_rss_bytes = 5000, peak_rss_bytes = 6200)

  expect_equal(plan$planned_workspace_bytes, 1000)
  expect_true(plan$fits_budget)
  expect_true(plan$measured_workspace_within_plan)
  expect_equal(plan$baseline_rss_bytes, 5000)
  expect_equal(plan$incremental_peak_rss_bytes, 1200)
  expect_equal(plan$absolute_peak_rss_bytes, 6200)
  expect_identical(plan$prediction_kind,
    "crossform_owned_workspace_upper_bound")
})

test_that("memory plans reject invalid and overflowing inputs", {
  expect_error(memory_plan(output_bytes = -1), "nonnegative",
    class = "effect_input_error")
  expect_error(memory_plan(workers = 1.5), "positive finite whole scalar",
    class = "effect_input_error")
  expect_error(memory_plan(workers = 2, n_active = 3), "cannot exceed",
    class = "effect_input_error")
  expect_error(memory_plan(safety_factor = 0.9), "greater than or equal",
    class = "effect_input_error")
  expect_error(memory_plan(budget_bytes = 0), "positive",
    class = "effect_input_error")
  expect_error(memory_plan(output_bytes = 0.5), "whole scalar",
    class = "effect_input_error")
  expect_error(memory_plan(source_block_bytes = 2^52, workers = 3,
    n_active = 3), "overflows", class = "effect_invariant_error")
  expect_error(memory_plan(baseline_rss_bytes = 100, peak_rss_bytes = 99),
    "cannot be smaller", class = "effect_input_error")
})

test_that("default plans contain no invented runtime RSS reserve", {
  plan <- memory_plan(safety_factor = 1)

  expect_equal(plan$modeled_workspace_bytes, 0)
  expect_null(plan$baseline_rss_bytes)
  expect_null(plan$absolute_peak_rss_bytes)
})

test_that("dense matrix byte planning is exact and allocation-free in size", {
  shapes <- list(c(0, 0), c(1, 1), c(4, 2), c(317, 41))
  for (shape in shapes) {
    expected <- as.double(utils::object.size(
      matrix(0, shape[[1L]], shape[[2L]])
    ))
    expect_identical(
      crossform:::.dense_double_matrix_bytes(shape[[1L]], shape[[2L]]),
      expected
    )
  }

  # This represents a 20 GB matrix.  The assertion is safe precisely because
  # the planner performs arithmetic and never constructs the matrix.
  expect_identical(
    crossform:::.dense_double_matrix_bytes(50000, 50000),
    20000000216
  )
  expect_error(
    crossform:::.dense_double_matrix_bytes(2^52, 2^52),
    "overflows"
  , class = "effect_invariant_error")
})

test_that("named array byte planning preserves metadata without values", {
  dimensions <- c(31, 4, 3)
  names <- list(NULL, paste0("effect_", seq_len(4)), paste0("run_", seq_len(3)))
  expected <- as.double(utils::object.size(array(
    0, dimensions, dimnames = names
  )))

  expect_identical(
    crossform:::.named_double_array_bytes(dimensions, names),
    expected
  )
})
