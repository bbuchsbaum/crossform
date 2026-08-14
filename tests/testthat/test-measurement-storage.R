test_that("memory measurement storage round-trips every logical block", {
  fixture <- measurement_kernel_fixture(seed = 2026081211)
  run <- effectagram:::.run_measurement_contraction(
    fixture$task, compute_policy(block_features = 2L), route = "pull_h"
  )
  plan <- effectagram:::.measurement_result_plan(fixture$task)
  store <- effectagram:::.memory_measurement_store(
    run$blocks, plan$block_index
  )

  expect_s3_class(store, "effect_measurement_store")
  expect_identical(store$codec, "measurement-block-double-v1")
  expect_true(store$manifest$complete)
  expect_identical(store$index$left, fixture$measurement_edges$left)
  expect_identical(store$index$right, fixture$measurement_edges$right)
  for (edge in seq_along(run$blocks)) {
    expect_equal(
      effectagram:::.read_measurement_store(store, edge),
      run$blocks[[edge]], tolerance = 0
    )
    expect_equal(
      effectagram:::.read_measurement_store(store, names(run$blocks)[[edge]]),
      run$blocks[[edge]], tolerance = 0
    )
  }
})

test_that("block-backed measurement storage is complete only after every write", {
  fixture <- measurement_kernel_fixture(seed = 2026081212)
  run <- effectagram:::.run_measurement_contraction(
    fixture$task, compute_policy(block_features = 3L), route = "pull_h"
  )
  plan <- effectagram:::.measurement_result_plan(fixture$task)
  path <- tempfile(fileext = ".emf")
  on.exit(unlink(c(path, paste0(path, ".manifest.rds"))), add = TRUE)
  store <- effectagram:::.file_measurement_store(
    path, plan$block_index, create = TRUE
  )

  expect_false(store$manifest$complete)
  expect_error(
    effectagram:::.validate_measurement_store(store),
    "incomplete"
  )
  written <- effectagram:::.write_measurement_block(
    store, 1L, run$blocks[[1L]]
  )
  expect_error(
    effectagram:::.validate_measurement_store(store, require_complete = FALSE),
    "stale"
  )
  expect_false(written$manifest$complete)
  for (edge in 2:length(run$blocks)) {
    written <- effectagram:::.write_measurement_block(
      written, edge, run$blocks[[edge]]
    )
  }

  expect_true(written$manifest$complete)
  expect_silent(effectagram:::.validate_measurement_store(written))
  expect_identical(file.info(path)$size,
    sum(plan$block_index$length_elements) * 8)
  for (edge in seq_along(run$blocks)) {
    expect_equal(
      effectagram:::.read_measurement_store(written, edge),
      run$blocks[[edge]], tolerance = 0
    )
  }
})

test_that("measurement stores reject overwrite, malformed blocks, and indices", {
  fixture <- measurement_kernel_fixture(seed = 2026081213)
  run <- effectagram:::.run_measurement_contraction(
    fixture$task, route = "pull_h"
  )
  plan <- effectagram:::.measurement_result_plan(fixture$task)
  path <- tempfile(fileext = ".emf")
  on.exit(unlink(c(path, paste0(path, ".manifest.rds"))), add = TRUE)
  store <- effectagram:::.file_measurement_store(
    path, plan$block_index, create = TRUE
  )

  expect_error(effectagram:::.file_measurement_store(
    path, plan$block_index, create = TRUE
  ), "overwrite")
  expect_error(effectagram:::.write_measurement_block(
    store, 1L, matrix(1, 1, 1)
  ), "dimensions")
  expect_error(effectagram:::.read_measurement_store(store, "missing"),
    "incomplete|not present")

  forged <- plan$block_index
  forged$left_axis[[1L]] <- forged$right_axis[[1L]]
  expect_error(effectagram:::.memory_measurement_store(run$blocks, forged),
    "index|noncanonical|identity")
})

test_that("storage codec never upgrades a partial edge result", {
  fixture <- measurement_kernel_fixture(seed = 2026081214)
  plan <- effectagram:::.measurement_result_plan(fixture$task)
  path <- tempfile(fileext = ".emf")
  on.exit(unlink(c(path, paste0(path, ".manifest.rds"))), add = TRUE)
  store <- effectagram:::.file_measurement_store(
    path, plan$block_index, create = TRUE
  )
  forged <- store
  forged$manifest$complete <- TRUE
  forged$manifest$written[] <- TRUE

  expect_error(effectagram:::.validate_measurement_store(forged),
    "stale|inconsistent")
  expect_identical(store$codec, "measurement-block-double-v1")
  expect_false(store$manifest$complete)
})
