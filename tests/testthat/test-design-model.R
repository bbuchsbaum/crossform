test_that("named row_ids bind to designs by partition name", {
  designs <- list(
    `run-1` = matrix(1, 2L, 1L),
    `run-2` = matrix(1, 3L, 1L)
  )
  rows <- crossform:::.normalize_design_rows(
    designs, list(`run-2` = 11:13, `run-1` = 1:2)
  )
  expect_identical(rows$row_ids, list(`run-1` = 1:2, `run-2` = 11:13))
  expect_identical(rownames(rows$designs[["run-2"]]), c("11", "12", "13"))

  positional <- crossform:::.normalize_design_rows(designs, list(1:2, 11:13))
  expect_identical(positional$row_ids, rows$row_ids)

  expect_error(crossform:::.normalize_design_rows(
    designs, list(`run-1` = 1:2, `run-3` = 11:13)
  ), "Named `row_ids`", class = "effect_input_error")
  expect_error(crossform:::.normalize_design_rows(
    designs, list(`run-2` = 11:13, 1:2)
  ), "Named `row_ids`", class = "effect_input_error")
})
