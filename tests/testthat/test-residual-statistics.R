direct_pair_values <- function(fit, partition, statistics) {
  residuals <- residual_block(
    fit, partition, seq_len(fit$relation$n_features)
  )
  product <- crossprod(residuals)
  product[cbind(statistics$pair_i, statistics$pair_j)]
}

test_that("fixed-width residual reads pad before residualization", {
  fixture <- residual_statistics_fixture()
  ordinary <- residual_block(fixture$fit, "run1", 17:19)
  padded <- crossform:::.residual_padded_block(
    fixture$fit, "run1", 17:19, 8L
  )

  expect_identical(dim(padded), c(37L, 8L))
  expect_equal(padded[, 1:3, drop = FALSE], ordinary, tolerance = 2e-13)
  expect_identical(padded[, 4:8, drop = FALSE], matrix(0, 37L, 5L))
})

test_that("residual pair statistics are bitwise invariant to workspace", {
  fixture <- residual_statistics_fixture()
  budgets <- residual_statistics_budgets(fixture)
  narrow <- residual_pair_statistics(
    fixture$fit, fixture$frame, workspace_bytes = budgets$minimum
  )
  narrow_reads <- fixture$reads()
  fixture$reset_reads()
  wide <- residual_pair_statistics(
    fixture$fit, fixture$frame, partitions = c("run3", "run1", "run2"),
    workspace_bytes = budgets$wider
  )
  wide_reads <- fixture$reads()

  expect_identical(narrow$pair_i, wide$pair_i)
  expect_identical(narrow$pair_j, wide$pair_j)
  expect_identical(narrow$partitions, wide$partitions)
  expect_identical(narrow$atomic, wide$atomic)
  expect_identical(narrow$signature, wide$signature)
  expect_gt(narrow$execution$memory$cache_capacity,
    0L)
  expect_gt(wide$execution$memory$cache_capacity,
    narrow$execution$memory$cache_capacity)
  expect_false(identical(narrow$execution, wide$execution))
  expect_true(all(narrow_reads >= wide_reads))
  expect_identical(
    unname(narrow_reads),
    unname(vapply(narrow$execution$atomic, `[[`, integer(1), "residual_reads"))
  )
  expect_identical(
    unname(wide_reads),
    unname(vapply(wide$execution$atomic, `[[`, integer(1), "residual_reads"))
  )
  expect_silent(crossform:::.validate_residual_pair_statistics(
    narrow, deep = TRUE
  ))
})

test_that("canonical pair statistics agree with a direct residual oracle", {
  fixture <- residual_statistics_fixture()
  budget <- residual_statistics_budgets(fixture)$wider
  statistics <- residual_pair_statistics(
    fixture$fit, fixture$frame, workspace_bytes = budget
  )

  for (partition in statistics$partitions) {
    expect_equal(
      statistics$atomic[[partition]]$cross_products,
      direct_pair_values(fixture$fit, partition, statistics),
      tolerance = 8e-13
    )
    expect_identical(
      statistics$atomic[[partition]]$residual_df,
      residual_df(fixture$fit, partition)
    )
  }
})

test_that("atomic statistics combine without rereading residuals", {
  fixture <- residual_statistics_fixture()
  budget <- residual_statistics_budgets(fixture)$wider
  statistics <- residual_pair_statistics(
    fixture$fit, fixture$frame, workspace_bytes = budget
  )
  fixture$reset_reads()
  scope <- crossform:::.residual_pair_scope(
    statistics, c("run3", "run1")
  )
  covariance <- crossform:::.residual_pair_scope_matrix(
    statistics, c("run3", "run1")
  )

  expect_identical(scope$partitions, c("run1", "run3"))
  expect_identical(
    scope$cross_products,
    statistics$atomic$run1$cross_products +
      statistics$atomic$run3$cross_products
  )
  expect_identical(scope$residual_df,
    statistics$atomic$run1$residual_df +
      statistics$atomic$run3$residual_df)
  expect_identical(fixture$reads(), c(run1 = 0L, run2 = 0L, run3 = 0L))
  direct <- (direct_pair_values(fixture$fit, "run1", statistics) +
    direct_pair_values(fixture$fit, "run3", statistics)) /
    scope$residual_df
  expect_equal(
    covariance[cbind(statistics$pair_i, statistics$pair_j)],
    direct, tolerance = 8e-13
  )
})

test_that("derived local precision inherits workspace bitwise invariance", {
  fixture <- residual_statistics_fixture()
  budgets <- residual_statistics_budgets(fixture)
  narrow <- residual_pair_statistics(
    fixture$fit, fixture$frame, workspace_bytes = budgets$minimum
  )
  wide <- residual_pair_statistics(
    fixture$fit, fixture$frame, workspace_bytes = budgets$wider
  )
  support <- crossform:::.support_index_support(
    fixture$frame$support_index, 10L
  )[[1L]]
  narrow_covariance <- as.matrix(
    crossform:::.residual_pair_scope_matrix(narrow)[support, support]
  )
  wide_covariance <- as.matrix(
    crossform:::.residual_pair_scope_matrix(wide)[support, support]
  )
  regularize <- function(value) {
    target <- diag(diag(value), nrow(value))
    solve(0.8 * value + 0.2 * target + diag(1e-8, nrow(value)))
  }

  expect_identical(narrow_covariance, wide_covariance)
  expect_identical(regularize(narrow_covariance),
    regularize(wide_covariance))
})

test_that("residual pair preflight refuses before reading", {
  fixture <- residual_statistics_fixture()
  minimum <- residual_statistics_budgets(fixture)$minimum
  expect_error(
    residual_pair_statistics(
      fixture$fit, fixture$frame, workspace_bytes = minimum - 1
    ),
    "requires at least.*workspace budget"
  )
  expect_identical(fixture$reads(), c(run1 = 0L, run2 = 0L, run3 = 0L))
})

test_that("residual pair statistics require an explicit error and support channel", {
  fixture <- residual_statistics_fixture()
  pure <- relation_fit(fixture$fit$relation)
  unsupported <- compile_frame(voxelwise(), fixture$domain)

  expect_error(
    residual_pair_statistics(pure, fixture$frame),
    "learned_metric_input"
  )
  expect_error(
    residual_pair_statistics(fixture$fit, unsupported),
    "explicit support index"
  )
})
