pg_crossfit_plan <- function(f, left, right, weight = NULL, directed = FALSE) {
  plan_geometry(f$rel, f$at, pairing(left, right, weight = weight, directed = directed,
    independence = "independent", generalizes_over = "run"))
}

test_that("[T26 T30 T48] cross-fit schedules use only declared disjoint products", {
  f <- pg_fixture()
  edges <- t(combn(letters[1:4], 2))
  complete <- pg_crossfit_plan(f, edges[, 1], edges[, 2], 1:6)
  schedule <- crossform:::.geometry_crossfit_schedule(complete, f$basis)
  expect_length(schedule$folds, 6L)
  expect_equal(schedule$edges$weight, (1:6) / 21)
  for (fold in schedule$folds) {
    expect_length(intersect(fold$training_support$observations, fold$evaluation_support$observations), 0L)
    expect_equal(nrow(fold$training$pairing), 1L)
    expected <- complete$pairing$weight[
      complete$pairing$left == fold$evaluation$pairing$left & complete$pairing$right == fold$evaluation$pairing$right]
    expect_equal(fold$weight, expected)
  }
  disconnected <- pg_crossfit_plan(f, c("a", "c"), c("b", "d"), c(2, 5))
  expect_length(crossform:::.geometry_crossfit_schedule(disconnected, f$basis)$folds, 2L)
  # A zero-weight disjoint edge supplies no independent training data.
  zero <- pg_crossfit_plan(f, c("a", "a", "c"), c("b", "c", "d"), c(1, 1, 0))
  star <- pg_crossfit_plan(f, rep("a", 3), c("b", "c", "d"))
  three <- pg_crossfit_plan(f, c("a", "a", "b"), c("b", "c", "c"))
  testthat::local_mocked_bindings(fit_geometry = function(...) stop("source execution began"), .package = "crossform")
  for (bad in list(zero, star, three)) {
    refusal <- catch_refusal(crossform:::.geometry_crossfit(bad, f$basis, rank = 1, penalty = .2))
    expect_identical(refusal$reasons[[1]], "insufficient_disjoint_edges")
  }
  directed <- pg_crossfit_plan(f, c("a", "c"), c("b", "d"), directed = TRUE)
  expect_identical(catch_refusal(crossform:::.geometry_crossfit(directed, f$basis, rank = 1, penalty = .2))$reasons[[1]],
    "directed_full_form_not_admitted")
})

test_that("[T49] manual folds match gain, unequal weights and separate squared offsets", {
  f <- pg_fixture()
  edges <- t(combn(letters[1:4], 2))
  plan <- pg_crossfit_plan(f, edges[, 1], edges[, 2], 1:6)
  result <- crossform:::.geometry_crossfit(plan, f$basis, rank = 2, penalty = .2)
  gains <- norms <- matrix(0, 2, 6)
  mean_forms <- replicate(2, matrix(0, 4, 4), simplify = FALSE)
  for (i in 1:6) {
    test_parts <- edges[i, ]; train_parts <- setdiff(letters[1:4], test_parts)
    train <- pg_plan(f, train_parts[1], train_parts[2])
    test <- pg_plan(f, test_parts[1], test_parts[2])
    fit <- fit_geometry(train, f$basis, rank = 2, penalty = .2)
    score <- score_geometry(fit, test)
    gains[, i] <- score$table$gain; norms[, i] <- score$table$prediction_norm_sq
    for (j in 1:2) mean_forms[[j]] <- mean_forms[[j]] + plan$pairing$weight[i] * pg_forms(fit)[[j]]
    expect_equal(result$scores[[i]]$table, score$table, tolerance = 1e-11)
    expect_identical(result$scores[[i]]$prediction_id, fit$prediction_id)
  }
  expect_equal(result$table$gain, as.numeric(gains %*% ((1:6)/21)), tolerance = 1e-11)
  expect_equal(result$table$prediction_norm_sq, as.numeric(norms %*% ((1:6)/21)), tolerance = 1e-11)
  expect_gt(max(abs(result$table$prediction_norm_sq - vapply(mean_forms, function(F) sum(F^2), 0))), .001)
  expect_equal(result$table$gain, 2 * result$table$inner_product - result$table$prediction_norm_sq)
  expect_match(result$estimand, "fold training sizes")
  reverse <- pg_crossfit_plan(f, edges[, 2], edges[, 1], 1:6)
  reversed <- crossform:::.geometry_crossfit(reverse, f$basis, rank = 2, penalty = .2)
  expect_equal(reversed$table, result$table, tolerance = 1e-11)
  wrong <- result; wrong$table$prediction_norm_sq <- vapply(mean_forms, function(F) sum(F^2), 0)
  expect_error(crossform:::.validate_geometry_crossfit(wrong), "squared norm")
})

test_that("[T44 T49] blocked cross-fit reopens and a failed fold has no complete aggregate", {
  f <- pg_fixture()
  plan <- pg_crossfit_plan(f, c("a", "c"), c("b", "d"), c(2, 5))
  reference <- crossform:::.geometry_crossfit(plan, f$basis, rank = 2, penalty = .2)
  path <- tempfile("crossfit-"); withr::defer(unlink(path, recursive = TRUE))
  blocked <- crossform:::.geometry_crossfit(plan, f$basis, rank = 2, penalty = .2,
    storage = "block", storage_path = path, row_block = 1)
  reopened <- readRDS(file.path(path, "crossfit.rds"))
  expect_equal(reopened$table, reference$table, tolerance = 1e-11)
  expect_s3_class(crossform:::.validate_geometry_crossfit(reopened), "effect_geometry_crossfit")
  expect_equal(sum(grepl("-fit$", list.dirs(path, recursive = FALSE))), 2L)
  expect_equal(pg_forms(reopened$fits[[1]]), pg_forms(reference$fits[[1]]), tolerance = 1e-11)
  expect_s3_class(crossform:::.validate_geometry_crossfit(reopened), "effect_geometry_crossfit")
  score <- score_geometry; called <- 0L
  testthat::local_mocked_bindings(score_geometry = function(...) {
    called <<- called + 1L
    if (called == 2L) stop("second fold failed")
    score(...)
  }, .package = "crossform")
  bad_path <- tempfile("crossfit-incomplete-"); withr::defer(unlink(bad_path, recursive = TRUE))
  expect_error(crossform:::.geometry_crossfit(plan, f$basis, rank = 2, penalty = .2,
    storage = "block", storage_path = bad_path), "second fold failed")
  expect_false(file.exists(bad_path))
  expect_identical(called, 2L)
})
