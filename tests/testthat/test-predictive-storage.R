pg_storage_fixture <- function(n = 9L) {
  f <- pg_fixture()
  f$at <- additive_frame(members = rep(list(1:2, 1L, 2L), length.out = n),
    domain = f$domain, measurements = paste0("node", seq_len(n)))
  f$plan <- pg_plan(f, "a", "b")
  f
}

test_that("[T25 T43] prediction and signed/mode stores are consumed only in bounded row blocks", {
  f <- pg_storage_fixture()
  reads <- list(); materialize <- crossform:::.geometry_fit_materialize
  guarded <- function(store) {
    original <- store$read
    store$read <- function(rows = NULL) {
      if (is.null(rows) || length(rows) > 2L) stop("unbounded predictive row read")
      reads[[length(reads) + 1L]] <<- rows
      original(rows)
    }
    store
  }
  testthat::local_mocked_bindings(.geometry_fit_materialize = function(...) {
    form <- materialize(...)
    form$total <- guarded(form$total); form$coherent <- guarded(form$coherent)
    form
  }, .package = "crossform")
  path <- tempfile("bounded-fit-"); score_path <- tempfile("bounded-score-")
  withr::defer(unlink(c(path, score_path), recursive = TRUE))
  fit <- fit_geometry(f$plan, f$basis, rank = 2, penalty = .2, row_block = 2,
    storage = "block", storage_path = path)
  fit$values <- guarded(fit$values)
  result <- score_geometry(fit, pg_plan(f), modes = TRUE, components = TRUE, row_block = 2,
    storage = "block", storage_path = score_path)
  result$mode_values <- guarded(result$mode_values)
  expect_equal(nrow(as.data.frame(result, view = "modes")), sum(result$table$rank))
  expect_true(length(reads) > nrow(result$table))
  expect_lte(max(lengths(reads)), 2L)
  expect_identical(fit$values$dim, c(9L, 10L))
  expect_identical(result$workspace$row_block, 2L)
})

test_that("[T43 T57] byte admission includes retained predictions and selects smaller row blocks", {
  f <- pg_storage_fixture()
  prepare <- crossform:::.geometry_fit_prepare(f$plan, f$basis)
  small <- crossform:::.geometry_prediction_workspace(prepare, 2, "memory", 1, "fit")$workspace
  large <- crossform:::.geometry_prediction_workspace(prepare, 2, "memory", 9, "fit")$workspace
  expect_gt(large$planned_workspace_bytes, small$planned_workspace_bytes)
  budget <- (small$planned_workspace_bytes + large$planned_workspace_bytes) / 2
  policy <- compute_policy(workspace_bytes = budget)
  limited <- pg_plan(f, "a", "b", compute = policy)
  fit <- fit_geometry(limited, f$basis, rank = 2, penalty = .2, row_block = 9)
  expect_lt(fit$workspace$row_block, 9L)
  expect_lte(fit$workspace$planned_workspace_bytes, budget)
  expect_identical(fit$workspace$workers, 1L)
  expect_equal(fit$workspace$planned_workspace_bytes,
    fit$workspace$extra_bytes + fit$receipt$memory$planned_workspace_bytes)
  reference <- fit_geometry(f$plan, f$basis, rank = 2, penalty = .2)
  expect_identical(fit$prediction_id, reference$prediction_id)
  expect_equal(pg_forms(fit), pg_forms(reference), tolerance = 1e-11)
  # A byte limit too small even for a single row refuses before execution.
  testthat::local_mocked_bindings(.geometry_fit_materialize = function(...) stop("early read"), .package = "crossform")
  refused <- catch_refusal(fit_geometry(pg_plan(f, "a", "b", compute = compute_policy(workspace_bytes = 1)),
    f$basis, rank = 2, penalty = .2))
  expect_identical(refused$reasons[[1]], "predictive_workspace_budget_exceeded")
  refused <- catch_refusal(score_geometry(fit, pg_plan(f, compute = compute_policy(workspace_bytes = 1))))
  expect_identical(refused$reasons[[1]], "predictive_workspace_budget_exceeded")
  changed <- fit; changed$workspace$planned_workspace_bytes <- 0
  expect_error(as.data.frame(changed), class = "effect_contract_error")
})

test_that("[T44 T57] admitted feature/row/storage routes agree and fixed-metric block output refuses", {
  f <- pg_storage_fixture(4)
  reference <- fit_geometry(f$plan, f$basis, rank = 2, penalty = .2)
  evidence <- score_geometry(reference, pg_plan(f), modes = TRUE)
  for (features in 1:2) for (storage in c("memory", "block")) {
    paths <- c(tempfile("route-fit-"), tempfile("route-score-"))
    withr::defer(unlink(paths, recursive = TRUE))
    fit <- fit_geometry(pg_plan(f, "a", "b", compute = compute_policy(block_features = features)),
      f$basis, rank = 2, penalty = .2, storage = storage,
      storage_path = if (storage == "block") paths[1] else NULL, row_block = features)
    score <- score_geometry(fit, pg_plan(f, compute = compute_policy(block_features = features)),
      modes = TRUE, storage = storage, storage_path = if (storage == "block") paths[2] else NULL,
      row_block = features)
    expect_equal(pg_forms(fit), pg_forms(reference), tolerance = 1e-11)
    expect_equal(as.data.frame(score), as.data.frame(evidence), tolerance = 1e-11)
    expect_equal(as.data.frame(score, view = "modes"), as.data.frame(evidence, view = "modes"), tolerance = 1e-11)
  }
  fixed <- neural_metric(matrix(c(2, .3, .3, 1), 2), f$domain)
  train <- pg_plan(f, "a", "b", metric = fixed)
  fit <- fit_geometry(train, f$basis, rank = 2, penalty = .2)
  score <- score_geometry(fit, pg_plan(f, metric = fixed))
  expect_true(all(is.finite(score$table$gain)))
  testthat::local_mocked_bindings(.geometry_fit_materialize = function(...) stop("early read"), .package = "crossform")
  path <- tempfile("unsupported-route-")
  expect_identical(catch_refusal(fit_geometry(train, f$basis, rank = 2, penalty = .2,
    storage = "block", storage_path = path))$reasons[[1]], "fixed_metric_block_storage_not_admitted")
  expect_false(file.exists(path))
  expect_identical(catch_refusal(compute_policy(workers = 2))$reasons[[1]], "worker_pool_not_implemented")
})

test_that("[T44] failures and interrupted writes leave no completed fit or score", {
  f <- pg_storage_fixture()
  reference <- fit_geometry(f$plan, f$basis, rank = 2, penalty = .2)
  # Failure after at least one prediction tile has been written.
  write_tile <- crossform:::.write_geometry_tile
  failures <- 0L
  testthat::local_mocked_bindings(.write_geometry_tile = function(store, rows, coordinates, value) {
    if (basename(store$path) %in% c("prediction.bin", "modes.bin") && min(rows) > 1L) {
      failures <<- failures + 1L
      stop("injected partial write")
    }
    write_tile(store, rows, coordinates, value)
  }, .package = "crossform")
  for (kind in c("fit", "score")) {
    path <- tempfile("incomplete-prediction-")
    withr::defer(unlink(path, recursive = TRUE))
    call <- if (kind == "fit") function() fit_geometry(f$plan, f$basis, rank = 2, penalty = .2,
      storage = "block", storage_path = path, row_block = 1) else function() score_geometry(reference,
        pg_plan(f), storage = "block", storage_path = path, row_block = 1, modes = TRUE)
    expect_error(call(), "injected partial write")
    expect_false(dir.exists(path))
    expect_false(file.exists(file.path(path, paste0(kind, ".rds"))))
  }
  expect_identical(failures, 2L)
})

test_that("[T44] read errors and R interrupts clean only the owned output directory", {
  f <- pg_storage_fixture(); fit <- fit_geometry(f$plan, f$basis, rank = 2, penalty = .2)
  materialize <- crossform:::.geometry_fit_materialize
  testthat::local_mocked_bindings(.geometry_fit_materialize = function(...) {
    form <- materialize(...)
    original <- form$total$read
    form$total$read <- function(rows) {
      if (min(rows) > 1) stop("injected read failure")
      original(rows)
    }
    form
  }, .package = "crossform")
  path <- tempfile("read-failure-"); withr::defer(unlink(path, recursive = TRUE))
  expect_error(score_geometry(fit, pg_plan(f), storage = "block", storage_path = path,
    row_block = 1), "injected read failure")
  expect_false(file.exists(path))
  existing <- tempfile("owned-by-caller-"); dir.create(existing)
  marker <- file.path(existing, "keep"); writeLines("keep", marker)
  withr::defer(unlink(existing, recursive = TRUE))
  expect_error(score_geometry(fit, pg_plan(f), storage = "block", storage_path = existing), "overwrite")
  expect_identical(readLines(marker), "keep")
  testthat::local_mocked_bindings(.geometry_score_from_form = function(...) {
    stop(structure(list(message = "injected interrupt", call = NULL), class = c("interrupt", "condition")))
  }, .package = "crossform")
  interrupted <- tryCatch(score_geometry(fit, pg_plan(f), storage = "block", storage_path = path),
    interrupt = function(e) e)
  expect_s3_class(interrupted, "interrupt")
  expect_false(file.exists(path))
})
