test_that("[T51] exact ties follow rank, penalty and canonical ID, independent of presentation order", {
  f <- pg_fixture()
  iso <- model_basis(kernels = list(m = tcrossprod(f$V)), conditions = f$rel$effect_space, normalize = "trace")
  inputs <- list(list(basis = f$basis, rank = 2L, penalty = 0),
    list(basis = f$basis, rank = 0L, penalty = .1),
    list(basis = f$basis, rank = 0L, penalty = .5),
    list(basis = iso, rank = 0L, penalty = .5))
  candidates <- crossform:::.geometry_candidates(inputs)
  values <- matrix(0, 2, length(candidates))
  chosen <- crossform:::.geometry_select_indices(values, candidates, "global", c(.25, .75))
  eligible <- which(vapply(candidates, function(x) x$rank == 0L && x$penalty == .5, FALSE))
  ids <- vapply(candidates, `[[`, "", "id")
  expect_identical(ids[chosen], rep(sort(ids[eligible], method = "radix")[1], 2))
  expect_identical(crossform:::.geometry_candidates(rev(inputs)), candidates)
  rank_two <- which(vapply(candidates, `[[`, 0L, "rank") == 2L)
  values[1, rank_two] <- 1e-14 # A nonzero gain difference is not an exact tie.
  per_node <- crossform:::.geometry_select_indices(values, candidates, "per_measurement")
  expect_identical(per_node[1], unname(rank_two))
  expect_true(per_node[2] %in% eligible)
  expect_error(crossform:::.geometry_candidates(c(inputs, inputs[1])), "distinct")
  expect_error(crossform:::.geometry_candidates(list(list(basis = f$basis, rank = 1))), "penalty")
})

test_that("[T50 T51] eight-run nested selection matches explicit weighted manual inner fits", {
  f <- pg_nested_fixture()
  for (scope in c("global", "per_measurement")) {
    result <- crossform:::.geometry_nested(f$train, f$test, f$candidates, f$inner, scope = scope)
    candidates <- result$selection$candidates
    expected <- matrix(0, 2, length(candidates))
    actual_train_edges <- as.data.frame(f$train$pairing)
    for (fold in f$inner) {
      plan_for <- function(parts) {
        selected <- actual_train_edges[actual_train_edges$left %in% parts & actual_train_edges$right %in% parts, ]
        plan_geometry(f$rel, f$at, pairing(selected$left, selected$right, selected$weight,
          independence = "independent", generalizes_over = "run"))
      }
      train <- plan_for(fold$train); validation <- plan_for(fold$validate)
      for (j in seq_along(candidates)) {
        candidate <- candidates[[j]]
        fit <- fit_geometry(train, candidate$basis, weights = candidate$weights,
          rank = candidate$rank, penalty = candidate$penalty)
        evidence <- score_geometry(fit, validation)
        expected[, j] <- expected[, j] + fold$weight * evidence$table$gain
      }
    }
    expect_equal(result$selection$mean_inner_gain, expected, tolerance = 1e-11)
    # Independent exact-tie ordering, outside the production selector.
    ranks <- vapply(candidates, `[[`, 0L, "rank"); penalties <- vapply(candidates, `[[`, 0, "penalty")
    ids <- vapply(candidates, `[[`, "", "id")
    select <- function(row) order(-row, ranks, -penalties, ids)[1]
    selected <- if (scope == "global") rep(select(colMeans(expected)), 2) else vapply(1:2, function(i) select(expected[i, ]), 0L)
    expect_identical(result$selection$selected, ids[selected])
    expect_identical(result$selection$scope, scope)
    for (i in 1:2) {
      candidate <- candidates[[selected[i]]]
      manual <- fit_geometry(f$train, candidate$basis, weights = candidate$weights,
        rank = candidate$rank, penalty = candidate$penalty)
      expect_equal(pg_forms(result$fits[[candidate$id]])[[i]], pg_forms(manual)[[i]], tolerance = 1e-11)
      expect_equal(result$table$gain[i], score_geometry(manual, f$test)$table$gain[i], tolerance = 1e-11)
      expect_length(intersect(result$fits[[candidate$id]]$training$observations, c("raw-g", "raw-h")), 0L)
      expect_length(result$fits[[candidate$id]]$training$upstream, 6L)
    }
  }
})

test_that("[T32 T36 T50] every fit and selection excludes outer-test values", {
  f <- pg_nested_fixture(); altered <- pg_nested_fixture(test_scale = -10)
  fitted <- list(); real_fit <- crossform:::.geometry_fit_candidate
  testthat::local_mocked_bindings(.geometry_fit_candidate = function(plan, candidate, component, row_block, upstream = list(), cache = NULL) {
    used <- names(crossform:::.geometry_training_support(plan)$endpoints)
    if (any(used %in% c("g", "h"))) stop("outer-test leakage")
    fitted[[length(fitted) + 1L]] <<- used
    real_fit(plan, candidate, component, row_block, upstream, cache = cache)
  }, .package = "crossform")
  a <- crossform:::.geometry_nested(f$train, f$test, f$candidates, f$inner, scope = "per_measurement")
  b <- crossform:::.geometry_nested(altered$train, altered$test, rev(altered$candidates),
    rev(altered$inner), scope = "per_measurement")
  expect_identical(a$selection$signature, b$selection$signature)
  expect_identical(a$selection$selected, b$selection$selected)
  expect_true(any(a$table$gain != b$table$gain))
  expect_true(all(lengths(fitted) %in% c(4L, 6L)))
  for (id in names(a$fits)) {
    expect_identical(a$fits[[id]]$prediction_id, b$fits[[id]]$prediction_id)
    expect_equal(pg_forms(a$fits[[id]]), pg_forms(b$fits[[id]]))
  }
  bad <- a$selection; bad$selected[] <- names(bad$candidates)[1]
  expect_error(crossform:::.validate_geometry_selection(bad), class = "effect_contract_error")
})

test_that("[T50 T51 T55] infeasible or undeclared selection choices refuse before fitting", {
  f <- pg_nested_fixture()
  testthat::local_mocked_bindings(.geometry_fit_candidate = function(...) stop("fitting began"), .package = "crossform")
  wrong <- f$inner; wrong[[3]]$validate <- c("g", "h")
  expect_identical(catch_refusal(crossform:::.geometry_nested(f$train, f$test, f$candidates, wrong))$reasons[[1]],
    "invalid_nested_partition_split")
  wrong <- f$inner; wrong[[2]]$train <- c("a", "b", "c", "d")
  expect_identical(catch_refusal(crossform:::.geometry_nested(f$train, f$test, f$candidates, wrong))$reasons[[1]],
    "invalid_nested_partition_split")
  wrong <- f$inner; wrong[[1]]$validate <- "e"
  expect_identical(catch_refusal(crossform:::.geometry_nested(f$train, f$test, f$candidates, wrong))$reasons[[1]],
    "inner_split_has_no_declared_products")
  wrong <- f$inner; wrong[[1]]$weight <- 2
  expect_error(crossform:::.geometry_nested(f$train, f$test, f$candidates, wrong), "sum to one")
  expect_error(crossform:::.geometry_nested(f$train, f$test, f$candidates, f$inner,
    measurement_weights = c(first = .8, second = .8)), "summing to one")
  expect_error(crossform:::.geometry_nested(f$train, f$test, f$candidates, f$inner,
    scope = "per_measurement", measurement_weights = c(first = .5, second = .5)), "does not pool")
  expect_identical(catch_refusal(crossform:::.geometry_nested(f$train, f$train, f$candidates, f$inner))$reasons[[1]],
    "overlapping_training_evaluation_origins")
})

test_that("[T52] the scheduled nested generator agrees with the full fit/select/refit/score route", {
  root <- testthat::test_path("..", "..", "benchmarks", "predictive-geometry")
  skip_if_not(file.exists(file.path(root, "selection-calibration.R")), "Source-tree calibration harness is not installed with test files")
  env <- new.env(parent = globalenv())
  for (file in c("calibration.R", "selection-calibration.R", "selection-config.R")) sys.source(file.path(root, file), env)
  withr::local_seed(2026090404)
  RNGkind("L'Ecuyer-CMRG"); set.seed(2026090404)
  for (arm in c("null", "aligned")) {
    design <- env$pgs_design(arm, env$predictive_selection_config)
    parts <- env$pgc_parts(design, .Random.seed)
    expected <- env$pgs_one(design, parts)
    names(parts) <- paste0("run", 1:8)
    domain <- abstract_domain(design$p, id = "selection-calibration")
    rel <- relation(parts, domain = domain, provenance = list(observation_origins = list(
      id = "eight-selection-origins", partitions = as.list(stats::setNames(names(parts), names(parts))),
      independence = "independent", assumption = "Independent Gaussian sampling and fixed model declarations.")))
    at <- compile_frame(whole_brain(normalization = "none"), domain)
    make_plan <- function(edges) plan_geometry(rel, at,
      pairing(names(parts)[edges$left], names(parts)[edges$right], edges$weight,
        independence = "independent", generalizes_over = "run"))
    inner <- lapply(design$inner, function(fold) list(
      train = names(parts)[unique(c(fold$train$left, fold$train$right))],
      validate = names(parts)[unique(c(fold$validate$left, fold$validate$right))], weight = fold$weight))
    candidates <- lapply(design$candidates, function(x) x[setdiff(names(x), "id")])
    actual <- crossform:::.geometry_nested(make_plan(design$train_edges), make_plan(design$test_edges), candidates, inner)
    chosen <- design$candidates[[as.integer(expected[["selected"]])]]
    expect_identical(unname(actual$selection$selected), chosen$id)
    expect_equal(actual$table$gain, expected[["gain"]], tolerance = 1e-10)
    F <- pg_forms(actual$fits[[chosen$id]])[[1]]
    expect_equal(actual$table$gain - (sum(design$truth^2) - sum((design$truth - F)^2)), expected[["error"]], tolerance = 1e-10)
    expect_s3_class(crossform:::.validate_geometry_nested(unserialize(serialize(actual, NULL))), "effect_geometry_nested")
    broken <- actual; broken$table$gain <- broken$table$gain + 1
    expect_error(crossform:::.validate_geometry_nested(broken), class = "effect_contract_error")
  }
})
