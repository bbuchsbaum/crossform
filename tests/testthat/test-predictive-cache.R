test_that("[T18 T53] candidate ranks share neural forms and shifted spectra within one training context", {
  f <- pg_nested_fixture()
  candidates <- crossform:::.geometry_candidates(f$candidates)
  cache <- crossform:::.new_geometry_training_cache(candidates)
  evaluate_cache <- crossform:::.new_geometry_training_cache(candidates, "evaluation")
  reads <- 0L; spectra <- 0L
  materialize <- crossform:::.geometry_fit_materialize; spectral <- crossform:::.geometry_fit_path
  testthat::local_mocked_bindings(.geometry_fit_materialize = function(...) {
    reads <<- reads + 1L; materialize(...)
  }, .geometry_fit_path = function(...) {
    spectra <<- spectra + 1L; spectral(...)
  }, .package = "crossform")
  fits <- lapply(candidates, function(candidate) crossform:::.geometry_fit_candidate(
    f$train, candidate, "total", 1L, cache = cache))
  expected_paths <- length(unique(vapply(candidates, `[[`, 0, "penalty"))) * f$train$measurements
  expect_identical(reads, 1L)
  expect_equal(spectra, expected_paths)
  stats <- crossform:::.geometry_cache_summary(cache)
  expect_equal(stats$counters[["form_misses"]], 1L)
  expect_equal(stats$counters[["form_hits"]], length(candidates) - 1L)
  expect_equal(stats$counters[["pool_misses"]], 1L)
  expect_equal(stats$counters[["path_hits"]], 2L * length(candidates) - expected_paths)
  scores <- lapply(fits, function(fit) crossform:::.geometry_score_candidate(fit, f$test, 1L, evaluate_cache))
  expect_identical(reads, 2L)
  expect_equal(evaluate_cache$counters[["form_hits"]], length(candidates) - 1L)
  expect_identical(evaluate_cache$counters[["path_misses"]], 0L)
  expect_error(crossform:::.geometry_fit_candidate(f$train, candidates[[1]], "total", 1L,
    cache = evaluate_cache), "matching training")
  # Cached and uncached physical paths produce identical predictions/gains.
  for (i in seq_along(candidates)) {
    reference <- crossform:::.geometry_fit_candidate(f$train, candidates[[i]], "total", 1L)
    expect_identical(fits[[i]]$prediction_id, reference$prediction_id)
    expect_equal(pg_forms(fits[[i]]), pg_forms(reference), tolerance = 1e-11)
    expect_equal(scores[[i]]$table, score_geometry(reference, f$test)$table, tolerance = 1e-11)
  }
})

test_that("[T29 T32 T53] actual source/fold/target/execution changes invalidate neural cache entries", {
  f <- pg_nested_fixture(); changed_test <- pg_nested_fixture(test_scale = -8)
  candidates <- crossform:::.geometry_candidates(f$candidates)
  candidate <- candidates[[which(vapply(candidates, `[[`, 0L, "rank") == 2L)[1]]]
  cache <- crossform:::.new_geometry_training_cache(candidates)
  run <- function(plan, component = "total", row_block = 1L) crossform:::.geometry_fit_candidate(
    plan, candidate, component, row_block, cache = cache)
  original <- run(f$train)
  same_training <- run(changed_test$train)
  expect_identical(original$prediction_id, same_training$prediction_id)
  expect_identical(cache$counters[["form_misses"]], 1L)
  expect_identical(cache$counters[["form_hits"]], 1L)
  parts <- lapply(f$rel$sources, function(source) source$read(1:2)); parts$a <- 1.1 * parts$a
  rel <- relation(parts, domain = f$domain, provenance = f$rel$provenance)
  source_changed <- plan_geometry(rel, f$at, f$train$pairing)
  shifted <- run(source_changed)
  expect_false(identical(original$prediction_id, shifted$prediction_id))
  expect_identical(cache$counters[["form_misses"]], 2L)
  subset <- crossform:::.geometry_prediction_pair_plan(f$train, as.data.frame(f$train$pairing)[1:3, ])
  reversed <- additive_frame(members = list(1L, 1:2), domain = f$domain, measurements = c("second", "first"))
  metric <- neural_metric(matrix(c(2, .3, .3, 1), 2), f$domain)
  variants <- list(subset,
    plan_geometry(f$rel, reversed, f$train$pairing),
    plan_geometry(f$rel, f$at, f$train$pairing, metric = metric),
    plan_geometry(f$rel, f$at, f$train$pairing, compute = compute_policy(block_features = 1)))
  for (plan in variants) {
    before <- cache$counters[["form_misses"]]
    actual <- run(plan)
    expected <- crossform:::.geometry_fit_candidate(plan, candidate, "total", 1L)
    expect_identical(cache$counters[["form_misses"]], before + 1L)
    expect_equal(pg_forms(actual), pg_forms(expected), tolerance = 1e-11)
    expect_identical(length(cache$forms), 1L)
  }
  coherent <- run(f$train, "coherent")
  expect_equal(pg_forms(coherent), pg_forms(crossform:::.geometry_fit_candidate(f$train, candidate, "coherent", 1L)))
  before <- cache$counters[["context_resets"]]
  run(f$train, "coherent", row_block = 2L)
  expect_identical(cache$counters[["context_resets"]], before + 1L)
  expect_identical(length(cache$forms), 1L)
})

test_that("[T09 T53] pool support, model ranks and normalization are distinct reusable preparations", {
  f <- pg_fixture()
  kernels <- list(a = tcrossprod(f$V[, 1]), b = tcrossprod(f$V[, 2]))
  basis <- model_basis(kernels = kernels, conditions = f$rel$effect_space, normalize = "trace")
  tuples <- lapply(list(c(a = 1, b = 0), c(a = .5, b = .5), c(a = 0, b = 1)), function(w)
    list(basis = basis, weights = w, rank = 1L, penalty = .2))
  # Same raw target, different model truncation/normalization declarations.
  tuples <- c(tuples, list(
    list(basis = model_basis(kernels = list(m = 10 * f$K), conditions = f$rel$effect_space,
      rank = 1, normalize = "trace"), rank = 1L, penalty = .2),
    list(basis = model_basis(kernels = list(m = 10 * f$K), conditions = f$rel$effect_space,
      normalize = "none"), rank = 1L, penalty = .2)))
  candidates <- crossform:::.geometry_candidates(tuples)
  cache <- crossform:::.new_geometry_training_cache(candidates)
  fits <- lapply(candidates, function(candidate) crossform:::.geometry_fit_candidate(f$plan,
    candidate, "total", 1L, cache = cache))
  expect_identical(cache$counters[["form_misses"]], 3L)
  expect_identical(cache$counters[["pool_misses"]], 5L)
  expect_identical(cache$counters[["path_misses"]], 10L)
  for (i in seq_along(candidates)) {
    reference <- crossform:::.geometry_fit_candidate(f$plan, candidates[[i]], "total", 1L)
    expect_equal(pg_forms(fits[[i]]), pg_forms(reference), tolerance = 1e-11)
  }
  unknown <- candidates[[1]]; unknown$penalty <- .333
  testthat::local_mocked_bindings(.geometry_fit_materialize = function(...) stop("undeclared recipe read"), .package = "crossform")
  cache <- crossform:::.new_geometry_training_cache(candidates)
  expect_error(crossform:::.geometry_fit_candidate(f$plan, unknown, "total", 1L, cache = cache), "not declared")
  expect_identical(cache$counters[["form_misses"]], 0L)
})

test_that("[T29 T31 T32 T53] effect meaning, extractor maps and upstream support invalidate cached forms", {
  f <- pg_fixture()
  spaces <- list(f$rel$effect_space,
    effect_space(f$rel$effects, basis_id = "alternate-estimand", units = "percent", scale = 2))
  candidates <- crossform:::.geometry_candidates(lapply(spaces, function(space)
    list(basis = model_basis(kernels = list(m = f$K), conditions = space, normalize = "trace"),
      rank = 1L, penalty = .2)))
  cache <- crossform:::.new_geometry_training_cache(candidates)
  parts <- lapply(f$rel$sources, function(source) source$read(1:2))
  for (candidate in candidates) {
    space <- spaces[[which(vapply(spaces, function(space)
      identical(space$signature, candidate$basis$centering$declaration$signature), FALSE))]]
    for (scale in c(1, 1.1)) {
      ext <- effect_extractor(diag(scale, 4), effects = space)
      rel <- relation(parts, extract = ext, domain = f$domain, provenance = f$rel$provenance)
      plan <- plan_geometry(rel, f$at, f$plan$pairing)
      before <- cache$counters[["form_misses"]]
      fit <- crossform:::.geometry_fit_candidate(plan, candidate, "total", 1L, cache = cache)
      reference <- crossform:::.geometry_fit_candidate(plan, candidate, "total", 1L)
      expect_identical(cache$counters[["form_misses"]], before + 1L)
      expect_equal(pg_forms(fit), pg_forms(reference), tolerance = 1e-11)
    }
    upstream <- list(crossform:::.geometry_training_support(pg_plan(f)))
    before <- cache$counters[["form_misses"]]
    fit <- crossform:::.geometry_fit_candidate(plan, candidate, "total", 1L, upstream, cache)
    expect_identical(cache$counters[["form_misses"]], before + 1L)
    expect_equal(pg_forms(fit), pg_forms(reference), tolerance = 1e-11)
    expect_false(identical(fit$prediction_id, reference$prediction_id))
  }
})

test_that("[T32 T49 T53 T57] nested reuse matches uncached choices and its cache memory is admitted", {
  f <- pg_nested_fixture()
  # One declared inner split suffices for this cache witness; the full three
  # fold oracle is in test-predictive-selection.R.
  inner <- f$inner[1]; inner[[1]]$weight <- 1
  cached <- crossform:::.geometry_nested(f$train, f$test, f$candidates, inner, reuse = TRUE)
  uncached <- crossform:::.geometry_nested(f$train, f$test, f$candidates, inner, reuse = FALSE)
  expect_identical(cached$selection$signature, uncached$selection$signature)
  expect_equal(cached$table, uncached$table, tolerance = 1e-11)
  expect_identical(cached$cache_statistics[[1]]$training$counters[["form_misses"]], 1L)
  expect_identical(cached$cache_statistics[[1]]$evaluation$counters[["form_misses"]], 1L)
  expect_true(all(vapply(cached$fits, function(fit) fit$workspace$categories[["training_cache"]] > 0, FALSE)))
  candidates <- crossform:::.geometry_candidates(f$candidates)
  cache <- crossform:::.new_geometry_training_cache(candidates)
  limited <- plan_geometry(f$rel, f$at, f$train$pairing, compute = compute_policy(workspace_bytes = 1))
  testthat::local_mocked_bindings(.geometry_fit_materialize = function(...) stop("unadmitted cache read"), .package = "crossform")
  expect_identical(catch_refusal(crossform:::.geometry_fit_candidate(limited, candidates[[1]], "total", 1L,
    cache = cache))$reasons[[1]], "predictive_workspace_budget_exceeded")
  expect_identical(cache$counters[["form_misses"]], 0L)
  expect_identical(cache$counters[["path_misses"]], 0L)
})
