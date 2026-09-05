test_that("[T14] public fitting compresses the signed full form before any PSD projection", {
  f <- pg_fixture(3L)
  G <- matrix(c(1, 2, 2, 1), 2)
  rel <- relation(list(a = f$V, b = f$V %*% G), domain = f$domain)
  at <- compile_frame(whole_brain(normalization = "none"), f$domain)
  plan <- plan_geometry(rel, at, pairing("a", "b", independence = "independent", generalizes_over = "run"))
  basis <- model_basis(features = f$V[, 1, drop = FALSE], conditions = rel$effect_space, normalize = "trace")
  fit <- fit_geometry(plan, basis, rank = 1, penalty = 0, retain_source = TRUE)
  row <- crossform:::.geometry_prediction_rows(fit, 1L)[[1L]]
  expect_equal(row$amplitudes, 1, tolerance = 1e-12)
  expect_equal(fit$diagnostics$prediction_norm_sq, 1, tolerance = 1e-12)
  expect_equal(as.numeric(geometry_component(fit$source, "total")), 1, tolerance = 1e-12)
  # The same full-form positive projection would instead yield amplitude 1.5.
  e <- eigen(G, symmetric = TRUE)
  wrong <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors)
  expect_gt(abs(row$amplitudes - wrong[1, 1]), .49)
})

test_that("[T53] a changed training edge measure invalidates reuse with unchanged endpoints and shapes", {
  f <- pg_fixture()
  candidates <- crossform:::.geometry_candidates(list(list(basis = f$basis, rank = 1L, penalty = .2)))
  cache <- crossform:::.new_geometry_training_cache(candidates)
  plan <- function(weights) plan_geometry(f$rel, f$at,
    pairing(c("a", "a", "b"), c("b", "c", "c"), weight = weights,
      independence = "independent", generalizes_over = "run"))
  first <- crossform:::.geometry_fit_candidate(plan(c(1, 2, 3)), candidates[[1]], "total", 1L, cache = cache)
  second_plan <- plan(c(3, 2, 1))
  second <- crossform:::.geometry_fit_candidate(second_plan, candidates[[1]], "total", 1L, cache = cache)
  oracle <- crossform:::.geometry_fit_candidate(second_plan, candidates[[1]], "total", 1L)
  expect_identical(cache$counters[["form_misses"]], 2L)
  expect_gt(max(abs(pg_forms(first)[[1]] - pg_forms(oracle)[[1]])), .01)
  expect_equal(pg_forms(second), pg_forms(oracle), tolerance = 1e-11)
})
