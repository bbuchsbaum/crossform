pg_preflight <- function(fit, plan, ...) crossform:::.geometry_score_prepare(fit, plan, ...)
pg_reason <- function(expr) catch_refusal(expr)$reasons[[1L]]

test_that("[T28 T31 T40] disjoint pairings match targets by stable measurement identity", {
  f <- pg_fixture()
  fit <- fit_geometry(f$plan, f$basis, rank = 2, penalty = .2)
  testthat::local_mocked_bindings(.geometry_fit_materialize = function(...) stop("neural read"), .package = "crossform")
  admitted <- pg_preflight(fit, pg_plan(f))
  expect_identical(admitted$mapping, 1:2)
  expect_true(admitted$validation$independent)
  expect_identical(admitted$prepared$support$observations, c("raw-c", "raw-d"))
  reversed <- additive_frame(members = list(1L, 1:2), domain = f$domain,
    measurements = c("second", "first"))
  expect_identical(pg_preflight(fit, pg_plan(f, at = reversed))$mapping, 2:1)
  # Different source bytes are not a requirement for independent origins.
  parts <- stats::setNames(rep(list(f$B), 4), letters[1:4])
  equal <- relation(parts, domain = f$domain, provenance = f$rel$provenance)
  expect_true(pg_preflight(fit, pg_plan(f, rel = equal))$validation$independent)
})

test_that("[T27 T28 T30] aliases, shared products and unknown ancestry refuse before execution", {
  f <- pg_fixture(); fit <- fit_geometry(f$plan, f$basis, rank = 1, penalty = .2)
  testthat::local_mocked_bindings(.geometry_fit_materialize = function(...) stop("neural read"), .package = "crossform")
  expect_identical(pg_reason(pg_preflight(fit, pg_plan(f, left = "a"))),
    "overlapping_training_evaluation_origins")
  copied <- crossform:::.relation_subset(f$rel, c("a", "b"), c("copy-c", "copy-d"))
  expect_identical(pg_reason(pg_preflight(fit, pg_plan(f, "copy-c", "copy-d", rel = copied))),
    "overlapping_training_evaluation_origins")
  expect_identical(pg_reason(pg_preflight(fit, pg_plan(f, independence = NULL))), "independence_undeclared")
  unknown <- pg_fixture(origins = FALSE)
  expect_identical(pg_reason(pg_preflight(fit, pg_plan(unknown))), "unknown_observation_origins")
  # Two evaluation endpoints can themselves have overlapping raw ancestry.
  raw_manifest <- list(id = "within-test", partitions = list(a = "a", b = "b", c = "c", d = "c"),
    independence = "independent", assumption = "Declared independent runs.")
  rel <- relation(stats::setNames(rep(list(f$B), 4), letters[1:4]), domain = f$domain,
    provenance = list(observation_origins = raw_manifest))
  train <- pg_plan(f, "a", "b", rel = rel)
  # Mock only materialization after obtaining the new fit in another scope.
  support <- crossform:::.geometry_training_support(train)
  altered <- fit; altered$training <- support
  altered$prediction_id <- crossform:::.geometry_fit_prediction_id(altered)
  altered$signature <- crossform:::.sha256_signature(crossform:::.geometry_fit_semantic(altered), "geometry-fit-sha256:")
  expect_identical(pg_reason(pg_preflight(altered, pg_plan(f, rel = rel))), "overlapping_neural_origins")
})

test_that("[T31 T40 T55] same dimensions do not admit different predictive targets", {
  f <- pg_fixture(); fit <- fit_geometry(f$plan, f$basis, rank = 1, penalty = .2)
  testthat::local_mocked_bindings(.geometry_fit_materialize = function(...) stop("neural read"), .package = "crossform")
  expect_identical(pg_reason(pg_preflight(fit, pg_plan(f), component = "configuration")),
    "prediction_component_mismatch")
  expect_identical(pg_reason(pg_preflight(fit, pg_plan(f, generalizes_over = "session"))),
    "prediction_generalizes_over_mismatch")
  swapped <- additive_frame(members = list(1L, 1:2), domain = f$domain, measurements = c("first", "second"))
  expect_identical(pg_reason(pg_preflight(fit, pg_plan(f, at = swapped))), "prediction_measurement_support_mismatch")
  missing <- additive_frame(members = list(1:2), domain = f$domain, measurements = "first")
  expect_identical(pg_reason(pg_preflight(fit, pg_plan(f, at = missing))), "prediction_measurements_mismatch")
  metric <- neural_metric(diag(c(2, 1)), domain = f$domain)
  expect_identical(pg_reason(pg_preflight(fit, pg_plan(f, metric = metric))), "prediction_metric_signature_mismatch")
  # Corrupt metadata is a contract error, not an opportunity to reindex by position.
  corrupt <- fit; corrupt$target$measurement_ids <- c("first", "first")
  expect_error(pg_preflight(corrupt, pg_plan(f)), class = "effect_contract_error")
  target <- fit$target; target$effect_space$units[] <- "other"
  expect_error(crossform:::.geometry_prediction_match(fit$target, target), class = "effect_contract_error")
})

test_that("[T29 T32] evaluation revisions change evidence provenance, never fitting dependencies", {
  f <- pg_fixture(); fit <- fit_geometry(f$plan, f$basis, rank = 1, penalty = .2)
  parts <- list(a = f$B, b = f$B, c = -100 * f$B, d = 3 * f$B)
  changed <- relation(parts, domain = f$domain, provenance = f$rel$provenance)
  refit <- fit_geometry(pg_plan(f, "a", "b", rel = changed), f$basis, rank = 1, penalty = .2)
  expect_identical(fit$prediction_id, refit$prediction_id)
  expect_equal(pg_forms(fit), pg_forms(refit))
  a <- pg_preflight(fit, pg_plan(f)); b <- pg_preflight(fit, pg_plan(f, rel = changed))
  expect_identical(a$validation$training, b$validation$training)
  expect_false(identical(a$validation$evaluation, b$validation$evaluation))
  changed_fit <- fit; changed_fit$training$endpoints$a$source_revision <- changed$sources$c$stable_revision
  expect_error(pg_preflight(changed_fit, pg_plan(f)), class = "effect_contract_error")
})
