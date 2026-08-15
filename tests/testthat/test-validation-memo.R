# The validation memo is an execution optimization: it re-admits an object only
# when that object is structurally identical to one already validated at the
# same depth or deeper. These tests hold it to that claim -- identical results
# whether cold or warm, and no weakening of any identity check.

memo_sweep_fixture <- function() {
  set.seed(4113L)
  conditions <- 4L
  partitions <- 4L
  observations <- 32L
  domain <- volume_domain(
    array(TRUE, c(4L, 4L, 3L)), spacing = c(3, 3, 3), id = "memo-sweep"
  )
  features <- domain$n_features
  frame <- compile_frame(searchlights(4), domain)
  condition <- factor(rep(seq_len(conditions), length.out = observations))
  design <- stats::model.matrix(~ 0 + condition)
  colnames(design) <- paste0("condition", seq_len(conditions))
  effects <- diag(conditions)
  dimnames(effects) <- list(colnames(design), colnames(design))
  signal <- matrix(rnorm(conditions * features, sd = 0.3), conditions, features)
  sources <- stats::setNames(lapply(seq_len(partitions), function(run) {
    design %*% signal + matrix(rnorm(observations * features), observations,
      features)
  }), paste0("run", seq_len(partitions)))
  fit <- lm_relation_fit(
    sources, design, effects, sampling_unit = "trial", domain = domain
  )
  metric <- noise_precision(
    diag(features), domain, covariance = diag(features),
    provenance = list(source = "memo-sweep-fixed-metric")
  )
  plan <- plan_geometry(
    fit$relation, at = frame, over = cross_partitions(fit$relation),
    metric = metric
  )
  list(plan = plan, fit = fit, nodes = nrow(frame$weights))
}

memo_node_result <- function(fixture, node, target) {
  covariance <- rdm_sampling_covariance(
    fixture$plan, fixture$fit, target = target, at = node
  )
  dimension <- covariance$dimension
  weights <- seq_len(dimension) / dimension
  map <- matrix(1 / seq_len(2L * dimension), 2L, dimension)
  list(
    diagonal = sampling_covariance(covariance),
    applied = sampling_covariance(covariance, "apply", weights),
    quadratic = sampling_covariance(covariance, "quadratic_form", weights),
    transport = sampling_covariance(covariance, "transport", map),
    dense = sampling_covariance(covariance, "materialize"),
    signal_factor = covariance$signal_factor,
    xi_factor = covariance$xi_factor,
    noise_trace = covariance$noise_trace,
    signature = covariance$signature,
    plan_id = covariance$plan$scientific_plan_id,
    plan_signature = covariance$plan$signature
  )
}

# Reproduces the validation behaviour this package had before the sampling
# sweep was optimized: every memo lookup misses, and a shallow geometry-plan
# check still rebuilds its neural metric in full. Results computed under it are
# the pre-optimization oracle.
with_unoptimized_validation <- function(code) {
  validate_metric <- effectagram:::.validate_neural_metric
  testthat::with_mocked_bindings(
    code,
    .validated_before = function(object, key, deep = TRUE) FALSE,
    .validate_neural_metric = function(x, deep = TRUE) {
      validate_metric(x, deep = TRUE)
    },
    .package = "effectagram"
  )
}

test_that("a node sweep is bit-identical to the pre-optimization path", {
  fixture <- memo_sweep_fixture()
  nodes <- unique(round(seq(1, fixture$nodes, length.out = 6L)))
  sweep <- function() {
    lapply(nodes, function(node) {
      lapply(c("plugin", "null"), function(target) {
        memo_node_result(fixture, node, target)
      })
    })
  }
  oracle <- with_unoptimized_validation(sweep())
  expect_identical(sweep(), oracle)
})

test_that("a warm memo leaves the point estimate and plan identity unchanged", {
  fixture <- memo_sweep_fixture()
  effectagram:::.reset_validation_memo()
  cold_plan_id <- fixture$plan$scientific_plan_id
  cold_rdm <- as.matrix(rdm(fixture$plan)$values)
  invisible(rdm_sampling_covariance(
    fixture$plan, fixture$fit, target = "plugin", at = 1L
  ))
  expect_identical(fixture$plan$scientific_plan_id, cold_plan_id)
  expect_identical(as.matrix(rdm(fixture$plan)$values), cold_rdm)
})

test_that("a fixed-metric geometry view matches the pre-optimization path", {
  fixture <- memo_sweep_fixture()
  bare <- plan_geometry(
    fixture$fit$relation, at = fixture$plan$frame,
    over = cross_partitions(fixture$fit$relation)
  )
  view <- function(plan) {
    value <- rdm(plan)
    list(values = as.matrix(value$values), index = value$index,
      plan_id = plan$scientific_plan_id, signature = plan$signature)
  }
  oracle <- with_unoptimized_validation(list(view(fixture$plan), view(bare)))
  expect_identical(list(view(fixture$plan), view(bare)), oracle)
  # A domain-wide identity metric must not move the point estimate either.
  expect_equal(oracle[[1L]]$values, oracle[[2L]]$values, tolerance = 1e-12)
})

test_that("the memo never re-admits a tampered neural metric", {
  fixture <- memo_sweep_fixture()
  metric <- fixture$plan$metric_schedule$metric
  expect_silent(effectagram:::.validate_neural_metric(metric))
  tampered <- metric
  tampered$value[1L, 1L] <- tampered$value[1L, 1L] + 1
  expect_error(
    effectagram:::.validate_neural_metric(tampered), "inconsistent"
  )
})

test_that("the memo never re-admits a tampered object", {
  fixture <- memo_sweep_fixture()
  covariance <- rdm_sampling_covariance(
    fixture$plan, fixture$fit, target = "plugin", at = 1L
  )
  expect_silent(effectagram:::.validate_sampling_covariance(covariance))
  tampered <- covariance
  tampered$noise_trace <- covariance$noise_trace * 2
  expect_error(
    effectagram:::.validate_sampling_covariance(tampered),
    "identity is inconsistent"
  )
  tampered_plan <- covariance$plan
  tampered_plan$spatial_scope <- "modeled"
  expect_error(
    effectagram:::.validate_evidence_sampling_plan(tampered_plan),
    "inconsistent"
  )
})

test_that("a shallow validation does not satisfy a later deep request", {
  fixture <- memo_sweep_fixture()
  forged <- fixture$plan
  forged$scientific_plan_id <- paste0("geometry-sha256:", strrep("a", 64L))
  effectagram:::.reset_validation_memo()
  expect_silent(effectagram:::.validate_geometry_plan(forged, deep = FALSE))
  expect_error(
    effectagram:::.validate_geometry_plan(forged, deep = TRUE),
    "identity is inconsistent"
  )
})

test_that("shallow geometry-plan validation still binds the metric schedule", {
  fixture <- memo_sweep_fixture()
  broken <- fixture$plan
  broken$metric_schedule$metric_signature <- paste0(
    "sha256:", strrep("b", 64L)
  )
  effectagram:::.reset_validation_memo()
  expect_error(
    effectagram:::.validate_geometry_plan(broken, deep = FALSE),
    "inconsistent"
  )
})

test_that("the memo store is bounded by its fixed key set", {
  fixture <- memo_sweep_fixture()
  effectagram:::.reset_validation_memo()
  for (node in seq_len(min(6L, fixture$nodes))) {
    invisible(rdm_sampling_covariance(
      fixture$plan, fixture$fit, target = "plugin", at = node
    ))
  }
  keys <- ls(effectagram:::.validation_memo, all.names = TRUE)
  expect_true(length(keys) <= 8L)
  expect_true(all(keys %in% c(
    "evidence_sampling_plan", "geometry_plan", "relation_fit",
    "sampling_covariance", "neural_metric"
  )))
  effectagram:::.reset_validation_memo()
  expect_identical(ls(effectagram:::.validation_memo, all.names = TRUE),
    character())
})
