# The validation memo is an execution optimization: it re-admits an object only
# when that object is structurally identical to one already validated at the
# same depth or deeper. These tests hold it to that claim -- identical results
# whether cold or warm, and no weakening of any identity check.

# The fixture is a pure function of a fixed seed, so every block in this file
# asked for byte-identical objects and paid to rebuild them. Building it once
# and handing out the same immutable value is not a shared-state hazard: R's
# copy-on-modify semantics mean the blocks below that forge or tamper with a
# plan still get their own copy, and the memo compares with `identical()`, so a
# reused object is admitted on exactly the same terms as a freshly built one.
#
# Its dimensions are the smallest that keep every structure these tests read
# non-degenerate: eighteen features carrying eighteen distinct overlapping
# searchlight supports of four to six voxels, four conditions so the sampling
# covariance keeps all six condition pairs, and two partitions so the
# cross-partition schedule is a genuine independent fold. Nothing asserted
# below is a claim about scale, so the sizes only have to stay above
# degeneracy, not reproduce a whole-brain sweep.
memo_sweep_fixture <- local({
  built <- NULL
  function() {
    if (!is.null(built)) {
      return(built)
    }
    set.seed(4113L)
    conditions <- 4L
    partitions <- 2L
    observations <- 8L
    domain <- volume_domain(
      array(TRUE, c(3L, 3L, 2L)), spacing = c(3, 3, 3), id = "memo-sweep"
    )
    features <- domain$n_features
    frame <- compile_frame(searchlights(4), domain)
    condition <- factor(rep(seq_len(conditions), length.out = observations))
    design <- stats::model.matrix(~ 0 + condition)
    colnames(design) <- paste0("condition", seq_len(conditions))
    effects <- diag(conditions)
    dimnames(effects) <- list(colnames(design), colnames(design))
    signal <- matrix(rnorm(conditions * features, sd = 0.3), conditions,
      features)
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
      fit$relation, at = frame,
      over = cross_partitions(fit$relation, independence = "independent"),
      metric = metric
    )
    built <<- list(plan = plan, fit = fit, nodes = nrow(frame$weights))
    built
  }
})

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
  validate_metric <- crossform:::.validate_neural_metric
  testthat::with_mocked_bindings(
    code,
    .validated_before = function(object, key, deep = TRUE) FALSE,
    .validate_neural_metric = function(x, deep = TRUE) {
      validate_metric(x, deep = TRUE)
    },
    .package = "crossform"
  )
}

# The oracle below re-validates every object in the plan graph on every call,
# and that costs about a second of CPU per (node, target) pair however small
# the fixture is -- the work is a fixed number of identifier and signature
# checks, not anything proportional to the data. So the sweep is kept to the
# smallest schedule that still covers what the claim is about: the first pair
# fills the memo cold, and the second is answered warm at a *different* node
# and a *different* target. Every reading of both pairs -- all five covariance
# operations, the factorizations, the noise trace, and both plan identities --
# is still compared field by field against the pre-optimization value.
test_that("a node sweep is bit-identical to the pre-optimization path", {
  fixture <- memo_sweep_fixture()
  nodes <- unique(round(seq(1, fixture$nodes, length.out = 2L)))
  targets <- c("plugin", "null")
  sweep <- function() {
    Map(function(node, target) {
      memo_node_result(fixture, node, target)
    }, nodes, targets)
  }
  oracle <- with_unoptimized_validation(sweep())
  expect_identical(sweep(), oracle)
})

test_that("a warm memo leaves the point estimate and plan identity unchanged", {
  fixture <- memo_sweep_fixture()
  crossform:::.reset_validation_memo()
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
    over = cross_partitions(
      fixture$fit$relation, independence = "independent"
    )
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
  expect_silent(crossform:::.validate_neural_metric(metric))
  tampered <- metric
  tampered$value[1L, 1L] <- tampered$value[1L, 1L] + 1
  expect_error(
    crossform:::.validate_neural_metric(tampered), "inconsistent"
  )
})

test_that("the memo never re-admits a tampered object", {
  fixture <- memo_sweep_fixture()
  covariance <- rdm_sampling_covariance(
    fixture$plan, fixture$fit, target = "plugin", at = 1L
  )
  expect_silent(crossform:::.validate_sampling_covariance(covariance))
  tampered <- covariance
  tampered$noise_trace <- covariance$noise_trace * 2
  expect_error(
    crossform:::.validate_sampling_covariance(tampered),
    "identity is inconsistent"
  )
  tampered_plan <- covariance$plan
  tampered_plan$spatial_scope <- "modeled"
  expect_error(
    crossform:::.validate_evidence_sampling_plan(tampered_plan),
    "inconsistent"
  )
})

test_that("a shallow validation does not satisfy a later deep request", {
  fixture <- memo_sweep_fixture()
  forged <- fixture$plan
  forged$scientific_plan_id <- paste0("geometry-sha256:", strrep("a", 64L))
  crossform:::.reset_validation_memo()
  expect_silent(crossform:::.validate_geometry_plan(forged, deep = FALSE))
  expect_error(
    crossform:::.validate_geometry_plan(forged, deep = TRUE),
    "identity is inconsistent"
  )
})

test_that("shallow geometry-plan validation still binds the metric schedule", {
  fixture <- memo_sweep_fixture()
  broken <- fixture$plan
  broken$metric_schedule$metric_signature <- paste0(
    "sha256:", strrep("b", 64L)
  )
  crossform:::.reset_validation_memo()
  expect_error(
    crossform:::.validate_geometry_plan(broken, deep = FALSE),
    "inconsistent"
  )
})

test_that("the memo store is bounded by its fixed key set", {
  fixture <- memo_sweep_fixture()
  crossform:::.reset_validation_memo()
  for (node in seq_len(min(6L, fixture$nodes))) {
    invisible(rdm_sampling_covariance(
      fixture$plan, fixture$fit, target = "plugin", at = node
    ))
  }
  keys <- ls(crossform:::.validation_memo, all.names = TRUE)
  expect_true(length(keys) <= 8L)
  expect_true(all(keys %in% c(
    "evidence_sampling_plan", "geometry_plan", "relation_fit",
    "sampling_covariance", "neural_metric"
  )))
  crossform:::.reset_validation_memo()
  expect_identical(ls(crossform:::.validation_memo, all.names = TRUE),
    character())
})
