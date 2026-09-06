# Transport conditioning is estimand-bearing provenance. These tests follow a
# functionally estimated operator from construction through every population
# inference layer and a serialization round trip.

ptc_fit <- function() {
  labels <- sprintf("s%02d", 1:6)
  sizes <- stats::setNames(7:12, labels)
  gains <- stats::setNames(c(0.8, 1.1, 0.9, 1.3, 1.0, 1.4), labels)
  subjects <- stats::setNames(lapply(labels, function(id) {
    ph_subject(id, sizes[[id]], gains[[id]])
  }), labels)
  transport <- stats::setNames(lapply(labels, function(id) {
    features <- sizes[[id]]
    block <- matrix(0, features, 2L)
    block[cbind(seq_len(features), 1L + (seq_len(features) %% 2L))] <- 1
    location_transport(
      block,
      native_index = paste0(id, "-f", seq_len(features)),
      group_index = c("group1", "group2"),
      semantics = "budget",
      provenance = list(
        method = "functional",
        details = "held-out response alignment v2",
        fitting_sample = c("session-A", "session-B"),
        cross_fit = c("fold-train-A", "fold-train-B"),
        cross_fit_folds = c("fold-A", "fold-B")
      )
    )
  }), labels)
  plan <- plan_population(subjects, transport)
  list(
    plan = plan,
    fit = estimate_population(plan, rbind(`face-house` = c(1, -1)))
  )
}

test_that("functional transport declares conditional inference completely", {
  carrier <- ptc_fit()$plan$transport[[1L]]
  conditioning <- carrier$provenance$conditioning

  expect_identical(conditioning$source, "held-out response alignment v2")
  expect_identical(conditioning$operator_status, "estimated")
  expect_identical(conditioning$fitting_sample, c("session-A", "session-B"))
  expect_identical(conditioning$cross_fit_folds, c("fold-A", "fold-B"))
  expect_identical(conditioning$circularity_control,
    "cross_fit_partitions_declared")
  expect_identical(conditioning$inference_scope,
    "conditional_on_realized_transport")
  expect_false(conditioning$uncertainty_propagated)
  expect_false(conditioning$marginal_over_transport)
  expect_true(all(c("transport_operator_estimation",
    "cross_fit_fold_assignment") %in% conditioning$excluded_uncertainty))
  expect_identical(conditioning$future$capability,
    "transport_uncertainty_propagation")
  expect_identical(conditioning$future$status, "not_implemented")
})

test_that("cross-fitting cannot be relabelled as marginal inference", {
  block <- diag(2)
  provenance <- list(
    method = "functional", details = "response alignment",
    cross_fit = c("train", "test"), marginal_over_transport = TRUE
  )
  expect_error(
    location_transport(block, c("n1", "n2"), c("g1", "g2"),
      semantics = "budget", provenance = provenance),
    "does not make inference marginal over transport estimation",
    class = "effect_input_error"
  )
})

test_that("conditioning survives plan execution serialization and views", {
  court <- ptc_fit()
  plan <- court$plan
  fit <- court$fit
  labels <- names(plan$subjects)

  expect_true(all(plan$subject_index$transport_status == "estimated"))
  expect_true(all(plan$subject_index$fitting_sample ==
    "session-A; session-B"))
  expect_true(all(plan$subject_index$cross_fit_folds == "fold-A; fold-B"))
  expect_false(any(plan$subject_index$transport_uncertainty_propagated))
  expect_false(any(plan$subject_index$marginal_over_transport))

  conditioning <- fit$coverage$conditioning
  expect_identical(names(conditioning$transport_by_subject), labels)
  expect_false(conditioning$uncertainty_propagated)
  expect_false(conditioning$marginal_over_transport)
  expect_true("transport_operator_estimation" %in%
    conditioning$excluded_uncertainty)

  robust <- population_uncertainty(fit, estimator = "HC3")
  expect_true("heteroskedasticity_robust_sandwich" %in%
    robust$between$assumptions)
  expect_identical(robust$between$conditioning, conditioning)
  # Concrete boundary: HC3 includes heteroskedastic subject residuals, while
  # the learned alignment and its fold assignment remain excluded.
  expect_true(all(c("transport_operator_estimation",
    "cross_fit_fold_assignment") %in%
      robust$between$conditioning$excluded_uncertainty))

  wild <- population_wild_bootstrap(
    fit, "(Intercept)", replicates = 99, seed = 812
  )
  expect_identical(wild$conditioning, conditioning)

  view <- contrast_energy(fit, c(1, -1))
  expect_identical(view$coverage$conditioning, conditioning)
  expect_identical(view$receipt$transport$conditioning, conditioning)

  copies <- lapply(list(plan, fit, robust, wild, view), function(x) {
    unserialize(serialize(x, NULL))
  })
  .validate_population_plan(copies[[1L]], deep = FALSE)
  .validate_population_result(copies[[2L]])
  .validate_population_uncertainty(copies[[3L]])
  .validate_population_bootstrap(copies[[4L]])
  .validate_population_view(copies[[5L]])
  expect_identical(copies[[2L]]$coverage$conditioning, conditioning)
  expect_identical(copies[[5L]]$receipt$transport$conditioning, conditioning)

  withr::local_options(list(width = 200))
  printed <- paste(utils::capture.output(print(view)), collapse = " ")
  expect_match(printed, "conditional_on_realized_transport", fixed = TRUE)
  expect_match(printed, "uncertainty not propa", fixed = TRUE)
  robust_print <- gsub("\\s+", " ", paste(
    utils::capture.output(print(robust)), collapse = " "
  ))
  expect_match(robust_print,
    "Cross-fitting limits circularity", fixed = TRUE)
})

test_that("fixed transports declare the distinct conditional boundary", {
  carrier <- anatomical_transport(
    cbind(c(0, 1)), cbind(c(0, 1)), semantics = "budget"
  )
  conditioning <- carrier$provenance$conditioning
  expect_identical(conditioning$operator_status, "fixed")
  expect_identical(conditioning$fitting_sample, character())
  expect_identical(conditioning$cross_fit_folds, character())
  expect_identical(conditioning$circularity_control, "fixed_operator")
  expect_identical(conditioning$excluded_uncertainty,
    "transport_source_choice")
})
