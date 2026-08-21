# Population coverage is part of the estimand, not a display statistic. These
# fixtures make absence covariate-linked: only participants whose `risk` is zero
# have native support near the third group node.

pcov_effects <- function() {
  effect_space(c("face", "house"), basis_id = "population-coverage:v1")
}

pcov_subject <- function(id, features) {
  domain <- abstract_domain(features,
    coordinates = cbind(x = seq_len(features) - 1),
    feature_ids = paste0("f", seq_len(features)), id = id)
  values <- matrix(seq_len(2L * features), 2L, features,
    dimnames = list(c("face", "house"), NULL)) / features
  relation <- relation(
    list(run1 = values, run2 = values * 1.1, run3 = values * 0.9),
    effects = pcov_effects(), domain = domain
  )
  plan_geometry(relation, compile_frame(voxelwise(), domain),
    cross_partitions(relation))
}

pcov_plan <- function(coverage_policy = "all_planned",
                      group_coords = c(0, 4, 9)) {
  sizes <- c(s01 = 6L, s02 = 8L, s03 = 10L, s04 = 12L)
  subjects <- stats::setNames(lapply(names(sizes), function(id) {
    pcov_subject(id, sizes[[id]])
  }), names(sizes))
  transport <- stats::setNames(lapply(names(sizes), function(id) {
    anatomical_transport(
      native_coords = cbind(seq_len(sizes[[id]]) - 1),
      group_coords = cbind(group_coords), semantics = "budget", radius = 1.5
    )
  }), names(sizes))
  covariates <- data.frame(
    risk = c(1, 2, 0, 0), row.names = names(sizes)
  )
  plan_population(subjects, transport, model = ~ risk, data = covariates,
    coverage_policy = coverage_policy)
}

pcov_bank <- function() rbind(`face-house` = c(1, -1))

test_that("all_planned makes an incompletely covered node unresolved", {
  plan <- pcov_plan()
  fit <- estimate_population(plan, pcov_bank())
  coverage <- fit$coverage

  expect_identical(coverage$contract_version, "population-estimand-v1")
  expect_identical(coverage$policy, "all_planned")
  expect_identical(coverage$planned_subjects, names(plan$subjects))
  expect_identical(coverage$subject_plan_id,
    stats::setNames(plan$subject_index$plan_id, names(plan$subjects)))
  expect_identical(coverage$transport_signature,
    stats::setNames(plan$subject_index$transport_signature,
      names(plan$subjects)))

  expect_identical(coverage$n["group3", "face-house"], 2L)
  expect_identical(
    coverage$status["group3", "face-house"],
    "planned_subject_unavailable"
  )
  expect_true(all(is.na(fit$coefficients["group3", "face-house", ])))
  expect_false(any(coverage$coefficient_estimable[
    "group3", "face-house", ]
  ))
  expect_true(all(coverage$exclusion_reason[
    "group3", "face-house", ] == "planned_subject_unavailable"
  ))

  key <- coverage$subject_set_id["group3", "face-house"]
  expect_identical(coverage$subject_sets[[key]], c("s03", "s04"))
  expect_identical(as.character(coverage$availability[
    "group3", "face-house", ]), c("FALSE", "FALSE", "TRUE", "TRUE"))
  expect_true(all(is.finite(fit$coefficients["group1", "face-house", ])))
  expect_identical(coverage$n["<sink>", "face-house"], 4L)
})

test_that("available_at_node records a changed target and local rank", {
  plan <- pcov_plan("available_at_node")
  fit <- estimate_population(plan, pcov_bank())
  coverage <- fit$coverage

  expect_identical(coverage$policy, "available_at_node")
  expect_identical(coverage$n["group3", "face-house"], 2L)
  expect_identical(coverage$design_rank["group3", "face-house"], 1L)
  expect_identical(coverage$residual_df["group3", "face-house"], 1L)
  expect_identical(coverage$status["group3", "face-house"], "rank_deficient")
  expect_identical(
    as.logical(coverage$coefficient_estimable["group3", "face-house", ]),
    c(TRUE, FALSE)
  )
  expect_true(is.finite(fit$coefficients[
    "group3", "face-house", "(Intercept)"
  ]))
  expect_true(is.na(fit$coefficients["group3", "face-house", "risk"]))
  expect_true(is.na(coverage$exclusion_reason[
    "group3", "face-house", "(Intercept)"
  ]))
  expect_identical(coverage$exclusion_reason[
    "group3", "face-house", "risk"
  ], "coefficient_not_estimable")

  expect_identical(coverage$n["group1", "face-house"], 4L)
  expect_identical(coverage$design_rank["group1", "face-house"], 2L)
  expect_identical(coverage$residual_df["group1", "face-house"], 2L)
  expect_identical(coverage$status["group1", "face-house"], "estimated")
})

test_that("cellwise uncertainty refuses a locally rank-deficient design", {
  fit <- estimate_population(pcov_plan("available_at_node"), pcov_bank())
  robust <- population_uncertainty(fit, estimator = "HC3")

  expect_identical(robust$between$status["group1", "face-house"],
    "estimated")
  expect_identical(robust$between$residual_df["group1", "face-house"], 2L)
  expect_identical(robust$between$status["group3", "face-house"],
    "refused")
  expect_identical(robust$between$reason["group3", "face-house"],
    "rank_deficient_design")
  expect_true(all(is.na(robust$between$covariance[
    "group3", "face-house", , ]
  )))
  expect_true(is.finite(fit$coefficients[
    "group3", "face-house", "(Intercept)"
  ]))
  expect_true(is.na(robust$between$se[
    "group3", "face-house", "(Intercept)"
  ]))
})

test_that("empty node coverage is distinct and subject sets are recoverable", {
  fit <- estimate_population(
    pcov_plan("available_at_node", group_coords = c(0, 4, 20)),
    pcov_bank()
  )
  coverage <- fit$coverage

  expect_identical(coverage$n["group3", "face-house"], 0L)
  expect_identical(coverage$design_rank["group3", "face-house"], 0L)
  expect_identical(coverage$residual_df["group3", "face-house"], 0L)
  expect_identical(coverage$status["group3", "face-house"],
    "empty_subject_set")
  key <- coverage$subject_set_id["group3", "face-house"]
  expect_identical(coverage$subject_sets[[key]], character())
  expect_true(all(is.na(fit$coefficients["group3", "face-house", ])))

  broken <- fit
  broken$coverage$subject_sets[[key]] <- "s01"
  expect_error(crossform:::.validate_population_result(broken),
    "not recoverable", class = "effect_contract_error")
})

test_that("coverage identity is stable across execution routes", {
  plan <- pcov_plan("available_at_node")
  first <- estimate_population(plan, pcov_bank())
  second <- estimate_population(plan, pcov_bank())
  complete <- materialize_population(plan)

  expect_identical(first$coverage$subject_set_id,
    second$coverage$subject_set_id)
  expect_identical(first$coverage$subject_sets,
    second$coverage$subject_sets)
  expect_identical(first$coverage$operator_coverage,
    complete$coverage$operator_coverage)
  expect_identical(first$coverage$subject_plan_id,
    complete$coverage$subject_plan_id)
  expect_true(all(complete$coverage$n[, ] ==
    first$coverage$n[, "face-house"]))
  expect_silent(crossform:::.validate_population_result(complete))
})

test_that("tables expose the exact target and exclusions", {
  fit <- estimate_population(pcov_plan("available_at_node"), pcov_bank())
  table <- as.data.frame(fit)
  row <- table$node == "group3" & table$term == "risk"

  expect_true(all(c(
    "coverage_policy", "planned_n", "n", "fraction", "n_eff",
    "mass_n_eff", "design_rank", "residual_df", "coverage_status",
    "coefficient_estimable", "exclusion_reason", "subject_set_id",
    "available_subjects", "excluded_subjects"
  ) %in% names(table)))
  expect_identical(table$coverage_policy[row], "available_at_node")
  expect_identical(table$planned_n[row], 4L)
  expect_identical(table$n[row], 2L)
  expect_false(table$coefficient_estimable[row])
  expect_identical(table$available_subjects[row], "s03,s04")
  expect_identical(table$excluded_subjects[row], "s01,s02")
})

test_that("views preserve exact coverage and refuse incompatible aggregation", {
  fit <- estimate_population(pcov_plan("available_at_node"), pcov_bank())
  selected <- contrast_energy(fit, c(face = 1, house = -1),
    term = "(Intercept)")

  expect_identical(selected$coverage$subject_set_id["group3", "energy"],
    fit$coverage$subject_set_id["group3", "face-house"])
  expect_identical(selected$coverage$subject_sets,
    fit$coverage$subject_sets)
  expect_identical(selected$coverage$n["group3", "energy"], 2L)
  expect_identical(as.data.frame(selected)$available_subjects[
    as.data.frame(selected)$node == "group3"
  ], "s03,s04")

  preserved <- contribution(selected, by = c("a", "a", "b"))
  expect_identical(preserved$coverage$n["a", "energy"], 4L)
  expect_identical(preserved$coverage$n["b", "energy"], 2L)
  expect_silent(crossform:::.validate_population_view(preserved))

  refusal <- catch_refusal(contribution(selected, by = c("a", "b", "b")))
  expect_identical(refusal$capability, "common_population_subject_set")
  expect_true(any(grepl("source_subject_sets_differ", refusal$reasons,
    fixed = TRUE)))
})
