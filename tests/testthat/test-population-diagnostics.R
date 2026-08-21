pdiag_fixture <- function() {
  effects <- effect_space(c("face", "house"), basis_id = "pdiag:v1")
  sizes <- c(s01 = 6L, s02 = 8L, s03 = 10L, s04 = 12L, s05 = 12L, s06 = 12L)
  gains <- c(s01 = 2.0, s02 = 1.7, s03 = 1.0, s04 = 0.8, s05 = 0.7, s06 = 0.6)
  subjects <- stats::setNames(lapply(names(sizes), function(id) {
    n <- sizes[[id]]
    domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
      feature_ids = paste0("f", seq_len(n)), id = id)
    values <- matrix(gains[[id]] * seq_len(2L * n), 2L, n,
      dimnames = list(c("face", "house"), NULL)) / n
    relation <- relation(list(run1 = values, run2 = values * 1.1),
      effects = effects, domain = domain)
    plan_geometry(relation, compile_frame(voxelwise(), domain),
                  cross_partitions(relation))
  }), names(sizes))
  transports <- stats::setNames(lapply(names(sizes), function(id) {
    anatomical_transport(cbind(seq_len(sizes[[id]]) - 1),
      cbind(c(0, 5, 11)), semantics = "budget", radius = 1.2)
  }), names(sizes))
  covariates <- data.frame(risk = c(2, 1.5, 0.5, 0, -0.5, -1),
                           row.names = names(sizes))
  plan <- plan_population(subjects, transports, model = ~ risk,
    data = covariates, coverage_policy = "available_at_node")
  estimate_population(plan, rbind(`face-house` = c(1, -1)))
}

test_that("diagnostics preserve provenance and detect injected associations", {
  fit <- pdiag_fixture()
  diagnostics <- population_diagnostics(fit, minimum_coverage = 0.9,
    minimum_transport_quality = 0.75, material_change = 0.15)

  expect_s3_class(diagnostics, "effect_population_diagnostics")
  expect_identical(diagnostics$parent_result_id, fit$scientific_plan_id)
  expect_identical(diagnostics$subject_provenance$subject,
                   fit$coverage$planned_subjects)
  expect_true(all(nzchar(diagnostics$subject_provenance$plan_id)))
  expect_true(any(abs(diagnostics$associations$availability_smd) > 0.5,
                  na.rm = TRUE))
  expect_true(any(abs(diagnostics$associations$outcome_correlation) > 0.5,
                  na.rm = TRUE))
  expect_true(any(diagnostics$warnings$warning ==
                    "coverage_below_declared_minimum"))
})

test_that("threshold sensitivity never changes the primary result", {
  fit <- pdiag_fixture()
  original <- serialize(fit, NULL)
  diagnostics <- population_diagnostics(fit,
    minimum_transport_quality = 0.95, material_change = 0.1)

  expect_identical(serialize(fit, NULL), original)
  expect_true(all(diagnostics$sensitivity$target_status ==
                    "sensitivity_descriptive_not_primary"))
  expect_true(any(diagnostics$sensitivity$removed_fraction > 0))
  expect_true(any(diagnostics$warnings$warning ==
                    "group_composition_changes_materially"))
  expect_match(diagnostics$interpretation, "does not replace")
})

test_that("diagnostic thresholds are explicit and bounded", {
  fit <- pdiag_fixture()
  expect_error(population_diagnostics(fit, minimum_coverage = 1.1),
               class = "effect_input_error")
  expect_error(population_diagnostics(fit, material_change = -0.1),
               class = "effect_input_error")
})
