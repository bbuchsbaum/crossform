test_that("every displayed effect is bound to exact support diagnostics", {
  fits <- pdec_results()
  decomposition <- do.call(population_decomposition, fits)
  effects <- population_scale_profile(decomposition, "age", interval = "HC3")
  diagnostics <- population_diagnostics(fits$total,
    minimum_coverage = 0.95, minimum_transport_quality = 0.95,
    material_change = 0.1)
  view <- population_diagnostic_view(effects, diagnostics)

  expect_s3_class(view, "effect_population_diagnostic_view")
  expect_false(anyNA(view$effect$support_key))
  expect_true(all(c("coverage_fraction", "n", "mass_n_eff",
                    "mean_sink_territory", "mean_transport_quality") %in%
                    names(view$effect)))
  expect_true(all(table(view$effect$support_key) == 3L))
  expect_setequal(names(view$panels),
                  c("coverage", "effective_n", "sink", "transport_quality"))
  expect_identical(length(unique(vapply(view$panels, `[[`, character(1),
                                               "scale_identity"))), 4L)
  expect_true(any(view$warnings$warning ==
                    "coverage_below_declared_minimum"))
  expect_true(any(view$warnings$warning ==
                    "sink_territory_above_declared_maximum"))
})

test_that("effect and diagnostic filters are synchronized and recorded", {
  fits <- pdec_results()
  decomposition <- do.call(population_decomposition, fits)
  effects <- population_component_view(decomposition, "(Intercept)")
  diagnostics <- population_diagnostics(fits$total)
  selected <- population_diagnostic_view(effects, diagnostics,
    node = unique(effects$data$node)[1:2], query = "face-house")

  expect_setequal(unique(selected$effect$node),
                  unique(effects$data$node)[1:2])
  expect_identical(unique(selected$effect$query), "face-house")
  expect_setequal(unique(selected$support$node),
                  unique(effects$data$node)[1:2])
  expect_identical(unique(selected$support$query), "face-house")
  expect_identical(selected$filters$operation,
                   "synchronized_exact_key_selection")
  expect_lt(selected$filters$effect_rows_after,
            selected$filters$effect_rows_before)
  expect_match(selected$interpretation, "separate units")
})

test_that("diagnostic views refuse effects from another population plan", {
  first <- pdec_results()
  second <- pdec_results(pdec_fixture(coverage_policy = "all_planned"))
  view <- population_component_view(
    do.call(population_decomposition, first), "age")
  diagnostics <- population_diagnostics(second$total)
  expect_error(population_diagnostic_view(view, diagnostics),
               class = "effect_capability_refusal")
})
