test_that("the population coefficient law covers terms nodes queries and missingness", {
  fits <- pdec_results()
  decomposition <- do.call(population_decomposition, fits)

  expect_s3_class(decomposition, "effect_population_decomposition")
  expect_lte(decomposition$max_coefficient_gap, 1e-10)
  expect_lte(decomposition$max_subject_value_gap, 1e-10)
  expect_lte(decomposition$direct_derived_total_covariance_gap, 1e-10)
  expect_identical(dimnames(decomposition$coefficient_gap)[[3L]],
                   c("(Intercept)", "age"))
  expect_identical(dimnames(decomposition$coefficient_gap)[[2L]],
                   c("face-house", "face-tool"))
  expect_true(any(fits$total$coverage$fraction < 1))
  expect_true(any(is.finite(decomposition$component_covariance)))
  expect_match(decomposition$interpretation, "not separate biological")
})

test_that("different plans and nonadditive normalization refuse", {
  fits <- pdec_results()
  other <- pdec_results(pdec_fixture(coverage_policy = "all_planned"))
  expect_error(population_decomposition(fits$total, fits$coherent,
    other$configuration), class = "effect_capability_refusal")

  normalized <- pdec_results(pdec_fixture(normalization = "unit_budget"))
  expect_error(do.call(population_decomposition, normalized),
               class = "effect_capability_refusal")
})

test_that("component roles cannot be relabelled by argument position", {
  fits <- pdec_results()
  expect_error(population_decomposition(fits$coherent, fits$total,
    fits$configuration), class = "effect_capability_refusal")
})
