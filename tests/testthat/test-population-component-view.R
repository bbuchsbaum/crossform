test_that("component views use one signed axis and preserve plotted data", {
  decomposition <- do.call(population_decomposition, pdec_results())
  view <- population_component_view(decomposition, term = "age")

  expect_s3_class(view, "effect_population_component_view")
  expect_setequal(unique(view$data$component),
                  c("total", "coherent", "configuration"))
  expect_identical(view$axis$type, "shared_symmetric")
  expect_equal(view$axis$limits[[1L]], -view$axis$limits[[2L]])
  expect_identical(view$axis$zero, 0)
  expect_identical(view$axis$sign, "positive_up_negative_down")
  expect_true(all(view$data$visual_magnitude ==
                    abs(view$data$estimate) / view$axis$limits[[2L]]))
  equal <- duplicated(abs(view$data$estimate)) |
    duplicated(abs(view$data$estimate), fromLast = TRUE)
  if (any(equal)) {
    expect_equal(view$data$visual_magnitude[equal],
                 abs(view$data$estimate[equal]) / view$axis$limits[[2L]])
  }
  expect_false(any(grepl("ratio", names(view$data), ignore.case = TRUE)))
  expect_identical(view$receipt$aggregation,
    "population_coefficient_not_participant_ratio_average")
})

test_that("component view exposes uncertainty coverage and scale mappings", {
  decomposition <- do.call(population_decomposition, pdec_results())
  view <- population_component_view(decomposition, "(Intercept)",
                                    query = "face-house")

  expect_true(all(c("estimate", "se", "lower", "upper", "sign", "units",
                    "node", "scale", "query", "term") %in% names(view$data)))
  expect_true(all(c("coverage_policy", "subject_set_id", "n", "fraction",
                    "n_eff", "mass_n_eff") %in% names(view$coverage)))
  expect_identical(unique(view$data$query), "face-house")
  expect_true(all(view$data$lower <= view$data$estimate |
                    is.na(view$data$lower)))
  expect_true(all(view$data$upper >= view$data$estimate |
                    is.na(view$data$upper)))
  expect_identical(nrow(view$data), 3L * nrow(view$coverage))
  expect_identical(view$receipt$coverage_preserved, TRUE)
})

test_that("component view selections refuse unknown terms and queries", {
  decomposition <- do.call(population_decomposition, pdec_results())
  expect_error(population_component_view(decomposition, "missing"),
               class = "effect_input_error")
  expect_error(population_component_view(decomposition, "age", "missing"),
               class = "effect_input_error")
})
