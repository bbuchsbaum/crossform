test_that("injected subject influence is localized and provenance is retained", {
  plan <- pdec_fixture(outlier_subject = "s04")
  decomposition <- do.call(population_decomposition, pdec_results(plan))
  split <- heterogeneity(plan, estimator = "cross_fit")
  influence <- population_influence(decomposition, heterogeneity = split)

  aggregate <- stats::aggregate(abs_delta ~ subject, influence$influence,
                                sum, na.rm = TRUE)
  expect_identical(aggregate$subject[which.max(aggregate$abs_delta)], "s04")
  expect_identical(influence$heterogeneity$estimator, "cross_fit")
  loading <- abs(influence$heterogeneity$leading_loading)
  expect_identical(names(which.max(loading)), "s04")
  expect_false(any(influence$subject_provenance$transport_self_fit))
  expect_true(all(nzchar(influence$influence$transport_signature)))
  expect_true(all(influence$influence$primary_subject_set_id != ""))
  expect_match(influence$interpretation, "descriptive")
})

test_that("influence links subject coverage transport and affected cells", {
  decomposition <- do.call(population_decomposition, pdec_results())
  influence <- population_influence(decomposition)

  expect_true(all(c("subject", "node", "scale", "query", "term", "component",
                    "contributed_to_primary", "primary_subject_set_id",
                    "transport_signature", "cross_fit", "sink_territory",
                    "delta", "abs_delta") %in% names(influence$influence)))
  absent <- !influence$influence$contributed_to_primary
  expect_true(any(absent))
  expect_true(all(is.na(influence$influence$leave_one_out_estimate[absent]) |
                    influence$influence$delta[absent] == 0))
  expect_identical(influence$mode, "bounded")
})

test_that("bounded defaults require an explicit deeper mode", {
  decomposition <- do.call(population_decomposition, pdec_results())
  expect_error(population_influence(decomposition, max_subjects = 4L),
               class = "effect_capability_refusal")
  deep <- population_influence(decomposition, mode = "deep", max_subjects = 4L)
  expect_identical(deep$mode, "deep")
  expect_gt(deep$bounds$requested_cells, 0)
})
