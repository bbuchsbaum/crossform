test_that("the 0.1 namespace exposes only implemented workflows", {
  exports <- getNamespaceExports("crossform")

  expect_true(all(c(
    "effect_space", "pair_query", "abstract_domain", "volume_domain", "relation",
    "compile_frame", "geometry", "evaluate_geometry", "contrast", "rdm",
    "rsa", "geometry_spectrum", "rdm_sampling_covariance",
    "sampling_covariance", "neuroim2_volume_domain",
    "neuroim2_searchlights"
  ) %in% exports))
  expect_false(any(c(
    "frame", "factor_frame", "nonlinear_query", "compile_lowering",
    "effect_geometry", "effect_view", "execution_receipt", "memory_plan"
  ) %in% exports))
})
