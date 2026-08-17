test_that("fmrireg independently reproduces the planned OLS point relation", {
  skip_if_not_installed("fmrireg", "0.1.2")
  if (!identical(as.character(utils::packageVersion("fmrireg")), "0.1.2")) {
    skip("The installed fmrireg version is outside the certified court.")
  }
  fixture <- relation_plan_fixture("svd", "treatment", "ols")
  native <- estimate_relation(fixture$plan)
  external <- fmrireg_relation(fixture$plan)

  for (partition in fixture$plan$partitions) {
    expect_equal(
      relation_block(external, partition, 1:5),
      relation_block(native, partition, 1:5),
      tolerance = 1e-12
    )
  }
  capabilities <- relation_fit_capabilities(external)
  expect_false(any(capabilities$residual_blocks))
  expect_false(any(capabilities$effect_covariance))
  expect_identical(external$provenance$relation_plan_id,
    fixture$plan$relation_plan_id)
})

test_that("fmrireg refuses observation models outside its certified slice", {
  skip_if_not_installed("fmrireg", "0.1.2")
  if (!identical(as.character(utils::packageVersion("fmrireg")), "0.1.2")) {
    skip("The installed fmrireg version is outside the certified court.")
  }
  fixture <- relation_plan_fixture("qr", "cell", "fixed_gls")
  refusal <- catch_refusal(fmrireg_relation(fixture$plan))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "supported_observation_model")
  expect_match(conditionMessage(refusal), "OLS")
})
