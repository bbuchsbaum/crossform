test_that("portable relation receipts satisfy the conformance court", {
  fixture <- relation_plan_fixture()
  court <- compiler_conformance(fixture$plan)

  expect_identical(court$partition, fixture$plan$partitions)
  fields <- setdiff(names(court), "partition")
  expect_true(all(vapply(court[fields], is.logical, logical(1))))
  expect_true(all(as.matrix(court[fields])))
})

test_that("compiler version admission fails with a named capability", {
  refusal <- catch_refusal(crossform:::.require_adapter_version(
    "fmridesign", "999.0.0"
  ))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "supported_compiler_version")
  expect_match(conditionMessage(refusal), "not been certified")
})
