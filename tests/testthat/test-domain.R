test_that("abstract domains retain explicit feature identity and coordinates", {
  domain <- abstract_domain(3, coordinates = cbind(x = 0:2, y = 1),
    feature_ids = c("a", "b", "c"), id = "native:s01")

  expect_s3_class(domain, "effect_domain")
  expect_identical(domain$n_features, 3L)
  expect_identical(domain$feature_ids, c("a", "b", "c"))
  expect_s3_class(domain$reference, "effect_domain_reference")
  expect_match(domain$reference$signature, "^sha256:[[:xdigit:]]{64}$")
  expect_error(abstract_domain(3, feature_ids = c("a", "a", "b")), "uniquely")
})

test_that("domain references include ordered features, units, and geometry", {
  first <- abstract_domain(3, coordinates = cbind(x = 0:2, y = 1),
    feature_ids = c("a", "b", "c"), id = "native:s01",
    coordinate_units = "mm")
  reversed <- abstract_domain(3, coordinates = cbind(x = 2:0, y = 1),
    feature_ids = c("c", "b", "a"), id = "native:s01",
    coordinate_units = "mm")
  other_units <- abstract_domain(3, coordinates = cbind(x = 0:2, y = 1),
    feature_ids = c("a", "b", "c"), id = "native:s01",
    coordinate_units = "voxel")

  expect_false(identical(first$reference$signature,
    reversed$reference$signature))
  expect_false(identical(first$reference$signature,
    other_units$reference$signature))
})

test_that("volume domains preserve mask order and physical coordinates", {
  mask <- array(FALSE, c(2, 3, 2))
  mask[c(1, 4, 12)] <- TRUE
  domain <- volume_domain(mask, spacing = c(2, 3, 4), id = "native-volume:s01")

  expect_identical(domain$n_features, 3L)
  expect_identical(domain$feature_ids, c(1L, 4L, 12L))
  expect_equal(domain$coordinates,
    sweep(domain$metadata$voxel - 1, 2, c(2, 3, 4), `*`), tolerance = 0)
  expect_error(volume_domain(array(FALSE, c(2, 2, 2))), "at least one")
})

test_that("mutated domains fail closed", {
  domain <- abstract_domain(2, coordinates = matrix(1:4, 2))
  domain$n_features <- 3L
  expect_error(crossform:::.validate_domain(domain), "identities|coordinates")
})
