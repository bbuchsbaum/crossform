test_that("effect spaces have stable semantic identity", {
  first <- effect_space(c("face", "house"), basis_id = "conditions:v1",
    units = "percent-signal", scale = c(1, 2),
    provenance = list(contrast_set = "stimulus"))
  same <- effect_space(c("face", "house"), basis_id = "conditions:v1",
    units = "percent-signal", scale = c(1, 2),
    provenance = list(contrast_set = "stimulus"))

  expect_s3_class(first, "effect_space")
  expect_match(first$signature, "^sha256:[[:xdigit:]]{64}$")
  expect_identical(first, same)
  expect_false(identical(first$signature,
    effect_space(c("face", "house"), basis_id = "conditions:v2")$signature))
  expect_false(identical(first$signature,
    effect_space(c("house", "face"), basis_id = "conditions:v1",
      units = "percent-signal", scale = c(2, 1),
      provenance = list(contrast_set = "stimulus"))$signature))
})

test_that("effect-space validators reject ambiguous and mutated identity", {
  expect_error(effect_space(character()), "nonempty character vector")
  expect_error(effect_space(c("a", "a")), "unique")
  expect_error(effect_space(c("a", "b"), units = c("z", "")), "unit")
  expect_error(effect_space(c("a", "b"), scale = c(1, 0)), "positive")
  expect_error(effect_space("a", provenance = list(fun = identity)),
    "nonportable")

  space <- effect_space(c("a", "b"), basis_id = "basis:v1")
  space$units[[1L]] <- "other"
  expect_error(crossform:::.validate_effect_space(space), "inconsistent")
})
