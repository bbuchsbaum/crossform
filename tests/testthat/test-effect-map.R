effect_map_fixture <- function() {
  conditions <- condition_space(
    c("face", "place", "object"),
    basis_id = "canonical-hrf-amplitude",
    units = "percent-signal-change",
    provenance = list(protocol = "fixture-v1")
  )
  weights <- rbind(
    place_minus_object = c(0, 1, -1),
    face_vs_others = c(1, -0.5, -0.5)
  )
  colnames(weights) <- conditions$coordinates
  effects <- effect_map(weights, conditions)

  cell_map <- cbind(
    face = c(1, 0, 0),
    place = c(0, 1, 0),
    object = c(0, 0, 1),
    drift = c(0, 0, 0)
  )
  rownames(cell_map) <- conditions$coordinates
  treatment_map <- cbind(
    intercept = c(1, 1, 1),
    place_minus_face = c(0, 1, 0),
    object_minus_face = c(0, 0, 1),
    drift = c(0, 0, 0)
  )
  rownames(treatment_map) <- conditions$coordinates

  list(
    conditions = conditions,
    effects = effects,
    cell = coefficient_parameterization(
      cell_map,
      conditions,
      coding_id = "cell-means",
      provenance = list(compiler = "fixture")
    ),
    treatment = coefficient_parameterization(
      treatment_map,
      conditions,
      coding_id = "treatment",
      provenance = list(compiler = "fixture")
    )
  )
}

test_that("condition spaces bind order basis units scale and provenance", {
  base <- condition_space(c("a", "b"), basis_id = "canonical")
  expect_s3_class(base, "effect_condition_space")
  expect_match(format(base), "2; canonical")

  variants <- list(
    condition_space(c("b", "a"), basis_id = "canonical"),
    condition_space(c("a", "b"), basis_id = "fir-bin-1"),
    condition_space(c("a", "b"), basis_id = "canonical", units = "psc"),
    condition_space(c("a", "b"), basis_id = "canonical", scale = 100),
    condition_space(c("a", "b"), basis_id = "canonical",
      provenance = list(model = "v2"))
  )
  expect_true(all(vapply(variants, function(value) {
    !identical(value$signature, base$signature)
  }, logical(1))))
})

test_that("effect maps are semantic and independent of coefficient coding", {
  fixture <- effect_map_fixture()
  cell <- lower_effect_map(fixture$effects, fixture$cell)
  treatment <- lower_effect_map(fixture$effects, fixture$treatment)

  expect_s3_class(fixture$effects, "effect_condition_map")
  expect_identical(cell$effect_map_id, treatment$effect_map_id)
  expect_identical(cell$effect_space, treatment$effect_space)
  expect_false(identical(cell$parameterization_id,
    treatment$parameterization_id))
  expect_false(identical(cell$lowering_id, treatment$lowering_id))
  expect_true(cell$capabilities$coding_invariant)
  expect_match(format(fixture$cell), "3 conditions x 4 coefficients")
  expect_match(format(cell), "2 effects x 4 coefficients")
  expect_identical(colnames(cell$target), fixture$cell$coefficients)
  expect_identical(colnames(treatment$target),
    fixture$treatment$coefficients)
})

test_that("semantic lowerings reproduce the same linear extractor", {
  fixture <- effect_map_fixture()
  cell <- lower_effect_map(fixture$effects, fixture$cell)
  treatment <- lower_effect_map(fixture$effects, fixture$treatment)

  set.seed(2026081503)
  condition <- rep(c("face", "place", "object"), length.out = 17L)
  semantic <- stats::model.matrix(~ condition - 1)
  colnames(semantic) <- sub("condition", "", colnames(semantic))
  semantic <- semantic[, fixture$conditions$coordinates, drop = FALSE]
  drift <- seq(-0.91, 1.13, length.out = nrow(semantic))
  cell_design <- cbind(semantic, drift = drift)
  coding <- rbind(
    c(1, 0, 0, 0),
    c(1, 1, 0, 0),
    c(1, 0, 1, 0),
    c(0, 0, 0, 1)
  )
  treatment_design <- cell_design %*% coding
  colnames(treatment_design) <- fixture$treatment$coefficients
  whitener <- diag(seq(0.79, 1.21, length.out = nrow(semantic)))

  cell_extractor <- lm_extractor(
    cell_design,
    cell$target,
    observation_whitener = whitener,
    effect_names = cell$effect_space
  )
  treatment_extractor <- lm_extractor(
    treatment_design,
    treatment$target,
    observation_whitener = whitener,
    effect_names = treatment$effect_space
  )
  response <- matrix(rnorm(nrow(semantic) * 9L), nrow(semantic), 9L)

  expect_equal(cell_extractor$map, treatment_extractor$map,
    tolerance = 1e-12)
  expect_equal(cell_extractor$map %*% response,
    treatment_extractor$map %*% response, tolerance = 1e-12)
})

test_that("functional units and HRF targets change semantic identity", {
  fixture <- effect_map_fixture()
  weights <- fixture$effects$weights
  changed_weights <- weights
  changed_weights[1, ] <- c(-1, 1, 0)

  variants <- list(
    effect_map(changed_weights, fixture$conditions),
    effect_map(weights, fixture$conditions, units = "z-score"),
    effect_map(weights, fixture$conditions, component = "fir-bin-2")
  )
  expect_true(all(vapply(variants, function(value) {
    !identical(value$effect_map_id, fixture$effects$effect_map_id)
  }, logical(1))))
})

test_that("condition-space mismatches refuse with typed remedies", {
  fixture <- effect_map_fixture()
  changed <- condition_space(
    fixture$conditions$coordinates,
    basis_id = "fir-bin-2",
    units = "percent-signal-change"
  )
  map <- diag(3)
  dimnames(map) <- list(changed$coordinates, c("face", "place", "object"))
  parameterization <- coefficient_parameterization(
    map, changed, coding_id = "cell-means"
  )
  refusal <- catch_refusal(lower_effect_map(fixture$effects, parameterization))

  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "valid_effect_lowering")
  expect_identical(refusal$namespace, "relation_plan")
  expect_length(refusal$reasons, 2L)
  expect_true(length(refusal$remedies) >= 1L)
})

test_that("rank-deficient parameterizations refuse before lowering", {
  conditions <- condition_space(c("a", "b", "c"))
  map <- matrix(c(
    1, 0,
    0, 1,
    1, 1
  ), 3L, 2L, byrow = TRUE,
  dimnames = list(conditions$coordinates, c("x", "y")))
  refusal <- catch_refusal(coefficient_parameterization(
    map, conditions, coding_id = "underidentified"
  ))

  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "valid_effect_lowering")
  expect_match(refusal$reasons, "semantic rank 2 for 3")
  expect_true(length(refusal$remedies) >= 1L)
})

test_that("raw targets bind parameterization and expose narrower capability", {
  fixture <- effect_map_fixture()
  cell <- lower_effect_map(fixture$effects, fixture$cell)
  treatment <- lower_effect_map(fixture$effects, fixture$treatment)
  raw_cell <- raw_effect_map(
    cell$target,
    effects = cell$effect_space,
    coefficients = colnames(cell$target),
    provenance = list(source = "external")
  )
  raw_treatment <- raw_effect_map(
    treatment$target,
    effects = treatment$effect_space,
    coefficients = colnames(treatment$target),
    provenance = list(source = "external")
  )

  expect_s3_class(raw_cell, "effect_raw_map")
  expect_false(raw_cell$capabilities$symbolic_effects)
  expect_false(raw_cell$capabilities$coding_invariant)
  expect_false(identical(raw_cell$effect_map_id,
    raw_treatment$effect_map_id))
})

test_that("semantic constructors reject positional and nonportable ambiguity", {
  expect_error(effect_map(matrix(1:4, 2)), "condition space")
  expect_error(condition_space(c("a", "a")), "unique")
  expect_error(condition_space("a", provenance = list(run = function() 1)),
    "nonportable")

  fixture <- effect_map_fixture()
  wrong_order <- fixture$effects$weights[, c("object", "place", "face")]
  expect_error(effect_map(wrong_order, fixture$conditions), "must follow")
})

test_that("semantic object validation detects identity tampering", {
  fixture <- effect_map_fixture()
  tampered <- fixture$effects
  tampered$weights[1, 1] <- tampered$weights[1, 1] + 1
  expect_error(crossform:::.validate_effect_map(tampered), "inconsistent")

  tampered_space <- fixture$conditions
  tampered_space$basis_id <- "changed"
  expect_error(crossform:::.validate_condition_space(tampered_space),
    "inconsistent")

  lowered <- lower_effect_map(fixture$effects, fixture$cell)
  lowered$target[1, 1] <- lowered$target[1, 1] + 1
  expect_error(crossform:::.validate_lowered_effect_map(lowered),
    "identity is inconsistent")
})
