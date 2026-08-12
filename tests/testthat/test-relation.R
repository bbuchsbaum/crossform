test_that("precomputed effect matrices form an identity relation", {
  effects <- c("face", "house")
  runs <- list(
    run1 = matrix(1:10, 2, 5, dimnames = list(effects, NULL)),
    run2 = matrix(11:20, 2, 5, dimnames = list(effects, NULL))
  )
  rel <- relation(runs, effects = effects, domain_id = "native:s01")

  expect_s3_class(rel, "effect_relation")
  expect_identical(rel$partitions, c("run1", "run2"))
  expect_identical(rel$n_features, 5L)
  expect_equal(relation_block(rel, "run2", c(2, 5)),
    runs$run2[, c(2, 5)], tolerance = 0)
})

test_that("named precomputed partitions align to the declared effect space", {
  space <- effect_space(c("face", "house"), basis_id = "conditions:v1",
    units = "percent-signal")
  canonical <- matrix(1:6, 2, 3,
    dimnames = list(c("face", "house"), NULL))
  reversed <- canonical[c("house", "face"), , drop = FALSE]
  rel <- relation(list(first = canonical, second = reversed), effects = space)

  expect_identical(rel$effect_space, space)
  expect_equal(relation_block(rel, "second", 1:3), canonical, tolerance = 0)
})

test_that("precomputed coordinate ambiguity fails during construction", {
  space <- effect_space(c("a", "b"), basis_id = "basis:v1")
  missing <- matrix(1:6, 2, 3, dimnames = list(c("a", "c"), NULL))
  partial <- matrix(1:6, 2, 3, dimnames = list(c("a", ""), NULL))

  expect_error(relation(list(run = missing), effects = space),
    "missing or extra")
  expect_error(relation(list(run = partial), effects = space),
    "complete and unique")
  expect_error(relation(list(run = unname(matrix(1:6, 2, 3)))),
    "complete row names")
})

test_that("extractor bases and units must match across partitions", {
  first <- effect_extractor(diag(2),
    effect_space(c("a", "b"), basis_id = "basis:v1", units = "z"))
  other_basis <- effect_extractor(diag(2),
    effect_space(c("a", "b"), basis_id = "basis:v2", units = "z"))
  other_units <- effect_extractor(diag(2),
    effect_space(c("a", "b"), basis_id = "basis:v1", units = "percent"))
  sources <- list(one = matrix(1, 2, 3), two = matrix(1, 2, 3))

  expect_error(relation(sources, extract = list(first, other_basis)),
    "identical effect space")
  expect_error(relation(sources, extract = list(first, other_units)),
    "identical effect space")
})

test_that("raw response and materialized effect relations are equivalent", {
  set.seed(61)
  design <- cbind(1, condition = rep(c(-0.5, 0.5), 6))
  target <- rbind(level = c(1, 0), condition = c(0, 1))
  extractor <- lm_extractor(design, target)
  raw <- list(run1 = matrix(rnorm(12 * 7), 12, 7))
  materialized <- list(run1 = extractor$map %*% raw$run1)
  raw_relation <- relation(raw, extract = extractor)
  effect_relation <- relation(materialized, effects = c("level", "condition"))

  for (features in list(1:7, c(1, 3, 7))) {
    expect_equal(
      relation_block(raw_relation, "run1", features),
      relation_block(effect_relation, "run1", features),
      tolerance = 1e-13
    )
  }
})

test_that("function-backed sources are lazy and validated", {
  matrix_source <- matrix(seq_len(30), 6, 5)
  reads <- 0L
  source <- function(features) {
    reads <<- reads + 1L
    matrix_source[, features, drop = FALSE]
  }
  extractor <- effect_extractor(matrix(c(1, 0, 0, 0, 0, -1), 1), "contrast")
  rel <- relation(list(run1 = source), extract = extractor,
    source_dims = list(c(6, 5)))

  expect_identical(reads, 0L)
  expect_error(relation_block(rel, "run1", 0), "valid neural feature")
  expect_identical(reads, 0L)
  got <- relation_block(rel, "run1", c(2, 4))
  expect_identical(reads, 1L)
  expect_equal(got, extractor$map %*% matrix_source[, c(2, 4)], tolerance = 0)
})

test_that("relation construction rejects incompatible partitions", {
  expect_error(relation(list(matrix(1, 2, 3), matrix(1, 2, 4)), effects = c("a", "b")),
    "feature dimension")

  extractor <- effect_extractor(diag(2), c("a", "b"))
  expect_error(relation(list(matrix(1, 3, 4)), extract = extractor),
    "observations")
})

test_that("mutated relation metadata fails before source reads", {
  reads <- 0L
  source <- function(features) {
    reads <<- reads + 1L
    matrix(1, 2, length(features))
  }
  rel <- relation(list(run1 = source), effects = c("a", "b"),
    source_dims = list(c(2, 3)))
  rel$n_features <- 4L

  expect_error(relation_block(rel, "run1", 1), "metadata")
  expect_identical(reads, 0L)
})

test_that("matrix sources receive strong deterministic revision identities", {
  source <- matrix(1:12, 3, 4)
  first <- relation(list(run1 = source), effects = c("a", "b", "c"))
  second <- relation(list(run1 = source), effects = c("a", "b", "c"))
  changed <- source
  changed[1, 1] <- 99
  third <- relation(list(run1 = changed), effects = c("a", "b", "c"))

  expect_match(first$capabilities$run1$stable_revision,
    "^sha256:[[:xdigit:]]{64}$")
  expect_identical(first$capabilities$run1$stable_revision,
    second$capabilities$run1$stable_revision)
  expect_false(identical(first$capabilities$run1$stable_revision,
    third$capabilities$run1$stable_revision))
})

test_that("a relation accepts a domain without exposing identity plumbing", {
  domain <- abstract_domain(4, id = "native:s01")
  rel <- relation(list(run1 = matrix(1, 2, 4)), effects = c("a", "b"),
    domain = domain)
  expect_identical(rel$domain_id, domain$id)
  expect_identical(rel$domain, domain$reference)
  expect_error(relation(list(run1 = matrix(1, 2, 3)), effects = c("a", "b"),
    domain = domain), "feature count")
  expect_error(relation(list(run1 = matrix(1, 2, 4)), effects = c("a", "b"),
    domain = domain, domain_id = "other"), "different neural domains")
})
