test_that("the documented first workflow runs through every public layer", {
  effects <- effect_space(c("face", "house", "object"),
    basis_id = "condition-means:v1", units = "percent-signal")
  domain <- abstract_domain(6,
    coordinates = cbind(x = 0:5, y = 0),
    feature_ids = paste0("feature", 1:6), id = "native:demo")
  run1 <- rbind(
    face = c(2, 1, 0.5, -0.5, 0, 1),
    house = c(0, 1, 2, 1, 0.5, 0),
    object = c(0.5, 0.5, 1, 1.5, 1, 0.5)
  )
  run2 <- rbind(
    face = c(2.2, 0.8, 0.4, -0.4, 0.1, 0.9),
    house = c(0.1, 1.2, 1.8, 0.9, 0.6, 0.1),
    object = c(0.4, 0.6, 1.1, 1.4, 0.9, 0.6)
  )
  rel <- relation(list(run1 = run1, run2 = run2), effects = effects,
    domain = domain)
  at <- compile_frame(searchlights(radius = 1.01, normalization = "local"), domain)
  over <- cross_partitions(rel)
  g <- geometry(rel, at, over)
  result <- contrast(g, c(face = 1, house = -1, object = 0))

  expected <- matrix(c(
    0.925, 0.850, 1.250, 2.100,
    0.133333333333333, 0.0166666666666667, 2.08333333333333, 2.100,
    -1.01666666666667, 1.03333333333333, 0.316666666666667, 1.350,
    -1.11666666666667, 1.24444444444444, 0.188888888888889, 1.43333333333333,
    -0.333333333333333, 0.111111111111111, 0.888888888888889, 1.000,
    0.200, 0.0375, 0.4875, 0.525
  ), 6, 4, byrow = TRUE)
  observed <- cbind(result$signed, result$coherent,
    result$configuration, result$total)
  expect_equal(observed, expected, tolerance = 1e-12)
  expect_equal(result$total, result$coherent + result$configuration,
    tolerance = 0)

  query <- bilinear_query(tcrossprod(c(1, -1, 0)), effects = effects)
  expect_equal(
    evaluate_geometry(rel, at, over, query)$values,
    query_geometry(g, query)$values,
    tolerance = 1e-12
  )
  expect_s3_class(rdm(g), "effect_rdm_view")
  expect_s3_class(geometry_spectrum(g), "effect_spectrum_view")
})

test_that("raw extraction and materialized effects compile to equal geometry", {
  design <- cbind(level = 1, condition = c(-1, -1, 1, 1))
  extractor <- lm_extractor(design, diag(2),
    effect_names = c("level", "condition"))
  raw <- list(
    run1 = matrix(seq_len(20) / 10, 4, 5),
    run2 = matrix(rev(seq_len(20)) / 11, 4, 5)
  )
  effects <- lapply(raw, function(value) extractor$map %*% value)
  domain <- abstract_domain(5, id = "raw-effect-law")
  raw_relation <- relation(raw, extract = extractor, domain = domain)
  effect_relation <- relation(effects, effects = extractor$effects, domain = domain)
  at <- compile_frame(whole_brain(), domain)
  over <- cross_partitions(raw_relation)

  from_raw <- geometry(raw_relation, at, over)
  from_effects <- geometry(effect_relation, at, over)
  for (component in c("total", "coherent", "configuration")) {
    expect_equal(geometry_component(from_raw, component),
      geometry_component(from_effects, component), tolerance = 1e-13)
  }
})

test_that("public frame laws preserve point decomposition and global total", {
  domain <- abstract_domain(5, coordinates = cbind(0:4, 0), id = "frame-laws")
  effects <- list(
    run1 = rbind(a = 1:5, b = c(2, 1, 0, -1, -2)),
    run2 = rbind(a = c(1.2, 1.8, 3.1, 3.9, 5.2),
      b = c(1.8, 1.1, 0.1, -1.2, -1.8))
  )
  rel <- relation(effects, domain = domain)
  over <- cross_partitions(rel)
  point <- geometry(rel, compile_frame(voxels(), domain), over)
  expect_equal(geometry_component(point, "configuration"),
    matrix(0, 5, 3), tolerance = 1e-14)

  local <- geometry(rel,
    compile_frame(searchlights(1.01, normalization = "conservative"), domain), over)
  global <- geometry(rel, compile_frame(whole_brain(normalization = "none"), domain), over)
  expect_equal(colSums(geometry_component(local, "total")),
    drop(geometry_component(global, "total")), tolerance = 1e-13)
})

test_that("public geometry is invariant to joint feature permutation and blocking", {
  set.seed(8102)
  matrices <- list(
    run1 = matrix(rnorm(3 * 7), 3, 7),
    run2 = matrix(rnorm(3 * 7), 3, 7)
  )
  rownames(matrices$run1) <- rownames(matrices$run2) <- c("a", "b", "c")
  weights <- matrix(runif(4 * 7), 4, 7)
  permutation <- sample(seq_len(7))
  original_domain <- abstract_domain(7, feature_ids = paste0("v", seq_len(7)),
    id = "permutation-law")
  permuted_domain <- abstract_domain(7,
    feature_ids = original_domain$feature_ids[permutation],
    id = "permutation-law")
  over <- cross_partitions(c("run1", "run2"))
  original <- geometry(
    relation(matrices, domain = original_domain),
    additive_frame(weights, domain = original_domain), over,
    compute = compute_policy(block_features = 1)
  )
  permuted <- geometry(
    relation(lapply(matrices, function(value) value[, permutation, drop = FALSE]),
      effects = c("a", "b", "c"), domain = permuted_domain),
    additive_frame(weights[, permutation, drop = FALSE],
      domain = permuted_domain), over,
    compute = compute_policy(block_features = 4)
  )
  expect_equal(geometry_component(permuted, "total"),
    geometry_component(original, "total"), tolerance = 1e-12)
  expect_equal(geometry_component(permuted, "coherent"),
    geometry_component(original, "coherent"), tolerance = 1e-12)
})

test_that("the cross-null sign distribution is exactly centered publicly", {
  set.seed(19)
  run1 <- matrix(rnorm(3 * 6), 3, 6)
  run2 <- matrix(rnorm(3 * 6), 3, 6)
  rownames(run1) <- rownames(run2) <- c("a", "b", "c")
  domain <- abstract_domain(6, id = "cross-null")
  at <- compile_frame(whole_brain(), domain)
  positive_relation <- relation(list(run1 = run1, run2 = run2), domain = domain)
  negative_relation <- relation(list(run1 = run1, run2 = -run2), domain = domain)
  positive <- geometry(positive_relation, at, cross_partitions(positive_relation))
  negative <- geometry(negative_relation, at, cross_partitions(negative_relation))

  expect_equal(
    geometry_component(positive, "total") +
      geometry_component(negative, "total"),
    matrix(0, 1, 6), tolerance = 1e-13
  )
})
