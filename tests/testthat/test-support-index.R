brute_ball_supports <- function(domain, radius) {
  lapply(seq_len(domain$n_features), function(center) {
    squared <- rowSums((domain$coordinates - matrix(
      domain$coordinates[center, ],
      nrow = domain$n_features,
      ncol = ncol(domain$coordinates),
      byrow = TRUE
    ))^2)
    which(squared <= radius^2 * (1 + 1e-12))
  })
}

brute_pair_pattern <- function(supports, features) {
  value <- matrix(FALSE, features, features)
  for (support in supports) value[support, support] <- TRUE
  value
}

test_that("coordinate-cell supports and pair patterns match brute force", {
  coordinates <- rbind(
    c(-1.1, 0.0), c(-0.2, 0.1), c(0.7, 0.0), c(1.5, 0.2),
    c(-0.9, 1.1), c(0.1, 1.0), c(1.2, 1.2), c(2.1, 1.0)
  )
  domain <- abstract_domain(
    nrow(coordinates), coordinates = coordinates,
    feature_ids = paste0("f", seq_len(nrow(coordinates))),
    id = "irregular-cells"
  )
  radius <- 1.05
  expected <- brute_ball_supports(domain, radius)
  index <- crossform:::.euclidean_support_index(domain, radius)

  expect_identical(
    crossform:::.support_index_support(index, seq_len(domain$n_features)),
    expected
  )
  expect_identical(
    as.matrix(index$pair_pattern),
    brute_pair_pattern(expected, domain$n_features)
  )
  expect_identical(index$construction$provider, "coordinate_cells")
  expect_identical(index$construction$lookup, "mixed_radix_vectorized")
  expect_silent(crossform:::.validate_support_index(index, deep = TRUE))
})

test_that("coordinate-cell hashing remains exact beyond mixed-radix range", {
  coordinates <- rbind(c(0, 0), c(1, 0), c(1e9, 1e9), c(1e9 + 1, 1e9))
  domain <- abstract_domain(
    nrow(coordinates), coordinates = coordinates, id = "wide-coordinate-cells"
  )
  radius <- 1.01
  expected <- brute_ball_supports(domain, radius)
  index <- crossform:::.euclidean_support_index(domain, radius)

  expect_identical(
    crossform:::.support_index_support(index, seq_len(domain$n_features)),
    expected
  )
  expect_identical(index$construction$lookup, "hashed_tuple")
})

test_that("volume stencils preserve holes, anisotropy, and brute-force order", {
  mask <- array(TRUE, c(5, 4, 3))
  mask[c(2, 7, 18, 35, 46)] <- FALSE
  domain <- volume_domain(mask, spacing = c(1, 1.5, 2.5),
    id = "anisotropic-mask")
  radius <- 2.55
  expected <- brute_ball_supports(domain, radius)
  index <- crossform:::.euclidean_support_index(domain, radius)

  expect_identical(
    crossform:::.support_index_support(index, seq_len(domain$n_features)),
    expected
  )
  expect_identical(
    as.matrix(index$pair_pattern),
    brute_pair_pattern(expected, domain$n_features)
  )
  expect_identical(index$construction$provider, "volume_stencil")
})

test_that("support identities, blocks, gathers, and costs are deterministic", {
  domain <- volume_domain(array(TRUE, c(4, 4, 3)), id = "support-tools")
  first <- crossform:::.euclidean_support_index(domain, 1.01)
  second <- crossform:::.euclidean_support_index(domain, 1.01)
  changed <- crossform:::.euclidean_support_index(domain, 1.5)

  expect_identical(first$signature, second$signature)
  expect_false(identical(first$signature, changed$signature))
  expect_identical(
    crossform:::.support_index_blocks(first, 19L),
    data.frame(
      block = 1:3,
      start = c(1L, 20L, 39L),
      end = c(19L, 38L, 48L),
      stringsAsFactors = FALSE
    )
  )
  support <- crossform:::.support_index_support(first, 7L)[[1L]]
  relation <- matrix(seq_len(3 * domain$n_features), 3)
  expect_identical(
    crossform:::.support_index_gather(first, relation, 7L),
    relation[, support, drop = FALSE]
  )
  expect_equal(first$cost$diagonal_metric_entries,
    length(first$members))
  expect_equal(first$cost$dense_metric_entries,
    sum(diff(first$ptr)^2))
  expect_equal(first$cost$one_dense_factorization_pass_units,
    sum(diff(first$ptr)^3))
  expect_gte(first$cost$pair_pattern_nnz, domain$n_features)
})

test_that("support preflight refuses persistent dense local schedules", {
  domain <- volume_domain(array(TRUE, c(12, 11, 10)), id = "preflight")
  index <- crossform:::.euclidean_support_index(domain, 2.1)
  budget <- index$cost$materialized_dense_metric_bytes * 2
  preflight <- crossform:::.support_index_preflight(
    index, budget, evaluation_edges = 3L
  )

  expect_true(preflight$largest_local_metric_fits)
  expect_true(preflight$materialized_dense_metric_per_edge_fits)
  expect_false(preflight$materialized_dense_metric_fits)
  expect_equal(preflight$materialized_dense_metric_bytes_total,
    3 * preflight$materialized_dense_metric_bytes_per_edge)
  expect_identical(preflight$dense_metric_lowering, "support_streamed_only")
  expect_identical(preflight$factorization_passes,
    "metric_capability_dependent")
  expect_equal(preflight$one_dense_factorization_pass_units_per_edge,
    sum(diff(index$ptr)^3))
  expect_equal(preflight$one_dense_factorization_pass_units_total,
    3 * sum(diff(index$ptr)^3))
  expect_identical(preflight$metric_schedule_storage,
    "derive_on_demand_from_frozen_recipe")
  expect_lt(index$cost$estimated_structural_bytes,
    index$cost$materialized_dense_metric_bytes)
})

test_that("support preflight requires explicit valid edge multiplicity", {
  domain <- volume_domain(array(TRUE, c(3, 3, 2)), id = "edge-preflight")
  index <- crossform:::.euclidean_support_index(domain, 1.01)

  expect_error(
    crossform:::.support_index_preflight(index, evaluation_edges = 0),
    "positive integer"
  )
  expect_error(
    crossform:::.support_index_preflight(index, evaluation_edges = 1.5),
    "positive integer"
  )
})

test_that("compiled searchlights retain their canonical support topology", {
  domain <- volume_domain(array(TRUE, c(5, 4, 3)), id = "compiled-support")
  frame <- compile_frame(
    searchlights(1.5, normalization = "local"), domain
  )

  expect_s3_class(frame$support_index, "effect_support_index")
  expect_identical(frame$support_index$node_ids, frame$index$measurement)
  expect_identical(
    as.matrix(crossform:::.support_index_membership(frame$support_index) != 0),
    as.matrix(frame$weights != 0)
  )
  expect_silent(crossform:::.validate_frame_for_compile(frame))
})

test_that("a 50k-volume topology builds without quadratic storage", {
  if (!identical(Sys.getenv("CROSSFORM_RUN_SCALE_TESTS"), "true")) {
    skip("Set CROSSFORM_RUN_SCALE_TESTS=true to run the 50k topology gate.")
  }
  domain <- volume_domain(array(TRUE, c(50, 40, 25)), id = "volume-50k")
  index <- crossform:::.euclidean_support_index(domain, 1.01)

  expect_identical(index$cost$nodes, 50000)
  expect_lte(index$cost$support_size[["max"]], 7)
  expect_lt(index$cost$pair_pattern_nnz, 30 * domain$n_features)
  expect_lt(index$cost$estimated_structural_bytes, 20 * 1024^2)
  expect_lt(index$cost$pair_pattern_nnz, domain$n_features^2 / 1000)
})

test_that("a 50k-coordinate topology uses vectorized cell lookup", {
  if (!identical(Sys.getenv("CROSSFORM_RUN_SCALE_TESTS"), "true")) {
    skip("Set CROSSFORM_RUN_SCALE_TESTS=true to run the 50k topology gate.")
  }
  coordinates <- as.matrix(expand.grid(
    x = seq_len(37L), y = seq_len(37L), z = seq_len(37L)
  ))
  domain <- abstract_domain(
    nrow(coordinates), coordinates = coordinates, id = "coordinates-50k"
  )
  elapsed <- system.time({
    index <- crossform:::.euclidean_support_index(domain, 1.1)
  })[["elapsed"]]

  expect_identical(index$cost$nodes, 50653)
  expect_identical(index$construction$lookup, "mixed_radix_vectorized")
  expect_lte(index$cost$support_size[["max"]], 7)
  expect_lt(index$cost$pair_pattern_nnz, 30 * domain$n_features)
  expect_lt(index$cost$estimated_structural_bytes, 20 * 1024^2)
  expect_lt(unname(elapsed), 10)
})

test_that("deep validation detects changed supports or pair topology", {
  domain <- volume_domain(array(TRUE, c(3, 3, 2)), id = "support-mutation")
  index <- crossform:::.euclidean_support_index(domain, 1.01)
  index$members[[1L]] <- if (index$members[[1L]] == 1L) 2L else 1L

  expect_error(
    crossform:::.validate_support_index(index, deep = TRUE),
    "strictly increasing|pair pattern|identity"
  )
})

test_that("deep validation rejects noncanonical support order", {
  domain <- volume_domain(array(TRUE, c(3, 3, 2)), id = "support-order")
  index <- crossform:::.euclidean_support_index(domain, 1.01)
  start <- index$ptr[[2L]] + 1
  end <- index$ptr[[3L]]
  index$members[start:end] <- rev(index$members[start:end])

  expect_error(
    crossform:::.validate_support_index(index, deep = TRUE),
    "strictly increasing"
  )
})
