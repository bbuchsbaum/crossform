# Regression: coupling tolerances are relative to the data's own magnitude.
#
# Every normalized coupling view is invariant to a global rescaling of the
# data (correlations, canonical correlations, CKA alignment, and Gaussian
# information are unit-free), and covariance blocks scale by the square of
# the factor. With absolute tolerances (1e-10 in the data's units) data
# rescaled by 1e-6 had every self-variance treated as zero, and data rescaled
# by 1e4 (raw BOLD magnitudes) had ordinary round-off asymmetry refused as a
# structurally asymmetric self-block.

scale_fixture_form <- function(scale, scalar = FALSE, partition = NULL,
                               flat = FALSE) {
  native <- abstract_domain(4, id = "scale:native:v1")
  times <- effect_space(paste0("time", 1:6), basis_id = "scale:time:v1")
  trend <- seq(-1.5, 1.5, length.out = 6)
  set.seed(20261011)
  noise <- list(matrix(rnorm(24), 6), matrix(rnorm(24), 6))
  session <- function(shift, i) scale * (cbind(
    trend + shift,
    -0.8 * trend + c(0.2, -0.1, 0.1, 0, -0.1, -0.1),
    sin(seq(shift, pi, length.out = 6)),
    cos(seq(shift, pi, length.out = 6))
  ) + noise[[i]])
  values <- list(session1 = session(0, 1), session2 = session(0.1, 2))
  if (flat) values <- lapply(values, function(value) {
    value[, 4L] <- 0
    value
  })
  if (!is.null(partition)) values <- values[partition]
  signals <- relation(values, effects = times, domain = native)
  operators <- list(
    anterior = matrix(c(1, 0.3, 0, 0), 1),
    posterior = matrix(c(0, 1, 0.7, 0), 1)
  )
  if (flat) operators$posterior <- matrix(c(0, 0, 0, 1), 1)
  if (!scalar) {
    operators$lateral <- rbind(c(0.2, 0, 1, 0), c(0, 0.1, 0, 1))
  }
  nodes <- measurement_frame(operators, domain = native,
    id = if (scalar) "scale:scalar-nodes:v1" else "scale:mixed-nodes:v1")
  pairs <- expand.grid(
    from = nodes$node_ids, to = nodes$node_ids, stringsAsFactors = FALSE
  )
  covariance <- variation_query(
    (diag(6) - matrix(1 / 6, 6, 6)) / 5, times, "time", "joint_covariance"
  )
  over <- pairing(
    signals$partitions, signals$partitions, directed = TRUE,
    self_pairs = "allow_biased", independence = "not_independent"
  )
  measurement_form(
    signals, edge_frame(pairs$from, pairs$to, nodes), covariance, over
  )
}

scale_views <- function(scale) {
  mixed <- scale_fixture_form(scale)
  scalar <- scale_fixture_form(scale, scalar = TRUE)
  # The ridge is a declared quantity in the data's units, so it is rescaled
  # with the data; everything downstream is then exactly scale-free.
  ridge <- 0.05 * scale^2
  model <- gaussian_covariance_model(list(assumption = "joint Gaussian"))
  list(
    correlation = connectivity(scalar, "correlation")$values$correlation,
    canonical = connectivity(mixed, "canonical",
      ridge = ridge)$values$canonical_correlation,
    geometry = connectivity(mixed,
      "geometry_alignment")$values$geometry_alignment,
    information = connectivity(mixed, "gaussian_information",
      ridge = ridge, model = model)$values$information,
    covariance = covariance_coupling(mixed)$values,
    partitioned = crossform:::.partitioned_pearson_coupling(
      list(
        scale_fixture_form(scale, scalar = TRUE, partition = "session1"),
        scale_fixture_form(scale, scalar = TRUE, partition = "session2")
      ),
      placement = "after_partition_aggregation"
    )$values$value
  )
}

test_that("normalized coupling views are invariant to the data's units", {
  reference <- scale_views(1)
  # The reference is a genuine, non-degenerate answer.
  expect_gt(max(reference$canonical), 0.5)
  expect_true(all(reference$information > 0))
  expect_true(any(abs(reference$correlation) < 0.999))

  for (scale in c(1e-6, 1e4)) {
    rescaled <- scale_views(scale)
    for (view in c("correlation", "canonical", "geometry", "information",
                   "partitioned")) {
      expect_equal(rescaled[[view]], reference[[view]], tolerance = 1e-8,
        info = sprintf("view %s at scale %g", view, scale))
    }
    # Covariance blocks are equivariant: they scale by the square.
    expect_equal(
      lapply(rescaled$covariance, function(block) block / scale^2),
      reference$covariance, tolerance = 1e-8,
      info = sprintf("covariance at scale %g", scale)
    )
  }
})

test_that("pseudo-inverse and self-block cutoffs are relative", {
  self_block <- matrix(c(2, 0.3, 0.3, 1.5), 2L, 2L)
  reference <- crossform:::.decomposition_inverse_sqrt(self_block, 1e-10)
  for (scale in c(1e-14, 1e-6, 1e4, 1e12)) {
    expect_equal(
      crossform:::.decomposition_inverse_sqrt(scale * self_block, 1e-10),
      reference / sqrt(scale), tolerance = 1e-10,
      info = sprintf("scale %g", scale)
    )
  }
  # A genuinely rank-deficient block still drops its null direction at any
  # scale, and an all-zero block whitens to zero rather than dividing by it.
  singular <- tcrossprod(c(1, 2))
  for (scale in c(1e-12, 1, 1e8)) {
    whitened <- crossform:::.decomposition_inverse_sqrt(
      scale * singular, 1e-10
    )
    expect_equal(qr(whitened)$rank, 1L)
  }
  expect_identical(
    crossform:::.decomposition_inverse_sqrt(matrix(0, 2, 2), 1e-10),
    matrix(0, 2, 2)
  )
  # Indefiniteness is also judged relative to the block's magnitude.
  expect_error(
    crossform:::.decomposition_inverse_sqrt(1e-12 * diag(c(1, -1)), 1e-10),
    "positive self-blocks"
  )
  expect_identical(crossform:::.relative_tolerance_floor(1e-10, 0),
    1e-10 * .Machine$double.xmin)
  expect_identical(crossform:::.relative_tolerance_floor(1e-10, c(-4, 2)),
    4e-10)
})

test_that("a node with no variation is refused at every data scale", {
  for (scale in c(1e-6, 1, 1e4)) {
    form <- scale_fixture_form(scale, scalar = TRUE, flat = TRUE)
    refusal <- catch_refusal(connectivity(form, "correlation"))
    expect_s3_class(refusal, "effect_capability_refusal")
    expect_identical(refusal$reasons, "self_variance_not_strictly_positive")
    expect_match(conditionMessage(refusal), "relative tolerance")
    expect_error(connectivity(form, "geometry_alignment"),
      "nonzero self-geometries")
  }
})
