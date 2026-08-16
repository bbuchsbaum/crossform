test_that("measurement legs bind source, output, basis, units, and support", {
  domain <- abstract_domain(3, feature_ids = c("v1", "v2", "v3"),
    id = "neural:v1")
  output <- crossform:::.measurement_axis(
    c("mode1", "mode2"), "roi:a:v1", basis_id = "gradients:v1",
    units = c("z", "z"), provenance = list(source = "fixed-atlas")
  )
  operator <- matrix(c(1, 0, 2, 0, -1, 3), 2, 3, byrow = TRUE)
  leg <- crossform:::.measurement_leg(
    operator, domain, output,
    provenance = list(estimator = "fixed-linear")
  )

  expect_s3_class(leg, "effect_measurement_leg")
  expect_identical(leg$source_domain, domain$reference)
  expect_identical(leg$output_space, output)
  expect_identical(leg$support, domain$feature_ids)
  expect_identical(leg$output_space$basis_id, "gradients:v1")
  expect_identical(leg$output_space$units, c("z", "z"))
  expect_silent(crossform:::.validate_measurement_leg(leg, domain))

  forged <- leg
  forged$support <- c("v1", "v2")
  expect_error(crossform:::.validate_measurement_leg(forged), "support")
})

test_that("equal metrics do not erase oriented measurement identity", {
  domain <- abstract_domain(3, id = "neural:v1")
  first_axis <- crossform:::.measurement_axis(
    c("a", "b"), "roi:v1", basis_id = "basis:original"
  )
  rotated_axis <- crossform:::.measurement_axis(
    c("u", "v"), "roi:v1", basis_id = "basis:rotated"
  )
  first_operator <- matrix(c(1, 2, 0, 0, 1, 3), 2, 3, byrow = TRUE)
  rotation <- matrix(c(0, -1, 1, 0), 2, 2)
  first <- crossform:::.measurement_leg(
    first_operator, domain, first_axis
  )
  rotated <- crossform:::.measurement_leg(
    rotation %*% first_operator, domain, rotated_axis
  )

  expect_equal(
    crossform:::.measurement_leg_metric(first),
    crossform:::.measurement_leg_metric(rotated),
    tolerance = 1e-14
  )
  expect_false(identical(first$signature, rotated$signature))
  q <- matrix(c(2, -1, 0, 3, 1, 4, -2, 5, 6), 3)
  expect_false(isTRUE(all.equal(
    first$operator %*% q %*% t(first$operator),
    first$operator %*% q %*% t(rotated$operator),
    tolerance = 1e-14
  )))
})

test_that("learned legs require frozen portable training provenance", {
  domain <- abstract_domain(2)
  output <- crossform:::.measurement_axis("mode", "learned:mode")
  training <- paste0("sha256:", paste(rep("a", 64), collapse = ""))

  expect_error(crossform:::.measurement_leg(
    matrix(c(1, 2), 1), domain, output, estimation = "learned_frozen"
  ), "frozen training provenance")
  learned <- crossform:::.measurement_leg(
    matrix(c(1, 2), 1), domain, output, estimation = "learned_frozen",
    provenance = list(
      frozen = TRUE,
      training_signature = training,
      split = "independent-training"
    )
  )
  expect_identical(learned$estimation, "learned_frozen")
  expect_silent(crossform:::.validate_measurement_leg(learned))
  expect_error(crossform:::.measurement_leg(
    matrix(c(1, 2), 1), domain, output, estimation = "learned_frozen",
    provenance = list(frozen = TRUE, training_signature = identity)
  ), "nonportable")
})

test_that("direct sums preserve named component identities and ranges", {
  domain <- abstract_domain(3, id = "neural:v1")
  coherent <- crossform:::.measurement_leg(
    matrix(c(1, 1, 1) / sqrt(3), 1), domain,
    crossform:::.measurement_axis("mean", "coherent:v1")
  )
  configuration <- crossform:::.measurement_leg(
    matrix(c(1, -1, 0, 1, 1, -2), 2, 3, byrow = TRUE), domain,
    crossform:::.measurement_axis(
      c("gradient1", "gradient2"), "configuration:v1",
      basis_id = "fixed-gradients:v1"
    )
  )
  combined <- crossform:::.direct_sum_measurement_leg(
    list(coherent = coherent, configuration = configuration),
    id = "roi:decomposition:v1"
  )

  expect_identical(combined$decomposition$components,
    c("coherent", "configuration"))
  expect_identical(combined$decomposition$ranges$coherent, 1L)
  expect_identical(combined$decomposition$ranges$configuration, 2:3)
  expect_identical(combined$operator,
    rbind(coherent$operator, configuration$operator))
  expect_identical(combined$output_space$coordinates,
    c("coherent::mean", "configuration::gradient1",
      "configuration::gradient2"))
})

test_that("measurement frames keep ordered legs and conservative dual metadata", {
  domain <- abstract_domain(3, id = "neural:v1")
  make_leg <- function(index) {
    operator <- matrix(0, 1, 3)
    operator[[index]] <- 1
    crossform:::.measurement_leg(
      operator, domain,
      crossform:::.measurement_axis(
        paste0("v", index), paste0("voxel:", index)
      )
    )
  }
  legs <- stats::setNames(lapply(1:3, make_leg), paste0("node", 1:3))
  frame <- crossform:::.measurement_frame(
    legs, injectivity = "positive_diagonal_coverage"
  )

  expect_identical(frame$node_ids, names(legs))
  expect_equal(frame$frame_operator, diag(3), tolerance = 0)
  expect_true(frame$coverage$complete)
  expect_true(frame$injectivity$guaranteed)
  expect_true(frame$dual$eligible)
  expect_false(frame$dual$materialized)
  expect_identical(frame$dual$stability, "not_certified")

  rank_deficient <- crossform:::.measurement_frame(legs[1:2])
  expect_false(rank_deficient$injectivity$guaranteed)
  expect_false(rank_deficient$dual$eligible)
  expect_error(crossform:::.measurement_frame(
    legs[1:2], injectivity = "positive_diagonal_coverage"
  ), "strictly positive diagonal")
})

test_that("measurement edges contain only explicit ordered requests", {
  domain <- abstract_domain(3)
  make_leg <- function(index) {
    operator <- matrix(0, 1, 3)
    operator[[index]] <- 1
    crossform:::.measurement_leg(
      operator, domain,
      crossform:::.measurement_axis(
        paste0("v", index), paste0("node-axis:", index)
      )
    )
  }
  frame <- crossform:::.measurement_frame(stats::setNames(
    lapply(1:3, make_leg), c("a", "b", "c")
  ))
  edges <- crossform:::.measurement_edges(
    c("a", "a"), c("b", "c"), frame, weight = c(1, -0.5)
  )

  expect_identical(nrow(edges$edges), 2L)
  expect_identical(edges$edges$left, c("a", "a"))
  expect_identical(edges$edges$right, c("b", "c"))
  expect_false("edges" %in% names(frame))
  reversed <- crossform:::.reverse_measurement_edges(edges, frame)
  expect_identical(reversed$edges$left, edges$edges$right)
  expect_identical(reversed$edges$right, edges$edges$left)
  expect_identical(
    crossform:::.reverse_measurement_edges(reversed, frame), edges
  )
  expect_error(crossform:::.measurement_edges(
    c("a", "a"), c("b", "b"), frame
  ), "unique explicit")
})

test_that("square-root additive legs exactly recover every local geometry", {
  set.seed(2026081203)
  domain <- abstract_domain(5, feature_ids = paste0("v", 1:5),
    id = "neural:v1")
  weights <- matrix(runif(4 * 5, 0.1, 2), 4, 5)
  relation <- matrix(rnorm(3 * 5), 3, 5)

  for (normalization in c("none", "local", "conservative")) {
    normalized <- switch(normalization,
      none = weights,
      local = weights / rowSums(weights),
      conservative = sweep(weights, 2L, colSums(weights), "/")
    )
    frame <- additive_frame(normalized, normalization, domain = domain)
    measured <- crossform:::.measurement_frame_from_additive(frame, "total")
    expect_true(measured$injectivity$guaranteed)
    for (node in seq_len(nrow(normalized))) {
      leg <- measured$legs[[node]]
      got <- relation %*% t(leg$operator) %*%
        (leg$operator %*% t(relation))
      expected <- relation %*% diag(normalized[node, ]) %*% t(relation)
      expect_equal(got, expected, tolerance = 2e-13,
        info = paste("normalization", normalization, "node", node))
    }
  }
})

test_that("coherent common-mode legs recover the certified coherent component", {
  set.seed(2026081204)
  domain <- abstract_domain(4, id = "neural:v1")
  weights <- matrix(runif(3 * 4, 0.2, 2), 3, 4)
  relation <- matrix(rnorm(2 * 4), 2, 4)
  frame <- additive_frame(weights, domain = domain)
  coherent <- crossform:::.measurement_frame_from_additive(frame, "coherent")

  expect_false(coherent$injectivity$guaranteed)
  for (node in seq_len(nrow(weights))) {
    leg <- coherent$legs[[node]]
    got <- relation %*% t(leg$operator) %*%
      (leg$operator %*% t(relation))
    first <- drop(relation %*% weights[node, ])
    expected <- tcrossprod(first) / sum(weights[node, ])
    expect_equal(got, expected, tolerance = 2e-13)
  }
})

test_that("additive adapters are canonical across dense and sparse storage", {
  domain <- abstract_domain(4, id = "neural:v1")
  weights <- matrix(c(
    1, 0, 2, 0,
    0, 3, 1, 4
  ), 2, 4, byrow = TRUE)
  dense <- additive_frame(weights, domain = domain)
  sparse <- additive_frame(Matrix::Matrix(weights, sparse = TRUE),
    domain = domain)
  dense_frame <- crossform:::.measurement_frame_from_additive(dense)
  sparse_frame <- crossform:::.measurement_frame_from_additive(sparse)

  expect_equal(dense_frame$frame_operator, sparse_frame$frame_operator,
    tolerance = 0)
  expect_identical(
    lapply(dense_frame$legs, `[[`, "operator"),
    lapply(sparse_frame$legs, `[[`, "operator")
  )
})

test_that("additive sparse identity is canonical without densification", {
  domain <- abstract_domain(4, id = "neural:sparse-signature")
  base <- Matrix::sparseMatrix(
    i = c(1L, 1L, 2L), j = c(1L, 3L, 2L), x = c(1, 2, 3),
    dims = c(2L, 4L)
  )
  stored_zero <- base
  stored_zero[2, 4] <- 1
  stored_zero@x[stored_zero@i == 1L &
    rep(seq_len(ncol(stored_zero)), diff(stored_zero@p)) == 4L] <- 0
  canonical <- Matrix::drop0(stored_zero)
  dense <- additive_frame(as.matrix(canonical), domain = domain)
  sparse_zero <- additive_frame(stored_zero, domain = domain)
  sparse <- additive_frame(canonical, domain = domain)

  expect_identical(
    crossform:::.additive_frame_signature(sparse_zero),
    crossform:::.additive_frame_signature(sparse)
  )
  expect_identical(
    crossform:::.additive_frame_signature(dense),
    crossform:::.additive_frame_signature(sparse)
  )
  expect_s4_class(
    crossform:::.canonical_additive_weights(stored_zero),
    "CsparseMatrix"
  )
  changed <- canonical
  changed[1, 1] <- changed[1, 1] + 1e-12
  expect_false(identical(
    crossform:::.additive_frame_signature(sparse),
    crossform:::.additive_frame_signature(
      additive_frame(changed, domain = domain)
    )
  ))
})
