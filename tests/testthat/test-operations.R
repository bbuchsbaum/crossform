# Contract tests for R/operations.R: the declarative operation-stage
# vocabulary (normalizers, edge transforms, partition reducers) and the plan
# that fixes their order and identity.

test_that("edge normalizers are declarations with an explicit zero policy", {
  plain <- inner_product()

  expect_s3_class(plain, "effect_edge_normalizer")
  expect_identical(names(plain), c("schema_version", "kind", "zero_policy"))
  expect_identical(plain$schema_version, 1L)
  expect_identical(plain$kind, "inner_product")
  # An inner product has no denominator, so there is no zero case to resolve.
  expect_identical(plain$zero_policy, "error")

  expect_identical(crossform:::covariance()$kind, "covariance")
  expect_identical(crossform:::covariance()$zero_policy, "error")
  expect_identical(crossform:::cosine()$kind, "cosine")
  expect_identical(crossform:::cosine("zero")$zero_policy, "zero")
  expect_identical(crossform:::correlation()$zero_policy, "error")
  expect_identical(crossform:::correlation("zero")$zero_policy, "zero")
  expect_error(crossform:::cosine("clip"), "'arg' should be one of")
  expect_error(crossform:::correlation("clip"), "'arg' should be one of")

  # The declaration computes nothing: identical declarations are identical
  # values, so they can key a plan signature.
  expect_identical(inner_product(), inner_product())
  expect_false(identical(inner_product(), crossform:::covariance()))
})

test_that("normalizer validation rejects every noncanonical shape", {
  valid <- inner_product()
  expect_identical(crossform:::.validate_edge_normalizer(valid), valid)

  reject <- function(mutate) {
    value <- valid
    value <- mutate(value)
    expect_error(crossform:::.validate_edge_normalizer(value),
      "Invalid edge normalizer")
  }
  reject(function(x) unclass(x))
  reject(function(x) {
    x$kind <- "spearman"
    x
  })
  reject(function(x) {
    x$zero_policy <- "clip"
    x
  })
  reject(function(x) {
    x$schema_version <- 2L
    x
  })
  reject(function(x) {
    structure(x[c("kind", "zero_policy", "schema_version")],
      class = "effect_edge_normalizer")
  })
  expect_error(
    crossform:::.validate_edge_normalizer(
      structure(list(schema_version = 1L, kind = "inner_product"),
        class = "effect_edge_normalizer")
    ),
    "Invalid edge normalizer"
  )
})

test_that("edge transforms carry exactly the policy their kind requires", {
  identity_transform <- crossform:::.identity_edge_transform()
  expect_identical(identity_transform$kind, "identity")
  expect_null(identity_transform$boundary)
  expect_null(identity_transform$delta)
  expect_null(identity_transform$ties)
  expect_identical(
    crossform:::.validate_edge_transform(identity_transform),
    identity_transform
  )

  strict <- fisher_z("error")
  expect_identical(strict$kind, "fisher_z")
  expect_identical(strict$boundary, "error")
  expect_null(strict$delta)
  clipped <- fisher_z("clip", delta = 1e-6)
  expect_identical(clipped$boundary, "clip")
  expect_identical(clipped$delta, 1e-6)

  ranked <- rank_edges()
  expect_identical(ranked$kind, "rank_edges")
  expect_identical(ranked$ties, "average")
  expect_identical(rank_edges("min")$ties, "min")
  expect_error(rank_edges("random"), "'arg' should be one of")

  # A boundary policy and its parameter must agree: `clip` needs a delta in
  # (0, 1), and `error` must not carry one.
  expect_error(fisher_z("clip"), "Invalid Fisher boundary policy")
  expect_error(fisher_z("clip", delta = 0), "Invalid Fisher boundary policy")
  expect_error(fisher_z("clip", delta = 1), "Invalid Fisher boundary policy")
  expect_error(fisher_z("error", delta = 0.1),
    "Invalid Fisher boundary policy")
  expect_error(fisher_z("shrink"), "'arg' should be one of")

  # Cross-kind parameters are refused rather than ignored.
  expect_error(
    crossform:::.validate_edge_transform(
      crossform:::.new_edge_transform("identity", boundary = "error")
    ),
    "cannot carry policy parameters"
  )
  expect_error(
    crossform:::.validate_edge_transform(
      crossform:::.new_edge_transform("fisher_z", boundary = "error",
        ties = "average")
    ),
    "Invalid Fisher boundary policy"
  )
  expect_error(
    crossform:::.validate_edge_transform(
      crossform:::.new_edge_transform("rank_edges")
    ),
    "Invalid edge-ranking tie policy"
  )
  expect_error(
    crossform:::.validate_edge_transform(
      crossform:::.new_edge_transform("winsorize")
    ),
    "Invalid edge transform"
  )
})

test_that("partition reducers name two distinct estimands", {
  edge_first <- reduce_partitions()
  aggregated <- aggregate_first()

  expect_s3_class(edge_first, "effect_partition_reducer")
  expect_identical(edge_first$kind, "weighted_sum")
  expect_identical(edge_first$weight_convention, "normalized_unit_mass")
  expect_identical(edge_first$order, "edge_first")
  expect_identical(aggregated$order, "aggregate_first")
  expect_false(identical(edge_first, aggregated))

  # `aggregate_first()` is the documented default of the measurement
  # pipeline, and `reduce_partitions()` is not.
  expect_identical(eval(formals(measurement_form)$reducer), aggregated)
  expect_identical(eval(formals(coupling)$reducer), aggregated)
  expect_false(identical(eval(formals(measurement_form)$reducer), edge_first))
})

test_that("the edge operation plan fixes the stage order and its lowering", {
  plan <- crossform:::.edge_operation_plan()

  expect_s3_class(plan, "effect_edge_operation_plan")
  expect_identical(plan$operation_order, c(
    "ordered_edge_product", "spatial_normalization", "edge_transform",
    "partition_reduction", "pair_query"
  ))
  expect_identical(plan$schema_version, 1L)
  expect_identical(plan$lowering, "bilinear")
  expect_identical(plan$normalizer, unclass(inner_product()))
  expect_identical(plan$transform,
    unclass(crossform:::.identity_edge_transform()))
  expect_identical(plan$reducer, unclass(reduce_partitions()))
  expect_match(plan$signature, "^sha256:[0-9a-f]{64}$")

  lowering <- function(normalizer, transform = NULL) {
    crossform:::.edge_operation_plan(normalizer, transform)$lowering
  }
  expect_identical(lowering(inner_product()), "bilinear")
  expect_identical(lowering(crossform:::covariance()), "centered_bilinear")
  expect_identical(lowering(crossform:::cosine()), "normalized_bilinear")
  expect_identical(lowering(crossform:::correlation()), "normalized_bilinear")
  # Any nonlinear transform forces the edges to be materialized, whatever the
  # normalizer would otherwise allow.
  expect_identical(
    lowering(crossform:::correlation(), fisher_z("error")),
    "required_edge_materialization"
  )
  expect_identical(
    lowering(inner_product(), rank_edges()),
    "required_edge_materialization"
  )

  # Fisher transformation is defined on correlations, so a non-correlation
  # normalizer is refused at plan time.
  expect_error(
    crossform:::.edge_operation_plan(inner_product(), fisher_z("error")),
    "requires correlation-valued edge input"
  )
  expect_error(
    crossform:::.edge_operation_plan(crossform:::cosine(), fisher_z("error")),
    "requires correlation-valued edge input"
  )
  # The one admitted pairing records both declarations verbatim.
  admitted <- crossform:::.edge_operation_plan(
    crossform:::correlation("zero"), fisher_z("clip", 1e-6)
  )
  expect_identical(admitted$normalizer, unclass(crossform:::correlation("zero")))
  expect_identical(admitted$transform, unclass(fisher_z("clip", 1e-6)))
  expect_identical(admitted$lowering, "required_edge_materialization")
})

test_that("every operation component enters the plan signature", {
  base <- crossform:::.edge_operation_plan()
  expect_identical(base$signature, crossform:::.edge_operation_plan()$signature)

  variants <- list(
    normalizer = crossform:::.edge_operation_plan(crossform:::covariance()),
    transform = crossform:::.edge_operation_plan(
      crossform:::correlation(), fisher_z("clip", 1e-6)
    ),
    delta = crossform:::.edge_operation_plan(
      crossform:::correlation(), fisher_z("clip", 1e-3)
    ),
    ties = crossform:::.edge_operation_plan(inner_product(), rank_edges("min")),
    reducer = crossform:::.edge_operation_plan(reducer = aggregate_first())
  )
  signatures <- c(base = base$signature,
    vapply(variants, `[[`, character(1), "signature"))
  expect_identical(anyDuplicated(signatures), 0L)

  expect_error(crossform:::.edge_operation_plan(reducer = "aggregate_first"),
    "reducer")
})

test_that("applying an edge transform equals its declared operation", {
  values <- matrix(c(0.3, -0.4, 0.8, 0.05), 2L, 2L)

  # Identity is exactly the identity, not a numerically-equal copy.
  expect_identical(
    crossform:::.apply_edge_transform(
      values, crossform:::.identity_edge_transform()
    ),
    values
  )

  expect_equal(
    crossform:::.apply_edge_transform(values, fisher_z("error")),
    atanh(values), tolerance = 1e-12
  )
  boundary <- matrix(c(1, -1, 0.5, 0.2), 2L, 2L)
  expect_error(
    crossform:::.apply_edge_transform(boundary, fisher_z("error")),
    "absolute-correlation boundary"
  )
  expect_equal(
    crossform:::.apply_edge_transform(boundary, fisher_z("clip", 1e-8)),
    atanh(pmin(1 - 1e-8, pmax(-1 + 1e-8, boundary))),
    tolerance = 1e-12
  )
  expect_true(all(is.finite(
    crossform:::.apply_edge_transform(boundary, fisher_z("clip", 1e-8))
  )))

  # Ranking is over the whole edge coordinate set, in column-major order,
  # and preserves the matrix shape.
  tied <- matrix(c(3, 1, 2, 1), 2L, 2L)
  expect_identical(
    crossform:::.apply_edge_transform(tied, rank_edges("average")),
    matrix(c(4, 1.5, 3, 1.5), 2L, 2L)
  )
  # The integer-valued tie methods inherit base R's integer storage mode;
  # only "average" can produce a fractional rank.
  expect_identical(
    crossform:::.apply_edge_transform(tied, rank_edges("min")),
    matrix(c(4L, 1L, 3L, 1L), 2L, 2L)
  )
  expect_identical(
    crossform:::.apply_edge_transform(tied, rank_edges("max")),
    matrix(c(4L, 2L, 3L, 2L), 2L, 2L)
  )
  expect_identical(
    crossform:::.apply_edge_transform(tied, rank_edges("first")),
    matrix(c(4L, 1L, 3L, 2L), 2L, 2L)
  )
  expect_identical(
    crossform:::.apply_edge_transform(tied, rank_edges("last")),
    matrix(c(4L, 2L, 3L, 1L), 2L, 2L)
  )
})
