# Cost-shape regression guards for the streamed first-moment pass.
#
# `contrast_energy()` needs the coherent component and the signed marginals, so
# it asks the streaming kernel to retain first moments while `rdm()` does not.
# That extra requirement used to cost more than the whole 28-column RDM
# contraction, because the kernel contracted the frame once per partition and
# then wrote the result into a three-way array, duplicating the durable array
# on every measurement tile. These tests pin the shape of the fixed pass:
# one product per partition family, one frame pass, and the same number of
# source reads whether or not first moments are retained.

first_moment_cost_fixture <- function() {
  set.seed(20260816)
  partitions <- c("run1", "run2", "run3")
  effects <- c("face", "house", "object")
  features <- 12L
  measurements <- 8L
  relation <- stats::setNames(lapply(partitions, function(partition) {
    matrix(stats::rnorm(length(effects) * features),
      length(effects), features)
  }), partitions)
  weights <- matrix(0, measurements, features)
  for (measurement in seq_len(measurements)) {
    support <- ((measurement - 1L) + seq_len(4L) - 1L) %% features + 1L
    weights[measurement, support] <- seq_along(support)
  }
  over <- pairing(
    c("run1", "run1", "run2"), c("run2", "run3", "run3"),
    weight = c(1, 2, 3)
  )
  list(
    partitions = partitions,
    effects = effects,
    relation = relation,
    frame = additive_frame(Matrix::Matrix(weights, sparse = TRUE)),
    edges = crossform:::.ordered_partition_edges(
      over, partitions, partitions, same_relation = TRUE
    )
  )
}

run_first_moment_cost <- function(fixture, retain, reads = NULL,
                                  row_tile = 2L, feature_block = 4L) {
  crossform:::.streamed_effect_form_contraction(
    frame = fixture$frame,
    read_left = function(partition, features) {
      if (!is.null(reads)) reads$count <- reads$count + 1L
      fixture$relation[[partition]][, features, drop = FALSE]
    },
    left_partitions = fixture$partitions,
    right_partitions = fixture$partitions,
    left_effects = fixture$effects,
    right_effects = fixture$effects,
    ordered_edges = fixture$edges,
    codec = "symmetric_packed",
    same_relation = TRUE,
    feature_block = feature_block,
    row_tile = row_tile,
    coordinate_tile = 2L,
    retain_first_moments = retain,
    form_total = TRUE
  )
}

test_that("retaining first moments reads every source block exactly once", {
  fixture <- first_moment_cost_fixture()
  plain_reads <- new.env(parent = emptyenv())
  plain_reads$count <- 0L
  retained_reads <- new.env(parent = emptyenv())
  retained_reads$count <- 0L

  plain <- run_first_moment_cost(fixture, FALSE, plain_reads)
  retained <- run_first_moment_cost(fixture, TRUE, retained_reads)

  expected_reads <- 3L * length(fixture$partitions)
  expect_identical(plain_reads$count, expected_reads)
  expect_identical(retained_reads$count, expected_reads)
  expect_identical(
    retained$diagnostics$relation_reads, plain$diagnostics$relation_reads
  )
  expect_identical(
    retained$diagnostics$feature_blocks, plain$diagnostics$feature_blocks
  )
  # The retained pass is a strict addition to the same reads: the total form
  # it also produces is unchanged.
  expect_identical(retained$value, plain$value)
})

test_that("first moments contract one partition family per product", {
  fixture <- first_moment_cost_fixture()
  retained <- run_first_moment_cost(fixture, TRUE)
  q <- length(fixture$effects)
  partitions <- length(fixture$partitions)
  row_tile <- 2L

  # One fused product per side, not one per partition: the live product is a
  # measurement tile by the whole (effect within partition) panel.
  expect_identical(
    retained$diagnostics$max_left_first_product_bytes,
    8 * row_tile * q * partitions
  )
  expect_identical(
    retained$diagnostics$max_left_first_replacement_bytes,
    retained$diagnostics$max_left_first_product_bytes
  )
  # The panel the frame multiplies is one transposed copy of the relation
  # blocks already accounted for, and it is inside the reported live bound.
  expect_gte(
    retained$diagnostics$max_live_temporary_bytes,
    retained$diagnostics$max_relation_bytes * 2 +
      retained$diagnostics$max_weight_slice_bytes +
      3 * retained$diagnostics$max_left_first_product_bytes
  )
})

test_that("the frame block a measurement tile reads is owned and accounted", {
  fixture <- first_moment_cost_fixture()

  # More than one measurement tile reads each feature block, so the block's
  # columns are taken once and sliced; the copy is reported.
  tiled <- run_first_moment_cost(fixture, TRUE, row_tile = 2L)
  expect_gt(tiled$diagnostics$max_frame_block_bytes, 0)
  expect_gte(
    tiled$diagnostics$max_live_temporary_bytes,
    tiled$diagnostics$max_frame_block_bytes
  )

  # One measurement tile reads it once, so no block is owned at all.
  single <- run_first_moment_cost(fixture, TRUE, row_tile = 1024L)
  expect_identical(single$diagnostics$max_frame_block_bytes, 0)

  # Slicing route is an implementation detail: both tilings agree exactly.
  expect_identical(single$value, tiled$value)
  expect_identical(single$first_moments$left, tiled$first_moments$left)
})

test_that("a contrast view costs no extra source pass over its plan", {
  set.seed(20260816)
  domain <- abstract_domain(24L, id = "first-moment-cost")
  effects <- c("face", "house", "object")
  partitions <- c("run1", "run2", "run3")
  sources <- stats::setNames(lapply(partitions, function(partition) {
    value <- matrix(stats::rnorm(length(effects) * 24L), length(effects), 24L)
    rownames(value) <- effects
    value
  }), partitions)
  relation <- relation(sources, domain = domain)
  frame <- compile_frame(
    regions(rep(c("a", "b", "c", "d"), each = 6L)), domain
  )
  plan <- plan_geometry(relation, frame, cross_partitions(
    relation, independence = "independent", generalizes_over = "run"
  ))
  weights <- c(face = 1, house = -1, object = 0)

  fused <- contrast_energy(plan, weights)
  projected <- contrast_energy(materialize_geometry(plan), weights)

  expect_equal(fused$total, projected$total, tolerance = 1e-12)
  expect_equal(fused$coherent, projected$coherent, tolerance = 1e-12)
  expect_equal(fused$signed, projected$signed, tolerance = 1e-12)
})
