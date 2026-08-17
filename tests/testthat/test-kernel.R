test_that("row-coordinate-feature tiling agrees with the dense oracle", {
  set.seed(20260812)
  weights <- matrix(runif(7 * 11), 7, 11)
  atoms <- matrix(rnorm(11 * 13), 11, 13)
  oracle <- weights %*% atoms

  tile_shapes <- list(c(1, 1, 1), c(3, 4, 5), c(20, 20, 20), c(2, 7, 3))
  for (shape in tile_shapes) {
    got <- crossform:::.tiled_contraction(
      weights, atoms,
      row_tile = shape[[1]], coordinate_tile = shape[[2]],
      feature_tile = shape[[3]]
    )
    expect_equal(got$value, oracle, tolerance = 1e-12)
    expect_lte(
      got$diagnostics$max_temporary_elements,
      min(shape[[1]], nrow(weights)) * min(shape[[2]], ncol(atoms))
    )
  }
})

test_that("writer-backed contraction never materializes the full output", {
  set.seed(42)
  weights <- matrix(runif(6 * 9), 6, 9)
  atoms <- matrix(rnorm(9 * 10), 9, 10)
  output <- matrix(NA_real_, nrow(weights), ncol(atoms))
  writer <- function(rows, coordinates, value) {
    output[rows, coordinates] <<- value
  }

  got <- crossform:::.tiled_contraction(
    weights, atoms, row_tile = 2, coordinate_tile = 3, feature_tile = 4,
    write_tile = writer
  )

  expect_null(got$value)
  expect_equal(output, weights %*% atoms, tolerance = 1e-12)
  expect_equal(got$diagnostics$max_temporary_elements, 6)
})

test_that("feature permutation preserves the contraction", {
  set.seed(17)
  weights <- matrix(runif(5 * 8), 5, 8)
  atoms <- matrix(rnorm(8 * 6), 8, 6)
  permutation <- sample(seq_len(8))

  original <- crossform:::.tiled_contraction(weights, atoms, 2, 3, 3)$value
  permuted <- crossform:::.tiled_contraction(
    weights[, permutation, drop = FALSE],
    atoms[permutation, , drop = FALSE], 2, 3, 3
  )$value

  expect_equal(permuted, original, tolerance = 1e-12)
})

test_that("tiled contraction rejects dimension and tile errors", {
  expect_error(
    crossform:::.tiled_contraction(matrix(1, 2, 3), matrix(1, 2, 2), 1, 1, 1),
    "feature dimension"
  , class = "effect_input_error")
  expect_error(
    crossform:::.tiled_contraction(matrix(1, 2, 2), matrix(1, 2, 2), 0, 1, 1),
    "positive integer"
  , class = "effect_input_error")
})

crossgram_fixture <- function(sparse = TRUE) {
  set.seed(918)
  partitions <- c("run1", "run2", "run3")
  effects <- c("face", "house", "object")
  relation <- setNames(lapply(partitions, function(x) {
    matrix(rnorm(length(effects) * 9), length(effects), 9)
  }), partitions)
  dense_weights <- matrix(c(
    1, 0, 0.5, 0, 0, 0, 0, 0, 0,
    0, 1, 0.5, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0.5, 1, 0, 0,
    0, 0, 0, 0, 0, 0.5, 0, 1, 1
  ), 4, 9, byrow = TRUE)
  weights <- if (sparse) Matrix::Matrix(dense_weights, sparse = TRUE) else dense_weights
  over <- pairing(
    c("run1", "run1", "run2"),
    c("run2", "run3", "run3"),
    weight = c(1, 2, 3)
  )
  list(
    partitions = partitions,
    effects = effects,
    relation = relation,
    frame = additive_frame(weights),
    over = over
  )
}

independent_crossgram_oracle <- function(fixture) {
  weights <- as.matrix(fixture$frame$weights)
  q <- length(fixture$effects)
  out <- matrix(0, nrow(weights), q * (q + 1) / 2)
  for (measurement in seq_len(nrow(weights))) {
    geometry <- matrix(0, q, q)
    for (feature in seq_len(ncol(weights))) {
      feature_geometry <- matrix(0, q, q)
      for (edge in seq_len(nrow(fixture$over))) {
        left <- fixture$relation[[fixture$over$left[[edge]]]][, feature]
        right <- fixture$relation[[fixture$over$right[[edge]]]][, feature]
        cross <- outer(left, right)
        feature_geometry <- feature_geometry + fixture$over$weight[[edge]] *
          0.5 * (cross + t(cross))
      }
      geometry <- geometry + weights[measurement, feature] * feature_geometry
    }
    packed <- numeric(ncol(out))
    k <- 0L
    for (column in seq_len(q)) {
      for (row in column:q) {
        k <- k + 1L
        packed[[k]] <- geometry[row, column] *
          if (row == column) 1 else sqrt(2)
      }
    }
    out[measurement, ] <- packed
  }
  out
}

run_streamed_fixture <- function(fixture, feature_block = 3L,
                                 row_tile = 2L, coordinate_tile = 2L,
                                 reader = NULL, accumulate_tile = NULL,
                                 retain_local_relations = FALSE) {
  if (is.null(reader)) {
    reader <- function(partition, features) {
      fixture$relation[[partition]][, features, drop = FALSE]
    }
  }
  crossform:::.streamed_crossgram_contraction(
    fixture$frame, reader, fixture$partitions, fixture$effects, fixture$over,
    feature_block = feature_block, row_tile = row_tile,
    coordinate_tile = coordinate_tile, accumulate_tile = accumulate_tile,
    retain_local_relations = retain_local_relations
  )
}

test_that("sparse streamed cross-Gram agrees with an independent relation oracle", {
  fixture <- crossgram_fixture(sparse = TRUE)
  got <- run_streamed_fixture(fixture)
  oracle <- independent_crossgram_oracle(fixture)

  expect_equal(got$value, oracle, tolerance = 1e-12)
  expect_identical(got$diagnostics$atom_count, 9L)
  expect_identical(got$diagnostics$measurement_kind,
    "static-owned-buffer-accounting")
})

test_that("sparse and dense frame representations are numerically equivalent", {
  sparse <- crossgram_fixture(sparse = TRUE)
  dense <- sparse
  dense$frame <- additive_frame(as.matrix(sparse$frame$weights))

  expect_equal(
    run_streamed_fixture(sparse, feature_block = 4)$value,
    run_streamed_fixture(dense, feature_block = 4)$value,
    tolerance = 1e-13
  )
})

test_that("each relation partition is read once per canonical feature block", {
  fixture <- crossgram_fixture()
  counts <- setNames(integer(length(fixture$partitions)), fixture$partitions)
  reader <- function(partition, features) {
    counts[[partition]] <<- counts[[partition]] + 1L
    fixture$relation[[partition]][, features, drop = FALSE]
  }
  got <- run_streamed_fixture(fixture, feature_block = 4, reader = reader)

  expect_identical(unname(counts), rep(3L, 3))
  expect_identical(got$diagnostics$feature_blocks, 3L)
  expect_identical(got$diagnostics$relation_reads, 9L)
})

test_that("the same relation reads can retain bounded coherent marginals", {
  fixture <- crossgram_fixture()
  got <- run_streamed_fixture(fixture, feature_block = 4,
    retain_local_relations = TRUE)
  weights <- as.matrix(fixture$frame$weights)

  for (partition in fixture$partitions) {
    expect_equal(
      unname(got$local_relations[, , partition]),
      weights %*% t(fixture$relation[[partition]]),
      tolerance = 1e-13
    )
  }
  expect_gt(got$diagnostics$durable_local_relation_bytes, 0)
  expect_gt(got$diagnostics$max_local_product_bytes, 0)
  expect_gt(got$diagnostics$max_local_existing_bytes, 0)
  expect_gt(got$diagnostics$max_local_replacement_bytes, 0)
})

test_that("coherent-only streaming omits atom formation and total output", {
  fixture <- crossgram_fixture()
  got <- crossform:::.streamed_crossgram_contraction(
    fixture$frame,
    function(partition, features) {
      fixture$relation[[partition]][, features, drop = FALSE]
    },
    fixture$partitions, fixture$effects, fixture$over,
    feature_block = 3, retain_local_relations = TRUE, form_total = FALSE
  )

  expect_null(got$value)
  expect_equal(got$diagnostics$atom_count, 0)
  expect_equal(got$diagnostics$max_atom_bytes, 0)
  expect_gt(got$diagnostics$durable_local_relation_bytes, 0)
})

test_that("coherent geometry from retained relations matches an explicit oracle", {
  fixture <- crossgram_fixture()
  streamed <- run_streamed_fixture(fixture, feature_block = 4,
    retain_local_relations = TRUE)
  mass <- Matrix::rowSums(fixture$frame$weights)
  got <- crossform:::.coherent_geometry_from_local(
    streamed$local_relations, fixture$over, mass, row_tile = 2
  )
  q <- length(fixture$effects)
  oracle <- matrix(0, nrow(fixture$frame$weights), q * (q + 1) / 2)
  for (measurement in seq_len(nrow(oracle))) {
    geometry <- matrix(0, q, q)
    for (edge in seq_len(nrow(fixture$over))) {
      left <- streamed$local_relations[measurement, , fixture$over$left[[edge]]]
      right <- streamed$local_relations[measurement, , fixture$over$right[[edge]]]
      cross <- outer(left, right)
      geometry <- geometry + fixture$over$weight[[edge]] *
        0.5 * (cross + t(cross)) / mass[[measurement]]
    }
    oracle[measurement, ] <- oracle_svec(geometry)
  }

  expect_equal(got$value, oracle, tolerance = 1e-12)
  expect_gt(got$diagnostics$max_tile_bytes, 0)
  expect_gt(got$diagnostics$max_work_bytes, 0)
})

test_that("canonical feature-block choices agree under the numeric contract", {
  fixture <- crossgram_fixture()
  oracle <- run_streamed_fixture(fixture, feature_block = 1)$value
  for (block in c(2L, 4L, 9L)) {
    got <- run_streamed_fixture(fixture, feature_block = block)$value
    agreement <- numerical_agreement(
      oracle, got, guarantee = "block_partition",
      contract = numerical_contract(atol = 1e-12, rtol = 1e-11)
    )
    expect_true(agreement$passed, info = paste("feature block", block))
  }
})

test_that("external tiled accumulation does not require an in-memory output", {
  fixture <- crossgram_fixture()
  accumulated <- matrix(0, nrow(fixture$frame$weights), 6)
  writer <- function(rows, coordinates, increment) {
    accumulated[rows, coordinates] <<-
      accumulated[rows, coordinates, drop = FALSE] + increment
  }
  got <- run_streamed_fixture(fixture, feature_block = 2,
    accumulate_tile = writer)

  expect_null(got$value)
  expect_false(got$diagnostics$external_accumulator_memory_measured)
  expect_equal(accumulated, independent_crossgram_oracle(fixture), tolerance = 1e-12)
})

test_that("streaming diagnostics statically bound every owned buffer category", {
  fixture <- crossgram_fixture()
  got <- run_streamed_fixture(fixture, feature_block = 3,
    row_tile = 2, coordinate_tile = 2)
  required <- c(
    "max_relation_bytes", "max_atom_bytes", "max_atom_work_bytes",
    "max_weight_slice_bytes", "max_atom_slice_bytes", "max_product_bytes",
    "max_existing_slice_bytes", "max_replacement_bytes"
  )
  expect_true(all(vapply(got$diagnostics[required], function(x) x > 0,
    logical(1))))
  expect_gte(
    got$diagnostics$max_live_temporary_bytes,
    sum(vapply(got$diagnostics[required], identity, numeric(1)))
  )
  expect_gt(got$diagnostics$durable_output_bytes, 0)
})

test_that("streamed cross-Gram fails closed on malformed relation blocks", {
  fixture <- crossgram_fixture()
  bad_reader <- function(partition, features) matrix(1, 1, length(features))
  expect_error(run_streamed_fixture(fixture, reader = bad_reader),
    "invalid effect-by-feature", class = "effect_input_error")

  forged_pairing <- fixture$over
  forged_pairing$weight[[1]] <- -1
  expect_error(crossform:::.streamed_crossgram_contraction(
    fixture$frame,
    function(partition, features) fixture$relation[[partition]][, features, drop = FALSE],
    fixture$partitions, fixture$effects, forged_pairing
  ), "nonnegative", class = "effect_input_error")
})

test_that("task observers record a failed task boundary", {
  fixture <- crossgram_fixture()
  events <- character()
  features <- list()
  observer <- function(event, feature_ids) {
    events <<- c(events, event)
    features[[length(features) + 1L]] <<- feature_ids
  }
  bad_reader <- function(partition, feature_ids) {
    if (feature_ids[[1L]] > 3L) stop("block failed")
    fixture$relation[[partition]][, feature_ids, drop = FALSE]
  }

  expect_error(crossform:::.streamed_crossgram_contraction(
    fixture$frame, bad_reader, fixture$partitions, fixture$effects,
    fixture$over, feature_block = 3, task_observer = observer), "block failed")
  expect_identical(events, c("started", "completed", "started", "failed"))
  expect_identical(features[[4L]], 4:6)
})
