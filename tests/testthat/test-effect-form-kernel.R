rectangular_kernel_fixture <- function(sparse = FALSE) {
  set.seed(20260813)
  left_partitions <- c("e1", "e2")
  right_partitions <- c("r1", "r2")
  left_effects <- c("encode_a", "encode_b")
  right_effects <- c("retrieve_a", "retrieve_b", "retrieve_c")
  left <- stats::setNames(lapply(left_partitions, function(partition) {
    matrix(rnorm(2 * 7), 2, 7)
  }), left_partitions)
  right <- stats::setNames(lapply(right_partitions, function(partition) {
    matrix(rnorm(3 * 7), 3, 7)
  }), right_partitions)
  weights <- matrix(runif(3 * 7), 3, 7)
  if (sparse) weights <- Matrix::Matrix(weights, sparse = TRUE)
  over <- pairing(
    c("e1", "e2"), c("r1", "r2"), weight = c(1, 3), directed = TRUE
  )
  edges <- crossform:::.ordered_partition_edges(
    over, left_partitions, right_partitions, FALSE
  )
  list(
    left = left,
    right = right,
    left_partitions = left_partitions,
    right_partitions = right_partitions,
    left_effects = left_effects,
    right_effects = right_effects,
    frame = additive_frame(weights),
    over = over,
    edges = edges
  )
}

rectangular_kernel_oracle <- function(fixture) {
  weights <- as.matrix(fixture$frame$weights)
  forms <- lapply(seq_len(nrow(weights)), function(measurement) {
    value <- matrix(0, length(fixture$left_effects),
      length(fixture$right_effects))
    for (feature in seq_len(ncol(weights))) {
      atom <- matrix(0, nrow(value), ncol(value))
      for (edge in seq_len(nrow(fixture$edges))) {
        atom <- atom + fixture$edges$weight[[edge]] * outer(
          fixture$left[[fixture$edges$left[[edge]]]][, feature],
          fixture$right[[fixture$edges$right[[edge]]]][, feature]
        )
      }
      value <- value + weights[measurement, feature] * atom
    }
    value
  })
  do.call(rbind, lapply(forms, as.vector))
}

run_rectangular_kernel <- function(fixture, feature_block = 3L,
                                   query = NULL, accumulate_tile = NULL,
                                   read_left = NULL, read_right = NULL) {
  if (is.null(read_left)) {
    read_left <- function(partition, features) {
      fixture$left[[partition]][, features, drop = FALSE]
    }
  }
  if (is.null(read_right)) {
    read_right <- function(partition, features) {
      fixture$right[[partition]][, features, drop = FALSE]
    }
  }
  crossform:::.streamed_effect_form_contraction(
    fixture$frame,
    read_left = read_left,
    read_right = read_right,
    left_partitions = fixture$left_partitions,
    right_partitions = fixture$right_partitions,
    left_effects = fixture$left_effects,
    right_effects = fixture$right_effects,
    ordered_edges = fixture$edges,
    codec = "rectangular",
    query = query,
    feature_block = feature_block,
    row_tile = 2L,
    coordinate_tile = 2L,
    accumulate_tile = accumulate_tile
  )
}

test_that("rectangular streamed forms match an independent dense-loop oracle", {
  fixture <- rectangular_kernel_fixture(sparse = TRUE)
  got <- run_rectangular_kernel(fixture)
  oracle <- rectangular_kernel_oracle(fixture)

  expect_identical(got$codec, "rectangular")
  expect_identical(got$logical_shape, c(2L, 3L))
  expect_equal(got$value, oracle, tolerance = 1e-12)
  expect_identical(got$diagnostics$atom_count, 7L)
  expect_gt(got$diagnostics$max_live_temporary_bytes, 0)
})

test_that("direct pair querying happens in feature tasks and equals late query", {
  fixture <- rectangular_kernel_fixture()
  H <- matrix(c(1, -2, 3, 0, 4, -1), 2, 3)
  physical_query <- matrix(as.vector(H), ncol = 1L)
  complete <- run_rectangular_kernel(fixture, feature_block = 2L)
  direct <- run_rectangular_kernel(
    fixture, feature_block = 2L, query = physical_query
  )
  task <- crossform:::.effect_form_feature_task(
    lapply(fixture$left, function(value) value[, 1:3, drop = FALSE]),
    lapply(fixture$right, function(value) value[, 1:3, drop = FALSE]),
    1:3,
    fixture$left_effects,
    fixture$right_effects,
    fixture$left_partitions,
    fixture$right_partitions,
    fixture$edges,
    codec = "rectangular",
    query = physical_query
  )

  expect_equal(
    direct$value,
    complete$value %*% physical_query,
    tolerance = 1e-12
  )
  expect_true(task$projected)
  expect_equal(task$diagnostics$atom_bytes, 0)
  expect_gt(task$diagnostics$query_atom_bytes, 0)

  feature_oracle <- numeric(3)
  for (feature in 1:3) {
    for (edge in seq_len(nrow(fixture$edges))) {
      left <- fixture$left[[fixture$edges$left[[edge]]]][, feature]
      right <- fixture$right[[fixture$edges$right[[edge]]]][, feature]
      feature_oracle[[feature]] <- feature_oracle[[feature]] +
        fixture$edges$weight[[edge]] * drop(t(left) %*% H %*% right)
    }
  }
  expect_equal(drop(task$atoms), feature_oracle, tolerance = 1e-14)
})

test_that("rectangular transpose and direct-sum laws hold", {
  fixture <- rectangular_kernel_fixture()
  forward <- run_rectangular_kernel(fixture, feature_block = 2L)
  reversed <- fixture
  reversed$left <- fixture$right
  reversed$right <- fixture$left
  reversed$left_partitions <- fixture$right_partitions
  reversed$right_partitions <- fixture$left_partitions
  reversed$left_effects <- fixture$right_effects
  reversed$right_effects <- fixture$left_effects
  reversed$edges <- crossform:::.ordered_partition_edges(
    pairing(
      fixture$over$right,
      fixture$over$left,
      fixture$over$weight,
      directed = TRUE
    ),
    reversed$left_partitions,
    reversed$right_partitions,
    FALSE
  )
  backward <- run_rectangular_kernel(reversed, feature_block = 2L)

  for (measurement in 1:3) {
    forward_form <- matrix(forward$value[measurement, ], 2, 3)
    backward_form <- matrix(backward$value[measurement, ], 3, 2)
    expect_equal(backward_form, t(forward_form), tolerance = 1e-12)
  }

  one <- run_rectangular_kernel(fixture, feature_block = 1L)$value
  whole <- run_rectangular_kernel(fixture, feature_block = 7L)$value
  expect_equal(one, whole, tolerance = 1e-12)
})

test_that("each required side partition is read once per feature block", {
  fixture <- rectangular_kernel_fixture()
  left_reads <- stats::setNames(integer(2), fixture$left_partitions)
  right_reads <- stats::setNames(integer(2), fixture$right_partitions)
  got <- run_rectangular_kernel(
    fixture,
    feature_block = 3L,
    read_left = function(partition, features) {
      left_reads[[partition]] <<- left_reads[[partition]] + 1L
      fixture$left[[partition]][, features, drop = FALSE]
    },
    read_right = function(partition, features) {
      right_reads[[partition]] <<- right_reads[[partition]] + 1L
      fixture$right[[partition]][, features, drop = FALSE]
    }
  )

  expect_identical(unname(left_reads), c(3L, 3L))
  expect_identical(unname(right_reads), c(3L, 3L))
  expect_identical(got$diagnostics$relation_reads, 12L)
})

test_that("rectangular streamed output supports block-backed accumulation", {
  fixture <- rectangular_kernel_fixture(sparse = TRUE)
  oracle <- rectangular_kernel_oracle(fixture)
  path <- tempfile(fileext = ".egm")
  on.exit(unlink(path), add = TRUE)
  store <- crossform:::.file_geometry_store(
    path, dim(oracle), create = TRUE, codec = "rectangular"
  )
  writer <- function(rows, coordinates, increment) {
    existing <- crossform:::.read_geometry_tile(store, rows, coordinates)
    crossform:::.write_geometry_tile(
      store, rows, coordinates, existing + increment
    )
  }
  got <- run_rectangular_kernel(
    fixture, feature_block = 2L, accumulate_tile = writer
  )

  expect_null(got$value)
  expect_equal(
    crossform:::.read_geometry_store(store),
    oracle,
    tolerance = 1e-12
  )
})

self_kernel_fixture <- function(sparse = FALSE) {
  set.seed(918)
  partitions <- c("run1", "run2", "run3")
  effects <- c("face", "house", "object")
  relation <- stats::setNames(lapply(partitions, function(partition) {
    matrix(rnorm(length(effects) * 9), length(effects), 9)
  }), partitions)
  weights <- matrix(c(
    1, 0, 0.5, 0, 0, 0, 0, 0, 0,
    0, 1, 0.5, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0.5, 1, 0, 0,
    0, 0, 0, 0, 0, 0.5, 0, 1, 1
  ), 4, 9, byrow = TRUE)
  if (sparse) weights <- Matrix::Matrix(weights, sparse = TRUE)
  list(
    partitions = partitions,
    effects = effects,
    relation = relation,
    frame = additive_frame(weights),
    over = pairing(
      c("run1", "run1", "run2"),
      c("run2", "run3", "run3"),
      weight = c(1, 2, 3)
    )
  )
}

run_self_kernel <- function(fixture, feature_block) {
  crossform:::.streamed_crossgram_contraction(
    fixture$frame,
    function(partition, features) {
      fixture$relation[[partition]][, features, drop = FALSE]
    },
    fixture$partitions, fixture$effects, fixture$over,
    feature_block = feature_block, row_tile = 2L, coordinate_tile = 2L
  )
}

test_that("legacy self streaming is the packed specialization of ordered atoms", {
  for (sparse in c(FALSE, TRUE)) {
    fixture <- self_kernel_fixture(sparse = sparse)
    edges <- crossform:::.ordered_partition_edges(
      fixture$over, fixture$partitions, fixture$partitions, TRUE
    )
    reader <- function(partition, features) {
      fixture$relation[[partition]][, features, drop = FALSE]
    }
    for (block in c(1L, 4L, 9L)) {
      old <- run_self_kernel(fixture, feature_block = block)$value
      lifted <- crossform:::.streamed_effect_form_contraction(
        fixture$frame,
        read_left = reader,
        left_partitions = fixture$partitions,
        right_partitions = fixture$partitions,
        left_effects = fixture$effects,
        right_effects = fixture$effects,
        ordered_edges = edges,
        codec = "symmetric_packed",
        same_relation = TRUE,
        feature_block = block,
        row_tile = 2L,
        coordinate_tile = 2L
      )$value
      expect_equal(lifted, old, tolerance = 1e-12)
    }
  }
})

test_that("self specialization is invariant to declared partition order", {
  fixture <- self_kernel_fixture(sparse = TRUE)
  reordered <- fixture
  reordered$partitions <- rev(fixture$partitions)
  reordered$relation <- fixture$relation[reordered$partitions]

  expect_equal(
    run_self_kernel(reordered, feature_block = 4L)$value,
    run_self_kernel(fixture, feature_block = 4L)$value,
    tolerance = 1e-12
  )
})

test_that("rectangular kernel memory plans bound named live objects", {
  fixture <- rectangular_kernel_fixture(sparse = TRUE)
  got <- run_rectangular_kernel(fixture, feature_block = 3L)
  plan <- crossform:::.effect_form_kernel_memory_plan(
    fixture$frame,
    fixture$left_effects,
    fixture$right_effects,
    fixture$left_partitions,
    fixture$right_partitions,
    codec = "rectangular",
    feature_block = 3L,
    row_tile = 2L,
    coordinate_tile = 2L,
    storage = "memory"
  )

  expect_gte(
    plan$planned_workspace_bytes,
    got$diagnostics$max_live_temporary_bytes +
      got$diagnostics$durable_output_bytes
  )
  expect_identical(plan$prediction_kind,
    "crossform_owned_workspace_upper_bound")
})

test_that("invalid rectangular kernel semantics fail before readers", {
  fixture <- rectangular_kernel_fixture()
  reads <- 0L
  reader <- function(partition, features) {
    reads <<- reads + 1L
    fixture$left[[partition]][, features, drop = FALSE]
  }

  expect_error(crossform:::.streamed_effect_form_contraction(
    fixture$frame,
    read_left = reader,
    read_right = reader,
    left_partitions = fixture$left_partitions,
    right_partitions = fixture$right_partitions,
    left_effects = fixture$left_effects,
    right_effects = fixture$right_effects,
    ordered_edges = fixture$edges,
    codec = "symmetric_packed"
  ), "self-adjoint", class = "effect_input_error")
  expect_error(run_rectangular_kernel(
    fixture, query = matrix(1, 5, 1), read_left = reader, read_right = reader
  ), "physical form", class = "effect_contract_error")
  expect_identical(reads, 0L)
})
