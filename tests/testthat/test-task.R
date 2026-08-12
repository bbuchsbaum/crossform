task_fixture <- function() {
  set.seed(731)
  partitions <- c("run1", "run2", "run3")
  effects <- c("face", "house", "object")
  relations <- setNames(lapply(partitions, function(partition) {
    matrix(rnorm(length(effects) * 7), nrow = length(effects),
      dimnames = list(effects, NULL))
  }), partitions)
  over <- pairing(
    left = c("run1", "run1", "run2"),
    right = c("run2", "run3", "run3"),
    weight = c(0.2, 0.3, 0.5)
  )
  weights <- Matrix::Matrix(matrix(c(
    1, 1, 0, 0, 0, 0, 0,
    0, 0.5, 1, 0.5, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 1
  ), nrow = 3, byrow = TRUE), sparse = TRUE)
  list(
    partitions = partitions,
    effects = effects,
    relations = relations,
    over = over,
    frame = additive_frame(weights)
  )
}

explicit_task_atoms <- function(fixture, feature_ids) {
  q <- length(fixture$effects)
  atoms <- matrix(0, length(feature_ids), q * (q + 1) / 2)
  for (feature_offset in seq_along(feature_ids)) {
    feature <- feature_ids[[feature_offset]]
    value <- matrix(0, q, q)
    for (edge in seq_len(nrow(fixture$over))) {
      left <- fixture$relations[[fixture$over$left[[edge]]]][, feature]
      right <- fixture$relations[[fixture$over$right[[edge]]]][, feature]
      cross <- outer(left, right)
      value <- value + fixture$over$weight[[edge]] *
        0.5 * (cross + t(cross))
    }
    atoms[feature_offset, ] <- effectagram:::.svec_symmetric(value)
  }
  atoms
}

make_task <- function(fixture, feature_ids, query = NULL) {
  relations <- lapply(fixture$relations, function(value) {
    value[, feature_ids, drop = FALSE]
  })
  effectagram:::.crossgram_feature_task(
    relations = relations,
    feature_ids = feature_ids,
    effects = fixture$effects,
    partitions = fixture$partitions,
    over = fixture$over,
    query = query
  )
}

test_that("the canonical feature task equals an explicit Gram oracle", {
  fixture <- task_fixture()
  task <- make_task(fixture, 2:5)

  expect_s3_class(task, "effect_feature_task_result")
  expect_equal(task$atoms, explicit_task_atoms(fixture, 2:5), tolerance = 1e-13)
  expect_identical(task$feature_ids, 2:5)
  expect_identical(task$partitions, fixture$partitions)
  expect_identical(task$effects, fixture$effects)
})

test_that("feature-task payloads are bounded and executor agnostic", {
  fixture <- task_fixture()
  task <- make_task(fixture, 1:3)

  expect_named(task, c(
    "feature_ids", "partitions", "effects", "relations", "atoms",
    "projected", "atoms_formed", "packed_width", "diagnostics"
  ))
  expect_false(any(vapply(task, is.function, logical(1))))
  expect_false(any(c("frame", "source", "reader", "output", "executor") %in%
    names(task)))
  component_bytes <- sum(vapply(
    c(task$relations, list(task$atoms)),
    function(value) as.double(utils::object.size(value)), numeric(1)
  ))
  expect_lt(
    as.double(utils::object.size(task)),
    component_bytes + 20000
  )
})

test_that("coherent-only feature tasks omit cross-Gram atoms", {
  fixture <- task_fixture()
  task <- effectagram:::.crossgram_feature_task(
    fixture$relations, 1:7, fixture$effects, fixture$partitions,
    fixture$over, form_atoms = FALSE
  )

  expect_null(task$atoms)
  expect_false(task$atoms_formed)
  expect_equal(task$diagnostics$atom_bytes, 0)
  expect_equal(task$diagnostics$max_atom_work_bytes, 0)
})

test_that("query projection is part of the same pure feature task", {
  fixture <- task_fixture()
  full <- make_task(fixture, 1:4)
  query <- matrix(seq_len(ncol(full$atoms) * 2) / 17,
    nrow = ncol(full$atoms), ncol = 2)
  projected <- make_task(fixture, 1:4, query = query)

  expect_true(projected$projected)
  expect_equal(projected$atoms, full$atoms %*% query, tolerance = 1e-13)
  expect_gt(projected$diagnostics$query_atom_bytes, 0)
  expect_equal(projected$diagnostics$atom_bytes,
    as.double(utils::object.size(full$atoms)))
})

test_that("ordered reduction is invariant to task completion order", {
  fixture <- task_fixture()
  tasks <- list(
    make_task(fixture, 1:2),
    make_task(fixture, 3:5),
    make_task(fixture, 6:7)
  )
  ordered <- effectagram:::.reduce_crossgram_tasks(
    tasks, fixture$frame, fixture$partitions, fixture$effects,
    output_width = ncol(tasks[[1]]$atoms), row_tile = 2,
    coordinate_tile = 2, retain_local_relations = TRUE
  )
  completed_out_of_order <- effectagram:::.reduce_crossgram_tasks(
    tasks[c(3, 1, 2)], fixture$frame, fixture$partitions, fixture$effects,
    output_width = ncol(tasks[[1]]$atoms), row_tile = 2,
    coordinate_tile = 2, retain_local_relations = TRUE
  )
  streamed <- effectagram:::.streamed_crossgram_contraction(
    fixture$frame,
    function(partition, features) {
      fixture$relations[[partition]][, features, drop = FALSE]
    },
    fixture$partitions, fixture$effects, fixture$over,
    feature_block = 2, row_tile = 2, coordinate_tile = 2,
    retain_local_relations = TRUE
  )

  expect_equal(completed_out_of_order$value, ordered$value, tolerance = 0)
  expect_equal(completed_out_of_order$local_relations,
    ordered$local_relations, tolerance = 0)
  expect_equal(ordered$value, streamed$value, tolerance = 1e-13)
  expect_equal(ordered$local_relations, streamed$local_relations,
    tolerance = 1e-13)
})

test_that("the canonical reducer rejects gaps, overlaps, and wrong identities", {
  fixture <- task_fixture()
  task_1 <- make_task(fixture, 1:2)
  task_2 <- make_task(fixture, 4:7)
  expect_error(
    effectagram:::.reduce_crossgram_tasks(
      list(task_1, task_2), fixture$frame, fixture$partitions,
      fixture$effects, output_width = ncol(task_1$atoms)
    ),
    "canonical contiguous order"
  )

  overlap <- make_task(fixture, 2:7)
  expect_error(
    effectagram:::.reduce_crossgram_tasks(
      list(task_1, overlap), fixture$frame, fixture$partitions,
      fixture$effects, output_width = ncol(task_1$atoms)
    ),
    "canonical contiguous order"
  )

  wrong <- task_1
  wrong$effects <- rev(wrong$effects)
  expect_error(
    effectagram:::.reduce_crossgram_tasks(
      list(wrong, make_task(fixture, 3:7)), fixture$frame,
      fixture$partitions, fixture$effects,
      output_width = ncol(task_1$atoms)
    ),
    "identities"
  )
})

test_that("feature-task construction fails before producing malformed payloads", {
  fixture <- task_fixture()
  bad <- fixture$relations
  bad$run1 <- bad$run1[, 1:2, drop = FALSE]
  expect_error(
    effectagram:::.crossgram_feature_task(
      bad, 1:3, fixture$effects, fixture$partitions, fixture$over
    ),
    "effect-by-feature"
  )

  expect_error(make_task(fixture, c(1, 1)), "strictly increasing")
  expect_error(make_task(fixture, c(2, 1)), "strictly increasing")
})
