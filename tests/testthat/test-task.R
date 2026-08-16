# Contract tests for R/task.R: the universal numerical primitive
# `.effect_form_feature_task()`, its input validator, and the physical query
# operator decoder. Every numeric assertion below is against a hand-written
# loop over ordered edges, not against another package path.

task_fixture <- local({
  built <- NULL
  function() {
    if (!is.null(built)) {
      return(built)
    }
    set.seed(20260820L)
    q_left <- 3L
    q_right <- 2L
    features <- 4L
    left_effects <- c("a", "b", "c")
    right_effects <- c("x", "y")
    make <- function(q) matrix(rnorm(q * features), q, features)
    left <- list(run1 = make(q_left), run2 = make(q_left))
    right <- list(s1 = make(q_right), s2 = make(q_right))
    self_edges <- crossform:::.ordered_partition_edges(
      pairing("run1", "run2", 1, directed = FALSE,
        independence = "independent"),
      names(left), names(left), TRUE
    )
    cross_edges <- crossform:::.ordered_partition_edges(
      pairing(c("run1", "run2"), c("s2", "s1"), c(1, 3), directed = TRUE,
        independence = "independent"),
      names(left), names(right), FALSE
    )
    declared_edges <- crossform:::.ordered_partition_edges(
      pairing(c("run1", "run2"), c("run2", "run1"), c(1, 1),
        directed = TRUE, independence = "independent"),
      names(left), names(left), TRUE
    )
    built <<- list(
      q_left = q_left, q_right = q_right, features = features,
      feature_ids = seq_len(features),
      left_effects = left_effects, right_effects = right_effects,
      left = left, right = right,
      self_edges = self_edges, cross_edges = cross_edges,
      declared_edges = declared_edges
    )
    built
  }
})

# An explicit loop over ordered edges: sum_e w_e * L_e' Op R_e per feature.
task_edge_oracle <- function(left, right, edges, feature, transform) {
  Reduce(`+`, lapply(seq_len(nrow(edges)), function(edge) {
    edges$weight[[edge]] * transform(
      left[[edges$left[[edge]]]][, feature],
      right[[edges$right[[edge]]]][, feature]
    )
  }))
}

test_that("an unqueried self task packs every symmetric form coordinate", {
  fixture <- task_fixture()
  got <- crossform:::.effect_form_feature_task(
    fixture$left, fixture$left, fixture$feature_ids,
    fixture$left_effects, fixture$left_effects,
    ordered_edges = fixture$self_edges,
    codec = "symmetric_packed", same_relation = TRUE
  )
  oracle <- t(vapply(fixture$feature_ids, function(feature) {
    oracle_svec(task_edge_oracle(
      fixture$left, fixture$left, fixture$self_edges, feature, outer
    ))
  }, numeric(6L)))

  expect_s3_class(got, "effect_form_feature_task_result")
  expect_identical(got$codec, "symmetric_packed")
  # q(q + 1)/2 is computed in double arithmetic, unlike the rectangular
  # width, which stays integer.
  expect_identical(got$physical_width, 6)
  expect_identical(got$logical_shape, c(3L, 3L))
  expect_false(got$projected)
  expect_true(got$atoms_formed)
  expect_identical(dim(got$atoms), c(4L, 6L))
  expect_equal(got$atoms, oracle, tolerance = 1e-12)
})

test_that("an unqueried rectangular task packs column-major vec coordinates", {
  fixture <- task_fixture()
  got <- crossform:::.effect_form_feature_task(
    fixture$left, fixture$right, fixture$feature_ids,
    fixture$left_effects, fixture$right_effects,
    ordered_edges = fixture$cross_edges,
    codec = "rectangular"
  )
  oracle <- t(vapply(fixture$feature_ids, function(feature) {
    as.vector(task_edge_oracle(
      fixture$left, fixture$right, fixture$cross_edges, feature, outer
    ))
  }, numeric(6L)))

  expect_identical(got$codec, "rectangular")
  expect_identical(got$physical_width, 6L)
  expect_identical(got$logical_shape, c(3L, 2L))
  expect_identical(dim(got$atoms), c(4L, 6L))
  expect_equal(got$atoms, oracle, tolerance = 1e-12)
  # The partition weights are the ones the pairing normalized, not the raw
  # `c(1, 3)` the caller wrote.
  expect_equal(sum(fixture$cross_edges$weight), 1, tolerance = 1e-12)
  expect_equal(fixture$cross_edges$weight, c(0.25, 0.75), tolerance = 1e-12)
})

test_that("a dense query is contracted directly rather than through atoms", {
  fixture <- task_fixture()
  query <- matrix(rnorm(6L * 2L), 6L, 2L)
  got <- crossform:::.effect_form_feature_task(
    fixture$left, fixture$right, fixture$feature_ids,
    fixture$left_effects, fixture$right_effects,
    ordered_edges = fixture$cross_edges,
    codec = "rectangular", query = query
  )
  operators <- crossform:::.physical_query_operators(query, 3L, 2L,
    "rectangular")
  oracle <- vapply(seq_len(2L), function(view) {
    vapply(fixture$feature_ids, function(feature) {
      sum(operators[[view]] * task_edge_oracle(
        fixture$left, fixture$right, fixture$cross_edges, feature, outer
      ))
    }, numeric(1))
  }, numeric(4L))

  expect_true(got$projected)
  expect_identical(dim(got$atoms), c(4L, 2L))
  expect_equal(got$atoms, oracle, tolerance = 1e-12)
  # The atoms are booked as query atoms, and the complete-form budget is zero.
  expect_identical(got$diagnostics$atom_bytes, 0)
  expect_identical(got$diagnostics$query_atom_bytes, 8 * 4 * 2)
  expect_identical(got$diagnostics$relation_bytes,
    8 * 4 * (2 * 3 + 2 * 2))
})

test_that("a packed query decodes to the symmetric operator it names", {
  fixture <- task_fixture()
  operator <- matrix(c(2, -1, 0.5, -1, 3, 1, 0.5, 1, -2), 3L, 3L)
  packed <- matrix(oracle_svec(operator), ncol = 1L)
  decoded <- crossform:::.physical_query_operators(packed, 3L, 3L,
    "symmetric_packed")

  expect_length(decoded, 1L)
  expect_equal(decoded[[1L]], operator, tolerance = 1e-12)
  expect_identical(
    crossform:::.physical_query_operators(
      matrix(as.numeric(1:6), 6L, 1L), 3L, 2L, "rectangular"
    )[[1L]],
    matrix(as.numeric(1:6), 3L, 2L)
  )

  got <- crossform:::.effect_form_feature_task(
    fixture$left, fixture$left, fixture$feature_ids,
    fixture$left_effects, fixture$left_effects,
    ordered_edges = fixture$self_edges,
    codec = "symmetric_packed", query = packed, same_relation = TRUE
  )
  oracle <- matrix(vapply(fixture$feature_ids, function(feature) {
    sum(operator * task_edge_oracle(
      fixture$left, fixture$left, fixture$self_edges, feature, outer
    ))
  }, numeric(1)), ncol = 1L)

  expect_equal(got$atoms, oracle, tolerance = 1e-12)
})

test_that("a structured pair-difference query tiles without a dense operator", {
  fixture <- task_fixture()
  query <- crossform:::.pair_difference_query(fixture$left_effects)
  got <- crossform:::.effect_form_feature_task(
    fixture$left, fixture$left, fixture$feature_ids,
    fixture$left_effects, fixture$left_effects,
    ordered_edges = fixture$self_edges,
    codec = "symmetric_packed", query = query, same_relation = TRUE
  )
  oracle <- vapply(seq_along(query$pair_left), function(pair) {
    direction <- numeric(3L)
    direction[[query$pair_left[[pair]]]] <- 1
    direction[[query$pair_right[[pair]]]] <- -1
    vapply(fixture$feature_ids, function(feature) {
      sum(tcrossprod(direction) * task_edge_oracle(
        fixture$left, fixture$left, fixture$self_edges, feature, outer
      ))
    }, numeric(1))
  }, numeric(4L))

  expect_identical(dim(got$atoms), c(4L, 3L))
  expect_equal(got$atoms, oracle, tolerance = 1e-12)

  # A pair-space coefficient map reduces the same edges to fewer views.
  coefficients <- matrix(c(1, -1, 0, 0, 0.5, 0.5), 2L, 3L, byrow = TRUE)
  reduced <- crossform:::.effect_form_feature_task(
    fixture$left, fixture$left, fixture$feature_ids,
    fixture$left_effects, fixture$left_effects,
    ordered_edges = fixture$self_edges,
    codec = "symmetric_packed", same_relation = TRUE,
    query = crossform:::.pair_difference_query(
      fixture$left_effects, coefficients = coefficients
    )
  )
  expect_identical(dim(reduced$atoms), c(4L, 2L))
  expect_equal(reduced$atoms, oracle %*% t(coefficients), tolerance = 1e-12)

  # The bounded pair tile is what the work budget accounts for; the full
  # pair-by-feature matrix is never a live allocation.
  expect_identical(reduced$diagnostics$max_atom_work_bytes,
    8 * 6 * min(64L, 3L) * 4)
  expect_lt(reduced$diagnostics$max_atom_work_bytes,
    8 * 6 * 64L * 4 + 1)
})

test_that("form_atoms = FALSE compiles the plan without forming any atom", {
  fixture <- task_fixture()
  got <- crossform:::.effect_form_feature_task(
    fixture$left, fixture$right, fixture$feature_ids,
    fixture$left_effects, fixture$right_effects,
    ordered_edges = fixture$cross_edges,
    codec = "rectangular", form_atoms = FALSE
  )

  expect_null(got$atoms)
  expect_false(got$atoms_formed)
  expect_false(got$projected)
  expect_identical(got$diagnostics$atom_bytes, 0)
  expect_identical(got$diagnostics$max_atom_work_bytes, 0)
  expect_identical(got$physical_width, 6L)
  expect_error(
    crossform:::.effect_form_feature_task(
      fixture$left, fixture$right, fixture$feature_ids,
      fixture$left_effects, fixture$right_effects,
      ordered_edges = fixture$cross_edges, codec = "rectangular",
      form_atoms = NA
    ),
    "`form_atoms` must be TRUE or FALSE"
  )
})

test_that("task inputs are validated before any product is formed", {
  fixture <- task_fixture()
  build <- function(...) {
    arguments <- list(
      left_relations = fixture$left, right_relations = fixture$right,
      feature_ids = fixture$feature_ids,
      left_effects = fixture$left_effects,
      right_effects = fixture$right_effects,
      ordered_edges = fixture$cross_edges, codec = "rectangular"
    )
    overrides <- list(...)
    # Replace whole arguments; `modifyList()` would recurse into the
    # relation families instead of substituting them.
    arguments[names(overrides)] <- overrides
    do.call(crossform:::.effect_form_feature_task, arguments)
  }

  expect_error(build(feature_ids = c(2, 1)),
    "strictly increasing positive integers")
  expect_error(build(feature_ids = c(0, 1)),
    "strictly increasing positive integers")
  expect_error(build(feature_ids = c(1, 1)),
    "strictly increasing positive integers")
  expect_error(build(feature_ids = integer()),
    "strictly increasing positive integers")

  expect_error(build(left_relations = unname(fixture$left)),
    "`left_relations` must match one named partition family")
  expect_error(
    build(right_partitions = "s1"),
    "`right_relations` must match one named partition family"
  )
  expect_error(
    build(left_partitions = c("run1", "run1")),
    "`left_relations` must match one named partition family"
  )
  # A consistently subset family is fine on its own; it is the ordered edges
  # that then no longer name a declared partition.
  expect_error(
    build(right_relations = fixture$right["s1"]),
    "Ordered partition edges are missing or noncanonical"
  )
  broken <- fixture$left
  broken$run1[1L, 1L] <- NA_real_
  expect_error(build(left_relations = broken),
    "Every left relation input must be a finite effect-by-feature matrix")
  wrong_width <- fixture$right
  wrong_width$s1 <- wrong_width$s1[, -1L, drop = FALSE]
  expect_error(build(right_relations = wrong_width),
    "Every right relation input must be a finite effect-by-feature matrix")

  expect_error(build(codec = "packed"), "'arg' should be one of")
  # A packed atom is only meaningful for a self-adjoint self-form task.
  expect_error(build(codec = "symmetric_packed"),
    "Symmetric-packed atoms require a self-adjoint self-form task")
  # A half-edge expansion is only admissible on a declared self relation.
  expect_error(
    crossform:::.effect_form_feature_task(
      fixture$left, fixture$left, fixture$feature_ids,
      fixture$left_effects, fixture$left_effects,
      ordered_edges = fixture$self_edges, codec = "symmetric_packed",
      same_relation = FALSE
    ),
    "Self-adjoint ordered-edge expansion is inconsistent"
  )

  expect_error(build(query = matrix(1, 5L, 1L)),
    "`query` must match the finite physical form coordinates")
  expect_error(build(query = matrix(NA_real_, 6L, 1L)),
    "`query` must match the finite physical form coordinates")
  expect_error(build(query = matrix(0, 6L, 0L)),
    "`query` must match the finite physical form coordinates")
  expect_error(
    build(query = crossform:::.pair_difference_query(fixture$left_effects)),
    "self-form queries"
  )
})

test_that("the same-relation flag only changes the accounted read volume", {
  fixture <- task_fixture()
  shared <- crossform:::.effect_form_feature_task(
    fixture$left, fixture$left, fixture$feature_ids,
    fixture$left_effects, fixture$left_effects,
    ordered_edges = fixture$declared_edges,
    codec = "rectangular", same_relation = TRUE
  )
  duplicated <- crossform:::.effect_form_feature_task(
    fixture$left, fixture$left, fixture$feature_ids,
    fixture$left_effects, fixture$left_effects,
    ordered_edges = fixture$declared_edges,
    codec = "rectangular", same_relation = FALSE
  )

  oracle <- t(vapply(fixture$feature_ids, function(feature) {
    as.vector(task_edge_oracle(
      fixture$left, fixture$left, fixture$declared_edges, feature, outer
    ))
  }, numeric(9L)))
  expect_equal(shared$atoms, oracle, tolerance = 1e-12)
  expect_equal(shared$atoms, duplicated$atoms, tolerance = 0)
  expect_gt(max(abs(shared$atoms)), 0)
  expect_identical(shared$diagnostics$relation_bytes, 8 * 4 * 2 * 3)
  expect_identical(duplicated$diagnostics$relation_bytes,
    2 * shared$diagnostics$relation_bytes)
})
