task_relation <- function(partitions, effects, domain, reads = NULL,
                          revision_character = "a") {
  effect_count <- if (inherits(effects, "effect_space")) {
    length(effects$coordinates)
  } else {
    length(effects)
  }
  if (is.null(reads)) {
    sources <- stats::setNames(lapply(seq_along(partitions), function(index) {
      matrix(index, effect_count, domain$reference$n_features)
    }), partitions)
    return(relation(sources, effects = effects, domain = domain))
  }
  source <- function(features) {
    reads$count <- reads$count + 1L
    matrix(1, effect_count, length(features))
  }
  revision <- paste0(
    "sha256:", paste(rep(revision_character, 64), collapse = "")
  )
  relation(
    stats::setNames(rep(list(source), length(partitions)), partitions),
    source_dims = rep(list(c(effect_count, domain$reference$n_features)),
      length(partitions)),
    effects = effects,
    domain = domain,
    capabilities = source_capabilities(TRUE, stable_revision = revision)
  )
}

test_that("undirected partition pairs compile to canonical ordered half edges", {
  domain <- abstract_domain(3, id = "shared:v1")
  rel <- task_relation(c("r1", "r2", "r3"), c("a", "b"), domain)
  over <- cross_partitions(rel)
  task <- crossform:::.compile_effect_task(rel, over)
  edges <- task$ordered_edges

  expect_s3_class(task, "effect_compiled_task")
  expect_identical(nrow(edges), 2L * nrow(over))
  expect_equal(sum(edges$weight), 1, tolerance = 1e-15)
  expect_identical(attr(edges, "expansion"), "self_adjoint_half_edges")
  for (position in seq(1L, nrow(edges), by = 2L)) {
    reverse <- position + 1L
    expect_identical(edges$left[[position]], edges$right[[reverse]])
    expect_identical(edges$right[[position]], edges$left[[reverse]])
    expect_identical(edges$weight[[position]], edges$weight[[reverse]])
    expect_identical(edges$weight[[position]],
      over$weight[[edges$input_edge[[position]]]] / 2)
  }

  flipped <- pairing(over$right, over$left, over$weight, directed = FALSE)
  flipped_task <- crossform:::.compile_effect_task(rel, flipped)
  expect_identical(flipped_task$ordered_edges, task$ordered_edges)
  expect_identical(flipped_task$task_id, task$task_id)
})

test_that("ordered half edges exactly retain the old self-adjoint estimator", {
  domain <- abstract_domain(2, id = "shared:v1")
  matrices <- list(
    r1 = matrix(c(1, 2, 3, -1), 2),
    r2 = matrix(c(2, 0, -2, 4), 2),
    r3 = matrix(c(5, -1, 1, 3), 2)
  )
  rel <- relation(matrices, effects = c("a", "b"), domain = domain)
  over <- cross_partitions(rel)
  task <- crossform:::.compile_effect_task(rel, over)

  old <- matrix(0, 2, 2)
  for (edge in seq_len(nrow(over))) {
    cross <- matrices[[over$left[[edge]]]] %*%
      t(matrices[[over$right[[edge]]]])
    old <- old + over$weight[[edge]] * 0.5 * (cross + t(cross))
  }
  ordered <- matrix(0, 2, 2)
  for (edge in seq_len(nrow(task$ordered_edges))) {
    ordered <- ordered + task$ordered_edges$weight[[edge]] *
      matrices[[task$ordered_edges$left[[edge]]]] %*%
        t(matrices[[task$ordered_edges$right[[edge]]]])
  }

  expect_equal(ordered, old, tolerance = 1e-15)
})

test_that("two homogeneous relation sides compile with structural direction", {
  domain <- abstract_domain(4, id = "shared:v1")
  left_space <- effect_space(c("encode_a", "encode_b"), basis_id = "encode:v1")
  right_space <- effect_space(c("retrieve_a", "retrieve_b", "retrieve_c"),
    basis_id = "retrieve:v1")
  left <- task_relation(c("e1", "e2"), left_space, domain)
  right <- task_relation(c("r1", "r2"), right_space, domain)
  over <- pairing(c("e1", "e2"), c("r1", "r2"), c(1, 3), directed = TRUE)
  H <- matrix(1:6, 2, 3)
  task <- crossform:::.compile_effect_task(
    left, over, right, pair_query(H, left_space, right_space)
  )

  expect_false(task$same_relation)
  expect_identical(task$left_space, left_space)
  expect_identical(task$right_space, right_space)
  expect_identical(task$ordered_edges$orientation, c("declared", "declared"))
  expect_identical(task$reducer$kind, "weighted_sum")
  expect_identical(task$bridge$kind, "identity")
  expect_identical(task$bridge$common_space, domain$reference)

  reversed <- crossform:::.compile_effect_task(
    right,
    pairing(over$right, over$left, over$weight, directed = TRUE),
    left,
    pair_query(t(H), right_space, left_space)
  )
  expect_false(identical(reversed$task_id, task$task_id))
  expect_identical(reversed$left_space, task$right_space)
  expect_identical(reversed$right_space, task$left_space)
  expect_identical(reversed$ordered_edges$left, task$ordered_edges$right)
  expect_identical(reversed$ordered_edges$right, task$ordered_edges$left)
})

test_that("task semantic failures occur before lazy relation reads", {
  reads <- new.env(parent = emptyenv())
  reads$count <- 0L
  domain <- abstract_domain(4, id = "shared:v1")
  changed_domain <- abstract_domain(4, feature_ids = 4:1, id = "shared:v1")
  left_space <- effect_space(c("l1", "l2"), basis_id = "left:v1")
  right_space <- effect_space(c("r1", "r2", "r3"), basis_id = "right:v1")
  left <- task_relation(c("e1", "e2"), left_space, domain, reads, "a")
  right <- task_relation(c("r1", "r2"), right_space, domain, reads, "b")
  wrong_domain <- task_relation(
    c("r1", "r2"), right_space, changed_domain, reads, "c"
  )
  over <- pairing("e1", "r1", directed = TRUE)

  expect_error(crossform:::.compile_effect_task(
    left, over, right,
    pair_query(matrix(1, 2, 3),
      effect_space(c("l1", "l2"), basis_id = "wrong"), right_space)
  ), "axes")
  expect_error(crossform:::.compile_effect_task(
    left, pairing("missing", "r1", directed = TRUE), right
  ), "endpoints")
  expect_error(crossform:::.compile_effect_task(left, over, wrong_domain),
    "Distinct neural spaces")
  expect_error(crossform:::.compile_effect_task(
    left, pairing("e1", "r1", directed = FALSE), right
  ), "undirected compatibility")
  expect_identical(reads$count, 0L)
})

test_that("task identity covers axes edges bridge reducer and query", {
  domain <- abstract_domain(3, id = "shared:v1")
  other_domain <- abstract_domain(3, id = "shared:v2")
  left_space <- effect_space(c("l1", "l2"), basis_id = "left:v1")
  right_space <- effect_space(c("r1", "r2"), basis_id = "right:v1")
  left <- task_relation(c("e1", "e2"), left_space, domain)
  right <- task_relation(c("r1", "r2"), right_space, domain)
  over <- pairing(c("e1", "e2"), c("r1", "r2"), c(1, 2), directed = TRUE)
  base <- crossform:::.compile_effect_task(
    left, over, right, pair_query(diag(2), left_space, right_space)
  )
  changed_query <- crossform:::.compile_effect_task(
    left, over, right, pair_query(matrix(c(1, 1, 0, 1), 2),
      left_space, right_space)
  )
  changed_edges <- crossform:::.compile_effect_task(
    left,
    pairing(rev(over$left), rev(over$right), rev(over$weight), directed = TRUE),
    right,
    pair_query(diag(2), left_space, right_space)
  )
  left_other <- task_relation(c("e1", "e2"), left_space, other_domain)
  right_other <- task_relation(c("r1", "r2"), right_space, other_domain)
  changed_bridge <- crossform:::.compile_effect_task(
    left_other, over, right_other, pair_query(diag(2), left_space, right_space)
  )

  expect_length(unique(c(
    base$task_id,
    changed_query$task_id,
    changed_edges$task_id,
    changed_bridge$task_id
  )), 4L)
  mutated <- base
  mutated$reducer$order <- 2L
  expect_error(crossform:::.validate_compiled_effect_task(mutated),
    "reducer")
})

test_that("legacy compiler plan hashes use canonical ordered task semantics", {
  domain <- abstract_domain(3, id = "shared:v1")
  rel <- task_relation(c("r1", "r2", "r3"), c("a", "b"), domain)
  frame <- additive_frame(diag(3), domain = domain)
  over <- cross_partitions(rel)
  flipped <- pairing(over$right, over$left, over$weight, directed = FALSE)
  reordered <- pairing(
    rev(over$left), rev(over$right), rev(over$weight), directed = FALSE
  )
  plan <- function(pairing) crossform:::.compiler_plan_id(
    rel, frame, pairing, "full_geometry", NULL, "full"
  )

  expect_identical(plan(flipped), plan(over))
  expect_false(identical(plan(reordered), plan(over)))
})

test_that("relation-family provenance is portable and changes identity", {
  domain <- abstract_domain(2)
  first <- relation(
    list(r1 = matrix(1, 2, 2), r2 = matrix(2, 2, 2)),
    effects = c("a", "b"),
    domain = domain,
    provenance = list(estimator = "v1")
  )
  changed <- relation(
    list(r1 = matrix(1, 2, 2), r2 = matrix(2, 2, 2)),
    effects = c("a", "b"),
    domain = domain,
    provenance = list(estimator = "v2")
  )

  expect_false(identical(
    crossform:::.relation_family_identity(first),
    crossform:::.relation_family_identity(changed)
  ))
  expect_error(relation(
    list(r1 = matrix(1, 2, 2)), effects = c("a", "b"),
    domain = domain, provenance = list(callback = identity)
  ), "nonportable")
})

test_that("same-relation task sessions deduplicate source handles", {
  values <- matrix(seq_len(8), 2, 4)
  path <- tempfile(fileext = ".bin")
  on.exit(unlink(path), add = TRUE)
  connection <- file(path, open = "wb")
  writeBin(as.double(values), connection, size = 8, endian = .Platform$endian)
  close(connection)
  descriptor <- file_matrix_source(path, dim(values))
  domain <- abstract_domain(4)
  rel <- relation(
    list(r1 = descriptor, r2 = descriptor),
    effects = c("a", "b"),
    domain = domain
  )
  task <- crossform:::.compile_effect_task(rel, cross_partitions(rel))
  opens <- 0L
  opener <- function(...) {
    opens <<- opens + 1L
    crossform:::.open_source_descriptor(...)
  }
  session <- crossform:::.open_effect_task_source_session(
    task, open_descriptor = opener
  )
  on.exit(session$close(), add = TRUE)

  expect_length(task$distinct_handle_keys, 1L)
  expect_identical(opens, 1L)
  expect_equal(session$read("left", "r1", 1:2), values[, 1:2], tolerance = 0)
  expect_equal(session$read("right", "r2", 3:4), values[, 3:4], tolerance = 0)
  expect_identical(opens, 1L)
  expect_identical(session$summary()$distinct_owned_handles, 1L)
})
