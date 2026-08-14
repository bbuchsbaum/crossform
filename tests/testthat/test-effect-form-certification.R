certification_form_fixture <- function(sparse = FALSE) {
  set.seed(2026081301)
  left_partitions <- c("encode_1", "encode_2")
  right_partitions <- c("retrieve_1", "retrieve_2")
  left_effects <- c("e1", "e2", "e3")
  right_effects <- c("r1", "r2", "r3", "r4")
  features <- 7L
  left <- stats::setNames(lapply(left_partitions, function(partition) {
    matrix(rnorm(length(left_effects) * features), length(left_effects))
  }), left_partitions)
  right <- stats::setNames(lapply(right_partitions, function(partition) {
    matrix(rnorm(length(right_effects) * features), length(right_effects))
  }), right_partitions)
  weights <- matrix(runif(3L * features, 0.2, 2), 3L, features)
  if (sparse) weights <- Matrix::Matrix(weights, sparse = TRUE)
  pairing <- pairing(
    left_partitions, right_partitions, c(0.25, 0.75), directed = TRUE
  )
  list(
    left = left,
    right = right,
    left_partitions = left_partitions,
    right_partitions = right_partitions,
    left_effects = left_effects,
    right_effects = right_effects,
    frame = additive_frame(weights),
    edges = effectagram:::.ordered_partition_edges(
      pairing, left_partitions, right_partitions, FALSE
    )
  )
}

certification_run_form <- function(fixture, feature_block, query = NULL,
                                   component = "complete") {
  effectagram:::.streamed_effect_form_components(
    fixture$frame,
    read_left = function(partition, features) {
      fixture$left[[partition]][, features, drop = FALSE]
    },
    read_right = function(partition, features) {
      fixture$right[[partition]][, features, drop = FALSE]
    },
    left_partitions = fixture$left_partitions,
    right_partitions = fixture$right_partitions,
    left_effects = fixture$left_effects,
    right_effects = fixture$right_effects,
    ordered_edges = fixture$edges,
    codec = "rectangular",
    query = query,
    component = component,
    feature_block = feature_block,
    row_tile = 2L,
    coordinate_tile = 3L
  )
}

certification_form_oracle <- function(fixture) {
  weights <- as.matrix(fixture$frame$weights)
  total <- coherent <- vector("list", nrow(weights))
  for (measurement in seq_len(nrow(weights))) {
    total_form <- coherent_form <- matrix(
      0, length(fixture$left_effects), length(fixture$right_effects)
    )
    mass <- sum(weights[measurement, ])
    for (edge in seq_len(nrow(fixture$edges))) {
      left <- fixture$left[[fixture$edges$left[[edge]]]]
      right <- fixture$right[[fixture$edges$right[[edge]]]]
      edge_weight <- fixture$edges$weight[[edge]]
      total_form <- total_form + edge_weight *
        oracle_effect_form(left, right, weights[measurement, ])
      coherent_form <- coherent_form + edge_weight * outer(
        drop(left %*% weights[measurement, ]),
        drop(right %*% weights[measurement, ])
      ) / mass
    }
    total[[measurement]] <- as.vector(total_form)
    coherent[[measurement]] <- as.vector(coherent_form)
  }
  list(total = do.call(rbind, total), coherent = do.call(rbind, coherent))
}

test_that("rectangular complete and direct paths pass one independent court", {
  dense <- certification_form_fixture(FALSE)
  sparse <- certification_form_fixture(TRUE)
  oracle <- certification_form_oracle(dense)
  complete_dense <- certification_run_form(dense, 1L)
  complete_sparse <- certification_run_form(sparse, 4L)

  expect_equal(complete_dense$total, oracle$total, tolerance = 1e-12)
  expect_equal(complete_dense$coherent, oracle$coherent, tolerance = 1e-12)
  expect_equal(complete_sparse$total, oracle$total, tolerance = 1e-12)
  expect_equal(complete_sparse$coherent, oracle$coherent, tolerance = 1e-12)

  query <- cbind(
    c(1, -1, 0, 2, 0, -2, 1, 0, 3, -1, 0, 2),
    c(0, 2, -1, 0, 1, 0, -2, 3, 0, 1, -1, 0)
  )
  for (component in c("total", "coherent", "configuration")) {
    direct <- certification_run_form(dense, 2L, query, component)
    complete <- switch(component,
      total = oracle$total,
      coherent = oracle$coherent,
      configuration = oracle$total - oracle$coherent
    )
    expect_equal(direct$value, complete %*% query, tolerance = 1e-12)
  }
})

test_that("execution changes preserve the scientific identity promised", {
  domain <- abstract_domain(4, id = "cert:self-neural")
  effects <- effect_space(c("a", "b"), basis_id = "cert:self-effects")
  rel <- relation(list(
    run1 = matrix(1:8, 2, 4),
    run2 = matrix(2:9, 2, 4)
  ), effects = effects, domain = domain)
  frame <- additive_frame(diag(4), domain = domain)
  over <- cross_partitions(rel)
  path <- tempfile("effect-form-certification-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)

  memory <- geometry(rel, frame, over,
    compute = compute_policy(block_features = 1L))
  blocked <- geometry(rel, frame, over, storage = "block",
    storage_path = path, compute = compute_policy(block_features = 4L))

  expect_identical(memory$receipt$scientific_plan_id,
    blocked$receipt$scientific_plan_id)
  expect_false(identical(memory$receipt, blocked$receipt))
  expect_equal(geometry_component(memory), geometry_component(blocked),
    tolerance = 0)
})

test_that("the installed namespace exposes the universal public vocabulary", {
  expected <- c(
    "pair_query", "match_coupling", "control_coupling", "pair_contrast",
    "pair_lm_query", "match_control", "measurement_space",
    "measurement_bridge", "reverse_bridge", "inner_product",
    "reduce_partitions", "aggregate_first"
  )
  expect_true(all(expected %in% getNamespaceExports("effectagram")))
})
