configuration_form_fixture <- function(singleton = FALSE) {
  set.seed(20260814)
  left_partitions <- c("encode_1", "encode_2")
  right_partitions <- c("retrieve_1", "retrieve_2")
  left_effects <- c("encode_a", "encode_b")
  right_effects <- c("retrieve_a", "retrieve_b", "retrieve_c")
  features <- 5L
  left <- stats::setNames(lapply(left_partitions, function(partition) {
    matrix(rnorm(length(left_effects) * features),
      length(left_effects), features)
  }), left_partitions)
  right <- stats::setNames(lapply(right_partitions, function(partition) {
    matrix(rnorm(length(right_effects) * features),
      length(right_effects), features)
  }), right_partitions)
  weights <- if (singleton) diag(features) else matrix(c(
    1, 2, 0, 3, 1,
    0.5, 1, 2, 0.5, 4,
    3, 0, 1, 2, 1
  ), 3, features, byrow = TRUE)
  over <- pairing(
    left_partitions, right_partitions,
    weight = c(1, 3), directed = TRUE
  )
  list(
    left = left,
    right = right,
    left_partitions = left_partitions,
    right_partitions = right_partitions,
    left_effects = left_effects,
    right_effects = right_effects,
    frame = additive_frame(Matrix::Matrix(weights, sparse = TRUE)),
    edges = crossform:::.ordered_partition_edges(
      over, left_partitions, right_partitions, FALSE
    )
  )
}

run_configuration_form <- function(fixture, query = NULL,
                                   form_total = TRUE) {
  crossform:::.streamed_effect_form_contraction(
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
    feature_block = 2L,
    row_tile = 2L,
    coordinate_tile = 2L,
    retain_first_moments = TRUE,
    form_total = form_total
  )
}

run_configuration_component <- function(fixture, component, query = NULL) {
  crossform:::.streamed_effect_form_components(
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
    feature_block = 2L,
    row_tile = 2L,
    coordinate_tile = 2L
  )
}

configuration_oracle <- function(fixture) {
  weights <- as.matrix(fixture$frame$weights)
  total <- coherent <- centered <- vector("list", nrow(weights))
  for (measurement in seq_len(nrow(weights))) {
    mass <- sum(weights[measurement, ])
    total_form <- coherent_form <- centered_form <- matrix(
      0, length(fixture$left_effects), length(fixture$right_effects)
    )
    for (edge in seq_len(nrow(fixture$edges))) {
      left <- fixture$left[[fixture$edges$left[[edge]]]]
      right <- fixture$right[[fixture$edges$right[[edge]]]]
      left_first <- drop(left %*% weights[measurement, ])
      right_first <- drop(right %*% weights[measurement, ])
      edge_total <- matrix(0, nrow(left), nrow(right))
      edge_centered <- matrix(0, nrow(left), nrow(right))
      for (feature in seq_len(ncol(left))) {
        edge_total <- edge_total + weights[measurement, feature] *
          outer(left[, feature], right[, feature])
        edge_centered <- edge_centered + weights[measurement, feature] *
          outer(
            left[, feature] - left_first / mass,
            right[, feature] - right_first / mass
          )
      }
      edge_weight <- fixture$edges$weight[[edge]]
      total_form <- total_form + edge_weight * edge_total
      coherent_form <- coherent_form + edge_weight *
        outer(left_first, right_first) / mass
      centered_form <- centered_form + edge_weight * edge_centered
    }
    total[[measurement]] <- total_form
    coherent[[measurement]] <- coherent_form
    centered[[measurement]] <- centered_form
  }
  list(
    total = do.call(rbind, lapply(total, as.vector)),
    coherent = do.call(rbind, lapply(coherent, as.vector)),
    centered = do.call(rbind, lapply(centered, as.vector))
  )
}

test_that("rectangular configuration equals explicit weighted centering", {
  fixture <- configuration_form_fixture()
  streamed <- run_configuration_form(fixture)
  coherent <- crossform:::.effect_form_coherent_from_first_moments(
    streamed$first_moments$left,
    streamed$first_moments$right,
    fixture$edges,
    streamed$mass,
    codec = "rectangular",
    row_tile = 2L
  )
  oracle <- configuration_oracle(fixture)

  expect_equal(streamed$value, oracle$total, tolerance = 1e-12)
  expect_equal(coherent$value, oracle$coherent, tolerance = 1e-12)
  expect_equal(
    streamed$value - coherent$value,
    oracle$centered,
    tolerance = 1e-12
  )
  expect_gt(streamed$diagnostics$durable_left_first_moment_bytes, 0)
  expect_gt(streamed$diagnostics$durable_right_first_moment_bytes, 0)
  expect_gt(streamed$diagnostics$max_left_first_product_bytes, 0)
  expect_gt(streamed$diagnostics$max_right_first_product_bytes, 0)
})

test_that("direct total, coherent, and configuration queries equal late queries", {
  fixture <- configuration_form_fixture()
  complete <- run_configuration_form(fixture)
  complete_coherent <- crossform:::.effect_form_coherent_from_first_moments(
    complete$first_moments$left,
    complete$first_moments$right,
    fixture$edges,
    complete$mass,
    codec = "rectangular"
  )$value
  operators <- list(
    matrix(c(1, -2, 3, 0, 4, -1), 2, 3),
    matrix(c(0, 2, -1, 3, 1, 4), 2, 3)
  )
  query <- do.call(cbind, lapply(operators, as.vector))
  direct <- run_configuration_form(fixture, query = query)
  direct_coherent <- crossform:::.effect_form_coherent_from_first_moments(
    direct$first_moments$left,
    direct$first_moments$right,
    fixture$edges,
    direct$mass,
    codec = "rectangular",
    query = query
  )$value

  expect_equal(direct$value, complete$value %*% query, tolerance = 1e-12)
  expect_equal(
    direct_coherent,
    complete_coherent %*% query,
    tolerance = 1e-12
  )
  expect_equal(
    direct$value - direct_coherent,
    (complete$value - complete_coherent) %*% query,
    tolerance = 1e-12
  )
  lowered_configuration <- run_configuration_component(
    fixture, "configuration", query
  )
  expect_equal(
    lowered_configuration$value,
    (complete$value - complete_coherent) %*% query,
    tolerance = 1e-12
  )
  expect_gt(
    lowered_configuration$diagnostics$stream[[
      "durable_left_first_moment_bytes"
    ]],
    0
  )
  expect_gt(
    lowered_configuration$diagnostics$stream[[
      "durable_right_first_moment_bytes"
    ]],
    0
  )
  expect_gte(
    lowered_configuration$memory$categories[["local_state"]],
    lowered_configuration$diagnostics$stream[[
      "durable_left_first_moment_bytes"
    ]] + lowered_configuration$diagnostics$stream[[
      "durable_right_first_moment_bytes"
    ]]
  )

  left_space <- effect_space(fixture$left_effects, basis_id = "encode:v1")
  right_space <- effect_space(fixture$right_effects, basis_id = "retrieve:v1")
  result <- crossform:::effect_form(
    complete$value, left_space, right_space, receipt_fixture(),
    codec = "rectangular", coherent = complete_coherent
  )
  expect_equal(
    unname(query_geometry(
      result,
      pair_query(operators[[1L]], left_space, right_space),
      component = "configuration"
    )$values),
    (complete$value - complete_coherent) %*% query[, 1L, drop = FALSE],
    tolerance = 1e-12
  )
})

test_that("singleton spatial frames have zero configuration", {
  fixture <- configuration_form_fixture(singleton = TRUE)
  streamed <- run_configuration_form(fixture)
  coherent <- crossform:::.effect_form_coherent_from_first_moments(
    streamed$first_moments$left,
    streamed$first_moments$right,
    fixture$edges,
    streamed$mass,
    codec = "rectangular"
  )$value

  expect_equal(streamed$value - coherent, matrix(0, nrow(streamed$value),
    ncol(streamed$value)), tolerance = 1e-14)
})

test_that("first-moment-only streaming omits total atoms and output", {
  fixture <- configuration_form_fixture()
  got <- run_configuration_form(fixture, form_total = FALSE)

  expect_null(got$value)
  expect_equal(got$diagnostics$atom_count, 0)
  expect_equal(got$diagnostics$max_atom_bytes, 0)
  expect_gt(got$diagnostics$durable_left_first_moment_bytes, 0)
  expect_gt(got$diagnostics$durable_right_first_moment_bytes, 0)
})

test_that("rectangular first-moment memory plans cover both side states", {
  fixture <- configuration_form_fixture()
  got <- run_configuration_form(fixture)
  plan <- crossform:::.effect_form_kernel_memory_plan(
    fixture$frame,
    fixture$left_effects,
    fixture$right_effects,
    fixture$left_partitions,
    fixture$right_partitions,
    codec = "rectangular",
    feature_block = 2L,
    row_tile = 2L,
    coordinate_tile = 2L,
    storage = "memory",
    retain_first_moments = TRUE
  )
  durable <- got$diagnostics$durable_output_bytes +
    got$diagnostics$durable_left_first_moment_bytes +
    got$diagnostics$durable_right_first_moment_bytes

  expect_gte(plan$categories[["local_state"]],
    got$diagnostics$durable_left_first_moment_bytes +
      got$diagnostics$durable_right_first_moment_bytes)
  expect_gte(plan$planned_workspace_bytes,
    durable + got$diagnostics$max_live_temporary_bytes)
})

test_that("packed self coherent specialization remains exact", {
  set.seed(511)
  partitions <- c("run1", "run2", "run3")
  effects <- c("a", "b", "c")
  local <- array(rnorm(4 * 3 * 3), c(4, 3, 3),
    dimnames = list(NULL, effects, partitions))
  over <- pairing(c("run1", "run1", "run2"),
    c("run2", "run3", "run3"), weight = c(1, 2, 3))
  mass <- c(2, 3, 4, 5)
  got <- crossform:::.coherent_geometry_from_local(local, over, mass)
  oracle <- matrix(0, 4, 6)
  for (measurement in 1:4) {
    form <- matrix(0, 3, 3)
    for (edge in seq_len(nrow(over))) {
      left <- local[measurement, , over$left[[edge]]]
      right <- local[measurement, , over$right[[edge]]]
      form <- form + over$weight[[edge]] *
        0.5 * (outer(left, right) + outer(right, left)) / mass[[measurement]]
    }
    oracle[measurement, ] <- crossform:::.svec_symmetric(form)
  }

  expect_equal(got$value, oracle, tolerance = 1e-14)
})
