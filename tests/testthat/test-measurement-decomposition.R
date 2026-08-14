decomposition_measurement_fixture <- function(seed = 2026081224) {
  set.seed(seed)
  q <- 4L
  p <- 5L
  domain <- abstract_domain(p, id = "decomposition:neural:v1")
  effects <- effect_space(paste0("condition", seq_len(q)),
    basis_id = "decomposition:effects:v1")
  values <- list(
    run1 = matrix(rnorm(q * p), q, p),
    run2 = matrix(rnorm(q * p), q, p)
  )
  rel <- relation(values, effects = effects, domain = domain)
  weights <- matrix(c(
    1, 2, 0, 3, 1,
    0, 1, 4, 2, 1
  ), 2, p, byrow = TRUE)
  additive <- additive_frame(weights, domain = domain)
  frame <- effectagram:::.measurement_frame_from_additive_decomposition(additive)
  pairs <- expand.grid(
    left = frame$node_ids,
    right = frame$node_ids,
    stringsAsFactors = FALSE
  )
  spatial_edges <- effectagram:::.measurement_edges(
    pairs$left, pairs$right, frame
  )
  over <- pairing(
    rel$partitions, rel$partitions,
    directed = TRUE,
    self_pairs = "allow_biased",
    independence = "not_independent"
  )
  partition_edges <- effectagram:::.ordered_partition_edges(
    over, rel$partitions, rel$partitions, TRUE
  )
  h <- crossprod(matrix(rnorm(q * q), q, q))
  query <- pair_query(h, effects, effects)
  task <- effectagram:::.new_evidence_task(
    rel, rel, TRUE, partition_edges,
    effectagram:::.closed_experimental_boundary(query),
    effectagram:::.open_neural_boundary(frame, frame, spatial_edges),
    effectagram:::.evidence_stage_plan(),
    effectagram:::.evidence_materialization(
      "measurement_form", "complete_form"
    )
  )
  run <- effectagram:::.run_measurement_contraction(
    task, compute_policy(block_features = 2L), route = "pull_h"
  )
  form <- effectagram:::.measurement_form_from_contraction(task, run)
  list(
    domain = domain,
    effects = effects,
    values = values,
    weights = weights,
    additive = additive,
    frame = frame,
    task = task,
    run = run,
    form = form,
    h = h
  )
}

test_that("coherent/configuration factors recompose the additive metric exactly", {
  fixture <- decomposition_measurement_fixture()
  for (node in seq_len(nrow(fixture$weights))) {
    leg <- fixture$frame$legs[[node]]
    decomposition <- leg$decomposition
    coherent_range <- decomposition$ranges$coherent
    configuration_range <- decomposition$ranges$configuration
    coherent <- leg$operator[coherent_range, , drop = FALSE]
    configuration <- leg$operator[configuration_range, , drop = FALSE]
    expected_total <- diag(fixture$weights[node, ])
    expected_coherent <- tcrossprod(fixture$weights[node, ]) /
      sum(fixture$weights[node, ])

    expect_equal(crossprod(leg$operator), expected_total, tolerance = 2e-13)
    expect_equal(crossprod(coherent), expected_coherent, tolerance = 2e-13)
    expect_equal(crossprod(configuration),
      expected_total - expected_coherent, tolerance = 3e-13)
    expect_identical(decomposition$orientation,
      c(coherent = "oriented", configuration = "subspace_only"))
  }
})

test_that("coherent and configuration node geometries retain certified values", {
  fixture <- decomposition_measurement_fixture(seed = 2026081225)
  relation <- fixture$values[[1L]]
  for (node in seq_len(nrow(fixture$weights))) {
    weight <- fixture$weights[node, ]
    mass <- sum(weight)
    leg <- fixture$frame$legs[[node]]
    coherent_operator <- leg$operator[
      leg$decomposition$ranges$coherent, , drop = FALSE
    ]
    configuration_operator <- leg$operator[
      leg$decomposition$ranges$configuration, , drop = FALSE
    ]
    total <- relation %*% diag(weight) %*% t(relation)
    first <- drop(relation %*% weight)
    coherent <- tcrossprod(first) / mass
    configuration <- total - coherent

    expect_equal(
      relation %*% t(coherent_operator) %*%
        (coherent_operator %*% t(relation)),
      coherent, tolerance = 3e-13
    )
    expect_equal(
      relation %*% t(configuration_operator) %*%
        (configuration_operator %*% t(relation)),
      configuration, tolerance = 4e-13
    )
  }
})

test_that("every crossed component block exactly reassembles its edge", {
  fixture <- decomposition_measurement_fixture(seed = 2026081226)
  for (edge in seq_len(nrow(fixture$form$block_index))) {
    lifted <- effectagram:::.lift_measurement_decomposition(
      fixture$form, edge
    )
    expect_s3_class(lifted, "effect_measurement_decomposition_view")
    expect_identical(nrow(lifted$index), 4L)
    expect_identical(names(lifted$blocks), c(
      "coherent::coherent", "configuration::coherent",
      "coherent::configuration", "configuration::configuration"
    ))
    expect_equal(
      effectagram:::.recompose_measurement_decomposition(lifted),
      effectagram:::.measurement_block(fixture$form, edge),
      tolerance = 0
    )
    expect_silent(effectagram:::.validate_measurement_decomposition_view(
      lifted
    ))
  }
})

test_that("raw entries require oriented bases while coordinate access remains exact", {
  fixture <- decomposition_measurement_fixture(seed = 2026081227)
  lifted <- effectagram:::.lift_measurement_decomposition(fixture$form, 2L)

  expect_silent(effectagram:::.measurement_component_block(
    lifted, "coherent", "coherent", interpretation = "raw"
  ))
  expect_error(effectagram:::.measurement_component_block(
    lifted, "configuration", "configuration", interpretation = "raw"
  ), "oriented bases|subspace")
  coordinate <- effectagram:::.measurement_component_block(
    lifted, "configuration", "configuration", interpretation = "coordinate"
  )
  position <- which(lifted$index$left_component == "configuration" &
    lifted$index$right_component == "configuration")
  expect_identical(coordinate, lifted$blocks[[position]])
})

test_that("reversal swaps crossed components and transposes every block", {
  fixture <- decomposition_measurement_fixture(seed = 2026081228)
  reversed <- effectagram:::.reverse_measurement_form(fixture$form)

  for (edge in seq_len(nrow(fixture$form$block_index))) {
    forward_lift <- effectagram:::.lift_measurement_decomposition(
      fixture$form, edge
    )
    reverse_lift <- effectagram:::.lift_measurement_decomposition(
      reversed, edge
    )
    for (row in seq_len(nrow(forward_lift$index))) {
      left <- forward_lift$index$left_component[[row]]
      right <- forward_lift$index$right_component[[row]]
      forward_block <- effectagram:::.measurement_component_block(
        forward_lift, left, right, "coordinate"
      )
      reverse_block <- effectagram:::.measurement_component_block(
        reverse_lift, right, left, "coordinate"
      )
      expect_equal(reverse_block, t(forward_block), tolerance = 0)
    }
  }
})

test_that("invariant summaries survive independent configuration rotations", {
  set.seed(2026081229)
  cross <- matrix(rnorm(12), 3, 4)
  left_seed <- matrix(rnorm(3 * 7), 3, 7)
  right_seed <- matrix(rnorm(4 * 8), 4, 8)
  left_self <- tcrossprod(left_seed) + diag(0.5, 3)
  right_self <- tcrossprod(right_seed) + diag(0.5, 4)
  left_rotation <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  right_rotation <- qr.Q(qr(matrix(rnorm(16), 4, 4)))
  rotated_cross <- left_rotation %*% cross %*% t(right_rotation)
  rotated_left <- left_rotation %*% left_self %*% t(left_rotation)
  rotated_right <- right_rotation %*% right_self %*% t(right_rotation)

  original <- effectagram:::.measurement_invariant_summary(
    cross, left_self, right_self
  )
  rotated <- effectagram:::.measurement_invariant_summary(
    rotated_cross, rotated_left, rotated_right
  )

  expect_false(isTRUE(all.equal(cross, rotated_cross, tolerance = 1e-12)))
  expect_equal(rotated$frobenius_norm, original$frobenius_norm,
    tolerance = 1e-12)
  expect_equal(rotated$singular_values, original$singular_values,
    tolerance = 1e-12)
  expect_equal(rotated$canonical_values, original$canonical_values,
    tolerance = 2e-11)
  expect_equal(rotated$subspace_angles, original$subspace_angles,
    tolerance = 2e-11)
  expect_equal(rotated$geometry_alignment, original$geometry_alignment,
    tolerance = 2e-12)
})

test_that("singleton nodes retain a zero-dimensional configuration component", {
  domain <- abstract_domain(3)
  frame <- additive_frame(matrix(c(0, 2, 0), 1), domain = domain)
  decomposition <- effectagram:::.measurement_frame_from_additive_decomposition(
    frame
  )
  leg <- decomposition$legs[[1L]]

  expect_identical(leg$decomposition$ranges$coherent, 1L)
  expect_identical(leg$decomposition$ranges$configuration, integer())
  expect_identical(dim(leg$operator), c(1L, 3L))
  expect_equal(crossprod(leg$operator), diag(c(0, 2, 0)), tolerance = 1e-15)
})

test_that("future named decompositions use the same crossed lifting law", {
  domain <- abstract_domain(3)
  low <- effectagram:::.measurement_leg(
    matrix(c(1, 0, 0), 1), domain,
    effectagram:::.measurement_axis("low", "frequency:low")
  )
  high <- effectagram:::.measurement_leg(
    matrix(c(0, 1, 1), 1), domain,
    effectagram:::.measurement_axis("high", "frequency:high")
  )
  combined <- effectagram:::.direct_sum_measurement_leg(
    list(low = low, high = high), "frequency:v1"
  )

  expect_identical(combined$decomposition$components, c("low", "high"))
  expect_identical(combined$decomposition$orientation,
    c(low = "oriented", high = "oriented"))
  expect_equal(crossprod(combined$operator),
    crossprod(low$operator) + crossprod(high$operator), tolerance = 0)
})
