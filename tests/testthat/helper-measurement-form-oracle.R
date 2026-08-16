# Independent dense oracle and fixtures for requested measurement blocks.

measurement_block_oracle <- function(left_relations, right_relations, h,
                                     partition_edges, measurement_edges,
                                     left_operators, right_operators) {
  blocks <- stats::setNames(lapply(seq_len(nrow(measurement_edges)), function(i) {
    matrix(0,
      nrow(left_operators[[measurement_edges$left[[i]]]]),
      nrow(right_operators[[measurement_edges$right[[i]]]])
    )
  }), paste0("edge_", sprintf("%06d", seq_len(nrow(measurement_edges)))))
  for (partition_edge in seq_len(nrow(partition_edges))) {
    left <- left_relations[[partition_edges$left[[partition_edge]]]]
    right <- right_relations[[partition_edges$right[[partition_edge]]]]
    for (edge in seq_len(nrow(measurement_edges))) {
      left_operator <- left_operators[[measurement_edges$left[[edge]]]]
      right_operator <- right_operators[[measurement_edges$right[[edge]]]]
      blocks[[edge]] <- blocks[[edge]] +
        partition_edges$weight[[partition_edge]] *
        left_operator %*% t(left) %*% h %*% right %*%
          t(right_operator)
    }
  }
  blocks
}

measurement_expect_blocks_equal <- function(got, expected, tolerance = 1e-11) {
  expect_identical(names(got), names(expected))
  expect_identical(vapply(got, dim, integer(2)),
    vapply(expected, dim, integer(2)))
  for (edge in seq_along(expected)) {
    expect_equal(got[[edge]], expected[[edge]], tolerance = tolerance,
      info = names(expected)[[edge]])
  }
}

measurement_kernel_fixture <- function(scalar = FALSE, low_rank_h = FALSE,
                                       seed = 2026081205,
                                       left_sources = NULL,
                                       right_sources = NULL) {
  set.seed(seed)
  q_left <- 4L
  q_right <- 3L
  p_left <- 5L
  p_right <- 6L
  left_partitions <- c("encode_1", "encode_2")
  right_partitions <- c("retrieve_1", "retrieve_2")
  left_values <- stats::setNames(lapply(left_partitions, function(...) {
    matrix(rnorm(q_left * p_left), q_left, p_left)
  }), left_partitions)
  right_values <- stats::setNames(lapply(right_partitions, function(...) {
    matrix(rnorm(q_right * p_right), q_right, p_right)
  }), right_partitions)
  if (is.null(left_sources)) left_sources <- left_values
  if (is.null(right_sources)) right_sources <- right_values
  left_domain <- abstract_domain(p_left, id = "kernel:left-neural:v1")
  right_domain <- abstract_domain(p_right, id = "kernel:right-neural:v1")
  left_space <- effect_space(paste0("l", seq_len(q_left)),
    basis_id = "kernel:left-effects:v1")
  right_space <- effect_space(paste0("r", seq_len(q_right)),
    basis_id = "kernel:right-effects:v1")
  left_relation <- relation(
    left_sources, effects = left_space, domain = left_domain
  )
  right_relation <- relation(
    right_sources, effects = right_space, domain = right_domain
  )
  left_operators <- if (scalar) {
    list(x1 = matrix(rnorm(p_left), 1), x2 = matrix(rnorm(p_left), 1),
      x3 = matrix(rnorm(p_left), 1))
  } else {
    list(x1 = matrix(rnorm(2 * p_left), 2, p_left),
      x2 = matrix(rnorm(3 * p_left), 3, p_left),
      x3 = matrix(rnorm(2 * p_left), 2, p_left))
  }
  right_operators <- if (scalar) {
    list(y1 = matrix(rnorm(p_right), 1), y2 = matrix(rnorm(p_right), 1),
      y3 = matrix(rnorm(p_right), 1))
  } else {
    list(y1 = matrix(rnorm(2 * p_right), 2, p_right),
      y2 = matrix(rnorm(2 * p_right), 2, p_right),
      y3 = matrix(rnorm(3 * p_right), 3, p_right))
  }
  make_legs <- function(operators, domain, prefix) {
    legs <- Map(function(operator, node) {
      crossform:::.measurement_leg(
        operator, domain,
        crossform:::.measurement_axis(
          paste0("mode", seq_len(nrow(operator))),
          paste0(prefix, ":", node),
          basis_id = paste0(prefix, ":fixed-basis")
        )
      )
    }, operators, names(operators))
    names(legs) <- names(operators)
    legs
  }
  left_frame <- crossform:::.measurement_frame(
    make_legs(left_operators, left_domain, "left-frame")
  )
  right_frame <- crossform:::.measurement_frame(
    make_legs(right_operators, right_domain, "right-frame")
  )
  measurement_edges <- crossform:::.measurement_edges(
    c("x1", "x2", "x3", "x1"),
    c("y2", "y1", "y3", "y3"),
    left_frame, right_frame,
    weight = c(1, 0.5, -0.25, 2)
  )
  over <- pairing(
    left_partitions, right_partitions, c(0.35, 0.65), directed = TRUE
  )
  partition_edges <- crossform:::.ordered_partition_edges(
    over, left_partitions, right_partitions, FALSE
  )
  h <- if (low_rank_h) {
    matrix(rnorm(q_left * 2L), q_left, 2L) %*%
      t(matrix(rnorm(q_right * 2L), q_right, 2L))
  } else {
    matrix(rnorm(q_left * q_right), q_left, q_right)
  }
  query <- pair_query(h, left_space, right_space)
  task <- crossform:::.new_evidence_task(
    left_relation, right_relation, FALSE, partition_edges,
    crossform:::.closed_experimental_boundary(query),
    crossform:::.open_neural_boundary(
      left_frame, right_frame, measurement_edges
    ),
    crossform:::.evidence_stage_plan(),
    crossform:::.evidence_materialization(
      "measurement_form", "complete_form"
    )
  )
  list(
    task = task,
    left_values = left_values,
    right_values = right_values,
    h = h,
    partition_edges = partition_edges,
    measurement_edges = measurement_edges$edges,
    left_operators = left_operators,
    right_operators = right_operators
  )
}
