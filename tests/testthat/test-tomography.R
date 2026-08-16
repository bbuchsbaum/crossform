tomography_frame_from_operator <- function(operator, domain, node_widths,
                                           prefix) {
  ends <- cumsum(node_widths)
  starts <- c(1L, head(ends, -1L) + 1L)
  ranges <- Map(seq.int, starts, ends)
  names(ranges) <- paste0("node", seq_along(ranges))
  legs <- Map(function(rows, node) {
    crossform:::.measurement_leg(
      operator[rows, , drop = FALSE], domain,
      crossform:::.measurement_axis(
        paste0(node, "_", seq_along(rows)),
        paste0(prefix, ":", node, ":v1"),
        basis_id = paste0(prefix, ":", node, ":basis:v1")
      )
    )
  }, ranges, names(ranges))
  names(legs) <- names(ranges)
  crossform:::.measurement_frame(legs)
}

tomography_self_fixture <- function(
    kind = c("parseval", "general", "deficient", "ill_conditioned"),
    complete = TRUE, seed = 2026081226) {
  kind <- match.arg(kind)
  set.seed(seed)
  p <- if (kind == "ill_conditioned") 2L else 3L
  q <- 7L
  domain <- abstract_domain(p, id = paste0("tomography:", kind, ":neural:v1"))
  effects <- effect_space(paste0("sample", seq_len(q)),
    basis_id = paste0("tomography:", kind, ":samples:v1"))
  b <- matrix(rnorm(q * p), q, p)
  values <- list(run1 = b)
  relation_value <- relation(values, effects = effects, domain = domain)
  h <- (diag(q) - matrix(1 / q, q, q)) / (q - 1)
  query <- crossform:::.variation_pair_query(
    h, effects, "trial", "joint_covariance",
    provenance = list(estimator = "centered-single-partition")
  )
  operator <- switch(kind,
    parseval = diag(3),
    general = matrix(c(
      1, 0, 0,
      0, 2, 0,
      1, 1, 0,
      0, 1, 1
    ), 4L, 3L, byrow = TRUE),
    deficient = matrix(c(
      1, 0, 0,
      0, 1, 0
    ), 2L, 3L, byrow = TRUE),
    ill_conditioned = diag(c(1, 1e-12))
  )
  widths <- switch(kind,
    parseval = c(2L, 1L),
    general = c(2L, 2L),
    deficient = c(1L, 1L),
    ill_conditioned = 2L
  )
  frame <- tomography_frame_from_operator(
    operator, domain, widths, paste0("tomography:", kind)
  )
  pairs <- if (complete) {
    expand.grid(
      left = frame$node_ids,
      right = frame$node_ids,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      left = frame$node_ids,
      right = frame$node_ids,
      stringsAsFactors = FALSE
    )
  }
  spatial_edges <- crossform:::.measurement_edges(
    pairs$left, pairs$right, frame
  )
  over <- pairing(
    "run1", "run1", directed = TRUE,
    self_pairs = "allow_biased", independence = "not_independent"
  )
  partition_edges <- crossform:::.ordered_partition_edges(
    over, "run1", "run1", TRUE
  )
  task <- crossform:::.new_evidence_task(
    relation_value, relation_value, TRUE, partition_edges,
    crossform:::.closed_experimental_boundary(
      query, "variation", "trial"
    ),
    crossform:::.open_neural_boundary(frame, frame, spatial_edges),
    crossform:::.evidence_stage_plan(),
    crossform:::.evidence_materialization(
      "measurement_form", "complete_form"
    )
  )
  contraction <- crossform:::.run_measurement_contraction(
    task, route = "pull_h"
  )
  form <- crossform:::.measurement_form_from_contraction(
    task, contraction,
    query_construction = "joint_covariance",
    edge_scope = if (complete) "frame_complete" else "requested"
  )
  list(
    form = form,
    frame = frame,
    stacked = operator,
    reference = crossprod(b, h %*% b)
  )
}

tomography_asymmetric_fixture <- function(seed = 2026081227) {
  set.seed(seed)
  q_left <- 6L
  q_right <- 5L
  p_left <- 3L
  p_right <- 2L
  left_domain <- abstract_domain(p_left, id = "tomography:asym:left:v1")
  right_domain <- abstract_domain(p_right, id = "tomography:asym:right:v1")
  left_effects <- effect_space(paste0("l", seq_len(q_left)),
    basis_id = "tomography:asym:left-effects:v1")
  right_effects <- effect_space(paste0("r", seq_len(q_right)),
    basis_id = "tomography:asym:right-effects:v1")
  left_b <- matrix(rnorm(q_left * p_left), q_left, p_left)
  right_b <- matrix(rnorm(q_right * p_right), q_right, p_right)
  left_relation <- relation(list(left_run = left_b),
    effects = left_effects, domain = left_domain)
  right_relation <- relation(list(right_run = right_b),
    effects = right_effects, domain = right_domain)
  h <- matrix(rnorm(q_left * q_right), q_left, q_right)
  query <- pair_query(h, left_effects, right_effects)
  left_operator <- matrix(c(
    1, 0, 0,
    0, 2, 0,
    1, 1, 1,
    0, 1, -1
  ), 4L, p_left, byrow = TRUE)
  right_operator <- matrix(c(
    1, 0,
    0, 2,
    1, -1
  ), 3L, p_right, byrow = TRUE)
  left_frame <- tomography_frame_from_operator(
    left_operator, left_domain, c(2L, 2L), "tomography:asym:left"
  )
  right_frame <- tomography_frame_from_operator(
    right_operator, right_domain, c(1L, 2L), "tomography:asym:right"
  )
  pairs <- expand.grid(
    left = left_frame$node_ids,
    right = right_frame$node_ids,
    stringsAsFactors = FALSE
  )
  spatial_edges <- crossform:::.measurement_edges(
    pairs$left, pairs$right, left_frame, right_frame
  )
  over <- pairing("left_run", "right_run", directed = TRUE)
  partition_edges <- crossform:::.ordered_partition_edges(
    over, "left_run", "right_run", FALSE
  )
  task <- crossform:::.new_evidence_task(
    left_relation, right_relation, FALSE, partition_edges,
    crossform:::.closed_experimental_boundary(query),
    crossform:::.open_neural_boundary(
      left_frame, right_frame, spatial_edges
    ),
    crossform:::.evidence_stage_plan(),
    crossform:::.evidence_materialization(
      "measurement_form", "complete_form"
    )
  )
  contraction <- crossform:::.run_measurement_contraction(
    task, route = "pull_h"
  )
  form <- crossform:::.measurement_form_from_contraction(
    task, contraction, edge_scope = "frame_complete"
  )
  list(
    form = form,
    left_frame = left_frame,
    right_frame = right_frame,
    reference = crossprod(left_b, h %*% right_b)
  )
}

test_that("Parseval frames assemble and reconstruct the global operator", {
  fixture <- tomography_self_fixture("parseval")
  assembled <- crossform:::.assemble_measurement_blocks(
    fixture$form, fixture$frame, fixture$frame
  )
  result <- crossform:::.reconstruct_neural_evidence(
    fixture$form, fixture$frame,
    reference_operator = fixture$reference
  )

  expect_equal(assembled$value, fixture$reference, tolerance = 1e-13)
  expect_identical(result$method, "parseval")
  expect_identical(result$status,
    "numerically_certified_exact_reconstruction")
  expect_true(result$lossless)
  expect_true(result$certified)
  expect_equal(result$operator, fixture$reference, tolerance = 1e-13)
  expect_lte(result$diagnostics$relative_reconstruction_residual, 1e-13)
  expect_silent(crossform:::.validate_measured_block_form(assembled))
  expect_silent(crossform:::.validate_tomography_result(result))

  algebraic <- crossform:::.reconstruct_neural_evidence(
    fixture$form, fixture$frame
  )
  expect_identical(algebraic$status, "exact_algebraic_reconstruction")
  expect_false(algebraic$certified)
})

test_that("general dual frames reconstruct and reverse independently", {
  fixture <- tomography_asymmetric_fixture()
  result <- crossform:::.reconstruct_neural_evidence(
    fixture$form, fixture$left_frame, fixture$right_frame,
    reference_operator = fixture$reference
  )
  expect_identical(result$method, "canonical_dual")
  expect_true(result$lossless)
  expect_equal(result$operator, fixture$reference, tolerance = 1e-10)

  reversed_form <- crossform:::.reverse_measurement_form(fixture$form)
  reversed <- crossform:::.reconstruct_neural_evidence(
    reversed_form, fixture$right_frame, fixture$left_frame,
    reference_operator = t(fixture$reference)
  )
  expect_identical(reversed$method, "canonical_dual")
  expect_equal(reversed$operator, t(result$operator), tolerance = 1e-10)
  expect_equal(reversed$left_projection, result$right_projection,
    tolerance = 1e-12)
  expect_equal(reversed$right_projection, result$left_projection,
    tolerance = 1e-12)
})

test_that("rank-deficient frames return only their projected operator", {
  fixture <- tomography_self_fixture("deficient")
  result <- crossform:::.reconstruct_neural_evidence(
    fixture$form, fixture$frame,
    reference_operator = fixture$reference
  )
  projection <- diag(c(1, 1, 0))
  expected <- projection %*% fixture$reference %*% projection

  expect_identical(result$method, "projected_pseudoinverse")
  expect_identical(result$status,
    "numerically_certified_projected_reconstruction")
  expect_false(result$lossless)
  expect_equal(result$left_projection, projection, tolerance = 1e-13)
  expect_equal(result$right_projection, projection, tolerance = 1e-13)
  expect_equal(result$operator, expected, tolerance = 1e-13)
  expect_identical(result$diagnostics$regularization$kind, "truncated_svd")
  expect_error(crossform:::.reconstruct_neural_evidence(
    fixture$form, fixture$frame, allow_projection = FALSE
  ), "projection is disabled", class = "effect_tomography_rejection")
})

test_that("diagonal-only blocks and incompatible bases cannot claim tomography", {
  incomplete <- tomography_self_fixture("deficient", complete = FALSE)
  expect_false(incomplete$form$capabilities$complete_edge_set)
  expect_error(crossform:::.assemble_measurement_blocks(
    incomplete$form, incomplete$frame, incomplete$frame
  ), "diagonal-only", class = "effect_capability_refusal")
  block_refusal <- catch_refusal(crossform:::.reconstruct_neural_evidence(
    incomplete$form, incomplete$frame
  ))
  expect_s3_class(block_refusal, "effect_capability_refusal")
  expect_identical(block_refusal$capability, "complete_edge_set")
  expect_identical(block_refusal$namespace, "tomography")
  expect_identical(block_refusal$reasons, "edge_set_is_not_frame_complete")

  q1 <- matrix(c(2, 0.75, 0.75, 3), 2L)
  q2 <- matrix(c(2, -0.75, -0.75, 3), 2L)
  expect_identical(diag(q1), diag(q2))
  expect_false(identical(q1, q2))

  complete <- tomography_self_fixture("parseval")
  altered <- tomography_frame_from_operator(
    complete$stacked,
    complete$frame$source_domain,
    c(2L, 1L),
    "tomography:altered-basis"
  )
  expect_error(crossform:::.reconstruct_neural_evidence(
    complete$form, altered
  ), "do not match.*bases")
})

test_that("ill-conditioning is rejected with frame diagnostics", {
  fixture <- tomography_self_fixture("ill_conditioned")
  error <- tryCatch(
    crossform:::.reconstruct_neural_evidence(
      fixture$form, fixture$frame,
      tolerance = 1e-15, max_condition = 1e8
    ),
    effect_tomography_rejection = identity
  )
  expect_s3_class(error, "effect_tomography_rejection")
  expect_match(conditionMessage(error), "ill-conditioned")
  expect_gt(error$diagnostics$left$condition_number, 1e8)
  expect_identical(error$diagnostics$regularization$ridge, 0)
})

test_that("resource preflight fails before any measurement block read", {
  fixture <- tomography_self_fixture("parseval")
  reads <- new.env(parent = emptyenv())
  reads$count <- 0L
  original_read <- fixture$form$store$read
  fixture$form$store$read <- function(edge) {
    reads$count <- reads$count + 1L
    original_read(edge)
  }
  plan <- crossform:::.tomography_resource_plan(
    fixture$form, fixture$frame, fixture$frame
  )
  expect_gt(plan$planned_workspace_bytes, 1)
  expect_silent(crossform:::.validate_tomography_resource_plan(plan))
  expect_error(crossform:::.reconstruct_neural_evidence(
    fixture$form, fixture$frame, workspace_bytes = 1
  ), "exceeding")
  expect_identical(reads$count, 0L)
})
