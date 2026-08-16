coupling_test_fixture <- function(multivariate = FALSE, rank_one = FALSE,
                                  crossvalidated = FALSE,
                                  reducer = reduce_partitions(),
                                  zero_second = FALSE,
                                  partition_subset = NULL,
                                  seed = 2026081225) {
  set.seed(seed)
  q <- 8L
  p <- 4L
  domain <- abstract_domain(p, id = "coupling:neural:v1")
  effects <- effect_space(paste0("sample", seq_len(q)),
    basis_id = "coupling:samples:v1")
  first <- matrix(rnorm(q * p), q, p)
  second <- matrix(rnorm(q * p), q, p)
  first[, 2L] <- -0.8 * first[, 1L] + 0.45 * first[, 2L]
  second[, 2L] <- -0.7 * second[, 1L] + 0.55 * second[, 2L]
  if (zero_second) {
    first[, 2L] <- 0
    second[, 2L] <- 0
  }
  values <- list(run1 = first, run2 = second)
  if (!is.null(partition_subset)) {
    values <- values[partition_subset]
  }
  rel <- relation(values, effects = effects, domain = domain)
  center <- diag(q) - matrix(1 / q, q, q)
  h <- if (rank_one) {
    contrast <- c(-1, -1, -1, -1, 1, 1, 1, 1)
    tcrossprod(contrast)
  } else {
    center / (q - 1)
  }
  query <- crossform:::.variation_pair_query(
    h, effects, sampling_axis = "trial",
    construction = "joint_covariance",
    provenance = list(estimator = if (rank_one) {
      "single-effect-direction"
    } else {
      "centered-within-partition"
    })
  )
  operators <- if (multivariate) {
    list(
      a = diag(p)[1:2, , drop = FALSE],
      b = diag(p)[3:4, , drop = FALSE]
    )
  } else {
    list(
      a = matrix(c(1, 0, 0, 0), 1L),
      b = matrix(c(0, 1, 0, 0), 1L)
    )
  }
  legs <- Map(function(operator, node) {
    crossform:::.measurement_leg(
      operator, domain,
      crossform:::.measurement_axis(
        paste0(node, seq_len(nrow(operator))),
        paste0("coupling:", node, ":v1"),
        basis_id = paste0("coupling:", node, ":fixed-basis:v1")
      )
    )
  }, operators, names(operators))
  names(legs) <- names(operators)
  frame <- crossform:::.measurement_frame(legs)
  pairs <- expand.grid(
    left = frame$node_ids,
    right = frame$node_ids,
    stringsAsFactors = FALSE
  )
  spatial_edges <- crossform:::.measurement_edges(
    pairs$left, pairs$right, frame
  )
  over <- if (crossvalidated) {
    cross_partitions(rel)
  } else {
    pairing(
      rel$partitions, rel$partitions,
      directed = TRUE,
      self_pairs = "allow_biased",
      independence = "not_independent"
    )
  }
  partition_edges <- crossform:::.ordered_partition_edges(
    over, rel$partitions, rel$partitions, TRUE
  )
  task <- crossform:::.new_evidence_task(
    rel, rel, TRUE, partition_edges,
    crossform:::.closed_experimental_boundary(
      query, role = "variation", sampling_axis = "trial"
    ),
    crossform:::.open_neural_boundary(frame, frame, spatial_edges),
    crossform:::.evidence_stage_plan(reducer = reducer),
    crossform:::.evidence_materialization(
      "measurement_form", "complete_form"
    )
  )
  run <- crossform:::.run_measurement_contraction(task, route = "pull_h")
  construction <- if (crossvalidated) "psd_variation" else
    "joint_covariance"
  form <- crossform:::.measurement_form_from_contraction(
    task, run,
    query_construction = construction,
    edge_scope = "frame_complete"
  )
  list(
    form = form,
    values = values,
    h = h,
    frame = frame,
    operators = operators,
    partition_edges = partition_edges
  )
}

coupling_expected_neural_form <- function(fixture) {
  Reduce(`+`, lapply(seq_len(nrow(fixture$partition_edges)), function(i) {
    row <- fixture$partition_edges[i, , drop = FALSE]
    left <- fixture$values[[row$left[[1L]]]]
    right <- fixture$values[[row$right[[1L]]]]
    row$weight[[1L]] * crossprod(left, fixture$h %*% right)
  }))
}

coupling_edge_position <- function(form, left, right) {
  position <- which(form$block_index$left == left &
    form$block_index$right == right)
  stopifnot(length(position) == 1L)
  position
}

test_that("effect coupling is the unrestricted algebraic measurement form", {
  fixture <- coupling_test_fixture(rank_one = TRUE)
  result <- crossform:::.effect_coupling(fixture$form)
  expected_q <- coupling_expected_neural_form(fixture)

  expect_s3_class(result, "effect_coupling_result")
  expect_identical(result$kind, "effect_coupling")
  expect_identical(result$normalization_axis, "none")
  expect_identical(result$summary_axis, "measurement_coordinates")
  expect_match(result$terminology, "no_covariance_claim")
  for (edge in seq_len(nrow(fixture$form$block_index))) {
    row <- fixture$form$block_index[edge, , drop = FALSE]
    expected <- fixture$operators[[row$left[[1L]]]] %*% expected_q %*%
      t(fixture$operators[[row$right[[1L]]]])
    expect_equal(result$values[[edge]], expected, tolerance = 1e-12)
  }
  expect_silent(crossform:::.validate_coupling_result(result))
  expect_error(crossform:::.pearson_coupling(fixture$form),
    "rank above one|rank-one effect direction")
})

test_that("scalar repeated variation recovers signed Pearson correlation", {
  fixture <- coupling_test_fixture()
  result <- crossform:::.pearson_coupling(fixture$form)
  q_form <- coupling_expected_neural_form(fixture)
  expected <- q_form[1L, 2L] / sqrt(q_form[1L, 1L] * q_form[2L, 2L])
  position <- coupling_edge_position(fixture$form, "a", "b")
  got <- result$values$correlation[
    result$values$edge_id == fixture$form$block_index$edge_id[[position]]
  ]

  expect_lt(expected, 0)
  expect_equal(got, expected, tolerance = 1e-12)
  expect_identical(result$normalization_axis, "experimental_samples")
  expect_identical(result$units, "correlation")
  expect_identical(
    crossform:::.connectivity_view(fixture$form, "correlation"), result
  )
  source_end <- length(fixture$form$plan$stages$order)
  expect_identical(result$stage_order[seq_len(source_end)],
    fixture$form$plan$stages$order)
  expect_gt(match("experimental_sample_normalization", result$stage_order),
    match("partition_reduction", result$stage_order))
  expect_silent(crossform:::.validate_coupling_result(result))
})

test_that("normalized coupling rejects missing joint covariance and zero self variance", {
  crossvalidated <- coupling_test_fixture(crossvalidated = TRUE)
  expect_silent(crossform:::.effect_coupling(crossvalidated$form))
  expect_error(crossform:::.covariance_coupling(crossvalidated$form),
    "joint covariance|positive self-blocks")
  expect_error(crossform:::.canonical_coupling(
    crossvalidated$form,
    crossform:::.measurement_regularization("ridge", 1e-4)
  ), "joint covariance|positive self-blocks")

  zero <- coupling_test_fixture(zero_second = TRUE)
  expect_error(crossform:::.pearson_coupling(zero$form),
    "strictly positive scalar self-variance")
})

test_that("CCA, geometry alignment, and Gaussian information match references", {
  fixture <- coupling_test_fixture(multivariate = TRUE)
  regularization <- crossform:::.measurement_regularization(
    "ridge", lambda_left = 0.05, lambda_right = 0.08
  )
  edge <- coupling_edge_position(fixture$form, "a", "b")
  edge_id <- fixture$form$block_index$edge_id[[edge]]
  cross <- crossform:::.measurement_block(fixture$form, edge)
  left_self <- crossform:::.measurement_block(
    fixture$form, coupling_edge_position(fixture$form, "a", "a")
  )
  right_self <- crossform:::.measurement_block(
    fixture$form, coupling_edge_position(fixture$form, "b", "b")
  )
  inverse_root <- function(value) {
    eig <- eigen((value + t(value)) / 2, symmetric = TRUE)
    eig$vectors %*% diag(1 / sqrt(eig$values), length(eig$values)) %*%
      t(eig$vectors)
  }
  normalized <- inverse_root(left_self + diag(0.05, nrow(left_self))) %*%
    cross %*% inverse_root(right_self + diag(0.08, nrow(right_self)))
  expected_rho <- svd(normalized, nu = 0L, nv = 0L)$d

  canonical <- crossform:::.canonical_coupling(
    fixture$form, regularization
  )
  got_rho <- canonical$values$canonical_correlation[
    canonical$values$edge_id == edge_id
  ]
  expect_equal(got_rho, expected_rho, tolerance = 1e-11)
  expect_identical(canonical$normalization_axis, "experimental_samples")

  alignment <- crossform:::.geometry_alignment(fixture$form)
  expected_alignment <- sum(cross^2) /
    sqrt(sum(left_self^2) * sum(right_self^2))
  got_alignment <- alignment$values$geometry_alignment[
    alignment$values$edge_id == edge_id
  ]
  expect_equal(got_alignment, expected_alignment, tolerance = 1e-12)
  expect_identical(alignment$normalization_axis, "form_entries")
  expect_match(alignment$terminology,
    "not_dynamic_informational_connectivity")

  model <- crossform:::.gaussian_covariance_model(
    list(assumption = "joint Gaussian observations")
  )
  information <- crossform:::.gaussian_information(
    fixture$form, regularization, model, units = "bits"
  )
  expected_information <- -0.5 * sum(log1p(-(expected_rho^2))) / log(2)
  got_information <- information$values$information[
    information$values$edge_id == edge_id
  ]
  expect_equal(got_information, expected_information, tolerance = 1e-11)
  expect_identical(information$units, "bits")
  expect_silent(crossform:::.validate_coupling_result(canonical))
  expect_silent(crossform:::.validate_coupling_result(alignment))
  expect_silent(crossform:::.validate_coupling_result(information))
})

test_that("stage order and axis identity are part of the coupling estimand", {
  edge_first <- coupling_test_fixture(reducer = reduce_partitions())
  aggregate <- coupling_test_fixture(reducer = aggregate_first())
  edge_result <- crossform:::.pearson_coupling(edge_first$form)
  aggregate_result <- crossform:::.pearson_coupling(aggregate$form)

  expect_false(identical(edge_first$form$plan$scientific_plan_id,
    aggregate$form$plan$scientific_plan_id))
  expect_false(identical(edge_result$signature, aggregate_result$signature))
  expect_false(identical(edge_result$stage_order,
    aggregate_result$stage_order))
  expect_identical(
    edge_first$form$plan$stages$normalization$axis, "neural_features"
  )
  expect_identical(edge_result$normalization_axis, "experimental_samples")
})

test_that("partition-pair and aggregate-first normalization are distinct estimands", {
  first <- coupling_test_fixture(partition_subset = "run1")
  second <- coupling_test_fixture(partition_subset = "run2")
  forms <- list(first$form, second$form)
  weights <- c(0.3, 0.7)
  transform <- fisher_z(boundary = "clip", delta = 1e-8)
  edge_first <- crossform:::.partitioned_pearson_coupling(
    forms, weights, "within_partition_pair", transform
  )
  aggregate_first_result <- crossform:::.partitioned_pearson_coupling(
    forms, weights, "after_partition_aggregation", transform
  )
  position <- coupling_edge_position(first$form, "a", "b")
  edge_id <- first$form$block_index$edge_id[[position]]
  q1 <- coupling_expected_neural_form(first)
  q2 <- coupling_expected_neural_form(second)
  correlation <- function(value) {
    value[1L, 2L] / sqrt(value[1L, 1L] * value[2L, 2L])
  }
  expected_edge <- sum(weights * atanh(c(correlation(q1), correlation(q2))))
  combined <- weights[[1L]] * q1 + weights[[2L]] * q2
  expected_aggregate <- atanh(correlation(combined))
  got_edge <- edge_first$values$value[edge_first$values$edge_id == edge_id]
  got_aggregate <- aggregate_first_result$values$value[
    aggregate_first_result$values$edge_id == edge_id
  ]

  expect_equal(got_edge, expected_edge, tolerance = 1e-12)
  expect_equal(got_aggregate, expected_aggregate, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(got_edge, got_aggregate, tolerance = 1e-8)))
  expect_false(identical(edge_first$signature,
    aggregate_first_result$signature))
  expect_identical(edge_first$partition_policy$placement,
    "within_partition_pair")
  expect_identical(aggregate_first_result$partition_policy$placement,
    "after_partition_aggregation")
  expect_lt(match("experimental_sample_normalization",
    edge_first$stage_order), match("partition_reduction",
    edge_first$stage_order))
  expect_gt(match("experimental_sample_normalization",
    aggregate_first_result$stage_order), match("partition_reduction",
    aggregate_first_result$stage_order))
  expect_silent(crossform:::.validate_coupling_result(edge_first))
  expect_silent(crossform:::.validate_coupling_result(
    aggregate_first_result
  ))
})

test_that("connectivity convenience validates view-specific declarations", {
  fixture <- coupling_test_fixture(multivariate = TRUE)
  expect_error(crossform:::.connectivity_view(fixture$form, "canonical"),
    "regularization")
  expect_error(crossform:::.connectivity_view(
    fixture$form, "gaussian_information",
    regularization = crossform:::.measurement_regularization("ridge", 0.1)
  ), "Gaussian model")
  expect_false(exists("informational_connectivity",
    envir = asNamespace("crossform"), inherits = FALSE))
})

test_that("coupling() takes the adjoint closure from the plan vocabulary", {
  set.seed(41219)
  domain <- abstract_domain(4L, id = "coupling-plan-domain")
  relation <- relation(
    list(
      run1 = matrix(rnorm(8), 2, 4, dimnames = list(c("a", "b"), NULL)),
      run2 = matrix(rnorm(8), 2, 4, dimnames = list(c("a", "b"), NULL))
    ),
    effects = effect_space(c("a", "b")), domain = domain
  )
  plan <- plan_geometry(
    relation, compile_frame(regions(c("r1", "r1", "r2", "r2")), domain),
    cross_partitions(relation)
  )
  h <- diag(2) - 0.5
  query <- variation_query(
    h, relation$effect_space, "trial", "psd_variation"
  )
  form <- coupling(plan, cbind("r1", "r2"), by = query)
  effect <- effect_coupling(form)
  expect_s3_class(effect, "effect_coupling_result")

  # The adapter is one vocabulary, not a new estimator: the same edges
  # declared through the measurement pipeline agree exactly.
  frame <- measurement_frame(plan$frame)
  manual <- measurement_form(
    relation, edge_frame("r1", "r2", frame), query, plan$pairing
  )
  expect_equal(
    effect$values[[1L]], effect_coupling(manual)$values[[1L]],
    tolerance = 1e-12
  )

  expect_error(coupling(plan, cbind("r1", "nowhere"), by = query),
    "identify measurements")
  other <- relation(
    list(other1 = matrix(rnorm(4), 1, 4, dimnames = list("z", NULL))),
    effects = effect_space("z"), domain = domain
  )
  rectangular <- plan_geometry(
    relation, compile_frame(voxelwise(), domain),
    pairing(c("run1", "run2"), c("other1", "other1"), directed = TRUE,
      independence = "independent"),
    right = other
  )
  expect_error(coupling(rectangular, cbind(1, 2), by = query),
    "self-form plan")
})
