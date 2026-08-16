joint_measurement_fixture <- function(seed = 2026081215,
                                      crossvalidated = FALSE) {
  set.seed(seed)
  q <- 5L
  p <- 4L
  domain <- abstract_domain(p, id = "joint:neural:v1")
  effects <- effect_space(paste0("sample", seq_len(q)),
    basis_id = "joint:samples:v1")
  values <- list(
    run1 = matrix(rnorm(q * p), q, p),
    run2 = matrix(rnorm(q * p), q, p)
  )
  rel <- relation(values, effects = effects, domain = domain)
  h <- diag(q) - matrix(1 / q, q, q)
  query <- crossform:::.variation_pair_query(
    h, effects, sampling_axis = "trial",
    construction = "joint_covariance",
    provenance = list(estimator = "centered-within-run")
  )
  legs <- list(
    a = crossform:::.measurement_leg(
      matrix(c(1, 0, 1, -1, 0, 2, 1, 0), 2, p, byrow = TRUE),
      domain,
      crossform:::.measurement_axis(c("a1", "a2"), "joint:a:v1")
    ),
    b = crossform:::.measurement_leg(
      matrix(c(0, 1, 1, 0, 2, -1, 0, 1), 2, p, byrow = TRUE),
      domain,
      crossform:::.measurement_axis(c("b1", "b2"), "joint:b:v1")
    )
  )
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
    crossform:::.evidence_stage_plan(),
    crossform:::.evidence_materialization(
      "measurement_form", "complete_form"
    )
  )
  list(task = task, values = values, frame = frame)
}

test_that("measurement_form is complete, durable, and distinct from effect_form", {
  fixture <- measurement_kernel_fixture(seed = 2026081216)
  run <- crossform:::.run_measurement_contraction(
    fixture$task, compute_policy(block_features = 2L), route = "pull_h"
  )
  form <- crossform:::.measurement_form_from_contraction(
    fixture$task, run
  )

  expect_s3_class(form, "effect_measurement_form")
  expect_false(inherits(form, "effect_form"))
  expect_identical(form$result_capability, "complete_form")
  expect_identical(form$completeness, "complete")
  expect_identical(form$edge_completeness, "requested_complete")
  expect_identical(form$plan$origin_task_id, fixture$task$task_id)
  expect_identical(form$plan$source_uses, fixture$task$source_uses)
  expect_identical(form$receipt$scientific_plan_id,
    form$plan$scientific_plan_id)
  expect_identical(form$receipt$execution$completion_status, "complete")
  expect_identical(form$receipt$execution$kernel_version,
    "measurement-form-v1")
  for (edge in seq_along(run$blocks)) {
    expect_equal(crossform:::.measurement_block(form, edge),
      run$blocks[[edge]], tolerance = 0)
  }
  expect_silent(crossform:::.validate_measurement_form(form))
})

test_that("memory and block codecs preserve one semantic result identity", {
  fixture <- measurement_kernel_fixture(seed = 2026081217)
  run <- crossform:::.run_measurement_contraction(
    fixture$task, compute_policy(block_features = 2L), route = "pull_h"
  )
  path <- tempfile(fileext = ".emf")
  on.exit(unlink(c(path, paste0(path, ".manifest.rds"))), add = TRUE)
  memory <- crossform:::.measurement_form_from_contraction(
    fixture$task, run, storage = "memory"
  )
  blocked <- crossform:::.measurement_form_from_contraction(
    fixture$task, run, storage = "block", path = path
  )

  expect_identical(memory$plan, blocked$plan)
  expect_identical(memory$capabilities, blocked$capabilities)
  expect_identical(memory$contract_signature, blocked$contract_signature)
  expect_identical(memory$codec, blocked$codec)
  expect_false(identical(memory$storage, blocked$storage))
  for (edge in seq_len(nrow(memory$block_index))) {
    expect_equal(crossform:::.measurement_block(memory, edge),
      crossform:::.measurement_block(blocked, edge), tolerance = 0)
  }
})

test_that("construction capabilities cannot be forged from diagnostics", {
  fixture <- measurement_kernel_fixture(seed = 2026081218)
  run <- crossform:::.run_measurement_contraction(fixture$task,
    route = "pull_h")
  form <- crossform:::.measurement_form_from_contraction(fixture$task, run)

  expect_false(form$capabilities$repeated_variation)
  expect_false(form$capabilities$joint_covariance)
  expect_false(form$capabilities$guaranteed_psd)
  expect_true(is.integer(form$diagnostics$blocks$effective_rank))
  expect_true(is.numeric(form$diagnostics$blocks$observed_min_eigenvalue))

  forged <- form
  forged$capabilities$guaranteed_psd <- TRUE
  expect_error(crossform:::.validate_measurement_form(forged),
    "cannot be forged")
  forged <- form
  forged$diagnostics$blocks$observed_min_eigenvalue[] <- 1
  expect_error(crossform:::.validate_measurement_form(forged),
    "diagnostics signature")
})

test_that("joint covariance requires same-partition self-products", {
  valid <- joint_measurement_fixture(crossvalidated = FALSE)
  run <- crossform:::.run_measurement_contraction(
    valid$task, compute_policy(block_features = 2L), route = "pull_h"
  )
  form <- crossform:::.measurement_form_from_contraction(
    valid$task, run,
    query_construction = "joint_covariance",
    edge_scope = "frame_complete"
  )

  expect_true(form$capabilities$self_form)
  expect_true(form$capabilities$symmetric)
  expect_true(form$capabilities$repeated_variation)
  expect_true(form$capabilities$joint_covariance)
  expect_true(form$capabilities$complete_edge_set)
  expect_true(form$capabilities$guaranteed_psd)
  expect_identical(form$capabilities$sampling_axis, "trial")

  blocks <- lapply(seq_len(nrow(form$block_index)), function(edge) {
    crossform:::.measurement_block(form, edge)
  })
  rows <- valid$frame$node_ids
  assembled <- do.call(rbind, lapply(rows, function(left) {
    do.call(cbind, lapply(rows, function(right) {
      position <- which(form$block_index$left == left &
        form$block_index$right == right)
      blocks[[position]]
    }))
  }))
  expect_gte(min(eigen(assembled, symmetric = TRUE,
    only.values = TRUE)$values), -1e-10)

  invalid <- joint_measurement_fixture(crossvalidated = TRUE)
  invalid_run <- crossform:::.run_measurement_contraction(
    invalid$task, route = "pull_h"
  )
  expect_error(crossform:::.measurement_form_from_contraction(
    invalid$task, invalid_run,
    query_construction = "joint_covariance",
    edge_scope = "frame_complete"
  ), "same-partition|crossvalidated")
  variation <- crossform:::.measurement_form_from_contraction(
    invalid$task, invalid_run,
    query_construction = "psd_variation",
    edge_scope = "frame_complete"
  )
  expect_true(variation$capabilities$repeated_variation)
  expect_false(variation$capabilities$joint_covariance)
  expect_false(variation$capabilities$guaranteed_psd)
})

test_that("reversal transposes blocks and swaps every measurement axis", {
  fixture <- measurement_kernel_fixture(seed = 2026081219)
  run <- crossform:::.run_measurement_contraction(fixture$task,
    route = "pull_h")
  form <- crossform:::.measurement_form_from_contraction(fixture$task, run)
  reversed <- crossform:::.reverse_measurement_form(form)

  expect_identical(reversed$left_frame, form$right_frame)
  expect_identical(reversed$right_frame, form$left_frame)
  expect_identical(reversed$block_index$left, form$block_index$right)
  expect_identical(reversed$block_index$right, form$block_index$left)
  expect_identical(reversed$receipt$derivation$kind, "reversal")
  expect_identical(reversed$receipt$derivation$parent,
    form$receipt$signature)
  for (edge in seq_len(nrow(form$block_index))) {
    expect_equal(crossform:::.measurement_block(reversed, edge),
      t(crossform:::.measurement_block(form, edge)), tolerance = 0)
  }
  expect_silent(crossform:::.validate_measurement_form(reversed))
})

test_that("swapped axes, incomplete stores, and forged plan identity fail", {
  fixture <- measurement_kernel_fixture(seed = 2026081220)
  run <- crossform:::.run_measurement_contraction(fixture$task,
    route = "pull_h")
  form <- crossform:::.measurement_form_from_contraction(fixture$task, run)

  swapped <- form
  swapped$left_frame <- form$right_frame
  expect_error(crossform:::.validate_measurement_form(swapped, probe = FALSE),
    "axes|identity|disagree")
  forged <- form
  forged$plan$edge_scope <- "frame_complete"
  expect_error(crossform:::.validate_measurement_form(forged, probe = FALSE),
    "complete|identity")

  path <- tempfile(fileext = ".emf")
  on.exit(unlink(c(path, paste0(path, ".manifest.rds"))), add = TRUE)
  incomplete <- crossform:::.file_measurement_store(
    path, form$block_index, create = TRUE
  )
  expect_error(crossform:::.new_measurement_form(
    incomplete, form$plan, form$capabilities, form$diagnostics, form$receipt
  ), "incomplete")
})

test_that("query-only measurement views cannot become complete forms", {
  fixture <- measurement_kernel_fixture(seed = 2026081221)
  run <- crossform:::.run_measurement_contraction(fixture$task,
    route = "pull_h")
  form <- crossform:::.measurement_form_from_contraction(fixture$task, run)
  values <- matrix(vapply(seq_len(nrow(form$block_index)), function(edge) {
    sum(crossform:::.measurement_block(form, edge)^2)
  }, numeric(1)), ncol = 1L)
  view <- crossform:::.measurement_view(
    values, form$plan, form$receipt,
    view = list(kind = "squared_energy")
  )

  expect_s3_class(view, "effect_measurement_view")
  expect_false(inherits(view, "effect_measurement_form"))
  expect_identical(view$completeness, "query_only")
  expect_identical(view$edge_completeness, "requested_complete")
  expect_silent(crossform:::.validate_measurement_view(view))
  upgraded <- view
  upgraded$completeness <- "complete"
  expect_error(crossform:::.validate_measurement_view(upgraded),
    "query-only")
  expect_error(crossform:::.measurement_block(view, 1L),
    "canonical complete")
})

test_that("typed view eligibility is capability-gated before block reads", {
  fixture <- measurement_kernel_fixture(seed = 2026081223)
  run <- crossform:::.run_measurement_contraction(fixture$task,
    route = "pull_h")
  arbitrary <- crossform:::.measurement_form_from_contraction(
    fixture$task, run
  )

  expect_silent(crossform:::.require_measurement_view_capabilities(
    arbitrary, "effect_coupling"
  ))
  # Capability gating is a classed refusal, so a caller can branch on the
  # cause rather than on the message prose.
  refusal <- catch_refusal(crossform:::.require_measurement_view_capabilities(
    arbitrary, "covariance_coupling"
  ))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "certified_repeated_variation")
  expect_identical(refusal$namespace, "coupling_views")
  expect_identical(refusal$reasons, "repeated_variation_not_certified")
  expect_match(conditionMessage(refusal), "repeated variation")

  joint <- joint_measurement_fixture(crossvalidated = FALSE)
  joint_run <- crossform:::.run_measurement_contraction(joint$task,
    route = "pull_h")
  covariance <- crossform:::.measurement_form_from_contraction(
    joint$task, joint_run,
    query_construction = "joint_covariance",
    edge_scope = "frame_complete"
  )
  expect_silent(crossform:::.require_measurement_view_capabilities(
    covariance, "covariance_coupling"
  ))
  unvalidated <- catch_refusal(
    crossform:::.require_measurement_view_capabilities(
      covariance, "canonical_coupling", positive_self_blocks = FALSE
    )
  )
  expect_s3_class(unvalidated, "effect_capability_refusal")
  expect_identical(unvalidated$capability, "coherent_joint_covariance")
  expect_identical(unvalidated$namespace, "coupling_views")
  # The joint covariance is certified here, so only the self-block reason is
  # unmet: the refusal reports every unmet reason, not the first one.
  expect_identical(unvalidated$reasons, "self_blocks_not_validated")
  expect_match(conditionMessage(unvalidated), "positive self-blocks")
  expect_silent(crossform:::.require_measurement_view_capabilities(
    covariance, "canonical_coupling", positive_self_blocks = TRUE
  ))
})

test_that("regularization identity is recorded but raw blocks cannot claim it applied", {
  fixture <- measurement_kernel_fixture(seed = 2026081222)
  run <- crossform:::.run_measurement_contraction(fixture$task,
    route = "pull_h")
  policy <- crossform:::.measurement_regularization(
    "ridge", lambda_left = 1e-6, lambda_right = 2e-6,
    applied = FALSE
  )
  form <- crossform:::.measurement_form_from_contraction(
    fixture$task, run, regularization = policy
  )

  expect_identical(form$plan$regularization, policy)
  expect_identical(form$diagnostics$regularization$policy, policy$signature)
  expect_false(form$diagnostics$regularization$applied)
  expect_error(crossform:::.measurement_form_from_contraction(
    fixture$task, run,
    regularization = crossform:::.measurement_regularization(
      "ridge", 1e-6, applied = TRUE
    )
  ), "Raw measurement blocks")
})
