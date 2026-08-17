evidence_ir_relation <- function(partitions, effects, domain, reads = NULL,
                                 revision_character = "a") {
  count <- length(effects$coordinates)
  if (is.null(reads)) {
    values <- stats::setNames(lapply(seq_along(partitions), function(index) {
      matrix(index, count, domain$reference$n_features)
    }), partitions)
    return(relation(values, effects = effects, domain = domain))
  }
  source <- function(features) {
    reads$count <- reads$count + 1L
    matrix(1, count, length(features))
  }
  revision <- paste0(
    "sha256:", paste(rep(revision_character, 64), collapse = "")
  )
  relation(
    stats::setNames(rep(list(source), length(partitions)), partitions),
    source_dims = rep(list(c(count, domain$reference$n_features)),
      length(partitions)),
    effects = effects,
    domain = domain,
    capabilities = source_capabilities(TRUE, stable_revision = revision)
  )
}

test_that("effect compilation is an exact adapter over the evidence-task IR", {
  domain <- abstract_domain(3, id = "neural:shared:v1")
  effects <- effect_space(c("a", "b"), basis_id = "conditions:v1")
  rel <- evidence_ir_relation(c("r1", "r2", "r3"), effects, domain)
  over <- cross_partitions(rel)

  evidence <- crossform:::.compile_effect_evidence_task(rel, over)
  legacy <- crossform:::.compile_effect_task(rel, over)

  expect_s3_class(evidence, "effect_evidence_task")
  expect_identical(crossform:::.as_compiled_effect_task(evidence), legacy)
  expect_identical(evidence$task_id, legacy$task_id)
  expect_identical(evidence$identity_schema, "effect-form-v1")
  expect_identical(evidence$experimental_boundary$state, "open")
  expect_identical(evidence$neural_boundary$state, "closed")
  expect_identical(evidence$materialization$kind, "effect_form")
  expect_identical(evidence$materialization$completeness, "complete_form")
  expect_identical(names(legacy), c(
    "left_relation", "right_relation", "same_relation", "left_relation_id",
    "right_relation_id", "left_space", "right_space", "ordered_edges",
    "bridge", "normalizer", "transform", "lowering", "reducer",
    "query_identity", "source_uses", "distinct_handle_keys", "task_id"
  ))
})

test_that("all four identified spaces and axis-labelled stages are explicit", {
  domain <- abstract_domain(4, id = "neural:shared:v1")
  left_space <- effect_space(c("encode1", "encode2"), basis_id = "encode:v1")
  right_space <- effect_space(c("retrieve1", "retrieve2", "retrieve3"),
    basis_id = "retrieve:v1")
  left <- evidence_ir_relation(c("e1", "e2"), left_space, domain)
  right <- evidence_ir_relation(c("r1", "r2"), right_space, domain)
  over <- pairing(c("e1", "e2"), c("r1", "r2"), c(1, 2), directed = TRUE)
  query <- pair_query(matrix(1:6, 2), left_space, right_space)
  task <- crossform:::.compile_effect_evidence_task(
    left, over, right, query,
    normalizer = correlation(),
    reducer = aggregate_first()
  )

  expect_identical(task$spaces, list(
    experimental_left = left_space,
    experimental_right = right_space,
    neural_left = domain$reference,
    neural_right = domain$reference
  ))
  expect_identical(task$stages$normalization$axis, "neural_features")
  expect_identical(task$stages$normalization$placement,
    "after_partition_aggregation")
  expect_identical(task$stages$reduction$axis, "partition_pairs")
  expect_identical(task$stages$transform$axis, "form_entries")
  expect_identical(task$materialization$completeness, "query_only")
  expect_identical(task$materialization$projection$query, query)
  expect_false(any(c("storage", "tile", "workers", "route") %in% names(task)))
})

test_that("legacy scientific identity is byte-for-byte the certified semantic", {
  domain <- abstract_domain(3, id = "neural:shared:v1")
  effects <- effect_space(c("a", "b"), basis_id = "conditions:v1")
  rel <- evidence_ir_relation(c("r1", "r2"), effects, domain)
  over <- cross_partitions(rel)
  task <- crossform:::.compile_effect_evidence_task(rel, over)
  expected_semantic <- crossform:::.effect_task_semantic(
    task$left_relation_id,
    task$right_relation_id,
    effects,
    effects,
    task$ordered_partition_products,
    task$neural_boundary$bridge,
    inner_product(),
    crossform:::.identity_edge_transform(),
    crossform:::.partition_reducer(),
    list(kind = "complete_form", query = NULL)
  )

  expect_identical(task$semantic, expected_semantic)
  expect_identical(task$task_id,
    crossform:::.effect_task_id(expected_semantic))
})

test_that("closed experimental and open neural boundaries form measurement IR", {
  domain <- abstract_domain(3, id = "neural:shared:v1")
  effects <- effect_space(c("a", "b", "c"), basis_id = "conditions:v1")
  rel <- evidence_ir_relation(c("r1", "r2"), effects, domain)
  partition_edges <- crossform:::.ordered_partition_edges(
    cross_partitions(rel), rel$partitions, rel$partitions, TRUE
  )
  additive <- additive_frame(diag(3), domain = domain)
  frame <- crossform:::.measurement_frame_from_additive(additive)
  spatial_edges <- crossform:::.measurement_edges(
    c("measurement_1", "measurement_2"),
    c("measurement_2", "measurement_3"),
    frame
  )
  query <- pair_query(diag(3), effects, effects)
  experimental <- crossform:::.closed_experimental_boundary(
    query, role = "variation", sampling_axis = "trial"
  )
  neural <- crossform:::.open_neural_boundary(frame, frame, spatial_edges)
  stages <- crossform:::.evidence_stage_plan()
  materialization <- crossform:::.evidence_materialization(
    "measurement_form", "complete_form"
  )
  task <- crossform:::.new_evidence_task(
    rel, rel, TRUE, partition_edges,
    experimental, neural, stages, materialization
  )

  expect_identical(task$identity_schema, "evidence-pairing-v1")
  expect_identical(task$experimental_boundary$state, "closed")
  expect_identical(task$neural_boundary$state, "open")
  expect_identical(task$materialization$kind, "measurement_form")
  expect_identical(nrow(task$neural_boundary$edges$edges), 2L)
  expect_silent(crossform:::.validate_evidence_task(task))
})

test_that("the neural boundary can close with an identified pair query", {
  left_domain <- abstract_domain(2, id = "neural:left:v1")
  right_domain <- abstract_domain(3, id = "neural:right:v1")
  left_space <- effect_space(c("l1", "l2"), basis_id = "left-effects:v1")
  right_space <- effect_space(c("r1", "r2"), basis_id = "right-effects:v1")
  left <- evidence_ir_relation(c("e1", "e2"), left_space, left_domain)
  right <- evidence_ir_relation(c("r1", "r2"), right_space, right_domain)
  over <- pairing(c("e1", "e2"), c("r1", "r2"), c(0.5, 0.5),
    directed = TRUE)
  partition_edges <- crossform:::.ordered_partition_edges(
    over, left$partitions, right$partitions, FALSE
  )
  k <- matrix(c(1, 0, -1, 2, 3, 1), 2, 3)
  neural_query <- crossform:::.neural_pair_query(
    k, left_domain, right_domain,
    provenance = list(bridge = "fixed-scientific-query")
  )
  task <- crossform:::.new_evidence_task(
    left, right, FALSE, partition_edges,
    crossform:::.open_experimental_boundary(left_space, right_space),
    crossform:::.closed_neural_query_boundary(neural_query),
    crossform:::.evidence_stage_plan(),
    crossform:::.evidence_materialization("effect_form", "complete_form")
  )

  expect_identical(task$neural_boundary$closure_kind, "query")
  expect_identical(task$neural_boundary$query$operator, k)
  reversed <- crossform:::.reverse_evidence_task(task)
  expect_identical(reversed$neural_boundary$query$operator, t(k))
  expect_identical(reversed$neural_boundary$query$left_domain,
    right_domain$reference)
  expect_identical(crossform:::.reverse_evidence_task(reversed), task)
})

test_that("typed boundary mismatches fail before any source read", {
  reads <- new.env(parent = emptyenv())
  reads$count <- 0L
  domain <- abstract_domain(3, id = "neural:shared:v1")
  other_domain <- abstract_domain(3, feature_ids = 3:1,
    id = "neural:shared:v1")
  effects <- effect_space(c("a", "b"), basis_id = "conditions:v1")
  other_effects <- effect_space(c("a", "b"), basis_id = "conditions:v2")
  rel <- evidence_ir_relation(c("r1", "r2"), effects, domain, reads)
  edges <- crossform:::.ordered_partition_edges(
    cross_partitions(rel), rel$partitions, rel$partitions, TRUE
  )
  bridge <- crossform:::.identity_measurement_bridge(rel, rel)
  neural <- crossform:::.closed_neural_boundary(bridge, rel, rel)
  stages <- crossform:::.evidence_stage_plan()
  materialization <- crossform:::.evidence_materialization("scalar_field",
    "query_only")

  expect_error(crossform:::.new_evidence_task(
    rel, rel, TRUE, edges,
    crossform:::.closed_experimental_boundary(
      pair_query(diag(2), other_effects, effects)
    ), neural, stages, materialization
  ), "Experimental boundary axes", class = "effect_contract_error")

  wrong_frame <- crossform:::.measurement_frame_from_additive(
    additive_frame(diag(3), domain = other_domain)
  )
  wrong_edges <- crossform:::.measurement_edges(
    "measurement_1", "measurement_2", wrong_frame
  )
  expect_error(crossform:::.new_evidence_task(
    rel, rel, TRUE, edges,
    crossform:::.closed_experimental_boundary(
      pair_query(diag(2), effects, effects)
    ),
    crossform:::.open_neural_boundary(
      wrong_frame, wrong_frame, wrong_edges
    ),
    stages,
    crossform:::.evidence_materialization(
      "measurement_form", "complete_form"
    )
  ), "Neural boundary axes", class = "effect_contract_error")
  expect_identical(reads$count, 0L)
})

test_that("reversing every typed side exactly matches reversed compilation", {
  domain <- abstract_domain(3, id = "neural:shared:v1")
  left_space <- effect_space(c("l1", "l2"), basis_id = "left:v1")
  right_space <- effect_space(c("r1", "r2", "r3"), basis_id = "right:v1")
  left <- evidence_ir_relation(c("e1", "e2"), left_space, domain)
  right <- evidence_ir_relation(c("r1", "r2"), right_space, domain)
  over <- pairing(c("e1", "e2"), c("r1", "r2"), c(0.4, 0.6),
    directed = TRUE)
  h <- matrix(1:6, 2, 3)
  forward <- crossform:::.compile_effect_evidence_task(
    left, over, right, pair_query(h, left_space, right_space)
  )
  reversed <- crossform:::.reverse_evidence_task(forward)
  independently_reversed <- crossform:::.compile_effect_evidence_task(
    right,
    pairing(over$right, over$left, over$weight, directed = TRUE),
    left,
    pair_query(t(h), right_space, left_space)
  )

  expect_identical(reversed, independently_reversed)
  expect_identical(crossform:::.reverse_evidence_task(reversed), forward)
  expect_identical(reversed$spaces$experimental_left,
    forward$spaces$experimental_right)
  expect_identical(reversed$ordered_partition_products$left,
    forward$ordered_partition_products$right)
})

test_that("materialization states are enforced rather than inferred by shape", {
  domain <- abstract_domain(2)
  effects <- effect_space(c("a", "b"))
  rel <- evidence_ir_relation(c("r1", "r2"), effects, domain)
  edges <- crossform:::.ordered_partition_edges(
    cross_partitions(rel), rel$partitions, rel$partitions, TRUE
  )
  experimental <- crossform:::.open_experimental_boundary(effects, effects)
  bridge <- crossform:::.identity_measurement_bridge(rel, rel)
  neural <- crossform:::.closed_neural_boundary(bridge, rel, rel)

  expect_error(crossform:::.new_evidence_task(
    rel, rel, TRUE, edges, experimental, neural,
    crossform:::.evidence_stage_plan(),
    crossform:::.evidence_materialization(
      "measurement_form", "complete_form"
    )
  ), "incompatible with its open boundaries", class = "effect_contract_error")
  expect_error(crossform:::.evidence_materialization(
    "scalar_field", "complete_form"
  ), "necessarily", class = "effect_input_error")
})

test_that("evidence-task validation detects semantic mutation", {
  domain <- abstract_domain(2)
  effects <- effect_space(c("a", "b"))
  rel <- evidence_ir_relation(c("r1", "r2"), effects, domain)
  task <- crossform:::.compile_effect_evidence_task(
    rel, cross_partitions(rel)
  )
  forged <- task
  forged$stages$normalization$axis <- "experimental_samples"

  expect_error(crossform:::.validate_evidence_task(forged),
    "stage identity", class = "effect_contract_error")
})
