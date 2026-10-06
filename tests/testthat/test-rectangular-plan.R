rectangular_fixture <- function(domain_id = "rectangular-plan-domain") {
  set.seed(30817)
  domain <- abstract_domain(6L, id = domain_id)
  encoding_effects <- c("enc_face", "enc_house", "enc_tool")
  retrieval_effects <- c("ret_face", "ret_house")
  encoding_sources <- list(
    enc_run1 = matrix(rnorm(3L * 6L), 3L, 6L,
      dimnames = list(encoding_effects, NULL)),
    enc_run2 = matrix(rnorm(3L * 6L), 3L, 6L,
      dimnames = list(encoding_effects, NULL))
  )
  retrieval_sources <- list(
    ret_run1 = matrix(rnorm(2L * 6L), 2L, 6L,
      dimnames = list(retrieval_effects, NULL)),
    ret_run2 = matrix(rnorm(2L * 6L), 2L, 6L,
      dimnames = list(retrieval_effects, NULL))
  )
  encoding <- relation(
    encoding_sources,
    effects = effect_space(encoding_effects, basis_id = "rect:encoding"),
    domain = domain
  )
  retrieval <- relation(
    retrieval_sources,
    effects = effect_space(retrieval_effects, basis_id = "rect:retrieval"),
    domain = domain
  )
  over <- pairing(
    c("enc_run1", "enc_run2"), c("ret_run2", "ret_run1"),
    directed = TRUE, independence = "independent",
    generalizes_over = "run"
  )
  list(
    domain = domain, encoding = encoding, retrieval = retrieval,
    over = over, frame = compile_frame(voxelwise(), domain),
    encoding_sources = encoding_sources,
    retrieval_sources = retrieval_sources
  )
}

test_that("a rectangular plan compiles with distinct experimental axes", {
  fixture <- rectangular_fixture()
  plan <- plan_geometry(
    fixture$encoding, fixture$frame, fixture$over,
    right = fixture$retrieval
  )
  expect_s3_class(plan, "effect_geometry_plan")
  expect_identical(plan$codec, "rectangular")
  expect_identical(plan$logical_shape, c(3L, 2L))
  expect_identical(plan$packed_width, 6L)
  expect_false(plan$task$same_relation)
})

test_that("rectangular pair queries execute query-first and match the oracle", {
  fixture <- rectangular_fixture()
  plan <- plan_geometry(
    fixture$encoding, fixture$frame, fixture$over,
    right = fixture$retrieval
  )
  operator <- matrix(0, 3L, 2L, dimnames = list(
    c("enc_face", "enc_house", "enc_tool"), c("ret_face", "ret_house")
  ))
  operator["enc_face", "ret_face"] <- 1
  operator["enc_house", "ret_house"] <- 1
  query <- pair_query(
    operator,
    fixture$encoding$effect_space,
    fixture$retrieval$effect_space
  )
  view <- evaluate_geometry(plan, query = query)
  expect_s3_class(view, "effect_view")

  # Direct oracle: per voxel, the weighted ordered edge sum of
  # tr(H' B_enc diag(e_v) B_ret').
  oracle <- vapply(seq_len(6L), function(voxel) {
    edge_values <- c(
      sum(operator * tcrossprod(
        fixture$encoding_sources$enc_run1[, voxel],
        fixture$retrieval_sources$ret_run2[, voxel]
      )),
      sum(operator * tcrossprod(
        fixture$encoding_sources$enc_run2[, voxel],
        fixture$retrieval_sources$ret_run1[, voxel]
      ))
    )
    mean(edge_values)
  }, numeric(1))
  expect_equal(drop(view$values), oracle, tolerance = 1e-12,
    ignore_attr = TRUE)
})

test_that("a rectangular plan materializes to a queryable rectangular form", {
  fixture <- rectangular_fixture()
  plan <- plan_geometry(
    fixture$encoding, fixture$frame, fixture$over,
    right = fixture$retrieval
  )
  form <- materialize_geometry(plan)
  expect_s3_class(form, "effect_form")
  expect_false(inherits(form, "effect_geometry"))
  expect_identical(form$codec, "rectangular")
  expect_false(form$capabilities$self_form)
  expect_identical(form$logical_shape, c(3L, 2L))

  operator <- matrix(rnorm(6L), 3L, 2L)
  query <- pair_query(
    operator,
    fixture$encoding$effect_space,
    fixture$retrieval$effect_space
  )
  direct <- evaluate_geometry(plan, query = query)
  projected <- query_geometry(form, query)
  expect_equal(direct$values, projected$values, tolerance = 1e-12)

  # Check each component against an independent first-principles oracle.
  # `configuration` is computed as `total - coherent`, so their sum is not
  # evidence about any of the three.
  total <- geometry_component(form, "total")
  coherent <- geometry_component(form, "coherent")
  configuration <- geometry_component(form, "configuration")
  oracle <- geometry_component_oracle(
    relation_values = fixture$encoding_sources,
    right_values = fixture$retrieval_sources,
    frame_weights = fixture$frame$weights,
    partition_edges = fixture$over,
    symmetrize = FALSE
  )
  pack <- function(component) {
    do.call(rbind, lapply(oracle[[component]], as.vector))
  }
  expect_equal(total, pack("total"), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(coherent, pack("coherent"), tolerance = 1e-12,
    ignore_attr = TRUE)
  expect_equal(configuration, pack("configuration"), tolerance = 1e-12,
    ignore_attr = TRUE)
})

test_that("rectangular plans refuse what their contract does not cover", {
  fixture <- rectangular_fixture()
  metric_refusal <- catch_refusal(plan_geometry(
    fixture$encoding, fixture$frame, fixture$over,
    metric = neural_metric(diag(6L), fixture$domain),
    right = fixture$retrieval
  ))
  expect_s3_class(metric_refusal, "effect_capability_refusal")
  expect_identical(metric_refusal$capability, "rectangular_fixed_metric")

  expect_error(
    plan_geometry(
      fixture$encoding, fixture$frame,
      pairing(c("enc_run1", "enc_run2"), c("ret_run2", "ret_run1"),
        independence = "independent"),
      right = fixture$retrieval
    ),
    "ordered endpoints"
  , class = "effect_input_error")

  plan <- plan_geometry(
    fixture$encoding, fixture$frame, fixture$over,
    right = fixture$retrieval
  )
  rdm_refusal <- catch_refusal(rdm(plan))
  expect_s3_class(rdm_refusal, "effect_capability_refusal")
  expect_identical(rdm_refusal$capability, "symmetric_self_form")
  expect_identical(rdm_refusal$namespace, "geometry_views")
  expect_identical(rdm_refusal$reasons, "rectangular_cross_axis_plan")
  expect_match(conditionMessage(rdm_refusal),
    "symmetric self form|self-form")
  expect_error(
    evaluate_geometry(plan, query = bilinear_query(diag(3L))),
    "pair_query"
  , class = "effect_input_error")
})

test_that("rectangular compiler memory plans size q_left * q_right coordinates", {
  fixture <- rectangular_fixture()
  left <- fixture$encoding
  # A wider right side (q_right = 3 > (q_left + 1) / 2) must cost more than
  # the packed symmetric q(q+1)/2 width.
  wide_right <- relation(
    fixture$encoding_sources,
    effects = effect_space(c("enc_face", "enc_house", "enc_tool"),
      basis_id = "rect:wide"),
    domain = fixture$domain
  )
  requirements <- crossform:::.component_requirements(NULL, "total")
  plan_for <- function(right_relation, width) {
    crossform:::.compiler_memory_plan(left, fixture$frame, compute_policy(),
      6L, 6L, 6L, width, "memory", requirements,
      right_relation = right_relation)
  }
  m <- nrow(fixture$frame$weights)
  symmetric <- plan_for(NULL, 6L)
  rectangular <- plan_for(wide_right, 9L)
  expect_equal(symmetric$categories[["atom_block"]], 6 * 6 * 8)
  expect_equal(rectangular$categories[["atom_block"]], 6 * 9 * 8)
  expect_equal(rectangular$categories[["local_state"]],
    (m * 3 * 2 + m * 3 * 2) * 8)
  expect_gt(rectangular$planned_workspace_bytes,
    symmetric$planned_workspace_bytes)
  narrow <- plan_for(fixture$retrieval, 6L)
  expect_equal(narrow$categories[["local_state"]],
    (m * 3 * 2 + m * 2 * 2) * 8)
})
