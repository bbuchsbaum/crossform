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
    over = over, frame = compile_frame(voxels(), domain),
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
  form <- geometry(plan)
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

  # Check the exact algebraic recomposition to floating-point tolerance.
  total <- geometry_component(form, "total")
  coherent <- geometry_component(form, "coherent")
  configuration <- geometry_component(form, "configuration")
  expect_equal(total, coherent + configuration, tolerance = 1e-12)
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
  )

  plan <- plan_geometry(
    fixture$encoding, fixture$frame, fixture$over,
    right = fixture$retrieval
  )
  expect_error(rdm(plan), "symmetric self form|self-form")
  expect_error(
    evaluate_geometry(plan, query = bilinear_query(diag(3L))),
    "pair_query"
  )
})
