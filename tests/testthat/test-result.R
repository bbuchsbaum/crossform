test_that("complete geometry and query-only view are distinct contracts", {
  geometry <- result_fixture()
  query <- matrix(c(1, 0, -1), ncol = 1)
  view <- query_geometry(geometry, query, component = "total")

  expect_s3_class(geometry, "effect_geometry")
  expect_s3_class(geometry$effect_space, "effect_space")
  expect_identical(geometry$completeness, "full")
  expect_s3_class(view, "effect_view")
  expect_identical(view$completeness, "query_only")
  expect_false(inherits(view, "effect_geometry"))
  expect_error(geometry_component(view), "complete effect_geometry",
    class = "effect_input_error")
})

test_that("geometry queries enforce experimental-space identity", {
  geometry <- result_fixture()
  incompatible <- bilinear_query(diag(2), effects = effect_space(
    c("a", "b"), basis_id = "different:basis"))

  expect_error(query_geometry(geometry, incompatible),
    "effect spaces are incompatible", class = "effect_contract_error")
})

test_that("full materialization and direct query agree", {
  geometry <- result_fixture()
  query <- matrix(c(2, -1, 0.5, 0, 1, 1), nrow = 3)

  materialized <- query_geometry(geometry, query, component = "configuration")
  direct <- effect_view(
    (geometry_component(geometry, "total") -
      geometry_component(geometry, "coherent")) %*% query,
    query = query,
    component = "configuration",
    receipt = geometry$receipt,
    index = geometry$index
  )
  colnames(direct$values) <- paste0("view", seq_len(ncol(direct$values)))

  expect_equal(materialized$values, direct$values, tolerance = 1e-14)
  expect_identical(materialized$completeness, "query_only")
  expect_s3_class(direct$effect_space, "effect_space")
})

test_that("bilinear operators compile only after their contracts are checked", {
  geometry <- result_fixture()
  operator <- matrix(c(2, -1, -1, 3), 2)
  packed <- matrix(oracle_svec(operator), ncol = 1)

  expect_equal(
    query_geometry(geometry, bilinear_query(operator))$values,
    query_geometry(geometry, packed)$values,
    tolerance = 1e-14
  )

  asymmetric <- bilinear_query(diag(2))
  asymmetric$operator[1, 2] <- 1
  expect_error(query_geometry(geometry, asymmetric), "symmetric",
    class = "effect_input_error")
  expect_error(
    query_geometry(geometry, bilinear_query(diag(3))),
    "experimental dimension"
  , class = "effect_contract_error")
})

test_that("invalid direct queries perform no store reads", {
  base <- result_fixture()
  total <- geometry_component(base, "total")
  coherent <- geometry_component(base, "coherent")
  reads <- new.env(parent = emptyenv())
  reads$total <- 0L
  reads$coherent <- 0L
  total_store <- crossform:::.block_geometry_store(dim(total), function(rows = NULL) {
    reads$total <- reads$total + 1L
    if (is.null(rows)) total else total[rows, , drop = FALSE]
  })
  coherent_store <- crossform:::.block_geometry_store(dim(coherent), function(rows = NULL) {
    reads$coherent <- reads$coherent + 1L
    if (is.null(rows)) coherent else coherent[rows, , drop = FALSE]
  })
  geometry <- effect_geometry(
    total_store, coherent_store, base$marginals,
    effects = base$effects, receipt = base$receipt, index = base$index
  )
  reads$total <- 0L
  reads$coherent <- 0L

  invalid <- list(
    function() query_geometry(geometry, c(1, 2, 3)),
    function() query_geometry(geometry, matrix(NA_real_, 3, 1)),
    function() query_geometry(geometry, matrix(numeric(), 3, 0)),
    function() query_geometry(geometry, matrix(1, 2, 1)),
    function() query_geometry(geometry, matrix(1, 3, 1), row_block = 0),
    function() query_geometry(geometry, matrix(1, 3, 1), component = "unknown")
  )
  for (attempt in invalid) {
    expect_error(attempt())
    expect_identical(reads$total, 0L)
    expect_identical(reads$coherent, 0L)
  }

  query <- matrix(c(1, -2, 0.5), ncol = 1)
  got <- query_geometry(geometry, query, component = "configuration", row_block = 1)
  expect_gt(reads$total, 0L)
  expect_gt(reads$coherent, 0L)
  expected <- (total - coherent) %*% query
  colnames(expected) <- "view1"
  expect_equal(got$values, expected, tolerance = 1e-14)
})

test_that("query labels survive bounded result allocation", {
  geometry <- result_fixture()
  query <- matrix(c(1, 0, 0, 0, 1, 0), 3, 2,
    dimnames = list(NULL, c("first", "second")))
  expect_identical(colnames(query_geometry(geometry, query)$values),
    c("first", "second"))
})

test_that("public geometry consumers reject mutated result contracts", {
  geometry <- result_fixture()

  mutated <- geometry
  mutated$completeness <- "query_only"
  expect_error(geometry_component(mutated), "canonical complete",
    class = "effect_input_error")

  mutated <- geometry
  mutated$effects <- rev(mutated$effects)
  expect_error(query_geometry(mutated, matrix(1, 3, 1)),
    "coordinate labels", class = "effect_contract_error")

  mutated <- geometry
  mutated$receipt$completion_status <- "invented"
  expect_error(rdm(mutated), "receipt|completion", class = "effect_input_error")

  mutated <- geometry
  mutated$storage <- "memory"
  mutated$total$representation <- "forged"
  expect_error(contrast_energy(mutated, c(a = 1, b = -1)),
    "storage metadata", class = "effect_contract_error")
})

test_that("storage representation cannot change semantic completeness", {
  geometry <- result_fixture()
  total <- geometry_component(geometry, "total")
  coherent <- geometry_component(geometry, "coherent")
  blocked_total <- crossform:::.block_geometry_store(
    dim(total), function(rows = NULL) if (is.null(rows)) total else total[rows, , drop = FALSE]
  )
  blocked_coherent <- crossform:::.block_geometry_store(
    dim(coherent),
    function(rows = NULL) if (is.null(rows)) coherent else coherent[rows, , drop = FALSE]
  )
  blocked <- effect_geometry(
    blocked_total, blocked_coherent, geometry$marginals,
    effects = geometry$effects, receipt = geometry$receipt, index = geometry$index
  )

  expect_identical(blocked$completeness, "full")
  expect_true("block_backed" %in% blocked$storage)
  expect_equal(geometry_component(blocked, "configuration"),
    geometry_component(geometry, "configuration"), tolerance = 1e-14)
  expect_equal(geometry_component(blocked, "total", rows = 2),
    total[2, , drop = FALSE], tolerance = 1e-14)
})

test_that("malformed complete and query-only results fail loudly", {
  geometry <- result_fixture()
  expect_error(
    effect_geometry(matrix(1, 1, 2), matrix(1, 2, 1), geometry$marginals,
      effects = geometry$effects, receipt = geometry$receipt),
    "identical dimensions"
  , class = "effect_input_error")
  expect_error(query_geometry(geometry, matrix(1, nrow = 2)), "input dimension",
    class = "effect_input_error")
  expect_error(effect_view(matrix(NA_real_, 1, 1), diag(1), "total",
    receipt = geometry$receipt), "finite", class = "effect_input_error")
})

test_that("full certification rejects malformed packed and marginal structure", {
  geometry <- result_fixture()
  expect_error(effect_geometry(
    matrix(1, 2, 2), matrix(1, 2, 2), geometry$marginals,
    effects = "effect", receipt = geometry$receipt
  ), "not triangular", class = "effect_input_error")

  empty <- structure(list(), semantics = "undirected_endpoint",
    class = c("effect_marginals", "list"))
  expect_error(effect_geometry(
    matrix(1, 2, 1), matrix(1, 2, 1), empty,
    effects = "effect", receipt = geometry$receipt
  ), "nonempty", class = "effect_input_error")

  wrong <- geometry$marginals
  colnames(wrong$endpoint) <- c("other", "b")
  expect_error(effect_geometry(
    matrix(1, 2, 1), matrix(1, 2, 1), wrong,
    effects = "effect", receipt = geometry$receipt
  ), "named effect columns", class = "effect_input_error")
})

test_that("full certification probes manifests, readers, indices, and receipts", {
  geometry <- result_fixture()
  forged <- geometry$total
  forged$manifest$complete <- FALSE
  expect_error(effect_geometry(
    forged, geometry$coherent, geometry$marginals,
    effects = geometry$effects, receipt = geometry$receipt,
    index = geometry$index
  ), "store manifest", class = "effect_input_error")

  broken_reader <- geometry$total
  broken_reader$read <- function(rows = NULL) matrix(1, 1, 1)
  expect_error(effect_geometry(
    broken_reader, geometry$coherent, geometry$marginals,
    effects = geometry$effects, receipt = geometry$receipt,
    index = geometry$index
  ), "reader cannot supply", class = "effect_input_error")

  expect_error(effect_geometry(
    geometry$total, geometry$coherent, geometry$marginals,
    effects = geometry$effects, receipt = geometry$receipt,
    index = c("same", "same")
  ), "uniquely identify", class = "effect_input_error")

  expect_error(effect_geometry(
    geometry$total, geometry$coherent, geometry$marginals,
    effects = geometry$effects, receipt = geometry$receipt,
    metadata = list(scientific_plan_id = "different")
  ), "different scientific plans", class = "effect_contract_error")
})
