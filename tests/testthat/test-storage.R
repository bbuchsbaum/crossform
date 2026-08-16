test_that("block-backed complete geometry matches in-memory geometry and queries", {
  set.seed(81)
  weights <- matrix(runif(9 * 12), 9, 12)
  total_atoms <- matrix(rnorm(12 * 10), 12, 10)
  coherent_atoms <- matrix(rnorm(12 * 10), 12, 10)
  total_path <- tempfile(fileext = ".egm")
  coherent_path <- tempfile(fileext = ".egm")
  on.exit(unlink(c(total_path, coherent_path)), add = TRUE)
  total_store <- crossform:::.file_geometry_store(total_path, c(9, 10), create = TRUE)
  coherent_store <- crossform:::.file_geometry_store(coherent_path, c(9, 10), create = TRUE)

  total_run <- crossform:::.tiled_contraction(
    weights, total_atoms, 3, 4, 5,
    write_tile = function(rows, coordinates, value) {
      crossform:::.write_geometry_tile(total_store, rows, coordinates, value)
    }
  )
  coherent_run <- crossform:::.tiled_contraction(
    weights, coherent_atoms, 3, 4, 5,
    write_tile = function(rows, coordinates, value) {
      crossform:::.write_geometry_tile(coherent_store, rows, coordinates, value)
    }
  )
  base <- result_fixture()
  effects <- paste0("e", 1:4)
  marginal_matrix <- matrix(0, 9, 4, dimnames = list(NULL, effects))
  marginals <- structure(list(endpoint = marginal_matrix),
    semantics = "undirected_endpoint", class = c("effect_marginals", "list"))
  blocked <- effect_geometry(total_store, coherent_store, marginals,
    effects = effects, receipt = base$receipt)
  memory <- effect_geometry(weights %*% total_atoms, weights %*% coherent_atoms,
    marginals, effects = effects, receipt = base$receipt)
  query <- matrix(rnorm(10 * 3), 10, 3)

  expect_null(total_run$value)
  expect_null(coherent_run$value)
  expect_equal(file.info(total_path)$size, 9 * 10 * 8)
  expect_equal(geometry_component(blocked, "configuration", rows = 2:5),
    geometry_component(memory, "configuration", rows = 2:5), tolerance = 1e-12)
  expect_equal(query_geometry(blocked, query, row_block = 2)$values,
    query_geometry(memory, query, row_block = 9)$values, tolerance = 1e-11)
})

test_that("block store refuses overwrite and malformed tiles", {
  path <- tempfile(fileext = ".egm")
  on.exit(unlink(path), add = TRUE)
  store <- crossform:::.file_geometry_store(path, c(3, 4), create = TRUE)

  expect_error(crossform:::.file_geometry_store(path, c(3, 4), create = TRUE),
    "Refusing to overwrite")
  expect_error(crossform:::.write_geometry_tile(
    store, c(1, 3), 1, matrix(1, 2, 1)
  ), "Invalid geometry tile")
})

test_that("selected block-store tiles round-trip without full-row reads", {
  path <- tempfile(fileext = ".egm")
  on.exit(unlink(path), add = TRUE)
  store <- crossform:::.file_geometry_store(path, c(5, 6), create = TRUE)
  value <- matrix(seq_len(6), 3, 2)
  crossform:::.write_geometry_tile(store, 2:4, c(2, 5), value)

  expect_equal(crossform:::.read_geometry_tile(store, 2:4, c(2, 5)), value)
  expect_error(crossform:::.read_geometry_tile(store, c(1, 3), 2),
    "Invalid geometry tile")
})
