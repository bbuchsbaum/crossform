test_that("every multivariate route agrees with an independent dense oracle", {
  fixture <- measurement_kernel_fixture()
  expected <- measurement_block_oracle(
    fixture$left_values, fixture$right_values, fixture$h,
    fixture$partition_edges, fixture$measurement_edges,
    fixture$left_operators, fixture$right_operators
  )
  production_reference <- crossform:::.measurement_form_dense_reference(
    fixture$task, fixture$left_values, fixture$right_values
  )
  measurement_expect_blocks_equal(production_reference, expected, 2e-12)

  for (route in c("forward_k", "pull_h", "multivariate_blocks")) {
    run <- crossform:::.run_measurement_contraction(
      fixture$task,
      compute = compute_policy(block_features = 2L),
      route = route,
      edge_tile = 2L
    )
    measurement_expect_blocks_equal(run$blocks, expected, 2e-11)
    expect_false(run$diagnostics$global_neural_operator_materialized)
    expect_false(run$diagnostics$kronecker_lift_materialized)
    expect_identical(run$route_plan$route, route)
  }
})

test_that("low-rank factorized contraction agrees without a neural operator", {
  fixture <- measurement_kernel_fixture(low_rank_h = TRUE)
  expected <- measurement_block_oracle(
    fixture$left_values, fixture$right_values, fixture$h,
    fixture$partition_edges, fixture$measurement_edges,
    fixture$left_operators, fixture$right_operators
  )
  plan <- crossform:::.measurement_route_plan(
    fixture$task, "factorized_h", feature_block = 2L, edge_tile = 3L
  )
  expect_identical(plan$factor_rank, 2L)
  expect_true(plan$eligible[["factorized_h"]])

  run <- crossform:::.run_measurement_contraction(
    fixture$task,
    compute = compute_policy(block_features = 2L),
    route = "factorized_h",
    edge_tile = 3L
  )
  measurement_expect_blocks_equal(run$blocks, expected, 3e-11)
})

test_that("scalar stacking computes only explicitly requested edges", {
  fixture <- measurement_kernel_fixture(scalar = TRUE)
  expected <- measurement_block_oracle(
    fixture$left_values, fixture$right_values, fixture$h,
    fixture$partition_edges, fixture$measurement_edges,
    fixture$left_operators, fixture$right_operators
  )
  plan <- crossform:::.measurement_route_plan(
    fixture$task, "scalar_stack", feature_block = 3L, edge_tile = 2L
  )
  expect_true(plan$eligible[["scalar_stack"]])
  run <- crossform:::.run_measurement_contraction(
    fixture$task,
    compute = compute_policy(block_features = 3L),
    route = "scalar_stack"
  )

  measurement_expect_blocks_equal(run$blocks, expected, 2e-11)
  expect_identical(length(run$blocks), nrow(fixture$measurement_edges))
  expect_identical(names(run$blocks),
    paste0("edge_", sprintf("%06d", seq_len(nrow(fixture$measurement_edges)))))
})

test_that("feature blocking, edge tiling, and eligible routes preserve values", {
  fixture <- measurement_kernel_fixture(low_rank_h = TRUE)
  runs <- list(
    tiny = crossform:::.run_measurement_contraction(
      fixture$task, compute_policy(block_features = 1L),
      route = "pull_h", edge_tile = 1L
    ),
    wide = crossform:::.run_measurement_contraction(
      fixture$task, compute_policy(block_features = 99L),
      route = "multivariate_blocks", edge_tile = 99L
    ),
    factored = crossform:::.run_measurement_contraction(
      fixture$task, compute_policy(block_features = 3L),
      route = "factorized_h", edge_tile = 2L
    )
  )
  measurement_expect_blocks_equal(runs$tiny$blocks, runs$wide$blocks, 3e-11)
  measurement_expect_blocks_equal(runs$factored$blocks, runs$wide$blocks, 3e-11)
  expect_identical(runs$tiny$task_id, runs$wide$task_id)
  expect_false(identical(runs$tiny$route_plan$signature,
    runs$wide$route_plan$signature))
})

test_that("block-backed relation sources obey revisions and bounded reads", {
  base <- measurement_kernel_fixture(seed = 2026081206)
  paths <- tempfile(pattern = "measurement-source-", fileext = ".bin",
    tmpdir = tempdir())
  paths <- c(paths, tempfile(pattern = "measurement-source-", fileext = ".bin",
    tmpdir = tempdir()), tempfile(pattern = "measurement-source-", fileext = ".bin",
    tmpdir = tempdir()), tempfile(pattern = "measurement-source-", fileext = ".bin",
    tmpdir = tempdir()))
  on.exit(unlink(paths), add = TRUE)
  values <- c(base$left_values, base$right_values)
  descriptors <- Map(function(value, path) {
    connection <- file(path, "wb")
    writeBin(as.double(value), connection, size = 8,
      endian = .Platform$endian)
    close(connection)
    file_matrix_source(path, dim(value))
  }, values, paths)
  left_sources <- stats::setNames(descriptors[1:2], names(base$left_values))
  right_sources <- stats::setNames(descriptors[3:4], names(base$right_values))
  fixture <- measurement_kernel_fixture(
    seed = 2026081206,
    left_sources = left_sources,
    right_sources = right_sources
  )
  run <- crossform:::.run_measurement_contraction(
    fixture$task, compute_policy(block_features = 2L), route = "pull_h"
  )
  expected <- measurement_block_oracle(
    fixture$left_values, fixture$right_values, fixture$h,
    fixture$partition_edges, fixture$measurement_edges,
    fixture$left_operators, fixture$right_operators
  )

  measurement_expect_blocks_equal(run$blocks, expected, 2e-11)
  expect_gt(run$diagnostics$source_bytes_read, 0)
  expect_lte(run$diagnostics$max_relation_block_bytes,
    as.double(utils::object.size(matrix(0, 4, 2))))

  connection <- file(paths[[1L]], "r+b")
  writeBin(999, connection, size = 8, endian = .Platform$endian)
  close(connection)
  failure <- tryCatch(
    crossform:::.run_measurement_contraction(
      fixture$task, compute_policy(block_features = 2L), route = "pull_h"
    ),
    error = identity
  )
  expect_s3_class(failure, "effect_measurement_kernel_error")
  expect_identical(failure$receipt$completion_status, "failed")
  expect_match(failure$receipt$message, "stale|revision", ignore.case = TRUE)
})

test_that("workspace rejection occurs before lazy source reads", {
  reads <- new.env(parent = emptyenv())
  reads$count <- 0L
  source <- function(features) {
    reads$count <- reads$count + 1L
    matrix(1, 4, length(features))
  }
  revision <- paste0("sha256:", paste(rep("b", 64), collapse = ""))
  capabilities <- source_capabilities(TRUE, stable_revision = revision)
  left_sources <- stats::setNames(rep(list(source), 2), c("encode_1", "encode_2"))
  right_source <- function(features) {
    reads$count <- reads$count + 1L
    matrix(1, 3, length(features))
  }
  right_sources <- stats::setNames(rep(list(right_source), 2),
    c("retrieve_1", "retrieve_2"))
  fixture <- measurement_kernel_fixture(seed = 2026081207)
  left_relation <- relation(
    left_sources,
    source_dims = rep(list(c(4L, 5L)), 2),
    effects = fixture$task$spaces$experimental_left,
    domain = fixture$task$spaces$neural_left,
    capabilities = capabilities
  )
  right_relation <- relation(
    right_sources,
    source_dims = rep(list(c(3L, 6L)), 2),
    effects = fixture$task$spaces$experimental_right,
    domain = fixture$task$spaces$neural_right,
    capabilities = capabilities
  )
  task <- crossform:::.new_evidence_task(
    left_relation, right_relation, FALSE,
    fixture$task$ordered_partition_products,
    fixture$task$experimental_boundary,
    fixture$task$neural_boundary,
    fixture$task$stages,
    fixture$task$materialization
  )
  failure <- tryCatch(
    crossform:::.run_measurement_contraction(
      task,
      compute_policy(block_features = 2L, workspace_bytes = 1),
      route = "pull_h"
    ),
    error = identity
  )

  expect_s3_class(failure, "effect_measurement_kernel_error")
  expect_match(failure$receipt$message, "workspace budget")
  expect_identical(reads$count, 0L)
  expect_false(failure$receipt$memory$fits_budget)
})

test_that("extreme reciprocal scales and near-degenerate H remain finite", {
  fixture <- measurement_kernel_fixture(low_rank_h = TRUE, seed = 2026081208)
  scaled_left <- lapply(fixture$left_values, `*`, 1e70)
  scaled_right <- lapply(fixture$right_values, `*`, 1e-70)
  left_relation <- relation(
    scaled_left,
    effects = fixture$task$spaces$experimental_left,
    domain = fixture$task$spaces$neural_left
  )
  right_relation <- relation(
    scaled_right,
    effects = fixture$task$spaces$experimental_right,
    domain = fixture$task$spaces$neural_right
  )
  h <- fixture$h
  h[, ncol(h)] <- h[, 1] + 1e-13 * h[, ncol(h)]
  query <- pair_query(
    h,
    fixture$task$spaces$experimental_left,
    fixture$task$spaces$experimental_right
  )
  task <- crossform:::.new_evidence_task(
    left_relation, right_relation, FALSE,
    fixture$task$ordered_partition_products,
    crossform:::.closed_experimental_boundary(query),
    fixture$task$neural_boundary,
    fixture$task$stages,
    fixture$task$materialization
  )
  expected <- measurement_block_oracle(
    scaled_left, scaled_right, h,
    fixture$partition_edges, fixture$measurement_edges,
    fixture$left_operators, fixture$right_operators
  )
  run <- crossform:::.run_measurement_contraction(
    task, compute_policy(block_features = 2L), route = "pull_h"
  )

  expect_true(all(vapply(run$blocks, function(x) all(is.finite(x)), logical(1))))
  measurement_expect_blocks_equal(run$blocks, expected, 5e-10)
})

test_that("route and memory plans account for every named live buffer", {
  fixture <- measurement_kernel_fixture()
  route <- crossform:::.measurement_route_plan(
    fixture$task, "pull_h", feature_block = 2L, edge_tile = 2L
  )
  memory <- crossform:::.measurement_kernel_memory_plan(
    fixture$task, route, "memory"
  )

  expect_identical(names(memory$buffers), c(
    "relation", "factor_slice", "projected_state", "query_or_factor",
    "edge", "output", "replacement"
  ))
  expect_true(all(memory$buffers > 0))
  expect_gte(memory$plan$planned_workspace_bytes,
    sum(memory$plan$categories))
  expect_true(all(is.finite(route$costs[route$eligible])))
})

test_that("measurement preflight uses arithmetic rather than buffer allocation", {
  fixture <- measurement_kernel_fixture()
  route <- crossform:::.measurement_route_plan(
    fixture$task, "pull_h", feature_block = 2L, edge_tile = 2L
  )

  memory <- crossform:::.measurement_kernel_memory_plan(
    fixture$task, route, "memory"
  )

  expect_s3_class(memory, "effect_measurement_memory_plan")
  expect_true(all(is.finite(memory$buffers)))
  expect_false(any(grepl(
    "object\\.size\\(matrix\\(|object\\.size\\(array\\(",
    paste(deparse(body(crossform:::.measurement_kernel_memory_plan)),
      collapse = "\n")
  )))
})

test_that("small benchmark receipts distinguish scalar and multivariate costs", {
  scalar <- measurement_kernel_fixture(scalar = TRUE, seed = 2026081209)
  multivariate <- measurement_kernel_fixture(scalar = FALSE, seed = 2026081210)
  scalar_record <- crossform:::.measurement_benchmark_record(
    scalar$task, iterations = 1L, feature_block = 5L, edge_tile = 4L
  )
  multivariate_record <- crossform:::.measurement_benchmark_record(
    multivariate$task, iterations = 1L, feature_block = 5L, edge_tile = 4L
  )

  expect_true(all(scalar_record$scalar_frame))
  expect_false(any(multivariate_record$scalar_frame))
  expect_true(all(scalar_record$median_elapsed_seconds >= 0))
  expect_true(all(multivariate_record$output_bytes > 0))
})
