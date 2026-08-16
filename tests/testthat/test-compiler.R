compiler_fixture <- function() {
  run1 <- matrix(c(2, 1, 0, -1, 0, 1, 3, 2), nrow = 2, byrow = TRUE,
    dimnames = list(c("face", "house"), NULL))
  run2 <- matrix(c(3, 0, 1, -2, 1, 2, 2, 1), nrow = 2, byrow = TRUE,
    dimnames = list(c("face", "house"), NULL))
  domain <- abstract_domain(4, coordinates = matrix(c(
    0, 0, 1, 0, 2, 0, 3, 0
  ), ncol = 2, byrow = TRUE), id = "line:v1")
  rel <- relation(list(run1 = run1, run2 = run2), domain = domain)
  at <- compile_frame(searchlights(radius = 1.01, normalization = "local"), domain)
  list(relation = rel, at = at, over = cross_partitions(
    rel$partitions, independence = "independent"
  ),
    matrices = list(run1 = run1, run2 = run2), domain = domain)
}

test_that("direct execution preserves compiled query labels", {
  fixture <- compiler_fixture()
  query <- matrix(c(1, 0, 0, 0, 1, 0), 3, 2,
    dimnames = list(NULL, c("first", "second")))
  got <- evaluate_geometry(fixture$relation, fixture$at, fixture$over, query)
  expect_identical(colnames(got$values), c("first", "second"))
})

compiler_oracle <- function(fixture) {
  weights <- as.matrix(fixture$at$weights)
  b1 <- fixture$matrices$run1
  b2 <- fixture$matrices$run2
  q <- nrow(b1)
  total <- coherent <- matrix(0, nrow(weights), q * (q + 1) / 2)
  local <- array(0, c(nrow(weights), q, 2),
    dimnames = list(NULL, rownames(b1), c("run1", "run2")))
  for (j in seq_len(nrow(weights))) {
    total_matrix <- matrix(0, q, q)
    for (v in seq_len(ncol(weights))) {
      cross <- outer(b1[, v], b2[, v])
      total_matrix <- total_matrix + weights[j, v] * 0.5 * (cross + t(cross))
    }
    total[j, ] <- oracle_svec(total_matrix)
    local[j, , 1] <- b1 %*% weights[j, ]
    local[j, , 2] <- b2 %*% weights[j, ]
    local_cross <- outer(local[j, , 1], local[j, , 2])
    coherent[j, ] <- oracle_svec(
      0.5 * (local_cross + t(local_cross)) / sum(weights[j, ])
    )
  }
  list(total = total, coherent = coherent)
}

test_that("public geometry compiles relation frame and pairing once", {
  fixture <- compiler_fixture()
  oracle <- compiler_oracle(fixture)
  got <- materialize_geometry(fixture$relation, fixture$at, fixture$over,
    compute = compute_policy(block_features = 2))

  expect_s3_class(got, "effect_geometry")
  expect_equal(geometry_component(got, "total"), oracle$total,
    tolerance = 1e-12)
  expect_equal(geometry_component(got, "coherent"), oracle$coherent,
    tolerance = 1e-12)
  expect_equal(geometry_component(got, "configuration"),
    oracle$total - oracle$coherent, tolerance = 1e-12)
  expect_identical(got$receipt$completion_status, "complete")
  expect_identical(got$receipt$task_count, 2)
  expect_identical(got$receipt$completed_task_count, 2)
  expect_identical(got$receipt$observed$task_counts,
    c(planned = 2, started = 2, completed = 2, failed = 0, retried = 0))
  expect_identical(got$receipt$observed$features_completed, 4)
  expect_gt(got$receipt$observed$bytes_read, 0)
  expect_identical(got$receipt$observed$tiles$feature_block, 2L)
  expect_true(all(c("source_admission", "feature_tasks", "coherent", "marginals") %in%
    names(got$receipt$observed$stage_seconds)))
  expect_true(got$receipt$observed$cleanup$success)
  expect_false(got$receipt$observed$checkpoint$enabled)
  expect_true(is.na(got$receipt$blas$observed_threads))
  expect_identical(got$metadata$pairing_estimate, "cross_generalized")
})

test_that("direct queries equal late queries without claiming full geometry", {
  fixture <- compiler_fixture()
  full <- materialize_geometry(fixture$relation, fixture$at, fixture$over,
    compute = compute_policy(block_features = 3))
  query <- bilinear_query(tcrossprod(c(1, -1)))

  for (component in c("total", "coherent", "configuration")) {
    direct <- evaluate_geometry(fixture$relation, fixture$at, fixture$over,
      query, component = component, compute = compute_policy(block_features = 2))
    late <- query_geometry(full, query, component = component)
    expect_s3_class(direct, "effect_view")
    expect_identical(direct$completeness, "query_only")
    expect_equal(direct$values, late$values, tolerance = 1e-12)
    expect_identical(direct$receipt$completion_status, "complete")
    expect_identical(direct$metadata$requirements$total,
      component %in% c("total", "configuration"))
    expect_identical(direct$metadata$requirements$coherent,
      component %in% c("coherent", "configuration"))
    expect_false(direct$metadata$requirements$marginals)
    if (component == "total") {
      expect_null(direct$metadata$diagnostics$coherent)
      expect_equal(direct$metadata$diagnostics$total$durable_local_relation_bytes,
        0)
      expect_equal(direct$receipt$memory$categories[["atom_block"]] > 0, TRUE)
    }
    if (component == "coherent") {
      expect_null(direct$metadata$diagnostics$total)
      expect_equal(direct$receipt$memory$categories[["atom_block"]], 0)
      expect_gt(direct$receipt$memory$categories[["output"]], 0)
    }
  }
})

test_that("compiler rejects a same-shaped query from another effect space", {
  fixture <- compiler_fixture()
  query <- bilinear_query(diag(2), effects = effect_space(
    c("face", "house"), basis_id = "different:basis"))

  expect_error(evaluate_geometry(fixture$relation, fixture$at, fixture$over,
    query), "effect spaces are incompatible")
})

test_that("experimental reparameterization preserves bound scalar queries", {
  fixture <- compiler_fixture()
  base_space <- effect_space(c("face", "house"), basis_id = "means:v1")
  transformed_space <- effect_space(c("sum", "difference"),
    basis_id = "sum-difference:v1")
  base <- relation(fixture$matrices, effects = base_space,
    domain = fixture$domain)
  transform <- matrix(c(1, 1, 1, -1), 2, byrow = TRUE)
  transformed <- relation(lapply(fixture$matrices, function(value) {
    unname(transform %*% value)
  }), effects = transformed_space, domain = fixture$domain)
  operator <- tcrossprod(c(1, -0.5))
  inverse <- solve(transform)
  transformed_operator <- t(inverse) %*% operator %*% inverse

  original <- evaluate_geometry(base, fixture$at, fixture$over,
    bilinear_query(operator, effects = base_space))
  changed <- evaluate_geometry(transformed, fixture$at, fixture$over,
    bilinear_query(transformed_operator, effects = transformed_space))

  expect_equal(original$values, changed$values, tolerance = 1e-12)
  expect_false(identical(original$receipt$scientific_plan_id,
    changed$receipt$scientific_plan_id))
})

test_that("block-backed public geometry remains complete and readable", {
  fixture <- compiler_fixture()
  path <- tempfile("crossform-store-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  memory <- materialize_geometry(fixture$relation, fixture$at, fixture$over)
  blocked <- materialize_geometry(fixture$relation, fixture$at, fixture$over,
    storage = "block", storage_path = path,
    compute = compute_policy(block_features = 2))

  expect_identical(blocked$completeness, "full")
  expect_true(all(file.exists(file.path(path, c("total.egm", "coherent.egm")))))
  expect_equal(geometry_component(blocked, "total"),
    geometry_component(memory, "total"), tolerance = 1e-12)
  expect_equal(geometry_component(blocked, "configuration"),
    geometry_component(memory, "configuration"), tolerance = 1e-12)
})

test_that("compiler reuses one admitted file handle across feature blocks", {
  values <- matrix(seq_len(16), 2, 8)
  path <- tempfile("crossform-session-", fileext = ".bin")
  on.exit(unlink(path), add = TRUE)
  connection <- file(path, open = "wb")
  writeBin(as.double(values), connection, size = 8, endian = .Platform$endian)
  close(connection)
  descriptor <- file_matrix_source(path, dim(values))
  domain <- abstract_domain(8)
  rel <- relation(list(run1 = descriptor, run2 = descriptor),
    effects = c("a", "b"), domain = domain)
  got <- materialize_geometry(rel, compile_frame(whole_brain(), domain),
    cross_partitions(rel, independence = "independent"),
    compute = compute_policy(block_features = 2))

  expect_identical(got$metadata$source_session$distinct_owned_handles, 1L)
  expect_identical(unname(got$metadata$source_session$read_count), c(4L, 4L))
  expect_identical(got$metadata$source_session$close_attempts, 1L)
  expect_true(got$metadata$source_session$closed)
})

test_that("compiler hashes a file source once at session admission", {
  values <- matrix(seq_len(12), 2, 6)
  path <- tempfile("crossform-hash-session-", fileext = ".bin")
  on.exit(unlink(path), add = TRUE)
  connection <- file(path, open = "wb")
  writeBin(as.double(values), connection, size = 8, endian = .Platform$endian)
  close(connection)
  descriptor <- file_matrix_source(path, dim(values))
  original_hash <- crossform:::.sha256_file
  hashes <- 0L
  testthat::local_mocked_bindings(
    .sha256_file = function(path) {
      hashes <<- hashes + 1L
      original_hash(path)
    },
    .package = "crossform"
  )
  domain <- abstract_domain(6)
  rel <- relation(list(run1 = descriptor, run2 = descriptor),
    effects = c("a", "b"), domain = domain)

  materialize_geometry(rel, compile_frame(whole_brain(), domain),
    cross_partitions(rel, independence = "independent"),
    compute = compute_policy(block_features = 1))

  expect_identical(hashes, 1L)
})

test_that("compile and budget failures occur before opaque source reads", {
  reads <- 0L
  source <- function(features) {
    reads <<- reads + 1L
    matrix(1, 2, length(features))
  }
  revision <- paste0("sha256:", paste(rep("a", 64), collapse = ""))
  rel <- relation(list(run1 = source, run2 = source),
    source_dims = list(c(2, 4), c(2, 4)), effects = c("a", "b"),
    capabilities = source_capabilities(TRUE, stable_revision = revision))
  at <- compile_frame(voxelwise(), abstract_domain(4))
  over <- cross_partitions(rel$partitions, independence = "independent")

  expect_error(evaluate_geometry(rel, at, over, matrix(1, 2, 1)),
    "packed-coordinate")
  expect_identical(reads, 0L)
  expect_error(materialize_geometry(rel, at, over,
    compute = compute_policy(workspace_bytes = 1)), "exceeding")
  expect_identical(reads, 0L)
})

test_that("a feasible workspace budget drives smaller legal tiles", {
  set.seed(44)
  p <- 64L
  m <- 32L
  domain <- abstract_domain(p)
  matrices <- list(
    run1 = matrix(rnorm(3 * p), 3, p),
    run2 = matrix(rnorm(3 * p), 3, p)
  )
  rownames(matrices$run1) <- rownames(matrices$run2) <- c("a", "b", "c")
  rel <- relation(matrices, domain = domain)
  weights <- matrix(runif(m * p), m, p)
  at <- additive_frame(weights, domain = domain)
  over <- cross_partitions(rel, independence = "independent")
  query <- matrix(1, 6, 1)
  requirements <- crossform:::.component_requirements(query, "total")
  unconstrained <- crossform:::.compiler_memory_plan(
    rel, at, compute_policy(), p, m, 1L, 1L, "memory", requirements
  )
  minimum <- crossform:::.compiler_memory_plan(
    rel, at, compute_policy(), 1L, 1L, 1L, 1L, "memory", requirements
  )
  budget <- floor((unconstrained$planned_workspace_bytes +
    minimum$planned_workspace_bytes) / 2)

  got <- evaluate_geometry(rel, at, over, query,
    compute = compute_policy(workspace_bytes = budget))

  expect_lte(got$receipt$memory$planned_workspace_bytes, budget)
  expect_true(got$metadata$tiles$feature_block < p ||
    got$metadata$tiles$row_tile < m)
  expect_identical(got$receipt$memory$budget_bytes, budget)
})

test_that("memory plans account for actual dense and sparse frame storage", {
  domain <- abstract_domain(20)
  matrices <- list(run1 = matrix(1, 2, 20), run2 = matrix(2, 2, 20))
  rel <- relation(matrices, effects = c("a", "b"), domain = domain)
  dense <- additive_frame(matrix(1, 10, 20), domain = domain)
  sparse <- additive_frame(Matrix::Matrix(diag(20)[1:10, ], sparse = TRUE),
    domain = domain)
  requirements <- crossform:::.component_requirements(matrix(1, 3, 1),
    "total")
  dense_plan <- crossform:::.compiler_memory_plan(rel, dense,
    compute_policy(), 10, 10, 1, 1, "memory", requirements)
  sparse_plan <- crossform:::.compiler_memory_plan(rel, sparse,
    compute_policy(), 10, 10, 1, 1, "memory", requirements)

  expect_equal(dense_plan$categories[["frame"]],
    as.double(utils::object.size(dense$weights)))
  expect_equal(sparse_plan$categories[["frame"]],
    as.double(utils::object.size(sparse$weights)))
  expect_false(identical(dense_plan$categories[["frame"]],
    sparse_plan$categories[["frame"]]))
})

test_that("exact domain mismatches fail before lazy source reads", {
  reads <- 0L
  source <- function(features) {
    reads <<- reads + 1L
    matrix(1, 2, length(features))
  }
  revision <- paste0("sha256:", paste(rep("b", 64), collapse = ""))
  relation_domain <- abstract_domain(4, feature_ids = letters[1:4],
    id = "same-label")
  reversed_domain <- abstract_domain(4, feature_ids = rev(letters[1:4]),
    id = "same-label")
  rel <- relation(list(run1 = source, run2 = source),
    source_dims = list(c(2, 4), c(2, 4)), effects = c("a", "b"),
    domain = relation_domain,
    capabilities = source_capabilities(TRUE, stable_revision = revision))
  at <- additive_frame(diag(4), domain = reversed_domain)

  expect_error(materialize_geometry(rel, at,
    cross_partitions(rel, independence = "independent")),
  "exact neural-domain")
  expect_identical(reads, 0L)
})

test_that("volume spacing participates in compiler domain identity", {
  mask <- array(TRUE, c(2, 2, 1))
  first <- volume_domain(mask, spacing = c(1, 1, 1), id = "volume")
  changed <- volume_domain(mask, spacing = c(2, 1, 1), id = "volume")
  matrices <- list(run1 = matrix(1, 2, 4), run2 = matrix(2, 2, 4))
  rel <- relation(matrices, effects = c("a", "b"), domain = first)

  expect_error(materialize_geometry(rel, compile_frame(voxelwise(), changed),
    cross_partitions(rel, independence = "independent")),
    "exact neural-domain")
})

test_that("opaque sources without revisions fail before reading", {
  reads <- 0L
  source <- function(features) {
    reads <<- reads + 1L
    matrix(1, 2, length(features))
  }
  rel <- relation(list(run1 = source, run2 = source),
    source_dims = list(c(2, 3), c(2, 3)), effects = c("a", "b"))
  at <- compile_frame(voxelwise(), abstract_domain(3))

  expect_error(materialize_geometry(rel, at,
    cross_partitions(rel$partitions, independence = "independent")),
    "explicit `source_capabilities")
  expect_identical(reads, 0L)
})

test_that("kernel failure removes only newly created block stores", {
  calls <- 0L
  failing <- function(features) {
    calls <<- calls + 1L
    if (calls > 2L) stop("source failed")
    matrix(1, 2, length(features))
  }
  revision <- paste0("sha256:", paste(rep("b", 64), collapse = ""))
  rel <- relation(list(run1 = failing, run2 = failing),
    source_dims = list(c(2, 4), c(2, 4)), effects = c("a", "b"),
    capabilities = source_capabilities(TRUE, stable_revision = revision))
  path <- tempfile("failed-geometry-")
  condition <- tryCatch(
    materialize_geometry(rel, compile_frame(voxelwise(), abstract_domain(4)),
      cross_partitions(rel$partitions, independence = "independent"),
      storage = "block", storage_path = path,
      compute = compute_policy(block_features = 2)),
    effect_execution_error = identity
  )

  expect_s3_class(condition, "effect_execution_error")
  expect_identical(condition$receipt$completion_status, "failed")
  expect_true(condition$cleanup_status$success)
  expect_false(dir.exists(path))
})

test_that("reporters remain nonsemantic at the public compiler boundary", {
  fixture <- compiler_fixture()
  got <- materialize_geometry(fixture$relation, fixture$at, fixture$over,
    reporter = function(event) stop("display failed"))
  expect_s3_class(got, "effect_geometry")
  expect_identical(got$receipt$completion_status, "complete")
  expect_false("reporter" %in% names(got$receipt))
  expect_true(length(got$receipt$observed$reporter_failures) >= 1L)
})
