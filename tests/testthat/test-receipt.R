test_that("receipt separates scientific and numerical execution identity", {
  receipt <- receipt_fixture()

  expect_s3_class(receipt, "effect_execution_receipt")
  expect_identical(receipt$scientific_plan_id, "plan-sha256:abc")
  expect_null(receipt$domain_signature)
  expect_identical(receipt$compute$process_backend, "sequential")
  expect_identical(receipt$sources[[1]]$stable_revision,
    paste0("sha256:", paste(rep("a", 64), collapse = "")))
  expect_identical(receipt$reduction_plan_id, "ascending-block-id")
  expect_identical(receipt$completion_status, "complete")
  expect_identical(receipt$task_count, receipt$completed_task_count)
  expect_identical(receipt$blas$requested_threads, 1L)
  expect_true(is.na(receipt$blas$observed_threads))
  expect_null(receipt$reporter)
})

test_that("compiler receipts expose exact neural-domain identity", {
  domain <- abstract_domain(3, feature_ids = c("left", "middle", "right"),
    id = "line:v1")
  matrices <- list(run1 = matrix(1, 2, 3), run2 = matrix(2, 2, 3))
  rel <- relation(matrices, effects = c("a", "b"), domain = domain)
  geometry <- geometry(rel, compile_frame(voxels(), domain), cross_partitions(rel))

  expect_identical(geometry$receipt$domain_signature,
    domain$reference$signature)
  expect_identical(geometry$metadata$frame$domain,
    domain$reference)
})

test_that("receipts preserve every compute field and completion fact", {
  compute <- compute_policy(block_features = 32, workspace_bytes = 4096)
  receipt <- execution_receipt(
    "plan", compute,
    list(source_capabilities(TRUE,
      stable_revision = paste0("sha256:", paste(rep("c", 64), collapse = "")))),
    memory_plan(workers = 1, budget_bytes = 4096),
    "kernel", "tasks-32", "ascending", precision = "double",
    completion_status = "complete", task_count = 9,
    completed_task_count = 9, elapsed_seconds = 1.25,
    blas = list(vendor = "Accelerate", requested_threads = 1,
      observed_threads = NA_integer_)
  )

  expect_identical(receipt$compute, compute)
  expect_identical(receipt$task_count, 9)
  expect_identical(receipt$completed_task_count, 9)
  expect_identical(receipt$elapsed_seconds, 1.25)
  expect_identical(receipt$blas, list(vendor = "Accelerate",
    requested_threads = 1L, observed_threads = NA_integer_))
})

test_that("receipt construction rejects contradictory execution claims", {
  source <- list(source_capabilities(TRUE,
    stable_revision = paste0("sha256:", paste(rep("d", 64), collapse = ""))))
  expect_error(execution_receipt(
    "plan", compute_policy(), source, memory_plan(workers = 2, n_active = 1),
    "kernel", "tasks", "reduction"
  ), "same worker count")
  expect_error(execution_receipt(
    "plan", compute_policy(workspace_bytes = 100), source,
    memory_plan(budget_bytes = 200), "kernel", "tasks", "reduction"
  ), "same memory budget")
  expect_error(execution_receipt(
    "plan", compute_policy(), source, memory_plan(),
    "kernel", "tasks", "reduction", precision = "single"
  ), "only.*double")
  expect_error(execution_receipt(
    "plan", compute_policy(), source, memory_plan(),
    "kernel", "tasks", "reduction", task_count = 3,
    completed_task_count = 2
  ), "every task")
  expect_error(execution_receipt(
    "plan", compute_policy(), source, memory_plan(),
    "kernel", "tasks", "reduction", completion_status = "planned",
    task_count = 3, completed_task_count = 1
  ), "cannot report")
  expect_error(execution_receipt(
    "plan", compute_policy(), source, memory_plan(),
    "kernel", "tasks", "reduction",
    blas = list(vendor = "OpenBLAS", requested_threads = 2,
      observed_threads = 2)
  ), "requested and observed")
  blocked_source <- list(source_capabilities(FALSE,
    stable_revision = paste0("sha256:", paste(rep("e", 64), collapse = ""))))
  expect_error(execution_receipt(
    "plan", compute_policy(), blocked_source, memory_plan(),
    "kernel", "tasks", "reduction"
  ), "block reads")
})

test_that("mutated receipt internals fail result certification", {
  geometry <- result_fixture()

  bad_compute <- geometry$receipt
  bad_compute$compute$workers <- 2
  expect_error(effect_geometry(
    geometry_component(geometry, "total"),
    geometry_component(geometry, "coherent"), geometry$marginals,
    effects = geometry$effects, receipt = bad_compute, index = geometry$index
  ), "workers")

  bad_memory <- geometry$receipt
  bad_memory$memory$modeled_workspace_bytes <- 999
  expect_error(effect_geometry(
    geometry_component(geometry, "total"),
    geometry_component(geometry, "coherent"), geometry$marginals,
    effects = geometry$effects, receipt = bad_memory, index = geometry$index
  ), "derived fields")

  bad_source <- geometry$receipt
  bad_source$sources[[1]]$stable_revision <- "path-and-mtime"
  expect_error(effect_geometry(
    geometry_component(geometry, "total"),
    geometry_component(geometry, "coherent"), geometry$marginals,
    effects = geometry$effects, receipt = bad_source, index = geometry$index
  ), "sha256")

  bad_status <- geometry$receipt
  bad_status$completed_task_count <- 0
  expect_error(effect_geometry(
    geometry_component(geometry, "total"),
    geometry_component(geometry, "coherent"), geometry$marginals,
    effects = geometry$effects, receipt = bad_status, index = geometry$index
  ), "completion counts")
})

test_that("source capabilities require a strong stable revision", {
  expect_error(
    source_capabilities(TRUE, stable_revision = ""),
    "sha256 identifier"
  )
  expect_error(
    source_capabilities(NA,
      stable_revision = paste0("sha256:", paste(rep("b", 64), collapse = ""))),
    "TRUE or FALSE"
  )
})

test_that("receipt rejects noncanonical execution objects", {
  expect_error(
    execution_receipt(
      "plan", compute_policy(), list(), memory_plan(),
      "kernel", "tasks", "reduction"
    ),
    "nonempty list"
  )
  expect_error(
    execution_receipt(
      "", compute_policy(),
      list(source_capabilities(TRUE,
        stable_revision = paste0("sha256:", paste(rep("b", 64), collapse = "")))),
      memory_plan(), "kernel", "tasks", "reduction"
    ),
    "nonempty character"
  )
})
