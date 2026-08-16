receipt_fixture <- function() {
  execution_receipt(
    scientific_plan_id = "plan-sha256:abc",
    compute = compute_policy(workers = 1, block_features = 64),
    sources = list(source_capabilities(
      block_read = TRUE,
      reopenable = FALSE,
      stable_revision = paste0("sha256:", paste(rep("a", 64), collapse = ""))
    )),
    memory = memory_plan(output_bytes = 128, contraction_bytes = 64),
    kernel_version = "crossgram-v1",
    task_partition_id = "features-64",
    reduction_plan_id = "ascending-block-id",
    precision = "double"
  )
}

result_fixture <- function() {
  local <- array(
    c(2, 4, 3, 5, 6, 8, 7, 9),
    dim = c(2, 2, 2),
    dimnames = list(c("m1", "m2"), c("a", "b"), c("r1", "r2"))
  )
  marginals <- crossform:::pairing_marginals(
    local,
    pairing("r1", "r2", directed = FALSE),
    mass = c(1, 2)
  )
  effect_geometry(
    total = matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE),
    coherent = matrix(c(0.25, 1, 2, 1, 2, 3), nrow = 2, byrow = TRUE),
    marginals = marginals,
    effects = c("a", "b"),
    receipt = receipt_fixture(),
    index = c("m1", "m2")
  )
}
