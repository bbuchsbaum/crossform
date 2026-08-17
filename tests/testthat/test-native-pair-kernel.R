# Native fused pair-difference kernel: oracle parity, signed weights,
# selected tiles, RSA contraction, and route-stable receipts.

native_pair_fixture <- function() {
  set.seed(20260816L)
  effects <- c("a", "b", "c", "d")
  left <- list(
    run1 = matrix(rnorm(4 * 5), 4L, 5L, dimnames = list(effects, NULL)),
    run2 = matrix(rnorm(4 * 5), 4L, 5L, dimnames = list(effects, NULL)),
    run3 = matrix(rnorm(4 * 5), 4L, 5L, dimnames = list(effects, NULL))
  )
  query <- crossform:::.pair_difference_query(effects)
  list(
    effects = effects,
    left = left,
    query = query,
    left_index = c(1L, 2L, 3L, 1L),
    right_index = c(2L, 3L, 1L, 3L),
    edge_weight = c(0.4, -0.25, 0.7, 0.15)
  )
}

native_geometry_fixture <- function(effects = 6L, features = 9L) {
  domain <- abstract_domain(
    features,
    coordinates = cbind(x = seq_len(features), y = 0),
    id = "native-pair-domain"
  )
  effect_names <- paste0("effect", seq_len(effects))
  relation <- relation(
    list(
      run1 = matrix(seq_len(effects * features) / 10, effects, features,
        dimnames = list(effect_names, NULL)),
      run2 = matrix(rev(seq_len(effects * features)) / 11, effects, features,
        dimnames = list(effect_names, NULL))
    ),
    domain = domain
  )
  list(
    relation = relation,
    frame = compile_frame(searchlights(radius = 1.01), domain),
    pairing = cross_partitions(relation),
    effect_names = effect_names
  )
}

test_that("the native kernel matches the retained two-pass R oracle", {
  fixture <- native_pair_fixture()
  native <- crossform:::.fused_pair_difference_atoms(
    fixture$left, fixture$left, fixture$left_index, fixture$right_index,
    fixture$edge_weight, fixture$query$pair_left, fixture$query$pair_right
  )
  oracle <- crossform:::.pair_difference_accumulate_oracle(
    fixture$left, fixture$left, fixture$left_index, fixture$right_index,
    fixture$edge_weight, fixture$query$pair_left, fixture$query$pair_right
  )

  expect_equal(unname(native), unname(oracle), tolerance = 1e-12)
  expect_identical(dim(native), c(5L, 6L))
})

test_that("selected pair tiles and RSA coefficients stay within 1e-12", {
  fixture <- native_pair_fixture()
  pair_indices <- c(6L, 1L, 4L)
  coefficients <- matrix(
    c(1, -1, 0.5, 0, 0, 2, 0.25, 0, -0.5, 1, 0, 0.75),
    2L, 6L, byrow = TRUE
  )

  native_tile <- crossform:::.fused_pair_difference_atoms(
    fixture$left, fixture$left, fixture$left_index, fixture$right_index,
    fixture$edge_weight, fixture$query$pair_left, fixture$query$pair_right,
    pair_indices = pair_indices
  )
  oracle_tile <- crossform:::.pair_difference_accumulate_oracle(
    fixture$left, fixture$left, fixture$left_index, fixture$right_index,
    fixture$edge_weight, fixture$query$pair_left, fixture$query$pair_right,
    pair_indices = pair_indices
  )
  expect_equal(unname(native_tile), unname(oracle_tile), tolerance = 1e-12)

  native_rsa <- crossform:::.fused_pair_difference_atoms(
    fixture$left, fixture$left, fixture$left_index, fixture$right_index,
    fixture$edge_weight, fixture$query$pair_left, fixture$query$pair_right,
    coefficients
  )
  oracle_rsa <- crossform:::.pair_difference_accumulate_oracle(
    fixture$left, fixture$left, fixture$left_index, fixture$right_index,
    fixture$edge_weight, fixture$query$pair_left, fixture$query$pair_right,
    coefficients
  )
  expect_equal(unname(native_rsa), unname(oracle_rsa), tolerance = 1e-12)
  expect_identical(dim(native_rsa), c(5L, 2L))
})

test_that("fused RDM and RSA agree with memory and block materialization", {
  fixture <- native_geometry_fixture()
  plan <- plan_geometry(
    fixture$relation, fixture$frame, fixture$pairing,
    compute = compute_policy(block_features = 4)
  )
  model <- as.matrix(stats::dist(seq_len(6L)))
  dimnames(model) <- list(fixture$effect_names, fixture$effect_names)
  selected <- cbind(
    fixture$effect_names[c(1L, 2L, 4L)],
    fixture$effect_names[c(3L, 6L, 5L)]
  )

  complete_memory <- materialize_geometry(plan)
  complete_block <- materialize_geometry(
    plan, storage = "block",
    storage_path = file.path(tempdir(), paste0("native-pair-", Sys.getpid()))
  )
  fused_rdm <- rdm(plan)
  selected_rdm <- rdm(plan, pairs = selected)
  fused_rsa <- rsa(plan, models = list(distance = model))
  late_rdm <- rdm(complete_memory)
  late_rsa <- rsa(complete_memory, models = list(distance = model))

  expect_equal(fused_rdm$values, late_rdm$values, tolerance = 1e-12)
  expect_equal(fused_rdm$values, rdm(complete_block)$values, tolerance = 1e-12)
  expect_equal(
    selected_rdm$values, rdm(complete_memory, pairs = selected)$values,
    tolerance = 1e-12
  )
  expect_equal(fused_rsa$coefficients, late_rsa$coefficients, tolerance = 1e-12)
  expect_identical(
    fused_rdm$receipt$scientific_plan_id,
    late_rdm$receipt$scientific_plan_id
  )
  expect_identical(
    fused_rdm$receipt$kernel_version, "additive-query-fused-native-v1"
  )
  expect_identical(
    late_rdm$receipt$kernel_version, "additive-packed-native-v1"
  )
  direct <- evaluate_geometry(
    plan,
    query = crossform:::.pair_difference_query(fixture$effect_names)
  )
  expect_identical(
    direct$metadata$execution_plan$lowering,
    "additive_query_fused_contraction"
  )
  expect_identical(
    direct$receipt$kernel_version, "additive-query-fused-native-v1"
  )
})
