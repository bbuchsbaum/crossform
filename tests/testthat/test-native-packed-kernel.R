# Native packed/coherent materialization: oracle parity, codecs, signed
# weights, and route-stable receipts.

packed_kernel_fixture <- function() {
  set.seed(20260816L)
  effects <- c("a", "b", "c")
  left <- list(
    run1 = matrix(rnorm(3 * 7), 3L, 7L, dimnames = list(effects, NULL)),
    run2 = matrix(rnorm(3 * 7), 3L, 7L, dimnames = list(effects, NULL))
  )
  right <- list(
    s1 = matrix(rnorm(2 * 7), 2L, 7L),
    s2 = matrix(rnorm(2 * 7), 2L, 7L)
  )
  list(
    effects = effects,
    left = left,
    right = right,
    self_index = c(1L, 2L, 1L),
    self_right_index = c(2L, 1L, 1L),
    self_weight = c(0.5, -0.25, 1.1),
    cross_left_index = c(1L, 2L),
    cross_right_index = c(2L, 1L),
    cross_weight = c(0.8, -0.3)
  )
}

test_that("packed and rectangular native atoms match the R oracle", {
  fixture <- packed_kernel_fixture()
  for (codec in c("symmetric_packed", "rectangular")) {
    if (identical(codec, "symmetric_packed")) {
      native <- crossform:::.packed_effect_form_atoms(
        fixture$left, fixture$left, fixture$self_index,
        fixture$self_right_index, fixture$self_weight, codec
      )
      oracle <- crossform:::.packed_effect_form_atoms_oracle(
        fixture$left, fixture$left, fixture$self_index,
        fixture$self_right_index, fixture$self_weight, codec
      )
    } else {
      native <- crossform:::.packed_effect_form_atoms(
        fixture$left, fixture$right, fixture$cross_left_index,
        fixture$cross_right_index, fixture$cross_weight, codec
      )
      oracle <- crossform:::.packed_effect_form_atoms_oracle(
        fixture$left, fixture$right, fixture$cross_left_index,
        fixture$cross_right_index, fixture$cross_weight, codec
      )
    }
    expect_equal(unname(native), unname(oracle), tolerance = 1e-12, info = codec)
  }
})

test_that("feature-task packed atoms stay on the native kernel", {
  fixture <- packed_kernel_fixture()
  edges <- crossform:::.ordered_partition_edges(
    pairing("run1", "run2", 1, directed = FALSE, independence = "independent"),
    names(fixture$left), names(fixture$left), TRUE
  )
  got <- crossform:::.effect_form_feature_task(
    fixture$left, fixture$left, seq_len(7L),
    fixture$effects, fixture$effects,
    ordered_edges = edges,
    codec = "symmetric_packed", same_relation = TRUE
  )
  oracle <- crossform:::.packed_effect_form_atoms_oracle(
    fixture$left, fixture$left,
    match(edges$left, names(fixture$left)),
    match(edges$right, names(fixture$left)),
    edges$weight, "symmetric_packed"
  )
  expect_equal(unname(got$atoms), unname(oracle), tolerance = 1e-12)
  expect_identical(got$diagnostics$max_atom_work_bytes, 0)
})

test_that("coherent native tiles match the retained first-moment R loop", {
  set.seed(17)
  measurements <- 11L
  q <- 4L
  partitions <- 3L
  left_first <- array(
    rnorm(measurements * q * partitions),
    dim = c(measurements, q, partitions),
    dimnames = list(NULL, paste0("e", seq_len(q)), paste0("p", seq_len(partitions)))
  )
  edges <- data.frame(
    left = c("p1", "p2", "p3"),
    right = c("p2", "p3", "p1"),
    weight = c(0.4, -0.2, 1.3),
    stringsAsFactors = FALSE
  )
  attr(edges, "expansion") <- "self_adjoint_half_edges"
  mass <- runif(measurements, 0.5, 2)
  native <- crossform:::.coherent_effect_form_atoms_cpp(
    left_first, left_first,
    match(edges$left, dimnames(left_first)[[3L]]),
    match(edges$right, dimnames(left_first)[[3L]]),
    edges$weight, mass, 3L, 5L, TRUE
  )
  oracle <- matrix(0, 5L, q * (q + 1L) / 2L)
  rows <- 3:7
  coordinate <- 0L
  for (column in seq_len(q)) {
    for (row in column:q) {
      coordinate <- coordinate + 1L
      work <- numeric(5L)
      for (edge in seq_len(nrow(edges))) {
        work <- work + edges$weight[[edge]] *
          left_first[rows, row, match(edges$left[[edge]], dimnames(left_first)[[3L]])] *
          left_first[rows, column, match(edges$right[[edge]], dimnames(left_first)[[3L]])]
      }
      work <- work / mass[rows]
      if (row != column) work <- sqrt(2) * work
      oracle[, coordinate] <- work
    }
  }
  expect_equal(unname(native), unname(oracle), tolerance = 1e-12)
})
