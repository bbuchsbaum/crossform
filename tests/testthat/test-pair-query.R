pair_query_spaces <- function() {
  list(
    left = effect_space(c("e1", "e2", "e3"), basis_id = "encoding:v1"),
    right = effect_space(c("r1", "r2", "r3", "r4"),
      basis_id = "retrieval:v1")
  )
}

full_pair_design <- function(spaces) {
  design <- expand.grid(
    left = spaces$left$coordinates,
    right = spaces$right$coordinates,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  design$match <- as.numeric(c(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
    0, 0, 0
  ))
  design$weight <- seq_len(nrow(design)) + 1
  design
}

explicit_pair_lm_operator <- function(design, coefficient, spaces,
                                      encoding_nuisance = FALSE,
                                      retrieval_nuisance = FALSE) {
  compiled <- crossform:::.pair_design_matrix(
    design, spaces$left, spaces$right,
    encoding_nuisance, retrieval_nuisance
  )
  X <- compiled$X
  weight <- if ("weight" %in% names(design)) design$weight else rep(1, nrow(X))
  contrast <- setNames(rep(0, ncol(X)), colnames(X))
  contrast[[coefficient]] <- 1
  influence <- drop(contrast %*% solve(
    crossprod(X, weight * X), sweep(t(X), 2L, weight, `*`)
  ))
  H <- matrix(0, length(spaces$left$coordinates),
    length(spaces$right$coordinates))
  for (row in seq_len(nrow(design))) {
    H[compiled$left_index[[row]], compiled$right_index[[row]]] <-
      H[compiled$left_index[[row]], compiled$right_index[[row]]] +
      influence[[row]]
  }
  H
}

test_that("couplings support unequal axes, missing matches, and eligibility", {
  spaces <- pair_query_spaces()
  eligible <- expand.grid(
    left = spaces$left$coordinates,
    right = c("r1", "r2", "r4"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  matches <- match_coupling(
    c("e1", "e2", "e2"), c("r1", "r2", "r2"),
    spaces$left, spaces$right, eligible
  )
  controls <- control_coupling(matches)
  query <- coupling_contrast(matches, controls)

  expect_identical(dim(matches$value), c(3L, 4L))
  expect_equal(matches$value["e2", "r2"], 2)
  expect_equal(sum(matches$value["e3", ]), 0)
  expect_true(all(controls$value[!matches$eligible] == 0))
  expect_s3_class(query, "effect_pair_query")
  expect_identical(query$metadata$claim, "observed_operator_balance_only")
  expect_equal(sum(query$operator), 0, tolerance = 1e-15)
  expect_error(match_coupling("e1", "r3", spaces$left, spaces$right,
    eligible), "eligible", class = "effect_input_error")
})

test_that("pair LM coefficients equal explicit weighted least squares", {
  spaces <- pair_query_spaces()
  design <- full_pair_design(spaces)
  query <- pair_lm_query(
    design, "match", spaces$left, spaces$right,
    encoding_nuisance = TRUE, retrieval_nuisance = TRUE
  )
  oracle <- explicit_pair_lm_operator(
    design, "match", spaces,
    encoding_nuisance = TRUE, retrieval_nuisance = TRUE
  )

  expect_equal(unname(query$operator), oracle, tolerance = 1e-14)
  expect_equal(query$metadata$diagnostics$rank,
    length(query$metadata$diagnostics$columns))
  expect_true(query$metadata$diagnostics$balance$zero_row_marginals)
  expect_true(query$metadata$diagnostics$balance$zero_column_marginals)
  expect_true(query$metadata$diagnostics$balance$additive_baseline_invariant)
})

test_that("zero marginal H cancels arbitrary additive pair baselines exactly", {
  spaces <- pair_query_spaces()
  design <- full_pair_design(spaces)
  query <- pair_lm_query(
    design, "match", spaces$left, spaces$right,
    encoding_nuisance = TRUE, retrieval_nuisance = TRUE
  )
  left_baseline <- c(10, -4, 7)
  right_baseline <- c(3, 11, -8, 5)
  additive <- outer(left_baseline, rep(1, 4)) +
    outer(rep(1, 3), right_baseline)

  expect_equal(sum(as.matrix(query$operator) * additive), 0, tolerance = 1e-13)
})

test_that("downstream diagonal weighting can leak a balanced column margin", {
  spaces <- pair_query_spaces()
  design <- full_pair_design(spaces)
  query <- pair_lm_query(
    design, "match", spaces$left, spaces$right,
    encoding_nuisance = TRUE, retrieval_nuisance = TRUE
  )
  Q <- as.matrix(query$operator)
  a <- c(1, 2, 5)
  leaked <- diag(a) %*% Q

  expect_true(all(abs(colSums(Q)) < 1e-12))
  expect_gt(max(abs(colSums(leaked))), 1e-4)
  expect_identical(query$metadata$claim,
    "balance_reported_after_final_weighting")
})

test_that("sparse and dense H compile and query identically", {
  spaces <- pair_query_spaces()
  design <- full_pair_design(spaces)
  dense <- pair_lm_query(design, "match", spaces$left, spaces$right)
  sparse <- pair_lm_query(
    design, "match", spaces$left, spaces$right, sparse = TRUE
  )

  expect_s4_class(sparse$operator, "Matrix")
  expect_equal(as.matrix(sparse$operator), dense$operator, tolerance = 0)
  expect_silent(crossform:::.validate_query_for_compile(sparse))

  values <- matrix(seq_len(24), 2, 12)
  result <- crossform:::effect_form(
    values, spaces$left, spaces$right, receipt_fixture(), codec = "rectangular"
  )
  expect_equal(
    query_geometry(result, sparse)$values,
    query_geometry(result, dense)$values,
    tolerance = 0
  )
})

test_that("rank-deficient and infeasible pair designs fail clearly", {
  spaces <- pair_query_spaces()
  design <- full_pair_design(spaces)
  design$duplicate <- design$match
  expect_error(pair_lm_query(
    design, "match", spaces$left, spaces$right
  ), "rank deficient", class = "effect_input_error")

  infeasible <- design[design$left == "e1", c("left", "right", "match")]
  expect_error(pair_lm_query(
    infeasible, "match", spaces$left, spaces$right,
    encoding_nuisance = TRUE
  ), "infeasible", class = "effect_input_error")
  expect_error(pair_lm_query(
    design[, c("left", "right")], "match", spaces$left, spaces$right
  ), "predictor", class = "effect_input_error")
  expect_error(pair_lm_query(
    full_pair_design(spaces), "missing", spaces$left, spaces$right
  ), "Unknown", class = "effect_input_error")
})

test_that("match-control compiles balanced nuisance-adjusted operators", {
  spaces <- pair_query_spaces()
  matches <- match_coupling(
    c("e1", "e2", "e3"), c("r1", "r2", "r3"),
    spaces$left, spaces$right
  )
  query <- match_control(matches)

  expect_identical(query$metadata$constructor, "match_control")
  expect_true(query$metadata$diagnostics$balance$zero_row_marginals)
  expect_true(query$metadata$diagnostics$balance$zero_column_marginals)
  expect_identical(query$metadata$claim,
    "observed_operator_balance_not_downstream_invariance")
})
