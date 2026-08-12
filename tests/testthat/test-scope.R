test_that("collapse is admitted only for additive fixed-bilinear work", {
  at <- additive_frame(matrix(c(1, 0, 0.5, 0.5), nrow = 2, byrow = TRUE))
  query <- bilinear_query(diag(2))

  plan <- compile_lowering(at, query)

  expect_true(plan$collapsed)
  expect_identical(plan$kind, "additive_contraction")
})

test_that("factor, local, adaptive, and nonlinear work have separate lowerings", {
  additive <- additive_frame(matrix(1, nrow = 1, ncol = 2))
  factor <- factor_frame(list(diag(2)))
  local_factor <- factor_frame(list(diag(2)), locally_estimated = TRUE)

  expect_identical(
    compile_lowering(factor, bilinear_query(diag(2)))$kind,
    "factor_contraction"
  )
  expect_identical(
    compile_lowering(local_factor, bilinear_query(diag(2)))$kind,
    "location_dependent_fit"
  )
  expect_identical(
    compile_lowering(additive, bilinear_query(diag(2), fixed = FALSE))$kind,
    "adaptive_query"
  )
  expect_identical(
    compile_lowering(additive, nonlinear_query(function(x) max(x)))$kind,
    "nonlinear_readout"
  )
})

test_that("additive frame rejects values outside theorem preconditions", {
  expect_error(additive_frame(matrix(c(1, -1), nrow = 1)), "nonnegative")
  expect_error(additive_frame(matrix(c(1, Inf), nrow = 1)), "finite")
})

test_that("compiler revalidates mutated and forged objects", {
  frame <- additive_frame(matrix(c(1, 0, 0, 1), 2), normalization = "local")
  frame$weights[1, 1] <- -1
  expect_error(compile_lowering(frame, bilinear_query(diag(2))), "nonnegative")

  frame <- additive_frame(diag(2), normalization = "local")
  frame$fixed <- FALSE
  expect_error(compile_lowering(frame, bilinear_query(diag(2))), "must be fixed")

  frame <- additive_frame(diag(2), normalization = "local")
  frame$domain_id <- ""
  expect_error(compile_lowering(frame, bilinear_query(diag(2))), "domain_id")

  factor <- factor_frame(list(diag(2)))
  factor$factors[[1]][1, 1] <- Inf
  expect_error(compile_lowering(factor, bilinear_query(diag(2))), "finite")

  query <- bilinear_query(diag(2))
  query$operator[1, 2] <- 1
  expect_error(compile_lowering(additive_frame(diag(2)), query), "symmetric")
})

test_that("sparse additive frames have an explicit collapse contract", {
  weights <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = 1, dims = c(3, 3))
  frame <- additive_frame(weights, normalization = "local", domain_id = "mesh:v1")
  plan <- compile_lowering(frame, bilinear_query(diag(2)))

  expect_true(inherits(frame$weights, "sparseMatrix"))
  expect_true(plan$collapsed)
})

test_that("normalization policies are revalidated at compilation", {
  expect_error(additive_frame(matrix(c(1, 1), 1), normalization = "local"),
    "sum to one")
  expect_error(additive_frame(matrix(c(1, 0), 1), normalization = "conservative"),
    "positive mass|columns must sum")
})
