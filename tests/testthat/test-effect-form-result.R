rectangular_form_fixture <- function(codec = "rectangular") {
  left <- effect_space(c("encoding_a", "encoding_b"),
    basis_id = "encoding:v1")
  right <- effect_space(c("retrieval_a", "retrieval_b", "retrieval_c"),
    basis_id = "retrieval:v1")
  forms <- list(
    matrix(c(1, 4, -2, 3, 5, 0), nrow = 2),
    matrix(c(2, -1, 6, 4, 0, 3), nrow = 2)
  )
  stored <- do.call(rbind, lapply(forms, function(value) as.vector(value)))
  list(
    left = left,
    right = right,
    forms = forms,
    result = crossform:::effect_form(
      total = stored,
      left_space = left,
      right_space = right,
      receipt = receipt_fixture(),
      index = c("m1", "m2"),
      codec = codec
    )
  )
}

test_that("pair queries bind independent ordered axes", {
  fixture <- rectangular_form_fixture()
  H <- matrix(c(1, -2, 3, 0, 4, -1), nrow = 2)
  query <- pair_query(H, fixture$left, fixture$right)

  expect_s3_class(query, "effect_pair_query")
  expect_identical(query$left_space, fixture$left)
  expect_identical(query$right_space, fixture$right)
  expect_identical(query$operator, H)

  wrong_left <- effect_space(fixture$left$coordinates, basis_id = "encoding:v2")
  wrong_query <- pair_query(H, wrong_left, fixture$right)
  expect_error(
    query_geometry(fixture$result, wrong_query),
    "axis identities"
  , class = "effect_contract_error")
  expect_error(
    query_geometry(
      fixture$result,
      pair_query(t(H), fixture$right, fixture$left)
    ),
    "axis identities"
  , class = "effect_contract_error")
  expect_error(pair_query(matrix(1, 3, 2), fixture$left, fixture$right),
    "dimension", class = "effect_contract_error")
})

test_that("complete rectangular forms answer pair queries by column-major vec", {
  fixture <- rectangular_form_fixture()
  H <- matrix(c(1, -2, 3, 0, 4, -1), nrow = 2)
  view <- query_geometry(
    fixture$result,
    pair_query(H, fixture$left, fixture$right)
  )
  oracle <- vapply(fixture$forms, function(value) sum(value * H), numeric(1))

  expect_s3_class(fixture$result, "effect_form")
  expect_identical(fixture$result$logical_shape, c(2L, 3L))
  expect_identical(fixture$result$codec, "rectangular")
  expect_false(fixture$result$capabilities$self_form)
  expect_false(fixture$result$capabilities$symmetric)
  expect_equal(drop(view$values), oracle, tolerance = 0)
  expect_identical(view$result_capability, "query_only")
  expect_null(view$effect_space)
  expect_identical(view$left_space, fixture$left)
  expect_identical(view$right_space, fixture$right)
})

test_that("independent left and right permutations are equivariant", {
  fixture <- rectangular_form_fixture()
  H <- matrix(c(1, -2, 3, 0, 4, -1), nrow = 2)
  left_order <- c(2, 1)
  right_order <- c(3, 1, 2)
  permuted_left <- effect_space(
    fixture$left$coordinates[left_order], basis_id = fixture$left$basis_id
  )
  permuted_right <- effect_space(
    fixture$right$coordinates[right_order], basis_id = fixture$right$basis_id
  )
  permuted_forms <- lapply(fixture$forms, function(value) {
    value[left_order, right_order, drop = FALSE]
  })
  permuted <- crossform:::effect_form(
    total = do.call(rbind, lapply(permuted_forms, as.vector)),
    left_space = permuted_left,
    right_space = permuted_right,
    receipt = receipt_fixture(),
    codec = "rectangular"
  )

  original_value <- query_geometry(
    fixture$result, pair_query(H, fixture$left, fixture$right)
  )$values
  permuted_value <- query_geometry(
    permuted,
    pair_query(H[left_order, right_order, drop = FALSE],
      permuted_left, permuted_right)
  )$values

  expect_equal(permuted_value, original_value, tolerance = 0)
})

test_that("packed and rectangular codecs contract symmetric forms identically", {
  space <- effect_space(c("a", "b", "c"), basis_id = "self:v1")
  forms <- list(
    crossprod(matrix(c(1, 2, 0, -1, 3, 4), nrow = 2)),
    crossprod(matrix(c(2, 0, 1, 5, -2, 3), nrow = 2))
  )
  packed_values <- do.call(rbind, lapply(forms, oracle_svec))
  rectangular_values <- do.call(rbind, lapply(forms, as.vector))
  zero_packed <- matrix(0, nrow(packed_values), ncol(packed_values))
  zero_rectangular <- matrix(0, nrow(rectangular_values), ncol(rectangular_values))
  packed <- crossform:::effect_form(
    packed_values, space, space, receipt_fixture(),
    codec = "symmetric_packed", symmetric = TRUE, coherent = zero_packed
  )
  rectangular <- crossform:::effect_form(
    rectangular_values, space, space, receipt_fixture(),
    codec = "rectangular", symmetric = TRUE, coherent = zero_rectangular
  )
  H <- matrix(c(1, 4, -2, 0, 3, 5, 2, -1, 6), nrow = 3)
  query <- pair_query(H, space, space)

  expect_equal(
    query_geometry(packed, query)$values,
    query_geometry(rectangular, query)$values,
    tolerance = 1e-14
  )
  expect_equal(
    rdm(packed)$values,
    rdm(rectangular)$values,
    tolerance = 1e-14
  )
  expect_equal(
    geometry_spectrum(packed)$values,
    geometry_spectrum(rectangular)$values,
    tolerance = 1e-14
  )
})

test_that("complete-form and query-only claims cannot be forged", {
  fixture <- rectangular_form_fixture()
  query <- pair_query(matrix(1, 2, 3), fixture$left, fixture$right)
  view <- query_geometry(fixture$result, query)

  forged <- fixture$result
  forged$result_capability <- "query_only"
  expect_error(crossform:::.validate_effect_form(forged), "complete effect_form",
    class = "effect_input_error")

  forged <- fixture$result
  forged$capabilities$symmetric <- TRUE
  expect_error(crossform:::.validate_effect_form(forged),
    "symmetric effect form|capabilities", class = "effect_input_error")

  forged <- fixture$result
  forged$total$manifest$format <- "packed-double-v1"
  expect_error(crossform:::.validate_effect_form(forged), "manifest",
    class = "effect_input_error")

  self <- result_fixture()
  forged_self <- self
  forged_self$capabilities$guaranteed_psd <- TRUE
  expect_error(crossform:::.validate_effect_form(forged_self),
    "contract signature", class = "effect_contract_error")

  forged_view <- view
  forged_view$result_capability <- "complete_form"
  class(forged_view) <- c("effect_form", class(forged_view))
  expect_error(crossform:::.validate_effect_view(forged_view), "query-only",
    class = "effect_input_error")

  forged_view <- view
  forged_view$query$operator[1, 1] <- 2
  expect_error(crossform:::.validate_effect_view(forged_view),
    "contract signature", class = "effect_contract_error")
})

test_that("self-only and component views have explicit capability guards", {
  fixture <- rectangular_form_fixture()

  rdm_refusal <- catch_refusal(rdm(fixture$result))
  expect_s3_class(rdm_refusal, "effect_capability_refusal")
  expect_identical(rdm_refusal$capability, "symmetric_self_form")
  expect_identical(rdm_refusal$namespace, "geometry_views")
  expect_identical(rdm_refusal$reasons, "form_is_not_a_symmetric_self_form")
  expect_match(conditionMessage(rdm_refusal), "symmetric self form")

  rsa_refusal <- catch_refusal(
    rsa(fixture$result, models = list(model = diag(3)))
  )
  expect_s3_class(rsa_refusal, "effect_capability_refusal")
  expect_identical(rsa_refusal$capability, "symmetric_self_form")
  expect_identical(rsa_refusal$namespace, "geometry_views")
  expect_identical(rsa_refusal$reasons, "form_is_not_a_symmetric_self_form")
  expect_match(conditionMessage(rsa_refusal), "symmetric self form")

  spectrum_refusal <- catch_refusal(geometry_spectrum(fixture$result))
  expect_s3_class(spectrum_refusal, "effect_capability_refusal")
  expect_identical(spectrum_refusal$capability, "symmetric_self_form")
  expect_identical(spectrum_refusal$namespace, "geometry_views")
  expect_identical(spectrum_refusal$reasons,
    "form_is_not_a_symmetric_self_form")
  expect_match(conditionMessage(spectrum_refusal), "symmetric self form")

  expect_error(
    geometry_component(fixture$result, "configuration"),
    "does not carry"
  , class = "effect_input_error")
})

test_that("rectangular block stores preserve codec and bounded queries", {
  fixture <- rectangular_form_fixture()
  stored <- geometry_component(fixture$result)
  path <- tempfile(fileext = ".egm")
  on.exit(unlink(path), add = TRUE)
  store <- crossform:::.file_geometry_store(
    path, dim(stored), create = TRUE, codec = "rectangular"
  )
  crossform:::.write_geometry_tile(
    store, seq_len(nrow(stored)), seq_len(ncol(stored)), stored
  )
  blocked <- crossform:::effect_form(
    store, fixture$left, fixture$right, receipt_fixture(),
    codec = "rectangular"
  )
  query <- pair_query(
    matrix(c(1, -2, 3, 0, 4, -1), nrow = 2),
    fixture$left,
    fixture$right
  )

  expect_identical(store$manifest$format, "rectangular-double-v1")
  expect_equal(
    query_geometry(blocked, query, row_block = 1)$values,
    query_geometry(fixture$result, query)$values,
    tolerance = 0
  )
})
