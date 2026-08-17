# Contract tests for R/query-structured.R: the structured rank-one
# pair-difference query that lets RDM and RSA views be evaluated without
# materializing a dense q(q+1)/2-by-q(q-1)/2 operator.

structured_effects <- c("face", "house", "tool")

test_that("the default query enumerates every unordered effect pair once", {
  query <- crossform:::.pair_difference_query(structured_effects)

  expect_s3_class(query, "effect_pair_difference_query")
  expect_true(crossform:::.is_pair_difference_query(query))
  expect_false(crossform:::.is_pair_difference_query(matrix(1, 2, 2)))
  expect_identical(query$kind, "pair_differences")
  expect_identical(query$effects, structured_effects)
  expect_identical(query$pair_left, c(1L, 1L, 2L))
  expect_identical(query$pair_right, c(2L, 3L, 3L))
  expect_null(query$coefficients)
  expect_identical(query$pair_labels,
    c("face - house", "face - tool", "house - tool"))
})

test_that("pairs may be named or indexed and are canonically ordered", {
  named <- crossform:::.pair_difference_query(
    structured_effects, cbind(c("tool", "house"), c("face", "tool"))
  )
  indexed <- crossform:::.pair_difference_query(
    structured_effects, cbind(c(3L, 2L), c(1L, 3L))
  )

  # A pair is unordered, so `(tool, face)` is stored as `(face, tool)`.
  expect_identical(named$pair_left, c(1L, 2L))
  expect_identical(named$pair_right, c(3L, 3L))
  expect_identical(named$pair_left, indexed$pair_left)
  expect_identical(named$pair_right, indexed$pair_right)
  expect_identical(named$pair_labels, c("face - tool", "house - tool"))

  custom <- crossform:::.pair_difference_query(
    structured_effects, cbind(1L, 2L), labels = "near"
  )
  expect_identical(custom$pair_labels, "near")
})

test_that("a malformed pair specification is refused with its reason", {
  expect_error(
    crossform:::.pair_difference_query("only"),
    "at least two distinct named effects"
  , class = "effect_input_error")
  expect_error(
    crossform:::.pair_difference_query(c("a", "a")),
    "at least two distinct named effects"
  , class = "effect_input_error")
  expect_error(
    crossform:::.pair_difference_query(c(1, 2)),
    "at least two distinct named effects"
  , class = "effect_input_error")
  expect_error(
    crossform:::.pair_difference_query(structured_effects,
      cbind("face", "lure")),
    "`lure`.*is not a declared experimental effect"
  , class = "effect_input_error")
  expect_error(
    crossform:::.pair_difference_query(structured_effects, "face"),
    "must be a two-column matrix naming the effect pairs"
  , class = "effect_input_error")
  expect_error(
    crossform:::.pair_difference_query(structured_effects, cbind(1L, 9L)),
    "two-column matrix of effect names or of indices in 1:3"
  , class = "effect_input_error")
  expect_error(
    crossform:::.pair_difference_query(structured_effects, cbind(1.5, 2)),
    "two-column matrix of effect names or of indices in 1:3"
  , class = "effect_input_error")
  expect_error(
    crossform:::.pair_difference_query(structured_effects, cbind(2L, 2L)),
    "distance from `house` to itself, which is zero by construction"
  , class = "effect_input_error")
  expect_error(
    crossform:::.pair_difference_query(structured_effects,
      rbind(c(1L, 2L), c(2L, 1L))),
    "duplicate effect pairs"
  , class = "effect_input_error")
  expect_error(
    crossform:::.pair_difference_query(structured_effects, cbind(1L, 2L),
      labels = c("a", "b")),
    "Pair labels must name every pair"
  , class = "effect_input_error")
  expect_error(
    crossform:::.pair_difference_query(structured_effects,
      coefficients = matrix(1, 1L, 2L)),
    "finite output-by-pair matrix aligned"
  , class = "effect_input_error")
  expect_error(
    crossform:::.pair_difference_query(structured_effects,
      coefficients = matrix(c(1, NA, 1), 1L, 3L)),
    "finite output-by-pair matrix aligned"
  , class = "effect_input_error")
})

test_that("output width, labels, and payload follow the query representation", {
  pairs_only <- crossform:::.pair_difference_query(structured_effects)
  coefficients <- matrix(c(1, -1, 0, 0, 1, -1), 2L, 3L,
    dimnames = list(c("linear", "quadratic"), NULL))
  reduced <- crossform:::.pair_difference_query(
    structured_effects, coefficients = coefficients
  )
  unnamed <- crossform:::.pair_difference_query(
    structured_effects, coefficients = unname(coefficients)
  )
  dense <- matrix(0, 6L, 2L, dimnames = list(NULL, c("first", "second")))

  expect_identical(crossform:::.query_output_width(pairs_only), 3L)
  expect_identical(crossform:::.query_output_width(reduced), 2L)
  expect_identical(crossform:::.query_output_width(dense), 2L)
  expect_error(crossform:::.query_output_width(NULL),
    "A query is required to report an output width",
    class = "effect_input_error")

  expect_identical(crossform:::.query_output_labels(pairs_only),
    pairs_only$pair_labels)
  expect_identical(crossform:::.query_output_labels(reduced),
    c("linear", "quadratic"))
  expect_identical(crossform:::.query_output_labels(unnamed),
    c("view1", "view2"))
  expect_identical(crossform:::.query_output_labels(dense),
    c("first", "second"))
  expect_identical(crossform:::.query_output_labels(unname(dense)),
    c("view1", "view2"))

  # The structured representation is what makes this cheap: three pairs cost
  # two index vectors, not a dense 6-by-3 operator.
  expect_identical(crossform:::.query_payload_bytes(pairs_only), 8 * 6)
  expect_identical(crossform:::.query_payload_bytes(reduced), 8 * (6 + 6))
  expect_identical(crossform:::.query_payload_bytes(dense), 8 * 12)
  expect_identical(crossform:::.query_payload_bytes(NULL), 0)
  expect_lt(
    crossform:::.query_payload_bytes(pairs_only),
    crossform:::.query_payload_bytes(matrix(0, 6L, 3L))
  )
})

test_that("query identity is the estimand, not its labels or dimnames", {
  query <- crossform:::.pair_difference_query(structured_effects)
  relabelled <- crossform:::.pair_difference_query(
    structured_effects, labels = c("one", "two", "three")
  )
  reordered <- crossform:::.pair_difference_query(
    structured_effects, cbind(c(2L, 1L, 1L), c(3L, 3L, 2L))
  )

  expect_identical(
    crossform:::.query_identity_semantic(query),
    list(kind = "pair_differences", pair_left = c(1L, 1L, 2L),
      pair_right = c(2L, 3L, 3L), coefficients = NULL)
  )
  # Labels are presentation; the same pairs are the same estimand.
  expect_identical(
    crossform:::.query_identity_semantic(relabelled),
    crossform:::.query_identity_semantic(query)
  )
  # A different pair order is a different readout, so it is a different
  # identity even though the set of pairs matches.
  expect_false(identical(
    crossform:::.query_identity_semantic(reordered),
    crossform:::.query_identity_semantic(query)
  ))

  coefficients <- matrix(c(1, -1, 0), 1L, 3L,
    dimnames = list("contrast", c("p1", "p2", "p3")))
  expect_identical(
    crossform:::.query_identity_semantic(
      crossform:::.pair_difference_query(structured_effects,
        coefficients = coefficients)
    )$coefficients,
    unname(coefficients)
  )
  dense <- matrix(1:12, 6L, 2L, dimnames = list(NULL, c("a", "b")))
  expect_identical(crossform:::.query_identity_semantic(dense), unname(dense))
})

test_that("structured edge products equal the explicit rank-one contraction", {
  query <- crossform:::.pair_difference_query(structured_effects)
  left <- matrix(c(1, 4, 2, 0, -1, 3, 5, 1, -2, 2, 0, 1), 3L, 4L)
  right <- matrix(c(0, 1, 2, 3, -1, 1, 2, 0, 1, -2, 4, 1), 3L, 4L)
  oracle <- t(vapply(seq_along(query$pair_left), function(pair) {
    direction <- numeric(3L)
    direction[[query$pair_left[[pair]]]] <- 1
    direction[[query$pair_right[[pair]]]] <- -1
    # tr(u u' L diag(e_v) R') = (u'L)_v (u'R)_v for each feature v.
    drop(direction %*% left) * drop(direction %*% right)
  }, numeric(4L)))

  expect_identical(dim(oracle), c(3L, 4L))
  expect_equal(
    crossform:::.pair_difference_edge_products(query, left, right),
    oracle, tolerance = 1e-12
  )
  # A pair tile is exactly the corresponding rows, so tiling is lossless.
  expect_equal(
    crossform:::.pair_difference_edge_products(query, left, right, c(3L, 1L)),
    oracle[c(3L, 1L), , drop = FALSE], tolerance = 1e-12
  )
  expect_identical(
    dim(crossform:::.pair_difference_edge_products(query, left, right, 2L)),
    c(1L, 4L)
  )
})

test_that("pair-difference queries are admitted only on one shared axis", {
  query <- crossform:::.pair_difference_query(structured_effects)

  expect_identical(
    crossform:::.validate_pair_difference_for_task(
      query, structured_effects, structured_effects, TRUE
    ),
    query
  )
  # A distance compares two effects of one axis, so a rectangular task, a
  # mismatched right axis, or a query over other effects is refused.
  expect_error(
    crossform:::.validate_pair_difference_for_task(
      query, structured_effects, structured_effects, FALSE
    ),
    "self-form queries"
  , class = "effect_input_error")
  expect_error(
    crossform:::.validate_pair_difference_for_task(
      query, structured_effects, c("x", "y", "z"), TRUE
    ),
    "self-form queries"
  , class = "effect_input_error")
  expect_error(
    crossform:::.validate_pair_difference_for_task(
      query, c("face", "house", "lure"), c("face", "house", "lure"), TRUE
    ),
    "self-form queries"
  , class = "effect_input_error")
})
