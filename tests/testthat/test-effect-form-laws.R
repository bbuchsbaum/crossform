# The executable law court for `effect-form-v1`
# (design/effect-form-contract.md). Every package value below is compared to
# an independent oracle in helper-effect-form-laws.R, which never calls
# production crossform code. Tolerances are explicit and never looser than
# 1e-10.

effect_form_law_tolerance <- 1e-10

# One relation over one weighted spatial node, so each package reading is a
# single oracle expression rather than a nested aggregate. `left` and `right`
# name distinct experimental axes over the same neural domain, giving both a
# self form (packed codec) and a rectangular form (vec codec) from the same
# numbers.
effect_form_law_fixture <- local({
  built <- NULL
  function() {
    if (!is.null(built)) {
      return(built)
    }
    set.seed(20260818L)
    p <- 6L
    domain <- abstract_domain(p, id = "laws:effect-form-neural:v1")
    left_space <- effect_space(c("a", "b", "c"),
      basis_id = "laws:effect-form-left:v1")
    right_space <- effect_space(c("x", "y"),
      basis_id = "laws:effect-form-right:v1")
    make <- function(q, space) {
      matrix(rnorm(q * p), q, p, dimnames = list(space$coordinates, NULL))
    }
    left_values <- list(run1 = make(3L, left_space), run2 = make(3L, left_space))
    right_values <- list(
      rrun1 = make(2L, right_space), rrun2 = make(2L, right_space)
    )
    weights <- c(0.5, 1, 1.5, 0.25, 2, 0.75)
    left <- relation(left_values, effects = left_space, domain = domain)
    right <- relation(right_values, effects = right_space, domain = domain)
    built <<- list(
      p = p, domain = domain,
      left_space = left_space, right_space = right_space,
      left_values = left_values, right_values = right_values,
      weights = weights,
      left = left, right = right,
      frame = additive_frame(matrix(weights, 1L, p), domain = domain),
      self_over = cross_partitions(left, independence = "independent"),
      cross_over = pairing("run1", "rrun1", 1,
        directed = TRUE, independence = "independent")
    )
    built
  }
})

effect_form_law_symmetrize <- function(value) (value + t(value)) / 2

# The self form the package materializes: one unordered partition pair,
# expanded into self-adjoint half edges, hence the symmetrization.
effect_form_law_self <- function(fixture) {
  materialize_geometry(plan_geometry(
    fixture$left, fixture$frame, fixture$self_over
  ))
}

effect_form_law_rectangular <- function(fixture) {
  materialize_geometry(plan_geometry(
    fixture$left, fixture$frame, fixture$cross_over, right = fixture$right
  ))
}


# Section 5: storage codecs -------------------------------------------------

test_that("the packed symmetric codec equals the oracle svec", {
  fixture <- effect_form_law_fixture()
  symmetric <- effect_form_law_symmetrize(oracle_effect_form(
    fixture$left_values$run1, fixture$left_values$run2, fixture$weights
  ))
  packed <- oracle_svec(symmetric)

  # The package codec and the oracle agree entry by entry, in the same order.
  expect_identical(length(packed), 6L)
  expect_equal(crossform:::.svec_symmetric(symmetric), packed,
    tolerance = effect_form_law_tolerance)
  # It is a lossless isometry: the round trip restores the matrix, and the
  # packed inner product is the Frobenius inner product.
  expect_equal(crossform:::.unsvec_symmetric(packed, 3L), symmetric,
    tolerance = effect_form_law_tolerance, ignore_attr = TRUE)
  other <- effect_form_law_symmetrize(oracle_effect_form(
    fixture$left_values$run2, fixture$left_values$run1, fixture$weights
  ))
  expect_equal(sum(packed * oracle_svec(other)), sum(symmetric * other),
    tolerance = effect_form_law_tolerance)
  # The off-diagonal sqrt(2) scaling is what makes that isometry hold, so a
  # naive lower-triangle read is a different vector.
  expect_false(isTRUE(all.equal(
    packed, symmetric[lower.tri(symmetric, diag = TRUE)], tolerance = 1e-6
  )))

  # And the materialized self form stores exactly that packed row.
  form <- effect_form_law_self(fixture)
  expect_identical(form$codec, "symmetric_packed")
  expect_equal(as.vector(geometry_component(form, "total")), packed,
    tolerance = effect_form_law_tolerance)
})

test_that("the rectangular codec is column-major vec of the oracle form", {
  fixture <- effect_form_law_fixture()
  form <- effect_form_law_rectangular(fixture)
  expected <- oracle_effect_form(
    fixture$left_values$run1, fixture$right_values$rrun1, fixture$weights
  )

  expect_identical(form$codec, "rectangular")
  expect_equal(as.vector(geometry_component(form, "total")),
    oracle_vec(expected), tolerance = effect_form_law_tolerance)
  expect_identical(
    length(oracle_vec(expected)),
    length(fixture$left_space$coordinates) *
      length(fixture$right_space$coordinates)
  )
})


# Section 7: the exact centered sufficient statistic -----------------------

test_that("total, coherent, and configuration each equal the oracle", {
  fixture <- effect_form_law_fixture()
  form <- effect_form_law_self(fixture)
  statistics <- oracle_centered_statistics(
    fixture$left_values$run1, fixture$left_values$run2, fixture$weights
  )
  packed <- function(value) oracle_svec(effect_form_law_symmetrize(value))

  # Each component is checked against its own first-principles oracle, so no
  # assertion here is satisfied by the identity configuration = total - coherent.
  expect_equal(as.vector(geometry_component(form, "total")),
    packed(statistics$total), tolerance = effect_form_law_tolerance)
  expect_equal(as.vector(geometry_component(form, "coherent")),
    packed(statistics$coherent), tolerance = effect_form_law_tolerance)
  expect_equal(as.vector(geometry_component(form, "configuration")),
    packed(statistics$centered), tolerance = effect_form_law_tolerance)

  # Coherent is the frame-weighted mean pattern outer product over the
  # frame mass, exactly as the contract defines it.
  expect_equal(statistics$mass, sum(fixture$weights),
    tolerance = effect_form_law_tolerance)
  expect_equal(statistics$coherent,
    outer(statistics$left_first, statistics$right_first) / statistics$mass,
    tolerance = effect_form_law_tolerance)
  expect_equal(statistics$left_first,
    as.vector(fixture$left_values$run1 %*% fixture$weights),
    tolerance = effect_form_law_tolerance)

  # And configuration is the weight-centered pattern cross-product computed
  # without ever forming `total`.
  expect_equal(as.vector(geometry_component(form, "configuration")),
    packed(oracle_explicit_center(
      fixture$left_values$run1, fixture$left_values$run2, fixture$weights
    )),
    tolerance = effect_form_law_tolerance)
  expect_equal(statistics$covariance, statistics$centered / statistics$mass,
    tolerance = effect_form_law_tolerance)
})

test_that("a singleton positive-mass frame has exactly zero configuration", {
  fixture <- effect_form_law_fixture()
  singleton <- additive_frame(
    matrix(c(0, 0, 1.5, 0, 0, 0), 1L, fixture$p), domain = fixture$domain
  )
  form <- materialize_geometry(plan_geometry(
    fixture$left, singleton, fixture$self_over
  ))
  statistics <- oracle_centered_statistics(
    fixture$left_values$run1, fixture$left_values$run2,
    c(0, 0, 1.5, 0, 0, 0)
  )

  expect_equal(max(abs(statistics$centered)), 0,
    tolerance = effect_form_law_tolerance)
  expect_equal(max(abs(geometry_component(form, "configuration"))), 0,
    tolerance = effect_form_law_tolerance)
  expect_equal(
    as.vector(geometry_component(form, "total")),
    as.vector(geometry_component(form, "coherent")),
    tolerance = effect_form_law_tolerance
  )
})


# Section 7: fixed linear queries ------------------------------------------

test_that("querying a form equals the oracle Frobenius pairing", {
  fixture <- effect_form_law_fixture()
  self_form <- effect_form_law_self(fixture)
  rectangular <- effect_form_law_rectangular(fixture)
  symmetric <- effect_form_law_symmetrize(oracle_effect_form(
    fixture$left_values$run1, fixture$left_values$run2, fixture$weights
  ))
  contrast <- c(1, -1, 0)
  operator <- tcrossprod(contrast)

  expect_equal(
    unname(drop(
      query_geometry(self_form, bilinear_query(operator), "total")$values
    )),
    oracle_query(symmetric, operator),
    tolerance = effect_form_law_tolerance
  )

  # The bound rectangular query carries axis identity, not just a shape.
  left_space <- list(
    id = fixture$left_space$basis_id,
    coordinates = fixture$left_space$coordinates
  )
  right_space <- list(
    id = fixture$right_space$basis_id,
    coordinates = fixture$right_space$coordinates
  )
  rectangular_operator <- matrix(c(1, -0.5, 0, 0, 0.25, -1), 3L, 2L)
  bound <- oracle_bind_query(rectangular_operator, left_space, right_space)
  expected <- oracle_apply_bound_query(
    oracle_effect_form(
      fixture$left_values$run1, fixture$right_values$rrun1, fixture$weights
    ),
    bound, left_space, right_space
  )
  expect_equal(
    unname(drop(query_geometry(rectangular, pair_query(
      rectangular_operator, fixture$left_space, fixture$right_space
    ))$values)),
    expected, tolerance = effect_form_law_tolerance
  )

  # Both the package and the oracle refuse a same-shaped query bound to a
  # different experimental identity.
  impostor <- effect_space(
    fixture$left_space$coordinates, basis_id = "laws:impostor:v1"
  )
  expect_error(
    query_geometry(rectangular, pair_query(
      rectangular_operator, impostor, fixture$right_space
    )),
    "The query and form axis identities are incompatible"
  , class = "effect_contract_error")
  expect_error(
    oracle_apply_bound_query(
      symmetric, bound, list(id = "laws:impostor:v1",
        coordinates = fixture$left_space$coordinates), right_space
    ),
    "axis identities differ"
  )
  expect_error(
    oracle_bind_query(rectangular_operator, right_space, left_space),
    "do not match their bound spaces"
  )
})


# Section 5, 10: capabilities and the result contract ----------------------

test_that("declared capabilities satisfy the oracle implication lattice", {
  fixture <- effect_form_law_fixture()
  self_form <- effect_form_law_self(fixture)
  rectangular <- effect_form_law_rectangular(fixture)
  read <- function(form) {
    oracle_capabilities(
      form$capabilities$self_form,
      form$capabilities$symmetric,
      form$capabilities$guaranteed_psd
    )
  }

  expect_identical(read(self_form),
    c(self_form = TRUE, symmetric = TRUE, guaranteed_psd = FALSE))
  expect_identical(read(rectangular),
    c(self_form = FALSE, symmetric = FALSE, guaranteed_psd = FALSE))
  # A crossvalidated self form is symmetric and indefinite, never PSD, so the
  # oracle's `guaranteed_psd implies symmetric` rule is not vacuous here.
  expect_false(self_form$capabilities$guaranteed_psd)
  expect_error(oracle_capabilities(TRUE, FALSE, TRUE),
    "guaranteed_psd implies symmetric")
  expect_error(oracle_capabilities(FALSE, TRUE, FALSE),
    "symmetric implies self_form")

  # The codec each form declares is the one the oracle contract permits.
  expect_identical(
    oracle_result_contract(
      read(self_form), "symmetric_packed", "complete_form"
    )$codec,
    self_form$codec
  )
  expect_identical(
    oracle_result_contract(
      read(rectangular), "rectangular", "complete_form"
    )$codec,
    rectangular$codec
  )
  expect_identical(self_form$result_capability, "complete_form")
  expect_error(
    oracle_result_contract(read(rectangular), "symmetric_packed",
      "complete_form"),
    "packed codec requires a symmetric guarantee"
  )

  # A query-only view is a different completeness, and has no form to read.
  view <- evaluate_geometry(
    plan_geometry(fixture$left, fixture$frame, fixture$self_over),
    query = bilinear_query(tcrossprod(c(1, -1, 0)))
  )
  # `$result_capability` is the field that carries the oracle's completeness
  # vocabulary on both result kinds; `$completeness` says "full" rather than
  # "complete_form" on a materialized form, so it is not the same word.
  expect_identical(
    view$result_capability,
    oracle_result_contract(read(self_form), "symmetric_packed",
      "query_only")$completeness
  )
  expect_identical(view$completeness, "query_only")
  expect_identical(
    self_form$result_capability,
    oracle_result_contract(read(self_form), "symmetric_packed",
      "complete_form")$completeness
  )
  expect_s3_class(view, "effect_view")
  expect_false(inherits(view, "effect_geometry"))
})


# Section 6: nonlinear edge transforms -------------------------------------

test_that("the Fisher transform equals the oracle under both boundary policies", {
  values <- matrix(c(0.1, -0.5, 0.87, 0.3), 2L, 2L)
  boundary <- matrix(c(1, -1, 0.5, 0.2), 2L, 2L)

  expect_equal(
    crossform:::.apply_edge_transform(values, fisher_z("error")),
    oracle_fisher(values, "error"),
    tolerance = effect_form_law_tolerance
  )
  expect_equal(
    crossform:::.apply_edge_transform(boundary, fisher_z("clip", 1e-6)),
    oracle_fisher(boundary, "clip", 1e-6),
    tolerance = effect_form_law_tolerance
  )
  # A correlation of exactly one is refused rather than silently clipped, and
  # both sides refuse it in the same place.
  expect_error(
    crossform:::.apply_edge_transform(boundary, fisher_z("error")),
    "absolute-correlation boundary"
  , class = "effect_input_error")
  expect_error(oracle_fisher(boundary, "error"), "Fisher boundary")
  expect_error(oracle_fisher(values, "clip", delta = 0), "clipping delta")
  expect_error(fisher_z("clip", delta = 0), "Fisher boundary policy",
    class = "effect_input_error")

  # Clipping is not the identity: the policy and its delta change the number,
  # so they belong in plan identity.
  expect_false(isTRUE(all.equal(
    oracle_fisher(boundary, "clip", 1e-6),
    oracle_fisher(boundary, "clip", 1e-3),
    tolerance = 1e-6
  )))
  expect_false(identical(
    crossform:::.edge_operation_plan(
      crossform:::correlation(), fisher_z("clip", 1e-6)
    )$signature,
    crossform:::.edge_operation_plan(
      crossform:::correlation(), fisher_z("clip", 1e-3)
    )$signature
  ))
  # Fisher is defined on correlations, so an uncentered inner product input
  # is refused at plan time rather than transformed anyway.
  expect_error(
    crossform:::.edge_operation_plan(inner_product(), fisher_z("error")),
    "requires correlation-valued edge input"
  , class = "effect_input_error")
})


# Section 9 of evidence-pairing-v1: the factorized bridge specialization ----

test_that("a measurement bridge factorizes the neural query exactly", {
  set.seed(20260818L + 1L)
  fixture <- effect_form_law_fixture()
  left_leg <- matrix(rnorm(2L * fixture$p), 2L, fixture$p)
  right_leg <- matrix(rnorm(2L * fixture$p), 2L, fixture$p)
  bridge <- measurement_bridge(
    left_leg, right_leg, fixture$domain, fixture$domain,
    measurement_space(2L, id = "laws:effect-form-common:v1")
  )
  applied <- crossform:::.apply_measurement_bridge(
    list(run1 = fixture$left_values$run1),
    list(rrun1 = fixture$right_values$rrun1), bridge
  )
  got <- applied$left$run1 %*% t(applied$right$rrun1)

  # Measuring both sides into the common space and contracting the induced
  # operator directly are the same form.
  expect_equal(
    got,
    oracle_bridge_form(
      fixture$left_values$run1, fixture$right_values$rrun1,
      left_leg, right_leg
    ),
    tolerance = effect_form_law_tolerance
  )
  # The second oracle contracts the induced operator K = L_left' L_right in
  # one step; the package value must equal that too, not just the two-step
  # measured form above.
  expect_equal(
    got,
    oracle_direct_bridge_form(
      fixture$left_values$run1, fixture$right_values$rrun1,
      left_leg, right_leg
    ),
    tolerance = effect_form_law_tolerance
  )
  expect_equal(
    t(left_leg) %*% right_leg, t(bridge$left_leg) %*% bridge$right_leg,
    tolerance = effect_form_law_tolerance
  )
})

test_that("distinct neural spaces require a bridge, as the oracle demands", {
  set.seed(20260818L + 2L)
  fixture <- effect_form_law_fixture()
  other_domain <- abstract_domain(4L, id = "laws:other-neural:v1")
  other <- relation(
    list(rrun1 = matrix(rnorm(8L), 2L, 4L,
      dimnames = list(c("x", "y"), NULL))),
    effects = fixture$right_space, domain = other_domain
  )
  over <- pairing("run1", "rrun1", 1,
    directed = TRUE, independence = "independent")

  expect_silent(oracle_require_bridge(fixture$domain$id, fixture$domain$id))
  expect_error(
    oracle_require_bridge(fixture$domain$id, other_domain$id),
    "distinct neural spaces require a bridge"
  )
  expect_error(
    crossform:::.compile_effect_task(fixture$left, over, other),
    "require an explicit"
  , class = "effect_contract_error")

  bridge <- measurement_bridge(
    matrix(rnorm(2L * fixture$p), 2L, fixture$p),
    matrix(rnorm(2L * 4L), 2L, 4L),
    fixture$domain, other_domain,
    measurement_space(2L, id = "laws:other-common:v1")
  )
  expect_silent(oracle_require_bridge(
    fixture$domain$id, other_domain$id, bridge
  ))
  task <- crossform:::.compile_effect_task(
    fixture$left, over, other, bridge = bridge
  )
  expect_identical(task$bridge$signature, bridge$signature)
})


# The accepted correlation-distance boundary -------------------------------

test_that("rdm() is squared distance and never a silent correlation distance", {
  fixture <- effect_form_law_fixture()
  form <- effect_form_law_self(fixture)
  geometry <- effect_form_law_symmetrize(oracle_effect_form(
    fixture$left_values$run1, fixture$left_values$run2, fixture$weights
  ))
  pairs <- utils::combn(3L, 2L)
  squared <- vapply(seq_len(ncol(pairs)), function(pair) {
    delta <- numeric(3L)
    delta[pairs[, pair]] <- c(1, -1)
    oracle_query(geometry, tcrossprod(delta))
  }, numeric(1))

  expect_equal(as.vector(rdm(form)$values), squared,
    tolerance = effect_form_law_tolerance)

  # The oracle correlation distance is a different, nonlinear quantity, and
  # `rdm()` must not be quietly reporting it (see
  # vignettes/correlation-distance-policy.Rmd).
  correlations <- oracle_weighted_correlation(
    fixture$left_values$run1, fixture$left_values$run1, fixture$weights
  )
  correlation_distance <- vapply(seq_len(ncol(pairs)), function(pair) {
    1 - correlations[pairs[1L, pair], pairs[2L, pair]]
  }, numeric(1))
  expect_equal(unname(diag(correlations)), rep(1, 3L),
    tolerance = effect_form_law_tolerance)
  expect_true(all(abs(correlations) <= 1 + effect_form_law_tolerance))
  expect_false(isTRUE(all.equal(
    squared, correlation_distance, tolerance = 1e-3
  )))

  # The oracle's zero-variance policy is the behavior a future named view
  # would have to implement; nothing in crossform 0.1 may produce it silently.
  flat <- fixture$left_values$run1
  flat[2L, ] <- 3
  expect_error(
    oracle_weighted_correlation(flat, flat, fixture$weights),
    "zero weighted variance"
  )
  tolerant <- oracle_weighted_correlation(
    flat, flat, fixture$weights, zero_variance = "zero"
  )
  expect_equal(unname(tolerant[2L, ]), c(0, 0, 0), tolerance = 0)
  expect_false("correlation_rdm" %in% getNamespaceExports("crossform"))
})
