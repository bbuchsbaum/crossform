# The executable law court for `evidence-pairing-v1`
# (design/evidence-pairing-contract.md, section 14). Every assertion here
# compares a production crossform surface against the independent oracles in
# helper-evidence-pairing-laws.R, which never call package code. Tolerances
# are explicit and never looser than 1e-10, except at one rank-one
# pseudo-inverse noted at its call site; fixtures are seeded and small enough
# that a failure reproduces from the printed dimensions alone.

pairing_law_tolerance <- 1e-10

# Four independently identified spaces, one partition per side, so every
# package value is one oracle expression rather than a weighted aggregate.
pairing_law_fixture <- local({
  built <- NULL
  function() {
    if (!is.null(built)) {
      return(built)
    }
    set.seed(20260816L)
    q_left <- 4L
    q_right <- 3L
    p_left <- 5L
    p_right <- 4L
    left_domain <- abstract_domain(p_left, id = "laws:left-neural:v1")
    right_domain <- abstract_domain(p_right, id = "laws:right-neural:v1")
    left_space <- effect_space(paste0("l", seq_len(q_left)),
      basis_id = "laws:left-effects:v1")
    right_space <- effect_space(paste0("r", seq_len(q_right)),
      basis_id = "laws:right-effects:v1")
    left_values <- matrix(rnorm(q_left * p_left), q_left, p_left,
      dimnames = list(left_space$coordinates, NULL))
    right_values <- matrix(rnorm(q_right * p_right), q_right, p_right,
      dimnames = list(right_space$coordinates, NULL))
    left_relation <- relation(list(left_run = left_values),
      effects = left_space, domain = left_domain)
    right_relation <- relation(list(right_run = right_values),
      effects = right_space, domain = right_domain)
    left_leg <- matrix(rnorm(2L * p_left), 2L, p_left)
    right_leg <- matrix(rnorm(3L * p_right), 3L, p_right)
    # A second right-hand relation over the *left* neural domain, for the
    # forward-transport laws, which close one shared neural boundary.
    shared_values <- matrix(rnorm(q_right * p_left), q_right, p_left,
      dimnames = list(right_space$coordinates, NULL))
    built <<- list(
      q_left = q_left, q_right = q_right,
      p_left = p_left, p_right = p_right,
      left_domain = left_domain, right_domain = right_domain,
      left_space = left_space, right_space = right_space,
      left_values = left_values, right_values = right_values,
      left_relation = left_relation, right_relation = right_relation,
      shared_values = shared_values,
      shared_relation = relation(list(right_run = shared_values),
        effects = right_space, domain = left_domain),
      left_leg = left_leg, right_leg = right_leg,
      h = matrix(rnorm(q_left * q_right), q_left, q_right),
      over = pairing("left_run", "right_run", 1,
        directed = TRUE, independence = "independent"),
      identity_left = measurement_frame(
        list(native_left = diag(p_left)), domain = left_domain,
        id = "laws:identity-left:v1"
      ),
      identity_right = measurement_frame(
        list(native_right = diag(p_right)), domain = right_domain,
        id = "laws:identity-right:v1"
      ),
      leg_left = measurement_frame(
        list(node_x = left_leg), domain = left_domain,
        id = "laws:left-legs:v1"
      ),
      leg_right = measurement_frame(
        list(node_y = right_leg), domain = right_domain,
        id = "laws:right-legs:v1"
      )
    )
    built
  }
})

# Q_H = B_L' H B_R, read through identity measurement legs so the package
# value is the native neural evidence form with nothing else applied.
pairing_law_native_form <- function(fixture, h = fixture$h, route = "auto") {
  form <- measurement_form(
    left = fixture$left_relation,
    between = edge_frame(
      "native_left", "native_right",
      fixture$identity_left, fixture$identity_right
    ),
    by = pair_query(h, fixture$left_space, fixture$right_space),
    over = fixture$over,
    right = fixture$right_relation,
    route = route
  )
  effect_coupling(form)$values[[1L]]
}

# The same evidence task the public constructor compiles, exposed so the
# contraction diagnostics can be read directly.
pairing_law_evidence_task <- function(fixture) {
  between <- edge_frame(
    "node_x", "node_y", fixture$leg_left, fixture$leg_right
  )
  crossform:::.new_evidence_task(
    fixture$left_relation, fixture$right_relation, FALSE,
    crossform:::.ordered_partition_edges(
      fixture$over, fixture$left_relation$partitions,
      fixture$right_relation$partitions, FALSE
    ),
    crossform:::.closed_experimental_boundary(
      pair_query(fixture$h, fixture$left_space, fixture$right_space)
    ),
    crossform:::.open_neural_boundary(
      between$from_frame, between$to_frame, between$edges
    ),
    crossform:::.evidence_stage_plan(),
    crossform:::.evidence_materialization(
      "measurement_form", "complete_form"
    )
  )
}

# The oracle's typed vocabulary, built once from the same fixture.
pairing_law_typed <- function(fixture) {
  experimental_left <- pair_oracle_space(
    fixture$left_space$basis_id, fixture$left_space$coordinates, "experimental"
  )
  experimental_right <- pair_oracle_space(
    fixture$right_space$basis_id, fixture$right_space$coordinates,
    "experimental"
  )
  neural_left <- pair_oracle_space(
    fixture$left_domain$id, paste0("nl", seq_len(fixture$p_left)), "neural"
  )
  neural_right <- pair_oracle_space(
    fixture$right_domain$id, paste0("nr", seq_len(fixture$p_right)), "neural"
  )
  list(
    experimental_left = experimental_left,
    experimental_right = experimental_right,
    neural_left = neural_left,
    neural_right = neural_right,
    left_relation = pair_oracle_relation(
      fixture$left_values, experimental_left, neural_left
    ),
    right_relation = pair_oracle_relation(
      fixture$right_values, experimental_right, neural_right
    )
  )
}


# Section 2, 3: the irreducible observable and its adjoint transports -------

test_that("the adjoint neural form equals the oracle transport", {
  fixture <- pairing_law_fixture()
  got <- pairing_law_native_form(fixture)

  expect_identical(dim(got), c(fixture$p_left, fixture$p_right))
  expect_equal(
    got,
    pair_oracle_adjoint(
      fixture$left_values, fixture$right_values, fixture$h
    ),
    tolerance = pairing_law_tolerance, ignore_attr = TRUE
  )
})

test_that("the forward transport equals the oracle on a weighted frame", {
  fixture <- pairing_law_fixture()
  # One relation pair over one shared neural domain: the rectangular geometry
  # plan closes the neural boundary with K = diag(frame weights).
  weights <- c(0.5, 1.25, 0.75, 2, 0.25)
  domain <- fixture$left_domain
  right_values <- fixture$shared_values
  plan <- plan_geometry(
    fixture$left_relation,
    additive_frame(matrix(weights, 1L, fixture$p_left), domain = domain),
    fixture$over,
    right = fixture$shared_relation
  )
  form <- materialize_geometry(plan)
  got <- matrix(
    geometry_component(form, "total")[1L, ], fixture$q_left, fixture$q_right
  )
  expected <- pair_oracle_forward(
    fixture$left_values, right_values, diag(weights)
  )

  expect_equal(got, expected, tolerance = pairing_law_tolerance,
    ignore_attr = TRUE)

  # Section 4: the packed rectangular codec is exactly column-major vec, and
  # the conceptual lift reproduces it on a bounded diagnostic materialization.
  lift <- pair_oracle_second_order_lift(
    fixture$left_values, right_values, materialize = TRUE
  )
  expect_identical(
    dim(lift),
    c(fixture$q_left * fixture$q_right, fixture$p_left * fixture$p_left)
  )
  expect_equal(
    as.vector(geometry_component(form, "total")[1L, ]),
    as.vector(lift %*% as.vector(diag(weights))),
    tolerance = pairing_law_tolerance
  )
})

test_that("forward, adjoint, and scalar evidence agree with the oracle", {
  set.seed(20260816L + 1L)
  fixture <- pairing_law_fixture()
  neural_query <- matrix(
    rnorm(fixture$p_left * fixture$p_right), fixture$p_left, fixture$p_right
  )
  native <- pairing_law_native_form(fixture)
  # <P*(H), K>_F read from the package's own adjoint materialization.
  got <- sum(native * neural_query)
  forward_side <- pair_oracle_evidence(
    fixture$left_values, fixture$right_values, fixture$h, neural_query
  )
  adjoint_side <- pair_oracle_evidence_adjoint(
    fixture$left_values, fixture$right_values, fixture$h, neural_query
  )

  expect_equal(forward_side, adjoint_side, tolerance = pairing_law_tolerance)
  expect_equal(got, forward_side, tolerance = pairing_law_tolerance)

  # The typed oracle carries the same number once every boundary identity is
  # checked rather than assumed from shape.
  typed <- pairing_law_typed(fixture)
  expect_equal(
    pair_oracle_typed_evidence(
      typed$left_relation, typed$right_relation,
      pair_oracle_bound_query(
        fixture$h, typed$experimental_left, typed$experimental_right, "effect"
      ),
      pair_oracle_bound_query(
        neural_query, typed$neural_left, typed$neural_right, "effect"
      )
    ),
    got, tolerance = pairing_law_tolerance
  )
})

test_that("the scalar pairing is bilinear in both arguments", {
  set.seed(20260816L + 2L)
  fixture <- pairing_law_fixture()
  second_h <- matrix(
    rnorm(fixture$q_left * fixture$q_right), fixture$q_left, fixture$q_right
  )
  neural_query <- matrix(
    rnorm(fixture$p_left * fixture$p_right), fixture$p_left, fixture$p_right
  )
  combined <- pairing_law_native_form(fixture, 2 * fixture$h - 3 * second_h)
  separate <- 2 * pairing_law_native_form(fixture) -
    3 * pairing_law_native_form(fixture, second_h)

  expect_equal(combined, separate, tolerance = pairing_law_tolerance)
  expect_equal(
    sum(combined * neural_query),
    2 * pair_oracle_evidence(
      fixture$left_values, fixture$right_values, fixture$h, neural_query
    ) - 3 * pair_oracle_evidence(
      fixture$left_values, fixture$right_values, second_h, neural_query
    ),
    tolerance = pairing_law_tolerance
  )
})


# Section 4: the lift is conceptual by default ------------------------------

test_that("the Kronecker lift is refused by default and bounded when asked", {
  fixture <- pairing_law_fixture()

  expect_error(
    pair_oracle_second_order_lift(fixture$left_values, fixture$right_values),
    "conceptual by default"
  )
  expect_error(
    pair_oracle_second_order_lift(
      fixture$left_values, fixture$right_values,
      materialize = TRUE, max_entries = 4L
    ),
    "bounded size guard"
  )

  # The production contraction never forms the lift on any admitted route,
  # and the lift is not even a candidate the route planner can price.
  task <- pairing_law_evidence_task(fixture)
  expected <- pair_oracle_measurement_form(
    fixture$left_values, fixture$right_values, fixture$h,
    fixture$left_leg, fixture$right_leg
  )
  for (route in c("forward_k", "pull_h", "multivariate_blocks")) {
    run <- crossform:::.run_measurement_contraction(task, route = route)
    expect_false(run$diagnostics$kronecker_lift_materialized, info = route)
    expect_false(
      run$diagnostics$global_neural_operator_materialized, info = route
    )
    expect_identical(run$route_plan$route, route)
    expect_equal(run$blocks[[1L]], expected,
      tolerance = pairing_law_tolerance, ignore_attr = TRUE, info = route)
  }
  expect_identical(
    names(crossform:::.measurement_route_plan(task)$costs),
    c("forward_k", "pull_h", "factorized_h", "scalar_stack",
      "multivariate_blocks")
  )
})


# Section 6: identified measurement legs and the four-leg law ---------------

test_that("the measured cross form equals the oracle leg contraction", {
  fixture <- pairing_law_fixture()
  form <- measurement_form(
    left = fixture$left_relation,
    between = edge_frame(
      "node_x", "node_y", fixture$leg_left, fixture$leg_right
    ),
    by = pair_query(fixture$h, fixture$left_space, fixture$right_space),
    over = fixture$over,
    right = fixture$right_relation
  )
  got <- effect_coupling(form)$values[[1L]]

  expect_identical(dim(got), c(2L, 3L))
  expect_equal(
    got,
    pair_oracle_measurement_form(
      fixture$left_values, fixture$right_values, fixture$h,
      fixture$left_leg, fixture$right_leg
    ),
    tolerance = pairing_law_tolerance, ignore_attr = TRUE
  )

  # The complete four-leg law: pushing a measured query forward and pulling
  # the experimental query backward are the same number.
  set.seed(20260816L + 3L)
  measured_query <- matrix(rnorm(2L * 3L), 2L, 3L)
  expect_equal(
    sum(fixture$h * pair_oracle_measured_effect_form(
      fixture$left_values, fixture$right_values,
      fixture$left_leg, fixture$right_leg, measured_query
    )),
    sum(measured_query * got),
    tolerance = pairing_law_tolerance
  )
})

test_that("two legs with the same metric give different oriented blocks", {
  fixture <- pairing_law_fixture()
  # A fixed 36-degree rotation: far enough from the identity that the
  # separation assertion below cannot depend on the RNG position.
  angle <- pi / 5
  rotation <- matrix(
    c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2L, 2L
  )
  rotated_leg <- rotation %*% fixture$left_leg
  rotated_frame <- measurement_frame(
    list(node_x = rotated_leg), domain = fixture$left_domain,
    id = "laws:left-legs-rotated:v1"
  )
  read <- function(frame) {
    effect_coupling(measurement_form(
      left = fixture$left_relation,
      between = edge_frame("node_x", "node_y", frame, fixture$leg_right),
      by = pair_query(fixture$h, fixture$left_space, fixture$right_space),
      over = fixture$over,
      right = fixture$right_relation
    ))$values[[1L]]
  }
  plain <- read(fixture$leg_left)
  rotated <- read(rotated_frame)

  # Same induced metric L'L, different oriented output coordinates.
  expect_equal(
    crossprod(rotated_leg), crossprod(fixture$left_leg),
    tolerance = pairing_law_tolerance
  )
  expect_false(isTRUE(all.equal(rotated, plain, tolerance = 1e-6)))
  expect_equal(rotated, rotation %*% plain,
    tolerance = pairing_law_tolerance)
  # Section 8: the Frobenius strength is the rotation-invariant summary.
  expect_equal(sum(rotated^2), sum(plain^2),
    tolerance = pairing_law_tolerance)
})


# Section 7: reversal ------------------------------------------------------

test_that("side reversal transposes every pairing object", {
  fixture <- pairing_law_fixture()
  forward <- effect_coupling(measurement_form(
    left = fixture$left_relation,
    between = edge_frame(
      "node_x", "node_y", fixture$leg_left, fixture$leg_right
    ),
    by = pair_query(fixture$h, fixture$left_space, fixture$right_space),
    over = fixture$over,
    right = fixture$right_relation
  ))$values[[1L]]
  reversed <- effect_coupling(measurement_form(
    left = fixture$right_relation,
    between = edge_frame(
      "node_y", "node_x", fixture$leg_right, fixture$leg_left
    ),
    by = pair_query(t(fixture$h), fixture$right_space, fixture$left_space),
    over = pairing("right_run", "left_run", 1,
      directed = TRUE, independence = "independent"),
    right = fixture$left_relation
  ))$values[[1L]]

  expect_equal(reversed, t(forward), tolerance = pairing_law_tolerance)
  expect_equal(
    reversed,
    t(pair_oracle_measurement_form(
      fixture$left_values, fixture$right_values, fixture$h,
      fixture$left_leg, fixture$right_leg
    )),
    tolerance = pairing_law_tolerance, ignore_attr = TRUE
  )

  set.seed(20260816L + 4L)
  neural_query <- matrix(
    rnorm(fixture$p_left * fixture$p_right), fixture$p_left, fixture$p_right
  )
  reversed_native <- effect_coupling(measurement_form(
    left = fixture$right_relation,
    between = edge_frame(
      "native_right", "native_left",
      fixture$identity_right, fixture$identity_left
    ),
    by = pair_query(t(fixture$h), fixture$right_space, fixture$left_space),
    over = pairing("right_run", "left_run", 1,
      directed = TRUE, independence = "independent"),
    right = fixture$left_relation
  ))$values[[1L]]
  # E_RL(H', K') = E_LR(H, K): both scalars are read from a package
  # materialization, and both are pinned to the oracle.
  expect_equal(
    sum(reversed_native * t(neural_query)),
    pair_oracle_evidence(
      fixture$left_values, fixture$right_values, fixture$h, neural_query
    ),
    tolerance = pairing_law_tolerance
  )
  expect_equal(
    sum(pairing_law_native_form(fixture) * neural_query),
    pair_oracle_evidence(
      fixture$right_values, fixture$left_values,
      t(fixture$h), t(neural_query)
    ),
    tolerance = pairing_law_tolerance
  )
  expect_equal(
    pairing_law_native_form(fixture),
    t(effect_coupling(measurement_form(
      left = fixture$right_relation,
      between = edge_frame(
        "native_right", "native_left",
        fixture$identity_right, fixture$identity_left
      ),
      by = pair_query(t(fixture$h), fixture$right_space, fixture$left_space),
      over = pairing("right_run", "left_run", 1,
        directed = TRUE, independence = "independent"),
      right = fixture$left_relation
    ))$values[[1L]]),
    tolerance = pairing_law_tolerance
  )
})


# Section 8: direct sums, decomposition lifting, low-rank contraction -------

test_that("forward transport is additive over a neural direct sum", {
  fixture <- pairing_law_fixture()
  domain <- fixture$left_domain
  right_values <- fixture$shared_values
  weights <- c(0.5, 1.25, 0.75, 2, 0.25)
  first <- c(weights[1:2], 0, 0, 0)
  second <- c(0, 0, weights[3:5])
  frame <- additive_frame(
    rbind(first, second, weights), domain = domain
  )
  form <- materialize_geometry(plan_geometry(
    fixture$left_relation, frame, fixture$over,
    right = fixture$shared_relation
  ))
  packed <- geometry_component(form, "total")
  unpack <- function(row) {
    matrix(packed[row, ], fixture$q_left, fixture$q_right)
  }

  expect_equal(unpack(3L), unpack(1L) + unpack(2L),
    tolerance = pairing_law_tolerance)
  expect_equal(
    unpack(1L),
    pair_oracle_forward(fixture$left_values, right_values, diag(first)),
    tolerance = pairing_law_tolerance, ignore_attr = TRUE
  )
  expect_equal(
    unpack(2L),
    pair_oracle_forward(fixture$left_values, right_values, diag(second)),
    tolerance = pairing_law_tolerance, ignore_attr = TRUE
  )
})

test_that("a low-rank query contracts without forming either full form", {
  set.seed(20260816L + 5L)
  fixture <- pairing_law_fixture()
  u <- matrix(rnorm(fixture$q_left * 2L), fixture$q_left, 2L)
  v <- matrix(rnorm(fixture$q_right * 2L), fixture$q_right, 2L)
  r <- matrix(rnorm(fixture$p_left * 2L), fixture$p_left, 2L)
  s <- matrix(rnorm(fixture$p_right * 2L), fixture$p_right, 2L)
  h <- tcrossprod(u, v)
  factorized <- pairing_law_native_form(fixture, h, route = "factorized_h")
  dense <- pairing_law_native_form(fixture, h, route = "pull_h")

  expect_equal(factorized, dense, tolerance = pairing_law_tolerance)
  expect_equal(
    sum(factorized * tcrossprod(r, s)),
    pair_oracle_factorized_evidence(
      fixture$left_values, fixture$right_values, u, v, r, s
    ),
    tolerance = pairing_law_tolerance
  )
  expect_equal(
    pair_oracle_factorized_evidence(
      fixture$left_values, fixture$right_values, u, v, r, s
    ),
    pair_oracle_evidence(
      fixture$left_values, fixture$right_values, h, tcrossprod(r, s)
    ),
    tolerance = pairing_law_tolerance
  )
})


# Section 1, 5: identity is not shape --------------------------------------

test_that("same-shaped but differently identified boundaries are rejected", {
  fixture <- pairing_law_fixture()
  impostor_space <- effect_space(
    fixture$left_space$coordinates, basis_id = "laws:impostor-effects:v1"
  )
  impostor_domain <- abstract_domain(
    fixture$p_left, id = "laws:impostor-neural:v1"
  )
  impostor_frame <- measurement_frame(
    list(native_left = diag(fixture$p_left)), domain = impostor_domain,
    id = "laws:impostor-identity:v1"
  )
  read <- function(by, between) {
    measurement_form(
      left = fixture$left_relation, between = between, by = by,
      over = fixture$over, right = fixture$right_relation
    )
  }
  native_edges <- edge_frame(
    "native_left", "native_right",
    fixture$identity_left, fixture$identity_right
  )

  expect_error(
    read(
      pair_query(fixture$h, impostor_space, fixture$right_space),
      native_edges
    ),
    "does not match the left and right experimental spaces"
  )
  expect_error(
    read(
      pair_query(fixture$h, fixture$left_space, fixture$right_space),
      edge_frame(
        "native_left", "native_right", impostor_frame, fixture$identity_right
      )
    ),
    "do not belong to their experimental relation sides"
  )

  # The oracle refuses the identical mismatch without looking at any value.
  typed <- pairing_law_typed(fixture)
  impostor <- pair_oracle_space(
    "laws:impostor-effects:v1", fixture$left_space$coordinates, "experimental"
  )
  expect_error(
    pair_oracle_assert_boundaries(
      typed$left_relation, typed$right_relation,
      pair_oracle_bound_query(
        fixture$h, impostor, typed$experimental_right, "effect"
      ),
      pair_oracle_bound_query(
        matrix(0, fixture$p_left, fixture$p_right),
        typed$neural_left, typed$neural_right, "effect"
      )
    ),
    "experimental query axis identities differ"
  )
  expect_error(
    pair_oracle_relation(
      fixture$left_values, typed$experimental_left, typed$neural_right
    ),
    "do not match its identified boundaries"
  )
})


# Section 9: exact effect-form specialization ------------------------------

test_that("a factorized bridge supplies K = L_left' L_right exactly", {
  fixture <- pairing_law_fixture()
  bridge <- measurement_bridge(
    fixture$left_leg,
    fixture$right_leg[seq_len(2L), , drop = FALSE],
    fixture$left_domain, fixture$right_domain,
    measurement_space(2L, id = "laws:common-modes:v1")
  )
  bridged <- crossform:::.apply_measurement_bridge(
    list(left_run = fixture$left_values),
    list(right_run = fixture$right_values), bridge
  )
  got <- bridged$left$left_run %*% t(bridged$right$right_run)
  neural_query <- t(bridge$left_leg) %*% bridge$right_leg

  expect_equal(
    got,
    pair_oracle_forward(
      fixture$left_values, fixture$right_values, neural_query
    ),
    tolerance = pairing_law_tolerance
  )
  expect_equal(
    sum(fixture$h * got),
    pair_oracle_evidence(
      fixture$left_values, fixture$right_values, fixture$h, neural_query
    ),
    tolerance = pairing_law_tolerance
  )

  # Reversal exchanges the legs and transposes the induced operator.
  reversed <- reverse_bridge(bridge)
  expect_equal(
    t(reversed$left_leg) %*% reversed$right_leg, t(neural_query),
    tolerance = pairing_law_tolerance
  )
  reversed_applied <- crossform:::.apply_measurement_bridge(
    list(right_run = fixture$right_values),
    list(left_run = fixture$left_values), reversed
  )
  expect_equal(
    reversed_applied$left$right_run %*% t(reversed_applied$right$left_run),
    t(pair_oracle_forward(
      fixture$left_values, fixture$right_values, neural_query
    )),
    tolerance = pairing_law_tolerance
  )
})


# Sections 10-13: capabilities, degeneracy, and normalized views -----------

# A self-form coupling fixture: two oriented scalar node measurements over one
# neural domain, closed by a declared repeated-variation query.
coupling_law_fixture <- local({
  cache <- list()
  function(rank_one = FALSE, partition = NULL, multivariate = FALSE) {
    key <- paste0(rank_one, "/", multivariate, "/",
      if (is.null(partition)) "both" else partition)
    if (!is.null(cache[[key]])) {
      return(cache[[key]])
    }
    set.seed(20260817L)
    q <- 8L
    p <- 4L
    domain <- abstract_domain(p, id = "laws:coupling-neural:v1")
    effects <- effect_space(paste0("trial", seq_len(q)),
      basis_id = "laws:coupling-trials:v1")
    run1 <- matrix(rnorm(q * p), q, p)
    run2 <- matrix(rnorm(q * p), q, p)
    run1[, 2L] <- -0.7 * run1[, 1L] + 0.5 * run1[, 2L]
    run2[, 2L] <- -0.6 * run2[, 1L] + 0.55 * run2[, 2L]
    values <- list(run1 = run1, run2 = run2)
    if (!is.null(partition)) values <- values[partition]
    rel <- relation(values, effects = effects, domain = domain)
    h <- if (rank_one) {
      tcrossprod(rep(c(-1, 1), each = q / 2L))
    } else {
      (diag(q) - matrix(1 / q, q, q)) / (q - 1L)
    }
    operators <- if (multivariate) {
      list(a = diag(p)[1:2, , drop = FALSE], b = diag(p)[3:4, , drop = FALSE])
    } else {
      list(a = matrix(c(1, 0, 0, 0), 1L), b = matrix(c(0, 1, 0, 0), 1L))
    }
    frame <- measurement_frame(
      operators, domain = domain,
      id = paste0("laws:coupling-frame:", if (multivariate) "mv" else "sc")
    )
    pairs <- expand.grid(
      from = c("a", "b"), to = c("a", "b"), stringsAsFactors = FALSE
    )
    form <- measurement_form(
      left = rel,
      between = edge_frame(pairs$from, pairs$to, frame),
      by = variation_query(
        h, effects, "trial", "joint_covariance",
        provenance = list(estimator = "centered-within-run")
      ),
      over = pairing(
        rel$partitions, rel$partitions, directed = TRUE,
        self_pairs = "allow_biased", independence = "not_independent"
      ),
      route = "pull_h"
    )
    value <- list(form = form, h = h, values = values, q = q)
    cache[[key]] <<- value
    value
  }
})

# Edge ids follow the row order of `expand.grid`: (a,a), (b,a), (a,b), (b,b).
coupling_law_blocks <- function(fixture) {
  blocks <- effect_coupling(fixture$form)$values
  list(
    cross = blocks[["edge_000003"]],
    left_self = blocks[["edge_000001"]],
    right_self = blocks[["edge_000004"]],
    edge_id = "edge_000003"
  )
}

test_that("construction capabilities are separate from observed diagnostics", {
  fixture <- coupling_law_fixture()
  capabilities <- fixture$form$capabilities
  diagnostics <- fixture$form$diagnostics
  oracle <- pair_oracle_query_capabilities(
    role = "variation", construction = "covariance",
    sampling_axis = "trial",
    effective_rank = diagnostics$experimental_effective_rank,
    observed_min_eigenvalue = min(eigen(fixture$h, symmetric = TRUE,
      only.values = TRUE)$values)
  )
  observed_eigenvalues <- fixture$form$diagnostics$blocks

  expect_true(capabilities$repeated_variation)
  expect_true(capabilities$joint_covariance)
  expect_identical(capabilities$sampling_axis,
    oracle$guarantees$sampling_axis)
  expect_true(oracle$guarantees$positive_query)
  expect_true(oracle$guarantees$joint_covariance)
  # The rank is a diagnostic of a finite realization, not a guarantee: the
  # centering operator on q trials has exactly q - 1 nonzero directions.
  expect_identical(oracle$diagnostics$effective_rank, fixture$q - 1L)
  expect_identical(oracle$provenance$construction, "covariance")
  # The smallest eigenvalue of the centering operator is zero to working
  # precision, and it is recorded as an observation, never as a guarantee.
  expect_equal(oracle$diagnostics$observed_min_eigenvalue, 0,
    tolerance = 1e-12)
  expect_false(isTRUE(oracle$guarantees$guaranteed_psd))
  expect_true(all(
    observed_eigenvalues$effective_rank <= fixture$q - 1L
  ))

  # An unlabelled effect query carries no sampling capability at all, even
  # though its matrix is the same size and also happens to be symmetric.
  plain <- pair_oracle_query_capabilities(role = "effect")
  expect_identical(plain$guarantees$role, "effect")
  expect_false(plain$guarantees$positive_query)
  expect_null(plain$guarantees$sampling_axis)
  expect_error(
    pair_oracle_query_capabilities(role = "variation"),
    "requires declared sampling provenance"
  )
  expect_error(
    pair_oracle_query_capabilities(role = "effect", construction = "covariance"),
    "require the variation role"
  )
})

test_that("the oracle view gate mirrors the package's admitted views", {
  fixture <- coupling_law_fixture()
  diagnostics <- fixture$form$diagnostics
  capable <- pair_oracle_query_capabilities(
    role = "variation", construction = "covariance", sampling_axis = "trial",
    effective_rank = diagnostics$experimental_effective_rank
  )
  bare <- pair_oracle_query_capabilities(role = "effect")

  expect_silent(pair_oracle_require_view(bare, "effect_coupling"))
  expect_silent(pair_oracle_require_view(capable, "covariance_coupling"))
  expect_silent(pair_oracle_require_view(
    capable, "canonical_coupling",
    positive_self_blocks = TRUE, regularization_policy = "ridge:0.05"
  ))
  expect_error(
    pair_oracle_require_view(bare, "covariance_coupling"),
    "requires valid repeated variation"
  )
  expect_error(
    pair_oracle_require_view(capable, "canonical_coupling"),
    "requires joint covariance and regularized self-blocks"
  )

  expect_s3_class(effect_coupling(fixture$form), "effect_coupling_result")
  expect_s3_class(covariance_coupling(fixture$form), "effect_coupling_result")
  expect_s3_class(
    canonical_coupling(fixture$form, ridge = 0.05), "effect_coupling_result"
  )
  expect_error(
    canonical_coupling(fixture$form),
    "argument \"ridge\" is missing"
  )
  refusal <- catch_refusal(
    connectivity(fixture$form, view = "canonical")
  )
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "declared_regularization")
})

test_that("normalized coupling views equal their oracles", {
  fixture <- coupling_law_fixture()
  blocks <- coupling_law_blocks(fixture)
  ridge <- 0.05

  expect_equal(
    connectivity(fixture$form, "geometry_alignment")$values$geometry_alignment[
      match(blocks$edge_id,
        connectivity(fixture$form, "geometry_alignment")$values$edge_id)
    ],
    pair_oracle_geometry_alignment(
      blocks$cross, blocks$left_self, blocks$right_self
    ),
    tolerance = pairing_law_tolerance
  )

  canonical <- connectivity(fixture$form, "canonical", ridge = ridge)$values
  expect_equal(
    canonical$canonical_correlation[canonical$edge_id == blocks$edge_id],
    pair_oracle_canonical_spectrum(
      blocks$cross,
      blocks$left_self + diag(ridge, nrow(blocks$left_self)),
      blocks$right_self + diag(ridge, nrow(blocks$right_self))
    ),
    tolerance = pairing_law_tolerance
  )

  correlation <- connectivity(fixture$form, "correlation")$values
  expect_equal(
    correlation$correlation[correlation$edge_id == blocks$edge_id],
    blocks$cross[[1L]] /
      sqrt(blocks$left_self[[1L]] * blocks$right_self[[1L]]),
    tolerance = pairing_law_tolerance
  )

  # The whitening step itself matches the oracle inverse square root.
  self_block <- matrix(c(2, 0.3, 0.3, 1.5), 2L, 2L)
  expect_equal(
    crossform:::.decomposition_inverse_sqrt(self_block, 1e-12),
    pair_oracle_inverse_sqrt(self_block),
    tolerance = pairing_law_tolerance
  )
  expect_error(
    pair_oracle_inverse_sqrt(matrix(c(1, 2, 3, 4), 2L)),
    "must be symmetric"
  )
  expect_error(
    pair_oracle_inverse_sqrt(diag(c(1, -1))),
    "not positive semidefinite"
  )
})

test_that("rank-one coupling degeneracy is algebraic and gains no capability", {
  # One partition and multivariate legs, so that H = c c' really does give
  # Q_H = u u' and the measured blocks are exactly a b', a a', b b'.
  fixture <- coupling_law_fixture(
    rank_one = TRUE, partition = "run1", multivariate = TRUE
  )
  blocks <- coupling_law_blocks(fixture)
  contrast <- rep(c(-1, 1), each = fixture$q / 2L)
  u <- drop(crossprod(fixture$values$run1, contrast))
  a <- (diag(4L)[1:2, , drop = FALSE]) %*% u
  b <- (diag(4L)[3:4, , drop = FALSE]) %*% u

  expect_equal(blocks$cross, tcrossprod(a, b),
    tolerance = pairing_law_tolerance)
  expect_equal(blocks$left_self, tcrossprod(a),
    tolerance = pairing_law_tolerance)
  expect_equal(blocks$right_self, tcrossprod(b),
    tolerance = pairing_law_tolerance)

  # Section 11: both perfect values are forced by the rank-one construction.
  expect_equal(
    pair_oracle_geometry_alignment(
      blocks$cross, blocks$left_self, blocks$right_self
    ),
    1, tolerance = pairing_law_tolerance
  )
  # The unregularized whitening inverts a rank-one block through an
  # eigenvalue threshold, so this one comparison carries the pseudo-inverse's
  # own conditioning rather than the file's 1e-10 floor.
  expect_equal(
    max(pair_oracle_canonical_spectrum(
      blocks$cross, blocks$left_self, blocks$right_self
    )),
    1, tolerance = 1e-8
  )

  # The block is still a valid effect coupling, and it still equals the
  # oracle contraction; only the normalized interpretation is refused.
  expect_s3_class(effect_coupling(fixture$form), "effect_coupling_result")
  for (view in c("correlation", "canonical", "geometry_alignment")) {
    refusal <- catch_refusal(
      connectivity(fixture$form, view = view, ridge = 0.05)
    )
    expect_s3_class(refusal, "effect_capability_refusal")
    expect_identical(refusal$capability, "nondegenerate_variation")
    expect_identical(refusal$namespace, "coupling_views")
  }
  expect_identical(fixture$form$diagnostics$experimental_effective_rank, 1L)
  expect_error(
    pair_oracle_require_view(
      pair_oracle_query_capabilities(
        role = "variation", construction = "covariance",
        sampling_axis = "trial", effective_rank = 1L
      ),
      "covariance_coupling"
    ),
    "rank above one"
  )
})

test_that("normalization placement changes the estimand, as the oracle says", {
  first <- coupling_law_fixture(partition = "run1")
  second <- coupling_law_fixture(partition = "run2")
  forms <- list(first$form, second$form)
  weights <- c(0.3, 0.7)
  blocks <- lapply(forms, function(form) {
    values <- effect_coupling(form)$values
    list(
      cross = values[["edge_000003"]][[1L]],
      left_self = values[["edge_000001"]][[1L]],
      right_self = values[["edge_000004"]][[1L]]
    )
  })
  cross <- vapply(blocks, `[[`, numeric(1), "cross")
  left_self <- vapply(blocks, `[[`, numeric(1), "left_self")
  right_self <- vapply(blocks, `[[`, numeric(1), "right_self")

  edge_first <- crossform:::.partitioned_pearson_coupling(
    forms, weights, "within_partition_pair"
  )
  aggregate_first_result <- crossform:::.partitioned_pearson_coupling(
    forms, weights, "after_partition_aggregation"
  )
  pick <- function(result) {
    result$values$value[result$values$edge_id == "edge_000003"]
  }

  expect_equal(
    pick(edge_first),
    pair_oracle_scalar_normalization(
      cross, left_self, right_self, weights, "within_then_reduce"
    ),
    tolerance = pairing_law_tolerance
  )
  expect_equal(
    pick(aggregate_first_result),
    pair_oracle_scalar_normalization(
      cross, left_self, right_self, weights, "reduce_then_normalize"
    ),
    tolerance = pairing_law_tolerance
  )
  # Section 12: normalize(sum) != sum(normalize) in general, so the two
  # placements are different estimands with different identities.
  expect_false(isTRUE(all.equal(
    pick(edge_first), pick(aggregate_first_result), tolerance = 1e-6
  )))
  expect_false(identical(
    edge_first$signature, aggregate_first_result$signature
  ))
  expect_error(
    pair_oracle_scalar_normalization(
      cross, c(left_self[[1L]], 0), right_self, weights
    ),
    "matching positive self moments"
  )
})

test_that("the stage ledger's recorded order and placements agree", {
  edge_first <- crossform:::.evidence_stage_plan(reducer = reduce_partitions())
  aggregated <- crossform:::.evidence_stage_plan(reducer = aggregate_first())

  # The package names the reduction stage differently in the two orders, and
  # records an explicit placement only on the normalization stage. Everything
  # the oracle needs is therefore derived from the recorded `order`: a stage
  # that appears before the reduction happens within a partition pair, and a
  # stage at or after it happens on the aggregate.
  reduction_name <- function(stages) {
    intersect(c("partition_reduction", "partition_sufficient_reduction"),
      stages$order)
  }
  placement_from_order <- function(stages, stage_name) {
    reduction <- match(reduction_name(stages), stages$order)
    if (match(stage_name, stages$order) < reduction) {
      "within_pair"
    } else {
      "after_aggregation"
    }
  }
  as_oracle_stages <- function(stages) {
    list(
      pair_oracle_stage(
        stages$partition_product$operation, stages$partition_product$axis,
        placement_from_order(stages, "ordered_partition_product")
      ),
      pair_oracle_stage(
        stages$normalization$operation$kind, stages$normalization$axis,
        placement_from_order(stages, "neural_feature_normalization")
      ),
      pair_oracle_stage(
        stages$transform$operation$kind, stages$transform$axis,
        placement_from_order(stages, "edge_transform")
      ),
      pair_oracle_stage(
        stages$reduction$operation$kind, stages$reduction$axis,
        placement_from_order(stages, reduction_name(stages))
      )
    )
  }

  expect_identical(
    pair_oracle_stage_identity(as_oracle_stages(edge_first)),
    c(
      "ordered_partition_product@partition_pairs@within_pair",
      "inner_product@neural_features@within_pair",
      "identity@form_entries@within_pair",
      "weighted_sum@partition_pairs@after_aggregation"
    )
  )
  expect_identical(
    pair_oracle_stage_identity(as_oracle_stages(aggregated)),
    c(
      "ordered_partition_product@partition_pairs@within_pair",
      "inner_product@neural_features@after_aggregation",
      "identity@form_entries@after_aggregation",
      "weighted_sum@partition_pairs@after_aggregation"
    )
  )

  # The one placement the package does record must agree with the order it
  # also records; a plan whose two statements disagreed would fail here.
  vocabulary <- c(within_partition_pair = "within_pair",
    after_partition_aggregation = "after_aggregation")
  for (stages in list(edge_first, aggregated)) {
    expect_identical(
      unname(vocabulary[[stages$normalization$placement]]),
      placement_from_order(stages, "neural_feature_normalization")
    )
  }
  # These are the fields the ledger actually carries; the oracle's demand for
  # a placement on every stage is met by the recorded order, not by a stored
  # `placement` field on the product, reduction, and projection stages.
  expect_null(edge_first$partition_product$placement)
  expect_null(edge_first$reduction$placement)
  expect_null(edge_first$projection$placement)
  expect_identical(edge_first$transform$placement, "before_partition_reduction")
  expect_identical(edge_first$projection$axis, "experimental_coordinates")

  expect_false(identical(edge_first$signature, aggregated$signature))
  expect_false(identical(edge_first$order, aggregated$order))
  expect_lt(
    match("neural_feature_normalization", edge_first$order),
    match("partition_reduction", edge_first$order)
  )
  expect_gt(
    match("neural_feature_normalization", aggregated$order),
    match("partition_sufficient_reduction", aggregated$order)
  )

  # A stage that omits its placement is not an admissible ledger entry.
  expect_error(
    pair_oracle_stage_identity(list(
      list(operation = "ordered_partition_product", axis = "partition_pairs")
    )),
    "must identify operation, axis, and placement"
  )
  expect_error(
    pair_oracle_stage("", "partition_pairs", "within_pair"),
    "requires one explicit name"
  )
  expect_error(
    pair_oracle_stage("ordered_partition_product", "voxels", "within_pair"),
    "'arg' should be one of"
  )
})
