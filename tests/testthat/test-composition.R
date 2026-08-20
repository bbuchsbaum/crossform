# The metric composition law (D6).
#
# `design/conservative-geometry-contract.md` section 5 states two compositions
# of a fixed metric with a spatial frame:
#
#   native    K_x = D(sqrt(w_x)) Q D(sqrt(w_x))
#   whitened  K_x = Q^(1/2) D(w_x) Q^(1/2)
#
# They agree for a diagonal `Q` and disagree for a dense one, where only the
# whitened composition closes the conservative ledger. This file is the
# executable half of that section: hand-computed dense algebra for both
# compositions, the conservation numbers the contract's section 9 tolerances
# name, the estimand separation in plan identity, and -- the gap G5 the
# 2026-08-20 review left open -- the pin on *which* root `"whitened"` means.

composition_fixture <- function(seed = 20260817L, features = 9L) {
  set.seed(seed)
  q <- 3L
  domain <- abstract_domain(
    features,
    coordinates = cbind(seq_len(features) - 1, 0),
    id = "composition-fixture"
  )
  block <- function() {
    value <- matrix(stats::rnorm(q * features), q, features)
    rownames(value) <- c("a", "b", "c")
    value
  }
  effects <- list(run1 = block(), run2 = block())
  relation <- relation(effects, domain = domain)
  raw <- matrix(stats::rnorm(features * features), features, features)
  list(
    domain = domain,
    relation = relation,
    over = cross_partitions(relation),
    left = effects$run1,
    right = effects$run2,
    q = q,
    conservative = compile_frame(
      searchlights(1.01, "conservative"), domain
    ),
    global = compile_frame(whole_brain("none"), domain),
    dense = crossprod(raw) / features + diag(features),
    diagonal = diag(stats::runif(features, 0.5, 2.5))
  )
}

# The oracle side: pack a symmetric matrix exactly as `symmetric_packed` does,
# and form the cross-partition estimator sym(B1 K B2^T) by hand.
composition_svec <- function(value) {
  value <- 0.5 * (value + t(value))
  q <- nrow(value)
  out <- numeric(q * (q + 1L) / 2L)
  slot <- 1L
  for (column in seq_len(q)) {
    for (row in column:q) {
      out[[slot]] <- if (row == column) {
        value[row, column]
      } else {
        sqrt(2) * value[row, column]
      }
      slot <- slot + 1L
    }
  }
  out
}

composition_oracle_nodes <- function(fixture, operator) {
  weights <- as.matrix(fixture$conservative$weights)
  left <- fixture$left
  right <- fixture$right
  t(vapply(seq_len(nrow(weights)), function(node) {
    kernel <- operator(weights[node, ])
    composition_svec(
      0.5 * (left %*% kernel %*% t(right) + right %*% kernel %*% t(left))
    )
  }, numeric(fixture$q * (fixture$q + 1L) / 2L)))
}

composition_total <- function(fixture, frame, metric = NULL,
                              composition = "native") {
  plan <- plan_geometry(
    fixture$relation, frame, fixture$over,
    metric = metric, composition = composition
  )
  list(
    plan = plan,
    total = geometry_component(materialize_geometry(plan), "total")
  )
}

# The symmetric positive-definite root, computed here from first principles so
# the test does not certify the implementation against itself.
composition_symmetric_root <- function(value) {
  decomposition <- eigen(value, symmetric = TRUE)
  root <- decomposition$vectors %*%
    (sqrt(decomposition$values) * t(decomposition$vectors))
  (root + t(root)) / 2
}

test_that("both compositions equal their hand-computed dense algebra", {
  fixture <- composition_fixture()
  metric <- neural_metric(fixture$dense, fixture$domain)
  root <- composition_symmetric_root(fixture$dense)

  native <- composition_total(fixture, fixture$conservative, metric, "native")
  whitened <- composition_total(
    fixture, fixture$conservative, metric, "whitened"
  )

  expected_native <- composition_oracle_nodes(fixture, function(weight) {
    fixture$dense * tcrossprod(sqrt(weight))
  })
  expected_whitened <- composition_oracle_nodes(fixture, function(weight) {
    root %*% diag(weight) %*% root
  })

  expect_equal(unname(native$total), expected_native, tolerance = 1e-12)
  expect_equal(unname(whitened$total), expected_whitened, tolerance = 1e-12)

  # The two are different estimands, not two routes to one number. The
  # contract measures the gap at roughly 13.5% of the largest node value on
  # this fixture; the invariant claim is only that it is not small.
  gap <- max(abs(expected_whitened - expected_native)) /
    max(abs(expected_native))
  expect_gt(gap, 0.01)
})

test_that("a dense metric conserves under whitening and not natively", {
  fixture <- composition_fixture()
  metric <- neural_metric(fixture$dense, fixture$domain)

  global <- drop(composition_total(fixture, fixture$global, metric)$total)
  native <- composition_total(fixture, fixture$conservative, metric)$total
  whitened <- composition_total(
    fixture, fixture$conservative, metric, "whitened"
  )$total

  # Contract section 9: the whitened law is asserted at 1e-12, the native
  # failure only as "not small", because its size and sign are fixture
  # specific (section 5.1, and oracle O2.c').
  expect_lt(max(abs(colSums(whitened) - global)), 1e-12)
  expect_gt(max(abs(colSums(native) - global)) / max(abs(global)), 0.01)
})

test_that("whitening means the symmetric root, not any root that conserves", {
  # Gap G5 of the contract's section 11.4, promoted from oracle O2.d'.
  #
  # `sum_x R D(w_x) R^T = R I R^T = Q` for ANY R with `R R^T = Q`, so the
  # conservation certificate cannot tell two roots apart. The per-node values
  # can, by a wide margin. `composition = "whitened"` therefore has to name a
  # root, and this test pins the one it names.
  fixture <- composition_fixture()
  metric <- neural_metric(fixture$dense, fixture$domain)
  symmetric <- composition_symmetric_root(fixture$dense)
  cholesky <- t(chol(fixture$dense))

  expect_equal(symmetric %*% t(symmetric), fixture$dense, tolerance = 1e-12)
  expect_equal(cholesky %*% t(cholesky), fixture$dense, tolerance = 1e-12)

  global <- drop(composition_total(fixture, fixture$global, metric)$total)
  by_symmetric <- composition_oracle_nodes(fixture, function(weight) {
    symmetric %*% diag(weight) %*% symmetric
  })
  by_cholesky <- composition_oracle_nodes(fixture, function(weight) {
    cholesky %*% diag(weight) %*% t(cholesky)
  })

  # Both roots conserve, to the same tolerance.
  expect_lt(max(abs(colSums(by_symmetric) - global)), 1e-12)
  expect_lt(max(abs(colSums(by_cholesky) - global)), 1e-12)

  # And their node values disagree materially. The contract measures 15.7% on
  # this fixture and 28.6% on another draw; the claim asserted here is that it
  # is not zero, never a particular percentage.
  root_gap <- max(abs(by_symmetric - by_cholesky)) / max(abs(by_symmetric))
  expect_gt(root_gap, 0.01)

  # So the package's whitened values must match the symmetric root exactly,
  # not merely conserve. A Cholesky implementation would pass every
  # conservation assertion above and fail this one by `root_gap`.
  executed <- composition_total(
    fixture, fixture$conservative, metric, "whitened"
  )$total
  expect_equal(unname(executed), by_symmetric, tolerance = 1e-12)
  expect_gt(max(abs(executed - by_cholesky)) / max(abs(by_symmetric)), 0.01)

  # The pinned convention is recorded, not implied.
  plan <- plan_geometry(
    fixture$relation, fixture$conservative, fixture$over,
    metric = metric, composition = "whitened"
  )
  expect_identical(plan$metric_schedule$root, "symmetric_psd_root")
  expect_identical(plan$metric_schedule$composition, "whitened")
})

test_that("composition is estimand-bearing in the scientific plan id", {
  fixture <- composition_fixture()
  metric <- neural_metric(fixture$dense, fixture$domain)

  native <- plan_geometry(
    fixture$relation, fixture$conservative, fixture$over, metric = metric
  )
  explicit_native <- plan_geometry(
    fixture$relation, fixture$conservative, fixture$over,
    metric = metric, composition = "native"
  )
  whitened <- plan_geometry(
    fixture$relation, fixture$conservative, fixture$over,
    metric = metric, composition = "whitened"
  )

  # The default is the behaviour that existed before the argument did, and
  # naming it explicitly must not move any identity.
  expect_identical(
    native$scientific_plan_id, explicit_native$scientific_plan_id
  )
  expect_identical(native$signature, explicit_native$signature)
  expect_null(native$metric_schedule$composition)

  # Two plans differing only in composition are two estimands.
  expect_false(identical(
    native$scientific_plan_id, whitened$scientific_plan_id
  ))
  expect_false(identical(
    native$metric_schedule$signature, whitened$metric_schedule$signature
  ))

  # The metric signature and the root convention are both inside the schedule
  # identity, so a different metric or a future different root would be a
  # different plan rather than a silent substitution.
  forged <- whitened$metric_schedule
  forged$root <- "lower_cholesky_root"
  expect_error(
    crossform:::.validate_geometry_metric_schedule(forged),
    class = "effect_contract_error"
  )
  expect_identical(
    whitened$metric_schedule$metric_signature, metric$signature
  )

  # A view derives its identity from its parent estimand, so the composition
  # reaches every view of the plan without any view knowing about it.
  weights <- c(a = 1, b = -1, c = 0)
  native_view <- contrast_energy(native, weights)$metadata$scientific_plan_id
  whitened_view <-
    contrast_energy(whitened, weights)$metadata$scientific_plan_id
  expect_true(nzchar(native_view))
  expect_false(identical(native_view, whitened_view))
})

test_that("a conservative contrast closes its ledger under whitening", {
  # The scientific payoff, stated as the law rather than as an internal: a
  # conservative crossnobis-style contrast under a *dense* metric conserves
  # under the whitened composition and does not under the native one.
  fixture <- composition_fixture()
  metric <- neural_metric(fixture$dense, fixture$domain)
  weights <- c(a = 1, b = -1, c = 0)

  energy <- function(frame, composition) {
    contrast_energy(
      plan_geometry(
        fixture$relation, frame, fixture$over,
        metric = metric, composition = composition
      ),
      weights
    )$total
  }

  global_native <- energy(fixture$global, "native")
  global_whitened <- energy(fixture$global, "whitened")
  # One global comparator: the whole-support contrast is the same number under
  # either composition, because a single unit-weight node has nothing to
  # whiten against.
  expect_equal(global_whitened, global_native, tolerance = 1e-12)

  expect_lt(
    abs(sum(energy(fixture$conservative, "whitened")) - global_whitened),
    1e-12
  )
  expect_gt(
    abs(sum(energy(fixture$conservative, "native")) - global_native) /
      abs(global_native),
    0.01
  )
})

test_that("the conservation certificate reports the law in force", {
  fixture <- composition_fixture()
  metric <- neural_metric(fixture$dense, fixture$domain)
  diagonal <- noise_precision(fixture$diagonal, fixture$domain)
  weights <- as.matrix(fixture$conservative$weights)

  node_metric <- function(builder) {
    lapply(seq_len(nrow(weights)), function(node) {
      positions <- which(weights[node, ] > 0)
      neural_metric(
        builder(positions), fixture$domain,
        support = fixture$domain$feature_ids[positions]
      )
    })
  }

  # Dense, native: refused, exactly as before D6. The refusal now names the
  # alternative without claiming it repairs anything.
  dense_native <- crossform:::.metric_frame_conservation(
    fixture$conservative,
    node_metric(function(positions) {
      fixture$dense[positions, positions, drop = FALSE]
    })
  )
  expect_identical(dense_native$composition, "native")
  expect_false(dense_native$feature_additive)
  expect_false(dense_native$identity_conservation)
  expect_identical(dense_native$global_metric_kind, "support_pair_operator")
  expect_match(dense_native$reason, "do not conserve a non-diagonal")
  expect_match(dense_native$reason, "composition = \\\"whitened\\\"")
  expect_error(
    crossform:::.require_metric_conservation(dense_native, "identity"),
    "not certified", class = "effect_input_error"
  )

  # Dense, whitened: certified conserved at the contract's 1e-12.
  dense_whitened <- crossform:::.metric_frame_conservation(
    fixture$conservative, composition = "whitened", metric = metric
  )
  expect_identical(dense_whitened$composition, "whitened")
  expect_identical(dense_whitened$root, "symmetric_psd_root")
  expect_true(dense_whitened$identity_conservation)
  expect_lt(dense_whitened$max_deviation, 1e-12)
  expect_identical(
    dense_whitened$global_metric_kind, "whitened_domain_operator"
  )
  expect_equal(
    dense_whitened$reference_mass, diag(fixture$dense), tolerance = 1e-12
  )
  expect_equal(
    dense_whitened$global_diagonal, diag(fixture$dense), tolerance = 1e-12
  )
  expect_silent(
    crossform:::.require_metric_conservation(dense_whitened, "identity")
  )
  # Conserving is not being feature additive, and the two targets say so
  # separately rather than one standing in for the other.
  expect_false(dense_whitened$feature_additive)
  expect_error(
    crossform:::.require_metric_conservation(dense_whitened,
      "feature_additive"),
    "Feature additivity is not certified, though conservation is",
    class = "effect_input_error"
  )

  # A local frame does not conserve under whitening either: the law needs unit
  # column mass, and whitening supplies no substitute for it.
  local_whitened <- crossform:::.metric_frame_conservation(
    compile_frame(searchlights(1.01), fixture$domain),
    composition = "whitened", metric = metric
  )
  expect_false(local_whitened$identity_conservation)
  expect_gt(local_whitened$max_deviation, 1e-8)
  expect_match(local_whitened$reason, "column mass is one")

  # Diagonal, native: feature additive, so the composition survives the frame.
  # The certificate on the *declared* frame compares against unit mass and so
  # reports that the summed diagonal is the metric, not the identity -- which
  # is the honest reading of that comparator, and why B5 gave the fold its own
  # reference mass.
  diagonal_native <- crossform:::.metric_frame_conservation(
    fixture$conservative,
    node_metric(function(positions) {
      fixture$diagonal[positions, positions, drop = FALSE]
    })
  )
  expect_identical(diagonal_native$composition, "native")
  expect_true(diagonal_native$feature_additive)
  expect_equal(
    diagonal_native$global_diagonal, diag(fixture$diagonal), tolerance = 1e-12
  )

  # Against the comparator the fold declares, it conserves exactly, and
  # `frame_conservation()` certifies the same frame the same way.
  folded <- crossform:::.metric_additive_frame(
    fixture$conservative,
    crossform:::.geometry_metric_schedule(fixture$conservative, diagonal)
  )
  folded_certificate <- crossform:::.metric_frame_conservation(folded)
  expect_true(folded_certificate$identity_conservation)
  expect_lt(folded_certificate$max_deviation, 1e-12)
  report <- frame_conservation(folded)
  expect_true(report$conserved)
  expect_identical(report$declared_normalization, "conservative")
  expect_identical(folded$metric_folded$composition, "diagonal_metric_fold")
  # `frame_conservation()` names the law its `reference_mass` is a mass of.
  expect_identical(report$composition, "diagonal_metric_fold")
  expect_identical(
    frame_conservation(fixture$conservative)$composition, "none"
  )
  # A whitened plan composes nothing into the frame, so the frame it measures
  # at is a plain conservative one and says so.
  whitened_plan <- plan_geometry(
    fixture$relation, fixture$conservative, fixture$over,
    metric = metric, composition = "whitened"
  )
  expect_identical(
    frame_conservation(whitened_plan$frame)$composition, "none"
  )
  expect_identical(whitened_plan$metric_schedule$composition, "whitened")

  # And at the geometry level, which is the claim the contract's section 5
  # actually makes: a diagonal metric conserves natively, a dense one does not.
  global <- drop(composition_total(fixture, fixture$global, diagonal)$total)
  local <- composition_total(fixture, fixture$conservative, diagonal)$total
  expect_lt(max(abs(colSums(local) - global)), 1e-12)
})

test_that("a diagonal metric composes identically both ways, and is still two plans", {
  fixture <- composition_fixture()
  metric <- noise_precision(fixture$diagonal, fixture$domain)

  native <- composition_total(fixture, fixture$conservative, metric, "native")
  whitened <- composition_total(
    fixture, fixture$conservative, metric, "whitened"
  )

  # D(sqrt(q)) D(w) D(sqrt(q)) = D(w q) = D(sqrt(w)) D(q) D(sqrt(w)), so the
  # numbers coincide exactly for a diagonal metric.
  expect_equal(
    unname(whitened$total), unname(native$total), tolerance = 1e-12
  )
  # The estimands still differ as declarations, because the reader of a plan
  # cannot know the metric was diagonal without looking, and the composition
  # is what the analysis record has to carry.
  expect_false(identical(
    native$plan$scientific_plan_id, whitened$plan$scientific_plan_id
  ))
  # Whitening a diagonal metric is feature additive, and the certificate says
  # so rather than reporting the dense answer.
  certificate <- crossform:::.metric_frame_conservation(
    fixture$conservative, composition = "whitened", metric = metric
  )
  expect_true(certificate$feature_additive)
  expect_true(certificate$identity_conservation)
  expect_identical(certificate$global_metric_kind, "native_diagonal")
})

test_that("the whitened composition refuses what it cannot root", {
  fixture <- composition_fixture()

  # A singular metric is positive semidefinite, so `neural_metric()` accepts
  # it; whitening cannot, because the null direction is a coordinate the
  # whitened frame has no way to represent. The refusal names the eigenvalue.
  singular <- fixture$dense
  decomposition <- eigen(singular, symmetric = TRUE)
  values <- decomposition$values
  values[[length(values)]] <- 0
  singular <- decomposition$vectors %*% (values * t(decomposition$vectors))
  singular <- (singular + t(singular)) / 2
  metric <- neural_metric(singular, fixture$domain)
  expect_false(metric$capabilities$positive_definite)
  expect_error(
    plan_geometry(
      fixture$relation, fixture$conservative, fixture$over,
      metric = metric, composition = "whitened"
    ),
    "positive-definite metric.*Eigenvalue",
    class = "effect_input_error"
  )

  # A support-local metric whitens nothing outside its own support, so the
  # composition admits a domain-wide operator only.
  local_metric <- neural_metric(
    diag(3), fixture$domain,
    support = fixture$domain$feature_ids[1:3]
  )
  expect_error(
    plan_geometry(
      fixture$relation, fixture$conservative, fixture$over,
      metric = local_metric, composition = "whitened"
    ),
    "domain-wide metric", class = "effect_input_error"
  )

  # And it needs a metric at all: whitened coordinates are defined by one.
  expect_error(
    plan_geometry(
      fixture$relation, fixture$conservative, fixture$over,
      composition = "whitened"
    ),
    "needs one", class = "effect_input_error"
  )
})

test_that("a learned local metric refuses the whitened composition", {
  setup <- metric_learning_setup()
  refusal <- catch_refusal(plan_geometry(
    setup$fixture$fit, setup$fixture$frame, setup$over,
    metric = shrinkage_precision(0.2), composition = "whitened",
    residual_workspace_bytes = setup$budgets$wider
  ))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "whitened_metric_composition")
  expect_identical(refusal$namespace, "geometry_plans")
  expect_identical(
    refusal$reasons, "learned_local_metric_has_no_global_root"
  )
  # The refusal explains the mathematics rather than the implementation: a
  # per-support operator has no single global root to whiten by.
  expect_match(conditionMessage(refusal), "no single domain-wide operator")

  # The same recipe compiles under the default composition, so this is a
  # refusal of the pairing and not of the recipe.
  expect_s3_class(
    plan_geometry(
      setup$fixture$fit, setup$fixture$frame, setup$over,
      metric = shrinkage_precision(0.2),
      residual_workspace_bytes = setup$budgets$wider
    ),
    "effect_geometry_plan"
  )
})

test_that("a whitened plan declares the transform it performed", {
  fixture <- composition_fixture()
  metric <- neural_metric(fixture$dense, fixture$domain)
  plan <- plan_geometry(
    fixture$relation, fixture$conservative, fixture$over,
    metric = metric, composition = "whitened"
  )

  # Whitening is not support local, so it does not go through the frame; what
  # executes is the implicit-identity lowering on whitened coordinates.
  expect_identical(plan$lowering, "additive_contraction")
  expect_identical(
    plan$metric_schedule$materialization, "whitened_effect_coordinates"
  )
  expect_true(plan$metric_schedule$feature_additive)
  expect_error(
    crossform:::.metric_additive_frame(
      fixture$conservative, plan$metric_schedule
    ),
    "no diagonal to fold", class = "effect_contract_error"
  )

  # The transform is visible in the plan, and its cost is stated in bytes: one
  # whitened effect-by-feature matrix per partition, plus the root.
  hints <- plan$execution_hints
  expect_true(hints$plan_time_effect_whitening)
  expect_equal(
    hints$whitened_effect_bytes,
    8 * 2 * fixture$q * fixture$domain$n_features
  )
  expect_equal(
    hints$metric_root_bytes, 8 * fixture$domain$n_features^2
  )
  expect_identical(
    plan$task$left_relation$provenance$construction,
    "whitened_effect_coordinates"
  )
  expect_identical(
    plan$task$left_relation$provenance$metric, metric$signature
  )

  # A declared workspace budget is enforced against the resident cost before
  # the source is read, because a global congruence cannot be streamed away.
  expect_error(
    plan_geometry(
      fixture$relation, fixture$conservative, fixture$over,
      compute = compute_policy(workspace_bytes = 64),
      metric = metric, composition = "whitened"
    ),
    "resident bytes", class = "effect_input_error"
  )

  # Reading the source in feature blocks does not change the numbers it
  # produces: the block size bounds the transient, never the estimand.
  blocked <- plan_geometry(
    fixture$relation, fixture$conservative, fixture$over,
    compute = compute_policy(workspace_bytes = 32 * 1024^2),
    metric = metric, composition = "whitened"
  )
  expect_identical(blocked$scientific_plan_id, plan$scientific_plan_id)
  expect_equal(
    unname(geometry_component(materialize_geometry(blocked), "total")),
    unname(geometry_component(materialize_geometry(plan), "total")),
    tolerance = 1e-12
  )
})
