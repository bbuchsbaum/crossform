# Laws stated in `design/conservative-geometry-contract.md`.
#
# Existing coverage this file deliberately does not duplicate:
#   * conservation of `total` for one conservative frame -- test-workflow.R
#     "public frame laws preserve point decomposition and global total", and
#     test-integrity-guards.R "conservative frames conserve total evidence;
#     local frames do not" (which also asserts coherent non-conservation);
#   * column mass of compiled conservative frames -- test-frame.R
#     "conservative frames partition global feature mass";
#   * the non-diagonal metric-conservation certificate -- test-metric.R
#     "conservative identity conservation is capability-gated";
#   * Frobenius consistency of the packed codec -- test-effect-form-laws.R
#     "the packed symmetric codec equals the oracle svec";
#   * coherence-fraction masking -- test-views.R "contrast returns one exact
#     decomposition and signed marginal".
#
# What is new here is the *frame family* law (claim 3) and the end-to-end
# numerical behaviour of the metric composition law (claim 5), neither of
# which any existing test exercises.

contract_fixture <- function(n = 9L, seed = 4471L, q = 3L,
                             id = "conservative-contract") {
  domain <- abstract_domain(n, coordinates = cbind(seq_len(n) - 1, 0), id = id)
  set.seed(seed)
  mk <- function() {
    m <- matrix(stats::rnorm(q * n), q, n)
    rownames(m) <- letters[seq_len(q)]
    m
  }
  effects <- list(run1 = mk(), run2 = mk())
  relation_value <- relation(effects, domain = domain)
  list(
    domain = domain, effects = effects, relation = relation_value,
    over = cross_partitions(relation_value),
    left = effects$run1, right = effects$run2,
    width = q * (q + 1L) / 2L
  )
}

# The cross-partition estimator for two partitions is sym(B1 K B2^T), read in
# the packed codec the package stores.
contract_form <- function(fixture, K) {
  crossform:::.svec_symmetric(0.5 * (
    fixture$left %*% K %*% t(fixture$right) +
    fixture$right %*% K %*% t(fixture$left)
  ))
}

contract_total <- function(fixture, frame, metric = NULL) {
  geometry_component(
    materialize_geometry(
      plan_geometry(fixture$relation, frame, fixture$over, metric = metric)
    ),
    "total"
  )
}

test_that("an alpha-weighted frame family conserves total scale by scale", {
  # Claim 3. Because every scale is column normalized, stacking scales with
  # weights alpha summing to one gives sum_{x in scale s} G_x = alpha_s *
  # G_Omega for EVERY s. The per-scale total is therefore fixed by alpha and
  # carries no information about the data.
  fixture <- contract_fixture(12L, id = "conservative-contract-family")
  domain <- fixture$domain

  scales <- list(
    point = voxelwise("conservative"),
    narrow = searchlights(1.01, "conservative"),
    wide = searchlights(2.01, "conservative")
  )
  frames <- lapply(scales, compile_frame, domain = domain)
  alpha <- c(point = 0.2, narrow = 0.5, wide = 0.3)
  expect_equal(sum(alpha), 1, tolerance = 0)

  family <- additive_frame(
    do.call(rbind, lapply(names(frames), function(s) {
      alpha[[s]] * as.matrix(frames[[s]]$weights)
    })),
    normalization = "conservative", domain = domain
  )
  expect_true(frame_conservation(family)$conserved)

  global <- drop(contract_total(fixture, compile_frame(whole_brain("none"),
    domain)))
  total <- contract_total(fixture, family)

  offset <- 0L
  for (s in names(frames)) {
    rows <- offset + seq_len(nrow(frames[[s]]$weights))
    offset <- offset + length(rows)
    expect_equal(colSums(total[rows, , drop = FALSE]), alpha[[s]] * global,
      tolerance = 1e-12)
  }
  # The whole family conserves as well, since the alphas sum to one.
  expect_equal(colSums(total), global, tolerance = 1e-12)
})

test_that("a conservative total is the frame contraction of a voxel ledger", {
  # Claim 3, first half: total = sum_v w_xv b_v^(L) b_v^(R)T. The ledger is
  # built here without calling any crossform contraction, so agreement means
  # the multiscale panel really is a reallocation of one univariate map.
  fixture <- contract_fixture(10L, seed = 991L,
    id = "conservative-contract-ledger")
  ledger <- t(vapply(seq_len(fixture$domain$n_features), function(v) {
    crossform:::.svec_symmetric(0.5 * (
      tcrossprod(fixture$left[, v], fixture$right[, v]) +
      tcrossprod(fixture$right[, v], fixture$left[, v])
    ))
  }, numeric(fixture$width)))

  for (spec in list(voxelwise("conservative"),
                    searchlights(1.01, "conservative"),
                    searchlights(2.51, "conservative"))) {
    frame <- compile_frame(spec, fixture$domain)
    expect_equal(
      unname(contract_total(fixture, frame)),
      unname(as.matrix(frame$weights) %*% ledger),
      tolerance = 1e-12
    )
  }
})

test_that("conservation survives a diagonal metric and fails for a dense one", {
  # Claim 5. The native composition is D(sqrt(w)) K D(sqrt(w)). For diagonal
  # K this is D(w q), still feature additive, so sum_x w_xv = 1 closes the
  # ledger. For dense K the off-diagonal entries pick up
  # S_uv = sum_x sqrt(w_xu w_xv), which is one only on the diagonal.
  fixture <- contract_fixture(9L, seed = 20260817L,
    id = "conservative-contract-metric")
  domain <- fixture$domain
  frame <- compile_frame(searchlights(1.01, "conservative"), domain)
  global_frame <- compile_frame(whole_brain("none"), domain)
  weights <- as.matrix(frame$weights)

  q_diag <- stats::runif(domain$n_features, 0.5, 2.5)
  diagonal <- neural_metric(diag(q_diag), domain)
  expect_true(metric_capabilities(diagonal)$feature_additive)
  expect_equal(
    colSums(contract_total(fixture, frame, diagonal)),
    drop(contract_total(fixture, global_frame, diagonal)),
    tolerance = 1e-12
  )

  a <- matrix(stats::rnorm(domain$n_features^2), domain$n_features)
  dense_value <- crossprod(a) / domain$n_features + diag(domain$n_features)
  dense <- neural_metric(dense_value, domain)
  expect_false(metric_capabilities(dense)$feature_additive)

  local_total <- colSums(contract_total(fixture, frame, dense))
  global_total <- drop(contract_total(fixture, global_frame, dense))
  # Not a tolerance question: the deviation is a fixed fraction of the signal.
  expect_gt(
    max(abs(local_total - global_total)) / max(abs(global_total)), 0.01
  )

  # The exact algebraic law for the failure.
  overlap <- matrix(0, domain$n_features, domain$n_features)
  for (x in seq_len(nrow(weights))) {
    overlap <- overlap + tcrossprod(sqrt(weights[x, ]))
  }
  expect_equal(diag(overlap), rep(1, domain$n_features), tolerance = 1e-12)
  expect_gt(max(abs(overlap[upper.tri(overlap)] - 1)), 0.1)
  expect_equal(unname(local_total),
    contract_form(fixture, overlap * dense_value), tolerance = 1e-12)
})

test_that("the whitened composition conserves where the native one does not", {
  # Claim 5, second half. crossform has no code path for
  # Q^(1/2) D(w) Q^(1/2), so this is a first-principles statement of the
  # estimand a `composition = "whitened"` option would have to compute. It is
  # a DIFFERENT estimand: the per-node values disagree with the native ones.
  fixture <- contract_fixture(9L, seed = 3313L,
    id = "conservative-contract-whitened")
  frame <- compile_frame(searchlights(1.01, "conservative"), fixture$domain)
  weights <- as.matrix(frame$weights)

  a <- matrix(stats::rnorm(fixture$domain$n_features^2),
    fixture$domain$n_features)
  value <- crossprod(a) / fixture$domain$n_features +
    diag(fixture$domain$n_features)
  decomposition <- eigen(value, symmetric = TRUE)
  root <- decomposition$vectors %*%
    diag(sqrt(decomposition$values)) %*% t(decomposition$vectors)
  expect_equal(root %*% root, value, tolerance = 1e-12)

  native <- t(vapply(seq_len(nrow(weights)), function(x) {
    contract_form(fixture, value * tcrossprod(sqrt(weights[x, ])))
  }, numeric(fixture$width)))
  whitened <- t(vapply(seq_len(nrow(weights)), function(x) {
    contract_form(fixture, root %*% diag(weights[x, ]) %*% root)
  }, numeric(fixture$width)))
  global <- contract_form(fixture, value)

  expect_equal(colSums(whitened), global, tolerance = 1e-12)
  expect_gt(max(abs(colSums(native) - global)) / max(abs(global)), 0.01)
  expect_gt(max(abs(whitened - native)) / max(abs(native)), 0.01)
})

test_that("contribution shares are undefined on the signed estimation layer", {
  # Claim 6. Crossvalidated node contributions are signed, so x / sum(x) is
  # not a share and PSD clipping is not a display choice: it changes the
  # conserved total.
  fixture <- contract_fixture(8L, seed = 303L,
    id = "conservative-contract-signed")
  frame <- compile_frame(searchlights(1.01, "conservative"), fixture$domain)
  view <- contrast_energy(
    plan_geometry(fixture$relation, frame, fixture$over),
    c(a = 1, b = -1, c = 0)
  )

  expect_true(any(view$total < 0))
  shares <- view$total / sum(view$total)
  expect_true(any(shares < 0 | shares > 1))
  expect_true(any(!view$coherence_fraction_valid))
  expect_true(all(is.na(view$coherence_fraction[!view$coherence_fraction_valid])))
  # Clipping breaks conservation, which is why it can never be silent.
  expect_gt(sum(pmax(view$total, 0)), sum(view$total))
})
