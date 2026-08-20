# Query banks over the sampling-covariance distance basis -- ticket D8.
#
# `sampling_covariance(x, queries = )` names a bank of K contrast-energy
# queries and returns their K-by-K analytic sampling covariance at one
# measurement. The estimand is the one `contrast_energy()` reports as `$total`:
# the crossvalidated quadratic form c'Gc.
#
# The claim the acceptance rests on is that this is not a second law but the
# existing one read in other coordinates. A centred contrast satisfies
#
#   c c' = - sum_{i<j} c_i c_j (e_i - e_j)(e_i - e_j)',
#
# so `svec(c c')` -- the physical query `contrast_energy()` lowers to -- lies
# in the span of the packed distance operators, and the bank's covariance is
# the transport of the distance-basis covariance through that lowering. Every
# agreement below is asserted against a lowering this file solves for itself,
# in packed coordinates, with base R only: the package's closed form is
# checked, not consulted.

query_bank_oracle_svec <- function(x) {
  q <- nrow(x)
  out <- numeric(q * (q + 1L) / 2L)
  k <- 0L
  for (column in seq_len(q)) {
    for (row in column:q) {
      k <- k + 1L
      out[[k]] <- x[row, column] * if (row == column) 1 else sqrt(2)
    }
  }
  out
}

# The lowering, solved rather than written down: express each query's packed
# operator in the packed distance basis. The system is square-free and of full
# column rank, so a solution exists exactly when the query is expressible, and
# the residual reported alongside is what "expressible" means numerically.
query_bank_oracle_lowering <- function(bank, effects) {
  pairs <- utils::combn(length(effects), 2L)
  packed_distance <- vapply(seq_len(ncol(pairs)), function(distance) {
    weights <- numeric(length(effects))
    weights[pairs[1L, distance]] <- 1
    weights[pairs[2L, distance]] <- -1
    query_bank_oracle_svec(tcrossprod(weights))
  }, numeric(length(effects) * (length(effects) + 1L) / 2L))
  map <- t(vapply(seq_len(nrow(bank)), function(query) {
    qr.solve(packed_distance, query_bank_oracle_svec(tcrossprod(bank[query, ])))
  }, numeric(ncol(pairs))))
  residual <- vapply(seq_len(nrow(bank)), function(query) {
    max(abs(
      packed_distance %*% map[query, ] -
        query_bank_oracle_svec(tcrossprod(bank[query, ]))
    ))
  }, numeric(1))
  dimnames(map) <- list(rownames(bank), paste(
    effects[pairs[1L, ]], effects[pairs[2L, ]], sep = " - "
  ))
  list(map = map, residual = residual)
}

query_bank_fixture <- function(conditions = 4L, features = 8L,
                               partitions = 4L) {
  set.seed(90211)
  observations <- 64L
  condition <- factor(rep(seq_len(conditions), length.out = observations))
  design <- stats::model.matrix(~ 0 + condition)
  colnames(design) <- paste0("condition", seq_len(conditions))
  effects <- diag(conditions)
  rownames(effects) <- colnames(design)
  # An AR-like residual covariance read under an identity metric, so the
  # whitened residual covariance stays anisotropic and the signal term's
  # dependence on it is identifiable rather than collapsing to a scalar.
  covariance <- stats::toeplitz(0.4^(0:(features - 1L)))
  factor <- chol(covariance)
  signal <- matrix(stats::rnorm(conditions * features, sd = 0.3),
    conditions, features)
  domain <- abstract_domain(
    features, coordinates = cbind(x = seq_len(features), y = 0, z = 0),
    id = "query-bank-domain", coordinate_units = "mm"
  )
  sources <- stats::setNames(lapply(seq_len(partitions), function(run) {
    design %*% signal +
      matrix(stats::rnorm(observations * features), observations, features) %*%
        factor
  }), paste0("run", seq_len(partitions)))
  fit <- lm_relation_fit(
    sources, design, effects, sampling_unit = "trial", domain = domain
  )
  metric <- noise_precision(
    diag(features), domain, covariance = diag(features),
    provenance = list(source = "query-bank-fixture")
  )
  list(
    fit = fit, domain = domain, metric = metric,
    effects = fit$relation$effects,
    plan = function(frame) {
      plan_geometry(
        fit$relation, frame,
        cross_partitions(fit$relation, independence = "independent"),
        metric = metric
      )
    }
  )
}

query_bank_matrix <- function(effects) {
  bank <- rbind(
    pairwise = c(1, -1, 0, 0),
    overlapping = c(1, 0, -1, 0),
    unbalanced = c(0.5, 0.5, -0.25, -0.75)
  )
  colnames(bank) <- effects
  bank
}

test_that("a query bank reproduces the transport it names", {
  fixture <- query_bank_fixture()
  plan <- fixture$plan(compile_frame(whole_brain(), fixture$domain))
  # `plugin` keeps the signal term alive: under `null` it vanishes and only
  # the quadratic noise term would be compared.
  distances <- rdm_sampling_covariance(
    plan, fixture$fit, target = "plugin", at = 1L
  )
  bank <- query_bank_matrix(fixture$effects)
  oracle <- query_bank_oracle_lowering(bank, fixture$effects)
  expect_lt(max(oracle$residual), 1e-12)

  energies <- sampling_covariance(distances, queries = bank)
  expect_s3_class(energies, "effect_sampling_covariance")
  expect_identical(energies$basis, "query_bank")
  expect_identical(energies$labels, rownames(bank))
  expect_identical(energies$dimension, 3L)

  # The lowering the package built is the lowering the packed system has.
  expect_equal(energies$source$query_lowering, oracle$map, tolerance = 1e-12)

  block <- sampling_covariance(energies, "materialize")
  transported <- sampling_covariance(
    distances, "transport", query = oracle$map
  )
  expect_identical(dim(block), c(3L, 3L))
  expect_identical(dimnames(block), list(rownames(bank), rownames(bank)))
  expect_lt(max(abs(block - transported)), 1e-12)
  # The off-diagonal is not incidentally zero, so the agreement is an
  # agreement about covariance and not only about variances.
  expect_gt(max(abs(block[upper.tri(block)])), 0)
})

test_that("a one-query bank is the distance it selects", {
  fixture <- query_bank_fixture()
  plan <- fixture$plan(compile_frame(whole_brain(), fixture$domain))
  distances <- rdm_sampling_covariance(
    plan, fixture$fit, target = "plugin", at = 1L
  )
  bank <- matrix(c(1, -1, 0, 0), nrow = 1L,
    dimnames = list("first pair", fixture$effects))
  oracle <- query_bank_oracle_lowering(bank, fixture$effects)

  energies <- sampling_covariance(distances, queries = bank)
  block <- sampling_covariance(energies, "materialize")
  expect_identical(dim(block), c(1L, 1L))
  expect_lt(max(abs(
    block - sampling_covariance(distances, "transport", query = oracle$map)
  )), 1e-12)

  # `e_i - e_j` is one of the basis coordinates, so its energy variance is
  # that distance's own variance, exactly.
  expect_equal(
    unname(block[1L, 1L]),
    unname(sampling_covariance(distances)[[distances$labels[[1L]]]]),
    tolerance = 1e-14
  )

  # A bare contrast vector is the same bank, numbered rather than named.
  bare <- sampling_covariance(distances, queries = c(1, -1, 0, 0))
  expect_identical(bare$labels, "query1")
  expect_equal(unname(sampling_covariance(bare, "materialize")),
    unname(block), tolerance = 1e-14)
})

test_that("the bank's coordinates are the energies contrast_energy reports", {
  # The covariance is of an estimand only if the estimand is the one the
  # package actually reports. The lowering applied to the crossvalidated
  # distances must reproduce `contrast_energy(plan, c)$total` exactly, or the
  # K-by-K block describes something the user never sees.
  fixture <- query_bank_fixture()
  plan <- fixture$plan(compile_frame(whole_brain(), fixture$domain))
  distances <- rdm_sampling_covariance(
    plan, fixture$fit, target = "plugin", at = 1L
  )
  bank <- query_bank_matrix(fixture$effects)
  oracle <- query_bank_oracle_lowering(bank, fixture$effects)

  observed <- rdm(plan)$values[1L, distances$labels]
  lowered <- drop(oracle$map %*% observed)
  reported <- vapply(seq_len(nrow(bank)), function(query) {
    contrast_energy(plan, bank[query, ])$total[[1L]]
  }, numeric(1))
  expect_equal(unname(lowered), unname(reported), tolerance = 1e-12)
})

test_that("each block is symmetric and positive semidefinite", {
  fixture <- query_bank_fixture()
  plan <- fixture$plan(compile_frame(whole_brain(), fixture$domain))
  bank <- query_bank_matrix(fixture$effects)
  for (target in c("null", "plugin")) {
    distances <- rdm_sampling_covariance(
      plan, fixture$fit, target = target, at = 1L
    )
    block <- sampling_covariance(
      distances, "materialize", queries = bank
    )
    expect_identical(max(abs(block - t(block))), 0)
    spectrum <- eigen(block, symmetric = TRUE, only.values = TRUE)$values
    expect_gt(min(spectrum), -1e-12 * max(spectrum))
  }
})

test_that("every operation reads the same block in the new basis", {
  fixture <- query_bank_fixture()
  plan <- fixture$plan(compile_frame(whole_brain(), fixture$domain))
  distances <- rdm_sampling_covariance(
    plan, fixture$fit, target = "plugin", at = 1L
  )
  bank <- query_bank_matrix(fixture$effects)
  energies <- sampling_covariance(distances, queries = bank)
  block <- sampling_covariance(energies, "materialize")

  variances <- sampling_covariance(energies)
  expect_identical(names(variances), rownames(bank))
  expect_equal(unname(variances), unname(diag(block)), tolerance = 1e-14)

  weights <- c(1, -1, 0.5)
  expect_equal(
    sampling_covariance(energies, "quadratic_form", query = weights),
    drop(weights %*% block %*% weights), tolerance = 1e-14
  )
  expect_equal(
    unname(sampling_covariance(energies, "apply", query = weights)),
    unname(drop(block %*% weights)), tolerance = 1e-14
  )
  # Supplying the operation alongside the bank is the same as supplying it
  # afterwards; the bank changes the basis, the operation reads it.
  expect_equal(
    sampling_covariance(distances, "materialize", queries = bank),
    block, tolerance = 1e-14
  )
})

test_that("named bank columns are aligned, not taken positionally", {
  fixture <- query_bank_fixture()
  plan <- fixture$plan(compile_frame(whole_brain(), fixture$domain))
  distances <- rdm_sampling_covariance(
    plan, fixture$fit, target = "plugin", at = 1L
  )
  bank <- query_bank_matrix(fixture$effects)
  shuffled <- bank[, rev(seq_len(ncol(bank))), drop = FALSE]

  expect_equal(
    sampling_covariance(distances, "materialize", queries = bank),
    sampling_covariance(distances, "materialize", queries = shuffled),
    tolerance = 1e-14
  )
  # A list of named contrast vectors is the same bank again.
  as_list <- lapply(seq_len(nrow(bank)), function(row) bank[row, ])
  names(as_list) <- rownames(bank)
  expect_equal(
    sampling_covariance(distances, "materialize", queries = as_list),
    sampling_covariance(distances, "materialize", queries = bank),
    tolerance = 1e-14
  )
})

test_that("a conservative frame carries the same agreement", {
  fixture <- query_bank_fixture()
  frame <- compile_frame(
    searchlights(1.01, normalization = "conservative"), fixture$domain
  )
  expect_identical(frame$normalization, "conservative")
  plan <- fixture$plan(frame)
  distances <- rdm_sampling_covariance(
    plan, fixture$fit, target = "plugin", at = 3L
  )
  bank <- query_bank_matrix(fixture$effects)
  oracle <- query_bank_oracle_lowering(bank, fixture$effects)

  block <- sampling_covariance(distances, "materialize", queries = bank)
  expect_lt(max(abs(
    block - sampling_covariance(distances, "transport", query = oracle$map)
  )), 1e-12)

  # Contract 7.6a: a conserved ledger buys no error bars. The block is one
  # measurement's marginal and says so.
  energies <- sampling_covariance(distances, queries = bank)
  expect_identical(energies$source$cross_node_covariance, "unavailable")
  expect_true(energies$source$local_only)
})

test_that("a frame family answers per row, named by measurement", {
  fixture <- query_bank_fixture()
  family <- frame_family(
    narrow = compile_frame(
      searchlights(1.01, normalization = "conservative"), fixture$domain
    ),
    wide = compile_frame(
      searchlights(2.01, normalization = "conservative"), fixture$domain
    ),
    alpha = c(narrow = 0.4, wide = 0.6)
  )
  plan <- fixture$plan(family)
  rows <- c(1L, 5L, 12L)
  distances <- rdm_sampling_covariance(
    plan, fixture$fit, target = "plugin", at = rows
  )
  bank <- query_bank_matrix(fixture$effects)
  oracle <- query_bank_oracle_lowering(bank, fixture$effects)

  energies <- sampling_covariance(distances, queries = bank)
  expect_s3_class(energies, "effect_sampling_covariance_batch")
  expect_identical(attr(energies, "basis"), "query_bank")
  expect_identical(names(energies),
    as.character(family$index$measurement[rows]))

  # The accessor the population layer reads: one K-by-K block per measurement,
  # keyed by the measurement identifier rather than by position.
  blocks <- sampling_covariance(distances, "materialize", queries = bank)
  expect_identical(names(blocks), names(energies))
  for (position in seq_along(rows)) {
    expect_identical(dim(blocks[[position]]), c(3L, 3L))
    expect_lt(max(abs(
      blocks[[position]] -
        sampling_covariance(
          distances[[position]], "transport", query = oracle$map
        )
    )), 1e-12)
  }
  expect_false(isTRUE(all.equal(blocks[[1L]], blocks[[3L]])))
})

test_that("an uncentred query is refused rather than approximated", {
  fixture <- query_bank_fixture()
  plan <- fixture$plan(compile_frame(whole_brain(), fixture$domain))
  distances <- rdm_sampling_covariance(
    plan, fixture$fit, target = "plugin", at = 1L
  )
  refusal <- catch_refusal(sampling_covariance(
    distances,
    queries = rbind(
      centred = c(1, -1, 0, 0), self = c(1, 0, 0, 0)
    )
  ))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "distance_basis_query_bank")
  expect_identical(refusal$namespace, "evidence_sampling")
  expect_identical(refusal$reasons,
    "uncentred_contrast_not_in_distance_basis")
  # It names the offending query and why the basis cannot hold it.
  expect_match(conditionMessage(refusal), "`self` sums to 1")
  expect_match(conditionMessage(refusal), "invariant to a shift")
  expect_match(conditionMessage(refusal), "will not manufacture")

  # The claim behind the refusal, checked: an uncentred query's packed
  # operator genuinely leaves the packed distance span.
  oracle <- query_bank_oracle_lowering(
    matrix(c(1, 0, 0, 0), nrow = 1L,
      dimnames = list("self", fixture$effects)),
    fixture$effects
  )
  expect_gt(oracle$residual, 1e-6)
})

test_that("a bank never reaches a learned metric or a whitened plan", {
  fixture <- query_bank_fixture()
  learned <- plan_crossnobis(
    fixture$fit, compile_frame(searchlights(100), fixture$domain),
    pairing("run1", "run2", independence = "independent"),
    metric = shrinkage_precision(0.2),
    training = metric_training_policy(
      "all_partitions_residual_orthogonality",
      justification = "test fixture; uncertainty remains unpropagated"
    )
  )
  learned_refusal <- catch_refusal(
    rdm_sampling_covariance(learned, fixture$fit, target = "null", at = 1L)
  )
  expect_identical(learned_refusal$capability, "fixed_metric_sampling_law")
  expect_identical(learned_refusal$reasons, "learned_metric_law_not_admitted")

  # Contract 5.2's owed refusal: a whitened plan refuses the residual channel,
  # and now says that the composition is what unbound it.
  whitened <- plan_geometry(
    fixture$fit$relation, compile_frame(whole_brain(), fixture$domain),
    cross_partitions(fixture$fit$relation, independence = "independent"),
    metric = fixture$metric, composition = "whitened"
  )
  whitened_refusal <- catch_refusal(
    rdm_sampling_covariance(whitened, fixture$fit, target = "null", at = 1L)
  )
  expect_s3_class(whitened_refusal, "effect_capability_refusal")
  expect_identical(whitened_refusal$capability,
    "whitened_composition_sampling_law")
  expect_identical(whitened_refusal$reasons,
    "whitened_composition_error_channel_not_bound")
  expect_match(conditionMessage(whitened_refusal), "Q\\^\\(1/2\\)")
  expect_match(conditionMessage(whitened_refusal), "refusal, not an approx")
  # The same answer when asked in advance rather than provoked.
  expect_identical(
    catch_refusal(sampling_capabilities(whitened, fixture$fit))$capability,
    "whitened_composition_sampling_law"
  )

  # The native composition of the same estimand is admitted, so the refusal
  # is about the composition and not about the fixture.
  native <- plan_geometry(
    fixture$fit$relation, compile_frame(whole_brain(), fixture$domain),
    cross_partitions(fixture$fit$relation, independence = "independent"),
    metric = fixture$metric
  )
  expect_s3_class(
    rdm_sampling_covariance(native, fixture$fit, target = "null", at = 1L),
    "effect_sampling_covariance"
  )
})

test_that("bank shapes are refused by the contrast aligner", {
  fixture <- query_bank_fixture()
  plan <- fixture$plan(compile_frame(whole_brain(), fixture$domain))
  distances <- rdm_sampling_covariance(
    plan, fixture$fit, target = "plugin", at = 1L
  )
  expect_error(
    sampling_covariance(distances, queries = rbind(a = c(1, -1, 0))),
    "has 3 values but the relation declares 4 effects",
    class = "effect_input_error"
  )
  expect_error(
    sampling_covariance(distances, queries = rbind(a = c(1, -1, NA, 0))),
    "must be finite", class = "effect_input_error"
  )
  expect_error(
    sampling_covariance(distances, queries = list()),
    "at least one contrast", class = "effect_input_error"
  )
  expect_error(
    sampling_covariance(distances, queries = "animacy"),
    "numeric contrast vector", class = "effect_input_error"
  )
  expect_error(
    sampling_covariance(distances, queries = list(
      a = c(1, -1, 0, 0), a = c(1, 0, -1, 0)
    )),
    "name every query exactly once", class = "effect_input_error"
  )
  expect_error(
    sampling_covariance(distances, queries = matrix(numeric(), 0L, 4L)),
    "at least one numeric row", class = "effect_input_error"
  )
  named <- rbind(a = c(nope = 1, condition2 = -1, condition3 = 0,
    condition4 = 0))
  expect_error(
    sampling_covariance(distances, queries = named),
    "is not a declared effect", class = "effect_input_error"
  )

  # A bank is attached to a covariance, not to the plan behind it, and the
  # message says which object to build first.
  expect_error(
    sampling_covariance(plan, queries = c(1, -1, 0, 0)),
    "build one with `rdm_sampling_covariance\\(\\)` first",
    class = "effect_input_error"
  )

  # A bank is lowered onto the distance basis, so a covariance already in the
  # query basis cannot take another one.
  energies <- sampling_covariance(
    distances, queries = query_bank_matrix(fixture$effects)
  )
  rebank <- catch_refusal(
    sampling_covariance(energies, queries = c(1, -1, 0, 0))
  )
  expect_identical(rebank$capability, "distance_basis_query_bank")
  expect_identical(rebank$reasons, "covariance_is_not_in_the_distance_basis")
})

test_that("cross-measurement covariance is refused with what is missing", {
  fixture <- query_bank_fixture()
  plan <- fixture$plan(
    compile_frame(
      searchlights(1.01, normalization = "conservative"), fixture$domain
    )
  )
  distances <- rdm_sampling_covariance(
    plan, fixture$fit, target = "plugin", at = c(1L, 2L)
  )
  bank <- query_bank_matrix(fixture$effects)
  refusal <- catch_refusal(sampling_covariance(
    distances, queries = bank, scope = "cross_measurement"
  ))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "cross_node_sampling_covariance")
  expect_identical(refusal$namespace, "evidence_sampling")
  expect_identical(refusal$reasons, c(
    "spatial_covariance_model_unavailable",
    "conservation_gives_no_uncertainty"
  ))
  expect_match(conditionMessage(refusal), "cross-support residual second")
  expect_match(conditionMessage(refusal), "conservation is a law about est")

  # The same refusal answers the distance basis, so no route quietly returns
  # per-measurement margins to a caller who asked for a joint covariance.
  expect_identical(
    catch_refusal(sampling_covariance(
      distances[[1L]], "materialize", scope = "cross_measurement"
    ))$capability,
    "cross_node_sampling_covariance"
  )
})

test_that("the bank is part of the artifact, and prints as one", {
  fixture <- query_bank_fixture()
  plan <- fixture$plan(compile_frame(whole_brain(), fixture$domain))
  distances <- rdm_sampling_covariance(
    plan, fixture$fit, target = "plugin", at = 1L
  )
  bank <- query_bank_matrix(fixture$effects)
  energies <- sampling_covariance(distances, queries = bank)

  expect_identical(energies$source$queries, bank)
  expect_identical(energies$source$construction,
    "query_bank_reduction_of_distance_basis")
  expect_identical(energies$source$basis_parent, "rdm")
  expect_identical(energies$source$distance_labels, distances$labels)
  # The packed query is the operator `contrast_energy()` lowers to.
  expect_equal(
    energies$source$packed_query["pairwise", ],
    query_bank_oracle_svec(tcrossprod(bank["pairwise", ])),
    tolerance = 1e-14
  )
  # The residual channel the block inherits is reported unchanged.
  expect_identical(energies$source$residual_df, distances$source$residual_df)
  expect_identical(energies$source$noise_trace_estimator,
    distances$source$noise_trace_estimator)

  # Two banks over one measurement are two artifacts.
  other <- sampling_covariance(
    distances, queries = bank[c("pairwise", "overlapping"), , drop = FALSE]
  )
  expect_false(identical(energies$signature, other$signature))

  output <- capture.output(print(energies))
  expect_true(any(grepl("basis:\\s+query_bank", output)))
  expect_true(any(grepl("3 contrast energies, lowered onto 6 distances",
    output)))
  expect_true(any(grepl("no cross-measurement covariance", output)))
  expect_match(format(energies), "3 queries")
})
