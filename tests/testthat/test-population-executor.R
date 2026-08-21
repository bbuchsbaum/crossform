# `estimate_population()` --- the query-first group executor of
# `population-form-v1`.
#
# The acceptance test of ticket E4 is an *oracle* rather than a regression
# snapshot: the naive route is written out in full below, in the other
# evaluation order, out of the package's public materialization verbs, and the
# executor has to agree with it. That is what §3's commutation claim asserts,
# and the only way to test a commutation claim is to compute both orders.
#
# The fixture has four participants on four different native frame sizes, a
# transport whose radius leaves real sink mass, and three experimental effects
# so the packed width (six) is meaningfully larger than the query bank (two).
# Heterogeneous frames are the case §3 says forces transport to precede the
# fit; a fixture where every participant shared a frame would never exercise
# it, and equal packed and query widths would hide the whole point of reading
# the query first.

pex_effects <- function() {
  effect_space(c("face", "house", "tool"), basis_id = "pop-exec:v1")
}

pex_subject <- function(id, features, gain = 1) {
  domain <- abstract_domain(features,
    coordinates = cbind(x = seq_len(features) - 1),
    feature_ids = paste0("f", seq_len(features)), id = id)
  # Deterministic, participant-specific, and genuinely signed: a conservative
  # ledger is signed by construction and a fixture of nonnegative numbers
  # would let a mass-law bug pass the budget certificate.
  values <- function(divisor) {
    set.seed(sum(as.integer(charToRaw(id))) + features +
      round(1000 * divisor))
    matrix(gain * stats::rnorm(3 * features) / divisor, 3, features,
      dimnames = list(c("face", "house", "tool"), NULL))
  }
  relation <- relation(
    list(run1 = values(1), run2 = values(1.3), run3 = values(0.8)),
    effects = pex_effects(), domain = domain
  )
  plan_geometry(relation, compile_frame(voxelwise(), domain),
    cross_partitions(relation))
}

# Group centres at 0, 4 and 9 with radius 1.5: a native node at x = 2, 6, 7 or
# 11 is two units from its nearest centre and therefore lands entirely in the
# sink, so the sink carries real mass rather than being an always-zero column
# no assertion can distinguish from an absent one. Group node 9 is out of
# reach of the two smaller participants under *any* radius --- nothing of
# theirs is closest to it --- which is what gives the density branch a group
# node reached by no native mass, and therefore an `NA` rather than a `0`.
pex_carrier <- function(features, semantics = "budget", radius = 1.5, ...) {
  anatomical_transport(
    native_coords = cbind(seq_len(features) - 1),
    group_coords = cbind(c(0, 4, 9)),
    semantics = semantics, radius = radius, ...
  )
}

pex_sizes <- c(s01 = 6L, s02 = 8L, s03 = 10L, s04 = 12L)
pex_gains <- c(s01 = 1, s02 = 1.6, s03 = 0.6, s04 = 1.2)

pex_subjects <- function(sizes = pex_sizes) {
  stats::setNames(lapply(names(sizes), function(id)
    pex_subject(id, sizes[[id]], pex_gains[[id]])), names(sizes))
}

pex_plan <- function(sizes = pex_sizes, semantics = "budget", ...) {
  plan_population(
    pex_subjects(sizes),
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      pex_carrier(sizes[[id]], semantics = semantics)),
    ...
  )
}

pex_bank <- function() {
  rbind(`face-house` = c(1, -1, 0), `face-tool` = c(1, 0, -1))
}

pex_packed <- function(bank) {
  width <- 3L * 4L / 2L
  matrix(vapply(seq_len(nrow(bank)), function(k) {
    crossform:::.svec_symmetric(tcrossprod(bank[k, ]))
  }, numeric(width)), nrow = width, dimnames = list(NULL, rownames(bank)))
}

# The naive dense route, written out rather than called: materialize each
# participant's *complete* packed geometry, transport all `p` packed
# coordinates, fit the group model on every one of them, and only then apply
# the query. This is `transport -> fit -> query`; the executor runs
# `query -> transport -> fit`. Nothing here touches `estimate_population()`.
pex_dense_oracle <- function(plan, bank, component = "total") {
  packed <- pex_packed(bank)
  labels <- names(plan$subjects)
  nodes <- nrow(plan$group_index) + 1L
  width <- nrow(packed)
  stack <- matrix(NA_real_, length(labels), nodes * width,
    dimnames = list(labels, NULL))
  for (label in labels) {
    geometry <- materialize_geometry(plan$subjects[[label]])
    rows <- geometry_component(geometry, component)
    stack[label, ] <- as.numeric(transport_values(plan$transport[[label]], rows))
  }
  terms <- plan$model$columns
  # Non-finite responses are solved around rather than through, matching the
  # driver: density returns `NA` at a group node reached by no native mass,
  # and dropping the participant instead would change the population.
  finite <- apply(is.finite(stack), 2L, all)
  coefficients <- matrix(NA_real_, length(terms), ncol(stack))
  coefficients[, finite] <- qr.coef(plan$model$qr, stack[, finite, drop = FALSE])
  dense <- array(t(coefficients), dim = c(nodes, width, length(terms)))
  # Query last: each term's `(m + 1)` by `p` coefficient block times the packed
  # query matrix.
  queried <- vapply(seq_along(terms), function(position) {
    dense[, , position] %*% packed
  }, matrix(0, nodes, nrow(bank)))
  array(queried, dim = c(nodes, nrow(bank), length(terms)),
    dimnames = list(
      c(as.character(plan$group_index$node), "<sink>"),
      rownames(bank), terms
    ))
}


# The acceptance oracle --------------------------------------------------------

test_that("the query-first executor equals the naive dense oracle", {
  plan <- pex_plan()
  bank <- pex_bank()
  fit <- estimate_population(plan, bank)
  oracle <- pex_dense_oracle(plan, bank)

  expect_identical(dim(fit$coefficients), dim(oracle))
  expect_identical(unname(dimnames(fit$coefficients)), dimnames(oracle))
  expect_identical(names(dimnames(fit$coefficients)),
    c("node", "query", "term"))
  expect_lt(max(abs(fit$coefficients - oracle)), 1e-12)

  # Not vacuously zero: the oracle carries real signal at every node and term.
  expect_gt(max(abs(oracle)), 1e-3)
})

test_that("the oracle agreement holds for the coherent and configuration ledgers", {
  plan <- pex_plan()
  bank <- pex_bank()
  for (component in c("coherent", "configuration")) {
    fit <- estimate_population(plan, bank, component = component)
    expect_lt(
      max(abs(fit$coefficients - pex_dense_oracle(plan, bank, component))),
      1e-12
    )
  }
})

test_that("the additive ledger identity survives transport and the group fit", {
  plan <- pex_plan()
  bank <- pex_bank()
  total <- estimate_population(plan, bank, component = "total")
  coherent <- estimate_population(plan, bank, component = "coherent")
  configuration <- estimate_population(plan, bank,
    component = "configuration")

  # `population-form-v1` §8.1: the two ledgers add back to the transported
  # total exactly. What fails is reading the coherent ledger as a group-node
  # common mode, which is why the names differ.
  expect_lt(max(abs(
    coherent$values + configuration$values - total$values
  )), 1e-12)
  expect_lt(max(abs(
    coherent$coefficients + configuration$coefficients - total$coefficients
  )), 1e-12)
  expect_identical(
    c(total$ledger, coherent$ledger, configuration$ledger),
    c("transported_total", "native_coherent_ledger",
      "native_configuration_ledger")
  )
})


# Commutation ------------------------------------------------------------------

test_that("query, transport and OLS commute in every definable order", {
  plan <- pex_plan()
  bank <- pex_bank()
  packed <- pex_packed(bank)
  labels <- names(plan$subjects)
  nodes <- nrow(plan$group_index) + 1L
  terms <- plan$model$columns

  # Order A: query, then transport, then fit.
  a_stack <- matrix(NA_real_, length(labels), nodes * nrow(bank),
    dimnames = list(labels, NULL))
  # Order B: transport, then fit, then query (the dense oracle above).
  b_stack <- matrix(NA_real_, length(labels), nodes * nrow(packed),
    dimnames = list(labels, NULL))
  # Order C: transport, then query, then fit.
  c_stack <- a_stack
  for (label in labels) {
    rows <- geometry_component(
      materialize_geometry(plan$subjects[[label]]), "total"
    )
    carrier <- plan$transport[[label]]
    a_stack[label, ] <- as.numeric(transport_values(carrier, rows %*% packed))
    carried <- transport_values(carrier, rows)
    b_stack[label, ] <- as.numeric(carried)
    c_stack[label, ] <- as.numeric(carried %*% packed)
  }
  reshape <- function(values, width) {
    array(t(values), dim = c(nodes, width, length(terms)))
  }
  order_a <- reshape(qr.coef(plan$model$qr, a_stack), nrow(bank))
  order_c <- reshape(qr.coef(plan$model$qr, c_stack), nrow(bank))
  dense <- reshape(qr.coef(plan$model$qr, b_stack), nrow(packed))
  order_b <- array(vapply(seq_along(terms), function(position) {
    dense[, , position] %*% packed
  }, matrix(0, nodes, nrow(bank))), dim = dim(order_a))

  # §11's tolerance for the commutation of query, transport and OLS under
  # subject-constant weights is `1e-12` absolute.
  expect_lt(max(abs(order_a - order_b)), 1e-12)
  expect_lt(max(abs(order_a - order_c)), 1e-12)

  # And the shipped executor is order A.
  expect_lt(
    max(abs(unname(estimate_population(plan, bank)$coefficients) - order_a)),
    1e-12
  )
})

test_that("a one-query bank reads the same physical query as contrast_energy()", {
  plan <- pex_plan()
  contrast <- c(1, -1, 0)
  fit <- estimate_population(plan, rbind(`face-house` = contrast))
  for (label in names(plan$subjects)) {
    native <- contrast_energy(plan$subjects[[label]], contrast)$total
    carried <- transport_values(plan$transport[[label]], native)
    expect_lt(max(abs(fit$values[, "face-house", label] - carried)), 1e-12)
  }
})


# Budget preservation ----------------------------------------------------------

test_that("transported totals over group nodes and the sink close on the native total", {
  plan <- pex_plan()
  bank <- pex_bank()
  fit <- estimate_population(plan, bank)

  for (label in names(plan$subjects)) {
    native <- geometry_component(
      materialize_geometry(plan$subjects[[label]]), "total"
    ) %*% pex_packed(bank)
    carried <- fit$values[, , label]
    # §2: the transported total *including* the sink is the native total, for
    # signed ledgers, with the tolerance scaled by the ledger's L1 norm.
    expect_lt(
      max(abs(colSums(carried) - colSums(native)) / colSums(abs(native))),
      1e-12
    )
    # And the group nodes alone are the native total minus the sink.
    sink <- carried[nrow(carried), ]
    expect_lt(
      max(abs(colSums(carried[-nrow(carried), , drop = FALSE]) -
        (colSums(native) - sink)) / colSums(abs(native))),
      1e-12
    )
    expect_gt(max(abs(sink)), 0)
  }
  expect_true(fit$receipt$budget$asserted)
  expect_lt(fit$receipt$budget$max_relative_deviation, 1e-12)
  expect_identical(fit$receipt$budget$scale, "relative_to_ledger_l1_norm")
})

test_that("the sink is a labelled row of the result and is fitted like any node", {
  plan <- pex_plan()
  fit <- estimate_population(plan, pex_bank())
  expect_identical(sum(fit$index$sink), 1L)
  expect_identical(as.character(fit$index$node[fit$index$sink]), "<sink>")
  expect_identical(nrow(fit$index), nrow(plan$group_index) + 1L)
  # Fitted rather than dropped: covariate-linked sink mass is how a
  # differentially failing transport becomes visible (§3.3).
  expect_true(all(is.finite(fit$coefficients[nrow(fit$index), , ])))
  expect_gt(max(abs(fit$coefficients[nrow(fit$index), , 1L])), 0)
  # The sink is reported in budget units under either semantics (§1.3).
  expect_identical(fit$index$units, c("budget", "budget", "budget", "budget"))
})


# The group fit ----------------------------------------------------------------

test_that("the fit uses the plan's prefactored QR and reports residuals and fitted", {
  plan <- pex_plan()
  fit <- estimate_population(plan, pex_bank())

  expect_identical(fit$residual_df, 4L - plan$model$rank)
  expect_lt(max(abs(fit$fitted + fit$residuals - fit$values)), 1e-12)
  # Fitted values are the design times the coefficients, at every node and
  # query, which is what makes this one OLS and not four.
  design <- plan$model$matrix
  nodes <- nrow(fit$index)
  for (query in dimnames(fit$values)[[2L]]) {
    estimates <- matrix(fit$coefficients[, query, , drop = FALSE], nrow = nodes)
    fitted <- matrix(fit$fitted[, query, , drop = FALSE], nrow = nodes)
    residuals <- matrix(fit$residuals[, query, , drop = FALSE], nrow = nodes)
    expect_lt(max(abs(design %*% t(estimates) - t(fitted))), 1e-12)
    # Residuals are orthogonal to the design, which is the defining property
    # of a least-squares fit and is not implied by the identity above.
    expect_lt(max(abs(crossprod(design, t(residuals)))), 1e-10)
  }
})

test_that("a group covariate is fitted at every node and query", {
  sizes <- pex_sizes
  covariates <- data.frame(age = c(24, 31, 47, 55), row.names = names(sizes))
  plan <- plan_population(
    pex_subjects(sizes),
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      pex_carrier(sizes[[id]])),
    model = ~ age, data = covariates
  )
  bank <- pex_bank()
  fit <- estimate_population(plan, bank)

  expect_identical(dimnames(fit$coefficients)[[3L]], c("(Intercept)", "age"))
  expect_identical(fit$residual_df, 2L)
  expect_lt(max(abs(fit$coefficients - pex_dense_oracle(plan, bank))), 1e-12)
})


# Normalization ----------------------------------------------------------------

test_that("`none` is the mean subject ledger in native evidence units", {
  plan <- pex_plan()
  bank <- pex_bank()
  fit <- estimate_population(plan, bank)

  expect_identical(fit$normalization, "none")
  # §4.1: `sum_j g_j = mean_i T_i` under the intercept-only model.
  totals <- fit$receipt$native_total
  expect_lt(max(abs(
    colSums(fit$coefficients[, , "(Intercept)"]) - colMeans(totals)
  )), 1e-12)
})

test_that("`unit_budget` is the unweighted mean of unit-normalized ledgers", {
  sizes <- pex_sizes
  plan <- plan_population(
    pex_subjects(sizes),
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      pex_carrier(sizes[[id]])),
    normalization = "unit_budget"
  )
  bank <- pex_bank()
  fit <- estimate_population(plan, bank)
  plain <- estimate_population(pex_plan(), bank)

  # §4.1: `sum_j g_j = 1` exactly, for every query.
  expect_lt(max(abs(colSums(fit$coefficients[, , "(Intercept)"]) - 1)), 1e-12)

  # It is a per-subject scalar rescaling of the ledger, and it commutes with
  # transport (§4.4) --- so rescaling the *transported* values reproduces it.
  totals <- plain$receipt$native_total
  for (label in names(plan$subjects)) {
    for (query in rownames(bank)) {
      expect_lt(max(abs(
        fit$values[, query, label] -
          plain$values[, query, label] / totals[label, query]
      )), 1e-12)
    }
  }

  # And it is a different estimand, not a rescaled display of the same one.
  expect_gt(
    max(abs(fit$coefficients - plain$coefficients)) /
      max(abs(plain$coefficients)),
    0.01
  )
  expect_identical(fit$receipt$normalization$budget_estimate,
    "same_data_ratio")
  expect_identical(fit$receipt$normalization$budget_floor, 0)
})

test_that("`unit_budget` marks a participant whose signed total is below the declared floor", {
  sizes <- pex_sizes
  plan <- plan_population(
    pex_subjects(sizes),
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      pex_carrier(sizes[[id]])),
    normalization = "unit_budget"
  )
  # §4.3: the divisor is a signed estimate, so the criterion is scale-free.
  # A floor of one is unreachable --- `|T_i| <= sum_x |c_ix|` always --- so
  # every participant is refused and nothing is estimated, which is the
  # refusal mechanism the contract requires rather than a divergent share.
  fit <- estimate_population(plan, pex_bank(), budget_floor = 1)
  expect_true(all(is.na(fit$coefficients)))
  expect_true(all(is.na(fit$values)))
  expect_false(any(fit$receipt$subjects$contributes))
  expect_identical(fit$receipt$unresolved_columns, 8L)
  expect_identical(fit$receipt$normalization$floor_criterion,
    "abs(native_total) > budget_floor * sum(abs(native_ledger))")
  # One admission per participant and query, reported as admissions rather
  # than exclusions so the ordinary case is a table of `TRUE` and not a gap.
  expect_identical(dim(fit$receipt$normalization$admitted), c(4L, 2L))
  expect_false(any(fit$receipt$normalization$admitted))
  expect_true(all(estimate_population(plan,
    pex_bank())$receipt$normalization$admitted))

  # The floor is part of the estimand, because it decides which participants
  # contributed at all.
  expect_false(identical(
    fit$scientific_plan_id,
    estimate_population(plan, pex_bank())$scientific_plan_id
  ))
})

test_that("`unit_budget` is refused on density semantics", {
  sizes <- pex_sizes
  plan <- plan_population(
    pex_subjects(sizes),
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      pex_carrier(sizes[[id]], semantics = "density")),
    normalization = "unit_budget"
  )
  refusal <- catch_refusal(estimate_population(plan, pex_bank()))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "budget_normalization_semantics")
  expect_true("unit_budget_requires_budget_semantics" %in% refusal$reasons)
})

test_that("`precision_weighted` is refused at plan construction, so it never reaches here", {
  sizes <- pex_sizes
  refusal <- catch_refusal(plan_population(
    pex_subjects(sizes),
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      pex_carrier(sizes[[id]])),
    normalization = "precision_weighted"
  ))
  expect_identical(refusal$capability, "precision_weighted_normalization")
  # No sealed plan can carry the mode, so the executor has no branch for it.
  expect_false("precision_weighted" %in%
    vapply(list(pex_plan()), function(plan) plan$normalization, character(1)))
})


# Density semantics ------------------------------------------------------------

test_that("density transports the same values through the declared ratio", {
  sizes <- pex_sizes
  plan <- plan_population(
    pex_subjects(sizes),
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      pex_carrier(sizes[[id]], semantics = "density"))
  )
  bank <- pex_bank()
  fit <- estimate_population(plan, bank)

  expect_identical(fit$semantics, "density")
  # Group rows are densities; the sink stays in budget units (§1.3).
  expect_identical(fit$index$units,
    c("density", "density", "density", "budget"))
  # Density conserves nothing, so no certificate is asserted and the receipt
  # says so rather than reporting a number it did not check.
  expect_false(fit$receipt$budget$asserted)
  expect_true(is.na(fit$receipt$budget$max_relative_deviation))

  # Group node 9 is reached by no native node of the two smaller
  # participants, so its density is `NA` and never `0`: zero is a
  # measurement and absence is not (section 1.3). The whole node-query cell
  # is unresolved, because a participant cannot be dropped to make one
  # estimable.
  expect_true(anyNA(fit$values))
  expect_true(all(is.na(fit$coefficients["group3", , ])))
  expect_true(all(is.finite(fit$coefficients["group1", , ])))
  expect_gt(fit$receipt$unresolved_columns, 0L)

  oracle <- pex_dense_oracle(plan, bank)
  expect_identical(as.vector(is.finite(fit$coefficients)),
    as.vector(is.finite(oracle)))
  finite <- is.finite(oracle)
  expect_lt(max(abs(fit$coefficients[finite] - oracle[finite])), 1e-12)
})


# Uncertainty ------------------------------------------------------------------

test_that("per-subject query-bank covariances are carried untransported with a refusal", {
  # A sampling covariance needs an error channel, so the participants here are
  # built through `estimate_relation()` rather than the bare relation the rest
  # of this file uses.
  build <- function(id, features) {
    set.seed(nchar(id) * 7L + features)
    domain <- abstract_domain(features,
      coordinates = cbind(x = seq_len(features) - 1),
      feature_ids = paste0("f", seq_len(features)), id = id)
    runs <- c("run-1", "run-2", "run-3")
    responses <- lapply(runs, function(run)
      matrix(stats::rnorm(6 * features), 6L, features))
    names(responses) <- runs
    indices <- lapply(runs, function(run)
      observation_index(paste0(run, "-scan-", 1:6), run))
    names(indices) <- runs
    facts <- study(observations(responses, indices, domain))
    design <- cbind(face = c(1, 0, 0, 1, 0, 0), house = c(0, 1, 0, 0, 1, 0),
      tool = c(0, 0, 1, 0, 0, 1))
    designs <- lapply(runs, function(run) {
      value <- design
      rownames(value) <- paste0(run, "-scan-", 1:6)
      value
    })
    names(designs) <- runs
    target <- diag(3)
    dimnames(target) <- list(c("face", "house", "tool"), colnames(design))
    fit <- estimate_relation(plan_relation(
      facts, raw_design_model(designs), raw_effect_map(target),
      observation_model("ols", sampling_unit = "scan")
    ))
    plan <- plan_geometry(fit$relation,
      compile_frame(voxelwise(), domain),
      cross_partitions(fit$relation, independence = "independent"))
    list(fit = fit, plan = plan, features = features)
  }
  sizes <- c(u01 = 5L, u02 = 6L)
  built <- lapply(stats::setNames(names(sizes), names(sizes)),
    function(id) build(id, sizes[[id]]))
  plan <- plan_population(
    lapply(built, `[[`, "plan"),
    lapply(built, function(value) pex_carrier(value$features, radius = 3))
  )
  # A centred bank: the crossvalidated distance basis carries no information
  # about an uncentred contrast, and D8 refuses one rather than approximating.
  bank <- rbind(`face-house` = c(1, -1, 0))
  uncertainty <- lapply(built, function(value) {
    rdm_sampling_covariance(value$plan, value$fit, target = "null",
      at = seq_len(value$features))
  })

  fit <- estimate_population(plan, bank, uncertainty = uncertainty)
  expect_identical(names(fit$uncertainty$native), c("u01", "u02"))
  expect_identical(fit$uncertainty$scope, "native_node_marginal")
  # One K-by-K block per *native* node, before transport --- which is what
  # makes it honest, and why the transported covariance is refused.
  expect_identical(length(fit$uncertainty$native$u01), 5L)
  expect_identical(
    fit$uncertainty$transported$capability, "transported_sampling_covariance"
  )
  expect_true("cross_node_sampling_covariance_unavailable" %in%
    fit$uncertainty$transported$reasons)

  # An uncentred contrast is a legitimate point estimate and an illegitimate
  # covariance, so the bank passes without `uncertainty` and refuses with it.
  expect_s3_class(
    estimate_population(plan, rbind(mean = c(1, 1, 1))),
    "effect_population_result"
  )
  refusal <- catch_refusal(estimate_population(plan, rbind(mean = c(1, 1, 1)),
    uncertainty = uncertainty))
  expect_identical(refusal$capability, "distance_basis_query_bank")

  misnamed <- uncertainty
  names(misnamed) <- c("u01", "u99")
  expect_identical(
    catch_refusal(estimate_population(plan, bank,
      uncertainty = misnamed))$capability,
    "aligned_subject_uncertainty"
  )
})


# Refusals and input guards ----------------------------------------------------

test_that("the executor refuses anything that is not a sealed population plan", {
  expect_error(estimate_population(list(), pex_bank()),
    "effect_population_plan", class = "effect_input_error")
  expect_error(
    estimate_population(pex_plan()$subjects$s01, pex_bank()),
    "effect_population_plan", class = "effect_input_error"
  )
})

test_that("a malformed query bank is refused rather than guessed at", {
  plan <- pex_plan()
  expect_error(estimate_population(plan), "queries",
    class = "effect_input_error")
  expect_error(estimate_population(plan, "face-house"),
    class = "effect_input_error")
  # Wrong width against the declared effect space.
  expect_error(estimate_population(plan, rbind(c(1, -1))),
    class = "effect_input_error")
  expect_error(estimate_population(plan, matrix(numeric(0), 0L, 3L)),
    class = "effect_input_error")
  # Duplicate query names cannot label a bank.
  expect_error(
    estimate_population(plan, rbind(a = c(1, -1, 0), a = c(1, 0, -1))),
    class = "effect_input_error"
  )
  expect_error(estimate_population(plan, pex_bank(), budget_floor = -1),
    class = "effect_input_error")
})

test_that("a group index carrying the result's own columns is refused", {
  sizes <- pex_sizes
  # The group index is an open table, so `sink` and `units` are only free by
  # convention; overwriting a caller's declared column would replace metadata
  # with bookkeeping and leave no trace of the substitution.
  group <- data.frame(node = c("g1", "g2", "g3"), units = c("mm", "mm", "mm"),
    stringsAsFactors = FALSE)
  plan <- plan_population(
    pex_subjects(sizes),
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      anatomical_transport(
        native_coords = cbind(seq_len(sizes[[id]]) - 1),
        group_coords = cbind(c(0, 4, 9)), semantics = "budget",
        radius = 1.5, group_index = group))
  )
  refusal <- catch_refusal(estimate_population(plan, pex_bank()))
  expect_identical(refusal$capability, "reserved_group_index_columns")
  expect_identical(refusal$reasons, "group_index_column_reserved:units")
})

test_that("a named bank names the queries and an unnamed one is numbered", {
  plan <- pex_plan()
  named <- estimate_population(plan, pex_bank())
  expect_identical(dimnames(named$coefficients)[[2L]],
    c("face-house", "face-tool"))
  unnamed <- estimate_population(plan, rbind(c(1, -1, 0), c(1, 0, -1)))
  expect_identical(dimnames(unnamed$coefficients)[[2L]],
    c("query1", "query2"))
  expect_lt(max(abs(unname(named$coefficients) -
    unname(unnamed$coefficients))), 1e-12)
})


# The record -------------------------------------------------------------------

test_that("the result is a sealed record with a receipt that names every read", {
  plan <- pex_plan()
  fit <- estimate_population(plan, pex_bank())

  expect_s3_class(fit, "effect_population_result")
  expect_silent(crossform:::.validate_population_result(fit))
  expect_match(fit$scientific_plan_id, "^population-sha256:")

  receipt <- fit$receipt
  expect_identical(receipt$contract_version, "population-form-v1")
  expect_identical(receipt$population_plan_id, plan$scientific_plan_id)
  expect_identical(receipt$evaluation_order, "transport_then_fit")
  expect_identical(nrow(receipt$subjects), 4L)
  expect_identical(names(receipt$subject_receipts), names(plan$subjects))
  for (value in receipt$subject_receipts) {
    expect_s3_class(value, "effect_execution_receipt")
  }
  expect_identical(receipt$subjects$transport_signature,
    plan$subject_index$transport_signature)
  expect_identical(dim(receipt$native_total), c(4L, 2L))
  expect_identical(dim(receipt$sink_budget), c(4L, 2L))

  # The identity moves with the bank and the component, because both name a
  # different group quantity.
  expect_false(identical(fit$scientific_plan_id,
    estimate_population(plan, pex_bank()[1L, , drop = FALSE])$scientific_plan_id))
  expect_false(identical(fit$scientific_plan_id,
    estimate_population(plan, pex_bank(),
      component = "coherent")$scientific_plan_id))
  # And it does not move when nothing about the estimand did.
  expect_identical(fit$scientific_plan_id,
    estimate_population(plan, pex_bank())$scientific_plan_id)
})

test_that("the result prints its ledger, its frame and its transport identity", {
  fit <- estimate_population(pex_plan(), pex_bank())
  output <- utils::capture.output(print(fit))
  expect_true(any(grepl("effect_population_result", output, fixed = TRUE)))
  # §8.1's required print line: the ledger name, the native frame the ledger
  # belongs to, and the transport identity with its cross-fit status.
  expect_true(any(grepl("transported_total", output, fixed = TRUE)))
  expect_true(any(grepl("conservative", output, fixed = TRUE)))
  expect_true(any(grepl("budget, anatomical", output, fixed = TRUE)))
  expect_true(any(grepl("cross-fit", output, fixed = TRUE)))
  expect_true(any(grepl("preserved", output, fixed = TRUE)))
  expect_match(format(fit), "^<effect_population_result: transported_total")

  # A transported component carries its own name and says what it is not.
  coherent <- utils::capture.output(
    print(estimate_population(pex_plan(), pex_bank(), component = "coherent"))
  )
  expect_true(any(grepl("native_coherent_ledger", coherent, fixed = TRUE)))
  expect_true(any(grepl("not a group-node common mode", coherent,
    fixed = TRUE)))
  # The bare names are forbidden as the ledger itself; naming the
  # `component` argument that selected it is not the same thing.
  expect_false(any(grepl("^  ledger: *coherent", coherent)))
  expect_false(any(grepl("^  ledger: *configuration", coherent)))
})

test_that("as.data.frame gives the coefficient table in long form", {
  fit <- estimate_population(pex_plan(), pex_bank())
  table <- as.data.frame(fit)
  expect_identical(nrow(table), 4L * 2L * 1L)
  expect_true(all(c("node", "sink", "units", "query", "term", "estimate") %in%
    names(table)))
  expect_identical(unique(table$query), c("face-house", "face-tool"))
  expect_identical(unique(table$term), "(Intercept)")
  expect_identical(table$estimate, as.numeric(fit$coefficients))
  expect_identical(sum(table$sink), 2L)
})
