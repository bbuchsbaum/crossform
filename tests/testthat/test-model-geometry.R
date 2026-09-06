# The model-coordinate reader (ticket M5 of the model-coordinate epic).
#
# `design/model-coordinate-geometry-contract.md` sections 4, 5 and 6, against
# `design/oracles/model-coordinate-geometry.R` sections O3, O5, O6 and O8.
# Every comparison on the right-hand side is base matrix algebra on the full
# form G read from the ORIGINAL plan and on the basis's Q and R; the reader
# never sees G, so agreement is a statement about the lowering, the packed
# query bank and the fits together.

mg_conditions <- paste0("c", 1:6)

mg_fixture <- function(seed = 20260905, p = 8L, runs = 3L, frame = "regions") {
  set.seed(seed)
  n <- length(mg_conditions)
  labels <- c(1, 1, 1, 0, 0, 0)
  category <- outer(labels, labels, function(a, b) (a - b)^2)
  dimnames(category) <- list(mg_conditions, mg_conditions)
  features <- matrix(rnorm(n * 2L), n, 2L, dimnames = list(mg_conditions, NULL))
  basis <- model_basis(list(category = category), "squared_euclidean",
    features = list(feature = features), conditions = mg_conditions)
  Q <- unname(basis$Q)
  R <- unname(basis$R)
  Tm <- Q %*% R
  H <- diag(n) - 1 / n
  orthogonal <- (H - Q %*% t(Q)) %*% rnorm(n)
  orthogonal <- orthogonal / sqrt(sum(orthogonal^2))
  truth <- Tm %*% matrix(rnorm(ncol(Tm) * p), ncol(Tm), p) +
    orthogonal %*% t(rnorm(p)) + rep(1, n) %*% t(rnorm(p))
  betas <- lapply(seq_len(runs), function(run) {
    b <- truth + matrix(rnorm(n * p, sd = 0.6), n, p)
    rownames(b) <- mg_conditions
    b
  })
  names(betas) <- paste0("run", seq_len(runs))
  domain <- abstract_domain(p, id = "model-geometry-fixture")
  relation <- relation(betas, effects = mg_conditions, domain = domain)
  at <- if (frame == "regions") {
    compile_frame(regions(rep(c("a", "b"), each = p / 2)), domain)
  } else {
    compile_frame(whole_brain(), domain)
  }
  over <- cross_partitions(relation, independence = "independent")
  plan <- plan_geometry(relation, at, over)
  full <- materialize_geometry(plan)
  packed <- geometry_component(full, "total")
  G <- lapply(seq_len(nrow(packed)), function(i) {
    crossform:::.unsvec_symmetric(packed[i, ], n)
  })
  list(basis = basis, Q = Q, R = R, T = Tm, H = H, P = Q %*% t(Q), plan = plan,
    relation = relation, domain = domain, at = at, over = over, G = G,
    betas = betas, n = n, p = p)
}

frob2 <- function(x) sum(x^2)

psd_rank <- function(S, s) {
  e <- eigen(S, symmetric = TRUE)
  keep <- seq_len(min(s, ncol(S)))
  lambda <- pmax(e$values[keep], 0)
  e$vectors[, keep, drop = FALSE] %*% (lambda * t(e$vectors[, keep, drop = FALSE]))
}

# KKT conditions of min 0.5 c'Mc - g'c over c >= 0.
kkt_holds <- function(c, M, g, tol = 1e-8) {
  gradient <- as.numeric(M %*% c) - g
  scale <- max(1, max(abs(g)))
  all(c >= -tol) && all(gradient >= -tol * scale) &&
    all(abs(c * gradient) <= tol * scale * max(1, max(abs(c))))
}

test_that("the trace decomposition is read from the original plan and closes", {
  fx <- mg_fixture()
  reading <- model_geometry(fx$plan, fx$basis)

  expect_s3_class(reading, "effect_model_geometry")
  expect_identical(length(reading$addressable), 2L)
  for (i in 1:2) {
    G <- fx$G[[i]]
    S <- t(fx$Q) %*% G %*% fx$Q
    expect_equal(reading$addressable[[i]], sum(diag(S)), tolerance = 1e-12)
    expect_equal(reading$centered_total[[i]], sum(diag(fx$H %*% G %*% fx$H)),
      tolerance = 1e-12)
    expect_equal(reading$orthogonal[[i]],
      sum(diag((fx$H - fx$P) %*% G)), tolerance = 1e-12)
  }
  expect_lt(reading$identity_error, 1e-12)
  expect_identical(reading$layer$estimation,
    c("addressable", "orthogonal", "centered_total"))
  expect_identical(reading$orthogonal_source, "original plan, operator H - P")
  # The orthogonal term exists because the fixture planted a direction
  # outside the model span; the split is not vacuous.
  expect_gt(min(abs(reading$orthogonal)), 1e-3)
  expect_identical(reading$source_scientific_plan_id,
    fx$plan$scientific_plan_id)
  expect_false(identical(reading$lowered_scientific_plan_id,
    fx$plan$scientific_plan_id))
  expect_identical(reading$basis_signature, fx$basis$signature)
  expect_identical(reading$reading, "latent descriptive layer; not for inference")
})

test_that("the nonnegative quadratic program satisfies its KKT conditions", {
  set.seed(4)
  for (trial in 1:20) {
    m <- sample(2:6, 1L)
    A <- matrix(rnorm(m * m), m, m)
    M <- crossprod(A) + diag(1e-6, m)
    g <- rnorm(m)
    fit <- crossform:::.nnls_quadratic(M, g)
    expect_true(fit$converged)
    expect_true(kkt_holds(fit$solution, M, g))
  }
  # Two variables against a brute-force grid.
  M <- matrix(c(2, 1, 1, 3), 2, 2)
  g <- c(1, -2)
  fit <- crossform:::.nnls_quadratic(M, g)
  objective <- function(c) 0.5 * drop(t(c) %*% M %*% c) - sum(g * c)
  grid <- expand.grid(a = seq(0, 2, by = 0.001), b = seq(0, 2, by = 0.001))
  best <- min(apply(grid, 1, objective))
  expect_lte(objective(fit$solution), best + 1e-6)
  expect_equal(fit$solution, c(0.5, 0), tolerance = 1e-10)
})

test_that("isotropic and diagonal fits are the quadratic programs of section 4.1", {
  fx <- mg_fixture()
  gram <- crossprod(fx$T)
  M <- gram^2
  columns <- fx$basis$columns
  isotropic <- model_geometry(fx$plan, fx$basis, structure = "isotropic")
  diagonal <- model_geometry(fx$plan, fx$basis, structure = "diagonal")
  expect_identical(colnames(isotropic$weights), c("category", "feature"))
  expect_identical(colnames(diagonal$weights),
    c("category.1", "feature.1", "feature.2"))
  expect_true(isotropic$tests_supplied_geometries)
  expect_false(diagonal$tests_supplied_geometries)
  for (i in 1:2) {
    G <- fx$G[[i]]
    S <- t(fx$Q) %*% G %*% fx$Q
    # g from the full form equals g from the compressed form (oracle O5).
    g <- vapply(seq_len(ncol(fx$T)), function(j) {
      drop(t(fx$T[, j]) %*% G %*% fx$T[, j])
    }, numeric(1))
    g_S <- vapply(seq_len(ncol(fx$R)), function(j) {
      drop(t(fx$R[, j]) %*% S %*% fx$R[, j])
    }, numeric(1))
    expect_equal(g, g_S, tolerance = 1e-10)
    c_diag <- diagonal$weights[i, ]
    expect_true(kkt_holds(c_diag, M, g))
    expect_equal(unname(diagonal$coefficients[, , i]), diag(c_diag, 3),
      tolerance = 1e-12)
    A <- fx$R %*% diag(c_diag, 3) %*% t(fx$R)
    expect_equal(diagonal$fit$residual[[i]], frob2(S - A), tolerance = 1e-10)
    expect_equal(diagonal$fit$fitted_energy[[i]], sum(diag(A)),
      tolerance = 1e-10)
    # The diagonal objective in full coordinates, per section 4.1, equals the
    # compressed residual plus the constant outside the span.
    constant <- frob2(G - fx$Q %*% S %*% t(fx$Q))
    expect_equal(frob2(G - fx$T %*% diag(c_diag, 3) %*% t(fx$T)),
      constant + diagonal$fit$residual[[i]], tolerance = 1e-8)

    k <- length(columns)
    g_iso <- vapply(columns, function(cols) sum(g[cols]), numeric(1))
    M_iso <- matrix(0, k, k)
    for (a in seq_len(k)) for (b in seq_len(k)) {
      M_iso[a, b] <- sum(M[columns[[a]], columns[[b]]])
    }
    alpha <- isotropic$weights[i, ]
    expect_true(kkt_holds(alpha, M_iso, g_iso))
    C_iso <- diag(rep(alpha, times = lengths(columns)), 3)
    expect_equal(unname(isotropic$coefficients[, , i]), C_iso, tolerance = 1e-12)
    # Isotropic is a restriction of diagonal, so it fits no better.
    expect_gte(isotropic$fit$residual[[i]], diagonal$fit$residual[[i]] - 1e-8)
  }
  expect_identical(isotropic$df, list(fit = 2, free = 12))
  expect_identical(diagonal$df, list(fit = 3, free = 12))
  expect_error(model_geometry(fx$plan, fx$basis, structure = "diagonal",
    rank = 1), "no rank budget", class = "effect_input_error")
})

test_that("the shared fit is the latent layer's rank budget and the block fit nests between", {
  fx <- mg_fixture()
  shared <- model_geometry(fx$plan, fx$basis, structure = "shared", rank = 1)
  shared_full <- model_geometry(fx$plan, fx$basis, structure = "shared")
  block <- model_geometry(fx$plan, fx$basis, structure = "block")
  diagonal <- model_geometry(fx$plan, fx$basis, structure = "diagonal")
  R_inv <- solve(fx$R)

  # The lowered geometry, materialized directly, and its latent layer.
  lowered <- relation(fx$betas, extract = fx$basis, domain = fx$domain)
  latent <- latent_geometry(materialize_geometry(plan_geometry(lowered, fx$at,
    cross_partitions(lowered, independence = "independent"))), rank = 1)

  for (i in 1:2) {
    S <- t(fx$Q) %*% fx$G[[i]] %*% fx$Q
    A_star <- psd_rank(S, 1L)
    expect_equal(unname(shared$coefficients[, , i]),
      R_inv %*% A_star %*% t(R_inv), tolerance = 1e-10)
    expect_equal(shared$fit$residual[[i]], frob2(S - A_star), tolerance = 1e-10)
    expect_equal(shared$fit$fitted_energy[[i]], sum(diag(A_star)),
      tolerance = 1e-10)
    expect_equal(shared$fit$fitted_energy[[i]], sum(latent$spectrum[i, ]),
      tolerance = 1e-10)
    expect_equal(shared$fit$clipped_negative_mass[[i]],
      latent$clipped_negative_mass[[i]], tolerance = 1e-10)
    expect_equal(shared$fit$truncated_positive_mass[[i]],
      latent$truncated_positive_mass[[i]], tolerance = 1e-10)
    # Nesting of the feasible sets (oracle O8): shared(q) <= block <= diagonal.
    expect_lte(shared_full$fit$residual[[i]], block$fit$residual[[i]] + 1e-8)
    expect_lte(block$fit$residual[[i]], diagonal$fit$residual[[i]] + 1e-8)
    # The block C is block diagonal with PSD blocks of the declared ranks.
    C <- unname(block$coefficients[, , i])
    expect_equal(C[1, 2:3], c(0, 0), tolerance = 1e-12)
    expect_gte(min(eigen(C[2:3, 2:3], symmetric = TRUE)$values), -1e-10)
    expect_gte(C[1, 1], -1e-12)
  }
  expect_true(all(block$fit$converged))
  expect_true(all(block$fit$iterations >= 1L))
  expect_identical(shared$rank, list(requested = 1L, effective = 1L,
    clamped = FALSE))
  expect_identical(shared$df, list(fit = 3, free = 5))
  expect_identical(shared_full$df, list(fit = 6, free = 12))
  expect_identical(block$rank$effective, c(category = 1L, feature = 2L))
  expect_identical(block$df, list(fit = 4, free = 12))
  expect_false(shared$tests_supplied_geometries)
  expect_null(block$fit$clipped_negative_mass)

  # Rank budgets: clamped and recorded, named per model, or refused.
  clamped <- model_geometry(fx$plan, fx$basis, structure = "shared", rank = 9)
  expect_true(clamped$rank$clamped)
  expect_identical(clamped$rank$effective, 3L)
  named <- model_geometry(fx$plan, fx$basis, structure = "block",
    rank = c(feature = 1, category = 1))
  expect_identical(named$rank$effective, c(category = 1L, feature = 1L))
  # One rank-1 PSD block on one dimension (1 parameter) and one on two
  # dimensions (2 parameters).
  expect_identical(named$df$fit, 3)
  expect_error(model_geometry(fx$plan, fx$basis, structure = "block",
    rank = c(other = 1, feature = 1)), "must match", class = "effect_input_error")
  # A shared budget is a plain latent identity change.
  expect_false(identical(shared$receipt$scientific_plan_id,
    shared_full$receipt$scientific_plan_id))
})

test_that("saturation refuses the shared fit and labels the others", {
  conditions <- c("face", "body", "house", "tool")
  set.seed(9)
  n <- 4L
  p <- 6L
  labels <- c(1, 1, 0, 0)
  category <- outer(labels, labels, function(a, b) (a - b)^2)
  dimnames(category) <- list(conditions, conditions)
  shape <- matrix(rnorm(8), 4, 2, dimnames = list(conditions, NULL))
  saturated <- model_basis(list(category = category), "squared_euclidean",
    features = list(shape = shape), conditions = conditions)
  expect_true(saturated$saturated)
  betas <- lapply(1:2, function(run) {
    b <- matrix(rnorm(n * p), n, p)
    rownames(b) <- conditions
    b
  })
  names(betas) <- c("run1", "run2")
  domain <- abstract_domain(p, id = "model-geometry-saturated")
  relation <- relation(betas, effects = conditions, domain = domain)
  plan <- plan_geometry(relation, compile_frame(whole_brain(), domain),
    cross_partitions(relation, independence = "independent"))

  refusal <- catch_refusal(model_geometry(plan, saturated, structure = "shared"))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "model_geometry_test")
  expect_identical(refusal$reasons[[1L]], "model_span_saturated")
  admitted <- model_geometry(plan, saturated, structure = "isotropic")
  expect_true(admitted$near_saturated)
  expect_equal(admitted$span_fraction, 1)
  # On a full-span basis the orthogonal term is zero at every measurement.
  expect_lt(max(abs(admitted$orthogonal)), 1e-10)
  # Below the warning fraction the label is off.
  proper <- model_basis(list(category = category), "squared_euclidean",
    conditions = conditions)
  expect_false(model_geometry(plan, proper)$near_saturated)
  expect_true(model_geometry(plan, proper, warning_fraction = 1 / 3)$near_saturated)
  expect_error(model_geometry(plan, proper, warning_fraction = 2),
    "at most one", class = "effect_input_error")
})

test_that("a block fit over one model is the shared fit and is gated as one", {
  conditions <- c("face", "body", "house", "tool")
  set.seed(10)
  n <- 4L
  p <- 6L
  shape <- matrix(rnorm(12), 4, 3, dimnames = list(conditions, NULL))
  saturated <- model_basis(features = list(shape = shape),
    conditions = conditions)
  expect_true(saturated$saturated)
  betas <- lapply(1:2, function(run) {
    b <- matrix(rnorm(n * p), n, p)
    rownames(b) <- conditions
    b
  })
  names(betas) <- c("run1", "run2")
  domain <- abstract_domain(p, id = "model-geometry-one-block")
  relation <- relation(betas, effects = conditions, domain = domain)
  plan <- plan_geometry(relation, compile_frame(whole_brain(), domain),
    cross_partitions(relation, independence = "independent"))
  refusal <- catch_refusal(model_geometry(plan, saturated, structure = "block"))
  expect_identical(refusal$capability, "model_geometry_test")
  expect_identical(refusal$reasons[[1L]], "model_span_saturated")

  # On a proper single-model basis the block fit coincides with the shared
  # fit, says so, and converges in one sweep.
  proper <- model_basis(features = list(shape = shape[, 1:2]),
    conditions = conditions)
  block <- model_geometry(plan, proper, structure = "block", max_sweeps = 1)
  shared <- model_geometry(plan, proper, structure = "shared")
  expect_equal(block$coefficients, shared$coefficients, tolerance = 1e-10)
  expect_equal(block$fit$residual, shared$fit$residual, tolerance = 1e-10)
  expect_true(all(block$fit$converged))
  expect_match(block$hypothesis, "the block fit is the shared fit", fixed = TRUE)
  expect_false(block$tests_supplied_geometries)
  expect_error(model_geometry(plan, proper, tolerance = 1), "below one",
    class = "effect_input_error")
  expect_true("converged" %in% names(as.data.frame(block)))
})

test_that("the reader refuses what it cannot lower", {
  fx <- mg_fixture()
  # A basis over the wrong condition order.
  reversed <- model_basis(list(category = outer(c(1, 1, 1, 0, 0, 0),
    c(1, 1, 1, 0, 0, 0), function(a, b) (a - b)^2) |>
    `dimnames<-`(list(mg_conditions, mg_conditions))), "squared_euclidean",
    conditions = rev(mg_conditions))
  expect_error(model_geometry(fx$plan, reversed), "plan's order",
    class = "effect_input_error")
  # An already-lowered plan.
  lowered <- relation(fx$betas, extract = fx$basis, domain = fx$domain)
  lowered_plan <- plan_geometry(lowered, fx$at,
    cross_partitions(lowered, independence = "independent"))
  q_basis <- model_basis(features = list(f = matrix(rnorm(3 * 2), 3, 2,
    dimnames = list(paste0("mc", 1:3), NULL))))
  expect_error(model_geometry(lowered_plan, q_basis), "already a lowered plan",
    class = "effect_input_error")
  # A whitened schedule cannot be re-planned on the lowered axis.
  metric <- neural_metric(diag(seq(1, 2, length.out = fx$p)), fx$domain)
  whitened <- plan_geometry(fx$relation, fx$at, fx$over, metric = metric,
    composition = "whitened")
  skip_if(whitened$metric_schedule$kind %in%
    c("implicit_identity_before_frame", "fixed_metric_before_frame"))
  refusal <- catch_refusal(model_geometry(whitened, fx$basis))
  expect_identical(refusal$capability, "model_coordinate_lowering")
  expect_identical(refusal$reasons[[1L]], "metric_schedule_not_lowerable")
  # A fixed native metric is admitted and read under that metric.
  native <- plan_geometry(fx$relation, fx$at, fx$over, metric = metric)
  under_metric <- model_geometry(native, fx$basis)
  expect_false(isTRUE(all.equal(under_metric$addressable,
    model_geometry(fx$plan, fx$basis)$addressable)))
  expect_lt(under_metric$identity_error, 1e-12)
})

test_that("the reader prints its layers apart and coerces to one row per measurement", {
  fx <- mg_fixture()
  reading <- model_geometry(fx$plan, fx$basis)
  printed <- capture.output(print(reading))
  expect_true(any(grepl("isotropic", printed, fixed = TRUE)))
  expect_true(any(grepl("tests the truncated model geometries", printed,
    fixed = TRUE)))
  expect_true(any(grepl("signed tr(S_x)", printed, fixed = TRUE)))
  expect_true(any(grepl("latent descriptive layer; not for inference", printed,
    fixed = TRUE)))
  expect_true(any(grepl("split closes", printed, fixed = TRUE)))
  expect_match(format(reading), "isotropic fit over 2 models", fixed = TRUE)

  frame <- as.data.frame(reading)
  expect_identical(names(frame), c("measurement", "addressable", "orthogonal",
    "centered_total", "fit_residual", "fitted_energy", "converged",
    "weight_category", "weight_feature"))
  expect_equal(frame$addressable + frame$orthogonal, frame$centered_total,
    tolerance = 1e-10)
  shared <- model_geometry(fx$plan, fx$basis, structure = "shared")
  expect_identical(names(as.data.frame(shared)), c("measurement", "addressable",
    "orthogonal", "centered_total", "fit_residual", "fitted_energy",
    "converged"))
  expect_true(any(grepl("rank", capture.output(print(shared)), fixed = TRUE)))
})

test_that("a forged reading fails closed", {
  fx <- mg_fixture()
  reading <- model_geometry(fx$plan, fx$basis)
  forged <- reading
  forged$addressable[[1L]] <- forged$addressable[[1L]] + 1
  expect_error(crossform:::.validate_effect_model_geometry(forged),
    "does not close", class = "effect_contract_error")
  forged <- reading
  forged$structure <- "shared"
  expect_error(crossform:::.validate_effect_model_geometry(forged),
    "misstates", class = "effect_contract_error")
  forged <- reading
  forged$fit$fitted_energy[[1L]] <- -1
  expect_error(crossform:::.validate_effect_model_geometry(forged),
    "negative", class = "effect_contract_error")
  forged <- reading
  forged$weights[1, 1] <- -0.1
  forged$coefficients[1, 1, 1] <- -0.1
  expect_error(crossform:::.validate_effect_model_geometry(forged),
    "negative weight", class = "effect_contract_error")
  # Every reported value is under the signature, and the derivable fields
  # are re-derived, so an edit to any of them is a contradiction.
  forge <- function(edit, pattern = "signature is inconsistent") {
    forged <- reading
    forged <- edit(forged)
    expect_error(crossform:::.validate_effect_model_geometry(forged), pattern,
      class = "effect_contract_error")
  }
  forge(function(r) { r$weights[1, 1] <- r$weights[1, 1] + 1; r },
    "not the diagonal")
  forge(function(r) { r$df$fit <- 99; r }, "degrees of freedom")
  forge(function(r) { r$near_saturated <- TRUE; r }, "saturation label")
  forge(function(r) { r$span_fraction <- 0.99; r }, "saturation label")
  forge(function(r) { r$coefficients[1, 1, 2] <- r$coefficients[1, 1, 2] + 1; r },
    "not the diagonal")
  forge(function(r) { r$fit$residual[] <- 0; r })
  forge(function(r) { r$fit$converged[[1L]] <- FALSE; r })
  forge(function(r) { r$orthogonal_source <- "elsewhere"; r }, "misstates")
  forge(function(r) { r$model_sizes[[1L]] <- 5L; r }, "saturation label")
  forged <- reading
  forged$component <- "coherent"
  expect_error(crossform:::.validate_effect_model_geometry(forged),
    "signature is inconsistent", class = "effect_contract_error")
  forged <- reading
  forged$extra <- 1
  expect_error(crossform:::.validate_effect_model_geometry(forged),
    "canonical", class = "effect_input_error")
})

# The edge-disjoint cross-fit (ticket M6) --------------------------------------
#
# Contract section 7: the learned readout of the shared structure is held
# out over partition edges under the exclude-evaluation discipline, the same
# unit a learned neural metric is trained on. The first test reproduces the
# package's cross-fitted energy from the ORIGINAL plan's own executor on the
# restricted pairings, so the lowering and the cross-fit are checked against
# each other through different code paths; the Monte Carlo test is the
# contract's O7 law on the package route.

cross_fit_oracle <- function(fx, rank) {
  partitions <- fx$relation$partitions
  edges <- as.data.frame(fx$plan$pairing)
  form_on <- function(subset, i) {
    plan <- plan_geometry(fx$relation, fx$at,
      cross_partitions(subset, independence = "independent"))
    packed <- geometry_component(materialize_geometry(plan), "total")
    t(fx$Q) %*% crossform:::.unsvec_symmetric(packed[i, ], fx$n) %*% fx$Q
  }
  m <- length(fx$G)
  vapply(seq_len(m), function(i) {
    sum(vapply(seq_len(nrow(edges)), function(e) {
      evaluation <- c(edges$left[[e]], edges$right[[e]])
      S_train <- form_on(setdiff(partitions, evaluation), i)
      S_eval <- form_on(evaluation, i)
      V <- eigen(S_train, symmetric = TRUE)$vectors[, seq_len(rank), drop = FALSE]
      edges$weight[[e]] * sum(diag(t(V) %*% S_eval %*% V))
    }, numeric(1))) / sum(edges$weight)
  }, numeric(1))
}

test_that("the cross-fitted energy is learned on disjoint edges and reduced with the pairing's weights", {
  fx <- mg_fixture(runs = 4L)
  policy <- metric_training_policy("exclude_evaluation")
  crossed <- model_geometry(fx$plan, fx$basis, structure = "shared", rank = 1,
    training = policy)
  plain <- model_geometry(fx$plan, fx$basis, structure = "shared", rank = 1)

  expect_identical(crossed$training, "exclude_evaluation")
  expect_true("cross_fit" %in% crossed$layer$estimation)
  expect_identical(crossed$cross_fit$rank, 1L)
  expect_identical(nrow(crossed$cross_fit$edges), 6L)
  expect_identical(dim(crossed$cross_fit$by_edge), c(2L, 6L))
  expect_equal(sum(crossed$cross_fit$edges$weight), 1, tolerance = 1e-12)
  expect_true(all(vapply(seq_len(6L), function(e) {
    trained <- strsplit(crossed$cross_fit$edges$trained_on[[e]], ",")[[1L]]
    !any(c(crossed$cross_fit$edges$left[[e]],
      crossed$cross_fit$edges$right[[e]]) %in% trained) && length(trained) == 2L
  }, logical(1))))
  expect_equal(crossed$cross_fit$energy, cross_fit_oracle(fx, 1L),
    tolerance = 1e-10)
  # The plug-in fit is unchanged by the cross-fit and is a different number.
  expect_identical(crossed$fit, plain$fit)
  expect_false(isTRUE(all.equal(crossed$cross_fit$energy,
    crossed$fit$fitted_energy)))
  expect_false(identical(crossed$receipt$scientific_plan_id,
    plain$receipt$scientific_plan_id))
  expect_null(plain$cross_fit)
  expect_null(plain$training)

  frame <- as.data.frame(crossed)
  expect_true("cross_fit_energy" %in% names(frame))
  expect_equal(frame$cross_fit_energy, crossed$cross_fit$energy)
  printed <- capture.output(print(crossed))
  expect_true(any(grepl("readout learned on edges disjoint", printed, fixed = TRUE)))
  expect_true(any(grepl("estimation layer, not the plug-in fit", printed, fixed = TRUE)))

  # Sealed: the reduction and the policy are checked.
  forged <- crossed
  forged$cross_fit$energy[[1L]] <- forged$cross_fit$energy[[1L]] + 1
  expect_error(crossform:::.validate_effect_model_geometry(forged),
    "weighted reduction", class = "effect_contract_error")
  forged <- crossed
  forged$training <- NULL
  expect_error(crossform:::.validate_effect_model_geometry(forged),
    "canonical", class = "effect_input_error")
  forged <- crossed
  forged$cross_fit <- NULL
  expect_error(crossform:::.validate_effect_model_geometry(forged),
    "canonical", class = "effect_input_error")
})

test_that("the cross-fit takes its edges from the declared pairing, not from the relation", {
  fx <- mg_fixture(runs = 5L)
  policy <- metric_training_policy("exclude_evaluation")

  # A pairing over three of the five runs: the floor counts the pairing's
  # partitions, so the two runs the plan excludes never train a readout.
  three <- plan_geometry(fx$relation, fx$at,
    cross_partitions(c("run1", "run2", "run3"), independence = "independent"))
  refusal <- catch_refusal(model_geometry(three, fx$basis, structure = "shared",
    rank = 1, training = policy))
  expect_identical(refusal$reasons[[1L]], "insufficient_disjoint_edges")
  expect_match(refusal$reasons[[2L]], "partitions: 3", fixed = TRUE)

  # A star has no edge disjoint from any evaluation edge.
  star <- plan_geometry(fx$relation, fx$at, pairing(rep("run1", 4L),
    paste0("run", 2:5), independence = "independent"))
  refusal <- catch_refusal(model_geometry(star, fx$basis, structure = "shared",
    rank = 1, training = policy))
  expect_identical(refusal$reasons[[1L]], "insufficient_disjoint_edges")
  expect_match(refusal$reasons[[2L]], "no disjoint training edge", fixed = TRUE)

  # A self-pair plan declares no cross-partition edge at all.
  self_pairs <- plan_geometry(fx$relation, fx$at,
    pairing(fx$relation$partitions, fx$relation$partitions,
      self_pairs = "allow_biased", independence = "not_independent"))
  refusal <- catch_refusal(model_geometry(self_pairs, fx$basis,
    structure = "shared", rank = 1, training = policy))
  expect_identical(refusal$reasons[[1L]], "insufficient_disjoint_edges")
  expect_match(refusal$reasons[[2L]], "only self products", fixed = TRUE)

  # A five-cycle with unequal weights: each readout trains on the two cycle
  # edges disjoint from its evaluation edge, with those edges' weights, and
  # the reduction uses the pairing's weights. The oracle re-plans the
  # ORIGINAL relation on the same restricted pairings, so the lowering and
  # the cross-fit are checked through different code paths.
  cycle <- pairing(paste0("run", 1:5), paste0("run", c(2:5, 1)),
    weight = 1:5, independence = "independent")
  cycle_plan <- plan_geometry(fx$relation, fx$at, cycle)
  crossed <- model_geometry(cycle_plan, fx$basis, structure = "shared", rank = 1,
    training = policy)
  edges <- as.data.frame(cycle)
  expect_identical(nrow(crossed$cross_fit$edges), 5L)
  expect_equal(crossed$cross_fit$edges$weight, edges$weight / sum(edges$weight),
    tolerance = 1e-12)
  expect_identical(crossed$cross_fit$edges$training_edges, rep(2L, 5L))
  expect_identical(crossed$cross_fit$edges$trained_on[[1L]], "run3,run4,run5")
  expect_identical(crossed$cross_fit$independence, "independent")
  form_on <- function(over, i) {
    packed <- geometry_component(materialize_geometry(plan_geometry(
      fx$relation, fx$at, over)), "total")
    t(fx$Q) %*% crossform:::.unsvec_symmetric(packed[i, ], fx$n) %*% fx$Q
  }
  oracle <- vapply(seq_along(fx$G), function(i) {
    sum(vapply(seq_len(nrow(edges)), function(e) {
      evaluation <- c(edges$left[[e]], edges$right[[e]])
      disjoint <- !(edges$left %in% evaluation) & !(edges$right %in% evaluation)
      train <- edges[disjoint, ]
      S_train <- form_on(pairing(train$left, train$right, weight = train$weight,
        independence = "independent"), i)
      S_eval <- form_on(pairing(evaluation[[1L]], evaluation[[2L]],
        independence = "independent"), i)
      v <- eigen(S_train, symmetric = TRUE)$vectors[, 1L]
      edges$weight[[e]] * drop(t(v) %*% S_eval %*% v)
    }, numeric(1))) / sum(edges$weight)
  }, numeric(1))
  expect_equal(crossed$cross_fit$energy, oracle, tolerance = 1e-10)

  # Bookkeeping is sealed: a readout trained on its own edge, or a rank the
  # object does not carry, is a contradiction.
  forged <- crossed
  forged$cross_fit$edges$trained_on[[1L]] <- "run1,run2"
  expect_error(crossform:::.validate_effect_model_geometry(forged),
    "disjoint-edge discipline", class = "effect_contract_error")
  forged <- crossed
  forged$cross_fit$rank <- 2L
  expect_error(crossform:::.validate_effect_model_geometry(forged),
    "disjoint-edge discipline", class = "effect_contract_error")
})

test_that("the cross-fit refuses below four partitions and outside the shared structure", {
  three <- mg_fixture(runs = 3L)
  policy <- metric_training_policy("exclude_evaluation")
  refusal <- catch_refusal(model_geometry(three$plan, three$basis,
    structure = "shared", rank = 1, training = policy))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "cross_fitted_model_energy")
  expect_identical(refusal$namespace, "model_coordinate")
  expect_identical(refusal$reasons[[1L]], "insufficient_disjoint_edges")
  four <- mg_fixture(runs = 4L)
  expect_error(model_geometry(four$plan, four$basis, structure = "isotropic",
    training = policy), "shared structure only", class = "effect_input_error")
  reuse <- metric_training_policy("all_partitions_residual_orthogonality",
    justification = "test")
  expect_error(model_geometry(four$plan, four$basis, structure = "shared",
    training = reuse), "exclude_evaluation", class = "effect_input_error")
  expect_error(model_geometry(four$plan, four$basis, structure = "shared",
    training = "exclude_evaluation"), class = "effect_input_error")
})

test_that("under pure noise the cross-fitted energy is centered on zero and the plug-in energy is not", {
  skip_on_cran()
  set.seed(20260906)
  n <- 6L
  p <- 4L
  reps <- 30L
  labels <- c(1, 1, 1, 0, 0, 0)
  category <- outer(labels, labels, function(a, b) (a - b)^2)
  dimnames(category) <- list(mg_conditions, mg_conditions)
  features <- matrix(rnorm(n * 2L), n, 2L, dimnames = list(mg_conditions, NULL))
  basis <- model_basis(list(category = category), "squared_euclidean",
    features = list(feature = features), conditions = mg_conditions)
  domain <- abstract_domain(p, id = "model-geometry-null")
  at <- compile_frame(whole_brain(), domain)
  policy <- metric_training_policy("exclude_evaluation")
  draws <- t(vapply(seq_len(reps), function(r) {
    betas <- lapply(1:4, function(run) {
      b <- matrix(rnorm(n * p), n, p)
      rownames(b) <- mg_conditions
      b
    })
    names(betas) <- paste0("run", 1:4)
    relation <- relation(betas, effects = mg_conditions, domain = domain)
    plan <- plan_geometry(relation, at,
      cross_partitions(relation, independence = "independent"))
    reading <- model_geometry(plan, basis, structure = "shared", rank = 1,
      training = policy)
    c(addressable = reading$addressable, latent = reading$fit$fitted_energy,
      cross_fit = reading$cross_fit$energy)
  }, numeric(3)))
  means <- colMeans(draws)
  se <- apply(draws, 2, stats::sd) / sqrt(reps)
  expect_lt(abs(means[["addressable"]]), 4 * se[["addressable"]])
  expect_lt(abs(means[["cross_fit"]]), 4 * se[["cross_fit"]])
  expect_gt(means[["latent"]], 4 * se[["latent"]])
})
