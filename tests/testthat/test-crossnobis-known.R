crossnobis_direct_oracle <- function(relation, frame, over, contrast, K) {
  weights <- as.matrix(frame$weights)
  relation_values <- lapply(relation$partitions, function(partition) {
    relation_block(relation, partition, seq_len(relation$n_features))
  })
  names(relation_values) <- relation$partitions
  vapply(seq_len(nrow(weights)), function(node) {
    support <- which(weights[node, ] > 0)
    root <- sqrt(weights[node, support])
    local_K <- K[support, support, drop = FALSE] * tcrossprod(root)
    sum(vapply(seq_len(nrow(over)), function(edge) {
      left <- drop(contrast %*%
        relation_values[[over$left[[edge]]]][, support, drop = FALSE])
      right <- drop(contrast %*%
        relation_values[[over$right[[edge]]]][, support, drop = FALSE])
      over$weight[[edge]] * drop(left %*% local_K %*% right)
    }, numeric(1)))
  }, numeric(1))
}

crossnobis_signal_fixture <- function(seed = 8401, signal = 0.7,
                                      noise_scale = 1, features = 12L,
                                      observation_whitener = NULL) {
  set.seed(seed)
  observations <- 80L
  domain <- abstract_domain(
    features,
    coordinates = cbind(x = seq_len(features), y = 0, z = 0),
    id = paste0("crossnobis-domain-", seed),
    coordinate_units = "mm"
  )
  frame <- compile_frame(searchlights(1.01), domain)
  condition <- rep(c(-0.5, 0.5), each = observations / 2)
  design <- cbind(intercept = 1, condition = condition)
  effects <- rbind(
    baseline = c(1, 0),
    condition = c(0, 1)
  )
  covariance <- toeplitz(0.55^(0:(features - 1L))) * noise_scale
  factor <- chol(covariance)
  planted <- numeric(features)
  planted[4:7] <- signal * c(1, -0.75, 0.5, 1.25)
  raw <- lapply(seq_len(3L), function(run) {
    noise <- matrix(rnorm(observations * features), observations, features) %*%
      factor
    design %*% rbind(rep(0, features), planted) + noise
  })
  names(raw) <- paste0("run", seq_along(raw))
  fit <- lm_relation_fit(
    raw, design, effects, domain = domain,
    observation_whitener = observation_whitener
  )
  relation <- fit$relation
  over <- cross_partitions(relation, independence = "independent")
  metric <- noise_precision(
    solve(covariance), domain, covariance = covariance,
    provenance = list(source = "simulation_truth")
  )
  list(
    fit = fit,
    relation = relation,
    domain = domain,
    frame = frame,
    over = over,
    metric = metric,
    covariance = covariance,
    planted = planted,
    contrast = c(baseline = 0, condition = 1)
  )
}

test_that("known-metric crossnobis is the exported evidence pairing", {
  fixture <- crossnobis_signal_fixture()
  plan <- plan_geometry(
    fixture$relation, fixture$frame, fixture$over,
    metric = fixture$metric,
    compute = compute_policy(workspace_bytes = 64 * 1024^2)
  )
  observed <- crossnobis(plan, fixture$contrast)
  expected <- crossnobis_direct_oracle(
    fixture$relation, fixture$frame, fixture$over, fixture$contrast,
    fixture$metric$value
  )

  expect_s3_class(observed, "effect_crossnobis_view")
  expect_identical(observed$estimand,
    "crossvalidated_squared_mahalanobis_contrast")
  expect_equal(observed$values, expected, tolerance = 2e-11)
  expect_identical(observed$metric, fixture$metric$signature)
  expect_identical(
    observed$metadata$execution_plan$lowering,
    "support_streamed_metric_query_contraction"
  )
  expect_identical(observed$receipt$kernel_version,
    "support-streamed-metric-v1")
  expect_false(observed$metadata$diagnostics$total$pair_atoms_materialized)
  expect_false(observed$metadata$diagnostics$total$pair_frame_materialized)
})

test_that("raw-response and precomputed-effect crossnobis paths agree", {
  fixture <- crossnobis_signal_fixture(
    observation_whitener = diag(seq(0.8, 1.2, length.out = 80L))
  )
  precomputed <- relation(
    lapply(fixture$relation$partitions, function(partition) {
      relation_block(
        fixture$fit, partition, seq_len(fixture$relation$n_features)
      )
    }),
    effects = fixture$relation$effect_space,
    partitions = fixture$relation$partitions,
    domain = fixture$domain
  )
  raw_plan <- plan_geometry(
    fixture$fit$relation, fixture$frame, fixture$over,
    metric = fixture$metric
  )
  precomputed_plan <- plan_geometry(
    precomputed, fixture$frame, fixture$over,
    metric = fixture$metric
  )

  expect_equal(
    crossnobis(raw_plan, fixture$contrast)$values,
    crossnobis(precomputed_plan, fixture$contrast)$values,
    tolerance = 2e-13
  )
  expect_identical(
    fixture$fit$error_models$run1$observation_whitener$kind,
    "explicit"
  )
})

test_that("diagonal and support-streamed routes agree on point supports", {
  fixture <- crossnobis_signal_fixture(features = 12L)
  frame <- compile_frame(voxelwise(), fixture$domain)
  diagonal <- seq(0.7, 1.8, length.out = fixture$relation$n_features)
  diagonal_metric <- noise_precision(diag(diagonal), fixture$domain)
  dense_value <- diag(diagonal)
  dense_value[1, 2] <- dense_value[2, 1] <- 0.05
  dense_metric <- noise_precision(dense_value, fixture$domain)
  diagonal_plan <- plan_geometry(
    fixture$relation, frame, fixture$over, metric = diagonal_metric
  )
  dense_plan <- plan_geometry(
    fixture$relation, frame, fixture$over, metric = dense_metric
  )
  diagonal_result <- crossnobis(diagonal_plan, fixture$contrast)
  dense_result <- crossnobis(dense_plan, fixture$contrast)

  expect_identical(diagonal_plan$lowering, "additive_contraction")
  expect_identical(dense_plan$lowering,
    "support_streamed_pair_contraction")
  expect_identical(
    diagonal_result$metadata$execution_plan$lowering,
    "additive_metric_query_fused_contraction"
  )
  expect_identical(
    dense_result$metadata$execution_plan$lowering,
    "support_streamed_metric_query_contraction"
  )
  expect_equal(diagonal_result$values, dense_result$values,
    tolerance = 2e-13)
})

test_that("crossnobis retains signed negative finite estimates", {
  domain <- abstract_domain(3, id = "signed-crossnobis")
  first <- rbind(baseline = 0, condition = c(1, -2, 0.5))
  second <- rbind(baseline = 0, condition = -c(1, -2, 0.5))
  relation <- relation(list(run1 = first, run2 = second), domain = domain)
  frame <- compile_frame(whole_brain(), domain)
  plan <- plan_geometry(
    relation, frame,
    cross_partitions(relation, independence = "independent"),
    metric = noise_precision(diag(3), domain)
  )
  value <- crossnobis(plan, c(baseline = 0, condition = 1))$values

  expect_length(value, 1L)
  expect_true(is.finite(value))
  expect_lt(value, 0)
  expect_equal(unname(value), -sum(c(1, -2, 0.5)^2) / 3,
    tolerance = 2e-14)
})

test_that("known-metric evidence is invariant to consistent neural scaling", {
  set.seed(8402)
  features <- 5L
  domain <- abstract_domain(features, id = "crossnobis-gauge-original")
  transformed_domain <- abstract_domain(
    features, id = "crossnobis-gauge-transformed"
  )
  scale <- c(0.5, 2, 1.5, 0.8, 3)
  T <- diag(scale)
  raw <- matrix(rnorm(features^2), features)
  covariance <- crossprod(raw) + diag(0.5, features)
  K <- solve(covariance)
  values <- list(
    run1 = rbind(a = rnorm(features), b = rnorm(features)),
    run2 = rbind(a = rnorm(features), b = rnorm(features))
  )
  transformed_values <- lapply(values, function(value) value %*% T)
  original <- relation(values, domain = domain)
  transformed <- relation(transformed_values, domain = transformed_domain)
  original_plan <- plan_geometry(
    original, compile_frame(whole_brain(), domain),
    cross_partitions(original, independence = "independent"),
    metric = noise_precision(K, domain, covariance = covariance)
  )
  transformed_covariance <- T %*% covariance %*% T
  transformed_plan <- plan_geometry(
    transformed, compile_frame(whole_brain(), transformed_domain),
    cross_partitions(transformed, independence = "independent"),
    metric = noise_precision(
      solve(transformed_covariance), transformed_domain,
      covariance = transformed_covariance
    )
  )

  expect_equal(
    crossnobis(original_plan, c(a = 1, b = -1))$values,
    crossnobis(transformed_plan, c(a = 1, b = -1))$values,
    tolerance = 3e-12
  )
})

test_that("known-metric evidence is invariant to an identified permutation", {
  set.seed(8404)
  features <- 6L
  permutation <- c(4L, 1L, 6L, 2L, 5L, 3L)
  domain <- abstract_domain(features, id = "crossnobis-permutation-original")
  permuted_domain <- abstract_domain(
    features, feature_ids = domain$feature_ids[permutation],
    id = "crossnobis-permutation-changed"
  )
  raw <- matrix(rnorm(features^2), features)
  covariance <- crossprod(raw) + diag(0.4, features)
  K <- solve(covariance)
  values <- list(
    run1 = rbind(a = rnorm(features), b = rnorm(features)),
    run2 = rbind(a = rnorm(features), b = rnorm(features))
  )
  permuted_values <- lapply(values, function(value) {
    value[, permutation, drop = FALSE]
  })
  original <- relation(values, domain = domain)
  changed <- relation(permuted_values, domain = permuted_domain)
  original_plan <- plan_geometry(
    original, compile_frame(whole_brain(), domain),
    cross_partitions(original, independence = "independent"),
    metric = noise_precision(K, domain, covariance = covariance)
  )
  permuted_covariance <- covariance[permutation, permutation, drop = FALSE]
  changed_plan <- plan_geometry(
    changed, compile_frame(whole_brain(), permuted_domain),
    cross_partitions(changed, independence = "independent"),
    metric = noise_precision(
      K[permutation, permutation, drop = FALSE], permuted_domain,
      covariance = permuted_covariance
    )
  )

  expect_equal(
    crossnobis(original_plan, c(a = 1, b = -1))$values,
    crossnobis(changed_plan, c(a = 1, b = -1))$values,
    tolerance = 3e-12
  )
})

test_that("near-singular fixed precision remains explicit and oracle-correct", {
  fixture <- crossnobis_signal_fixture()
  eigenvalues <- c(1, rep(1e-7, fixture$relation$n_features - 1L))
  covariance <- diag(eigenvalues)
  metric <- noise_precision(
    diag(1 / eigenvalues), fixture$domain, covariance = covariance,
    tolerance = 1e-12
  )
  plan <- plan_geometry(
    fixture$relation, fixture$frame, fixture$over, metric = metric
  )
  observed <- crossnobis(plan, fixture$contrast)$values
  expected <- crossnobis_direct_oracle(
    fixture$relation, fixture$frame, fixture$over,
    fixture$contrast, metric$value
  )

  expect_true(all(is.finite(observed)))
  expect_equal(observed, expected, tolerance = 1e-7)
})

test_that("crossnobis refuses undeclared metrics and coordinate identity errors", {
  fixture <- crossnobis_signal_fixture()
  euclidean_plan <- plan_geometry(
    fixture$relation, fixture$frame, fixture$over
  )
  generic_plan <- plan_geometry(
    fixture$relation, fixture$frame, fixture$over,
    metric = neural_metric(diag(fixture$relation$n_features), fixture$domain)
  )
  other_domain <- abstract_domain(
    fixture$relation$n_features,
    coordinates = fixture$domain$coordinates,
    id = "wrong-crossnobis-domain",
    coordinate_units = "mm"
  )

  implicit <- catch_refusal(crossnobis(euclidean_plan, fixture$contrast))
  expect_s3_class(implicit, "effect_capability_refusal")
  expect_identical(implicit$capability, "declared_noise_metric")
  expect_identical(implicit$namespace, "geometry_views")
  expect_identical(implicit$reasons,
    "implicit_identity_metric_is_not_a_noise_model")
  expect_match(conditionMessage(implicit), "explicit noise-precision")

  undeclared <- catch_refusal(crossnobis(generic_plan, fixture$contrast))
  expect_s3_class(undeclared, "effect_capability_refusal")
  expect_identical(undeclared$capability, "declared_noise_metric")
  expect_identical(undeclared$namespace, "geometry_views")
  # The two refusals share a capability but not a reason: one plan has no
  # metric at all, the other has one that was never declared as noise.
  expect_identical(undeclared$reasons, "metric_role_is_not_noise_precision")
  expect_match(conditionMessage(undeclared), "noise_precision", fixed = TRUE)
  expect_error(
    plan_geometry(
      fixture$relation, fixture$frame, fixture$over,
      metric = noise_precision(
        diag(fixture$relation$n_features), other_domain
      )
    ),
    "share one exact domain"
  , class = "effect_contract_error")
  biased <- pairing("run1", "run1", self_pairs = "allow_biased",
    independence = "not_independent")
  biased_plan <- plan_geometry(
    fixture$relation, fixture$frame, biased, metric = fixture$metric
  )
  expect_error(crossnobis(biased_plan, fixture$contrast),
    "requires cross-partition edges declared independent",
    class = "effect_input_error")
})

test_that("known-metric Monte Carlo recovers null and planted targets", {
  set.seed(8403)
  replications <- 64L
  observations <- 36L
  features <- 4L
  domain <- abstract_domain(features, id = "crossnobis-monte-carlo")
  frame <- compile_frame(whole_brain(), domain)
  null_condition <- rep(c(-0.5, 0.5), each = observations / 2)
  signal_condition <- rep(c(-0.5, 0.5), times = observations / 2)
  design <- cbind(
    intercept = 1,
    null = null_condition,
    signal = signal_condition
  )
  effects <- rbind(
    baseline = c(1, 0, 0),
    null = c(0, 1, 0),
    signal = c(0, 0, 1)
  )
  covariance <- toeplitz(0.4^(0:(features - 1L)))
  factor <- chol(covariance)
  K <- solve(covariance)
  metric <- noise_precision(K, domain, covariance = covariance)
  contrasts <- list(
    null = c(baseline = 0, null = 1, signal = 0),
    signal = c(baseline = 0, null = 0, signal = 1)
  )
  signal <- c(0.55, -0.35, 0.2, 0.4)
  target <- drop(signal %*% K %*% signal) / features
  estimates <- matrix(NA_real_, replications, 2L,
    dimnames = list(NULL, c("null", "signal")))
  # Prespecified per-replication SDs from an independent 120-replication
  # reference run for this exact design. Derive expected Monte Carlo SE from
  # the live replication count so changing test duration cannot silently
  # change either the bias or variance criterion.
  expected_standard_deviations <- c(
    null = 0.0058 * sqrt(120),
    signal = 0.0119 * sqrt(120)
  )
  expected_standard_errors <- expected_standard_deviations /
    sqrt(replications)

  for (replication in seq_len(replications)) {
    run_noise <- lapply(seq_len(2L), function(run) {
      matrix(rnorm(observations * features), observations, features) %*%
        factor
    })
    names(run_noise) <- c("run1", "run2")
    raw <- lapply(run_noise, function(noise) {
      design %*% rbind(
        rep(0, features), rep(0, features), signal
      ) + noise
    })
    fit <- lm_relation_fit(raw, design, effects, domain = domain)
    plan <- plan_geometry(
      fit$relation, frame,
      cross_partitions(fit$relation, independence = "independent"),
      metric = metric
    )
    estimates[replication, ] <- vapply(contrasts, function(contrast) {
      unname(crossnobis(plan, contrast)$values)
    }, numeric(1))
  }
  means <- colMeans(estimates)
  standard_errors <- apply(estimates, 2L, stats::sd) / sqrt(replications)

  expect_lte(abs(means[["null"]]),
    4 * expected_standard_errors[["null"]])
  expect_lte(abs(means[["signal"]] - target),
    4 * expected_standard_errors[["signal"]])
  expect_lte(standard_errors[["null"]],
    1.5 * expected_standard_errors[["null"]])
  expect_lte(standard_errors[["signal"]],
    1.5 * expected_standard_errors[["signal"]])
  expect_true(any(estimates[, "null"] < 0))
})

test_that("crossnobis is the named total of the metric-carrying contrast", {
  set.seed(72901)
  domain <- abstract_domain(4L, id = "crossnobis-alias-domain")
  sources <- list(
    run1 = matrix(rnorm(12), 3, 4, dimnames = list(c("a", "b", "c"), NULL)),
    run2 = matrix(rnorm(12), 3, 4, dimnames = list(c("a", "b", "c"), NULL))
  )
  relation <- relation(
    sources, effects = effect_space(c("a", "b", "c")), domain = domain
  )
  precision <- diag(c(1, 2, 0.5, 1.5))
  plan <- plan_geometry(
    relation, compile_frame(whole_brain(), domain),
    cross_partitions(relation, independence = "independent"),
    metric = noise_precision(precision, domain, covariance = solve(precision))
  )
  weights <- c(a = 1, b = -1, c = 0)
  named <- crossnobis(plan, weights)
  family <- contrast_energy(plan, weights)
  expect_equal(unname(named$values), unname(family$total),
    tolerance = 1e-12)

  # The decomposition of the same estimand is available from
  # contrast_energy(). Each of the three components is checked against an
  # independent first-principles oracle rather than against each other: the
  # implementation computes `configuration` as `total - coherent`, so
  # asserting their sum would prove nothing about either one.
  frame_weights <- as.matrix(compile_frame(whole_brain(), domain)$weights)
  over <- cross_partitions(relation, independence = "independent")
  oracle <- geometry_contrast_oracle(
    weights = unname(weights),
    relation_values = sources,
    frame_weights = frame_weights,
    partition_edges = over,
    metric = precision
  )
  expect_equal(unname(family$total), oracle$total, tolerance = 1e-12)
  expect_equal(unname(family$coherent), oracle$coherent, tolerance = 1e-12)
  expect_equal(unname(family$configuration), oracle$configuration,
    tolerance = 1e-12)
  # The crossnobis alias is the same named total, so it inherits the check.
  expect_equal(unname(named$values), oracle$total, tolerance = 1e-12)

  # The oracle's metric-general coherent operator is a rank-one projection
  # along the Riesz representative of the frame-weighted mean. Check the two
  # properties that pin it down: it annihilates nothing else (configuration
  # kills that direction exactly), and it reduces to the identity-metric
  # formula the contract writes, `(B_L w)(B_R w)' / sum(w)`.
  support_weight <- frame_weights[1L, ]
  local_metric <- precision * tcrossprod(sqrt(support_weight))
  profile <- support_weight / sum(support_weight)
  representative <- drop(solve(local_metric, profile))
  coherent_metric <- local_metric %*%
    (tcrossprod(representative) /
      drop(crossprod(representative, local_metric %*% representative))) %*%
    local_metric
  expect_equal(
    drop((local_metric - coherent_metric) %*% representative),
    numeric(length(representative)), tolerance = 1e-12
  )
  expect_identical(qr(coherent_metric)$rank, 1L)

  euclidean <- contrast_energy(
    plan_geometry(relation, compile_frame(whole_brain(), domain), over),
    weights
  )
  contract_coherent <- drop(
    unname(weights) %*% Reduce(`+`, lapply(seq_len(nrow(over)), function(edge) {
      left <- sources[[over$left[[edge]]]] %*% support_weight
      right <- sources[[over$right[[edge]]]] %*% support_weight
      cross <- tcrossprod(drop(left), drop(right)) / sum(support_weight)
      over$weight[[edge]] * (cross + t(cross)) / 2
    })) %*% unname(weights)
  )
  expect_equal(unname(euclidean$coherent), contract_coherent,
    tolerance = 1e-12)
})
