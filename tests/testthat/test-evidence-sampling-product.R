sampling_product_fixture <- function(conditions = 4L, features = 6L,
                                     partitions = 4L) {
  set.seed(81331)
  observations <- 64L
  condition <- factor(rep(seq_len(conditions), length.out = observations))
  design <- stats::model.matrix(~ 0 + condition)
  colnames(design) <- paste0("condition", seq_len(conditions))
  effects <- diag(conditions)
  rownames(effects) <- colnames(design)
  covariance <- toeplitz(0.35^(0:(features - 1L)))
  factor <- chol(covariance)
  signal <- matrix(rnorm(conditions * features, sd = 0.25),
    conditions, features)
  domain <- abstract_domain(
    features, coordinates = cbind(x = seq_len(features), y = 0, z = 0),
    id = "sampling-product-domain", coordinate_units = "mm"
  )
  sources <- stats::setNames(lapply(seq_len(partitions), function(run) {
    design %*% signal +
      matrix(rnorm(observations * features), observations, features) %*%
        factor
  }), paste0("run", seq_len(partitions)))
  fit <- lm_relation_fit(
    sources, design, effects, sampling_unit = "trial", domain = domain
  )
  evidence <- plan_geometry(
    fit$relation, compile_frame(whole_brain(), domain),
    cross_partitions(fit$relation, independence = "independent"),
    metric = noise_precision(
      solve(covariance), domain, covariance = covariance,
      provenance = list(source = "simulation_truth")
    )
  )
  list(
    fit = fit, evidence = evidence, domain = domain, covariance = covariance,
    signal = signal
  )
}

test_that("relation-fit product path supplies exact fixed-metric covariance", {
  fixture <- sampling_product_fixture()
  plan <- crossform:::.compile_fixed_metric_rdm_sampling(
    fixture$evidence, fixture$fit,
    target = crossform:::.sampling_target("null")
  )
  covariance <- crossform:::.fixed_metric_rdm_sampling_covariance(plan)
  observed <- crossform:::.execute_fixed_metric_rdm_sampling(plan)

  expect_s3_class(covariance, "effect_sampling_covariance")
  expect_identical(names(observed), covariance$labels)
  expect_true(all(is.finite(observed)))
  expect_true(all(observed >= 0))
  expect_identical(covariance$source$local_only, TRUE)
  expect_identical(covariance$source$specialization,
    "fixed_metric_equal_partition_rdm")
})

test_that("relation-fit product agrees with an independent direct oracle", {
  fixture <- sampling_product_fixture()
  # Two metrics. The exact precision metric makes Sigma_R spherical, which is
  # the historically flattering case; the identity metric leaves the AR-like
  # toeplitz residual covariance anisotropic, so the signal term's dependence
  # on Sigma_R is identifiable here and not only in the oracle.
  metrics <- list(
    exact_precision = solve(fixture$covariance),
    identity = diag(nrow(fixture$covariance))
  )
  for (name in names(metrics)) {
    evidence <- plan_geometry(
      fixture$fit$relation,
      compile_frame(whole_brain(), fixture$domain),
      cross_partitions(fixture$fit$relation, independence = "independent"),
      metric = noise_precision(
        metrics[[name]], fixture$domain,
        covariance = solve(metrics[[name]]),
        provenance = list(source = paste0("simulation_", name))
      )
    )
    covariance <- rdm_sampling_covariance(
      evidence, fixture$fit, target = "plugin", at = 1L
    )
    observed <- sampling_covariance(covariance, "materialize")
    estimates <- lapply(fixture$fit$relation$partitions, function(partition) {
      relation_block(fixture$fit, partition,
        seq_len(fixture$fit$relation$n_features))
    })
    mean_patterns <- Reduce(`+`, estimates) / length(estimates)
    # whole_brain(normalization = "local") composes K with weights 1 / P.
    # Keeping Sigma_R in the base whitened coordinates and letting Eq. 13
    # carry its explicit P normalization is the same law, so the two must
    # agree exactly rather than coincidentally.
    effective_root <- chol(metrics[[name]])
    whitened_signal <- mean_patterns %*% t(effective_root)
    residual_products <- lapply(fixture$fit$relation$partitions,
      function(partition) {
        crossprod(residual_block(
          fixture$fit, partition,
          seq_len(fixture$fit$relation$n_features)
        ))
      })
    total_df <- sum(vapply(
      fixture$fit$relation$partitions, function(partition) {
        residual_df(fixture$fit, partition)
      }, integer(1)
    ))
    residual <- Reduce(`+`, residual_products) / total_df
    residual_whitened <- effective_root %*% residual %*% t(effective_root)
    sigma_k <- effect_covariance(fixture$fit, 1L)
    components <- sampling_oracle_components(
      whitened_signal, sigma_k, residual_whitened
    )
    # `residual_whitened` is the pooled plug-in with `total_df` degrees of
    # freedom, so the oracle corrects the quadratic noise term exactly as the
    # package does; the signal term is linear and is left alone.
    expected <- sampling_oracle_eq13(
      components$differences, components$xi, residual_whitened,
      length(fixture$fit$relation$partitions),
      residual_df = total_df
    )$covariance

    expect_equal(unname(observed), unname(expected), tolerance = 4e-12,
      info = name)
    expect_identical(dim(observed), dim(expected))
  }
})

test_that("beta-only and learned-metric product paths refuse honestly", {
  fixture <- sampling_product_fixture()
  beta_only <- catch_refusal(
    crossform:::.compile_fixed_metric_rdm_sampling(
      fixture$evidence, fixture$fit$relation
    )
  )
  expect_s3_class(beta_only, "effect_capability_refusal")
  expect_identical(beta_only$capability, "sampling_covariance")
  expect_identical(beta_only$namespace, "evidence_sampling")
  expect_true("missing_error_channel" %in% beta_only$reasons)
  expect_match(conditionMessage(beta_only),
    "lm_relation_fit.*beta matrices alone")

  learned <- plan_crossnobis(
    fixture$fit,
    compile_frame(searchlights(100), fixture$domain),
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
  expect_s3_class(learned_refusal, "effect_capability_refusal")
  expect_identical(learned_refusal$capability, "fixed_metric_sampling_law")
  expect_identical(learned_refusal$namespace, "evidence_sampling")
  expect_identical(learned_refusal$reasons, "learned_metric_law_not_admitted")
  expect_match(conditionMessage(learned_refusal), "metric was learned")
})

test_that("overlapping supports reuse one exact sparse residual statistic", {
  fixture <- sampling_product_fixture(features = 9L)
  frame <- compile_frame(searchlights(2.01), fixture$domain)
  evidence <- plan_geometry(
    fixture$fit$relation, frame,
    cross_partitions(fixture$fit$relation, independence = "independent"),
    metric = noise_precision(
      solve(fixture$covariance), fixture$domain,
      covariance = fixture$covariance,
      provenance = list(source = "simulation_truth")
    )
  )
  plan <- crossform:::.compile_fixed_metric_rdm_sampling(
    evidence, fixture$fit,
    target = crossform:::.sampling_target("null")
  )
  resources <- crossform:::.fixed_metric_rdm_sampling_resources(
    plan, strategy = "shared_pair_statistics"
  )

  expect_s3_class(resources$residual_statistics,
    "effect_residual_pair_statistics")
  for (node in c(1L, 5L, 9L)) {
    sparse <- crossform:::.execute_fixed_metric_rdm_sampling(
      plan, node, resources = resources
    )
    node_value <- crossform:::.frame_metric_node_accessor(frame)(node)
    direct_residual <- Reduce(`+`, lapply(
      fixture$fit$relation$partitions, function(partition) {
        crossprod(residual_block(
          fixture$fit, partition, node_value$support_positions
        ))
      }
    )) / sum(vapply(fixture$fit$relation$partitions, function(partition) {
      residual_df(fixture$fit, partition)
    }, integer(1)))
    gathered <- crossform:::.sampling_node_residual_covariance(
      plan, node_value, resources
    )

    expect_equal(gathered, direct_residual, tolerance = 3e-13)
    expect_true(all(is.finite(sparse)))
  }
})

test_that("one-node sampling queries do not compile whole-frame residual state", {
  fixture <- sampling_product_fixture(features = 9L)
  frame <- compile_frame(searchlights(2.01), fixture$domain)
  evidence <- plan_geometry(
    fixture$fit$relation, frame,
    cross_partitions(fixture$fit$relation, independence = "independent"),
    metric = noise_precision(
      solve(fixture$covariance), fixture$domain,
      covariance = fixture$covariance,
      provenance = list(source = "simulation_truth")
    )
  )
  plan <- crossform:::.compile_fixed_metric_rdm_sampling(
    evidence, fixture$fit,
    target = crossform:::.sampling_target("null")
  )
  resources <- crossform:::.fixed_metric_rdm_sampling_resources(plan)

  expect_identical(resources$strategy, "node_local")
  expect_null(resources$residual_statistics)
  expect_true(all(is.finite(
    crossform:::.execute_fixed_metric_rdm_sampling(
      plan, node = 5L, resources = resources
    )
  )))
})

test_that("batched sampling covariance matches repeated single-node calls", {
  fixture <- sampling_product_fixture(features = 9L)
  frame <- compile_frame(searchlights(2.01), fixture$domain)
  evidence <- plan_geometry(
    fixture$fit$relation, frame,
    cross_partitions(fixture$fit$relation, independence = "independent"),
    metric = noise_precision(
      solve(fixture$covariance), fixture$domain,
      covariance = fixture$covariance,
      provenance = list(source = "simulation_truth")
    )
  )
  nodes <- c(1L, 4L, 7L, 9L)
  query_vector <- c(1, -1, 0, 0, 0, 0)
  selected <- cbind(1L, 2L)
  transport <- rbind(first = query_vector, second = rev(query_vector))

  repeated <- lapply(nodes, function(node) {
    rdm_sampling_covariance(
      evidence, fixture$fit, target = "null", at = node
    )
  })
  batched <- rdm_sampling_covariance(
    evidence, fixture$fit, target = "null", at = nodes
  )
  explicit_local <- rdm_sampling_covariance(
    evidence, fixture$fit, target = "null", at = nodes,
    residual_strategy = "node_local"
  )

  expect_s3_class(batched, "effect_sampling_covariance_batch")
  expect_identical(attr(batched, "basis"), "rdm")
  expect_length(batched, length(nodes))
  expect_identical(batched[[1L]]$plan$scientific_plan_id,
    repeated[[1L]]$plan$scientific_plan_id)
  expect_identical(batched[[1L]]$source$execution$route, "batched")
  expect_identical(
    batched[[1L]]$source$execution$residual_strategy,
    "shared_pair_statistics"
  )
  expect_true(batched[[1L]]$source$execution$shared_residual_statistics)
  expect_null(repeated[[1L]]$source$execution)
  expect_identical(
    explicit_local[[1L]]$source$execution$residual_strategy,
    "node_local"
  )
  expect_false(explicit_local[[1L]]$source$execution$shared_residual_statistics)

  for (index in seq_along(nodes)) {
    expect_equal(
      sampling_covariance(batched[[index]]),
      sampling_covariance(repeated[[index]]),
      tolerance = 1e-12
    )
    expect_equal(
      sampling_covariance(
        batched[[index]], "selected_entries", query = selected
      ),
      sampling_covariance(
        repeated[[index]], "selected_entries", query = selected
      ),
      tolerance = 1e-12
    )
    expect_equal(
      sampling_covariance(
        batched[[index]], "quadratic_form", query = query_vector
      ),
      sampling_covariance(
        repeated[[index]], "quadratic_form", query = query_vector
      ),
      tolerance = 1e-12
    )
    expect_equal(
      sampling_covariance(
        batched[[index]], "transport", query = transport
      ),
      sampling_covariance(
        repeated[[index]], "transport", query = transport
      ),
      tolerance = 1e-12
    )
  }
  expect_equal(
    sampling_covariance(batched, "quadratic_form", query = query_vector),
    lapply(repeated, sampling_covariance, "quadratic_form",
      query = query_vector),
    tolerance = 1e-12
  )
})

test_that("public RDM sampling covariance is explicit and exactly queryable", {
  fixture <- sampling_product_fixture()
  covariance <- rdm_sampling_covariance(
    fixture$evidence, fixture$fit, target = "plugin", at = 1L
  )
  dense <- sampling_covariance(covariance, "materialize")
  vector <- seq_len(nrow(dense)) / 10
  map <- rbind(first = vector, second = rev(vector))

  expect_s3_class(covariance, "effect_sampling_covariance")
  expect_identical(covariance$basis, "rdm")
  expect_match(format(covariance), "factorized")
  expect_equal(
    sampling_covariance(covariance),
    stats::setNames(diag(dense), covariance$labels),
    tolerance = 3e-13
  )
  expect_equal(
    sampling_covariance(
      covariance, "quadratic_form", query = vector
    ),
    drop(crossprod(vector, dense %*% vector)), tolerance = 3e-13
  )
  expect_equal(
    sampling_covariance(covariance, "transport", query = map),
    map %*% dense %*% t(map), tolerance = 4e-13
  )
  target_refusal <- catch_refusal(
    rdm_sampling_covariance(fixture$evidence, fixture$fit)
  )
  expect_s3_class(target_refusal, "effect_capability_refusal")
  expect_identical(target_refusal$capability, "calibration_target_declared")
  expect_identical(target_refusal$namespace, "evidence_sampling")
  expect_identical(target_refusal$reasons, "calibration_target_not_declared")
  expect_match(conditionMessage(target_refusal),
    "signal-dependent covariance target is not inferred")

  beta_only <- catch_refusal(
    rdm_sampling_covariance(
      fixture$evidence, fixture$fit$relation, target = "null", at = 1L
    )
  )
  expect_s3_class(beta_only, "effect_capability_refusal")
  expect_identical(beta_only$capability, "sampling_covariance")
  expect_identical(beta_only$namespace, "evidence_sampling")
  expect_true("missing_error_channel" %in% beta_only$reasons)
  expect_match(conditionMessage(beta_only),
    "lm_relation_fit.*beta matrices alone")
})

# Rendering ------------------------------------------------------------------
#
# The RDM basis of a sampling covariance is a field, not a subclass, so what a
# reader sees for that basis is pinned here rather than left to whichever
# method happens to dispatch. The content digest is masked: it hashes
# double-precision row factors, so it is a platform fact rather than a
# rendering choice.

mask_sampling_digest <- function(lines) {
  sub("(signature:\\s+sha256:)[0-9a-f]+", "\\1<digest>", lines)
}

test_that("an RDM sampling covariance and its batch render for a reader", {
  fixture <- sampling_product_fixture(features = 9L)
  frame <- compile_frame(searchlights(2.01), fixture$domain)
  evidence <- plan_geometry(
    fixture$fit$relation, frame,
    cross_partitions(fixture$fit$relation, independence = "independent"),
    metric = noise_precision(
      solve(fixture$covariance), fixture$domain,
      covariance = fixture$covariance,
      provenance = list(source = "simulation_truth")
    )
  )
  single <- rdm_sampling_covariance(
    evidence, fixture$fit, target = "null", at = 1L
  )
  batched <- rdm_sampling_covariance(
    evidence, fixture$fit, target = "null", at = c(1L, 4L, 7L)
  )

  expect_snapshot(print(single), transform = mask_sampling_digest)
  expect_snapshot(format(single))
  expect_snapshot(print(batched), transform = mask_sampling_digest)
})
