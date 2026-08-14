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
    cross_partitions(fit$relation),
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
  plan <- effectagram:::.compile_fixed_metric_rdm_sampling(
    fixture$evidence, fixture$fit,
    target = effectagram:::.sampling_target("null")
  )
  covariance <- effectagram:::.fixed_metric_rdm_sampling_covariance(plan)
  observed <- effectagram:::.execute_fixed_metric_rdm_sampling(plan)

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
  covariance <- rdm_sampling_covariance(
    fixture$evidence, fixture$fit, target = "plugin"
  )
  observed <- sampling_covariance(covariance, "materialize")
  effects <- fixture$fit$relation$effects
  contrasts <- sampling_oracle_condition_contrasts(length(effects))
  estimates <- lapply(fixture$fit$relation$partitions, function(partition) {
    relation_block(fixture$fit, partition,
      seq_len(fixture$fit$relation$n_features))
  })
  mean_patterns <- Reduce(`+`, estimates) / length(estimates)
  metric_root <- chol(solve(fixture$covariance))
  # whole_brain(normalization = "local") composes K with weights 1 / P.
  # Here Sigma_R is kept in the base whitened coordinates and Eq. 13 retains
  # its explicit P normalization, matching the product specialization.
  effective_root <- metric_root
  whitened_signal <- mean_patterns %*% t(effective_root)
  residual_products <- lapply(fixture$fit$relation$partitions,
    function(partition) {
      crossprod(residual_block(
        fixture$fit, partition,
        seq_len(fixture$fit$relation$n_features)
      ))
    })
  total_df <- sum(vapply(fixture$fit$relation$partitions, function(partition) {
    residual_df(fixture$fit, partition)
  }, integer(1)))
  residual <- Reduce(`+`, residual_products) / total_df
  residual_whitened <- effective_root %*% residual %*% t(effective_root)
  sigma_k <- effect_covariance(fixture$fit, 1L)
  components <- sampling_oracle_components(
    whitened_signal, sigma_k, residual_whitened
  )
  expected <- sampling_oracle_eq13(
    components$delta, components$xi, residual_whitened,
    length(fixture$fit$relation$partitions)
  )$covariance

  expect_equal(unname(observed), unname(expected), tolerance = 4e-12)
  expect_identical(dim(observed), dim(expected))
})

test_that("beta-only and learned-metric product paths refuse honestly", {
  fixture <- sampling_product_fixture()
  expect_error(
    effectagram:::.compile_fixed_metric_rdm_sampling(
      fixture$evidence, fixture$fit$relation
    ),
    "lm_relation_fit.*beta matrices alone"
  )
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
  expect_error(
    rdm_sampling_covariance(learned, fixture$fit, target = "null"),
    "metric was learned"
  )
})

test_that("overlapping supports reuse one exact sparse residual statistic", {
  fixture <- sampling_product_fixture(features = 9L)
  frame <- compile_frame(searchlights(2.01), fixture$domain)
  evidence <- plan_geometry(
    fixture$fit$relation, frame, cross_partitions(fixture$fit$relation),
    metric = noise_precision(
      solve(fixture$covariance), fixture$domain,
      covariance = fixture$covariance,
      provenance = list(source = "simulation_truth")
    )
  )
  plan <- effectagram:::.compile_fixed_metric_rdm_sampling(
    evidence, fixture$fit,
    target = effectagram:::.sampling_target("null")
  )
  resources <- effectagram:::.fixed_metric_rdm_sampling_resources(
    plan, strategy = "shared_pair_statistics"
  )

  expect_s3_class(resources$residual_statistics,
    "effect_residual_pair_statistics")
  for (node in c(1L, 5L, 9L)) {
    sparse <- effectagram:::.execute_fixed_metric_rdm_sampling(
      plan, node, resources = resources
    )
    node_value <- effectagram:::.frame_metric_node_accessor(frame)(node)
    direct_residual <- Reduce(`+`, lapply(
      fixture$fit$relation$partitions, function(partition) {
        crossprod(residual_block(
          fixture$fit, partition, node_value$support_positions
        ))
      }
    )) / sum(vapply(fixture$fit$relation$partitions, function(partition) {
      residual_df(fixture$fit, partition)
    }, integer(1)))
    gathered <- effectagram:::.sampling_node_residual_covariance(
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
    fixture$fit$relation, frame, cross_partitions(fixture$fit$relation),
    metric = noise_precision(
      solve(fixture$covariance), fixture$domain,
      covariance = fixture$covariance,
      provenance = list(source = "simulation_truth")
    )
  )
  plan <- effectagram:::.compile_fixed_metric_rdm_sampling(
    evidence, fixture$fit,
    target = effectagram:::.sampling_target("null")
  )
  resources <- effectagram:::.fixed_metric_rdm_sampling_resources(plan)

  expect_identical(resources$strategy, "node_local")
  expect_null(resources$residual_statistics)
  expect_true(all(is.finite(
    effectagram:::.execute_fixed_metric_rdm_sampling(
      plan, node = 5L, resources = resources
    )
  )))
})

test_that("public RDM sampling covariance is explicit and exactly queryable", {
  fixture <- sampling_product_fixture()
  covariance <- rdm_sampling_covariance(
    fixture$evidence, fixture$fit, target = "plugin"
  )
  dense <- sampling_covariance(covariance, "materialize")
  vector <- seq_len(nrow(dense)) / 10
  map <- rbind(first = vector, second = rev(vector))

  expect_s3_class(covariance, "effect_rdm_sampling_covariance")
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
  expect_error(
    rdm_sampling_covariance(fixture$evidence, fixture$fit),
    "signal-dependent covariance target is not inferred"
  )
  expect_error(
    rdm_sampling_covariance(
      fixture$evidence, fixture$fit$relation, target = "null"
    ),
    "lm_relation_fit.*beta matrices alone"
  )
})
