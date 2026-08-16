# The fixed-metric sampling law under a NON-spherical whitened residual
# covariance. Every configuration here has Sigma_w != c I, which is where an
# isotropic surrogate for the signal term and the exact law part company.

sampling_nonspherical_cache <- new.env(parent = emptyenv())

sampling_nonspherical_memo <- function(key, value) {
  if (exists(key, envir = sampling_nonspherical_cache, inherits = FALSE)) {
    return(get(key, envir = sampling_nonspherical_cache, inherits = FALSE))
  }
  value <- force(value)
  assign(key, value, envir = sampling_nonspherical_cache)
  value
}

sampling_nonspherical_truth <- function() {
  sampling_nonspherical_memo("truth", {
    conditions <- 3L
    features <- 6L
    partitions <- 5L
    per_condition <- 10L
    labels <- paste0("condition", seq_len(conditions))
    design <- stats::model.matrix(
      ~ 0 + factor(rep(labels, each = per_condition))
    )
    colnames(design) <- labels
    set.seed(81380)
    list(
      conditions = conditions,
      features = features,
      partitions = partitions,
      labels = labels,
      design = design,
      contrasts = sampling_oracle_condition_contrasts(conditions),
      sigma_k = solve(crossprod(design)),
      # AR(1), rho = 0.75: strongly non-spherical, and not diagonal.
      sigma_r = toeplitz(0.75^(0:(features - 1L))),
      patterns = matrix(
        rnorm(conditions * features, sd = 0.6), conditions, features
      )
    )
  })
}

# One crossform relation fit drawn from that truth. `index` selects an
# independent dataset; index 1 is the reference fit used by the exact tests.
sampling_nonspherical_fit <- function(index = 1L) {
  truth <- sampling_nonspherical_truth()
  key <- paste0("fit:", index)
  sampling_nonspherical_memo(key, {
    set.seed(81390L + index)
    root <- chol(truth$sigma_r)
    domain <- abstract_domain(truth$features, id = "sampling-nonspherical")
    observations <- nrow(truth$design)
    sources <- stats::setNames(
      lapply(seq_len(truth$partitions), function(partition) {
        truth$design %*% truth$patterns +
          matrix(
            rnorm(observations * truth$features), observations, truth$features
          ) %*% root
      }), paste0("run", seq_len(truth$partitions))
    )
    effect_map <- diag(truth$conditions)
    dimnames(effect_map) <- list(truth$labels, truth$labels)
    fit <- lm_relation_fit(
      sources, truth$design, effect_map, effect_names = truth$labels,
      sampling_unit = "trial", domain = domain
    )
    list(fit = fit, domain = domain)
  })
}

# The four configurations. `metric` names the base neural metric; the frame
# normalization supplies the weights w, so the effective metric crossform
# whitens with is D(sqrt(w)) K D(sqrt(w)) and the distance normalization is
# nu = 1 (the frame carries the averaging, not an explicit divisor).
sampling_nonspherical_configurations <- function() {
  list(
    list(
      label = "no frame normalization, identity metric",
      normalization = "none", metric = "identity", spherical = FALSE
    ),
    list(
      label = "no frame normalization, exact precision metric",
      normalization = "none", metric = "precision", spherical = FALSE
    ),
    list(
      label = "local frame normalization, identity metric",
      normalization = "local", metric = "identity", spherical = FALSE
    ),
    # The one historical coincidence: local normalization plus the exact
    # precision metric makes Sigma_w = I / P, the single configuration in
    # which an isotropic surrogate for the signal term is also correct.
    list(
      label = "local frame normalization, exact precision metric",
      normalization = "local", metric = "precision", spherical = TRUE
    )
  )
}

sampling_nonspherical_setup <- function(configuration, index = 1L) {
  truth <- sampling_nonspherical_truth()
  built <- sampling_nonspherical_fit(index)
  metric <- if (identical(configuration$metric, "identity")) {
    noise_precision(
      diag(truth$features), built$domain, covariance = diag(truth$features),
      provenance = list(source = "nonspherical-identity")
    )
  } else {
    noise_precision(
      solve(truth$sigma_r), built$domain, covariance = truth$sigma_r,
      provenance = list(source = "nonspherical-true-precision")
    )
  }
  evidence <- plan_geometry(
    built$fit$relation,
    compile_frame(whole_brain(configuration$normalization), built$domain),
    cross_partitions(built$fit$relation, independence = "independent"),
    metric = metric
  )
  weight <- if (identical(configuration$normalization, "local")) {
    rep(1 / truth$features, truth$features)
  } else {
    rep(1, truth$features)
  }
  # Independently reconstructed effective metric and whitening root.
  effective <- metric$value * tcrossprod(sqrt(weight))
  root <- chol(effective)
  list(
    truth = truth, fit = built$fit, domain = built$domain, metric = metric,
    evidence = evidence, effective = effective, root = root,
    whitened_patterns = truth$patterns %*% t(root),
    whitened_residual = root %*% truth$sigma_r %*% t(root),
    xi = truth$contrasts %*% truth$sigma_k %*% t(truth$contrasts)
  )
}

# What crossform computed before 2026-08-15: the whitened residual covariance
# replaced by (tr(Sigma_w^2) / nu) I inside the signal term only.
sampling_nonspherical_isotropic_surrogate <- function(differences, xi, sigma_w,
                                                      partitions,
                                                      normalization) {
  noise_trace <- sum(sigma_w * sigma_w) / normalization^2
  4 * xi * (noise_trace * tcrossprod(differences) / normalization) /
    partitions +
    2 * xi^2 * noise_trace / (partitions * (partitions - 1))
}

sampling_nonspherical_draws <- function(configuration, replications = 12000L,
                                        seed = 81400L, null = FALSE) {
  key <- paste("draws", configuration$normalization, configuration$metric,
    replications, seed, null, sep = ":")
  sampling_nonspherical_memo(key, {
    setup <- sampling_nonspherical_setup(configuration)
    differences <- setup$truth$contrasts %*% setup$whitened_patterns
    if (isTRUE(null)) differences[] <- 0
    sampling_oracle_distance_draws(
      differences, setup$xi, setup$whitened_residual,
      setup$truth$partitions, normalization = 1, replications = replications,
      seed = seed
    )
  })
}

# Pre-registered criterion, stated once. For N independent replications the
# Monte Carlo standard error of entry (r, s) of the empirical covariance is
# sd((x_r - xbar_r)(x_s - xbar_s)) / sqrt(N), computed from the same draws.
# Every analytic entry must sit within 4 of those standard errors. The
# multiplier and the replication count are fixed here, so changing test
# duration cannot silently relax the criterion.
sampling_nonspherical_standard_errors <- function(draws) {
  centered <- sweep(draws, 2L, colMeans(draws))
  dimension <- ncol(draws)
  value <- matrix(0, dimension, dimension)
  for (row in seq_len(dimension)) {
    for (column in seq_len(dimension)) {
      products <- centered[, row] * centered[, column]
      value[row, column] <- stats::sd(products) / sqrt(nrow(draws))
    }
  }
  value
}

test_that("crossform's estimator is the whitened all-pairs functional", {
  for (configuration in sampling_nonspherical_configurations()) {
    setup <- sampling_nonspherical_setup(configuration)
    whitened <- lapply(setup$fit$relation$partitions, function(partition) {
      relation_block(
        setup$fit, partition, seq_len(setup$truth$features)
      ) %*% t(setup$root)
    })
    oracle <- rowMeans(sampling_oracle_partition_edge_distances(
      whitened, setup$truth$contrasts, normalization = 1
    ))
    package <- drop(rdm(setup$evidence)$values)

    expect_equal(unname(package), unname(oracle), tolerance = 1e-12,
      info = configuration$label)
  }
})

test_that("the analytic law matches Monte Carlo under non-spherical noise", {
  for (configuration in sampling_nonspherical_configurations()) {
    setup <- sampling_nonspherical_setup(configuration)
    plan <- crossform:::.compile_fixed_metric_rdm_sampling(
      setup$evidence, setup$fit
    )
    analytic <- crossform:::.sampling_covariance_materialize(
      crossform:::.sampling_covariance_from_components(
        plan, setup$truth$contrasts, setup$whitened_patterns,
        setup$truth$sigma_k, setup$whitened_residual, normalization = 1,
        labels = rownames(setup$truth$contrasts)
      )
    )
    draws <- sampling_nonspherical_draws(configuration)
    empirical <- stats::cov(draws)
    tolerance <- 4 * sampling_nonspherical_standard_errors(draws)

    expect_lt(max(abs(unname(analytic) - empirical) / tolerance), 1,
      label = paste("analytic vs Monte Carlo,", configuration$label)
    )
    # Monte Carlo also confirms the estimator is unbiased for the whitened
    # squared distances the law is a covariance of.
    expect_equal(
      unname(colMeans(draws)),
      unname(diag(tcrossprod(setup$truth$contrasts %*% setup$whitened_patterns))),
      tolerance = 0.05, info = configuration$label
    )
  }
})

test_that("the pre-2026-08-15 isotropic surrogate fails that same court", {
  results <- vapply(sampling_nonspherical_configurations(), function(configuration) {
    setup <- sampling_nonspherical_setup(configuration)
    differences <- setup$truth$contrasts %*% setup$whitened_patterns
    surrogate <- sampling_nonspherical_isotropic_surrogate(
      differences, setup$xi, setup$whitened_residual,
      setup$truth$partitions, normalization = 1
    )
    draws <- sampling_nonspherical_draws(configuration)
    empirical <- stats::cov(draws)
    tolerance <- 4 * sampling_nonspherical_standard_errors(draws)
    max(abs(surrogate - empirical) / tolerance)
  }, numeric(1))
  spherical <- vapply(sampling_nonspherical_configurations(), function(x) {
    isTRUE(x$spherical)
  }, logical(1))

  # Three configurations reject the surrogate outright; the fourth is the
  # Sigma_w = I / P coincidence under which it happens to be correct, which
  # is exactly why the defect survived so long.
  expect_true(all(results[!spherical] > 10))
  expect_true(all(results[spherical] < 1))
})

test_that("the null target is exact and the plug-in target is algebraic", {
  for (configuration in sampling_nonspherical_configurations()) {
    setup <- sampling_nonspherical_setup(configuration)
    estimates <- lapply(setup$fit$relation$partitions, function(partition) {
      relation_block(setup$fit, partition, seq_len(setup$truth$features))
    })
    partition_mean <- Reduce(`+`, estimates) / length(estimates)
    residual_df_total <- sum(vapply(
      setup$fit$relation$partitions,
      function(partition) residual_df(setup$fit, partition), integer(1)
    ))
    residual <- Reduce(`+`, lapply(
      setup$fit$relation$partitions, function(partition) {
        crossprod(residual_block(
          setup$fit, partition, seq_len(setup$truth$features)
        ))
      }
    )) / residual_df_total
    whitened_residual <- setup$root %*% residual %*% t(setup$root)
    differences <- setup$truth$contrasts %*% partition_mean %*% t(setup$root)

    null_law <- sampling_oracle_scalar_law(
      0 * differences, setup$xi, whitened_residual,
      setup$truth$partitions, normalization = 1
    )
    plugin_law <- sampling_oracle_scalar_law(
      differences, setup$xi, whitened_residual,
      setup$truth$partitions, normalization = 1
    )
    null_package <- sampling_covariance(
      rdm_sampling_covariance(setup$evidence, setup$fit, target = "null"),
      "materialize"
    )
    plugin_package <- sampling_covariance(
      rdm_sampling_covariance(setup$evidence, setup$fit, target = "plugin"),
      "materialize"
    )

    expect_equal(unname(null_package), null_law, tolerance = 1e-12,
      info = configuration$label)
    expect_equal(unname(plugin_package), plugin_law, tolerance = 1e-12,
      info = configuration$label)
  }
})

test_that("the null target reproduces the Monte Carlo null covariance", {
  for (configuration in sampling_nonspherical_configurations()) {
    setup <- sampling_nonspherical_setup(configuration)
    plan <- crossform:::.compile_fixed_metric_rdm_sampling(
      setup$evidence, setup$fit,
      target = crossform:::.sampling_target("null")
    )
    analytic <- crossform:::.sampling_covariance_materialize(
      crossform:::.sampling_covariance_from_components(
        plan, setup$truth$contrasts, 0 * setup$whitened_patterns,
        setup$truth$sigma_k, setup$whitened_residual, normalization = 1,
        labels = rownames(setup$truth$contrasts)
      )
    )
    draws <- sampling_nonspherical_draws(configuration, null = TRUE)
    empirical <- stats::cov(draws)
    tolerance <- 4 * sampling_nonspherical_standard_errors(draws)

    expect_lt(max(abs(unname(analytic) - empirical) / tolerance), 1,
      label = paste("null target vs Monte Carlo,", configuration$label)
    )
  }
})

test_that("the plug-in target is biased upward by a stated amount", {
  # Substituting the partition mean of the ESTIMATES for the signal inflates
  # the signal term, because for M partitions
  #
  #   E[mu-hat_r Sigma_w mu-hat_s'] = mu_r Sigma_w mu_s'
  #                                   + Xi_rs tr(Sigma_w^2) / M,
  #
  # and the pooled residual covariance inflates the noise trace by
  # E[tr(S-hat^2)] = tr(Sigma_w^2)(1 + 1/df) + tr(Sigma_w)^2 / df. Both terms
  # are exactly predictable, so the policy is disclosed rather than hidden.
  configuration <- sampling_nonspherical_configurations()[[1L]]
  setup <- sampling_nonspherical_setup(configuration)
  partitions <- setup$truth$partitions
  differences <- setup$truth$contrasts %*% setup$whitened_patterns
  sigma_w <- setup$whitened_residual
  residual_df_total <- sum(vapply(
    setup$fit$relation$partitions,
    function(partition) residual_df(setup$fit, partition), integer(1)
  ))
  noise_trace <- sum(sigma_w * sigma_w)
  expected_noise_trace <- noise_trace * (1 + 1 / residual_df_total) +
    sum(diag(sigma_w))^2 / residual_df_total
  exact <- sampling_oracle_scalar_law(
    differences, setup$xi, sigma_w, partitions, normalization = 1
  )
  predicted <- 4 * setup$xi *
    (differences %*% sigma_w %*% t(differences) +
      setup$xi * noise_trace / partitions) / partitions +
    2 * setup$xi^2 * expected_noise_trace /
      (partitions * (partitions - 1))

  replications <- 20L
  observed <- vapply(seq_len(replications), function(index) {
    replicate_setup <- sampling_nonspherical_setup(configuration, index)
    diag(sampling_covariance(
      rdm_sampling_covariance(
        replicate_setup$evidence, replicate_setup$fit, target = "plugin"
      ),
      "materialize"
    ))
  }, numeric(nrow(setup$truth$contrasts)))
  observed_mean <- rowMeans(observed)
  standard_error <- apply(observed, 1L, stats::sd) / sqrt(replications)

  # The inflation is real, not a rounding artifact.
  expect_true(all(diag(predicted) > diag(exact) * 1.05))
  # And it is exactly the stated amount: 4 standard errors of the mean over
  # 20 independent datasets, with the replication count fixed here.
  expect_lt(
    max(abs(observed_mean - diag(predicted)) / (4 * standard_error)), 1
  )
})
