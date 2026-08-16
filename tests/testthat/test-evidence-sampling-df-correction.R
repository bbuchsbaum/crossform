# The residual covariance the product path uses is a PLUG-IN with nu degrees
# of freedom, and it enters the sampling law's signal-independent term through
# the quadratic functional tr(Sigma_w^2). Taking tr(S_w^2) at face value
# overstates that functional, and therefore the reported standard error, by
# sqrt(1 + (1 + P_eff) / nu). This file is the end-to-end court for the
# correction: it refits `lm_relation_fit()` per replication, so the residual
# covariance is estimated exactly as a user's would be, and compares the
# empirical spread of the distance estimator against the standard error the
# package reports.
#
# Pre-registered criterion, stated once and fixed here:
#
#   P = 120 features, nu = 168 residual df (6 runs x (32 - 4)), M = 6
#   partitions, whole_brain("local") frame, identity metric, null truth.
#   N = 400 replications, so the Monte Carlo standard error of a ratio of
#   standard deviations is 1 / sqrt(2 (N - 1)) = 0.0354. The band is
#   1 +/- 4 x 0.0354 = [0.8584, 1.1416].
#
#   * empirical sd / mean corrected SE must fall INSIDE that band;
#   * empirical sd / mean UNCORRECTED SE must fall OUTSIDE it, low.
#
# The uncorrected ratio is predicted to be 1 / sqrt(1 + (1 + P) / nu) = 0.762
# under isotropic residual noise, where P_eff = P. Observed on 2026-08-16:
# corrected 0.972 to 1.025 over the six distances, uncorrected 0.741 to 0.781.
# Neither the band, the replication count, nor the seed may be widened to make
# a run pass.

sampling_df_correction_cache <- new.env(parent = emptyenv())

sampling_df_correction_memo <- function(key, value) {
  if (exists(key, envir = sampling_df_correction_cache, inherits = FALSE)) {
    return(get(key, envir = sampling_df_correction_cache, inherits = FALSE))
  }
  value <- force(value)
  assign(key, value, envir = sampling_df_correction_cache)
  value
}

sampling_df_correction_experiment <- function(features = 120L,
                                              conditions = 4L,
                                              per_condition = 8L,
                                              partitions = 6L,
                                              replications = 400L,
                                              seed = 20260816L) {
  key <- paste("df-correction", features, conditions, per_condition,
    partitions, replications, seed, sep = ":")
  sampling_df_correction_memo(key, {
    labels <- letters[seq_len(conditions)]
    set.seed(seed)
    order <- sample(rep(labels, each = per_condition))
    design <- stats::model.matrix(~ 0 + factor(order, levels = labels))
    colnames(design) <- labels
    effect_map <- diag(conditions)
    dimnames(effect_map) <- list(labels, labels)
    contrasts <- sampling_oracle_condition_contrasts(conditions)
    domain <- abstract_domain(features, id = "df-correction")
    frame <- compile_frame(whole_brain("local"), domain)
    observations <- nrow(design)
    positions <- seq_len(features)

    draws <- vapply(seq_len(replications), function(replication) {
      sources <- stats::setNames(
        lapply(seq_len(partitions), function(partition) {
          matrix(rnorm(observations * features), observations, features)
        }), paste0("run", seq_len(partitions))
      )
      fit <- lm_relation_fit(
        sources, design, effect_map, sampling_unit = "trial", domain = domain
      )
      plan <- plan_geometry(
        fit$relation, frame,
        cross_partitions(fit$relation, independence = "independent",
          generalizes_over = "run")
      )
      covariance <- rdm_sampling_covariance(
        plan, fit, target = "null", at = 1L
      )
      corrected <- sqrt(sampling_covariance(covariance))

      # The estimator the law describes, computed by the oracle rather than by
      # `rdm()`; the two are pinned equal in
      # test-evidence-sampling-nonspherical.R. The frame weight 1 / P is the
      # normalization.
      blocks <- lapply(fit$relation$partitions, function(partition) {
        relation_block(fit, partition, positions)
      })
      distances <- rowMeans(sampling_oracle_partition_edge_distances(
        blocks, contrasts, normalization = features
      ))

      # Independent reconstruction of the plug-in whitened residual
      # covariance, and of the inflation factor the correction removes. Under
      # `target = "null"` the signal term is exactly zero, so the reported
      # standard error is proportional to the square root of the noise trace
      # and the uncorrected standard error is this multiple of it.
      pooled <- Reduce(`+`, lapply(fit$relation$partitions,
        function(partition) crossprod(residual_block(fit, partition, positions))
      ))
      degrees <- sum(vapply(fit$relation$partitions,
        function(partition) residual_df(fit, partition), integer(1)))
      plug_in <- pooled / degrees / features
      raw_trace <- sum(plug_in * plug_in)
      corrected_trace <- (degrees^2 / ((degrees - 1) * (degrees + 2))) *
        (raw_trace - sum(diag(plug_in))^2 / degrees)

      c(distances, corrected, sqrt(raw_trace / corrected_trace),
        degrees, covariance$source$residual_effective_dimension)
    }, numeric(2L * nrow(contrasts) + 3L))

    dimension <- nrow(contrasts)
    list(
      labels = rownames(contrasts),
      estimates = t(draws[seq_len(dimension), , drop = FALSE]),
      corrected_se = t(draws[dimension + seq_len(dimension), , drop = FALSE]),
      inflation = draws[2L * dimension + 1L, ],
      residual_df = draws[2L * dimension + 2L, ],
      effective_dimension = draws[2L * dimension + 3L, ],
      features = features, partitions = partitions,
      replications = replications
    )
  })
}

test_that("the plug-in residual df correction calibrates the reported SE", {
  experiment <- sampling_df_correction_experiment()
  replications <- experiment$replications
  monte_carlo_se <- 1 / sqrt(2 * (replications - 1))
  band <- 4 * monte_carlo_se

  empirical <- apply(experiment$estimates, 2L, stats::sd)
  corrected <- colMeans(experiment$corrected_se)
  uncorrected <- colMeans(
    experiment$corrected_se * experiment$inflation
  )

  expect_identical(unique(experiment$residual_df), 168)
  # Isotropic residual noise: the effective dimension the package reports is
  # the support size, so the predicted inflation is exact.
  expect_equal(mean(experiment$effective_dimension), experiment$features,
    tolerance = 0.05)
  expect_equal(
    mean(experiment$inflation),
    sqrt(1 + (1 + experiment$features) / 168), tolerance = 0.02
  )

  # The correction lands inside the pre-registered band ...
  expect_lt(max(abs(empirical / corrected - 1)), band)
  # ... and the estimator the package used before 2026-08-16 does not. It is
  # low, meaning it overstated the standard error.
  expect_true(all(empirical / uncorrected < 1 - band))
})

test_that("the correction is the Wishart-unbiased quadratic estimator", {
  # Algebraic identity rather than simulation: for S ~ W_P(nu, Sigma) / nu,
  # E[tr(S^2) - tr(S)^2 / nu] = tr(Sigma^2) (nu - 1)(nu + 2) / nu^2, so the
  # estimator below is unbiased. Checked here against a direct Wishart
  # expectation computed by Monte Carlo in the residual channel alone, which
  # involves no crossform code path.
  features <- 12L
  degrees <- 30L
  sigma <- toeplitz(0.6^(0:(features - 1L)))
  root <- chol(sigma)
  set.seed(4242)
  raw <- numeric(4000L)
  unbiased <- numeric(4000L)
  for (draw in seq_along(raw)) {
    noise <- matrix(rnorm(degrees * features), degrees, features) %*% root
    sample_covariance <- crossprod(noise) / degrees
    raw[[draw]] <- sum(sample_covariance * sample_covariance)
    unbiased[[draw]] <- (degrees^2 / ((degrees - 1) * (degrees + 2))) *
      (raw[[draw]] - sum(diag(sample_covariance))^2 / degrees)
  }
  truth <- sum(sigma * sigma)

  expect_equal(mean(raw) / truth,
    1 + (1 + sum(diag(sigma))^2 / truth) / degrees, tolerance = 0.02)
  expect_equal(mean(unbiased) / truth, 1, tolerance = 0.02)

  # And the package's internal helper computes exactly that estimator.
  sample_covariance <- crossprod(
    matrix(rnorm(degrees * features), degrees, features) %*% root
  ) / degrees
  quadratic <- crossform:::.sampling_unbiased_noise_trace(
    sample_covariance, degrees
  )
  expect_true(quadratic$corrected)
  expect_equal(quadratic$value,
    (degrees^2 / ((degrees - 1) * (degrees + 2))) *
      (sum(sample_covariance * sample_covariance) -
        sum(diag(sample_covariance))^2 / degrees))
  expect_equal(quadratic$effective_dimension,
    sum(diag(sample_covariance))^2 / quadratic$value)

  # `NULL` degrees of freedom means "this is Sigma, not an estimate of it",
  # which is what an oracle or a known-covariance caller supplies.
  known <- crossform:::.sampling_unbiased_noise_trace(sigma, NULL)
  expect_false(known$corrected)
  expect_equal(known$value, truth)
})

test_that("too few residual df for the support is a capability refusal", {
  # Two runs of a short design over a large support: the residual covariance
  # spreads over more effective directions than the residual df can estimate.
  features <- 60L
  conditions <- 3L
  labels <- letters[seq_len(conditions)]
  set.seed(20260817L)
  order <- rep(labels, each = 4L)
  design <- stats::model.matrix(~ 0 + factor(order, levels = labels))
  colnames(design) <- labels
  effect_map <- diag(conditions)
  dimnames(effect_map) <- list(labels, labels)
  domain <- abstract_domain(features, id = "df-starved")
  sources <- stats::setNames(lapply(seq_len(3L), function(partition) {
    matrix(rnorm(nrow(design) * features), nrow(design), features)
  }), paste0("run", seq_len(3L)))
  fit <- lm_relation_fit(
    sources, design, effect_map, sampling_unit = "trial", domain = domain
  )
  plan <- plan_geometry(
    fit$relation, compile_frame(whole_brain("local"), domain),
    cross_partitions(fit$relation, independence = "independent")
  )

  refusal <- catch_refusal(
    rdm_sampling_covariance(plan, fit, target = "null", at = 1L)
  )

  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "sufficient_residual_df")
  expect_identical(refusal$namespace, "evidence_sampling")
  expect_identical(refusal$reasons, "residual_df_below_effective_dimension")
  expect_match(conditionMessage(refusal), "residual degrees of freedom")
  expect_match(refusal$remedies, "smaller support", all = FALSE)
})
