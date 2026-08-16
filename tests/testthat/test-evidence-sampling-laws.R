# The oracle helpers in helper-evidence-sampling-laws.R are base-R only. Each
# block below states a law with them and then holds crossform's own code to
# that same law, so the file is a court and not a private conversation.

sampling_law_relation <- function(conditions, features, partitions,
                                  seed = 51501L) {
  set.seed(seed)
  labels <- paste0("condition", seq_len(conditions))
  design <- kronecker(rep(1, 4L), diag(conditions))
  colnames(design) <- labels
  effect_map <- diag(conditions)
  dimnames(effect_map) <- list(labels, labels)
  domain <- abstract_domain(
    features,
    id = sprintf("sampling-law-%d-%d-%d", conditions, features, partitions)
  )
  sources <- stats::setNames(lapply(seq_len(partitions), function(index) {
    matrix(rnorm(4L * conditions * features), 4L * conditions, features)
  }), paste0("run", seq_len(partitions)))
  fit <- lm_relation_fit(
    sources, design, effect_map, sampling_unit = "trial", domain = domain
  )
  list(fit = fit, domain = domain, labels = labels, features = features)
}

sampling_law_geometry <- function(built, normalization = "local",
                                  metric = NULL) {
  plan_geometry(
    built$fit$relation,
    compile_frame(whole_brain(normalization), built$domain),
    cross_partitions(built$fit$relation, independence = "independent"),
    metric = metric
  )
}

sampling_law_plan <- function(built, ...) {
  crossform:::.compile_evidence_sampling_plan(
    sampling_law_geometry(built, ...), built$fit
  )
}

sampling_law_package_covariance <- function(built, contrasts, patterns,
                                            sigma_k, sigma_r,
                                            normalization = ncol(patterns),
                                            ...) {
  plan <- sampling_law_plan(built, ...)
  crossform:::.sampling_covariance_materialize(
    crossform:::.sampling_covariance_from_components(
      plan, contrasts, patterns, sigma_k, sigma_r,
      normalization = normalization, labels = rownames(contrasts)
    )
  )
}

test_that("Diedrichsen Eq. 10 is exactly the all-unordered-pairs estimator", {
  set.seed(81311)
  observed <- vapply(2:12, function(partitions) {
    values <- matrix(rnorm(partitions * 17), partitions, 17)
    abs(
      sampling_oracle_leave_one_out(values) -
        sampling_oracle_all_pairs(values)
    )
  }, numeric(1))

  expect_lt(max(observed), 5e-16)

  # crossform's own estimator is that same equal-weight functional. A
  # whole-brain frame with local normalization puts weight 1 / P on every
  # feature, so the composed identity metric divides by P exactly as the
  # oracle does.
  built <- sampling_law_relation(conditions = 4L, features = 9L,
    partitions = 5L)
  contrasts <- sampling_oracle_condition_contrasts(4L)
  estimates <- lapply(built$fit$relation$partitions, function(partition) {
    relation_block(built$fit, partition, seq_len(built$features))
  })
  package <- drop(rdm(sampling_law_geometry(built))$values)
  oracle_all_pairs <- vapply(seq_len(nrow(contrasts)), function(distance) {
    sampling_oracle_all_pairs(t(vapply(estimates, function(value) {
      drop(contrasts[distance, ] %*% value)
    }, numeric(built$features))))
  }, numeric(1))
  oracle_leave_one_out <- vapply(seq_len(nrow(contrasts)), function(distance) {
    sampling_oracle_leave_one_out(t(vapply(estimates, function(value) {
      drop(contrasts[distance, ] %*% value)
    }, numeric(built$features))))
  }, numeric(1))

  expect_equal(unname(package), oracle_all_pairs, tolerance = 1e-13)
  expect_equal(unname(package), oracle_leave_one_out, tolerance = 1e-13)
})

test_that("Eq. 13 equals an independent endpoint-overlap enumeration", {
  set.seed(81312)
  conditions <- 5L
  features <- 13L
  patterns <- matrix(rnorm(conditions * features), conditions, features)
  sigma_k_raw <- matrix(rnorm(conditions * conditions), conditions, conditions)
  sigma_k <- tcrossprod(sigma_k_raw) / conditions + diag(0.3, conditions)
  sigma_r <- toeplitz(0.42^(0:(features - 1L)))
  components <- sampling_oracle_components(patterns, sigma_k, sigma_r)

  for (partitions in c(2L, 3L, 4L, 7L, 12L)) {
    equation <- sampling_oracle_eq13(
      components$differences, components$xi, sigma_r, partitions
    )
    enumerated <- sampling_oracle_endpoint_enumeration(
      components$differences, components$xi, sigma_r, partitions
    )
    built <- sampling_law_relation(conditions, features, partitions)
    package <- sampling_law_package_covariance(
      built, components$contrasts, patterns, sigma_k, sigma_r
    )

    expect_equal(
      enumerated$signal_coefficient, 4 / partitions,
      tolerance = 2e-15
    )
    expect_equal(
      enumerated$noise_coefficient,
      2 / (partitions * (partitions - 1L)),
      tolerance = 2e-15
    )
    expect_equal(
      enumerated$covariance, equation$covariance,
      tolerance = 3e-14
    )
    expect_equal(
      unname(package), unname(equation$covariance), tolerance = 3e-13
    )
  }
})

test_that("Eq. 13 is one complete distance-by-distance covariance", {
  set.seed(81313)
  conditions <- 4L
  features <- 9L
  patterns <- matrix(rnorm(conditions * features), conditions, features)
  sigma_k <- matrix(c(
    1.3, 0.2, 0.1, 0.0,
    0.2, 1.1, 0.3, 0.1,
    0.1, 0.3, 1.4, 0.2,
    0.0, 0.1, 0.2, 1.2
  ), 4, 4, byrow = TRUE)
  sigma_r <- toeplitz(0.3^(0:(features - 1L)))
  components <- sampling_oracle_components(patterns, sigma_k, sigma_r)
  result <- sampling_oracle_eq13(
    components$differences, components$xi, sigma_r, partitions = 6L
  )
  built <- sampling_law_relation(conditions, features, partitions = 6L)
  package <- sampling_law_package_covariance(
    built, components$contrasts, patterns, sigma_k, sigma_r
  )

  expect_identical(dim(result$covariance), c(6L, 6L))
  expect_equal(result$covariance, t(result$covariance), tolerance = 2e-15)
  expect_gt(max(abs(result$covariance[row(result$covariance) !=
    col(result$covariance)])), 0)
  expect_equal(result$covariance, result$signal + result$noise,
    tolerance = 0)
  expect_gt(sum(abs(result$signal)), 0)
  expect_gt(sum(abs(result$noise)), 0)
  expect_gte(min(eigen(result$covariance, symmetric = TRUE,
    only.values = TRUE)$values), -1e-12)

  # The package delivers the same complete off-diagonal law, not a diagonal
  # of variances with the cross terms quietly dropped.
  expect_identical(dim(package), c(6L, 6L))
  expect_equal(unname(package), unname(result$covariance), tolerance = 3e-13)
  expect_gt(max(abs(package[row(package) != col(package)])), 0)
})

test_that("signal and noise terms have their distinct partition scaling", {
  set.seed(81314)
  conditions <- 3L
  features <- 11L
  patterns <- matrix(rnorm(conditions * features), conditions, features)
  sigma_k <- matrix(c(1.2, 0.4, 0.1, 0.4, 0.9, 0.2, 0.1, 0.2, 1.1), 3, 3)
  sigma_r <- diag(features)
  components <- sampling_oracle_components(patterns, sigma_k, sigma_r)

  small <- sampling_oracle_eq13(
    components$differences, components$xi, sigma_r, partitions = 3L
  )
  large <- sampling_oracle_eq13(
    components$differences, components$xi, sigma_r, partitions = 12L
  )

  expect_equal(unname(small$signal / large$signal),
    matrix(4, nrow(components$xi), nrow(components$xi)), tolerance = 2e-15)
  expect_equal(
    unname(small$noise / large$noise),
    matrix((12 * 11) / (3 * 2), nrow(components$xi), nrow(components$xi)),
    tolerance = 2e-15
  )

  # The same two scalings must be visible in package output. Evaluating the
  # package at zero patterns isolates its noise term, and the difference from
  # the full covariance isolates its signal term.
  package_terms <- function(partitions) {
    built <- sampling_law_relation(conditions, features, partitions)
    full <- sampling_law_package_covariance(
      built, components$contrasts, patterns, sigma_k, sigma_r
    )
    noise <- sampling_law_package_covariance(
      built, components$contrasts, 0 * patterns, sigma_k, sigma_r
    )
    list(signal = full - noise, noise = noise)
  }
  package_small <- package_terms(3L)
  package_large <- package_terms(12L)

  expect_equal(unname(package_small$signal / package_large$signal),
    matrix(4, nrow(components$xi), nrow(components$xi)), tolerance = 1e-10)
  expect_equal(
    unname(package_small$noise / package_large$noise),
    matrix((12 * 11) / (3 * 2), nrow(components$xi), nrow(components$xi)),
    tolerance = 1e-10
  )
  expect_equal(unname(package_small$signal), unname(small$signal),
    tolerance = 3e-13)
  expect_equal(unname(package_small$noise), unname(small$noise),
    tolerance = 3e-13)
})

test_that("the two naive signal-dominated SE factors are not conflated", {
  for (partitions in c(3L, 4L, 6L, 8L, 12L)) {
    factors <- sampling_oracle_naive_signal_factors(partitions)
    expect_equal(
      unname(factors[["sample_sd"]]), sqrt(partitions + 1),
      tolerance = 2e-15
    )
    expect_equal(
      unname(factors[["marginal_edge"]]), sqrt(partitions - 1),
      tolerance = 2e-15
    )
  }

  # Both factors are statements about the package's own analytic variance in
  # the signal-dominated limit. Scaling the residual covariance by t sends
  # the signal term to O(t) and the noise term to O(t^2), so at t = 1e-8 the
  # package variance is its signal term to twelve digits. One unit of signal
  # is u = Xi_rr mu_r Sigma_R mu_r' / nu^2; the exact variance of the
  # all-pairs mean is 4u / M, one edge product has marginal variance 2u, and
  # edges sharing an endpoint contribute u each.
  set.seed(81315)
  conditions <- 3L
  features <- 7L
  patterns <- matrix(rnorm(conditions * features), conditions, features)
  sigma_k <- diag(c(1.1, 0.7, 1.3))
  sigma_r <- 1e-8 * toeplitz(0.6^(0:(features - 1L)))
  components <- sampling_oracle_components(patterns, sigma_k, sigma_r)
  unit <- diag(components$xi) *
    diag(components$differences %*% sigma_r %*% t(components$differences)) /
    features^2

  for (partitions in c(3L, 4L, 6L, 8L, 12L)) {
    built <- sampling_law_relation(conditions, features, partitions)
    package <- diag(sampling_law_package_covariance(
      built, components$contrasts, patterns, sigma_k, sigma_r
    ))
    factors <- sampling_oracle_naive_signal_factors(partitions)
    edges <- choose(partitions, 2L)
    marginal_edge_variance <- 2 * unit / edges
    shared <- sampling_oracle_endpoint_enumeration(
      components$differences, components$xi, sigma_r, partitions
    )
    mean_off_diagonal <- (shared$shared_endpoints - 2 * edges) /
      (edges * (edges - 1L))
    sample_sd_variance <- (2 - mean_off_diagonal) * unit / edges

    expect_equal(unname(package / (4 * unit / partitions)),
      rep(1, length(unit)), tolerance = 1e-6)
    expect_equal(
      unname(sqrt(package / marginal_edge_variance)),
      rep(unname(factors[["marginal_edge"]]), length(unit)),
      tolerance = 1e-6
    )
    expect_equal(
      unname(sqrt(package / sample_sd_variance)),
      rep(unname(factors[["sample_sd"]]), length(unit)),
      tolerance = 1e-6
    )
  }
})
