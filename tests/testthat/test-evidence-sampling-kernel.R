sampling_kernel_fixture <- function(dimension = 7L) {
  set.seed(81320)
  design <- cbind(1, condition = rep(c(-0.5, 0.5), 8))
  effects <- rbind(level = c(1, 0), condition = c(0, 1))
  domain <- abstract_domain(7L, id = "sampling-kernel-domain")
  sources <- stats::setNames(lapply(seq_len(4L), function(index) {
    matrix(rnorm(16 * 7), 16, 7)
  }), paste0("run", seq_len(4L)))
  fit <- lm_relation_fit(
    sources, design, effects, sampling_unit = "trial", domain = domain
  )
  evidence <- plan_geometry(
    fit$relation, compile_frame(whole_brain(), domain),
    cross_partitions(fit$relation, independence = "independent")
  )
  base <- list(fit = fit, evidence = evidence)
  set.seed(81321)
  delta_raw <- matrix(rnorm(dimension^2), dimension)
  xi_raw <- matrix(rnorm(dimension^2), dimension)
  delta <- tcrossprod(delta_raw) / dimension
  xi <- tcrossprod(xi_raw) / dimension + diag(0.2, dimension)
  plan <- crossform:::.compile_evidence_sampling_plan(
    base$evidence, base$fit
  )
  covariance <- crossform:::.sampling_covariance_form(
    plan, sqrt(0.37) * delta_raw / sqrt(dimension),
    t(chol(xi)), noise_trace = 0.37,
    partitions = 4L, labels = paste0("d", seq_len(dimension)),
    source = list(oracle = "test-fixture")
  )
  dense <- sampling_oracle_eq13(
    delta, xi,
    sigma_r = diag(sqrt(0.37 * 13), 13L),
    partitions = 4L
  )$covariance
  list(base = base, plan = plan, covariance = covariance, dense = dense)
}

test_that("structured covariance equals the independent Eq. 13 oracle", {
  fixture <- sampling_kernel_fixture()
  observed <- crossform:::.sampling_covariance_materialize(
    fixture$covariance
  )

  expect_equal(unname(observed), unname(fixture$dense), tolerance = 2e-14)
  expect_identical(dimnames(observed), list(
    fixture$covariance$labels, fixture$covariance$labels
  ))
})

test_that("component constructor preserves Eq. 13 without dense factors", {
  fixture <- sampling_kernel_fixture()
  set.seed(81323)
  conditions <- 5L
  features <- 11L
  # The component constructor now binds the contrasts to the plan's effect
  # space, so the oracle needs a plan whose experimental dimension is the
  # oracle's own five conditions.
  labels <- paste0("condition", seq_len(conditions))
  design <- kronecker(rep(1, 4L), diag(conditions))
  colnames(design) <- labels
  effect_map <- diag(conditions)
  dimnames(effect_map) <- list(labels, labels)
  domain <- abstract_domain(3L, id = "sampling-kernel-five-domain")
  sources <- stats::setNames(lapply(seq_len(4L), function(index) {
    matrix(rnorm(4L * conditions * 3L), 4L * conditions, 3L)
  }), paste0("run", seq_len(4L)))
  five_fit <- lm_relation_fit(
    sources, design, effect_map, sampling_unit = "trial", domain = domain
  )
  five_plan <- crossform:::.compile_evidence_sampling_plan(
    plan_geometry(
      five_fit$relation, compile_frame(whole_brain(), domain),
      cross_partitions(five_fit$relation, independence = "independent")
    ),
    five_fit
  )
  contrasts <- sampling_oracle_condition_contrasts(conditions)
  patterns <- matrix(rnorm(conditions * features), conditions, features)
  effect_raw <- matrix(rnorm(conditions^2), conditions, conditions)
  effect_covariance <- tcrossprod(effect_raw) / conditions +
    diag(0.2, conditions)
  residual_covariance <- toeplitz(0.25^(0:(features - 1L)))
  components <- sampling_oracle_components(
    patterns, effect_covariance, residual_covariance
  )
  expected <- sampling_oracle_eq13(
    components$delta, components$xi, residual_covariance,
    five_plan$partition$count
  )$covariance
  observed <- crossform:::.sampling_covariance_from_components(
    five_plan, contrasts, patterns, effect_covariance,
    residual_covariance, labels = rownames(contrasts)
  )

  expect_equal(
    unname(crossform:::.sampling_covariance_materialize(observed)),
    unname(expected), tolerance = 3e-13
  )
  expect_identical(dim(observed$signal_factor),
    c(nrow(contrasts), features))
  expect_identical(nrow(observed$xi_factor), nrow(contrasts))
  expect_false(any(c("delta", "xi", "covariance") %in% names(observed)))
})

test_that("every exact query route agrees with dense covariance", {
  fixture <- sampling_kernel_fixture()
  dense <- crossform:::.sampling_covariance_materialize(
    fixture$covariance
  )
  set.seed(81322)
  vector <- rnorm(nrow(dense))
  matrix_value <- matrix(rnorm(nrow(dense) * 3L), nrow(dense), 3L)
  transport <- matrix(rnorm(4L * nrow(dense)), 4L, nrow(dense))
  selected <- matrix(c(1L, 1L, 2L, 5L, 7L, 3L), byrow = TRUE, ncol = 2L)

  expect_equal(
    crossform:::.sampling_covariance_diagonal(fixture$covariance),
    stats::setNames(diag(dense), fixture$covariance$labels),
    tolerance = 2e-14
  )
  expect_equal(
    crossform:::.sampling_covariance_entries(
      fixture$covariance, selected[, 1L], selected[, 2L]
    ),
    dense[selected], tolerance = 2e-14
  )
  expect_equal(
    crossform:::.sampling_covariance_apply(fixture$covariance, vector),
    drop(dense %*% vector), tolerance = 3e-14
  )
  expect_equal(
    crossform:::.sampling_covariance_apply(
      fixture$covariance, matrix_value
    ),
    dense %*% matrix_value, tolerance = 3e-14
  )
  expect_equal(
    crossform:::.sampling_covariance_quadratic(
      fixture$covariance, vector
    ),
    drop(crossprod(vector, dense %*% vector)), tolerance = 3e-14
  )
  expect_equal(
    crossform:::.sampling_covariance_transport(
      fixture$covariance, transport
    ),
    transport %*% dense %*% t(transport), tolerance = 5e-14
  )
  named_vector <- stats::setNames(vector, rev(fixture$covariance$labels))
  aligned_vector <- named_vector[fixture$covariance$labels]
  named_transport <- transport[, rev(seq_len(ncol(transport))), drop = FALSE]
  colnames(named_transport) <- rev(fixture$covariance$labels)
  expect_equal(
    crossform:::.sampling_covariance_apply(
      fixture$covariance, named_vector
    ),
    drop(dense %*% aligned_vector), tolerance = 3e-14
  )
  expect_equal(
    crossform:::.sampling_covariance_transport(
      fixture$covariance, named_transport
    ),
    transport %*% dense %*% t(transport), tolerance = 5e-14
  )
  expect_error(
    crossform:::.sampling_covariance_apply(
      fixture$covariance,
      stats::setNames(vector, paste0("wrong", seq_along(vector)))
    ),
    "identify every sampling-covariance coordinate"
  )
})

test_that("compiled operation plans execute the requested exact route", {
  fixture <- sampling_kernel_fixture()
  vector <- seq_len(fixture$covariance$dimension) / 10
  operation <- crossform:::.sampling_operation("quadratic_form", vector)
  plan <- crossform:::.rebind_evidence_sampling_plan(
    fixture$plan, operation,
    crossform:::.sampling_materialization("query_only")
  )
  observed <- crossform:::.execute_evidence_sampling_plan(
    plan, fixture$covariance
  )

  expect_equal(
    observed,
    crossform:::.sampling_covariance_quadratic(
      fixture$covariance, vector
    ),
    tolerance = 0
  )
  expect_false(identical(plan$scientific_plan_id,
    fixture$plan$scientific_plan_id))
})

test_that("dense materialization is explicit and size-preflighted", {
  fixture <- sampling_kernel_fixture(dimension = 20L)
  materialize <- crossform:::.rebind_evidence_sampling_plan(
    fixture$plan,
    crossform:::.sampling_operation("materialize"),
    crossform:::.sampling_materialization("dense_covariance")
  )

  expect_error(
    crossform:::.execute_evidence_sampling_plan(
      materialize, fixture$covariance, max_bytes = 8
    ),
    "exceeding.*materialization budget"
  )
  expect_identical(
    dim(crossform:::.execute_evidence_sampling_plan(
      materialize, fixture$covariance,
      max_bytes = 8 * 20^2 + 64 * 20
    )),
    c(20L, 20L)
  )
})

test_that("sampling covariance detects plan and source mutation", {
  fixture <- sampling_kernel_fixture()
  forged <- fixture$covariance
  forged$xi_factor[1L, 1L] <- forged$xi_factor[1L, 1L] + 1

  expect_error(
    crossform:::.validate_sampling_covariance(forged),
    "identity"
  )
  different_target <- crossform:::.compile_evidence_sampling_plan(
    fixture$base$evidence, fixture$base$fit,
    target = crossform:::.sampling_target("null")
  )
  expect_error(
    crossform:::.execute_evidence_sampling_plan(
      different_target, fixture$covariance
    ),
    "different scientific identities"
  )
})
