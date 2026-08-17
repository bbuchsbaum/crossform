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
  # The fixture builds its row factors by hand, so the oracle is applied to
  # the same signal Gram (0.37 Delta) and noise trace (0.37) directly.
  dense <- sampling_oracle_eq13_terms(
    signal_gram = 0.37 * delta, xi = xi, noise_trace = 0.37,
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
    components$differences, components$xi, residual_covariance,
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

test_that("the component constructor equals a scalar-loop law exactly", {
  # An algebraic court for the constructor: random small inputs, an
  # entry-by-entry base-R evaluation of
  #   Cov(d_r, d_s) = 4 Xi_rs (mu_r Sigma_R mu_s') / (M nu^2)
  #                 + 2 Xi_rs^2 tr(Sigma_R^2) / (M (M - 1) nu^2),
  # and no shared code path with the implementation.
  set.seed(81325)
  for (case in seq_len(6L)) {
    conditions <- sample(2:5, 1L)
    features <- sample(2:7, 1L)
    partitions <- sample(2:5, 1L)
    normalization <- sample(c(1, features, 2.5), 1L)
    labels <- paste0("condition", seq_len(conditions))
    design <- kronecker(rep(1, 4L), diag(conditions))
    colnames(design) <- labels
    effect_map <- diag(conditions)
    dimnames(effect_map) <- list(labels, labels)
    domain <- abstract_domain(
      features, id = sprintf("sampling-scalar-law-%d", case)
    )
    sources <- stats::setNames(lapply(seq_len(partitions), function(index) {
      matrix(rnorm(4L * conditions * features), 4L * conditions, features)
    }), paste0("run", seq_len(partitions)))
    fit <- lm_relation_fit(
      sources, design, effect_map, sampling_unit = "trial", domain = domain
    )
    plan <- crossform:::.compile_evidence_sampling_plan(
      plan_geometry(
        fit$relation, compile_frame(whole_brain(), domain),
        cross_partitions(fit$relation, independence = "independent")
      ),
      fit
    )
    contrasts <- sampling_oracle_condition_contrasts(conditions)
    patterns <- matrix(
      rnorm(conditions * features), conditions, features
    )
    effect_raw <- matrix(rnorm(conditions^2), conditions, conditions)
    effect_covariance <- tcrossprod(effect_raw) / conditions +
      diag(0.15, conditions)
    residual_raw <- matrix(
      rnorm(features * (features + 2L)), features, features + 2L
    )
    residual_covariance <- tcrossprod(residual_raw) / (features + 2L)
    observed <- crossform:::.sampling_covariance_materialize(
      crossform:::.sampling_covariance_from_components(
        plan, contrasts, patterns, effect_covariance, residual_covariance,
        normalization = normalization, labels = rownames(contrasts)
      )
    )
    expected <- sampling_oracle_scalar_law(
      contrasts %*% patterns,
      contrasts %*% effect_covariance %*% t(contrasts),
      residual_covariance, partitions, normalization
    )

    expect_equal(unname(observed), expected, tolerance = 1e-12,
      info = sprintf(
        "case %d: %d conditions, %d features, %d partitions, nu = %g",
        case, conditions, features, partitions, normalization
      ))
  }
})

test_that("a rank-deficient residual covariance keeps a valid factor", {
  set.seed(81326)
  conditions <- 4L
  features <- 5L
  partitions <- 3L
  labels <- paste0("condition", seq_len(conditions))
  design <- kronecker(rep(1, 4L), diag(conditions))
  colnames(design) <- labels
  effect_map <- diag(conditions)
  dimnames(effect_map) <- list(labels, labels)
  domain <- abstract_domain(features, id = "sampling-rank-deficient")
  sources <- stats::setNames(lapply(seq_len(partitions), function(index) {
    matrix(rnorm(4L * conditions * features), 4L * conditions, features)
  }), paste0("run", seq_len(partitions)))
  fit <- lm_relation_fit(
    sources, design, effect_map, sampling_unit = "trial", domain = domain
  )
  plan <- crossform:::.compile_evidence_sampling_plan(
    plan_geometry(
      fit$relation, compile_frame(whole_brain(), domain),
      cross_partitions(fit$relation, independence = "independent")
    ),
    fit
  )
  contrasts <- sampling_oracle_condition_contrasts(conditions)
  patterns <- matrix(rnorm(conditions * features), conditions, features)
  effect_covariance <- diag(0.4, conditions)
  deficient_raw <- matrix(rnorm(features * 2L), features, 2L)
  deficient <- tcrossprod(deficient_raw)
  degenerate <- matrix(0, features, features)

  rank_two <- crossform:::.sampling_covariance_from_components(
    plan, contrasts, patterns, effect_covariance, deficient,
    normalization = 1, labels = rownames(contrasts)
  )
  empty <- crossform:::.sampling_covariance_from_components(
    plan, contrasts, patterns, effect_covariance, degenerate,
    normalization = 1, labels = rownames(contrasts)
  )

  expect_identical(ncol(rank_two$signal_factor), 2L)
  expect_equal(
    unname(crossform:::.sampling_covariance_materialize(rank_two)),
    sampling_oracle_scalar_law(
      contrasts %*% patterns,
      contrasts %*% effect_covariance %*% t(contrasts),
      deficient, partitions, 1
    ), tolerance = 1e-12
  )
  expect_identical(ncol(empty$signal_factor), 1L)
  expect_true(all(
    crossform:::.sampling_covariance_materialize(empty) == 0
  ))
  expect_error(
    crossform:::.sampling_covariance_from_components(
      plan, contrasts, patterns, effect_covariance,
      diag(c(1, 1, 1, 1, -1)), normalization = 1,
      labels = rownames(contrasts)
    ),
    "Residual covariance must be positive semidefinite"
  , class = "effect_input_error")
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
  , class = "effect_input_error")
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
  , class = "effect_input_error")
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
  , class = "effect_contract_error")
  different_target <- crossform:::.compile_evidence_sampling_plan(
    fixture$base$evidence, fixture$base$fit,
    target = crossform:::.sampling_target("null")
  )
  expect_error(
    crossform:::.execute_evidence_sampling_plan(
      different_target, fixture$covariance
    ),
    "different scientific identities"
  , class = "effect_contract_error")
})

test_that("a sampling covariance in the general evidence basis renders", {
  fixture <- sampling_kernel_fixture()

  # The digest hashes double-precision row factors, so it is masked; every
  # other line is a rendering decision this snapshot is here to hold.
  expect_snapshot(
    print(fixture$covariance),
    transform = function(lines) {
      sub("(signature:\\s+sha256:)[0-9a-f]+", "\\1<digest>", lines)
    }
  )
  expect_snapshot(format(fixture$covariance))
})
