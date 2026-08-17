oracle_metric_components <- function(K, a) {
  denominator <- drop(crossprod(a, solve(K, a)))
  coherent <- tcrossprod(a) / denominator
  list(
    coherent = coherent,
    configuration = K - coherent,
    denominator = denominator
  )
}

oracle_frame_metric <- function(weight, K) {
  root <- sqrt(weight)
  K * tcrossprod(root)
}

test_that("same-space metric capabilities gate exact lowering classes", {
  domain <- abstract_domain(3, id = "metric-capabilities")
  identity <- neural_metric(diag(3), domain)
  diagonal <- neural_metric(diag(c(1, 2, 3)), domain)
  almost_diagonal <- diag(3)
  almost_diagonal[1, 2] <- almost_diagonal[2, 1] <- 1e-16
  dense <- neural_metric(almost_diagonal, domain)
  recipe <- crossform:::.metric_recipe(
    "shrinkage_precision", domain,
    native_diagonal = FALSE,
    positive_definite = TRUE,
    inverse_quadratic_recipe = TRUE,
    hyperparameters = list(shrinkage = "oas")
  )

  expect_true(metric_capabilities(identity)$identity)
  expect_true(metric_capabilities(diagonal)$feature_additive)
  expect_false(metric_capabilities(dense)$feature_additive)
  expect_true(metric_capabilities(dense)$support_dense)
  expect_true(metric_capabilities(recipe)$learned_recipe)
  expect_false(metric_capabilities(recipe)$materialized)
  expect_true(metric_capabilities(recipe)$inverse_quadratic_recipe)
  expect_identical(crossform:::.metric_lowering(identity),
    "additive_contraction")
  expect_identical(crossform:::.metric_lowering(dense),
    "support_streamed_pair_contraction")
  expect_identical(crossform:::.metric_lowering(recipe),
    "derive_then_support_streamed_pair_contraction")
})

test_that("metric and bridge roles cannot be confused", {
  domain <- abstract_domain(2, id = "metric-vs-bridge")
  bridge <- measurement_bridge(
    diag(2), diag(2), domain, domain,
    measurement_space(2, "shared:metric-vs-bridge")
  )

  expect_error(metric_capabilities(bridge), "bridge.*not.*metric",
    class = "effect_input_error")
  expect_error(crossform:::metric_components(bridge, c(0.5, 0.5)),
    "bridges do not admit", class = "effect_input_error")
})

test_that("support-local metrics never require a domain-wide dense matrix", {
  domain <- abstract_domain(50000, id = "support-local-metric")
  support <- c(101L, 9001L, 49999L)
  metric <- neural_metric(diag(c(2, 3, 4)), domain, support = support)

  expect_identical(dim(metric$value), c(3L, 3L))
  expect_identical(metric$support, support)
  expect_identical(metric$positions, support)
  expect_lt(as.double(utils::object.size(metric$value)), 1024)
})

test_that("frame weights and metrics compose by square-root congruence", {
  domain <- abstract_domain(4, id = "frame-metric-law")
  weight <- c(0, 0.2, 0.3, 0.5)
  frame <- additive_frame(matrix(weight, 1), domain = domain)
  support <- 2:4
  base_value <- matrix(c(
    2.0, 0.3, -0.1,
    0.3, 1.5, 0.2,
    -0.1, 0.2, 1.2
  ), 3, 3, byrow = TRUE)
  base <- neural_metric(base_value, domain, support = support)
  composed <- crossform:::.compose_frame_metric(frame, base, 1)

  expect_equal(composed$metric$value,
    oracle_frame_metric(weight[support], base_value), tolerance = 1e-14)
  expect_equal(composed$coherent$value,
    weight[support] / sum(weight[support]), tolerance = 0)
  expect_identical(composed$composition, "sqrt_weight_congruence")
  expect_identical(crossform:::.metric_lowering(composed$metric),
    "support_streamed_pair_contraction")
})

test_that("identity base metrics recover the existing additive metric exactly", {
  domain <- abstract_domain(4, id = "identity-frame-metric")
  weight <- c(0.1, 0.2, 0.3, 0.4)
  frame <- additive_frame(matrix(weight, 1), normalization = "local",
    domain = domain)
  composed <- crossform:::.compose_frame_metric(frame, NULL, 1)

  expect_equal(composed$metric$value, diag(weight), tolerance = 0)
  expect_true(metric_capabilities(composed$metric)$feature_additive)
  expect_identical(crossform:::.metric_lowering(composed$metric),
    "additive_contraction")
})

test_that("conservative identity conservation is capability-gated", {
  domain <- abstract_domain(3, id = "metric-conservation")
  weights <- rbind(c(0.5, 1, 0), c(0.5, 0, 1))
  frame <- additive_frame(weights, normalization = "conservative",
    domain = domain)
  identity <- crossform:::.metric_frame_conservation(frame)
  dense_metrics <- list(
    neural_metric(matrix(c(1, 0.2, 0.2, 1), 2), domain, c(1L, 2L)),
    neural_metric(matrix(c(1, -0.1, -0.1, 1), 2), domain, c(1L, 3L))
  )
  dense <- crossform:::.metric_frame_conservation(frame, dense_metrics)

  expect_true(identity$feature_additive)
  expect_true(identity$identity_conservation)
  expect_equal(identity$global_diagonal, rep(1, 3), tolerance = 0)
  expect_silent(crossform:::.require_metric_conservation(identity, "identity"))
  expect_false(dense$feature_additive)
  expect_false(dense$identity_conservation)
  expect_null(dense$global_diagonal)
  expect_error(crossform:::.require_metric_conservation(dense, "identity"),
    "not certified.*non-diagonal", class = "effect_input_error")
})

test_that("rank-one coherent and PSD configuration laws hold", {
  set.seed(8201)
  domain <- abstract_domain(8, id = "component-laws")
  raw <- matrix(rnorm(64), 8, 8)
  K <- crossprod(raw) + diag(0.75, 8)
  a <- seq_len(8) / sum(seq_len(8))
  metric <- neural_metric(K, domain)
  functional <- crossform:::coherent_functional(a, domain, label = "weighted_mean")
  got <- crossform:::metric_components(metric, functional)
  expected <- oracle_metric_components(K, a)

  expect_identical(got$coherent_rank, 1L)
  expect_true(got$configuration_psd)
  expect_equal(got$coherent, expected$coherent, tolerance = 3e-13)
  expect_equal(got$configuration, expected$configuration,
    tolerance = 3e-13)
  expect_equal(got$coherent + got$configuration, K, tolerance = 2e-15)
  expect_gte(min(eigen(got$configuration, symmetric = TRUE,
    only.values = TRUE)$values), -1e-10)
})

test_that("the metric split reproduces signed cross-partition components", {
  weight <- c(1, 2, 4, 3)
  mass <- sum(weight)
  domain <- abstract_domain(4, id = "signed-component-law")
  metric <- neural_metric(diag(weight), domain)
  functional <- crossform:::coherent_functional(weight / mass, domain)
  components <- crossform:::metric_components(metric, functional)
  left <- matrix(c(
    1, -2, 0.5, 3,
    -1, 0.5, 2, -0.25
  ), 2, 4, byrow = TRUE)
  right <- -left
  total <- left %*% diag(weight) %*% t(right)
  coherent <- tcrossprod(
    drop(left %*% weight), drop(right %*% weight)
  ) / mass
  observed_coherent <- left %*% components$coherent %*% t(right)
  observed_configuration <-
    left %*% components$configuration %*% t(right)

  expect_equal(observed_coherent, coherent, tolerance = 2e-14)
  expect_equal(observed_coherent + observed_configuration, total,
    tolerance = 2e-14)
  expect_true(any(total < 0))
})

test_that("precision metrics keep the raw mean and reuse retained covariance", {
  set.seed(8202)
  dimension <- 6L
  domain <- abstract_domain(dimension, id = "precision-components")
  raw <- matrix(rnorm(dimension^2), dimension)
  Sigma <- crossprod(raw) + diag(0.5, dimension)
  K <- solve(Sigma)
  a <- rep(1 / dimension, dimension)
  relation <- matrix(rnorm(3 * dimension), 3, dimension)
  retained <- neural_metric(K, domain, inverse = Sigma,
    provenance = list(kind = "precision_family"))
  solved <- neural_metric(K, domain,
    provenance = list(kind = "precision_family"))
  retained_components <- crossform:::metric_components(retained, a)
  solved_components <- crossform:::metric_components(solved, a)
  amplitude <- drop(relation %*% a)
  expected <- tcrossprod(amplitude) /
    drop(crossprod(a, Sigma %*% a))

  expect_equal(
    relation %*% retained_components$coherent %*% t(relation),
    expected, tolerance = 4e-13
  )
  expect_equal(retained_components$coherent, solved_components$coherent,
    tolerance = 4e-13)
  expect_identical(retained_components$inverse_quadratic_mode,
    "retained_inverse_metric")
  expect_identical(retained_components$factorization_count, 0L)
  expect_identical(solved_components$inverse_quadratic_mode,
    "cholesky_solve")
  expect_identical(solved_components$factorization_count, 1L)
  expect_identical(retained$signature, solved$signature)
})

test_that("metric-aware components are covariant under neural reparameterization", {
  set.seed(8203)
  dimension <- 5L
  domain <- abstract_domain(dimension, id = "metric-gauge")
  raw <- matrix(rnorm(dimension^2), dimension)
  K <- crossprod(raw) + diag(1, dimension)
  transform <- matrix(rnorm(dimension^2), dimension) + diag(2, dimension)
  while (abs(det(transform)) < 0.1) {
    transform <- transform + diag(0.25, dimension)
  }
  inverse_transform <- solve(transform)
  transformed_K <- t(inverse_transform) %*% K %*% inverse_transform
  a <- seq_len(dimension) / sum(seq_len(dimension))
  transformed_a <- solve(t(transform), a)
  relation <- matrix(rnorm(4 * dimension), 4, dimension)
  transformed_relation <- relation %*% t(transform)
  original <- crossform:::metric_components(neural_metric(K, domain), a)
  changed <- crossform:::metric_components(
    neural_metric(transformed_K, domain), transformed_a
  )

  expect_equal(
    transformed_relation %*% transformed_K %*% t(transformed_relation),
    relation %*% K %*% t(relation), tolerance = 1e-11
  )
  expect_equal(
    transformed_relation %*% changed$coherent %*% t(transformed_relation),
    relation %*% original$coherent %*% t(relation), tolerance = 1e-11
  )
  expect_equal(
    transformed_relation %*% changed$configuration %*%
      t(transformed_relation),
    relation %*% original$configuration %*% t(relation),
    tolerance = 1e-11
  )
})

test_that("singular and unfrozen metrics refuse coherent decomposition", {
  domain <- abstract_domain(3, id = "metric-refusal")
  singular <- neural_metric(diag(c(1, 1, 0)), domain)
  recipe <- crossform:::.metric_recipe(
    "unfrozen_precision", domain,
    inverse_quadratic_recipe = TRUE
  )

  expect_false(metric_capabilities(singular)$positive_definite)
  expect_error(crossform:::metric_components(singular, rep(1 / 3, 3)),
    "requires an SPD metric", class = "effect_input_error")
  expect_error(crossform:::metric_components(recipe, rep(1 / 3, 3)),
    "derived and frozen", class = "effect_input_error")
})

test_that("component identity binds both metric and coherent functional", {
  domain <- abstract_domain(3, id = "component-identity")
  first_metric <- neural_metric(diag(c(1, 2, 3)), domain)
  second_metric <- neural_metric(diag(c(1, 2, 4)), domain)
  first_functional <- crossform:::coherent_functional(c(0.2, 0.3, 0.5), domain)
  second_functional <- crossform:::coherent_functional(c(0.3, 0.2, 0.5), domain)
  baseline <- crossform:::metric_components(first_metric, first_functional)

  expect_false(identical(
    baseline$signature,
    crossform:::metric_components(second_metric, first_functional)$signature
  ))
  expect_false(identical(
    baseline$signature,
    crossform:::metric_components(first_metric, second_functional)$signature
  ))
})

test_that("mutated metric identities fail at the deep boundary", {
  domain <- abstract_domain(2, id = "metric-mutation")
  metric <- neural_metric(diag(2), domain)
  metric$value[1, 1] <- 2

  expect_error(crossform:::.validate_neural_metric(metric, deep = TRUE),
    "identity", class = "effect_contract_error")
})

test_that("a folded diagonal metric keeps the declared normalization", {
  domain <- abstract_domain(4, id = "metric-fold-provenance")
  precision <- noise_precision(diag(c(0.5, 2, 1.25, 3)), domain)
  labels <- c("a", "a", "b", "b")

  conservative <- compile_frame(
    regions(labels, normalization = "conservative"), domain
  )
  schedule <- crossform:::.geometry_metric_schedule(conservative, precision)
  folded <- crossform:::.metric_additive_frame(conservative, schedule)

  # The composed weights are no longer column-normalized, so the frame cannot
  # claim "conservative" without lying about the operator it carries. The
  # declared normalization survives as provenance instead.
  expect_identical(folded$normalization, "none")
  expect_true(folded$metric_folded$folded)
  expect_identical(folded$metric_folded$declared_normalization, "conservative")
  expect_identical(folded$metric_folded$metric_kind, "native_diagonal")
  expect_identical(folded$metric_folded$composition, "diagonal_metric_fold")
  expect_identical(folded$metric_folded$metric_signature, precision$signature)
  expect_equal(folded$metric_folded$reference_mass, c(0.5, 2, 1.25, 3),
    tolerance = 0)

  # Conservation is unchanged by the fold: the conserving per-feature mass is
  # the metric diagonal, because the global comparator is read under the same
  # metric.
  report <- frame_conservation(folded)
  expect_true(report$conserved)
  expect_true(report$metric_folded)
  expect_identical(report$declared_normalization, "conservative")
  expect_identical(report$normalization, "none")
  expect_equal(report$max_deviation, 0, tolerance = 1e-12)

  certificate <- crossform:::.metric_frame_conservation(folded)
  expect_true(certificate$identity_conservation)
  expect_true(certificate$metric_folded)
  expect_identical(certificate$declared_normalization, "conservative")
  expect_match(certificate$reason, "folded metric diagonal")

  # A local frame records "local" and is no more conserved after the fold
  # than before it.
  local_frame <- compile_frame(regions(labels), domain)
  local_folded <- crossform:::.metric_additive_frame(
    local_frame, crossform:::.geometry_metric_schedule(local_frame, precision)
  )
  expect_identical(local_folded$metric_folded$declared_normalization, "local")
  expect_true(local_folded$metric_folded$folded)
  expect_false(frame_conservation(local_folded)$conserved)
  expect_identical(
    frame_conservation(local_folded)$declared_normalization, "local"
  )

  # The implicit identity schedule folds nothing, so the frame is returned
  # untouched and carries no fold provenance.
  untouched <- crossform:::.metric_additive_frame(
    conservative, crossform:::.geometry_metric_schedule(conservative)
  )
  expect_null(untouched$metric_folded)
  expect_identical(untouched$normalization, "conservative")

  # This is why the declared normalization is provenance rather than a claim
  # on `$normalization`: the frame validator checks the column-mass law
  # against the weights it is handed, so a folded frame that declared
  # "conservative" would be refused outright.
  mislabeled <- folded
  mislabeled$normalization <- "conservative"
  expect_error(crossform:::.validate_frame_for_compile(mislabeled),
    "columns must sum to one", class = "effect_input_error")

  forged <- folded
  forged$metric_folded$declared_normalization <- "invented"
  expect_error(crossform:::.validate_frame_for_compile(forged),
    "metric-fold provenance", class = "effect_input_error")
})

test_that("a folded diagonal metric conserves local totals globally", {
  domain <- abstract_domain(5, coordinates = cbind(0:4, 0),
    id = "metric-fold-conservation")
  effects <- list(
    run1 = rbind(a = c(1, 0, 2, 1, -1), b = c(0, 1, 1, 0, 2)),
    run2 = rbind(a = c(1.1, 0.1, 1.9, 0.8, -1.2), b = c(0.1, 0.9, 1.2, 0.2, 1.8))
  )
  rel <- relation(effects, domain = domain)
  over <- cross_partitions(rel, independence = "independent")
  precision <- noise_precision(diag(c(0.5, 2, 1.25, 3, 0.75)), domain)
  query <- bilinear_query(tcrossprod(c(1, -1)))

  local_plan <- plan_geometry(
    rel, compile_frame(searchlights(1.01, normalization = "conservative"),
      domain), over, metric = precision
  )
  global_plan <- plan_geometry(
    rel, compile_frame(whole_brain(normalization = "none"), domain), over,
    metric = precision
  )
  # Pin the route: this conservation law is the one the metric-folded frame
  # implements, so the test must fail if the plan stops using that lowering.
  expect_identical(
    crossform:::.compile_geometry_execution_plan(
      local_plan, query = query, component = "total"
    )$lowering,
    "additive_metric_query_fused_contraction"
  )

  local_view <- evaluate_geometry(local_plan, query = query,
    component = "total")
  global_view <- evaluate_geometry(global_plan, query = query,
    component = "total")

  expect_equal(sum(local_view$values), unname(drop(global_view$values)),
    tolerance = 1e-12)
})
