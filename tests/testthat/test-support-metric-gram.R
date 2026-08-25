# The whitened Gram route is an execution route for the support-streamed
# pair-difference lowering, not a second estimand. What has to be shown is
# that it reproduces the declared numbers, that it declines wherever it is not
# the cheaper association, and that declining leaves the reference route --
# its values, its diagnostics and its observations -- exactly as it was.

gram_route_fixture <- function(effects = 6L, side = 3L, partitions = 3L,
                               radius = 1.5, seed = 20240824) {
  set.seed(seed)
  n_features <- side^3
  coordinates <- as.matrix(expand.grid(
    x = seq_len(side), y = seq_len(side), z = seq_len(side)
  ))
  domain <- abstract_domain(n_features, coordinates = coordinates,
    id = "gram-route", coordinate_units = "mm")
  scale <- seq(0.8, 1.3, length.out = n_features)
  covariance <- toeplitz(0.5^(seq_len(n_features) - 1L))
  covariance <- diag(scale) %*% covariance %*% diag(scale)
  names <- sprintf("e%02d", seq_len(effects))
  truth <- matrix(stats::rnorm(effects * n_features, sd = 0.4),
    effects, n_features, dimnames = list(names, NULL))
  blocks <- stats::setNames(lapply(seq_len(partitions), function(partition) {
    truth + matrix(stats::rnorm(effects * n_features), effects, n_features)
  }), paste0("run", seq_len(partitions)))
  relation <- relation(blocks, domain = domain)
  list(
    domain = domain,
    relation = relation,
    blocks = blocks,
    effects = names,
    covariance = covariance,
    metric = noise_precision(solve(covariance), domain,
      covariance = covariance),
    frame = compile_frame(
      searchlights(radius = radius, normalization = "local"), domain
    ),
    over = cross_partitions(relation, independence = "independent")
  )
}

# An independent reference: the pair-difference RDM written straight from the
# definition, out of the dense frame weights, the dense metric and the
# relation blocks. It shares no code with either execution route.
gram_route_oracle <- function(fixture, pairs = NULL) {
  weights <- as.matrix(fixture$frame$weights)
  metric <- fixture$metric$value
  effects <- fixture$effects
  if (is.null(pairs)) pairs <- t(utils::combn(seq_along(effects), 2L))
  values <- vapply(seq_len(nrow(weights)), function(node) {
    support <- which(weights[node, ] > 0)
    root <- sqrt(weights[node, support])
    local_metric <- metric[support, support, drop = FALSE] * tcrossprod(root)
    accumulated <- matrix(0, length(effects), length(effects))
    for (edge in seq_len(nrow(fixture$over))) {
      left <- fixture$blocks[[fixture$over$left[[edge]]]][, support, drop = FALSE]
      right <- fixture$blocks[[fixture$over$right[[edge]]]][, support, drop = FALSE]
      accumulated <- accumulated +
        fixture$over$weight[[edge]] * (left %*% local_metric %*% t(right))
    }
    vapply(seq_len(nrow(pairs)), function(pair) {
      i <- pairs[pair, 1L]
      j <- pairs[pair, 2L]
      accumulated[i, i] + accumulated[j, j] -
        accumulated[i, j] - accumulated[j, i]
    }, numeric(1))
  }, numeric(nrow(pairs)))
  matrix(values, nrow(weights), nrow(pairs), byrow = TRUE)
}

# Counts how often the route accepted, and returns it to its normal self when
# the calling test finishes.
local_gram_route_spy <- function(accepted, envir = parent.frame()) {
  route <- crossform:::.support_metric_gram_route
  testthat::local_mocked_bindings(
    .support_metric_gram_route = function(...) {
      value <- route(...)
      if (!is.null(value)) accepted$count <- accepted$count + 1L
      value
    },
    .package = "crossform",
    .env = envir
  )
}

local_gram_route_declined <- function(envir = parent.frame()) {
  testthat::local_mocked_bindings(
    .support_metric_gram_route = function(...) NULL,
    .package = "crossform",
    .env = envir
  )
}

test_that("the Gram route reproduces an independent pair-difference oracle", {
  fixture <- gram_route_fixture()
  plan <- plan_geometry(fixture$relation, at = fixture$frame,
    over = fixture$over, metric = fixture$metric)
  accepted <- new.env(parent = emptyenv())
  accepted$count <- 0L
  local_gram_route_spy(accepted)
  observed <- rdm(plan)

  expect_identical(accepted$count, 1L)
  expect_equal(unname(as.matrix(observed$values)),
    gram_route_oracle(fixture), tolerance = 1e-12)
})

test_that("the Gram route and the reference route agree to double precision", {
  fixture <- gram_route_fixture()
  plan <- plan_geometry(fixture$relation, at = fixture$frame,
    over = fixture$over, metric = fixture$metric)
  native <- rdm(plan)
  reference <- local({
    local_gram_route_declined()
    rdm(plan)
  })

  expect_equal(as.matrix(native$values), as.matrix(reference$values),
    tolerance = 1e-13)
  expect_identical(native$pairs, reference$pairs)
  expect_identical(native$index, reference$index)
  expect_identical(
    native$receipt$observed$task_counts,
    reference$receipt$observed$task_counts
  )
  expect_identical(
    native$receipt$observed$features_completed,
    reference$receipt$observed$features_completed
  )
})

test_that("the Gram route reports the same support tasks it performed", {
  fixture <- gram_route_fixture()
  plan <- plan_geometry(fixture$relation, at = fixture$frame,
    over = fixture$over, metric = fixture$metric)
  view <- evaluate_geometry(plan,
    query = crossform:::.pair_difference_query(fixture$effects))
  diagnostics <- view$metadata$diagnostics$total

  expect_identical(diagnostics$support_tasks,
    nrow(fixture$frame$weights))
  expect_identical(diagnostics$relation_reads,
    length(fixture$relation$partitions))
  expect_false(diagnostics$pair_atoms_materialized)
  expect_false(diagnostics$pair_frame_materialized)
  expect_identical(diagnostics$measurement_kind,
    "static-owned-buffer-accounting")
  expect_gt(diagnostics$max_support_size, 1L)
  expect_identical(diagnostics$metric_factorizations, 0L)
})

test_that("a query naming one pair is answered by the same route", {
  # The Gram route forms every pair whether or not the query asks for them, so
  # a selective query looks like the case where the difference route should
  # win. Measured, it is not: on this frame the Gram route is faster for a
  # single pair at every q tested up to 1600. What the route must do is report
  # exactly the pair asked for.
  fixture <- gram_route_fixture()
  plan <- plan_geometry(fixture$relation, at = fixture$frame,
    over = fixture$over, metric = fixture$metric)
  selected <- cbind(fixture$effects[1L], fixture$effects[2L])
  accepted <- new.env(parent = emptyenv())
  accepted$count <- 0L
  local_gram_route_spy(accepted)
  observed <- rdm(plan, pairs = selected)
  reference <- local({
    local_gram_route_declined()
    rdm(plan, pairs = selected)
  })

  expect_identical(accepted$count, 1L)
  expect_identical(ncol(as.matrix(observed$values)), 1L)
  expect_equal(unname(as.matrix(observed$values)),
    gram_route_oracle(fixture, pairs = cbind(1L, 2L)), tolerance = 1e-12)
  expect_equal(as.matrix(observed$values), as.matrix(reference$values),
    tolerance = 1e-13)
})

test_that("the route declines when its accumulator would not fit its budget", {
  # The q-by-q accumulator is the one resource in which the Gram route is
  # worse than the difference route, and the budget is what bounds it. Below
  # the budget nothing else may decline a query this size.
  fixture <- gram_route_fixture()
  plan <- plan_geometry(fixture$relation, at = fixture$frame,
    over = fixture$over, metric = fixture$metric)
  accepted <- new.env(parent = emptyenv())
  accepted$count <- 0L
  local_gram_route_spy(accepted)
  testthat::local_mocked_bindings(
    .support_gram_accumulator_budget_bytes = function() 8 * 6^2 - 1,
    .package = "crossform"
  )
  declined <- rdm(plan)

  expect_identical(accepted$count, 0L)
  expect_equal(unname(as.matrix(declined$values)), gram_route_oracle(fixture),
    tolerance = 1e-12)
})

test_that("the coherent component keeps the reference route", {
  fixture <- gram_route_fixture()
  plan <- plan_geometry(fixture$relation, at = fixture$frame,
    over = fixture$over, metric = fixture$metric)
  accepted <- new.env(parent = emptyenv())
  accepted$count <- 0L
  local_gram_route_spy(accepted)
  observed <- rdm(plan, component = "coherent")

  expect_identical(accepted$count, 0L)
  expect_identical(ncol(as.matrix(observed$values)), 15L)
})

# The kernel arguments, for the two tests that call it directly rather than
# through a plan. Partition endpoints are zero-based and ascending; a pair of
# effect index vectors selects the reported columns.
gram_route_arguments <- function(fixture, metric = fixture$metric$value,
                                 pairs = NULL, path = 1L) {
  rows <- methods::as(fixture$frame$weights, "RsparseMatrix")
  n_features <- ncol(fixture$frame$weights)
  if (is.null(pairs)) pairs <- t(utils::combn(length(fixture$effects), 2L))
  blocks <- unname(lapply(fixture$relation$partitions, function(partition) {
    fixture$blocks[[partition]]
  }))
  edges <- t(utils::combn(length(blocks), 2L)) - 1L
  list(
    rows@p, rows@j, rows@x, metric,
    as.integer(seq_len(n_features) - 1L), blocks,
    as.integer(edges[, 1L]), as.integer(edges[, 2L]),
    rep(1 / (2 * nrow(edges)), nrow(edges)),   # the half-edge weight
    as.integer(pairs[, 1L] - 1L), as.integer(pairs[, 2L] - 1L), path
  )
}

# The same numbers straight from the definition, for a metric supplied
# directly rather than through `neural_metric()`.
gram_route_metric_oracle <- function(fixture, metric, weight) {
  weights <- as.matrix(fixture$frame$weights)
  effects <- fixture$effects
  blocks <- fixture$blocks[fixture$relation$partitions]
  pairs <- t(utils::combn(seq_along(effects), 2L))
  values <- vapply(seq_len(nrow(weights)), function(node) {
    support <- which(weights[node, ] > 0)
    root <- sqrt(weights[node, support])
    local_metric <- metric[support, support, drop = FALSE] * tcrossprod(root)
    accumulated <- matrix(0, length(effects), length(effects))
    for (a in seq_along(blocks)) for (b in seq_along(blocks)) {
      if (a == b) next
      accumulated <- accumulated + weight *
        (blocks[[a]][, support, drop = FALSE] %*% local_metric %*%
           t(blocks[[b]][, support, drop = FALSE]))
    }
    vapply(seq_len(nrow(pairs)), function(pair) {
      i <- pairs[pair, 1L]
      j <- pairs[pair, 2L]
      accumulated[i, i] + accumulated[j, j] -
        accumulated[i, j] - accumulated[j, i]
    }, numeric(1))
  }, numeric(nrow(pairs)))
  matrix(values, nrow(weights), nrow(pairs), byrow = TRUE)
}

test_that("an indefinite local metric is contracted, not refused", {
  # The reduction factors no metric, so the route carries the same admission
  # as the reference loop: a metric that no Cholesky would accept is still a
  # bilinear form, and the estimand it defines is still the declared one.
  fixture <- gram_route_fixture()
  indefinite <- fixture$metric$value
  indefinite[1L, 1L] <- -abs(indefinite[1L, 1L])
  indefinite[2L, 3L] <- indefinite[3L, 2L] <- 5
  expect_lt(min(eigen(indefinite, symmetric = TRUE, only.values = TRUE)$values), 0)

  observed <- do.call(crossform:::.support_metric_gram_pairs_cpp,
    gram_route_arguments(fixture, metric = indefinite))
  edges <- nrow(t(utils::combn(length(fixture$relation$partitions), 2L)))

  expect_equal(observed$value,
    gram_route_metric_oracle(fixture, indefinite, 1 / (2 * edges)),
    tolerance = 1e-12)
  expect_true(all(is.finite(observed$value)))
})

test_that("the two Gram layouts compute the same numbers", {
  fixture <- gram_route_fixture()
  sweep <- do.call(crossform:::.support_metric_gram_pairs_cpp,
    gram_route_arguments(fixture, path = 1L))
  paired <- do.call(crossform:::.support_metric_gram_pairs_cpp,
    gram_route_arguments(fixture, path = 2L))

  expect_equal(sweep$value, paired$value, tolerance = 1e-13)
  expect_identical(sweep$mass, paired$mass)
  expect_identical(sweep$max_support, paired$max_support)
  expect_equal(unname(sweep$value), gram_route_oracle(fixture),
    tolerance = 1e-12)
})

test_that("a subset of pairs reports exactly the columns it names", {
  fixture <- gram_route_fixture()
  chosen <- rbind(c(2L, 5L), c(1L, 3L))
  full <- do.call(crossform:::.support_metric_gram_pairs_cpp,
    gram_route_arguments(fixture))
  subset <- do.call(crossform:::.support_metric_gram_pairs_cpp,
    gram_route_arguments(fixture, pairs = chosen))
  all_pairs <- t(utils::combn(length(fixture$effects), 2L))
  position <- match(paste(chosen[, 1L], chosen[, 2L]),
    paste(all_pairs[, 1L], all_pairs[, 2L]))

  expect_identical(ncol(subset$value), 2L)
  expect_equal(subset$value, full$value[, position, drop = FALSE],
    tolerance = 1e-13)
})
