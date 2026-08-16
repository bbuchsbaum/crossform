public_measurement_fixture <- function(rank_one = FALSE,
                                       mode = "scalar",
                                       seed = 2026081228) {
  set.seed(seed)
  q <- 8L
  p <- 4L
  domain <- abstract_domain(p, id = paste0("public:", mode, ":neural:v1"))
  effects <- effect_space(paste0("trial", seq_len(q)),
    basis_id = paste0("public:", mode, ":trials:v1"))
  run1 <- matrix(rnorm(q * p), q, p)
  run2 <- matrix(rnorm(q * p), q, p)
  run1[, 2L] <- -0.75 * run1[, 1L] + 0.5 * run1[, 2L]
  run2[, 2L] <- -0.65 * run2[, 1L] + 0.55 * run2[, 2L]
  rel <- relation(list(run1 = run1, run2 = run2),
    effects = effects, domain = domain)
  operators <- if (mode == "scalar") {
    list(
      anterior = matrix(c(1, 0, 0, 0), 1L,
        dimnames = list("mean", NULL)),
      posterior = matrix(c(0, 1, 0, 0), 1L,
        dimnames = list("mean", NULL))
    )
  } else {
    list(
      anterior = diag(p)[1:2, , drop = FALSE],
      posterior = diag(p)[3:4, , drop = FALSE]
    )
  }
  frame <- measurement_frame(
    operators, domain = domain, id = paste0("public:", mode, ":frame:v1")
  )
  pairs <- expand.grid(
    from = c("anterior", "posterior"),
    to = c("anterior", "posterior"),
    stringsAsFactors = FALSE
  )
  between <- edge_frame(pairs$from, pairs$to, frame)
  h <- if (rank_one) {
    contrast <- c(-1, -1, -1, -1, 1, 1, 1, 1)
    tcrossprod(contrast)
  } else {
    (diag(q) - matrix(1 / q, q, q)) / (q - 1)
  }
  by <- variation_query(
    h, effects, "trial", "joint_covariance",
    provenance = list(estimator = if (rank_one) {
      "single-effect-direction"
    } else {
      "centered-within-run"
    })
  )
  over <- pairing(
    rel$partitions, rel$partitions,
    directed = TRUE,
    self_pairs = "allow_biased",
    independence = "not_independent"
  )
  form <- measurement_form(
    left = rel, between = between, by = by, over = over,
    route = "pull_h"
  )
  list(
    form = form,
    frame = frame,
    between = between,
    by = by,
    rel = rel,
    values = list(run1 = run1, run2 = run2),
    h = h,
    operators = operators
  )
}

test_that("the public boundary vocabulary is narrow and axis-explicit", {
  arguments <- names(formals(measurement_form))
  expect_identical(arguments[1:5],
    c("left", "between", "by", "over", "right"))
  expect_true(all(c(
    "measurement_frame", "edge_frame", "variation_query",
    "measurement_form", "effect_coupling", "covariance_coupling",
    "canonical_coupling", "geometry_alignment", "connectivity",
    "gaussian_covariance_model", "measurement_components",
    "reconstruct_evidence"
  ) %in% getNamespaceExports("crossform")))
  expect_false(any(c(
    "functional_connectivity_model", "informational_connectivity_model",
    "searchlight_connectivity_engine", "roi_connectivity_result"
  ) %in% getNamespaceExports("crossform")))
})

test_that("public measurement frames fail before brain-scale dense conversion", {
  features <- 6000L
  domain <- abstract_domain(features, id = "public:small-node-gate:v1")
  weights <- Matrix::sparseMatrix(
    i = 1L, j = 1L, x = 1, dims = c(1L, features)
  )
  frame <- additive_frame(weights, domain = domain)

  expect_error(
    measurement_frame(frame),
    "small-node limit.*support-local geometry plan"
  )
  expect_equal(
    eval(formals(reconstruct_evidence)$workspace_bytes),
    512 * 1024^2,
    tolerance = 0
  )
})

test_that("measurement and edge frames preserve explicit oriented legs", {
  domain <- abstract_domain(3, id = "public:frame-domain:v1")
  operators <- list(
    seed = matrix(c(1, 0, 0), 1L,
      dimnames = list("common", NULL)),
    target = matrix(c(0, 1, 0, 0, 0, 1), 2L, 3L,
      dimnames = list(c("gradient", "residual"), NULL))
  )
  frame <- measurement_frame(operators, domain,
    id = "public:oriented-frame:v1")
  between <- edge_frame("seed", "target", frame)

  expect_s3_class(frame, "effect_measurement_frame")
  expect_identical(frame$legs$target$operator, unname(operators$target))
  expect_identical(frame$legs$target$output_space$coordinates,
    c("gradient", "residual"))
  expect_s3_class(between, "effect_edge_frame")
  expect_identical(between$edges$edges$left, "seed")
  expect_identical(between$edges$edges$right, "target")
  expect_error(edge_frame(frame = frame), "argument.*missing")

  additive <- compile_frame(
    regions(c("A", "A", "B")), domain
  )
  total <- measurement_frame(additive, mode = "total")
  decomposed <- measurement_frame(
    additive, mode = "coherent_configuration"
  )
  expect_true(all(vapply(total$legs, function(leg) {
    is.null(leg$decomposition)
  }, logical(1))))
  expect_true(all(vapply(decomposed$legs, function(leg) {
    identical(leg$decomposition$components,
      c("coherent", "configuration"))
  }, logical(1))))
})

test_that("measurement_form recovers signed scalar connectivity", {
  fixture <- public_measurement_fixture()
  result <- connectivity(fixture$form, "correlation")
  expected_q <- Reduce(`+`, lapply(fixture$values, function(value) {
    crossprod(value, fixture$h %*% value) / length(fixture$values)
  }))
  expected <- expected_q[1L, 2L] /
    sqrt(expected_q[1L, 1L] * expected_q[2L, 2L])
  edge <- which(fixture$form$block_index$left == "anterior" &
    fixture$form$block_index$right == "posterior")
  got <- result$values$correlation[
    result$values$edge_id == fixture$form$block_index$edge_id[[edge]]
  ]

  expect_s3_class(fixture$form, "effect_measurement_form")
  expect_identical(fixture$form$plan$query_role, "variation")
  expect_identical(fixture$form$plan$sampling_axis, "trial")
  expect_identical(fixture$form$plan$query_construction,
    "joint_covariance")
  expect_identical(fixture$form$plan$stages$reduction$operation$order,
    "aggregate_first")
  expect_lt(expected, 0)
  expect_equal(got, expected, tolerance = 1e-12)
  expect_s3_class(covariance_coupling(fixture$form),
    "effect_coupling_result")
})

test_that("a rank-one query remains effect coupling but not connectivity", {
  fixture <- public_measurement_fixture(rank_one = TRUE)
  effect <- effect_coupling(fixture$form)
  expect_identical(effect$kind, "effect_coupling")
  refusal <- catch_refusal(connectivity(fixture$form, "correlation"))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "nondegenerate_variation")
  expect_identical(refusal$namespace, "coupling_views")
  expect_identical(refusal$reasons, "rank_one_variation_axis")
  expect_match(conditionMessage(refusal), "effective rank 1")
  expect_error(canonical_coupling(fixture$form, ridge = 1e-4),
    class = "effect_capability_refusal")
})

test_that("public forms cross independent experimental and neural sides", {
  set.seed(2026081229)
  left_domain <- abstract_domain(3, id = "public:cross:left-neural:v1")
  right_domain <- abstract_domain(4, id = "public:cross:right-neural:v1")
  left_effects <- effect_space(c("encode_a", "encode_b"),
    basis_id = "public:encoding:v1")
  right_effects <- effect_space(c("retrieve_a", "retrieve_b", "lure"),
    basis_id = "public:retrieval:v1")
  left_b <- matrix(rnorm(6), 2L, 3L)
  right_b <- matrix(rnorm(12), 3L, 4L)
  left <- relation(list(encode = left_b), effects = left_effects,
    domain = left_domain)
  right <- relation(list(retrieve = right_b), effects = right_effects,
    domain = right_domain)
  left_frame <- measurement_frame(
    list(seed = matrix(c(1, 0, 1), 1L)), left_domain,
    id = "public:cross:left-frame:v1"
  )
  right_frame <- measurement_frame(
    list(target = matrix(c(1, 0, 0, 1, 0, 1, 0, 0), 2L, 4L)),
    right_domain, id = "public:cross:right-frame:v1"
  )
  between <- edge_frame(
    "seed", "target", left_frame, to_frame = right_frame
  )
  h <- matrix(c(1, -0.5, 0, 0, 0.5, -1), 2L, 3L, byrow = TRUE)
  by <- pair_query(h, left_effects, right_effects)
  over <- pairing("encode", "retrieve", directed = TRUE)
  form <- measurement_form(
    left, between, by, over, right = right, route = "pull_h"
  )
  effect <- effect_coupling(form)
  expected <- left_frame$legs$seed$operator %*%
    crossprod(left_b, h %*% right_b) %*%
    t(right_frame$legs$target$operator)

  expect_false(form$capabilities$self_form)
  expect_false(form$capabilities$symmetric)
  expect_equal(effect$values[[1L]], expected, tolerance = 1e-12)
  expect_error(connectivity(form, "canonical", ridge = 1e-4),
    "repeated variation")
})

test_that("public canonical, alignment, and Gaussian views stay capability gated", {
  fixture <- public_measurement_fixture(mode = "multivariate")
  canonical <- canonical_coupling(fixture$form, ridge = 0.05)
  alignment <- geometry_alignment(fixture$form)
  model <- gaussian_covariance_model(
    list(assumption = "joint Gaussian trial variation")
  )
  information <- connectivity(
    fixture$form, "gaussian_information",
    ridge = 0.05, model = model, units = "bits"
  )

  expect_identical(canonical$kind, "canonical_coupling")
  expect_identical(alignment$kind, "geometry_alignment")
  expect_identical(alignment$units, "linear_cka")
  expect_identical(information$kind, "gaussian_mutual_information")
  expect_identical(information$units, "bits")
  expect_error(connectivity(fixture$form, "canonical"), "regularization")
  expect_error(connectivity(
    fixture$form, "gaussian_information", ridge = 0.05
  ), "Gaussian model")
})

test_that("coherent and configuration components lift across an edge", {
  set.seed(2026081230)
  q <- 5L
  p <- 5L
  domain <- abstract_domain(p, id = "public:decomposition:neural:v1")
  effects <- effect_space(paste0("sample", seq_len(q)),
    basis_id = "public:decomposition:samples:v1")
  b <- matrix(rnorm(q * p), q, p)
  rel <- relation(list(run = b), effects = effects, domain = domain)
  additive <- additive_frame(matrix(c(
    1, 2, 1, 0, 0,
    0, 0, 1, 2, 1
  ), 2L, p, byrow = TRUE), domain = domain)
  frame <- measurement_frame(additive, mode = "coherent_configuration")
  between <- edge_frame(
    rep(frame$node_ids, each = 2L),
    rep(frame$node_ids, times = 2L), frame
  )
  h <- diag(q) - matrix(1 / q, q, q)
  form <- measurement_form(
    rel, between,
    variation_query(h, effects, "trial", "joint_covariance"),
    pairing("run", "run", directed = TRUE,
      self_pairs = "allow_biased", independence = "not_independent"),
    route = "pull_h"
  )
  edge <- which(form$block_index$left == frame$node_ids[[1L]] &
    form$block_index$right == frame$node_ids[[2L]])
  components <- measurement_components(form, edge)

  expect_identical(
    paste(components$left_component, components$right_component, sep = "::"),
    c("coherent::coherent", "configuration::coherent",
      "coherent::configuration", "configuration::configuration")
  )
  expect_identical(components$raw_entries_meaningful,
    c(TRUE, FALSE, FALSE, FALSE))
  expect_true(all(is.finite(components$frobenius_strength)))
  expect_true(all(components$frobenius_strength >= 0))
})

test_that("public reconstruction distinguishes exact and projected results", {
  set.seed(2026081231)
  q <- 6L
  p <- 3L
  domain <- abstract_domain(p, id = "public:tomography:neural:v1")
  effects <- effect_space(paste0("sample", seq_len(q)),
    basis_id = "public:tomography:samples:v1")
  b <- matrix(rnorm(q * p), q, p)
  rel <- relation(list(run = b), effects = effects, domain = domain)
  frame <- measurement_frame(list(
    first = diag(p)[1:2, , drop = FALSE],
    second = diag(p)[3, , drop = FALSE]
  ), domain, id = "public:tomography:parseval:v1")
  pairs <- expand.grid(from = frame$node_ids, to = frame$node_ids,
    stringsAsFactors = FALSE)
  between <- edge_frame(pairs$from, pairs$to, frame)
  h <- (diag(q) - matrix(1 / q, q, q)) / (q - 1)
  form <- measurement_form(
    rel, between,
    variation_query(h, effects, "trial", "joint_covariance"),
    pairing("run", "run", directed = TRUE,
      self_pairs = "allow_biased", independence = "not_independent"),
    route = "pull_h"
  )
  reference <- crossprod(b, h %*% b)
  result <- reconstruct_evidence(
    form, between, reference_operator = reference
  )

  expect_identical(result$method, "parseval")
  expect_identical(result$status,
    "numerically_certified_exact_reconstruction")
  expect_true(result$lossless)
  expect_equal(result$operator, reference, tolerance = 1e-12)
  expect_error(reconstruct_evidence(
    form, between, workspace_bytes = 1
  ), "exceeding")
})
