if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

test_that("near-singular PSD total geometry agrees at high contrast leverage", {
  oracle <- cg_oracle()
  set.seed(2101)
  q <- 4L
  d <- 4L
  effects <- paste0("e", seq_len(q))
  blocks <- stats::setNames(lapply(1:3, function(run) {
    matrix(stats::rnorm(q * d), q, d, dimnames = list(effects, NULL))
  }), paste0("run", 1:3))
  metric <- diag(c(1, 1e-4, 1e-8, 1e-12))
  plan <- cg_plan(blocks, metric, "cg-near-singular")
  gamma <- matrix(1 / 6, 3L, 3L,
    dimnames = list(names(blocks), names(blocks)))
  diag(gamma) <- 0
  contrast <- c(e1 = 1e6, e2 = -1e6, e3 = 3, e4 = -3)
  expected <- oracle$oracle_energy(
    tcrossprod(contrast), blocks, gamma, metric / d
  )
  observed <- drop(evaluate_geometry(
    plan, bilinear_query(tcrossprod(contrast)), component = "total"
  )$values)
  court <- oracle$oracle_numeric_comparison(
    expected, observed, matrices = list(metric / d), operations = 1728L
  )

  expect_true(court$pass,
    info = sprintf("kappa=%g error=%g bound=%g",
      court$condition, court$error, court$bound))
  expect_equal(unname(crossnobis(plan, contrast)$values), unname(observed),
    tolerance = court$bound)
  expect_error(contrast_energy(plan, contrast), "requires an SPD metric",
    class = "effect_execution_error")
})

test_that("unequal partition weights agree with the explicit Gamma", {
  oracle <- cg_oracle()
  set.seed(2201)
  q <- 3L
  d <- 3L
  effects <- paste0("e", 1:q)
  runs <- paste0("run", 1:4)
  blocks <- stats::setNames(lapply(runs, function(run) {
    matrix(stats::rnorm(q * d), q, d, dimnames = list(effects, NULL))
  }), runs)
  edges <- utils::combn(runs, 2L)
  over <- pairing(edges[1L, ], edges[2L, ],
    weight = c(1, 2, 5, 11, 17, 23), directed = FALSE,
    independence = "independent")
  domain <- abstract_domain(d, id = "cg-unbalanced")
  relation <- relation(blocks, domain = domain)
  metric <- diag(c(2, 0.5, 3))
  plan <- plan_geometry(relation, compile_frame(whole_brain(), domain), over,
    metric = noise_precision(metric, domain))
  gamma <- matrix(0, 4L, 4L, dimnames = list(runs, runs))
  for (edge in seq_len(nrow(over))) {
    left <- match(over$left[[edge]], runs)
    right <- match(over$right[[edge]], runs)
    gamma[left, right] <- over$weight[[edge]] / 2
    gamma[right, left] <- over$weight[[edge]] / 2
  }
  contrast <- c(2, -3, 1)
  expected <- oracle$oracle_energy(
    tcrossprod(contrast), blocks, gamma, metric / d
  )
  observed <- contrast_energy(plan, contrast)$total
  expect_equal(observed, expected, tolerance = 2e-12)
  expect_identical(sum(gamma), 1)
})

test_that("zero effects and extreme but finite scales are deterministic", {
  effects <- c("a", "b", "c")
  zero <- stats::setNames(rep(list(matrix(0, 3L, 3L,
    dimnames = list(effects, NULL))), 3L), paste0("run", 1:3))
  zero_plan <- cg_plan(zero, diag(3), "cg-zero")
  expect_identical(contrast_energy(zero_plan, c(1, -1, 0))$total, 0)
  expect_true(all(rdm(zero_plan)$values == 0))

  set.seed(2301)
  base <- stats::setNames(lapply(1:2, function(run) {
    matrix(stats::rnorm(9), 3L, 3L, dimnames = list(effects, NULL))
  }), c("run1", "run2"))
  scale <- 1e80
  scaled <- lapply(base, `*`, scale)
  plan <- cg_plan(scaled, diag(3) / scale, "cg-extreme")
  value <- drop(evaluate_geometry(plan,
    bilinear_query(tcrossprod(c(1, -1, 0))), component = "total")$values)
  expect_true(is.finite(value))
  expect_equal(unname(value),
    unname(scale * contrast_energy(cg_plan(base, diag(3), "cg-base"),
      c(1, -1, 0))$total),
    tolerance = abs(value) * 2e-12)
})

test_that("invalid cells axes normalization and RSA rank fail early", {
  case <- cg_case(2401)
  broken <- case$blocks
  broken[[1L]][1, 1] <- NA_real_
  expect_error(relation(broken, domain = abstract_domain(case$d)),
    "finite", class = "effect_input_error")

  wrong_metric <- diag(case$d + 1L)
  metric_domain <- abstract_domain(case$d,
    feature_ids = paste0("f", seq_len(case$d)), id = "wrong-metric")
  expect_error(noise_precision(wrong_metric, metric_domain,
    support = paste0("f", seq_len(case$d + 1L))),
    "support|match|domain", class = "effect_input_error")

  expect_error(evaluate_geometry(case$plan,
    bilinear_query(diag(case$q + 1L))),
    "effect dimension", class = "effect_input_error")

  normalized <- catch_refusal(rdm(case$plan, normalize = "correlation"))
  expect_identical(normalized$capability, "guaranteed_psd")
  expect_identical(normalized$reasons,
    "signed_cross_generalized_diagonals")

  constant <- matrix(1, case$q, case$q,
    dimnames = list(case$effects, case$effects))
  diag(constant) <- 0
  expect_error(rsa(case$plan, models = list(constant = constant)),
    "rank deficient", class = "effect_input_error")
})

test_that("condition-aware courts do not accept visible perturbations", {
  oracle <- cg_oracle()
  reference <- c(1e-12, 1, 1e12)
  roundoff <- reference + c(1e-15, 1e-13, 1e-3)
  stable <- oracle$oracle_numeric_comparison(reference, roundoff,
    matrices = list(diag(c(1, 1e-4))), operations = 10L)
  expect_true(stable$pass)

  perturbed <- reference
  perturbed[[2L]] <- perturbed[[2L]] + 1e-3
  detected <- oracle$oracle_numeric_comparison(reference, perturbed,
    matrices = list(diag(2)), operations = 10L)
  expect_false(detected$pass)
  expect_gt(detected$error[[2L]], detected$bound[[2L]])
})
