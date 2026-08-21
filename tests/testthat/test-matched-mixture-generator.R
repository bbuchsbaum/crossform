if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

mixture_generator_environment <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(testthat::test_path(
    "..", "..", "benchmarks", "matched-interpretability",
    "00-mixture-generator.R"
  ), envir = environment)
  environment
}

test_that("fixed-total mixtures obey the analytic component path", {
  generator <- mixture_generator_environment()
  theta <- c(0, pi / 6, pi / 4, pi / 3, pi / 2)
  fixtures <- generator$fixed_total_mixture_grid(
    theta, n_features = 12L, total_magnitude = 3.5, seed = 4101L
  )
  truth <- do.call(rbind, lapply(fixtures, function(x) unlist(x$truth[c(
    "total", "coherent", "configuration", "coherent_share",
    "configuration_share"
  )])))

  expect_equal(truth[, "total"], rep(3.5, length(theta)), tolerance = 1e-12)
  expect_equal(truth[, "coherent"], 3.5 * cos(theta)^2,
               tolerance = 1e-12)
  expect_equal(truth[, "configuration"], 3.5 * sin(theta)^2,
               tolerance = 1e-12)
  expect_equal(truth[, "coherent"] + truth[, "configuration"],
               truth[, "total"], tolerance = 1e-12)
  expect_equal(truth[, "coherent_share"] + truth[, "configuration_share"],
               rep(1, length(theta)), tolerance = 1e-12)
  expect_equal(unname(truth[1L, c("coherent_share", "configuration_share")]),
               c(1, 0), tolerance = 1e-12)
  expect_equal(unname(truth[length(theta),
                            c("coherent_share", "configuration_share")]),
               c(0, 1), tolerance = 1e-12)
})

test_that("mixture metadata is complete and relation ready", {
  generator <- mixture_generator_environment()
  fixture <- generator$fixed_total_mixture(
    n_features = 8L, total_magnitude = 2, theta = pi / 3, seed = 4102L
  )
  metadata <- fixture$metadata

  expect_s3_class(fixture, "crossform_matched_mixture")
  expect_identical(metadata$schema_version, "matched-mixture-v1")
  expect_match(metadata$basis_id, "p8:seed4102", fixed = TRUE)
  expect_identical(metadata$normalization$direction, "euclidean_unit")
  expect_identical(metadata$normalization$frame, "none")
  expect_identical(metadata$normalization$coherent_projector,
                   "11T/n_features")
  expect_identical(metadata$total_magnitude, 2)
  expect_identical(metadata$theta, pi / 3)
  expect_equal(metadata$theta_degrees, 60, tolerance = 1e-14)
  expect_identical(metadata$seed, 4102L)
  expect_length(metadata$expected_spectrum, 5L)
  expect_true(metadata$strict_basis)
  expect_equal(unname(metadata$basis$gram), diag(2), tolerance = 1e-12)
  expect_lt(abs(metadata$basis$configuration_mean), 1e-12)
  expect_equal(drop(fixture$contrast %*% fixture$effect_matrix),
               fixture$effect_pattern, tolerance = 1e-15)
  expect_equal(fixture$truth$neural_spectrum,
               c(2, rep(0, 7L)), tolerance = 1e-12)
})

test_that("seeded configuration is reproducible without changing caller RNG", {
  generator <- mixture_generator_environment()
  set.seed(991L)
  before <- .Random.seed
  first <- generator$fixed_total_mixture(n_features = 10L, seed = 4103L)
  expect_identical(.Random.seed, before)
  second <- generator$fixed_total_mixture(n_features = 10L, seed = 4103L)
  other <- generator$fixed_total_mixture(n_features = 10L, seed = 4104L)

  expect_identical(first, second)
  expect_false(isTRUE(all.equal(
    first$metadata$basis$configuration,
    other$metadata$basis$configuration
  )))
})

test_that("invalid bases refuse unless the exact relaxation is declared", {
  generator <- mixture_generator_environment()
  mean <- rep(0.5, 4L)
  configuration <- c(2, 0, 0, 0)
  expect_error(
    generator$fixed_total_mixture(
      n_features = 4L, mean_direction = mean,
      configuration_direction = configuration
    ),
    "norm 1"
  )

  configuration <- c(1, -1, 0, 0) / sqrt(2)
  nonorthogonal <- (mean + configuration) / sqrt(2)
  expect_error(
    generator$fixed_total_mixture(
      n_features = 4L, mean_direction = mean,
      configuration_direction = nonorthogonal
    ),
    "orthogonal"
  )

  noncanonical <- c(1, -1, 0, 0) / sqrt(2)
  configuration <- c(0, 0, 1, -1) / sqrt(2)
  expect_error(
    generator$fixed_total_mixture(
      n_features = 4L, mean_direction = noncanonical,
      configuration_direction = configuration
    ),
    "normalized constant mode"
  )

  relaxed <- generator$fixed_total_mixture(
    n_features = 4L, mean_direction = 2 * mean,
    configuration_direction = c(1, 0, 0, 0),
    allow_unnormalized = TRUE,
    allow_nonorthogonal = TRUE,
    allow_noncanonical_mean = TRUE
  )
  expect_false(relaxed$metadata$strict_basis)
  expect_true(all(relaxed$metadata$explicit_relaxations))
  expect_equal(relaxed$truth$coherent + relaxed$truth$configuration,
               relaxed$truth$total, tolerance = 1e-12)
})

test_that("noiseless public decomposition recovers generator truth", {
  generator <- mixture_generator_environment()
  fixture <- generator$fixed_total_mixture(
    n_features = 9L, total_magnitude = 2.75, theta = pi / 5,
    seed = 4105L
  )
  domain <- abstract_domain(9L, id = "matched-mixture-production")
  blocks <- list(run1 = fixture$effect_matrix, run2 = fixture$effect_matrix)
  rel <- relation(blocks, domain = domain)
  plan <- plan_geometry(
    rel,
    compile_frame(whole_brain(normalization = "none"), domain),
    cross_partitions(rel, independence = "independent")
  )
  observed <- contrast_energy(plan, fixture$contrast)

  expect_equal(unname(observed$total), fixture$truth$total, tolerance = 1e-12)
  expect_equal(unname(observed$coherent), fixture$truth$coherent,
               tolerance = 1e-12)
  expect_equal(unname(observed$configuration), fixture$truth$configuration,
               tolerance = 1e-12)
  expect_equal(observed$coherent + observed$configuration, observed$total,
               tolerance = 1e-12)
})
