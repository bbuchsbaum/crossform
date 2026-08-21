if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

test_that("contrast crossnobis RDM and fixed RSA share the randomized court", {
  oracle <- cg_oracle()
  seeds <- cg_property_seeds()
  for (position in seq_along(seeds)) {
    case <- cg_case(seeds[[position]], rank_deficient = position %% 2L == 0L)
    info <- paste(case$replay,
      "mapping=H=ccT; D(G)=uTGu; beta=(XtX)^-1 Xt D(G)")

    contrast <- if (case$rank_deficient) {
      evaluate_geometry(
        case$plan, bilinear_query(tcrossprod(case$contrast)),
        component = "total"
      )
    } else {
      contrast_energy(case$plan, case$contrast)
    }
    cross <- crossnobis(case$plan, case$contrast)
    distances <- rdm(case$plan)
    regression <- rsa(case$plan, models = case$models)
    expected_contrast <- as.numeric(crossprod(
      case$contrast, case$geometry %*% case$contrast
    ))

    comparison <- oracle$oracle_numeric_comparison(
      expected_contrast,
      if (case$rank_deficient) drop(contrast$values) else contrast$total,
      matrices = list(case$metric / case$d),
      operations = case$partitions^2 * case$q^2 * case$d^2
    )
    expect_true(comparison$pass, info = info)
    observed_contrast <- if (case$rank_deficient) {
      drop(contrast$values)
    } else {
      contrast$total
    }
    expect_equal(unname(cross$values), unname(observed_contrast),
      tolerance = comparison$bound,
      info = info)
    expect_equal(as.numeric(distances$values), as.numeric(case$rdm),
      tolerance = max(comparison$bound), info = info)
    expect_equal(as.numeric(regression$coefficients), unname(case$rsa),
      tolerance = 5e-11, info = info)
    if (case$rank_deficient) {
      expect_error(contrast_energy(case$plan, case$contrast),
        "requires an SPD metric", class = "effect_execution_error",
        info = info)
    } else {
      expect_equal(contrast$total,
        contrast$coherent + contrast$configuration,
        tolerance = 2e-12, info = info)
    }
  }
})

test_that("linearity and scale are production-route laws", {
  case <- cg_case(1301)
  h1 <- tcrossprod(case$contrast)
  h2 <- tcrossprod(case$contrast2)
  first <- contrast_energy(case$plan, case$contrast)$total
  second <- contrast_energy(case$plan, case$contrast2)$total
  combined <- drop(evaluate_geometry(
    case$plan, bilinear_query(h1 + h2)
  )$values)
  expect_equal(unname(combined), unname(first + second), tolerance = 2e-12,
    info = case$replay)

  factor <- -1.7
  scaled_blocks <- lapply(case$blocks, `*`, factor)
  scaled <- cg_plan(scaled_blocks, case$metric, "cg-scale")
  expect_equal(
    contrast_energy(scaled, case$contrast)$total,
    factor^2 * first, tolerance = 5e-12, info = case$replay
  )

  metric_factor <- 2.3
  metric_scaled <- cg_plan(case$blocks, metric_factor * case$metric,
    "cg-metric-scale")
  expect_equal(
    contrast_energy(metric_scaled, case$contrast)$total,
    metric_factor * first, tolerance = 5e-12, info = case$replay
  )
})

test_that("effect relabeling is permutation equivariant", {
  case <- cg_case(1401, rank_deficient = FALSE)
  permutation <- rev(seq_len(case$q))
  permuted_blocks <- lapply(case$blocks, function(value) {
    value[permutation, , drop = FALSE]
  })
  permuted <- cg_plan(permuted_blocks, case$metric, "cg-permutation")

  original_contrast <- contrast_energy(case$plan, case$contrast)$total
  permuted_contrast <- contrast_energy(permuted, case$contrast)$total
  expect_equal(permuted_contrast, original_contrast, tolerance = 2e-12,
    info = case$replay)

  original_rdm <- as.data.frame(rdm(case$plan))
  permuted_rdm <- as.data.frame(rdm(permuted))
  key <- function(x) paste(pmin(x$left, x$right), pmax(x$left, x$right),
    sep = "::")
  expect_equal(
    permuted_rdm$estimate[match(key(original_rdm), key(permuted_rdm))],
    original_rdm$estimate, tolerance = 2e-12, info = case$replay
  )

  original_rsa <- rsa(case$plan, models = case$models)$coefficients
  permuted_rsa <- rsa(permuted, models = case$models)$coefficients
  expect_equal(permuted_rsa, original_rsa, tolerance = 2e-12,
    info = case$replay)
})

test_that("the bounded and deep courts expose replayable failure context", {
  deep <- identical(Sys.getenv("CROSSFORM_DEEP_EQUIVALENCE"), "true")
  expect_identical(cg_property_seeds(), if (deep) 1201:1250 else 1201:1206)
  case <- cg_case(1501, rank_deficient = TRUE)
  expect_match(case$replay, "seed=1501", fixed = TRUE)
  expect_match(case$replay, "metric_rank=", fixed = TRUE)
  expect_match(case$replay, "kappa=", fixed = TRUE)
  if (deep) {
    expect_identical(length(cg_property_seeds()), 50L)
  } else {
    expect_lt(length(cg_property_seeds()), 10L)
  }
})
