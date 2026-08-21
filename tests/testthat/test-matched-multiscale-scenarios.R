if (!file.exists(testthat::test_path("..", "..", "design", "unification-contract.md")))
  testthat::skip("source-checkout scientific court")

matched_multiscale_environment <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(testthat::test_path(
    "..", "..", "benchmarks", "matched-interpretability",
    "00-mixture-generator.R"
  ), envir = environment)
  sys.source(testthat::test_path(
    "..", "..", "benchmarks", "matched-interpretability",
    "01-multiscale-scenarios.R"
  ), envir = environment)
  environment
}

test_that("every scale carries its declared share of one fixed total", {
  generator <- matched_multiscale_environment()
  bundle <- generator$matched_multiscale_scenarios()
  expected <- bundle$expected

  expect_identical(bundle$schema_version, "matched-multiscale-v1")
  expect_equal(as.numeric(tapply(expected$total, expected$scenario, sum)),
               rep(4, 3L), tolerance = 1e-12)
  expect_equal(
    expected$total,
    expected$alpha * bundle$metadata$total_magnitude,
    tolerance = 1e-12
  )
  expect_equal(expected$coherent + expected$configuration,
               expected$total, tolerance = 1e-12)
  expect_equal(expected$coherent_share + expected$configuration_share,
               rep(1, nrow(expected)), tolerance = 1e-12)
})

test_that("broad fine and mixed scenarios meet prespecified scale targets", {
  generator <- matched_multiscale_environment()
  bundle <- generator$matched_multiscale_scenarios()
  expected <- bundle$expected

  broad <- expected[expected$scenario == "broad_coherent", ]
  fine <- expected[expected$scenario == "fine_configuration", ]
  mixed <- expected[expected$scenario == "mixed_broad_fine", ]
  expect_equal(broad$coherent_share, rep(1, nrow(broad)), tolerance = 1e-12)
  expect_gt(fine$configuration_share[[2L]],
            bundle$metadata$fine_min_configuration_share_at_first_nonpoint)
  expect_identical(
    fine$scale[[which.max(fine$configuration_share)]],
    max(bundle$metadata$radii)
  )
  crossing <- mixed$scale[
    mixed$configuration_share >= bundle$metadata$mixed_transition_threshold
  ][[1L]]
  expect_identical(crossing, bundle$metadata$mixed_expected_first_crossing)
  expect_true(all(diff(fine$configuration_share) > 0))
  expect_true(all(diff(mixed$configuration_share) > 0))
})

test_that("public coherence spectra equal independent frame formulas", {
  generator <- matched_multiscale_environment()
  bundle <- generator$matched_multiscale_scenarios()
  for (scenario in names(bundle$scenarios)) {
    fixture <- generator$matched_multiscale_plan(bundle, scenario)
    observed <- as.data.frame(coherence_spectrum(
      fixture$plan, fixture$contrast
    ))
    expected <- bundle$expected[bundle$expected$scenario == scenario, ]
    expected <- expected[match(observed$scale, expected$scale), ]
    info <- paste("scenario", scenario)

    expect_equal(observed$total, expected$total, tolerance = 2e-12,
                 info = info)
    expect_equal(observed$coherent, expected$coherent, tolerance = 2e-12,
                 info = info)
    expect_equal(observed$configuration, expected$configuration,
                 tolerance = 2e-12, info = info)
    expect_equal(observed$coherence_fraction, expected$coherent_share,
                 tolerance = 2e-12, info = info)
  }
})

test_that("frame order and relabeling do not change scale values", {
  generator <- matched_multiscale_environment()
  bundle <- generator$matched_multiscale_scenarios()
  original <- generator$matched_multiscale_plan(bundle, "mixed_broad_fine")
  original_spectrum <- as.data.frame(coherence_spectrum(
    original$plan, original$contrast
  ))

  order <- rev(seq_along(bundle$frames))
  reordered_frames <- bundle$frames[order]
  names(reordered_frames) <- paste0("renamed-", seq_along(reordered_frames))
  reordered_alpha <- bundle$metadata$alpha[order]
  names(reordered_alpha) <- names(reordered_frames)
  reordered_family <- do.call(
    frame_family, c(reordered_frames, list(alpha = reordered_alpha))
  )
  reordered <- generator$matched_multiscale_plan(
    bundle, "mixed_broad_fine", family = reordered_family,
    id_suffix = "reordered"
  )
  reordered_spectrum <- as.data.frame(coherence_spectrum(
    reordered$plan, reordered$contrast
  ))
  original_spectrum <- original_spectrum[order(original_spectrum$scale), ]
  reordered_spectrum <- reordered_spectrum[order(reordered_spectrum$scale), ]

  expect_equal(reordered_spectrum$total, original_spectrum$total,
               tolerance = 1e-12)
  expect_equal(reordered_spectrum$coherent, original_spectrum$coherent,
               tolerance = 1e-12)
  expect_equal(reordered_spectrum$configuration,
               original_spectrum$configuration, tolerance = 1e-12)
  expect_equal(reordered_spectrum$coherence_fraction,
               original_spectrum$coherence_fraction, tolerance = 1e-12)
})

test_that("boundary overlap and sparse-support cases are explicit", {
  generator <- matched_multiscale_environment()
  bundle <- generator$matched_multiscale_scenarios()
  structure <- unique(bundle$expected[c(
    "family", "scale", "n_nodes", "min_support", "max_support",
    "boundary_nodes", "overlapping_features", "sparse_weights"
  )])
  structure <- structure[order(structure$scale), ]

  expect_true(all(structure$sparse_weights))
  expect_true(all(structure$n_nodes == bundle$metadata$n_features))
  expect_identical(structure$min_support[[1L]], 1)
  expect_identical(structure$max_support[[1L]], 1)
  expect_true(all(structure$boundary_nodes[-1L] > 0L))
  expect_true(all(structure$min_support[-1L] < structure$max_support[-1L]))
  expect_true(all(structure$overlapping_features[-1L] ==
                    bundle$metadata$n_features))
  expect_true(frame_conservation(bundle$family)$conserved)
})
