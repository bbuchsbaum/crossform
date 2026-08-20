# Multiscale searchlights -- `searchlights(radius = c(...))`, ticket D3.
#
# `searchlights()` with several radii is a request for a conservative frame
# family: one column-normalized member frame per radius, stacked by
# `frame_family()` under weights that sum to one. The law it exists to make
# reachable is `design/conservative-geometry-contract.md` section 3.1,
#
#     sum over the rows of scale s of G  =  alpha_s * G_Omega,
#
# which holds block by block whatever the data say. Two consequences are what
# this file pins:
#
#   * per-scale *energy* is the analyst's own weight vector and is never a
#     finding (section 3.1's normative consequence), and
#   * the coherent *share* of each block's fixed budget is invariant to those
#     weights, which is what makes the coherence spectrum -- the thing a
#     multiscale family is actually for -- well posed (section 3.2).
#
# The single-radius path is unchanged by the new argument, and the file says so
# against an independent oracle rather than against the implementation.

multiscale_domain <- function(id = "multiscale-volume") {
  volume_domain(array(TRUE, c(3L, 3L, 2L)), spacing = c(2, 2, 2), id = id)
}

multiscale_fixture <- function(seed = 91021L, q = 3L, id = "multiscale") {
  domain <- multiscale_domain(id)
  n <- domain$n_features
  set.seed(seed)
  mk <- function() {
    value <- matrix(stats::rnorm(q * n), q, n)
    rownames(value) <- letters[seq_len(q)]
    value
  }
  relation_value <- relation(list(run1 = mk(), run2 = mk()), domain = domain)
  list(domain = domain, relation = relation_value,
    over = cross_partitions(relation_value))
}

multiscale_total <- function(fixture, frame) {
  geometry_component(
    materialize_geometry(fixture$relation, frame, fixture$over), "total"
  )
}

# An independent Euclidean-ball membership operator, built from the domain's
# own coordinates rather than from the package's support index, so the
# single-radius assertions below do not check the implementation against
# itself.
oracle_searchlight_weights <- function(domain, radius, normalization) {
  distances <- as.matrix(stats::dist(domain$coordinates))
  membership <- matrix(as.numeric(distances <= radius), nrow(distances))
  if (identical(normalization, "local")) {
    membership <- membership / rowSums(membership)
  } else if (identical(normalization, "conservative")) {
    membership <- t(t(membership) / colSums(membership))
  }
  membership
}

# Construction ---------------------------------------------------------------

test_that("several radii compile to one conservative frame family", {
  domain <- multiscale_domain("multiscale-construct")
  specification <- searchlights(c(2.1, 4.1), "conservative")

  expect_s3_class(specification, "effect_frame_spec")
  expect_identical(specification$kind, "searchlight_family")
  expect_identical(specification$normalization, "conservative")
  expect_identical(specification$radius, c(2.1, 4.1))

  family <- compile_frame(specification, domain)

  expect_s3_class(family, "effect_frame")
  expect_identical(family$representation, "additive_diagonal")
  expect_identical(family$normalization, "conservative")
  expect_identical(nrow(family$weights), 2L * domain$n_features)
  expect_identical(ncol(family$weights), domain$n_features)

  # The family is exactly the stack an analyst could have assembled by hand
  # from two single-radius frames, which is the point: the constructor adds
  # the checking and the metadata, not a different operator.
  members <- lapply(c(2.1, 4.1), function(radius) {
    compile_frame(searchlights(radius, "conservative"), domain)
  })
  expect_equal(
    as.matrix(family$weights),
    rbind(0.5 * as.matrix(members[[1L]]$weights),
      0.5 * as.matrix(members[[2L]]$weights)),
    tolerance = 0, ignore_attr = TRUE
  )
})

test_that("family rows carry the radius they came from", {
  domain <- multiscale_domain("multiscale-metadata")
  family <- compile_frame(searchlights(c(2.1, 4.1), "conservative"), domain)
  index <- family$index

  expect_identical(names(index),
    c("measurement", "family", "node", "scale", "center", "alpha"))
  expect_identical(nrow(index), nrow(family$weights))
  expect_identical(unique(index$family), c("radius-2.1", "radius-4.1"))

  # `$index$scale` is each member's own radius, per row rather than frame
  # wide, which is what a stacked frame's `$specification` can no longer say.
  expect_identical(unique(index$scale[index$family == "radius-2.1"]), 2.1)
  expect_identical(unique(index$scale[index$family == "radius-4.1"]), 4.1)

  # The node label is the domain feature the neighborhood is centered on, and
  # the identifier that reaches a result is unique across scales even though
  # the node label is not.
  expect_identical(index$node[index$family == "radius-2.1"],
    as.character(domain$feature_ids))
  expect_true(any(duplicated(index$node)))
  expect_false(anyDuplicated(index$measurement) > 0L)
  expect_identical(index$measurement, paste0(index$family, "::", index$node))
  expect_identical(index$center, index$node)

  # The originating request survives on the stacked specification, which
  # otherwise remembers only the members it was handed.
  request <- family$specification$request
  expect_identical(request$kind, "searchlight_family")
  expect_identical(request$radius, c(2.1, 4.1))
  expect_identical(family$specification$kind, "frame_family")
  expect_identical(names(family$specification$members),
    c("radius-2.1", "radius-4.1"))
})

test_that("weights default to equal and are applied exactly as passed", {
  domain <- multiscale_domain("multiscale-weights")

  equal <- compile_frame(searchlights(c(2.1, 3.1, 4.1), "conservative"),
    domain)
  expect_identical(unique(equal$index$alpha), 1 / 3)
  expect_identical(equal$specification$alpha,
    stats::setNames(rep(1 / 3, 3L),
      c("radius-2.1", "radius-3.1", "radius-4.1")))

  weighted <- compile_frame(
    searchlights(c(2.1, 3.1, 4.1), "conservative",
      weights = c(0.2, 0.5, 0.3)),
    domain
  )
  alpha <- weighted$index$alpha
  expect_identical(unique(alpha[weighted$index$family == "radius-2.1"]), 0.2)
  expect_identical(unique(alpha[weighted$index$family == "radius-3.1"]), 0.5)
  expect_identical(unique(alpha[weighted$index$family == "radius-4.1"]), 0.3)

  # The weight the row records is the weight in the operator, not a label.
  rows <- which(weighted$index$family == "radius-3.1")
  member <- compile_frame(searchlights(3.1, "conservative"), domain)
  expect_equal(as.matrix(weighted$weights[rows, , drop = FALSE]),
    0.5 * as.matrix(member$weights), tolerance = 0, ignore_attr = TRUE)
})

test_that("named weights are matched to scale names, not to argument order", {
  domain <- multiscale_domain("multiscale-named")
  specification <- searchlights(c(4, 8), "conservative",
    weights = c(`radius-8` = 0.75, `radius-4` = 0.25))

  expect_identical(specification$weights,
    c(`radius-4` = 0.25, `radius-8` = 0.75))
  family <- compile_frame(specification, domain)
  expect_identical(unique(family$index$alpha[
    family$index$family == "radius-4"]), 0.25)
})

# One radius is what it always was ------------------------------------------

test_that("a single radius compiles to exactly the frame it did before", {
  domain <- multiscale_domain("multiscale-single")

  # The specification still has its three canonical fields, in order, with the
  # radius stored as the caller passed it.
  specification <- searchlights(2.1)
  expect_identical(names(unclass(specification)),
    c("kind", "normalization", "radius"))
  expect_identical(specification$kind, "searchlights")
  expect_identical(specification$normalization, "local")
  expect_identical(specification$radius, 2.1)
  expect_identical(typeof(searchlights(4L)$radius), "integer")

  for (normalization in c("none", "local", "conservative")) {
    frame <- compile_frame(searchlights(2.1, normalization), domain)

    expect_identical(frame$normalization, normalization)
    expect_identical(frame$specification$kind, "searchlights")
    expect_identical(names(frame$index), "measurement")
    expect_identical(frame$index$measurement, domain$feature_ids)
    expect_equal(as.matrix(frame$weights),
      oracle_searchlight_weights(domain, 2.1, normalization),
      tolerance = 1e-12, ignore_attr = TRUE)
  }

  # And the compiled object is the same object whichever way the arguments
  # are written, `weights` having no default that could leak into it.
  expect_identical(compile_frame(searchlights(2.1), domain),
    compile_frame(searchlights(radius = 2.1, normalization = "local"), domain))
})

# Refusals -------------------------------------------------------------------

test_that("only conservative normalization admits several radii", {
  expect_error(searchlights(c(4, 8, 12)),
    "needs `normalization = \"conservative\"`", class = "effect_input_error")
  expect_error(searchlights(c(4, 8, 12)),
    "conservative-geometry-contract", class = "effect_input_error")
  expect_error(searchlights(c(4, 8), "none"),
    "2 radii asked for \"none\"", class = "effect_input_error")
  expect_error(searchlights(c(4, 8), "local"),
    "no conserved budget for `weights` to divide",
    class = "effect_input_error")
})

test_that("multiscale weights are refused rather than repaired", {
  expect_error(searchlights(c(4, 8), "conservative", weights = c(1, 1, 1)),
    "must be 2 finite numbers, one per radius", class = "effect_input_error")
  expect_error(searchlights(c(4, 8), "conservative", weights = c(0.5, 0.2)),
    "must sum to one", class = "effect_input_error")
  expect_error(searchlights(c(4, 8), "conservative", weights = c(0.5, 0.2)),
    "never renormalized for you", class = "effect_input_error")
  expect_error(searchlights(c(4, 8), "conservative", weights = c(0, 1)),
    "must be positive", class = "effect_input_error")
  expect_error(
    searchlights(c(4, 8), "conservative", weights = c(a = 0.5, b = 0.5)),
    "must name every radius exactly once", class = "effect_input_error")
})

test_that("weights alongside a single radius are refused, not ignored", {
  # Silently dropping them would let a mistyped multiscale request succeed as
  # a single-scale one, which is the failure the metadata exists to prevent.
  expect_error(searchlights(4, "conservative", weights = 1),
    "applies only when `radius` names more than one scale",
    class = "effect_input_error")
  expect_error(searchlights(4, "conservative", weights = c(0.5, 0.5)),
    "received one radius and 2 weights", class = "effect_input_error")
})

test_that("radii that cannot be told apart are refused", {
  expect_error(searchlights(c(4, 4), "conservative"),
    "must name a distinct scale", class = "effect_input_error")
  expect_error(searchlights(c(3, -4), "conservative"),
    "several of them for a multiscale family", class = "effect_input_error")
  expect_error(searchlights(c(3, NA), "conservative"),
    "several of them for a multiscale family", class = "effect_input_error")
  expect_error(searchlights(numeric(0), "conservative"),
    "several of them for a multiscale family", class = "effect_input_error")
})

test_that("a noncanonical multiscale specification is refused at compile", {
  domain <- multiscale_domain("multiscale-noncanonical")
  specification <- searchlights(c(4, 8), "conservative")
  specification$weights <- unname(specification$weights)

  expect_error(compile_frame(specification, domain),
    "missing or noncanonical", class = "effect_input_error")
})

# The law (contract sections 3.1 and 3.2) ------------------------------------

test_that("a multiscale family conserves overall and scale by scale", {
  fixture <- multiscale_fixture(id = "multiscale-law")
  domain <- fixture$domain
  weights <- c(`radius-2.1` = 0.2, `radius-3.1` = 0.5, `radius-4.1` = 0.3)
  family <- compile_frame(
    searchlights(c(2.1, 3.1, 4.1), "conservative", weights = weights), domain
  )

  global <- drop(multiscale_total(fixture,
    compile_frame(whole_brain("none"), domain)))
  total <- multiscale_total(fixture, family)

  # Whole family: the rows sum to the global geometry, because sum(alpha) = 1.
  expect_equal(colSums(total), global, tolerance = 1e-12)
  expect_lt(max(abs(colSums(total) - global)), 1e-12)

  # Scale by scale: each block carries exactly its own weight of that same
  # budget. The rows are found through `$index$scale`, which is the metadata
  # this constructor exists to attach.
  worst <- 0
  for (radius in c(2.1, 3.1, 4.1)) {
    rows <- which(family$index$scale == radius)
    expect_identical(length(rows), domain$n_features)
    alpha <- weights[[paste0("radius-", radius)]]
    observed <- colSums(total[rows, , drop = FALSE])
    expect_equal(observed, alpha * global, tolerance = 1e-12)
    worst <- max(worst, max(abs(observed - alpha * global)))
  }
  expect_lt(worst, 1e-12)

  # `frame_conservation()` certifies the same thing without any geometry.
  report <- frame_conservation(family)
  expect_true(report$conserved)
  expect_identical(report$members$family, names(weights))
  expect_identical(report$members$alpha, unname(weights))
  expect_lt(max(report$members$max_deviation), 1e-12)
})

test_that("per-scale energy is the weight vector and the share is not", {
  # Section 3.1's normative consequence, and section 3.2's reason it is not
  # the end of the story: the energy column is alpha times a constant, while
  # the coherent share is invariant to alpha altogether.
  fixture <- multiscale_fixture(id = "multiscale-spectrum")
  domain <- fixture$domain
  radii <- c(2.1, 3.1, 4.1)
  contrast <- c(a = 1, b = -1, c = 0)

  spectrum <- function(weights) {
    family <- compile_frame(
      searchlights(radii, "conservative", weights = weights), domain
    )
    values <- as.data.frame(contrast_energy(
      plan_geometry(fixture$relation, family, fixture$over), contrast
    ))
    values$scale <- family$index$scale
    values
  }

  first <- spectrum(c(0.2, 0.5, 0.3))
  second <- spectrum(c(0.6, 0.1, 0.3))

  # Energy per scale follows the weights exactly: the same ratio, whichever
  # weights were chosen, so the panel reports the analyst, not the brain.
  energy <- function(values) vapply(split(values$total, values$scale), sum,
    numeric(1))
  budget <- sum(energy(first))
  expect_equal(unname(energy(first)), c(0.2, 0.5, 0.3) * budget,
    tolerance = 1e-12)
  expect_equal(unname(energy(second)), c(0.6, 0.1, 0.3) * budget,
    tolerance = 1e-12)

  # The coherent share is unchanged by that same reweighting, row for row.
  expect_equal(first$coherence_fraction, second$coherence_fraction,
    tolerance = 1e-10)

  # And it does vary with scale, which is what makes it worth reporting.
  share <- vapply(split(first$coherent, first$scale), sum, numeric(1)) /
    energy(first)
  expect_gt(diff(range(share)), 0.01)
})

# The neuroim2 provider ------------------------------------------------------

test_that("neuroim2_searchlights builds the same kind of family", {
  testthat::skip_if_not_installed("neuroim2", minimum_version = "0.19.0")
  values <- array(FALSE, c(5L, 5L, 4L))
  values[2:4, 2:4, 2:3] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(values,
    neuroim2::NeuroSpace(c(5L, 5L, 4L), spacing = c(3, 3, 3)))
  domain <- neuroim2_volume_domain(mask)

  family <- neuroim2_searchlights(mask, c(4, 6), domain = domain,
    normalization = "conservative", weights = c(0.4, 0.6))

  expect_s3_class(family, "effect_frame")
  expect_identical(family$normalization, "conservative")
  expect_identical(nrow(family$weights), 2L * domain$n_features)
  expect_identical(unique(family$index$family), c("radius-4", "radius-6"))
  expect_identical(unique(family$index$scale[
    family$index$family == "radius-6"]), 6)
  expect_identical(unique(family$index$alpha[
    family$index$family == "radius-4"]), 0.4)
  expect_identical(family$index$node[family$index$family == "radius-4"],
    as.character(domain$feature_ids))

  # The members are the same single-radius frames the adapter has always
  # built, only column normalized and weighted.
  member <- neuroim2_searchlights(mask, 6, domain = domain,
    normalization = "conservative")
  rows <- which(family$index$family == "radius-6")
  expect_equal(as.matrix(family$weights[rows, , drop = FALSE]),
    0.6 * as.matrix(member$weights), tolerance = 0, ignore_attr = TRUE)

  # The pinned upstream geometry travels with the request.
  request <- family$specification$request
  expect_identical(request$kind, "neuroim2_searchlight_family")
  expect_identical(request$radius, c(4, 6))
  expect_identical(request$upstream_commit, "77b1ddb")

  report <- frame_conservation(family)
  expect_true(report$conserved)
  expect_identical(report$members$alpha, c(0.4, 0.6))
  expect_lt(max(report$members$max_deviation), 1e-12)
})

test_that("the neuroim2 provider refuses what searchlights() refuses", {
  testthat::skip_if_not_installed("neuroim2", minimum_version = "0.19.0")
  values <- array(FALSE, c(5L, 5L, 4L))
  values[2:4, 2:4, 2:3] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(values,
    neuroim2::NeuroSpace(c(5L, 5L, 4L), spacing = c(3, 3, 3)))

  expect_error(neuroim2_searchlights(mask, c(4, 6)),
    "needs `normalization = \"conservative\"`", class = "effect_input_error")
  expect_error(neuroim2_searchlights(mask, c(4, 6),
    normalization = "conservative", weights = c(0.5, 0.2)),
    "must sum to one", class = "effect_input_error")
  expect_error(neuroim2_searchlights(mask, 4, weights = 1),
    "applies only when `radius` names more than one scale",
    class = "effect_input_error")

  # One radius is still one frame, with no family metadata attached.
  single <- neuroim2_searchlights(mask, 4, normalization = "none")
  expect_identical(names(single$index), "measurement")
  expect_identical(single$specification$kind, "neuroim2_searchlights")
})
