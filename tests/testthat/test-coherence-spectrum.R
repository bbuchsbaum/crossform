# `coherence_spectrum()` -- coherent share versus scale, ticket D5.
#
# `design/conservative-geometry-contract.md` is the source of every assertion
# here:
#
#   section 3.1  per scale, sum_{x in s} G_{s,x} = alpha_s G_Omega. Per-scale
#                ENERGY is therefore the analyst's own weight vector times a
#                constant, and the contract makes it normative that no panel
#                of it may be presented as evidence about spatial scale.
#   section 3.2  the coherent SHARE is what is left when that is taken away,
#                and the reviewer's amendment strengthens the claim from "not
#                proportional to alpha" to *exactly invariant* to it. That is
#                what makes the spectrum well-posed, and it is the test this
#                file exists for.
#   section 4    coherent does not conserve, so a coherent energy is
#                frame-relative, and a share is masked, never clamped.
#   section 11.4 gap G3 -- the share is a function of (location, scale) and
#                not a number, so a location-wise collapse must be a declared
#                reduction. None is offered; `by_location = TRUE` returns the
#                table and stops.
#
# The per-node laws underneath are already covered by
# `test-conservative-geometry-contract.R`, `test-frame-family.R` and
# `test-multiscale-searchlights.R`, and this file does not restate them. What
# is new is the reduction: the grouping, the aggregated share, and the
# alpha-invariance the aggregation is supposed to expose.

spectrum_domain <- function(n = 15L, id = "d5") {
  abstract_domain(n, coordinates = cbind(seq_len(n) - 1, 0),
    feature_ids = paste0("v", seq_len(n)), id = id)
}

# One effect carries the whole contrast, the other is flat zero, so
# `c(face = 1, house = -1)` reads exactly the spatial pattern handed in. That
# makes every number below predictable in closed form rather than merely
# reproducible, which is what lets the point fixture assert 1, 1/3 and 1/5.
spectrum_pattern <- function(signal) {
  rbind(face = signal, house = rep(0, length(signal)))
}

spectrum_plan <- function(run1, run2 = run1, radii = c(0.5, 1.01, 2.01),
                          weights = NULL, normalization = "conservative",
                          id = "d5") {
  domain <- spectrum_domain(length(run1), id)
  relation_value <- relation(
    list(run1 = spectrum_pattern(run1), run2 = spectrum_pattern(run2)),
    domain = domain
  )
  frame <- if (length(radii) == 1L) {
    compile_frame(searchlights(radii, normalization), domain)
  } else {
    compile_frame(searchlights(radii, normalization, weights = weights), domain)
  }
  list(
    domain = domain,
    relation = relation_value,
    frame = frame,
    plan = plan_geometry(relation_value, frame,
      cross_partitions(relation_value, independence = "independent"))
  )
}

spectrum_contrast <- c(face = 1, house = -1)

spectrum_point_signal <- function(n = 15L, at = 8L) replace(rep(0, n), at, 1)

# A signal that is the same vector everywhere, with a small smooth ripple so
# the shares are not the degenerate exact 1 a perfectly constant field gives.
spectrum_smooth_signal <- function(n = 15L, phase = 0) {
  1 + 0.10 * sin(seq_len(n) / 3 + phase)
}

# Shape and the two routes ---------------------------------------------------

test_that("a spectrum is one row per scale carrying the weight that fixed it", {
  fixture <- spectrum_plan(spectrum_point_signal(), id = "d5-shape")
  spectrum <- coherence_spectrum(fixture$plan, spectrum_contrast)

  # The record is `contribution()`'s, because the spectrum is that
  # aggregation under a scale-resolved grouping. No new sealed class.
  expect_s3_class(spectrum, "effect_contrast_view")
  expect_identical(length(spectrum$total), 3L)

  table <- as.data.frame(spectrum)
  expect_identical(as.character(table$measurement), c("0.5", "1.01", "2.01"))
  expect_identical(table$scale, c(0.5, 1.01, 2.01))
  expect_identical(table$n_rows, rep(15L, 3L))
  expect_identical(sum(table$n_rows), length(fixture$frame$index$measurement))

  # `family` and `alpha` are carried through because they take one value
  # inside every group; `center` is not, because it does not.
  expect_identical(table$family,
    c("radius-0.5", "radius-1.01", "radius-2.01"))
  expect_equal(table$alpha, rep(1 / 3, 3L), tolerance = 1e-12)
  expect_false("center" %in% names(table))

  # A signed marginal is a local weighted mean, so it is masked exactly as
  # `contribution()` masks it rather than summed into a density times three.
  expect_true(all(is.na(spectrum$signed)))
  expect_identical(spectrum$metadata$aggregation$masked, "signed")
})

test_that("the plan route and the evaluated-view route agree exactly", {
  fixture <- spectrum_plan(spectrum_smooth_signal(), id = "d5-routes")
  from_plan <- coherence_spectrum(fixture$plan, spectrum_contrast)

  view <- contrast_energy(fixture$plan, spectrum_contrast)
  from_view <- coherence_spectrum(view, using = fixture$frame$index)

  # Bit-identical, not merely close: the view route aggregates the same
  # numbers the plan route evaluates, and the derived identity is a function
  # of the parent identity and the grouping, which are the same on both.
  expect_identical(from_view$total, from_plan$total)
  expect_identical(from_view$coherent, from_plan$coherent)
  expect_identical(from_view$coherence_fraction, from_plan$coherence_fraction)
  expect_identical(from_view$receipt$scientific_plan_id,
    from_plan$receipt$scientific_plan_id)

  # A result's `$index` is the identifier vector alone, so the view route
  # cannot recover the scales by itself and says so rather than guessing.
  bare <- try(coherence_spectrum(view), silent = TRUE)
  expect_s3_class(bare, "try-error")
  expect_match(conditionMessage(attr(bare, "condition")),
    "frame family", fixed = TRUE)
})

# (a) A point effect: the share falls as the neighborhood outgrows it --------

test_that("a single-voxel effect loses coherent share as the radius grows", {
  # One voxel carries the contrast. At a scale that sees only that voxel the
  # node is a singleton, which has exactly zero configuration
  # (`effect-form-v1` section 7), so the share is 1. At a wider scale the same
  # unit of evidence is spread over nodes that also see silent voxels, and the
  # common mode those nodes share is a smaller part of what they hold.
  fixture <- spectrum_plan(spectrum_point_signal(), id = "d5-point")
  spectrum <- coherence_spectrum(fixture$plan, spectrum_contrast)
  share <- spectrum$coherence_fraction

  expect_true(all(spectrum$coherence_fraction_valid))

  # Closed form. A radius covering k features gives every node containing the
  # signal voxel weight alpha/k on it and total frame mass alpha, so
  # phi_s = sum_x w_x8^2 / (sum_v w_xv) / alpha_s = k (alpha/k)^2 / alpha^2 =
  # 1/k. The three radii cover 1, 3 and 5 features.
  expect_equal(share, c(1, 1 / 3, 1 / 5), tolerance = 1e-12)

  # The assertion the closed form exists to make: monotone decrease over three
  # scales, strictly, with room to spare.
  expect_true(all(diff(share) < 0))
  expect_gt(share[[1L]] - share[[3L]], 0.5)

  # And the column that does NOT vary: each scale carries exactly its weight
  # of the whole-domain total, which is the law that makes an energy panel a
  # picture of `weights` (contract section 3.1).
  global <- contrast_energy(
    plan_geometry(fixture$relation,
      compile_frame(whole_brain("none"), fixture$domain),
      cross_partitions(fixture$relation, independence = "independent")),
    spectrum_contrast
  )
  expect_equal(spectrum$total, rep(global$total / 3, 3L), tolerance = 1e-12)
  expect_equal(sum(spectrum$total), global$total, tolerance = 1e-12)
})

# (b) A spatially smooth effect: the share stays high ------------------------

test_that("a spatially coherent effect keeps its share across scales", {
  # The same signal vector at every voxel, up to a small smooth ripple. A
  # node's weighted common mode is then almost the whole of what the node
  # holds, at any radius, so the share stays near one instead of falling.
  #
  # The radii are all larger than one feature deliberately. A singleton scale
  # has exactly zero configuration, so its aggregated configuration sits on
  # the mask boundary and is pinned separately below; this test is about the
  # scales where the split is a real quantity.
  fixture <- spectrum_plan(
    spectrum_smooth_signal(), spectrum_smooth_signal(phase = 0.4),
    radii = c(1.01, 2.01, 3.01), id = "d5-smooth"
  )
  spectrum <- coherence_spectrum(fixture$plan, spectrum_contrast)
  share <- spectrum$coherence_fraction

  expect_true(all(spectrum$coherence_fraction_valid))
  expect_true(all(share > 0.99))

  # No strong decrease: the point fixture above loses more than half its
  # share over three scales, and this one loses less than a hundredth of that.
  expect_lt(share[[1L]] - share[[3L]], 0.005)

  # The comparison is the point of the pair, so it is asserted rather than
  # left to a reader of two test names.
  point <- coherence_spectrum(
    spectrum_plan(spectrum_point_signal(), radii = c(1.01, 2.01, 3.01),
      id = "d5-smooth-point")$plan,
    spectrum_contrast
  )
  expect_true(all(share > point$coherence_fraction))
})

# (c) Alpha invariance -------------------------------------------------------

test_that("two families differing only in alpha give identical shares", {
  # Contract section 3.2, as strengthened by the 2026-08-20 review: both
  # components are homogeneous of degree one under a row rescaling, so the
  # ratio cancels alpha exactly. This is what lets a spectrum be reported
  # without disclosing the weights.
  #
  # Non-singleton radii, for the reason given in the smooth test: the exact
  # zero of a singleton scale's configuration is a separate question and is
  # pinned in its own test below.
  radii <- c(1.01, 2.01, 3.01)
  signal <- spectrum_smooth_signal()
  equal <- coherence_spectrum(
    spectrum_plan(signal, spectrum_smooth_signal(phase = 0.4), radii = radii,
      id = "d5-alpha-equal")$plan,
    spectrum_contrast
  )
  skewed <- coherence_spectrum(
    spectrum_plan(signal, spectrum_smooth_signal(phase = 0.4), radii = radii,
      weights = c(0.2, 0.5, 0.3), id = "d5-alpha-skew")$plan,
    spectrum_contrast
  )

  expect_equal(equal$coherence_fraction, skewed$coherence_fraction,
    tolerance = 1e-14)
  expect_lt(max(abs(equal$coherence_fraction - skewed$coherence_fraction)),
    1e-14)
  expect_identical(equal$coherence_fraction_valid,
    skewed$coherence_fraction_valid)

  # The energy columns must differ, or the invariance would be vacuous: they
  # are the two weight vectors times one constant.
  expect_equal(equal$total, rep(equal$total[[1L]], 3L), tolerance = 1e-12)
  expect_equal(skewed$total / sum(skewed$total), c(0.2, 0.5, 0.3),
    tolerance = 1e-12)
  expect_gt(max(abs(equal$total - skewed$total)), 0.1)

  # Coherent energy moves with alpha too -- it is degree one, not invariant --
  # so it is the ratio and only the ratio that is the finding.
  expect_gt(max(abs(equal$coherent - skewed$coherent)), 0.1)
  expect_equal(skewed$coherent / equal$coherent, c(0.2, 0.5, 0.3) * 3,
    tolerance = 1e-12)

  # The provenance says which column is which, so a reader of the object does
  # not have to have read the contract to know.
  expect_identical(skewed$metadata$aggregation$alpha_invariant,
    "coherence_fraction")
  expect_identical(skewed$metadata$aggregation$alpha_fixed, "total")
  expect_equal(unname(skewed$metadata$aggregation$alpha), c(0.2, 0.5, 0.3),
    tolerance = 1e-12)
})

# (d) Masking ----------------------------------------------------------------

test_that("a scale whose coherent energy is negative reports NA, not a clamp", {
  # Pointwise the two runs agree in sign, so every voxel's crossvalidated
  # contribution is positive and the whole-domain total is positive at every
  # scale. But inside a three-voxel window the two runs' weighted common modes
  # have opposite signs -- (1, 1, -3) sums to -1 while (3, 3, -1) sums to 5 --
  # so the coherent part of a wide node is negative. That is a real property
  # of a signed cross-generalized estimate, not a numerical accident, and the
  # share is undefined there: the components are not a nonnegative partition.
  n <- 15L
  fixture <- spectrum_plan(
    rep(c(1, 1, -3), length.out = n), rep(c(3, 3, -1), length.out = n),
    id = "d5-mask"
  )
  spectrum <- coherence_spectrum(fixture$plan, spectrum_contrast)

  expect_identical(spectrum$coherence_fraction_valid, c(TRUE, FALSE, FALSE))
  expect_identical(is.na(spectrum$coherence_fraction), c(FALSE, TRUE, TRUE))

  # Masked, never clamped: the offending components are reported as they are,
  # negative, rather than pushed to zero to make a fraction reportable.
  expect_true(all(spectrum$coherent[2:3] < -1))
  expect_true(all(spectrum$total > 0))
  expect_true(all(spectrum$configuration[2:3] > spectrum$total[2:3]))

  # The total is still budget-exact where the share is not defined: the two
  # laws are independent, and a masked share does not withdraw the ledger.
  expect_equal(spectrum$total, rep(spectrum$total[[1L]], 3L),
    tolerance = 1e-12)

  # The fraction is recomputed from the aggregated components, never averaged.
  # Averaging the per-node fractions of the widest scale would report a
  # number here, which is the error the mask exists to prevent.
  view <- contrast_energy(fixture$plan, spectrum_contrast)
  widest <- fixture$frame$index$scale == 2.01
  averaged <- mean(view$coherence_fraction[widest], na.rm = TRUE)
  expect_true(is.finite(averaged))
  expect_true(is.na(spectrum$coherence_fraction[[3L]]))
})

test_that("a singleton scale sits on the mask boundary and is never clamped", {
  # A singleton positive-mass frame has exactly zero configuration
  # (`effect-form-v1` section 7), so a point scale's aggregated configuration
  # is an exact zero in exact arithmetic and lands within one unit in the last
  # place of it in floating point. Which side it lands on decides whether the
  # mask fires, so the reported value is either exactly 1 or `NA` -- and never
  # a clamped number in between, which is the property that actually matters.
  #
  # This is characterised rather than repaired: the guard is the contract's
  # (section 4, normative 3) and is shared with every per-node coherence
  # fraction in the package, so relaxing it here would be a silent divergence
  # from the rule `contribution()` and `contrast_energy()` both apply.
  for (weights in list(NULL, c(0.2, 0.5, 0.3), c(0.5, 0.25, 0.25))) {
    spectrum <- coherence_spectrum(
      spectrum_plan(spectrum_smooth_signal(),
        spectrum_smooth_signal(phase = 0.4), weights = weights,
        id = "d5-singleton")$plan,
      spectrum_contrast
    )
    label <- paste(if (is.null(weights)) "equal" else weights, collapse = ",")
    expect_lt(abs(spectrum$configuration[[1L]]), 1e-14)
    reported <- spectrum$coherence_fraction[[1L]]
    expect_true(is.na(reported) || identical(reported, 1), info = label)
  }
})

# (e) The (center, scale) table ----------------------------------------------

test_that("by_location returns the (scale, center) table and collapses nothing", {
  # Gap G3: the coherent share is a function of (location, scale), not a
  # number, because a location sits in several scales at once.
  fixture <- spectrum_plan(spectrum_point_signal(), id = "d5-location")
  located <- coherence_spectrum(fixture$plan, spectrum_contrast,
    by_location = TRUE)

  expect_s3_class(located, "effect_contrast_view")
  expect_identical(length(located$total), 45L)

  table <- as.data.frame(located)
  expect_true(all(c("scale", "center") %in% names(table)))
  expect_identical(sort(unique(table$scale)), c(0.5, 1.01, 2.01))
  expect_identical(length(unique(table$center)), 15L)
  expect_identical(table$n_rows, rep(1L, 45L))
  expect_identical(as.character(table$measurement[[1L]]),
    paste0(table$scale[[1L]], "::", table$center[[1L]]))

  # `by = c("scale", "center")` is the same request spelled out, and asking
  # for both spellings at once is an error rather than a precedence rule.
  spelled <- coherence_spectrum(fixture$plan, spectrum_contrast,
    by = c("scale", "center"))
  expect_identical(spelled$coherence_fraction, located$coherence_fraction)
  expect_error(
    coherence_spectrum(fixture$plan, spectrum_contrast, by_location = TRUE,
      by = "scale"),
    "same choice twice"
  )

  # A location's scale profile is a filter on this table, not a reduction.
  # The signal voxel's profile is the per-scale spectrum of the point fixture,
  # because at that location it is the only node carrying anything.
  profile <- table[table$center == "v8", ]
  profile <- profile[order(profile$scale), ]
  expect_identical(nrow(profile), 3L)
  expect_equal(profile$coherence_fraction, c(1, 1 / 3, 1 / 5),
    tolerance = 1e-12)

  # Nothing collapses the table for you, and the record says so: an
  # alpha-weighted mean over scales or an argmax-scale map would be a declared
  # reduction with its own certificate, and neither is offered.
  expect_identical(located$metadata$aggregation$location_collapse, "none")
  expect_true(located$metadata$aggregation$by_location)
  expect_identical(located$metadata$aggregation$resolved_by,
    c("scale", "center"))

  # Summing the located table by scale reproduces the per-scale spectrum
  # exactly, which is what makes the two tables one object read at two
  # granularities rather than two computations.
  by_scale <- coherence_spectrum(fixture$plan, spectrum_contrast)
  regrouped <- vapply(split(located$total, table$scale), sum, numeric(1))
  expect_equal(unname(regrouped[order(as.numeric(names(regrouped)))]),
    by_scale$total, tolerance = 1e-12)
})

# (f) Provenance -------------------------------------------------------------

test_that("the spectrum records what it is and what it may not be read as", {
  fixture <- spectrum_plan(spectrum_smooth_signal(),
    spectrum_smooth_signal(phase = 0.4), weights = c(0.2, 0.5, 0.3),
    id = "d5-provenance")
  spectrum <- coherence_spectrum(fixture$plan, spectrum_contrast)
  record <- spectrum$metadata$aggregation

  expect_identical(record$reduction, "coherence_spectrum")
  expect_identical(record$aggregated_by, "scale")
  expect_identical(record$resolved_by, "scale")
  expect_identical(record$groups, 3L)
  expect_identical(record$measurements, 45L)

  # Contract section 4: a coherent energy is a share of this family's own
  # coherent mass, and two families give two incomparable denominators.
  expect_true(record$frame_relative)
  expect_identical(record$frame_relative_components,
    c("coherent", "configuration"))
  expect_identical(record$budget_exact, "total")
  expect_identical(record$frame_normalization, "conservative")
  expect_false(record$overlap_split)

  # The family and its weights, so the fixed column can be read as the weight
  # vector it is without going back to the frame.
  expect_identical(unname(record$family),
    c("radius-0.5", "radius-1.01", "radius-2.01"))
  expect_identical(names(record$alpha), c("0.5", "1.01", "2.01"))
  expect_equal(unname(record$alpha), c(0.2, 0.5, 0.3), tolerance = 1e-12)

  # The receipt is the parent's under a derived identity: same run, different
  # scientific object.
  parent <- contrast_energy(fixture$plan, spectrum_contrast)
  expect_false(identical(spectrum$receipt$scientific_plan_id,
    parent$receipt$scientific_plan_id))
  expect_identical(record$frame_normalization, "conservative")
  expect_identical(spectrum$metadata$aggregated_from,
    parent$receipt$scientific_plan_id)

  # The grouping enters the derived identity, so the per-scale and the
  # located readings of one run are distinguishable objects.
  located <- coherence_spectrum(fixture$plan, spectrum_contrast,
    by_location = TRUE)
  expect_false(identical(spectrum$receipt$scientific_plan_id,
    located$receipt$scientific_plan_id))

  # The print says which column is the finding, because the table alone does
  # not, and reading the fixed column as a result is the specific error the
  # contract forbids.
  output <- gsub("\\s+", " ",
    paste(utils::capture.output(print(spectrum)), collapse = " "))
  expect_match(output, "coherence_spectrum", fixed = TRUE)
  expect_match(output, "not a finding", fixed = TRUE)
  expect_match(output, "frame_relative: TRUE", fixed = TRUE)
})

# Refusals -------------------------------------------------------------------

test_that("a locally normalized frame is refused before anything is grouped", {
  # Section 1.1: a local frame is a detection map, its overlapping nodes
  # double-count, and adding them up estimates nothing. That is true of a
  # scale group exactly as it is of a region, so the refusal is the same one
  # `contribution()` raises and it fires before the grouping is even resolved.
  domain <- spectrum_domain(id = "d5-local")
  relation_value <- relation(
    list(run1 = spectrum_pattern(spectrum_smooth_signal()),
      run2 = spectrum_pattern(spectrum_smooth_signal(phase = 0.4))),
    domain = domain
  )
  frame <- compile_frame(searchlights(1.01), domain)
  view <- contrast_energy(
    plan_geometry(relation_value, frame,
      cross_partitions(relation_value, independence = "independent")),
    spectrum_contrast
  )
  metadata <- data.frame(
    measurement = as.character(view$index),
    scale = rep(1.01, length(view$index)),
    stringsAsFactors = FALSE
  )

  refusal <- catch_refusal(coherence_spectrum(view, using = metadata))
  expect_identical(refusal$capability, "conservative_frame")
  expect_match(refusal$message, "detection map", fixed = TRUE)
})

test_that("readouts without a coherence decomposition are refused with reasons", {
  fixture <- spectrum_plan(spectrum_smooth_signal(), id = "d5-refusals")

  # A query-only view holds one component, and a share is a ratio of two.
  single <- evaluate_geometry(fixture$plan,
    query = bilinear_query(tcrossprod(c(1, -1))), component = "total")
  refusal <- catch_refusal(coherence_spectrum(single))
  expect_identical(refusal$capability, "coherence_decomposition")
  expect_match(refusal$message, "one component", fixed = TRUE)

  # A comparative readout carries dissimilarities, not a split of a budget.
  distances <- rdm(fixture$plan)
  expect_identical(catch_refusal(coherence_spectrum(distances))$capability,
    "coherence_decomposition")

  # Packed geometry names no contrast to take the share of.
  geometry <- materialize_geometry(fixture$plan)
  expect_identical(catch_refusal(coherence_spectrum(geometry))$capability,
    "coherence_decomposition")

  # A plan needs the contrast; a view already has one and cannot take another.
  expect_error(coherence_spectrum(fixture$plan), "`weights` is required")
  expect_error(
    coherence_spectrum(contrast_energy(fixture$plan, spectrum_contrast),
      spectrum_contrast, using = fixture$frame$index),
    "applies only to a plan"
  )
  expect_error(coherence_spectrum(fixture$plan, spectrum_contrast,
    by = character()), "one or more columns")
})

test_that("a grouping column that is missing or NA is named, not dropped", {
  domain <- spectrum_domain(9L, id = "d5-na")
  relation_value <- relation(
    list(run1 = spectrum_pattern(spectrum_smooth_signal(9L)),
      run2 = spectrum_pattern(spectrum_smooth_signal(9L, phase = 0.4))),
    domain = domain
  )

  # A family whose second member is a region frame: its rows have no scale
  # and no center, and a spectrum grouped by scale would drop them silently
  # if the NA were treated as a group of its own or filtered away.
  family <- frame_family(
    narrow = compile_frame(searchlights(1.01, "conservative"), domain),
    parcels = compile_frame(
      regions(rep(c("a", "b", "c"), each = 3L), "conservative"), domain
    ),
    alpha = c(narrow = 0.5, parcels = 0.5)
  )
  plan <- plan_geometry(relation_value, family,
    cross_partitions(relation_value, independence = "independent"))

  expect_error(coherence_spectrum(plan, spectrum_contrast),
    "grouping `scale` is missing")

  # Grouping by `family` is the remedy the message names, and it works.
  by_family <- coherence_spectrum(plan, spectrum_contrast, by = "family")
  expect_identical(as.character(as.data.frame(by_family)$measurement),
    c("narrow", "parcels"))
  expect_equal(by_family$total, rep(sum(by_family$total) / 2, 2L),
    tolerance = 1e-12)

  # A frame that is not a family has no scale column at all, and the message
  # says which constructors produce one rather than only what is absent.
  plain <- plan_geometry(relation_value,
    compile_frame(searchlights(1.01, "conservative"), domain),
    cross_partitions(relation_value, independence = "independent"))
  bare <- try(coherence_spectrum(plain, spectrum_contrast), silent = TRUE)
  expect_s3_class(bare, "try-error")
  expect_match(conditionMessage(attr(bare, "condition")),
    "frame_family()", fixed = TRUE)
})

# The forbidden panel --------------------------------------------------------

test_that("plot refuses the per-scale energy profile and draws the split", {
  # Contract section 3.1 is normative: a panel of energy against scale is a
  # panel of the analyst's own weights and may not be presented as evidence
  # about spatial scale. `plot()` would draw exactly that as its profile, so
  # an explicit request for it is refused and the default is the panel whose
  # content is alpha-invariant.
  fixture <- spectrum_plan(spectrum_smooth_signal(), id = "d5-plot")
  spectrum <- coherence_spectrum(fixture$plan, spectrum_contrast)

  refusal <- catch_refusal(plot(spectrum, which = "profile"))
  expect_identical(refusal$capability, "scale_energy_panel")
  expect_match(refusal$message, "picture of `weights`", fixed = TRUE)
  expect_identical(catch_refusal(plot(spectrum, which = "both"))$capability,
    "scale_energy_panel")

  # The unqualified call draws the decomposition instead of refusing.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(spectrum))
  expect_silent(plot(spectrum, which = "decomposition"))

  # Nothing changes for a contrast view that is not a spectrum.
  view <- contrast_energy(fixture$plan, spectrum_contrast)
  expect_silent(plot(view, which = "profile"))
  expect_silent(plot(view))
})
