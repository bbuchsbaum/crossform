# `population_prevalence()` --- the descriptive group summaries of ticket E9.
#
# Four things are held here.
#
#   1. The arithmetic. Both fractions are checked against oracles written out
#      in this file with explicit loops over `fit$values` --- a different
#      evaluation order from the implementation's array reductions and, for
#      the alignment reference, a different algebra: the oracle averages the
#      complement of each participant explicitly where the implementation uses
#      `(S - v_i) / (n - 1)`. Where the fixture makes a fraction determinate
#      (a planted contrast every participant carries) the literal is asserted
#      too, because an oracle that shares a bug with its subject agrees with
#      it.
#   2. The leave-one-out reference, which is the one design decision in the
#      alignment measure that a reader cannot check from the number. A
#      leave-one-in mean puts `||v_i||^2 / n` inside every product, so it
#      counts an outlier as agreeing with a group it is pulling; the arithmetic
#      is exercised directly on a hand-made profile where the two answers
#      differ, which no fixture built through the executor would guarantee.
#   3. The latent-layer discipline. The print carries the same sentence
#      `latent_geometry()` prints, the record refuses to carry a field named
#      like an error bar, and constructing a prevalence leaves the result's own
#      `$uncertainty` byte-identical. That last one is the acceptance
#      criterion: the descriptive and the inferential layers are separate
#      objects and neither reaches into the other.
#   4. The refusals, each with the remedy it names.
#
# The fixture plants one effect --- face over house --- in **both** runs of
# every participant, so the crossvalidated geometry carries it and every
# participant is positive on that contrast; the rest of each participant's
# pattern is a phase-shifted cosine that reproduces across runs in none of
# them, so the `face-tool` contrast sits near the `0.5` the layer exists to
# warn about. A fixture in which every value were positive would let a
# prevalence of `1` pass for arithmetic.

pv_effects <- function() {
  effect_space(c("face", "house", "tool"), basis_id = "pop-prev:v1")
}

pv_subject <- function(id, features, phase, planted = 0.8, gain = 1) {
  domain <- abstract_domain(features,
    coordinates = cbind(x = seq_len(features) - 1),
    feature_ids = paste0("f", seq_len(features)), id = id)
  # The planted direction is in *effect* space and identical across runs, so
  # it survives the cross-partition product; the cosine is phase-shifted
  # between the two runs, so its cross-run product is signed and carries no
  # reproducible direction. Deterministic, because a prevalence checked
  # against an oracle must not depend on an RNG stream.
  shared <- outer(c(planted, -planted, 0), rep(1, features))
  values <- function(shift) {
    matrix(gain * cos(seq_len(3 * features) + shift), 3, features,
      dimnames = list(c("face", "house", "tool"), NULL)) + shared
  }
  relation <- relation(list(run1 = values(0), run2 = values(phase)),
    effects = pv_effects(), domain = domain)
  plan_geometry(relation, compile_frame(voxelwise(), domain),
    cross_partitions(relation))
}

# Group centres at 0, 4 and 9 with radius 1.5. A native node at x = 2, 6 or 7
# is two units from its nearest centre and lands entirely in the sink, so the
# sink carries real mass; group node 9 is out of reach of the two smaller
# participants under any radius, which is what gives the fixture a group node
# with partial subject coverage and --- under density semantics --- an `NA`
# rather than a `0`.
pv_carrier <- function(features, semantics = "budget", radius = 1.5, ...) {
  anatomical_transport(
    native_coords = cbind(seq_len(features) - 1),
    group_coords = cbind(c(0, 4, 9)),
    semantics = semantics, radius = radius, ...
  )
}

pv_sizes <- c(s01 = 6L, s02 = 8L, s03 = 10L, s04 = 12L)
pv_phases <- c(s01 = 0.2, s02 = 1.1, s03 = 2.4, s04 = 0.6)

pv_subjects <- function(sizes = pv_sizes, planted = 0.8, gain = 1) {
  stats::setNames(lapply(names(sizes), function(id)
    pv_subject(id, sizes[[id]], pv_phases[[id]], planted = planted,
      gain = gain)), names(sizes))
}

pv_plan <- function(sizes = pv_sizes, semantics = "budget", planted = 0.8,
                    gain = 1, ...) {
  plan_population(
    pv_subjects(sizes, planted = planted, gain = gain),
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      pv_carrier(sizes[[id]], semantics = semantics)),
    ...
  )
}

pv_bank <- function() {
  rbind(`face-house` = c(1, -1, 0), `face-tool` = c(1, 0, -1))
}

pv_fit <- function(..., queries = pv_bank()) {
  estimate_population(pv_plan(...), queries)
}

# The oracles. Explicit loops over the shipped `$values` array, sharing no
# code with `R/population-prevalence.R`.

pv_sign_oracle <- function(values, threshold = 0) {
  keep <- seq_len(dim(values)[[1L]])
  out <- matrix(NA_real_, length(keep), dim(values)[[2L]],
    dimnames = dimnames(values)[1:2])
  for (u in keep) {
    for (k in seq_len(dim(values)[[2L]])) {
      v <- values[u, k, ]
      v <- v[is.finite(v)]
      out[u, k] <- if (length(v)) sum(v > threshold) / length(v) else NA_real_
    }
  }
  out
}

# The reference for participant `i` is the mean of the *other* participants,
# formed here by averaging the complement explicitly rather than by subtracting
# `v_i` from a total.
pv_alignment_oracle <- function(values) {
  out <- rep(NA_real_, dim(values)[[1L]])
  for (u in seq_len(dim(values)[[1L]])) {
    profile <- matrix(values[u, , ], dim(values)[[2L]], dim(values)[[3L]])
    ok <- which(colSums(!is.finite(profile)) == 0L)
    if (length(ok) < 2L) next
    positive <- 0L
    for (i in ok) {
      others <- setdiff(ok, i)
      reference <- rowMeans(profile[, others, drop = FALSE])
      if (sum(profile[, i] * reference) > 0) positive <- positive + 1L
    }
    out[[u]] <- positive / length(ok)
  }
  out
}

# 1. The arithmetic ------------------------------------------------------------

test_that("sign prevalence is the fraction of participants above the threshold", {
  fit <- pv_fit()
  prevalence <- population_prevalence(fit)

  keep <- !fit$index$sink
  expected <- pv_sign_oracle(fit$values[keep, , , drop = FALSE])
  expect_equal(unname(prevalence$sign$fraction), unname(expected))

  # The counts and the denominators are published beside the fraction, and
  # the fraction has to be the one they produce.
  expect_identical(prevalence$sign$count,
    structure(as.integer(round(prevalence$sign$fraction *
      prevalence$sign$resolved)), dim = dim(prevalence$sign$count),
      dimnames = dimnames(prevalence$sign$count)))
  expect_true(all(prevalence$sign$resolved <= length(prevalence$subjects)))

  # The planted contrast is in both runs of every participant, so every
  # participant that reached a group node at all is positive there: the
  # numerator is the coverage count exactly. This literal is a property of the
  # fixture, not of the oracle, and it is also the reason a bare fraction is
  # not enough --- at the group node only half the participants reach, the
  # fraction reads `0.5` while every participant present carries the effect.
  expect_identical(unname(prevalence$sign$count[, "face-house"]),
    unname(prevalence$coverage$contributing[, "face-house"]))
  full <- prevalence$coverage$contributing[, "face-house"] ==
    length(prevalence$subjects)
  expect_gt(sum(full), 0L)
  expect_identical(unname(prevalence$sign$fraction[full, "face-house"]),
    rep(1, sum(full)))

  # And the contrast that reproduces in nobody sits at the null reference.
  expect_identical(prevalence$reference, 0.5)
  expect_true(all(prevalence$sign$fraction[, "face-tool"] <= 0.75))
})

test_that("the threshold moves the count and never the denominator", {
  fit <- pv_fit()
  at_zero <- population_prevalence(fit)
  keep <- !fit$index$sink
  values <- fit$values[keep, , , drop = FALSE]
  high <- max(values[is.finite(values)])

  raised <- population_prevalence(fit, threshold = high / 2)
  expect_equal(unname(raised$sign$fraction),
    unname(pv_sign_oracle(values, high / 2)))
  expect_identical(raised$sign$resolved, at_zero$sign$resolved)
  expect_true(all(raised$sign$count <= at_zero$sign$count))
  expect_true(any(raised$sign$count < at_zero$sign$count))

  # Above every value nothing is counted; below every value everything is.
  above <- population_prevalence(fit, threshold = high + 1)
  below <- population_prevalence(fit,
    threshold = min(values[is.finite(values)]) - 1)
  expect_true(all(above$sign$fraction == 0))
  expect_true(all(below$sign$fraction == 1))

  # The comparison is strictly greater and takes no tolerance: a threshold set
  # exactly at a participant's own value excludes that participant.
  one <- values[1L, 1L, 1L]
  exact <- population_prevalence(fit, threshold = one)
  expect_identical(exact$sign$count[[1L, 1L]],
    at_zero$sign$count[[1L, 1L]] - sum(values[1L, 1L, ] > 0 &
      values[1L, 1L, ] <= one))

  # The threshold names a different estimand, so it is in the identity.
  expect_false(identical(raised$scientific_plan_id,
    at_zero$scientific_plan_id))
})

test_that("alignment prevalence scores against a leave-one-out reference", {
  fit <- pv_fit()
  prevalence <- population_prevalence(fit)
  keep <- !fit$index$sink

  expect_equal(unname(prevalence$alignment$fraction),
    pv_alignment_oracle(fit$values[keep, , , drop = FALSE]))
  expect_identical(prevalence$alignment$reference,
    "leave_one_out_participant_mean")
  expect_identical(prevalence$receipt$prevalence$alignment_reference,
    "leave_one_out_participant_mean")
  expect_true(all(prevalence$alignment$resolved <= length(prevalence$subjects)))
})

test_that("a leave-one-in reference would count an outlier as agreeing", {
  # Four participants at one group node, profiles `3, 1, 1, -2` along the
  # first query coordinate. Against the plain group mean `0.75` the first
  # three are positive, so a leave-one-in count reports `0.75`. Against the
  # mean of the *others*, participant one's reference is exactly zero --- it
  # is the whole reason the group mean is positive at all --- so it drops out,
  # and the honest fraction is `0.5`. The gap is the self term
  # `||v_i||^2 / n` that leaving the participant in puts inside its own
  # product.
  profile <- rbind(c(3, 1, 1, -2), c(0, 0, 0, 0))
  values <- array(profile, c(1L, 2L, 4L))

  loo <- crossform:::.population_prevalence_alignment(values)
  expect_identical(loo$count, 2L)
  expect_equal(loo$fraction, 0.5)

  leave_one_in <- sum(colSums(profile * rowMeans(profile)) > 0) / ncol(profile)
  expect_equal(leave_one_in, 0.75)
  expect_gt(leave_one_in, loo$fraction)
  # The participant the two answers disagree about is the one that carries the
  # group mean.
  expect_true(colSums(profile * rowMeans(profile))[[1L]] > 0)

  # With one participant there is no cross-participant product at all, so the
  # fraction is withheld rather than reported as the `1` a self-product would
  # always produce.
  alone <- crossform:::.population_prevalence_alignment(
    array(c(3, 4), c(1L, 2L, 1L)))
  expect_identical(alone$resolved, 1L)
  expect_true(is.na(alone$fraction))
  expect_true(is.na(alone$count))
})

test_that("the readout Gram says whether alignment is a Frobenius product", {
  # A contrast bank is not orthonormal in packed coordinates: `svec(w w')` has
  # norm `||w||^2 = 2` for a two-condition contrast, so the Gram is 4 on the
  # diagonal and the departure from the identity is 3 --- an O(1) fact about
  # the bank, not a tolerance.
  contrasts <- population_prevalence(pv_fit())
  expect_false(contrasts$alignment$frobenius_equivalent)
  expect_equal(contrasts$alignment$readout_gram_deviation, 3)
  expect_identical(contrasts$alignment$inner_product,
    "query_readout_euclidean")

  # A bank of single-condition energies is orthonormal, and there the readout
  # inner product *is* the Frobenius inner product on the bank's span.
  energies <- population_prevalence(pv_fit(
    queries = rbind(face = c(1, 0, 0), house = c(0, 1, 0))))
  expect_true(energies$alignment$frobenius_equivalent)
  expect_equal(energies$alignment$readout_gram_deviation, 0)
  expect_true(energies$receipt$prevalence$alignment_frobenius_equivalent)
})

# 2. Denominators, the sink, and coverage --------------------------------------

test_that("the sink is excluded from every fraction and stays on the receipt", {
  fit <- pv_fit()
  prevalence <- population_prevalence(fit)

  # The precedent is `heterogeneity()`'s: the sink is unmapped territory in
  # budget units at no location, so it is not a row of a geometry summary. It
  # is excluded, the exclusion is recorded, and the budget stays readable.
  expect_false(any(prevalence$index$sink))
  expect_identical(nrow(prevalence$index), nrow(fit$index) - 1L)
  expect_false(crossform:::.transport_sink_label %in%
    rownames(prevalence$sign$fraction))
  expect_true(prevalence$receipt$prevalence$sink_excluded)
  expect_identical(prevalence$receipt$prevalence$sink_reason,
    "sink_is_budget_units_at_no_location")
  expect_identical(prevalence$receipt$sink_budget, fit$receipt$sink_budget)
  expect_true(any(fit$receipt$sink_budget != 0))
})

test_that("an unresolved cell leaves the denominator rather than dividing by N", {
  # Density semantics returns `NA` at a group node no native mass reached, and
  # group node 9 is out of reach of the two smaller participants.
  fit <- pv_fit(semantics = "density")
  values <- fit$values[!fit$index$sink, , , drop = FALSE]
  expect_true(anyNA(values))

  prevalence <- population_prevalence(fit)
  expect_equal(unname(prevalence$sign$fraction), unname(pv_sign_oracle(values)))
  # The row that lost participants lost them from the denominator, not from
  # the numerator only.
  unresolved <- apply(is.na(values), c(1L, 2L), sum)
  expect_identical(prevalence$sign$resolved,
    structure(as.integer(length(prevalence$subjects) - unresolved),
      dim = dim(unresolved), dimnames = dimnames(prevalence$sign$resolved)))
  expect_true(any(prevalence$sign$resolved < length(prevalence$subjects)))
  # A participant with no finite profile at a node is out of the alignment
  # count at that node too.
  expect_true(any(prevalence$alignment$resolved <
    length(prevalence$subjects)))
})

test_that("coverage counts contributing participants and marks a declared floor", {
  fit <- pv_fit()
  prevalence <- population_prevalence(fit)
  values <- fit$values[!fit$index$sink, , , drop = FALSE]

  contributing <- apply(is.finite(values) & values != 0, c(1L, 2L), sum)
  expect_identical(prevalence$coverage$contributing,
    structure(as.integer(contributing), dim = dim(contributing),
      dimnames = dimnames(prevalence$coverage$contributing)))
  expect_identical(prevalence$coverage$definition,
    "participants_with_a_finite_nonzero_transported_value")
  # Under budget semantics an unreached group node receives budget zero rather
  # than `NA`, so the coverage count is the number that separates a low
  # prevalence from a thin one.
  expect_true(any(prevalence$coverage$contributing <
    prevalence$sign$resolved))

  # No floor by default: `population-form-v1` section 14.3 records the
  # threshold as an open maintainer decision, so the number is reported and
  # nothing is marked.
  expect_true(is.na(prevalence$coverage$floor))
  expect_identical(prevalence$coverage$below_floor, character(0))

  marked <- population_prevalence(fit, coverage_floor = 4L)
  expect_identical(marked$coverage$floor, 4L)
  expect_identical(marked$coverage$below_floor,
    rownames(marked$sign$fraction)[apply(marked$coverage$contributing, 1L,
      min) < 4L])
  expect_gt(length(marked$coverage$below_floor), 0L)

  # The floor decides what is *marked* and not what was counted, so it stays
  # out of the identity, exactly as `heterogeneity()`'s node selection does.
  expect_identical(marked$scientific_plan_id, prevalence$scientific_plan_id)
  expect_identical(marked$sign$fraction, prevalence$sign$fraction)
})

test_that("query selection restricts both measures and enters the identity", {
  fit <- pv_fit()
  full <- population_prevalence(fit)
  one <- population_prevalence(fit, query = "face-house")

  expect_identical(one$query_labels, "face-house")
  expect_identical(dim(one$sign$fraction), c(nrow(one$index), 1L))
  expect_identical(one$sign$fraction[, "face-house"],
    full$sign$fraction[, "face-house"])
  # The alignment inner product is taken over exactly the selected
  # coordinates, so restricting the bank is a different estimand and not a
  # subset of the same one.
  expect_false(identical(one$scientific_plan_id, full$scientific_plan_id))
  expect_identical(one$alignment$fraction,
    stats::setNames(
      pv_alignment_oracle(fit$values[!fit$index$sink, "face-house", ,
        drop = FALSE]),
      rownames(one$sign$fraction)))

  expect_identical(population_prevalence(fit, query = 2L)$query_labels,
    "face-tool")

  # The receipt's own `$queries` names the bank the *run* was read through and
  # is carried unchanged; the subset this record counted over is beside the
  # threshold it belongs with, so neither field has to be read as the other.
  expect_identical(one$receipt$queries, fit$receipt$queries)
  expect_identical(one$receipt$prevalence$queries, "face-house")
  expect_identical(one$receipt$readout$labels, fit$receipt$readout$labels)
  expect_error(population_prevalence(fit, query = "face-nose"),
    class = "effect_input_error")
  expect_error(population_prevalence(fit, query = 7L),
    class = "effect_input_error")
})

# 3. The latent-layer discipline -----------------------------------------------

test_that("the record is a latent-layer record and says so in the D7 words", {
  prevalence <- population_prevalence(pv_fit())

  expect_identical(prevalence$layer, "latent_descriptive")
  expect_identical(prevalence$receipt$prevalence$layer, "latent_descriptive")
  # The sentence is the one `latent_geometry()` prints, shared as a constant
  # so the two descriptive layers cannot drift apart on it.
  expect_identical(prevalence$reading,
    "latent descriptive layer; not for inference")
  expect_identical(prevalence$reading, crossform:::.latent_reading_line)

  printed <- utils::capture.output(print(prevalence))
  expect_true(any(grepl("latent descriptive layer; not for inference", printed,
    fixed = TRUE)))
  expect_true(any(grepl("latent descriptive", printed, fixed = TRUE)))
  # The fractions are labelled descriptive on the printed page, not only in
  # the record.
  expect_true(any(grepl("descriptive only", printed, fixed = TRUE)))
  expect_true(any(grepl("near 0.5", printed, fixed = TRUE)))
})

test_that("the descriptive layer carries no inferential quantity and builds none", {
  plan <- pv_plan()
  fit <- estimate_population(plan, pv_bank())
  before <- fit$uncertainty
  prevalence <- population_prevalence(fit)

  # Nothing on the record may look like an error bar, at any depth.
  flatten <- function(x) {
    if (!is.list(x)) return(character())
    c(names(x), unlist(lapply(x, flatten), use.names = FALSE))
  }
  named <- flatten(prevalence[c("sign", "alignment", "coverage")])
  expect_identical(
    intersect(named, c("se", "t", "lower", "upper", "p_value", "p",
      "statistic", "level", "confidence")),
    character(0)
  )
  expect_identical(prevalence$receipt$prevalence$inference, "none_derivable")

  # And constructing one does not reach into the inferential layer: the
  # result's own uncertainty block is untouched, and the two verbs stay
  # separate objects with separate identities.
  expect_identical(fit$uncertainty, before)
  bars <- population_uncertainty(fit)
  expect_identical(fit$uncertainty, before)
  expect_false(identical(bars$scientific_plan_id,
    prevalence$scientific_plan_id))
  expect_false(inherits(prevalence, "effect_population_uncertainty"))

  # A validator that would accept a forged fraction is not a validator: the
  # fraction is re-derived from the count and denominator the record
  # publishes.
  forged <- prevalence
  forged$sign$fraction[[1L, 1L]] <- 0.5 * forged$sign$fraction[[1L, 1L]] + 0.1
  expect_error(crossform:::.validate_population_prevalence(forged),
    class = "effect_contract_error")

  smuggled <- prevalence
  smuggled$sign$se <- prevalence$sign$fraction
  expect_error(crossform:::.validate_population_prevalence(smuggled),
    class = "effect_contract_error")
})

test_that("the ledger name and the transport identity are on the printed page", {
  # `population-form-v1` section 8.1 makes this line normative: a transported
  # component summary states the native frame family it is a ledger of and the
  # transport that carried it, and a coherent ledger is never printed under
  # the bare name `coherent`.
  plan <- pv_plan()
  coherent <- population_prevalence(
    estimate_population(plan, pv_bank(), component = "coherent"))
  total <- population_prevalence(estimate_population(plan, pv_bank()))

  expect_identical(coherent$ledger, "native_coherent_ledger")
  expect_identical(total$ledger, "transported_total")
  # Two components of one plan are two results, so their prevalences are two
  # objects: the identity descends from the *result* and not from the plan,
  # which a plan id would have collided.
  expect_false(identical(coherent$scientific_plan_id,
    total$scientific_plan_id))
  expect_identical(coherent$receipt$prevalence$parent,
    estimate_population(plan, pv_bank(),
      component = "coherent")$scientific_plan_id)
  printed <- utils::capture.output(print(coherent))
  expect_true(any(grepl("native_coherent_ledger", printed, fixed = TRUE)))
  expect_true(any(grepl("frame:", printed, fixed = TRUE)))
  expect_true(any(grepl("transport:", printed, fixed = TRUE)))
  expect_true(any(grepl("not a group-node common mode", printed, fixed = TRUE)))
  expect_false(any(grepl("^  ledger: *coherent$", printed)))
})

test_that("the printed record is stable", {
  prevalence <- population_prevalence(pv_fit())
  # The estimand line shows a digest over BLAS-computed content, which the
  # package's own numerical contract does not promise bitwise across
  # platforms (`numerical_contract()$bitwise_across_platforms`); scrub it the
  # way every other snapshot in the suite scrubs digests.
  scrub <- function(lines) {
    sub("(population-sha256:)[0-9a-f]+", "\\1<digest>", lines)
  }
  expect_snapshot(print(prevalence), transform = scrub)
  expect_snapshot(format(prevalence), transform = scrub)
})

test_that("one measure at a time, because they are not one table", {
  prevalence <- population_prevalence(pv_fit())

  sign <- as.data.frame(prevalence)
  expect_identical(names(sign),
    c("node", "query", "prevalence", "participants", "resolved",
      "contributing"))
  expect_identical(nrow(sign),
    nrow(prevalence$index) * length(prevalence$query_labels))
  expect_equal(sign$prevalence, as.numeric(prevalence$sign$fraction))

  aligned <- as.data.frame(prevalence, measure = "alignment")
  expect_identical(names(aligned),
    c("node", "prevalence", "participants", "resolved"))
  expect_identical(nrow(aligned), nrow(prevalence$index))
  expect_error(as.data.frame(prevalence, measure = "both"))
})

# 4. The refusals --------------------------------------------------------------

test_that("a complete-form result has no participant axis to count over", {
  form <- materialize_population(pv_plan())
  refusal <- catch_refusal(population_prevalence(form))

  expect_identical(refusal$capability, "population_participant_values")
  expect_identical(refusal$namespace, "population_prevalence")
  expect_identical(refusal$reasons,
    "streamed_route_retains_no_participant_values")
  expect_match(refusal$remedies, "estimate_population", fixed = TRUE)
  expect_null(form$values)
})

test_that("a group-model term is refused, with the inferential layer as remedy", {
  fit <- pv_fit()
  refusal <- catch_refusal(population_prevalence(fit, term = "(Intercept)"))

  expect_identical(refusal$capability, "participant_term_decomposition")
  expect_identical(refusal$reasons, "participant_values_carry_no_term_axis")
  expect_true(any(grepl("population_uncertainty", refusal$remedies,
    fixed = TRUE)))
  # The refusal fires before anything else is checked, so a caller who passed
  # both `term` and an impossible query still learns about `term`.
  expect_identical(
    catch_refusal(population_prevalence(fit, query = 1L, term = 1L))$capability,
    "participant_term_decomposition")
})

test_that("a sign-flipped unit-budget participant is refused, not marked", {
  # Section 4.3's divisor is a signed native total. The fixture's phase-shifted
  # cosines put at least one participant's total below zero, and that
  # participant's whole ledger is then sign-flipped relative to the others.
  plan <- pv_plan(normalization = "unit_budget")
  fit <- estimate_population(plan, pv_bank())
  expect_true(any(fit$receipt$native_total < 0))

  refusal <- catch_refusal(population_prevalence(fit))
  expect_identical(refusal$capability, "comparable_ledger_orientation")
  expect_identical(refusal$namespace, "population_prevalence")
  expect_true("unit_budget_divisor_negative" %in% refusal$reasons)
  expect_true(any(grepl("normalization = \\\"none\\\"", refusal$remedies)))
  # The refusal names the participants, so which ones is readable without
  # re-running anything.
  flipped <- rownames(fit$receipt$native_total)[
    rowSums(fit$receipt$native_total < 0) > 0L]
  expect_identical(sort(sub("^participant:", "",
    grep("^participant:", refusal$reasons, value = TRUE))), sort(flipped))

  # A query on which every divisor is positive is not refused: the gate is per
  # query, because the divisor is.
  positive <- colnames(fit$receipt$native_total)[
    colSums(fit$receipt$native_total < 0) == 0L]
  skip_if(!length(positive), "no query has a uniformly positive divisor")
  expect_s3_class(population_prevalence(fit, query = positive[[1L]]),
    "effect_population_prevalence")
})

test_that("the form-space reading travels as a record, not as an error", {
  prevalence <- population_prevalence(pv_fit())
  refusal <- prevalence$receipt$form_prevalence_refusal

  # `heterogeneity()`'s `$receipt$latent_refusal` is the precedent: a reader of
  # any prevalence record learns that the Frobenius reading exists and why
  # this record does not carry it.
  expect_identical(refusal$capability, "participant_form_prevalence")
  expect_identical(refusal$namespace, "population_prevalence")
  expect_true("participant_transported_forms_not_retained" %in% refusal$reasons)
  expect_true("basis:query_bank" %in% refusal$reasons)
  expect_true(any(grepl("heterogeneity", refusal$remedies, fixed = TRUE)))
  expect_match(refusal$message, "Frobenius", fixed = TRUE)

  # It is a record and it did not stop the verb.
  expect_s3_class(prevalence, "effect_population_prevalence")
  expect_false(inherits(refusal, "condition"))
})

test_that("the verb refuses what is not a population result", {
  expect_error(population_prevalence(42), class = "effect_input_error")
  expect_error(population_prevalence(pv_plan()), class = "effect_input_error")
  expect_error(population_prevalence(pv_fit(), threshold = c(0, 1)),
    class = "effect_input_error")
  expect_error(population_prevalence(pv_fit(), threshold = NA_real_),
    class = "effect_input_error")
  expect_error(population_prevalence(pv_fit(), coverage_floor = 0L),
    class = "effect_input_error")
  expect_error(population_prevalence(pv_fit(), coverage_floor = 99L),
    class = "effect_input_error")
})
