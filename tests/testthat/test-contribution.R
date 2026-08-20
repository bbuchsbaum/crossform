# `contribution()` -- additive aggregation over a territory, ticket D4.
#
# `design/conservative-geometry-contract.md` is the source of every assertion
# here:
#
#   section 1.2  a conservative frame partitions one fixed global budget, so a
#                node's value is a share and shares add up. Section 1.1 says
#                the opposite of a local frame, whose overlapping nodes
#                double-count: summing a detection map estimates nothing.
#   section 2    claim 2 -- sum_x G_x = G_Omega, with the comparator taken
#                under `whole_brain("none")`. This is the identity a ledger
#                has to reproduce, and it is what the conservation tests below
#                measure rather than assume.
#   section 3.1  per scale, sum_{x in s} G_{s,x} = alpha_s G_Omega. Aggregating
#                a `frame_family()` by `family` is the readout of that law.
#   section 4    claim 4 -- coherent does NOT conserve, so an aggregated
#                coherent budget is frame-relative and a coherence fraction is
#                masked, never clamped.
#   section 11.4 gap G4 -- grouping is by row, budget-exact; overlap-splitting
#                is not offered.
#
# What is new here is the aggregation itself. The per-node laws it rests on are
# already covered by `test-conservative-geometry-contract.R` and
# `test-frame-family.R`, and this file does not restate them.

contribution_relation <- function(n = 12L, q = 3L, seed = 4471L, id = "d4") {
  domain <- abstract_domain(n, coordinates = cbind(seq_len(n) - 1, 0),
    feature_ids = paste0("vox", seq_len(n)), id = id)
  set.seed(seed)
  mk <- function() {
    value <- matrix(stats::rnorm(q * n), q, n)
    rownames(value) <- letters[seq_len(q)]
    value
  }
  relation_value <- relation(list(run1 = mk(), run2 = mk()), domain = domain)
  list(domain = domain, relation = relation_value,
    over = cross_partitions(relation_value),
    contrast = stats::setNames(c(1, -1, rep(0, q - 2L)), letters[seq_len(q)]))
}

contribution_view <- function(fixture, frame) {
  contrast_energy(
    plan_geometry(fixture$relation, frame, fixture$over), fixture$contrast
  )
}

# The comparator of claim 2: the same contrast under the unnormalized
# whole-domain operator. `whole_brain()` defaults to `"local"`, which divides
# by the feature count; `"none"` is the one that answers "what is the total".
contribution_global <- function(fixture) {
  contribution_view(
    fixture, compile_frame(whole_brain("none"), fixture$domain)
  )
}

# Conservation ---------------------------------------------------------------

test_that("a region ledger grouped by an external label adds back exactly", {
  fixture <- contribution_relation(id = "d4-regions")
  labels <- rep(paste0("r", 1:4), each = 3L)
  frame <- compile_frame(regions(labels, "conservative"), fixture$domain)
  view <- contribution_view(fixture, frame)

  # Four regions, two networks. The grouping arrives as a plain vector: a
  # network label is the user's own science, not something a frame carries.
  network <- c("dorsal", "dorsal", "ventral", "ventral")
  ledger <- contribution(view, by = network)

  expect_s3_class(ledger, "effect_contrast_view")
  expect_identical(length(ledger$total), 2L)
  expect_identical(as.character(ledger$index$measurement),
    c("dorsal", "ventral"))
  expect_identical(ledger$index$n_rows, c(2L, 2L))
  expect_identical(sum(ledger$index$n_rows), length(view$total))

  # Budget-exact: the groups partition the rows, and the rows partition the
  # whole-domain total (claim 2). Both halves are measured, not assumed.
  expect_equal(sum(ledger$total), sum(view$total), tolerance = 1e-12)
  global <- contribution_global(fixture)
  expect_equal(sum(ledger$total), global$total, tolerance = 1e-12)

  # `signed` is masked, not summed. It is the local weighted *mean* contrast,
  # so it is a density: summing it over four regions gives four times the
  # whole-domain value, which is the arithmetic section 1.1 forbids.
  expect_true(all(is.na(ledger$signed)))
  expect_identical(length(ledger$signed), 2L)
  expect_identical(ledger$metadata$aggregation$masked, "signed")
  expect_identical(ledger$metadata$aggregation$budget_exact, "total")
  expect_gt(abs(sum(view$signed) - global$signed), 1e-6)

  # Each group is the plain sum of its own rows -- no reweighting, no split.
  for (position in seq_along(network)) {
    group <- network[[position]]
    rows <- which(network == group)
    expect_equal(
      ledger$total[[match(group, as.character(ledger$index$measurement))]],
      sum(view$total[rows]), tolerance = 1e-12
    )
  }
})

test_that("a frame family aggregates by family and by scale under the alpha law", {
  fixture <- contribution_relation(id = "d4-family")
  alpha <- c(narrow = 0.4, wide = 0.6)
  family <- frame_family(
    narrow = compile_frame(searchlights(1.01, "conservative"), fixture$domain),
    wide = compile_frame(searchlights(2.01, "conservative"), fixture$domain),
    alpha = alpha
  )
  view <- contribution_view(fixture, family)
  global <- contribution_global(fixture)

  # `by` names a column of the family's own per-row metadata, joined through
  # `using`. A view carries measurement identifiers, not the table behind them.
  by_family <- contribution(view, by = "family", using = family$index)
  expect_identical(as.character(by_family$index$measurement),
    c("narrow", "wide"))
  expect_identical(by_family$index$n_rows,
    c(fixture$domain$n_features, fixture$domain$n_features))
  expect_equal(sum(by_family$total), global$total, tolerance = 1e-12)

  # Section 3.1: each block carries exactly its alpha share of the budget, so
  # a per-scale energy panel is a plot of the analyst's own alpha vector.
  expect_equal(unname(by_family$total),
    unname(alpha[c("narrow", "wide")] * global$total), tolerance = 1e-12)

  # The same rows grouped by the scale itself. A numeric key would be
  # flattened by the `measurement` column, so it also gets a typed column.
  by_scale <- contribution(view, by = "scale", using = family$index)
  expect_identical(by_scale$index$scale, c(1.01, 2.01))
  expect_identical(as.character(by_scale$index$measurement), c("1.01", "2.01"))
  expect_equal(by_scale$total, by_family$total, tolerance = 0)
  expect_equal(sum(by_scale$total), global$total, tolerance = 1e-12)

  # The family's `$index` may be given in any order: the join is by
  # `measurement`, not by position.
  shuffled <- family$index[rev(seq_len(nrow(family$index))), , drop = FALSE]
  expect_equal(contribution(view, by = "family", using = shuffled)$total,
    by_family$total, tolerance = 0)
})

test_that("a query-only effect view aggregates on the same rule", {
  fixture <- contribution_relation(id = "d4-query")
  frame <- compile_frame(searchlights(1.01, "conservative"), fixture$domain)
  plan <- plan_geometry(fixture$relation, frame, fixture$over)
  query <- bilinear_query(tcrossprod(unname(fixture$contrast)))
  view <- evaluate_geometry(plan, query = query)
  labels <- rep(c("anterior", "posterior"), each = 6L)

  ledger <- contribution(view, by = labels)
  expect_s3_class(ledger, "effect_view")
  expect_identical(dim(ledger$values), c(2L, 1L))
  expect_identical(as.character(ledger$index$measurement),
    c("anterior", "posterior"))
  expect_equal(sum(ledger$values), sum(view$values), tolerance = 1e-12)
  expect_equal(sum(ledger$values), drop(contribution_global(fixture)$total),
    tolerance = 1e-12)

  # The record stays canonical, so every reader of an `effect_view` still
  # works on the aggregate.
  expect_silent(crossform:::.validate_effect_view(ledger))
  frame_data <- as.data.frame(ledger)
  expect_identical(nrow(frame_data), 2L)
  expect_true(all(c("measurement", "n_rows") %in% names(frame_data)))

  # `total` is the only component with a global comparator, so it is the only
  # one whose aggregate is not labelled frame-relative.
  expect_false(ledger$metadata$aggregation$frame_relative)
  coherent <- contribution(
    evaluate_geometry(plan, query = query, component = "coherent"),
    by = labels
  )
  expect_true(coherent$metadata$aggregation$frame_relative)
  expect_identical(coherent$metadata$aggregation$frame_relative_components,
    "coherent")
})

# Frame-relative budgets and masked fractions --------------------------------

test_that("coherence fractions are recomputed from the aggregate and masked", {
  # Claim 4. A fraction of sums is not a sum of fractions, and a group whose
  # aggregated components are not a nonnegative partition reports NA rather
  # than a clamped number. This fixture is chosen so that both ways of failing
  # the mask occur: one group has a positive total but a negative coherent
  # part, another has a negative total outright.
  fixture <- contribution_relation(seed = 777L, id = "d4-mask")
  frame <- compile_frame(searchlights(1.01, "conservative"), fixture$domain)
  view <- contribution_view(fixture, frame)
  groups <- rep(c("g1", "g2", "g3"), each = 4L)
  ledger <- contribution(view, by = groups)

  expect_identical(ledger$coherence_fraction_valid, c(FALSE, TRUE, FALSE))
  expect_true(all(is.na(ledger$coherence_fraction[
    !ledger$coherence_fraction_valid])))

  # `g1` fails on the coherent part alone, so the mask is not a sign test on
  # the total; `g3` fails on the total.
  expect_gt(ledger$total[[1L]], 0)
  expect_lt(ledger$coherent[[1L]], 0)
  expect_lt(ledger$total[[3L]], 0)

  # Where it is defined, the fraction is the aggregated coherent over the
  # aggregated total -- and that is not the average of the nodes' fractions.
  expect_equal(ledger$coherence_fraction[[2L]],
    ledger$coherent[[2L]] / ledger$total[[2L]], tolerance = 1e-12)
  node_valid <- groups == "g2" & view$coherence_fraction_valid
  expect_true(any(node_valid))
  expect_gt(
    abs(ledger$coherence_fraction[[2L]] -
      mean(view$coherence_fraction[node_valid])),
    1e-6
  )

  # The exact partition survives aggregation: total = coherent + configuration
  # holds group by group, whatever the mask says about the fraction.
  expect_equal(ledger$total, ledger$coherent + ledger$configuration,
    tolerance = 1e-12)

  # Coherent budgets are this frame's own (claim 4), and the object says so
  # rather than leaving the reader to remember it.
  record <- ledger$metadata$aggregation
  expect_true(record$frame_relative)
  expect_identical(record$frame_relative_components,
    c("coherent", "configuration"))
  expect_identical(record$budget_exact, "total")
  printed <- paste(utils::capture.output(print(ledger)), collapse = " ")
  expect_match(printed, "frame_relative: TRUE")
  # A grouping passed as a bare symbol names itself, so the print says which
  # variable the territories came from.
  expect_identical(record$aggregated_by, "groups")
  expect_match(printed, "aggregated_by: `groups`")
  expect_match(printed, "masked: `signed`")
})

# Provenance -----------------------------------------------------------------

test_that("the aggregate records what it aggregated and by what", {
  fixture <- contribution_relation(id = "d4-provenance")
  labels <- rep(paste0("r", 1:4), each = 3L)
  frame <- compile_frame(regions(labels, "conservative"), fixture$domain)
  view <- contribution_view(fixture, frame)
  network <- c("dorsal", "dorsal", "ventral", "ventral")
  ledger <- contribution(view, by = network)
  record <- ledger$metadata$aggregation

  expect_identical(record$aggregated_by, "network")
  expect_identical(record$groups, 2L)
  expect_identical(record$measurements, 4L)
  expect_identical(record$frame_normalization, "conservative")
  expect_identical(record$group_keys, c("dorsal", "ventral"))
  # Gap G4's resolution, recorded rather than implied: no node's mass was
  # divided between territories.
  expect_false(record$overlap_split)
  expect_identical(ledger$metadata$aggregated_from,
    view$receipt$scientific_plan_id)

  # The aggregate is a different scientific object from the view it came from,
  # and the identity is a function of the grouping alone.
  expect_false(identical(ledger$receipt$scientific_plan_id,
    view$receipt$scientific_plan_id))
  expect_identical(contribution(view, by = network)$receipt$scientific_plan_id,
    ledger$receipt$scientific_plan_id)
  expect_false(identical(
    contribution(view, by = c("a", "a", "b", "b"))$receipt$scientific_plan_id,
    ledger$receipt$scientific_plan_id
  ))
  expect_match(ledger$receipt$task_partition_id, "projected_view")
  expect_identical(ledger$receipt$sources, view$receipt$sources)
})

# Refusals -------------------------------------------------------------------

test_that("a detection map is refused rather than summed", {
  # Section 1.1. `regions()` and `searchlights()` default to `"local"`, so the
  # refusal is what a first attempt actually meets.
  fixture <- contribution_relation(id = "d4-local")
  labels <- rep(paste0("r", 1:4), each = 3L)
  view <- contribution_view(fixture,
    compile_frame(regions(labels), fixture$domain))

  refusal <- catch_refusal(
    contribution(view, by = c("dorsal", "dorsal", "ventral", "ventral"))
  )
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "conservative_frame")
  expect_identical(refusal$namespace, "geometry_views")
  expect_identical(refusal$reasons, "local_frame_is_a_detection_map")
  expect_match(refusal$message, "double-count")
  expect_match(refusal$remedies, "conservative")

  # A frame that records no normalization at all is refused for its own
  # reason: the guard is not allowed to assume the safe case.
  bare <- view
  bare$metadata <- list()
  expect_identical(
    catch_refusal(contribution(bare, by = c("a", "a", "b", "b")))$reasons,
    "frame_normalization_not_recorded"
  )
})

test_that("views whose values are not contributions are refused by kind", {
  fixture <- contribution_relation(id = "d4-kinds")
  frame <- compile_frame(searchlights(1.01, "conservative"), fixture$domain)
  plan <- plan_geometry(fixture$relation, frame, fixture$over)
  geometry <- materialize_geometry(plan)
  labels <- rep(c("a", "b"), each = 6L)

  # Eigenvalues do not add across measurements at all.
  spectrum <- catch_refusal(
    contribution(geometry_spectrum(geometry), by = labels)
  )
  expect_identical(spectrum$capability, "additive_contribution")
  expect_identical(spectrum$reasons, "eigenvalues_are_not_additive")

  # A distance and a regression coefficient are comparisons, not shares of a
  # declared whole. The remedy is the additive route to the same numbers.
  for (view in list(rdm(plan), rsa(plan, models = list(
    m = matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), 3, 3)
  )))) {
    refusal <- catch_refusal(contribution(view, by = labels))
    expect_identical(refusal$capability, "additive_contribution")
    expect_identical(refusal$reasons, "comparative_readout_is_not_a_budget")
    expect_match(refusal$remedies, "evaluate_geometry")
  }

  # Packed geometry has no single per-measurement value, and a plan has none
  # at all.
  expect_identical(
    catch_refusal(contribution(geometry, by = labels))$reasons,
    "packed_geometry_is_not_a_readout"
  )
  expect_identical(
    catch_refusal(contribution(plan, by = labels))$reasons,
    "plan_holds_no_values"
  )
  expect_error(contribution(fixture$relation, by = labels),
    "effect_contrast_view", class = "effect_input_error")
})

test_that("a grouping that does not cover the rows exactly is refused", {
  fixture <- contribution_relation(id = "d4-grouping")
  labels <- rep(paste0("r", 1:4), each = 3L)
  frame <- compile_frame(regions(labels, "conservative"), fixture$domain)
  view <- contribution_view(fixture, frame)
  network <- c("dorsal", "dorsal", "ventral", "ventral")

  expect_error(contribution(view), "`by` is required",
    class = "effect_input_error")

  # A column name that names nothing says what is available instead of
  # failing on a NULL.
  family <- frame_family(
    narrow = compile_frame(searchlights(1.01, "conservative"), fixture$domain),
    alpha = c(narrow = 1)
  )
  family_view <- contribution_view(fixture, family)
  expect_error(
    contribution(family_view, by = "region", using = family$index),
    "names no column", class = "effect_input_error")
  expect_error(
    contribution(family_view, by = "region", using = family$index),
    "`family`", class = "effect_input_error")
  # With no table at all the message says that, rather than listing nothing.
  expect_error(contribution(view, by = "network"),
    "no `using` table was supplied", class = "effect_input_error")

  # A grouping vector must have one entry per measurement, and the message
  # names both counts.
  expect_error(contribution(view, by = c("dorsal", "ventral")),
    "view has 4 measurements and `by` has 2 entries",
    class = "effect_input_error")

  # A missing label would take part of the budget with it.
  expect_error(contribution(view, by = c("dorsal", NA, "ventral", "ventral")),
    "must belong to exactly one group", class = "effect_input_error")
  point_family <- frame_family(
    point = compile_frame(voxelwise("conservative"), fixture$domain),
    narrow = compile_frame(searchlights(1.01, "conservative"), fixture$domain),
    alpha = c(point = 0.5, narrow = 0.5)
  )
  expect_error(
    contribution(contribution_view(fixture, point_family), by = "scale",
      using = point_family$index),
    "The grouping `scale` is missing for 12 measurements",
    class = "effect_input_error")

  # `using` is metadata for a named column; alongside a grouping vector it
  # would be silently unused.
  expect_error(contribution(view, by = network, using = family$index),
    "only meaningful when `by` names", class = "effect_input_error")
})

test_that("a metadata table must cover every measurement exactly once", {
  fixture <- contribution_relation(id = "d4-using")
  family <- frame_family(
    narrow = compile_frame(searchlights(1.01, "conservative"), fixture$domain),
    wide = compile_frame(searchlights(2.01, "conservative"), fixture$domain),
    alpha = c(narrow = 0.5, wide = 0.5)
  )
  view <- contribution_view(fixture, family)

  expect_error(
    contribution(view, by = "family",
      using = family$index[-1L, , drop = FALSE]),
    "has no row for 1 measurement", class = "effect_input_error")

  doubled <- rbind(family$index, family$index[1L, , drop = FALSE])
  expect_error(contribution(view, by = "family", using = doubled),
    "repeats", class = "effect_input_error")

  # A table with no `measurement` column is positional, so its length is the
  # only thing that can be checked.
  positional <- data.frame(lobe = rep(c("front", "back"), length.out = 12L),
    stringsAsFactors = FALSE)
  expect_error(contribution(view, by = "lobe", using = positional),
    "carries no `measurement` column to join on",
    class = "effect_input_error")
  expect_error(contribution(view, by = "family", using = "family"),
    "must be a nonempty data frame", class = "effect_input_error")
})
