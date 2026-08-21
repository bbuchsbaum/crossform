# `heterogeneity()` --- the subject Gram of `population-form-v1` sections 5
# and 6.
#
# Four claims, each with its own test below.
#
#   1. The `N`-by-`N` Gram *is* the `h`-by-`h` geometry-space covariance. Its
#      nonzero spectrum matches to `1e-10`, its modes span the same subspace,
#      and the covariance the comparison is against is built here from the
#      package's public single-participant verbs --- materialize each
#      participant's complete packed geometry, transport it, residualize it ---
#      touching nothing this file is testing. That oracle is also the drift
#      test that binds the streamed route to the dense one.
#   2. The `sqrt(2)` packing is load-bearing. A naive `upper.tri` packing of
#      the same forms gives a different Gram and a different spectrum, and the
#      test states by how much: it is an `O(1)` disagreement, not a tolerance.
#   3. `sum_u K_u` is the global Gram. The per-node Grams are built densely
#      here and summed; section 5 measures the additivity at `3.55e-15` and
#      this pins it at `1e-12`.
#   4. The plug-in split is biased upward and the cross-fitted one is not. The
#      direction is the oracle's (`design/oracles/population-geometry-split.R`
#      section P6.b, `+62.7%`), and it is checked on a synthetic population
#      built like the oracle's: real between-participant structure plus
#      independent within-participant noise, so the plug-in Gram's trace has
#      somewhere inflated to come from.
#
# The fixture carries **four** runs, not three: a cross-fitted Gram needs two
# disjoint halves and each half needs a cross-generalized estimate of its own,
# so four partitions is the floor and a three-run fixture cannot test the
# contract's own estimator.

hx_effects <- function() {
  effect_space(c("face", "house", "tool"), basis_id = "pop-het:v1")
}

hx_subject <- function(id, features, gain = 1, runs = 4L,
                       effects = hx_effects(), amplitude = 0) {
  q <- length(effects$coordinates)
  domain <- abstract_domain(features,
    coordinates = cbind(x = seq_len(features) - 1),
    feature_ids = paste0("f", seq_len(features)), id = id)
  # Signed by construction, and one independent draw per run so that the
  # within-participant noise the cross-fit removes actually exists.
  #
  # `amplitude` plants shared between-participant structure: one fixed
  # direction in *effect* space, present in every run of a participant and
  # scaled by that participant's own amplitude. Effect space is the one axis
  # every participant shares -- native frames differ in size -- so this is the
  # only place a low-rank heterogeneity mode can be planted such that the
  # transported forms carry it at every group node.
  planted <- amplitude * outer(c(0.9, -0.6, 0.3)[seq_len(q)],
    rep(1, features))
  values <- function(k) {
    set.seed(sum(as.integer(charToRaw(id))) + features + 97L * k)
    matrix(gain * stats::rnorm(q * features) / (1 + k / 10), q,
      features, dimnames = list(effects$coordinates, NULL)) + planted
  }
  sources <- stats::setNames(lapply(seq_len(runs), values),
    paste0("run", seq_len(runs)))
  relation <- relation(sources, effects = effects, domain = domain)
  plan_geometry(relation, compile_frame(voxelwise(), domain),
    cross_partitions(relation))
}

hx_carrier <- function(features, semantics = "budget", radius = 1.5, ...) {
  anatomical_transport(
    native_coords = cbind(seq_len(features) - 1),
    group_coords = cbind(c(0, 4, 9)),
    semantics = semantics, radius = radius, ...
  )
}

# The general Gram oracle uses one common participant set at every ordinary
# node. Coverage-varying targets are tested at the population-result layer and
# must not become an NA-filled covariance geometry by accident.
hx_sizes <- c(s01 = 10L, s02 = 11L, s03 = 12L, s04 = 13L, s05 = 14L)
hx_gains <- c(s01 = 1, s02 = 1.6, s03 = 0.6, s04 = 1.2, s05 = 0.9)

hx_plan <- function(sizes = hx_sizes, runs = 4L, semantics = "budget", ...) {
  subjects <- stats::setNames(lapply(names(sizes), function(id)
    hx_subject(id, sizes[[id]], hx_gains[[id]], runs = runs)), names(sizes))
  plan_population(subjects,
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      hx_carrier(sizes[[id]], semantics = semantics)), ...)
}

# The dense oracle: every participant's complete packed geometry, transported,
# with the sink dropped and the group model residualized out. This is the
# `N`-by-`m`-by-`p` stack the streamed route refuses to build, built here on
# purpose so there is something to be equal to. It calls no population
# executor.
hx_dense_deviations <- function(plan, component = "total") {
  labels <- names(plan$subjects)
  full_nodes <- nrow(plan$group_index) + 1L
  width <- plan$subjects[[1L]]$packed_width
  stack <- matrix(NA_real_, length(labels), full_nodes * width,
    dimnames = list(labels, NULL))
  for (label in labels) {
    rows <- geometry_component(
      materialize_geometry(plan$subjects[[label]]), component
    )
    stack[label, ] <- as.numeric(transport_values(plan$transport[[label]], rows))
  }
  # The transport stacked node fastest, so the sink is every `full_nodes`th
  # row of each coordinate block.
  keep <- rep(c(rep(TRUE, full_nodes - 1L), FALSE), times = width)
  kept <- stack[, keep, drop = FALSE]
  # A column carrying a non-finite response is a measurement: density
  # semantics returns `NA` at a group node reached by no native mass. It is
  # solved around rather than through, and comes back absent.
  finite <- apply(is.finite(kept), 2L, all)
  deviations <- matrix(NA_real_, nrow(kept), ncol(kept),
    dimnames = dimnames(kept))
  deviations[, finite] <- qr.resid(plan$model$qr, kept[, finite, drop = FALSE])
  deviations
}


# 1. The acceptance oracle -----------------------------------------------------

test_that("the streamed subject Gram is the geometry-space covariance", {
  plan <- hx_plan()
  split <- heterogeneity(plan, estimator = "plug_in")
  deviations <- hx_dense_deviations(plan)
  residual_df <- length(plan$subjects) - plan$model$rank

  expect_identical(split$estimator, "plug_in")
  expect_identical(split$space, "packed_form")
  expect_identical(split$residual_df, as.integer(residual_df))

  # The Gram itself, streamed against dense.
  gram <- tcrossprod(deviations) / residual_df
  expect_lt(max(abs(split$gram - gram)), 1e-12)

  # The `h`-by-`h` covariance, whose nonzero spectrum the `N`-by-`N` Gram
  # reproduces. `h` is the packed width times the retained group nodes: 18
  # here against a 5-by-5 Gram.
  covariance <- crossprod(deviations) / residual_df
  expect_gt(ncol(covariance), nrow(gram))
  dense <- sort(eigen(covariance, symmetric = TRUE,
    only.values = TRUE)$values, decreasing = TRUE)
  streamed <- sort(unname(split$spectrum), decreasing = TRUE)
  expect_lt(max(abs(streamed[seq_len(residual_df)] -
    dense[seq_len(residual_df)])), 1e-10)
  # Everything past the residual degrees of freedom is a numerical zero in
  # both: rank is bounded by `N - rank(X)`, which is section 5's `N - 1` for
  # the intercept-only design.
  expect_lt(max(abs(dense[-seq_len(residual_df)])), 1e-10)
  expect_lt(max(abs(streamed[-seq_len(residual_df)])), 1e-10)

  # Section 5's mode recovery: `v = C^T u` is an eigenvector of the covariance
  # with the same eigenvalue, so the two decompositions name the same
  # subspace. Measured as the principal angle between the leading planes.
  modes <- crossprod(deviations, split$loadings[, 1:2, drop = FALSE])
  modes <- apply(modes, 2L, function(column) column / sqrt(sum(column^2)))
  reference <- eigen(covariance, symmetric = TRUE)$vectors[, 1:2, drop = FALSE]
  angles <- svd(crossprod(modes, reference))$d
  expect_lt(max(abs(angles - 1)), 1e-8)
})

test_that("a mode form unpacks to the symmetric geometry it names", {
  plan <- hx_plan()
  split <- heterogeneity(plan, estimator = "plug_in", nodes = c(1L, 2L),
    modes = 2L)
  deviations <- hx_dense_deviations(plan)
  loadings <- split$loadings[, 1:2, drop = FALSE]
  weighted <- crossprod(deviations, loadings)
  scale <- sqrt(colSums(weighted^2))

  expect_identical(dim(split$mode_forms), c(2L, 2L, 6L))
  expect_identical(names(dimnames(split$mode_forms)),
    c("node", "mode", "coordinate"))
  expect_identical(split$nodes, c("group1", "group2"))

  # The tile stacked node fastest, so node `u` of coordinate `c` is row
  # `(c - 1) * m + u` of the dense contraction.
  m <- nrow(split$index)
  for (mode in 1:2) {
    expected <- matrix(weighted[, mode], m, 6L)[1:2, ] / scale[[mode]]
    expect_lt(max(abs(split$mode_forms[, mode, ] - expected)), 1e-11)
  }

  # Every slice is one packed form and unpacks to a symmetric matrix.
  form <- crossform:::.unsvec_symmetric(split$mode_forms["group1", "mode1", ],
    3L)
  expect_equal(form, t(form))
  # The mode is one unit-norm direction over the whole retained block, so a
  # single node's slice carries a share of that norm and every node's slices
  # together carry all of it. That is what makes the reported slices
  # comparable across nodes.
  expect_lt(sum(split$mode_forms["group1", "mode1", ]^2), 1)
  everywhere <- heterogeneity(plan, estimator = "plug_in",
    nodes = seq_len(m), modes = 1L)
  expect_equal(sum(everywhere$mode_forms[, "mode1", ]^2), 1)
})


# 2. Frobenius consistency -----------------------------------------------------

test_that("the sqrt(2) packing is load-bearing for the Gram", {
  plan <- hx_plan()
  split <- heterogeneity(plan, estimator = "plug_in")
  deviations <- hx_dense_deviations(plan)
  residual_df <- length(plan$subjects) - plan$model$rank

  # Undo the codec's off-diagonal weight to get the naive `upper.tri` packing
  # of the same forms, then take the same Gram of it. The scale vector is the
  # packed coordinate table's own -- one entry per coordinate, repeated over
  # the node axis the stack varies fastest along -- so this is exactly the
  # packing the contract measures failing.
  naive <- sweep(deviations, 2L,
    rep(split$coordinates$scale, each = nrow(split$index)), "/")
  naive_gram <- tcrossprod(naive) / residual_df

  # An `O(1)` disagreement, not a tolerance: the test would pass with a
  # tolerance of 0.1 and fail with any packing correction, which is the point.
  relative <- max(abs(naive_gram - split$gram)) / max(abs(split$gram))
  expect_gt(relative, 0.05)
  naive_spectrum <- sort(eigen(naive_gram, symmetric = TRUE,
    only.values = TRUE)$values, decreasing = TRUE)
  expect_gt(abs(naive_spectrum[[1L]] - max(split$spectrum)) /
    max(split$spectrum), 0.05)

  # The claim the Gram actually rests on, asserted against the *forward* codec
  # rather than against the scale table: `K[i, j]` is a sum of Frobenius inner
  # products of `q`-by-`q` forms, and it equals the packed inner product only
  # because `.svec_symmetric()` carries the root two. Strip the root two from
  # the codec and this breaks, because `unsvec` would still reconstruct the
  # right forms while the packed inner product stopped being their Frobenius
  # one. That is the failure mode the contract says a future storage change
  # must not be able to cause silently.
  m <- nrow(split$index)
  frobenius <- function(i, j) {
    sum(vapply(seq_len(m), function(node) {
      columns <- (seq_len(nrow(split$coordinates)) - 1L) * m + node
      sum(crossform:::.unsvec_symmetric(deviations[i, columns], 3L) *
        crossform:::.unsvec_symmetric(deviations[j, columns], 3L))
    }, numeric(1)))
  }
  for (pair in list(c(1L, 1L), c(1L, 2L), c(3L, 5L))) {
    expect_equal(split$gram[pair[[1L]], pair[[2L]]] * residual_df,
      frobenius(pair[[1L]], pair[[2L]]))
  }

  # And the forward codec itself, on forms this file draws rather than on
  # anything the population layer produced: it is Frobenius consistent and the
  # naive packing is not.
  set.seed(20260820)
  worst_codec <- 0
  worst_naive <- 0
  for (draw in seq_len(50L)) {
    a <- matrix(stats::rnorm(9L), 3L)
    a <- (a + t(a)) / 2
    b <- matrix(stats::rnorm(9L), 3L)
    b <- (b + t(b)) / 2
    target <- sum(a * b)
    worst_codec <- max(worst_codec, abs(sum(
      crossform:::.svec_symmetric(a) * crossform:::.svec_symmetric(b)
    ) - target) / max(1, abs(target)))
    worst_naive <- max(worst_naive, abs(sum(
      a[lower.tri(a, diag = TRUE)] * b[lower.tri(b, diag = TRUE)]
    ) - target) / max(1, abs(target)))
  }
  expect_lt(worst_codec, 1e-12)
  expect_gt(worst_naive, 0.1)

  expect_identical(split$receipt$gram$packing,
    "symmetric_packed_sqrt2_off_diagonal")
})


# 3. Node additivity -----------------------------------------------------------

test_that("the global Gram is the sum of the per-node Grams", {
  plan <- hx_plan()
  split <- heterogeneity(plan, estimator = "plug_in")
  deviations <- hx_dense_deviations(plan)
  residual_df <- length(plan$subjects) - plan$model$rank
  m <- nrow(split$index)
  width <- nrow(split$coordinates)

  accumulated <- matrix(0, nrow(split$gram), ncol(split$gram))
  for (node in seq_len(m)) {
    columns <- (seq_len(width) - 1L) * m + node
    accumulated <- accumulated +
      tcrossprod(deviations[, columns, drop = FALSE]) / residual_df
  }
  expect_lt(max(abs(accumulated - split$gram)), 1e-12)
  expect_true(split$receipt$gram$node_additive)

  # The tile is a memory bound and nothing else: one coordinate at a time is
  # the same Gram and the same estimand identity.
  tiled <- heterogeneity(plan, estimator = "plug_in", coordinate_tile = 1L)
  expect_lt(max(abs(tiled$gram - split$gram)), 1e-12)
  expect_identical(tiled$scientific_plan_id, split$scientific_plan_id)
  expect_identical(tiled$receipt$streaming$passes_per_subject, 6L)
})

test_that("a node nobody measured is held out, not zeroed", {
  # Density semantics returns `NA` at a group node reached by no native mass,
  # and it is `NA` rather than `0` because zero is a measurement and absence is
  # not. The third group centre sits at x = 9 with radius 1.5, so the three
  # smaller participants reach it with nothing at all.
  partial_sizes <- c(s01 = 6L, s02 = 8L, s03 = 10L, s04 = 12L, s05 = 7L)
  plan <- hx_plan(sizes = partial_sizes, semantics = "density")
  split <- heterogeneity(plan, estimator = "plug_in", nodes = 1:3, modes = 2L)
  width <- nrow(split$coordinates)

  # One unreached node times the packed width, held out of the accumulation.
  expect_identical(split$receipt$unresolved_columns, as.integer(width))
  expect_true(all(is.finite(split$gram)))

  # The Gram is the Gram of the resolved cells only, and equals a dense
  # accumulation over the first two nodes.
  deviations <- hx_dense_deviations(plan)
  resolved <- apply(is.finite(deviations), 2L, all)
  expect_identical(sum(!resolved), as.integer(width))
  expect_lt(max(abs(split$gram -
    tcrossprod(deviations[, resolved, drop = FALSE]) / split$residual_df)),
    1e-12)

  # And the unmeasured node's mode geometry is absent, not a zero form: a zero
  # there would read as "this mode has no geometry here" rather than as "this
  # node was never measured".
  expect_true(all(is.na(split$mode_forms["group3", , ])))
  expect_false(anyNA(split$mode_forms[c("group1", "group2"), , ]))

  lines <- capture.output(print(split))
  expect_true(any(grepl("held out of the Gram", lines)))
})

test_that("the sink is excluded from the Gram and said to be", {
  plan <- hx_plan()
  split <- heterogeneity(plan, estimator = "plug_in")
  expect_false(any(split$index$sink))
  expect_identical(nrow(split$index), nrow(plan$group_index))
  expect_true(split$receipt$gram$sink_excluded)
  # The sink budget is still readable where a failing transport is visible.
  expect_identical(dim(split$receipt$sink_budget),
    c(length(plan$subjects), nrow(split$coordinates)))
  expect_true(any(abs(split$receipt$sink_budget) > 0))
})


# 4. The cross-fitted estimator ------------------------------------------------

test_that("the cross-fitted Gram reads across two disjoint partition halves", {
  plan <- hx_plan()
  split <- heterogeneity(plan)

  expect_identical(split$estimator, "cross_fit")
  record <- split$receipt$cross_fit
  expect_identical(record$split, "interleaved_partitions")
  halves <- record$partitions[["s01"]]
  expect_identical(halves$a, c("run1", "run3"))
  expect_identical(halves$b, c("run2", "run4"))
  expect_length(intersect(halves$a, halves$b), 0L)
  expect_true(all(vapply(record$half_plan_ids, function(id)
    grepl("^population-sha256:", id), logical(1))))
  expect_false(identical(record$half_plan_ids[["a"]],
    record$half_plan_ids[["b"]]))

  # It is the same matrix section 6.4 writes down, built here from the two
  # halves read separately through the public planner.
  halves_gram <- lapply(c("a", "b"), function(side) {
    subjects <- stats::setNames(lapply(names(plan$subjects), function(label) {
      subject <- plan$subjects[[label]]
      plan_geometry(subject$task$left_relation, subject$frame,
        cross_partitions(record$partitions[[label]][[side]]),
        compute = subject$compute)
    }), names(plan$subjects))
    hx_dense_deviations(plan_population(subjects, plan$transport))
  })
  cross <- tcrossprod(halves_gram[[1L]], halves_gram[[2L]])
  expected <- (cross + t(cross)) / 2 / split$residual_df
  expect_lt(max(abs(split$gram - expected)), 1e-12)
})

test_that("the plug-in Gram is larger than the cross-fitted one", {
  # Built like the oracle's fixture: planted between-participant structure
  # that both halves see, plus independent per-run noise that only the plug-in
  # estimator books as heterogeneity.
  effects <- hx_effects()
  sizes <- hx_sizes
  amplitudes <- c(s01 = -1.2, s02 = 0.9, s03 = 1.5, s04 = -0.4, s05 = 0.2)
  subjects <- stats::setNames(lapply(names(sizes), function(id)
    hx_subject(id, sizes[[id]], gain = 1, runs = 4L, effects = effects,
      amplitude = amplitudes[[id]])), names(sizes))
  plan <- plan_population(subjects,
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      hx_carrier(sizes[[id]])))

  plug <- heterogeneity(plan, estimator = "plug_in")
  cross <- heterogeneity(plan)

  # The direction the oracle records: within-participant sampling noise lands
  # entirely in the plug-in heterogeneity trace and nowhere in the
  # cross-fitted one.
  expect_gt(sum(diag(plug$gram)), sum(diag(cross$gram)))
  expect_gt(sum(plug$spectrum), sum(cross$spectrum))
  # And it is not a rounding difference.
  expect_gt((sum(diag(plug$gram)) - sum(diag(cross$gram))) /
    abs(sum(diag(cross$gram))), 0.1)

  # The cross-fitted Gram is indefinite by construction; the plug-in one is a
  # Gram of real vectors and is not.
  expect_gt(sum(cross$spectrum < -1e-8 * max(abs(cross$spectrum))), 0L)
  expect_gte(min(plug$spectrum), -1e-8 * max(plug$spectrum))

  # Section 6.3: the loading directions survive the inflation, so the two
  # estimators agree on the leading mode up to sign even where their
  # eigenvalues do not.
  expect_gt(abs(sum(plug$loadings[, "mode1"] * cross$loadings[, "mode1"])),
    0.7)
})

test_that("a declared partition split is honoured and validated", {
  plan <- hx_plan()
  declared <- heterogeneity(plan,
    partitions = list(c("run1", "run2"), c("run3", "run4")))
  expect_identical(declared$receipt$cross_fit$split, "declared")
  expect_identical(declared$receipt$cross_fit$partitions[["s03"]]$a,
    c("run1", "run2"))
  expect_false(identical(declared$gram, heterogeneity(plan)$gram))
  # The split is part of the estimand: two splits are two records.
  expect_false(identical(declared$scientific_plan_id,
    heterogeneity(plan)$scientific_plan_id))

  shared <- catch_refusal(heterogeneity(plan,
    partitions = list(c("run1", "run2"), c("run2", "run3"))))
  expect_identical(shared$capability, "cross_fitted_subject_gram")
  expect_true(any(grepl("partition_in_both_halves", shared$reasons)))

  absent <- catch_refusal(heterogeneity(plan,
    partitions = list(c("run1", "run2"), c("run3", "run9"))))
  expect_true(any(grepl("partition_absent_for_subject", absent$reasons)))

  short <- catch_refusal(heterogeneity(plan,
    partitions = list("run1", c("run2", "run3"))))
  expect_true(any(grepl("half_carries_fewer_than_two", short$reasons)))
})

test_that("a transport that saw the split's partitions is reported", {
  # A functional transport learned from the responses has seen some runs. If
  # those overlap a cross-fit half then that half's inner product is not fully
  # held out. crossform cannot repair it -- the transport is a declared input
  # -- so it is measured, recorded and printed rather than ignored.
  plan <- hx_plan()
  # The same operator, redeclared as functional. `anatomical_transport()`
  # refuses a declared method, and `location_transport()` appends the sink
  # itself, so the group columns are handed over without it.
  circular <- lapply(stats::setNames(names(hx_sizes), names(hx_sizes)),
    function(id) {
      carrier <- plan$transport[[id]]
      external_transport(
        carrier$matrix[, seq_len(ncol(carrier$matrix) - 1L), drop = FALSE],
        semantics = "budget",
        native_index = carrier$native_index,
        group_index = carrier$group_index,
        provenance = list(method = "functional",
          details = "response clustering", cross_fit = c("run1", "run3")))
    })
  circular_plan <- plan_population(plan$subjects, circular)
  split <- heterogeneity(circular_plan)

  overlap <- split$receipt$cross_fit$transport_partition_overlap
  expect_length(overlap, 2L * length(hx_sizes))
  expect_true(all(grepl("^s0[1-5]:run[13]$", overlap)))
  lines <- capture.output(print(split))
  expect_true(any(grepl("not fully held out", lines)))

  # An anatomical transport never saw the responses and owes no such record.
  expect_length(
    heterogeneity(plan)$receipt$cross_fit$transport_partition_overlap, 0L)
})

test_that("a whitened metric schedule cannot be split into halves", {
  # Whitening freezes one set of effect coordinates from the whole partition
  # set, so a half of it is not the same measurement and the two halves would
  # not be comparable. The refusal is named rather than silently re-planning
  # onto a different schedule.
  effects <- hx_effects()
  sizes <- hx_sizes[1:2]
  subjects <- stats::setNames(lapply(names(sizes), function(id) {
    features <- sizes[[id]]
    domain <- abstract_domain(features,
      coordinates = cbind(x = seq_len(features) - 1),
      feature_ids = paste0("f", seq_len(features)), id = id)
    plain <- hx_subject(id, features, effects = effects)
    plan_geometry(plain$task$left_relation, plain$frame,
      plain$pairing,
      metric = neural_metric(diag(features) + 0.1, domain),
      composition = "whitened")
  }), names(sizes))
  plan <- plan_population(subjects,
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      hx_carrier(sizes[[id]])))
  expect_identical(plan$subjects[[1L]]$metric_schedule$kind,
    "whitened_metric_before_frame")

  refusal <- catch_refusal(heterogeneity(plan))
  expect_identical(refusal$capability, "cross_fitted_subject_gram")
  expect_true(any(grepl("metric_schedule_not_repartitionable", refusal$reasons)))
  expect_true(any(grepl("whitened_metric_before_frame", refusal$reasons)))
  # The plug-in route reads one pass and needs no re-plan, so it still runs.
  expect_s3_class(heterogeneity(plan, estimator = "plug_in"),
    "effect_population_heterogeneity")
})

test_that("three partitions cannot be cross-fitted, and say so", {
  plan <- hx_plan(runs = 3L)
  refusal <- catch_refusal(heterogeneity(plan))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "cross_fitted_subject_gram")
  expect_true(any(grepl("half_carries_fewer_than_two", refusal$reasons)))
  expect_true(any(grepl("plug_in", refusal$remedies)))
  # And the plug-in route still runs on the same plan.
  expect_s3_class(heterogeneity(plan, estimator = "plug_in"),
    "effect_population_heterogeneity")
})


# The latent layer and its refusals ---------------------------------------------

test_that("nonnegative functionals live only on the cross-fitted projection", {
  plan <- hx_plan()
  cross <- heterogeneity(plan)
  plug <- heterogeneity(plan, estimator = "plug_in")

  expect_null(plug$latent)
  refusal <- plug$receipt$latent_refusal
  expect_identical(refusal$capability, "plug_in_spectrum_functionals")
  expect_true(any(grepl("within_subject_noise", refusal$reasons)))
  expect_match(refusal$message, "62.7%", fixed = TRUE)
  expect_identical(plug$receipt$plug_in_reference$measured_inflation, 0.627)

  expect_identical(cross$latent$method, "psd_projection")
  expect_identical(cross$latent$operator, "eigenvalue_truncation_at_zero")
  expect_true(all(cross$latent$spectrum >= 0))
  expect_equal(cross$latent$spectrum, pmax(unname(cross$spectrum), 0))
  # The moved mass is the magnitude of what was clipped, recorded rather than
  # absorbed (section 6.5, gap G7's discipline one level up).
  expect_equal(cross$latent$moved_mass, sum(pmax(-cross$spectrum, 0)))
  expect_equal(cross$latent$moved_share,
    cross$latent$moved_mass / sum(abs(cross$spectrum)))
  expect_gt(cross$latent$n_eff, 0)
  expect_lte(cross$latent$n_eff, length(cross$subjects))
  expect_equal(cross$latent$cumulative[[length(cross$subjects)]], 1)
  expect_identical(cross$latent$negative_modes,
    crossform:::.heterogeneity_negative_modes(cross$spectrum))

  # `nearest_psd` is a declared member of the closed set and is not
  # implemented; it refuses here for the same reason it refuses in
  # `latent_geometry()`, rather than falling back to a different operator.
  reserved <- catch_refusal(heterogeneity(plan, projection = "nearest_psd"))
  expect_identical(reserved$capability, "nearest_psd_projection")
})


# The query-space route ---------------------------------------------------------

test_that("an estimated result gives the Gram in its own query space", {
  plan <- hx_plan()
  bank <- rbind(`face-house` = c(1, -1, 0), `face-tool` = c(1, 0, -1))
  fit <- estimate_population(plan, bank)
  split <- heterogeneity(fit)

  expect_identical(split$estimator, "plug_in")
  expect_identical(split$space, "query_bank")
  expect_identical(split$receipt$gram$packing, "declared_query_bank")
  expect_identical(as.character(split$coordinates$coordinate), rownames(bank))
  expect_null(split$effects)
  expect_null(split$latent)

  # It is the Gram of the residuals the result already carries, with the sink
  # dropped.
  keep <- !fit$index$sink
  block <- t(matrix(fit$residuals[keep, , , drop = FALSE],
    sum(keep) * nrow(bank), length(plan$subjects)))
  expect_lt(max(abs(split$gram - tcrossprod(block) / split$residual_df)), 1e-12)

  # Cross-fitting is a property of the plan, not of a finished result.
  refusal <- catch_refusal(heterogeneity(fit, estimator = "cross_fit"))
  expect_identical(refusal$capability, "query_space_cross_fit")
  expect_true(any(grepl("heterogeneity", refusal$remedies)))

  # A complete-form result carries no participant values at all.
  form <- catch_refusal(heterogeneity(materialize_population(plan)))
  expect_identical(form$capability, "heterogeneity_readout")
})

test_that("query-space mode values agree with the dense contraction", {
  plan <- hx_plan()
  fit <- estimate_population(plan, rbind(`face-house` = c(1, -1, 0),
    `face-tool` = c(1, 0, -1)))
  split <- heterogeneity(fit, nodes = "group2", modes = 2L)
  keep <- !fit$index$sink
  block <- t(matrix(fit$residuals[keep, , , drop = FALSE], sum(keep) * 2L,
    length(plan$subjects)))
  weighted <- crossprod(block, split$loadings[, 1:2, drop = FALSE])
  scale <- sqrt(colSums(weighted^2))
  for (mode in 1:2) {
    expected <- matrix(weighted[, mode], sum(keep), 2L)[2L, ] / scale[[mode]]
    expect_lt(max(abs(split$mode_forms["group2", mode, ] - expected)), 1e-11)
  }
})


# The record, its identity and its refusals -------------------------------------

test_that("the record is sealed, re-derivable and named by its estimator", {
  plan <- hx_plan()
  split <- heterogeneity(plan)

  expect_s3_class(split, "effect_population_heterogeneity")
  expect_identical(names(split), crossform:::.population_heterogeneity_fields)
  expect_identical(split$ledger, "transported_total")
  expect_identical(split$component, "total")
  expect_identical(split$semantics, "budget")

  # The spectrum is a function of the Gram, and the validator re-derives it
  # rather than trusting it.
  forged <- split
  forged$spectrum[[1L]] <- forged$spectrum[[1L]] + 1
  expect_error(crossform:::.validate_population_heterogeneity(forged),
    class = "effect_contract_error")

  # The estimator is part of the identity: two estimators are two estimands.
  plug <- heterogeneity(plan, estimator = "plug_in")
  expect_false(identical(plug$scientific_plan_id, split$scientific_plan_id))
  # The node selection is not: it decides what is displayed.
  shown <- heterogeneity(plan, nodes = "group1")
  expect_identical(shown$scientific_plan_id, split$scientific_plan_id)

  # A latent layer on a plug-in spectrum is a contract error, not an option.
  smuggled <- plug
  smuggled$latent <- split$latent
  expect_error(crossform:::.validate_population_heterogeneity(smuggled),
    class = "effect_contract_error")
})

test_that("the components carry their section 8.1 ledger names", {
  plan <- hx_plan()
  total <- heterogeneity(plan, estimator = "plug_in")
  for (component in c("coherent", "configuration")) {
    split <- heterogeneity(plan, estimator = "plug_in", component = component)
    expect_identical(split$component, component)
    expect_identical(split$ledger,
      unname(crossform:::.population_ledger_names[[component]]))
    expect_false(identical(split$scientific_plan_id, total$scientific_plan_id))
    # The component is carried and not merely named: each Gram is the Gram of
    # that component's own transported deviations.
    expect_lt(max(abs(split$gram -
      tcrossprod(hx_dense_deviations(plan, component)) / split$residual_df)),
      1e-12)
  }

  # On a voxelwise frame the ledger identity is degenerate and that is what a
  # reader should see: a singleton support has no configuration, so the
  # coherent Gram *is* the total one and the configuration Gram is identically
  # zero. Reported rather than hidden, and it is why this fixture cannot be
  # used to show the three Grams apart.
  coherent <- heterogeneity(plan, estimator = "plug_in", component = "coherent")
  configuration <- heterogeneity(plan, estimator = "plug_in",
    component = "configuration")
  expect_lt(max(abs(coherent$gram - total$gram)), 1e-12)
  expect_lt(max(abs(configuration$gram)), 1e-20)
})

test_that("the group model, not the mean, defines the deviations", {
  plan <- hx_plan(data = data.frame(
    subject = names(hx_sizes),
    age = c(20, 32, 45, 51, 27), stringsAsFactors = FALSE
  ), model = ~ age)
  split <- heterogeneity(plan, estimator = "plug_in")
  expect_identical(split$residual_df, 3L)
  expect_identical(split$receipt$gram$divisor, "residual_df")
  expect_identical(split$receipt$gram$divisor_value, 3)

  deviations <- hx_dense_deviations(plan)
  expect_lt(max(abs(split$gram - tcrossprod(deviations) / 3)), 1e-12)
  # Two of the five eigenvalues are numerical zeros, because the design
  # annihilates a two-dimensional subspace of the participant axis.
  spectrum <- sort(unname(split$spectrum), decreasing = TRUE)
  expect_lt(max(abs(spectrum[4:5])), 1e-10)
})

test_that("a saturated group design is refused rather than answered", {
  sizes <- hx_sizes[1:2]
  subjects <- stats::setNames(lapply(names(sizes), function(id)
    hx_subject(id, sizes[[id]], hx_gains[[id]])), names(sizes))
  plan <- plan_population(subjects,
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      hx_carrier(sizes[[id]])),
    model = ~ group, data = data.frame(subject = names(sizes),
      group = c("a", "b"), stringsAsFactors = FALSE))
  refusal <- catch_refusal(heterogeneity(plan, estimator = "plug_in"))
  expect_identical(refusal$capability, "identified_heterogeneity")
  expect_match(refusal$message, "saturated")
})

test_that("unmatched arguments and unknown nodes are refused", {
  plan <- hx_plan()
  expect_error(heterogeneity(plan, estimator = "plug_in", normalize = TRUE),
    class = "effect_input_error")
  refusal <- catch_refusal(heterogeneity(plan, estimator = "plug_in",
    nodes = "<sink>"))
  expect_identical(refusal$capability, "selected_group_nodes")
  expect_true(any(grepl("group_node_absent", refusal$reasons)))
  expect_error(heterogeneity(plan, estimator = "plug_in", nodes = 9L),
    class = "effect_input_error")
  expect_error(heterogeneity(plan, estimator = "plug_in", nodes = "group1",
    modes = 9L), class = "effect_input_error")
})

test_that("unit_budget is refused on the form route, as it is for the executor", {
  plan <- hx_plan(normalization = "unit_budget")
  refusal <- catch_refusal(heterogeneity(plan, estimator = "plug_in"))
  expect_identical(refusal$capability, "complete_form_normalization")
})


# Presentation -----------------------------------------------------------------

test_that("printing names the estimator and prints the caveat", {
  plan <- hx_plan()
  plug <- heterogeneity(plan, estimator = "plug_in")
  cross <- heterogeneity(plan, nodes = "group1", modes = 1L)

  plug_lines <- capture.output(print(plug))
  expect_identical(plug_lines[[1L]], "<effect_population_heterogeneity>")
  expect_lte(length(plug_lines), 20L)
  expect_true(all(nchar(plug_lines) <= 80L))
  expect_true(any(grepl("plug_in", plug_lines)))
  expect_true(any(grepl("62.7%", plug_lines, fixed = TRUE)))
  expect_true(any(grepl("refused on the plug-in spectrum", plug_lines)))

  cross_lines <- capture.output(print(cross))
  expect_lte(length(cross_lines), 20L)
  expect_true(all(nchar(cross_lines) <= 80L))
  expect_true(any(grepl("cross_fit", cross_lines)))
  expect_true(any(grepl("indefinite by construction", cross_lines)))
  expect_true(any(grepl("interleaved: \\[run1, run3\\]", cross_lines)))
  expect_true(any(grepl("1 at group1", cross_lines)))

  expect_identical(print(plug), plug)
  expect_match(format(cross), "^<effect_population_heterogeneity")

  table <- as.data.frame(cross)
  expect_identical(nrow(table), 25L)
  expect_identical(names(table),
    c("subject", "mode", "eigenvalue", "loading"))
  expect_identical(table$loading[1:5], unname(cross$loadings[, "mode1"]))
})
