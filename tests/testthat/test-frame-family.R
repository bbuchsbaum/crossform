# Frame families -- `frame_family()`, ticket D2.
#
# `design/conservative-geometry-contract.md` section 3.1 states the law this
# constructor exists to make usable: stacking column-normalized frames with
# weights alpha summing to one gives a conservative frame that conserves
# block by block, sum over the rows of scale s of G = alpha_s * G_Omega. The
# contract's section 11.4 records four gaps this file is the evidence for:
#
#   G1  alpha normalization is a decision, not an assumption -- `frame_family()`
#       refuses weights that do not sum to one and never renormalizes.
#   G2  each member must be column normalized ON ITS OWN; a conserving stack
#       does not imply it, so the check is per block and at construction.
#   G11 node labels must survive stacking -- `$index` and `$specification` did
#       not survive the `rbind()` + `additive_frame()` route.
#   G12 oracle-only claims have no regression protection. Claim 3b
#       (`design/oracles/conservative-multiscale-ledger.R` section O1.c) and
#       claims 7b-7e (`design/oracles/conservative-transport-readiness.R`
#       sections O3.b-O3.f) are promoted here rather than into
#       `test-conservative-geometry-contract.R`, because every one of them is
#       a statement about the family route this ticket owns.
#
# `test-conservative-geometry-contract.R` already asserts claim 3b for a family
# built by hand with `rbind()`. What is new here is the same law through the
# named constructor, plus the per-row metadata that route could not carry.

family_fixture <- function(n = 12L, seed = 4471L, q = 3L,
                           id = "frame-family") {
  domain <- abstract_domain(n, coordinates = cbind(seq_len(n) - 1, 0),
    feature_ids = paste0("vox", seq_len(n)), id = id)
  set.seed(seed)
  mk <- function() {
    value <- matrix(stats::rnorm(q * n), q, n)
    rownames(value) <- letters[seq_len(q)]
    value
  }
  relation_value <- relation(list(run1 = mk(), run2 = mk()), domain = domain)
  list(
    domain = domain, relation = relation_value,
    over = cross_partitions(relation_value)
  )
}

family_scales <- function(domain) {
  list(
    point = compile_frame(voxelwise("conservative"), domain),
    narrow = compile_frame(searchlights(1.01, "conservative"), domain),
    wide = compile_frame(searchlights(2.01, "conservative"), domain)
  )
}

family_total <- function(fixture, frame) {
  geometry_component(
    materialize_geometry(fixture$relation, frame, fixture$over), "total"
  )
}

# Construction ---------------------------------------------------------------

test_that("a frame family is the alpha-scaled row bind of its members", {
  fixture <- family_fixture(id = "frame-family-construct")
  frames <- family_scales(fixture$domain)
  alpha <- c(point = 0.2, narrow = 0.5, wide = 0.3)

  family <- frame_family(point = frames$point, narrow = frames$narrow,
    wide = frames$wide, alpha = alpha)

  expect_s3_class(family, "effect_frame")
  expect_identical(family$representation, "additive_diagonal")
  expect_identical(family$normalization, "conservative")
  expect_identical(nrow(family$weights),
    sum(vapply(frames, function(f) nrow(f$weights), integer(1))))
  expect_identical(ncol(family$weights), fixture$domain$n_features)

  # The rows are the members' rows, scaled by the weight the row records.
  offset <- 0L
  for (name in names(frames)) {
    rows <- offset + seq_len(nrow(frames[[name]]$weights))
    offset <- offset + length(rows)
    expect_equal(
      as.matrix(family$weights[rows, , drop = FALSE]),
      alpha[[name]] * as.matrix(frames[[name]]$weights),
      tolerance = 0, ignore_attr = TRUE
    )
  }

  # Conservative because every member is, and the weights sum to one.
  expect_equal(as.numeric(Matrix::colSums(family$weights)),
    rep(1, fixture$domain$n_features), tolerance = 1e-12)
  expect_true(frame_conservation(family)$conserved)
})

test_that("family members default to equal weights and positional names", {
  domain <- abstract_domain(6L, coordinates = cbind(seq_len(6L) - 1, 0),
    id = "frame-family-defaults")
  family <- frame_family(
    compile_frame(voxelwise("conservative"), domain),
    compile_frame(searchlights(1.01, "conservative"), domain)
  )
  expect_identical(unique(family$index$alpha), 0.5)
  expect_identical(unique(family$index$family), c("frame1", "frame2"))
  expect_identical(names(family$specification$members),
    c("frame1", "frame2"))

  # Named `alpha` is matched to member names, not to argument order.
  reordered <- frame_family(
    point = compile_frame(voxelwise("conservative"), domain),
    narrow = compile_frame(searchlights(1.01, "conservative"), domain),
    alpha = c(narrow = 0.75, point = 0.25)
  )
  expect_identical(reordered$specification$alpha,
    c(point = 0.25, narrow = 0.75))
  expect_identical(unique(reordered$index$alpha[
    reordered$index$family == "point"]), 0.25)
})

test_that("a one-member family is admitted and is the member, relabelled", {
  domain <- abstract_domain(5L, id = "frame-family-single")
  member <- compile_frame(voxelwise("conservative"), domain)
  family <- frame_family(point = member)

  expect_identical(family$specification$alpha, c(point = 1))
  expect_equal(as.matrix(family$weights), as.matrix(member$weights),
    tolerance = 0, ignore_attr = TRUE)
  expect_identical(family$index$node, as.character(member$index$measurement))
})

# Per-row metadata (gap G11, contract section 7.1) ----------------------------

test_that("family rows carry family, node, scale, center and alpha", {
  fixture <- family_fixture(id = "frame-family-metadata")
  domain <- fixture$domain
  labels <- rep(c("roi1", "roi2", "roi3"), each = 4L)
  frames <- list(
    point = compile_frame(voxelwise("conservative"), domain),
    narrow = compile_frame(searchlights(1.01, "conservative"), domain),
    atlas = compile_frame(regions(labels, "conservative"), domain)
  )
  alpha <- c(point = 0.25, narrow = 0.5, atlas = 0.25)
  family <- frame_family(point = frames$point, narrow = frames$narrow,
    atlas = frames$atlas, alpha = alpha)

  index <- family$index
  expect_identical(names(index),
    c("measurement", "family", "node", "scale", "center", "alpha"))
  expect_identical(nrow(index), nrow(family$weights))

  # The identity is unique across the family, which the node label alone is
  # not: `vox1` is a row of both the point and the searchlight member.
  expect_false(anyDuplicated(index$measurement) > 0L)
  expect_true(any(duplicated(index$node)))
  expect_identical(index$measurement,
    paste0(index$family, "::", index$node))

  point <- index[index$family == "point", ]
  narrow <- index[index$family == "narrow", ]
  atlas <- index[index$family == "atlas", ]
  expect_identical(nrow(point), nrow(frames$point$weights))
  expect_identical(nrow(atlas), 3L)

  # Node labels survive from the member's own `$index$measurement`.
  expect_identical(point$node, domain$feature_ids)
  expect_identical(narrow$node, domain$feature_ids)
  expect_identical(atlas$node, c("roi1", "roi2", "roi3"))

  # Scale is the member's own scale parameter, per row rather than frame wide.
  expect_identical(unique(narrow$scale), 1.01)
  expect_true(all(is.na(point$scale)))
  expect_true(all(is.na(atlas$scale)))

  # A center exists exactly where the member anchors its rows at a feature.
  expect_identical(point$center, domain$feature_ids)
  expect_identical(narrow$center,
    as.character(frames$narrow$support_index$node_ids))
  expect_true(all(is.na(atlas$center)))

  # The applied weight travels with the row.
  expect_identical(unique(point$alpha), 0.25)
  expect_identical(unique(narrow$alpha), 0.5)
  expect_identical(unique(atlas$alpha), 0.25)

  # The member specifications survive on the family specification.
  specification <- family$specification
  expect_identical(specification$kind, "frame_family")
  expect_identical(specification$normalization, "conservative")
  expect_identical(names(specification$members), names(frames))
  expect_identical(specification$members$narrow$specification,
    frames$narrow$specification)
  expect_identical(specification$members$atlas$declared_normalization,
    "conservative")
  expect_identical(specification$alpha, alpha)
})

test_that("a declared member without an index falls back to row positions", {
  domain <- abstract_domain(4L, id = "frame-family-declared")
  declared <- additive_frame(
    matrix(c(1, 0, 0, 0, 0, 1, 1, 1), 2L, 4L, byrow = TRUE),
    normalization = "conservative", domain = domain
  )
  family <- frame_family(declared = declared,
    point = compile_frame(voxelwise("conservative"), domain),
    alpha = c(0.5, 0.5))

  rows <- family$index[family$index$family == "declared", ]
  expect_identical(rows$node, c("1", "2"))
  expect_true(all(is.na(rows$scale)))
  expect_true(all(is.na(rows$center)))
})

# Refusals -------------------------------------------------------------------

test_that("family weights that do not sum to one are refused, not rescaled", {
  # Gap G1. The per-block law reads the weight that was actually applied, so
  # a silent renormalization would change the estimand of every row.
  domain <- abstract_domain(6L, coordinates = cbind(seq_len(6L) - 1, 0),
    id = "frame-family-alpha")
  point <- compile_frame(voxelwise("conservative"), domain)
  narrow <- compile_frame(searchlights(1.01, "conservative"), domain)

  expect_error(
    frame_family(point = point, narrow = narrow, alpha = c(1, 1)),
    "must sum to one", class = "effect_input_error"
  )
  expect_error(
    frame_family(point = point, narrow = narrow, alpha = c(1, 1)),
    "never renormalized for you", class = "effect_input_error"
  )
  expect_error(
    frame_family(point = point, narrow = narrow, alpha = c(0.5, 0.4)),
    "sums to 0.9", class = "effect_input_error"
  )
  # A nonpositive weight is not a way to drop a member: it makes rows of zero
  # mass, which is not a share of anything.
  expect_error(
    frame_family(point = point, narrow = narrow, alpha = c(0, 1)),
    "must be positive", class = "effect_input_error"
  )
  expect_error(
    frame_family(point = point, narrow = narrow, alpha = c(-0.5, 1.5)),
    "must be positive", class = "effect_input_error"
  )
  expect_error(
    frame_family(point = point, narrow = narrow, alpha = 1),
    "one per family member", class = "effect_input_error"
  )
  expect_error(
    frame_family(point = point, narrow = narrow,
      alpha = c(point = 0.5, other = 0.5)),
    "must name every family member exactly once", class = "effect_input_error"
  )

  # Weights that do sum to one are accepted at the stated tolerance.
  expect_s3_class(
    frame_family(point = point, narrow = narrow, alpha = c(1, 2) / 3),
    "effect_frame"
  )
})

test_that("a member that is not column normalized on its own is refused", {
  # Gap G2. The stacked frame can conserve while no single block does, so the
  # check has to be per block: `sum over x in s of w_xv = 1` for every feature
  # separately is strictly stronger than the stack summing to one.
  domain <- abstract_domain(6L, coordinates = cbind(seq_len(6L) - 1, 0),
    id = "frame-family-blocks")
  point <- compile_frame(voxelwise("conservative"), domain)

  expect_error(
    frame_family(point = point,
      narrow = compile_frame(searchlights(1.01), domain), alpha = c(0.5, 0.5)),
    "not column normalized on its own", class = "effect_input_error"
  )

  # A member that leaves a feature uncovered cannot be column normalized at
  # all, and the message says so rather than reporting a bare deviation.
  uncovered <- additive_frame(
    matrix(c(1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0), 2L, 6L, byrow = TRUE),
    normalization = "none", domain = domain
  )
  expect_error(
    frame_family(point = point, partial = uncovered, alpha = c(0.5, 0.5)),
    "1 feature carry no mass at all", class = "effect_input_error"
  )

  # The case the stack cannot detect: two members whose alpha-weighted sum has
  # unit columns while neither member does. The stacked frame conserves, so a
  # check on the stack passes it, and the per-scale law is nonetheless false
  # for both blocks.
  heavy <- additive_frame(rbind(c(1.5, 1.5, 1.5, 0.5, 0.5, 0.5)),
    normalization = "none", domain = domain)
  light <- additive_frame(rbind(c(0.5, 0.5, 0.5, 1.5, 1.5, 1.5)),
    normalization = "none", domain = domain)
  expect_equal(as.numeric(colSums(rbind(
    0.5 * as.matrix(heavy$weights), 0.5 * as.matrix(light$weights)
  ))), rep(1, 6L), tolerance = 1e-12)
  expect_true(frame_conservation(additive_frame(
    rbind(0.5 * as.matrix(heavy$weights), 0.5 * as.matrix(light$weights)),
    normalization = "conservative", domain = domain
  ))$conserved)
  expect_gt(max(abs(as.numeric(Matrix::colSums(heavy$weights)) - 1)), 0.1)
  expect_error(
    frame_family(heavy = heavy, light = light, alpha = c(0.5, 0.5)),
    "not column normalized on its own", class = "effect_input_error"
  )

  # A metric-folded member carries the metric diagonal in its columns, so it
  # has no fixed unit budget for alpha to divide.
  small <- abstract_domain(4L, id = "frame-family-folded")
  precision <- noise_precision(diag(c(0.5, 2, 1.25, 3)), small)
  conservative <- compile_frame(
    regions(c("a", "a", "b", "b"), normalization = "conservative"), small
  )
  folded <- crossform:::.metric_additive_frame(
    conservative,
    crossform:::.geometry_metric_schedule(conservative, precision)
  )
  expect_error(
    frame_family(folded = folded,
      point = compile_frame(voxelwise("conservative"), small),
      alpha = c(0.5, 0.5)),
    "diagonal metric folded into its weights", class = "effect_input_error"
  )
})

test_that("a family refuses members bound to different neural domains", {
  first <- abstract_domain(6L, coordinates = cbind(seq_len(6L) - 1, 0),
    id = "frame-family-domain-a")
  second <- abstract_domain(6L, coordinates = cbind(seq_len(6L) - 1, 0),
    id = "frame-family-domain-b")

  expect_error(
    frame_family(
      point = compile_frame(voxelwise("conservative"), first),
      narrow = compile_frame(searchlights(1.01, "conservative"), second),
      alpha = c(0.5, 0.5)
    ),
    "must be compiled against the same one", class = "effect_contract_error"
  )

  # Equal feature counts do not make two domains one domain.
  wider <- abstract_domain(7L, coordinates = cbind(seq_len(7L) - 1, 0),
    id = "frame-family-domain-a")
  expect_error(
    frame_family(
      point = compile_frame(voxelwise("conservative"), first),
      wide = compile_frame(voxelwise("conservative"), wider),
      alpha = c(0.5, 0.5)
    ),
    class = "effect_contract_error"
  )
})

test_that("a family refuses malformed members and duplicate identities", {
  domain <- abstract_domain(5L, id = "frame-family-refusals")
  point <- compile_frame(voxelwise("conservative"), domain)

  expect_error(frame_family(), "at least one compiled frame",
    class = "effect_input_error")
  expect_error(frame_family(point, "not a frame"),
    "must be a compiled `effect_frame`", class = "effect_input_error")
  expect_error(frame_family(point = point, point = point, alpha = c(0.5, 0.5)),
    "names more than one member", class = "effect_input_error")
  expect_error(
    frame_family(point = point, other = point, alpha = c(0.5, 0.5),
      normalization = "local"),
    "must be \"conservative\"", class = "effect_input_error"
  )
  expect_error(frame_family(point = point, tolerance = -1),
    class = "effect_input_error")
})

# The conservation law, end to end -------------------------------------------

test_that("an alpha-weighted family conserves total, overall and by block", {
  # Contract claim 3b, promoted from `conservative-multiscale-ledger.R` section
  # O1.c (gap G12). The comparator is the UNNORMALIZED whole-brain operator:
  # `whole_brain()` defaults to "local", which divides by the feature count.
  fixture <- family_fixture(id = "frame-family-conservation")
  domain <- fixture$domain
  frames <- family_scales(domain)
  alpha <- c(point = 0.2, narrow = 0.5, wide = 0.3)
  family <- frame_family(point = frames$point, narrow = frames$narrow,
    wide = frames$wide, alpha = alpha)

  global <- drop(family_total(fixture,
    compile_frame(whole_brain("none"), domain)))
  total <- family_total(fixture, family)

  # Whole family: sum over ALL rows is the global geometry, since sum(alpha) = 1.
  expect_equal(colSums(total), global, tolerance = 1e-12)
  expect_lt(max(abs(colSums(total) - global)), 1e-12)

  # Block by block: each scale's rows carry exactly alpha_s of the budget. The
  # rows are found through `$index$family`, not through a positional offset --
  # that is what the metadata is for.
  worst <- 0
  for (name in names(frames)) {
    rows <- which(family$index$family == name)
    expect_identical(length(rows), nrow(frames[[name]]$weights))
    observed <- colSums(total[rows, , drop = FALSE])
    expect_equal(observed, alpha[[name]] * global, tolerance = 1e-12)
    worst <- max(worst, max(abs(observed - alpha[[name]] * global)))
  }
  expect_lt(worst, 1e-12)

  # A fixed query reads the same statement as one number per scale, which is
  # the panel the contract's section 3.1 forbids presenting as a finding: the
  # column is alpha times a constant.
  query <- crossform:::.svec_symmetric(diag(3L))
  budget <- sum(query * global)
  energy <- vapply(names(frames), function(name) {
    rows <- which(family$index$family == name)
    sum(total[rows, , drop = FALSE] %*% query)
  }, numeric(1))
  expect_equal(unname(energy), unname(alpha * budget), tolerance = 1e-12)
})

test_that("frame_conservation certifies a family block by block", {
  fixture <- family_fixture(id = "frame-family-certificate")
  frames <- family_scales(fixture$domain)
  alpha <- c(point = 0.2, narrow = 0.5, wide = 0.3)
  family <- frame_family(point = frames$point, narrow = frames$narrow,
    wide = frames$wide, alpha = alpha)

  report <- frame_conservation(family)
  expect_true(report$conserved)
  expect_identical(report$component, "total")
  expect_identical(report$normalization, "conservative")

  members <- report$members
  expect_s3_class(members, "data.frame")
  expect_identical(names(members),
    c("family", "alpha", "measurements", "max_deviation", "conserved"))
  expect_identical(members$family, names(frames))
  expect_identical(members$alpha, unname(alpha))
  expect_identical(members$measurements,
    vapply(frames, function(f) nrow(f$weights), integer(1)),
    ignore_attr = TRUE)
  expect_true(all(members$conserved))
  expect_lt(max(members$max_deviation), 1e-12)

  # An ordinary compiled frame is not a family and reports no member block.
  expect_null(frame_conservation(frames$narrow)$members)
})

# Per-row metadata reaches the result (gap G11, contract section 7.1) ---------

test_that("family metadata survives a geometry evaluation", {
  fixture <- family_fixture(id = "frame-family-result")
  domain <- fixture$domain
  frames <- family_scales(domain)
  family <- frame_family(point = frames$point, narrow = frames$narrow,
    wide = frames$wide, alpha = c(0.2, 0.5, 0.3))

  view <- contrast_energy(
    plan_geometry(fixture$relation, family, fixture$over),
    c(a = 1, b = -1, c = 0)
  )
  values <- as.data.frame(view)

  # The result indexes by the family's own identifiers, not by row position.
  expect_identical(as.character(values$measurement), family$index$measurement)
  expect_false(all(values$measurement == as.character(seq_len(nrow(values)))))

  # Which makes the result joinable back to the per-row metadata, so a view
  # can be grouped by scale, by center or by member.
  joined <- merge(values, family$index, by = "measurement", sort = FALSE)
  expect_identical(nrow(joined), nrow(values))
  expect_identical(sort(unique(joined$family)), sort(names(frames)))
  expect_identical(sort(unique(joined$scale[joined$family == "wide"])), 2.01)

  by_scale <- vapply(split(joined$total, joined$family), sum, numeric(1))
  expect_identical(sort(names(by_scale)), sort(names(frames)))
  expect_equal(sum(by_scale), sum(values$total), tolerance = 1e-12)
})

test_that("the bare rbind route is still provenance blind", {
  # Claim 7e, promoted from `conservative-transport-readiness.R` section O3.f
  # (gap G12), in its post-D2 form: the hand-stacked route remains numerically
  # sound and remains provenance blind, which is the whole reason
  # `frame_family()` exists. If `additive_frame()` ever starts carrying an
  # index of its own, this test is the one that says so.
  fixture <- family_fixture(8L, id = "frame-family-provenance")
  domain <- fixture$domain
  point <- compile_frame(voxelwise("conservative"), domain)
  narrow <- compile_frame(searchlights(1.01, "conservative"), domain)

  hand_stacked <- additive_frame(
    rbind(0.5 * as.matrix(point$weights), 0.5 * as.matrix(narrow$weights)),
    normalization = "conservative", domain = domain
  )
  expect_true(frame_conservation(hand_stacked)$conserved)
  expect_null(hand_stacked$index)
  expect_null(hand_stacked$specification)

  hand_values <- as.data.frame(contrast_energy(
    plan_geometry(fixture$relation, hand_stacked, fixture$over),
    c(a = 1, b = -1, c = 0)
  ))
  expect_identical(as.character(hand_values$measurement),
    as.character(seq_len(nrow(hand_values))))

  # The constructor gives the same operator with the provenance attached.
  family <- frame_family(point = point, narrow = narrow, alpha = c(0.5, 0.5))
  expect_equal(as.matrix(family$weights), as.matrix(hand_stacked$weights),
    tolerance = 0, ignore_attr = TRUE)
  expect_false(is.null(family$index))
  expect_false(is.null(family$specification))
})

# Promoted transport-readiness oracle claims (gap G12) ------------------------

test_that("geometry covariance has rank at most N-1 and the subject Gram recovers it", {
  # Claim 7b, promoted from `conservative-transport-readiness.R` section O3.b.
  # The packed codec is a Frobenius isometry, so a population covariance over
  # N subjects can be read from the N x N subject Gram instead of the P x P
  # covariance. For realistic q that is the difference between a feasible
  # population layer and an infeasible one.
  set.seed(20260817)
  q <- 5L
  n_subjects <- 6L
  width <- q * (q + 1L) / 2L
  symmetric <- function() {
    value <- matrix(stats::rnorm(q * q), q, q)
    0.5 * (value + t(value))
  }
  packed <- t(vapply(seq_len(n_subjects),
    function(s) crossform:::.svec_symmetric(symmetric()), numeric(width)))
  centered <- sweep(packed, 2L, colMeans(packed))

  covariance <- crossprod(centered) / (n_subjects - 1L)
  gram <- tcrossprod(centered) / (n_subjects - 1L)
  expect_identical(dim(covariance), c(as.integer(width), as.integer(width)))
  expect_identical(dim(gram), c(n_subjects, n_subjects))
  expect_identical(qr(covariance)$rank, n_subjects - 1L)

  spectrum <- function(x) {
    sort(eigen(x, symmetric = TRUE, only.values = TRUE)$values,
      decreasing = TRUE)
  }
  top <- seq_len(n_subjects - 1L)
  expect_equal(spectrum(covariance)[top], spectrum(gram)[top],
    tolerance = 1e-12)
  expect_lt(max(abs(spectrum(covariance)[n_subjects:width])), 1e-12)
})

test_that("row-stochastic transport preserves the budget and a missing sink does not", {
  # Claim 7c, promoted from `conservative-transport-readiness.R` section O3.c.
  # Conservation is what makes partial coverage detectable: a detection map has
  # no budget to check against, so the same omission leaves no numerical trace.
  set.seed(20260817)
  assignment <- c(1, 1, 2, 2, 2, 3, 3, 4, 4, 4)
  n_native <- length(assignment)
  n_group <- 4L
  transport <- matrix(0, n_native, n_group)
  transport[cbind(seq_len(n_native), assignment)] <- 1
  expect_equal(rowSums(transport), rep(1, n_native), tolerance = 0)

  contribution <- stats::rnorm(n_native)
  budget <- sum(contribution)
  expect_equal(sum(crossprod(transport, contribution)), budget,
    tolerance = 1e-12)

  # Partial coverage without a sink: the rows stop being stochastic and the
  # mass is gone rather than accounted for.
  partial <- transport
  partial[assignment == 4L, ] <- 0
  lost <- budget - sum(crossprod(partial, contribution))
  expect_gt(abs(lost / budget), 0.01)

  # With an explicit sink column the budget closes and the unmapped mass is
  # visible in the sink.
  sunk <- cbind(partial, sink = 0)
  sunk[rowSums(partial) == 0, n_group + 1L] <- 1
  expect_equal(rowSums(sunk), rep(1, n_native), tolerance = 0)
  expect_equal(sum(crossprod(sunk, contribution)), budget, tolerance = 1e-12)
  expect_equal(drop(crossprod(sunk, contribution))[[n_group + 1L]], lost,
    tolerance = 1e-12)
})

test_that("budget and density transport semantics differ, and neither is a default", {
  # Claim 7d, promoted from `conservative-transport-readiness.R` section O3.d.
  # Two subjects carrying the same total evidence contribute unequally per
  # territory under budget semantics, purely because their native frames have
  # different resolutions.
  subjects <- list(
    fine = list(contribution = rep(1, 10) / 10,
      assignment = c(1, 1, 2, 2, 2, 3, 3, 4, 4, 4)),
    coarse = list(contribution = rep(1, 4) / 4, assignment = 1:4)
  )
  n_group <- 4L
  read <- function(subject) {
    transport <- matrix(0, length(subject$contribution), n_group)
    transport[cbind(seq_along(subject$contribution), subject$assignment)] <- 1
    budget <- as.numeric(crossprod(transport, subject$contribution))
    list(budget = budget, density = budget / colSums(transport))
  }
  fine <- read(subjects$fine)
  coarse <- read(subjects$coarse)

  # Both subjects carry budget one, and both semantics preserve that total.
  expect_equal(sum(fine$budget), 1, tolerance = 1e-12)
  expect_equal(sum(coarse$budget), 1, tolerance = 1e-12)

  # The per-territory numbers nonetheless disagree between the two readings.
  expect_equal(fine$budget, c(0.2, 0.3, 0.2, 0.3), tolerance = 1e-12)
  expect_equal(coarse$budget, rep(0.25, 4), tolerance = 1e-12)
  expect_equal(fine$density, rep(0.1, 4), tolerance = 1e-12)
  expect_equal(coarse$density, rep(0.25, 4), tolerance = 1e-12)
  expect_gt(max(abs(fine$budget - coarse$budget)), 0.01)
  expect_gt(max(abs(fine$density - coarse$density)), 0.01)
})
