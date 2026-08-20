# `plan_population()` --- the group-level estimand of `population-form-v1`.
#
# The fixtures are deliberately tiny and fully deterministic: a population
# plan reads no data, so nothing here needs a realistic dataset, and every
# assertion below is about which combinations are admitted, what identity they
# seal, and what the refusals say. The two participants have *different*
# native frame sizes on purpose --- heterogeneous native frames are the case
# `population-form-v1` §3 says forces transport to precede the fit, and a
# fixture where every subject had the same frame would never exercise it.

population_effects <- function() {
  effect_space(c("face", "house"), basis_id = "population-conditions:v1")
}

# The cheapest compiled conservative geometry plan: a deterministic relation
# over two runs, a `voxelwise()` frame (conservative by default), and the
# cross-partition pairing.
population_subject_plan <- function(id, features, scale = 1,
                                    normalization = "conservative") {
  domain <- abstract_domain(features,
    coordinates = cbind(x = seq_len(features) - 1),
    feature_ids = paste0("f", seq_len(features)), id = id)
  values <- function(divisor) {
    matrix(scale * seq_len(2 * features) / (features * divisor), 2, features,
      dimnames = list(c("face", "house"), NULL))
  }
  relation <- relation(
    list(run1 = values(1), run2 = values(2)),
    effects = population_effects(), domain = domain
  )
  plan_geometry(
    relation,
    compile_frame(voxelwise(normalization = normalization), domain),
    cross_partitions(relation)
  )
}

# Two group nodes at x = 0 and x = 3, reached by every native node.
population_carrier <- function(features, group = c(0, 3), ...) {
  anatomical_transport(
    native_coords = cbind(seq_len(features) - 1),
    group_coords = cbind(group), semantics = "budget", ...
  )
}

population_fixture <- function(features = c(s01 = 5L, s02 = 6L), ...) {
  labels <- names(features)
  list(
    subjects = stats::setNames(lapply(labels, function(label)
      population_subject_plan(label, features[[label]])), labels),
    transport = stats::setNames(lapply(labels, function(label)
      population_carrier(features[[label]], ...)), labels)
  )
}

population_plan_fixture <- function(...) {
  fixture <- population_fixture()
  plan_population(fixture$subjects, fixture$transport, ...)
}


# Construction ----------------------------------------------------------------

test_that("a population plan seals the subjects, transports and group nodes", {
  fixture <- population_fixture()
  plan <- plan_population(fixture$subjects, fixture$transport)

  expect_s3_class(plan, "effect_population_plan")
  expect_identical(names(plan$subjects), c("s01", "s02"))
  expect_identical(names(plan$transport), c("s01", "s02"))
  expect_identical(nrow(plan$group_index), 2L)
  expect_identical(plan$semantics, "budget")
  expect_identical(plan$normalization, "none")
  expect_false(plan$allow_nonconservative)

  # The evaluation order is recorded, not chosen: with heterogeneous native
  # frames there is no common node axis, so fit-then-transport is not even
  # definable (`population-form-v1` §3).
  expect_identical(plan$fit, list(
    kind = "ols", weights = "subject_constant", commuting = TRUE,
    evaluation_order = "transport_then_fit"
  ))

  # The default group model is the group mean, and its matrix has one row per
  # participant carrying that participant's name.
  expect_identical(plan$model$columns, "(Intercept)")
  expect_identical(rownames(plan$model$matrix), c("s01", "s02"))
  expect_identical(plan$model$rank, 1L)
  expect_identical(nrow(plan$data), 2L)

  expect_true(all(plan$subject_index$conserved))
  expect_identical(plan$subject_index$declared_normalization,
    c("conservative", "conservative"))
  expect_identical(plan$subject_index$measurements, c(5L, 6L))
  expect_identical(plan$subject_index$plan_id,
    c(fixture$subjects$s01$scientific_plan_id,
      fixture$subjects$s02$scientific_plan_id))

  expect_match(plan$scientific_plan_id, "^population-sha256:")
  expect_silent(crossform:::.validate_population_plan(plan, deep = TRUE))
})

test_that("subjects are stored in a canonical order whatever order they arrive in", {
  fixture <- population_fixture()
  forward <- plan_population(fixture$subjects, fixture$transport)
  reversed <- plan_population(
    fixture$subjects[c("s02", "s01")], fixture$transport[c("s02", "s01")]
  )
  expect_identical(names(reversed$subjects), c("s01", "s02"))
  expect_identical(forward$scientific_plan_id, reversed$scientific_plan_id)
  expect_identical(forward$signature, reversed$signature)
})

test_that("a group model is built once, with a pivoted QR on the plan", {
  fixture <- population_fixture(c(s01 = 5L, s02 = 6L, s03 = 7L, s04 = 4L))
  covariates <- data.frame(
    age = c(24, 31, 27, 44), group = c("a", "b", "a", "b"),
    row.names = c("s01", "s02", "s03", "s04")
  )
  plan <- plan_population(fixture$subjects, fixture$transport,
    model = ~ age + group, data = covariates)

  design <- plan$model$matrix
  expect_identical(plan$model$columns,
    c("(Intercept)", "age", "groupb"))
  expect_identical(plan$model$terms, c("age", "group"))
  expect_identical(plan$model$intercept, 1L)
  expect_identical(plan$model$rank, 3L)

  # The stored factorization is the group fit's whole numerical content under
  # OLS: R'R reproduces the pivoted cross-product exactly.
  factorization <- plan$model$qr
  expect_s3_class(factorization, "qr")
  pivoted <- design[, factorization$pivot, drop = FALSE]
  expect_equal(crossprod(qr.R(factorization)), crossprod(pivoted),
    tolerance = 1e-12)
  expect_equal(unname(qr.Q(factorization) %*% qr.R(factorization)),
    unname(pivoted), tolerance = 1e-12)

  # The design rows follow the canonical subject order, and the covariates
  # follow their participant rather than their position.
  expect_identical(rownames(design), c("s01", "s02", "s03", "s04"))
  expect_identical(unname(design[, "age"]), c(24, 31, 27, 44))
})

test_that("design rows bind to participants by name, in any row order", {
  fixture <- population_fixture(c(s01 = 5L, s02 = 6L))
  shuffled <- data.frame(subject = c("s02", "s01"), age = c(31, 24))
  plan <- plan_population(fixture$subjects, fixture$transport,
    model = ~ age, data = shuffled)
  expect_identical(unname(plan$model$matrix[, "age"]), c(24, 31))
  expect_identical(rownames(plan$data), c("s01", "s02"))
})

test_that("unit_budget is admitted and recorded as a plan-identity field", {
  fixture <- population_fixture()
  none <- plan_population(fixture$subjects, fixture$transport)
  share <- plan_population(fixture$subjects, fixture$transport,
    normalization = "unit_budget")
  expect_identical(share$normalization, "unit_budget")
  expect_false(identical(none$scientific_plan_id, share$scientific_plan_id))
})


# Identity --------------------------------------------------------------------

test_that("the scientific identity is stable across constructions", {
  first <- population_plan_fixture()
  second <- population_plan_fixture()
  expect_identical(first$scientific_plan_id, second$scientific_plan_id)
  expect_identical(first$signature, second$signature)

  # Stable across sessions, not merely within one: the digest is a SHA-256
  # over a pinned serialization, so a plan that was saved and reloaded --
  # which is what a second session holds -- carries and re-derives the same
  # identity rather than merely remembering it.
  reloaded <- unserialize(serialize(first, NULL))
  expect_identical(reloaded$scientific_plan_id, first$scientific_plan_id)
  crossform:::.reset_validation_memo()
  expect_silent(crossform:::.validate_population_plan(reloaded, deep = TRUE))
})

test_that("every identity ingredient moves the identity", {
  base <- population_plan_fixture()
  fixture <- population_fixture()

  # 1. A subject's own geometry estimand.
  altered <- fixture
  altered$subjects$s02 <- population_subject_plan("s02", 6L, scale = 3)
  moved_subject <- plan_population(altered$subjects, altered$transport)
  expect_false(identical(base$scientific_plan_id,
    moved_subject$scientific_plan_id))

  # 2. The transport. Two transports are two estimands (§1.5).
  altered <- fixture
  altered$transport$s01 <- population_carrier(5L, group = c(1, 3))
  altered$transport$s02 <- population_carrier(6L, group = c(1, 3))
  moved_transport <- plan_population(altered$subjects, altered$transport)
  expect_false(identical(base$scientific_plan_id,
    moved_transport$scientific_plan_id))

  # 3. The group model formula.
  covariates <- data.frame(age = c(24, 31), row.names = c("s01", "s02"))
  with_age <- plan_population(fixture$subjects, fixture$transport,
    model = ~ age, data = covariates)
  expect_false(identical(base$scientific_plan_id, with_age$scientific_plan_id))

  # 4. The covariate *content*, at a fixed formula.
  older <- plan_population(fixture$subjects, fixture$transport, model = ~ age,
    data = data.frame(age = c(24, 32), row.names = c("s01", "s02")))
  expect_false(identical(with_age$scientific_plan_id,
    older$scientific_plan_id))

  # A covariate `data` carries but the model does not use is not part of the
  # estimand, so it leaves the scientific identity alone -- and still moves
  # the physical signature, because the plan stores the table.
  extra <- plan_population(fixture$subjects, fixture$transport, model = ~ age,
    data = data.frame(age = c(24, 31), handedness = c("l", "r"),
      row.names = c("s01", "s02")))
  expect_identical(with_age$scientific_plan_id, extra$scientific_plan_id)
  expect_false(identical(with_age$signature, extra$signature))

  # 5. The budget normalization.
  share <- plan_population(fixture$subjects, fixture$transport,
    normalization = "unit_budget")
  expect_false(identical(base$scientific_plan_id, share$scientific_plan_id))

  # 6. The non-conservative admission, which is an estimand change and not a
  # relaxed setting on the same one -- so it moves the identity even when the
  # frames it would have admitted are all conservative anyway.
  admitted <- plan_population(fixture$subjects, fixture$transport,
    allow_nonconservative = TRUE)
  expect_false(identical(base$scientific_plan_id,
    admitted$scientific_plan_id))
})

test_that("the compute policy moves the signature and not the estimand", {
  fixture <- population_fixture()
  base <- plan_population(fixture$subjects, fixture$transport)
  blocked <- plan_population(fixture$subjects, fixture$transport,
    compute = compute_policy(block_features = 8L))
  expect_identical(base$scientific_plan_id, blocked$scientific_plan_id)
  expect_false(identical(base$signature, blocked$signature))
})

test_that("a tampered plan fails its own signature", {
  plan <- population_plan_fixture()

  forged <- plan
  forged$normalization <- "unit_budget"
  expect_error(crossform:::.validate_population_plan(forged, deep = FALSE),
    class = "effect_contract_error")

  forged <- plan
  forged$subject_index$conserved <- c(FALSE, FALSE)
  expect_error(crossform:::.validate_population_plan(forged, deep = FALSE),
    class = "effect_contract_error")

  forged <- plan
  forged$semantics <- NULL
  expect_error(crossform:::.validate_population_plan(forged, deep = FALSE),
    class = "effect_input_error")

  expect_error(crossform:::.validate_population_plan(list()),
    class = "effect_input_error")
})


# Refusals --------------------------------------------------------------------

test_that("participants must share one experimental space", {
  fixture <- population_fixture()
  domain <- abstract_domain(6L, coordinates = cbind(x = seq_len(6) - 1),
    feature_ids = paste0("f", seq_len(6)), id = "s02")
  values <- matrix(seq_len(12) / 6, 2, 6,
    dimnames = list(c("face", "object"), NULL))
  other <- relation(list(run1 = values, run2 = values / 2),
    effects = effect_space(c("face", "object"), basis_id = "other:v1"),
    domain = domain)
  fixture$subjects$s02 <- plan_geometry(other,
    compile_frame(voxelwise(), domain), cross_partitions(other))

  refusal <- catch_refusal(
    plan_population(fixture$subjects, fixture$transport)
  )
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "common_experimental_space")
  expect_identical(refusal$namespace, "population_plans")
  expect_true(any(grepl("^disagreeing_subject:s02:", refusal$reasons)))
  expect_match(refusal$remedies, "effect_space", all = FALSE)
})

test_that("non-conservative subject frames are refused unless declared", {
  fixture <- population_fixture()
  fixture$subjects$s01 <- population_subject_plan("s01", 5L,
    normalization = "local")

  refusal <- catch_refusal(
    plan_population(fixture$subjects, fixture$transport)
  )
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "conservative_subject_geometry")
  expect_true("normalization_not_conservative:s01:local" %in% refusal$reasons)
  expect_match(refusal$remedies, "allow_nonconservative", all = FALSE)

  # The escape hatch exists, and is recorded rather than silent.
  admitted <- plan_population(fixture$subjects, fixture$transport,
    allow_nonconservative = TRUE)
  expect_true(admitted$allow_nonconservative)
  expect_identical(admitted$subject_index$declared_normalization,
    c("local", "conservative"))
})

test_that("a metric-folded conservative frame counts as conservative", {
  # Folding a diagonal metric into a conservative frame leaves weights that
  # are no longer column-normalized, so the frame declares
  # `normalization = "none"` -- the truth about the operator it carries --
  # and records the declared normalization under `$metric_folded`. Reading
  # `$normalization` here would refuse a conservative frame for carrying a
  # metric, so the plan reads `frame_conservation()$declared_normalization`.
  features <- 5L
  domain <- abstract_domain(features, coordinates = cbind(x = 0:4),
    feature_ids = paste0("f", seq_len(features)), id = "s01")
  base <- compile_frame(voxelwise(), domain)
  precision <- noise_precision(diag(c(0.5, 2, 1.25, 3, 1.5)), domain)
  folded <- crossform:::.metric_additive_frame(
    base, crossform:::.geometry_metric_schedule(base, precision)
  )
  expect_identical(folded$normalization, "none")
  report <- frame_conservation(folded)
  expect_identical(report$declared_normalization, "conservative")
  expect_true(report$metric_folded)

  values <- function(divisor) {
    matrix(seq_len(2 * features) / (features * divisor), 2, features,
      dimnames = list(c("face", "house"), NULL))
  }
  relation <- relation(list(run1 = values(1), run2 = values(2)),
    effects = population_effects(), domain = domain)
  fixture <- population_fixture()
  fixture$subjects$s01 <- plan_geometry(relation, folded,
    cross_partitions(relation))

  plan <- plan_population(fixture$subjects, fixture$transport)
  expect_identical(plan$subject_index$metric_folded, c(TRUE, FALSE))
  expect_true(all(plan$subject_index$conserved))
  expect_identical(plan$subject_index$declared_normalization,
    c("conservative", "conservative"))
})

test_that("learned-metric subject plans are refused for now", {
  setup <- metric_learning_setup()
  learned <- plan_geometry(
    setup$fixture$fit, setup$fixture$frame, setup$over,
    metric = shrinkage_precision(0.2),
    residual_workspace_bytes = setup$budgets$wider
  )
  expect_identical(learned$metric_schedule$kind, "learned_local_before_frame")

  fixture <- population_fixture()
  fixture$subjects$s01 <- learned
  fixture$transport$s01 <- population_carrier(learned$measurements)

  refusal <- catch_refusal(
    plan_population(fixture$subjects, fixture$transport)
  )
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "fixed_metric_population")
  expect_true("learned_metric_subject_plan:s01" %in% refusal$reasons)
})

test_that("subject and transport names must match, one for one", {
  fixture <- population_fixture()

  renamed <- fixture$transport
  names(renamed) <- c("s01", "s99")
  refusal <- catch_refusal(plan_population(fixture$subjects, renamed))
  expect_identical(refusal$capability, "aligned_subject_transports")
  expect_true("subject_without_transport:s02" %in% refusal$reasons)
  expect_true("transport_without_subject:s99" %in% refusal$reasons)

  expect_error(plan_population(fixture$subjects, unname(fixture$transport)),
    class = "effect_input_error")
  expect_error(
    plan_population(fixture$subjects, list(s01 = 1, s02 = 2)),
    class = "effect_input_error"
  )
})

test_that("a transport's native rows must count the subject's measurements", {
  fixture <- population_fixture()
  fixture$transport$s02 <- population_carrier(9L)

  refusal <- catch_refusal(
    plan_population(fixture$subjects, fixture$transport)
  )
  expect_identical(refusal$capability, "aligned_subject_transports")
  expect_true("native_size_mismatch:s02:9_rows_for_6_measurements" %in%
    refusal$reasons)
})

test_that("transports must land on one group node set under one semantics", {
  fixture <- population_fixture()
  fixture$transport$s02 <- population_carrier(6L, group = c(0, 3, 5))
  refusal <- catch_refusal(
    plan_population(fixture$subjects, fixture$transport)
  )
  expect_identical(refusal$capability, "shared_group_nodes")
  expect_true("group_index_disagrees:s02" %in% refusal$reasons)

  mixed <- population_fixture()
  mixed$transport$s02 <- anatomical_transport(
    native_coords = cbind(seq_len(6) - 1), group_coords = cbind(c(0, 3)),
    semantics = "density"
  )
  refusal <- catch_refusal(plan_population(mixed$subjects, mixed$transport))
  expect_identical(refusal$capability, "shared_group_nodes")
  expect_true("semantics_disagrees:s02:density" %in% refusal$reasons)
})

test_that("the group design must be identified, and the refusal names the columns", {
  fixture <- population_fixture(c(s01 = 5L, s02 = 6L, s03 = 7L))
  covariates <- data.frame(
    age = c(24, 31, 27), months = c(288, 372, 324),
    row.names = c("s01", "s02", "s03")
  )
  refusal <- catch_refusal(
    plan_population(fixture$subjects, fixture$transport,
      model = ~ age + months, data = covariates)
  )
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "identified_group_model")
  expect_true("model_rank_deficient" %in% refusal$reasons)
  expect_true("aliased_column:months" %in% refusal$reasons)
  expect_match(refusal$message, "rank deficient")

  # Fewer participants than columns is its own reason, reported alongside the
  # aliased names rather than instead of them.
  pair <- population_fixture(c(s01 = 5L, s02 = 6L))
  short <- catch_refusal(
    plan_population(pair$subjects, pair$transport, model = ~ age + group,
      data = data.frame(age = c(24, 31), group = c("a", "b"),
        row.names = c("s01", "s02")))
  )
  expect_identical(short$capability, "identified_group_model")
  expect_true("subjects_fewer_than_columns:2_of_3" %in% short$reasons)
  expect_true("aliased_column:groupb" %in% short$reasons)
})

test_that("group covariates must be present, complete, and one-sided", {
  fixture <- population_fixture()

  absent <- catch_refusal(
    plan_population(fixture$subjects, fixture$transport, model = ~ age)
  )
  expect_identical(absent$capability, "identified_group_model")
  expect_true("model_variable_absent:age" %in% absent$reasons)

  incomplete <- catch_refusal(
    plan_population(fixture$subjects, fixture$transport, model = ~ age,
      data = data.frame(age = c(24, NA), row.names = c("s01", "s02")))
  )
  expect_true("model_variable_incomplete:age" %in% incomplete$reasons)
  expect_match(incomplete$message, "will not drop a participant")

  two_sided <- catch_refusal(
    plan_population(fixture$subjects, fixture$transport, model = age ~ 1,
      data = data.frame(age = c(24, 31), row.names = c("s01", "s02")))
  )
  expect_true("model_formula_has_response" %in% two_sided$reasons)
})

test_that("design rows must name their participant rather than trusting order", {
  fixture <- population_fixture()
  refusal <- catch_refusal(
    plan_population(fixture$subjects, fixture$transport, model = ~ age,
      data = data.frame(age = c(24, 31)))
  )
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "aligned_subject_rows")
  expect_true("data_rows:unlabelled" %in% refusal$reasons)

  wrong <- catch_refusal(
    plan_population(fixture$subjects, fixture$transport, model = ~ age,
      data = data.frame(age = c(24, 31), row.names = c("s01", "s99")))
  )
  expect_identical(wrong$capability, "aligned_subject_rows")
})

test_that("precision_weighted is refused until the cross-node covariance exists", {
  fixture <- population_fixture()
  refusal <- catch_refusal(
    plan_population(fixture$subjects, fixture$transport,
      normalization = "precision_weighted")
  )
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "precision_weighted_normalization")
  expect_identical(refusal$namespace, "population_plans")
  expect_true("per_subject_budget_variance_unavailable" %in% refusal$reasons)
  expect_true("cross_node_sampling_covariance_route_absent" %in%
    refusal$reasons)
  expect_match(refusal$remedies, "gap G8", all = FALSE)
  expect_match(refusal$message, "cross_measurement")

  # The cause this gate names is live, not historical, and the two records
  # agree. `sampling_covariance()`'s query bank landed while this file was
  # written and is *not* the missing object: the cross-measurement scope a
  # conserved-budget variance would have to take is itself refused, for want
  # of a spatial model and because conservation is a law about estimates
  # rather than about their uncertainty. Asserted against the scope guard
  # directly -- it is a pure function of the scope string, so this grounds the
  # gate without standing up a sampling fixture and without asserting anything
  # about `sampling_covariance()`'s own arguments.
  grounding <- catch_refusal(
    crossform:::.require_local_sampling_scope("cross_measurement")
  )
  expect_s3_class(grounding, "effect_capability_refusal")
  expect_identical(grounding$capability, "cross_node_sampling_covariance")
  expect_true(all(c("spatial_covariance_model_unavailable",
    "conservation_gives_no_uncertainty") %in% grounding$reasons))
  expect_true(crossform:::.require_local_sampling_scope("measurement"))

  # It stays in the closed set --- the set *is* the plan identity, so an
  # argument that silently dropped the mode would be a different contract.
  expect_identical(
    eval(formals(plan_population)$normalization),
    c("none", "unit_budget", "precision_weighted")
  )
})

test_that("rectangular cross-axis subject plans have no symmetric form to pack", {
  domain <- abstract_domain(5L, coordinates = cbind(x = seq_len(5) - 1),
    feature_ids = paste0("f", seq_len(5)), id = "s01")
  left_values <- matrix(seq_len(10) / 5, 2, 5,
    dimnames = list(c("face", "house"), NULL))
  right_values <- matrix(seq_len(15) / 5, 3, 5,
    dimnames = list(c("a", "b", "c"), NULL))
  left <- relation(list(run1 = left_values, run2 = left_values / 2),
    effects = population_effects(), domain = domain)
  right <- relation(list(run1 = right_values, run2 = right_values / 2),
    effects = effect_space(c("a", "b", "c"), basis_id = "right:v1"),
    domain = domain)
  rectangular <- plan_geometry(left, compile_frame(voxelwise(), domain),
    pairing("run1", "run2", directed = TRUE), right = right)

  fixture <- population_fixture()
  fixture$subjects$s01 <- rectangular

  refusal <- catch_refusal(
    plan_population(fixture$subjects, fixture$transport)
  )
  expect_identical(refusal$capability, "symmetric_packed_population")
  expect_true("rectangular_subject_plan:s01" %in% refusal$reasons)
})

test_that("a population needs at least two participants and typed inputs", {
  fixture <- population_fixture()
  expect_error(
    plan_population(fixture$subjects["s01"], fixture$transport["s01"]),
    class = "effect_input_error"
  )
  expect_error(plan_population(unname(fixture$subjects), fixture$transport),
    class = "effect_input_error")
  expect_error(
    plan_population(list(s01 = 1, s02 = 2), fixture$transport),
    class = "effect_input_error"
  )
  expect_error(
    plan_population(fixture$subjects, fixture$transport, model = "~ 1"),
    class = "effect_input_error"
  )
  expect_error(
    plan_population(fixture$subjects, fixture$transport, data = 1:2),
    class = "effect_input_error"
  )
  expect_error(
    plan_population(fixture$subjects, fixture$transport,
      allow_nonconservative = NA),
    class = "effect_input_error"
  )
})


# Printing --------------------------------------------------------------------

test_that("a population plan prints its estimand-bearing choices", {
  fixture <- population_fixture(c(s01 = 5L, s02 = 6L, s03 = 7L))
  plan <- plan_population(fixture$subjects, fixture$transport,
    model = ~ age, normalization = "unit_budget",
    data = data.frame(age = c(24, 31, 27), row.names = c("s01", "s02", "s03")))
  # Digests are scrubbed: a printed identity is checked for *shape* here, and
  # its value is pinned by the identity tests above. A snapshot holding the
  # literal hash would break on any upstream change to a subject plan's own
  # estimand, which is not what this test is about.
  scrub <- function(lines) {
    sub("((?:estimand|signature):\\s+(?:population-)?sha256:)[0-9a-f]+",
      "\\1<digest>", lines)
  }
  expect_snapshot(print(plan), transform = scrub)
  expect_identical(format(plan),
    "<effect_population_plan: 3 subjects, 2 group nodes, ~age, unit_budget>")
})

test_that("the printed sink line reports partial coverage", {
  # A radius-limited transport leaves the far native nodes unmapped, and the
  # unmapped territory is a property of the operator alone -- printable from
  # a plan that has read no data.
  fixture <- population_fixture(c(s01 = 5L, s02 = 6L))
  fixture$transport$s01 <- population_carrier(5L, group = c(0, 1), radius = 1)
  fixture$transport$s02 <- population_carrier(6L, group = c(0, 1), radius = 1)
  plan <- plan_population(fixture$subjects, fixture$transport)
  expect_true(all(plan$subject_index$sink_territory > 0))
  expect_snapshot(print(plan), transform = function(lines) {
    sub("((?:estimand|signature):\\s+(?:population-)?sha256:)[0-9a-f]+",
      "\\1<digest>", lines)
  })
})
