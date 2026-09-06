# Pins the numeric claims of `vignettes/population-form.Rmd`.
#
# The vignette asserts every number it prints, so a broken identity stops it
# knitting -- but vignettes are not built by `R CMD check --no-build-vignettes`
# and are not built at all by a bare `devtools::test()`. This file is the
# safety net: it rebuilds the article's fixture verbatim and re-asserts the
# quantities the prose names, so a change that would silently move a printed
# number in the article fails here first.
#
# It deliberately pins *values and bounds*, never printed output. The article
# prints several objects, and the print methods are free to change without
# this file having an opinion.
#
# The fixture is built once and memoized. Every other population test file
# rebuilds its plan per `test_that()`; this one is a duplicate of work those
# files already do at full detail, so it buys its own place in the suite by
# being cheap.

vig_effects <- function() {
  effect_space(c("face", "house", "tool"),
    basis_id = "population-vignette:v1")
}

vig_subject <- function(id, features, gain = 1, tilt = 0, runs = 4L,
                        normalization = "conservative") {
  domain <- abstract_domain(features,
    coordinates = cbind(x = seq_len(features) - 1),
    feature_ids = paste0("f", seq_len(features)), id = id)
  consensus <- outer(c(0.9, -0.9, 0), rep(1, features))
  odd <- outer(c(0, 0.8, -0.8), rep(1, features))
  block <- function(k) {
    set.seed(1000L * k + sum(as.integer(charToRaw(id))) + features)
    values <- matrix(gain * stats::rnorm(3L * features), 3L, features,
      dimnames = list(c("face", "house", "tool"), NULL))
    values + consensus + tilt * odd
  }
  relation <- relation(
    stats::setNames(lapply(seq_len(runs), block), paste0("run", seq_len(runs))),
    effects = vig_effects(), domain = domain
  )
  plan_geometry(relation,
    compile_frame(voxelwise(normalization = normalization), domain),
    cross_partitions(relation))
}

vig_carrier <- function(features, semantics = "budget", radius = 2) {
  anatomical_transport(
    native_coords = cbind(seq_len(features) - 1),
    group_coords = cbind(c(0, 5, 11)),
    semantics = semantics, radius = radius,
    native_index = paste0("f", seq_len(features)))
}

# Full ordinary-node coverage keeps the article's OLS and commutation examples
# on one planned population; the radius still leaves genuine sink territory.
vig_sizes <- c(s01 = 12L, s02 = 14L, s03 = 16L, s04 = 13L, s05 = 15L,
  s06 = 17L)
vig_gains <- c(s01 = 1, s02 = 1.3, s03 = 0.8, s04 = 1.1, s05 = 0.9, s06 = 1.2)
vig_tilts <- c(s01 = 0, s02 = 0, s03 = 0, s04 = 0, s05 = 0, s06 = 4)

vig_bank <- function() {
  rbind(`face-house` = c(1, -1, 0), `house-tool` = c(0, 1, -1))
}

vig_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    subjects <- stats::setNames(lapply(names(vig_sizes), function(id)
      vig_subject(id, vig_sizes[[id]], vig_gains[[id]], vig_tilts[[id]])),
      names(vig_sizes))
    transports <- stats::setNames(lapply(names(vig_sizes), function(id)
      vig_carrier(vig_sizes[[id]])), names(vig_sizes))
    plan <- plan_population(subjects, transports)
    cache <<- list(subjects = subjects, transports = transports, plan = plan,
      fit = estimate_population(plan, vig_bank()))
    cache
  }
})


# Section 1 --- the three objects ----------------------------------------------

test_that("the article's transport carries a signed ledger and needs its sink", {
  P <- rbind(
    f1 = c(anterior = 1,   posterior = 0),
    f2 = c(anterior = 0.6, posterior = 0.4),
    f3 = c(anterior = 0,   posterior = 0.7),
    f4 = c(anterior = 0,   posterior = 0)
  )
  carrier <- external_transport(P, semantics = "budget",
    provenance = list(details = "partial-volume warp, atlas-tool 2.1"))
  expect_identical(dim(carrier$matrix), c(4L, 3L))
  expect_identical(nrow(carrier$group_index), 2L)

  ledger <- c(f1 = 1.5, f2 = -0.5, f3 = 2, f4 = 0.25)
  carried <- transport_values(carrier, ledger)
  expect_lt(abs(sum(carried) - sum(ledger)), 1e-12 * sum(abs(ledger)))
  # The prose names 0.85 as the signed mass deleting the sink would lose.
  expect_equal(sum(ledger) - sum(carried[c("anterior", "posterior")]), 0.85,
    tolerance = 1e-12)

  # Density is a different estimand and conserves nothing.
  dense <- external_transport(P, semantics = "density",
    provenance = list(details = "the same operator, read as a density"))
  expect_gt(abs(sum(transport_values(dense, ledger)) - sum(ledger)), 1e-6)

  # And a functional operator must name the partitions that built it.
  refusal <- catch_refusal(external_transport(P, semantics = "budget",
    provenance = list(method = "functional",
      details = "hyperalignment on the analysis runs")))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "cross_fit_provenance")
  expect_identical(refusal$namespace, "location_transport")
  expect_identical(refusal$reasons, "cross_fit_partitions_not_declared")
})

test_that("the article's plan is the estimand its prose describes", {
  fixture <- vig_fixture()
  plan <- fixture$plan
  expect_identical(plan$semantics, "budget")
  expect_identical(plan$normalization, "none")
  expect_identical(plan$fit$evaluation_order, "transport_then_fit")
  expect_true(plan$fit$commuting)
  expect_identical(nrow(plan$group_index), 3L)
  expect_true(all(plan$subject_index$conserved))
  expect_true(all(plan$subject_index$sink_territory > 0))

  # The conservative gate, by name.
  loose <- fixture$subjects
  loose$s01 <- vig_subject("s01", vig_sizes[["s01"]], vig_gains[["s01"]],
    normalization = "local")
  gate <- catch_refusal(plan_population(loose, fixture$transports))
  expect_identical(gate$capability, "conservative_subject_geometry")
  expect_true("normalization_not_conservative:s01:local" %in% gate$reasons)

  # The group model binds by name, not by position.
  aged <- plan_population(fixture$subjects, fixture$transports, model = ~ age,
    data = data.frame(age = c(24, 31, 27, 44, 38, 22),
      row.names = c("s02", "s01", "s04", "s03", "s06", "s05")))
  expect_identical(aged$model$columns, c("(Intercept)", "age"))
  expect_identical(aged$model$rank, 2L)
  expect_identical(rownames(aged$model$matrix), names(vig_sizes))
  expect_identical(unname(aged$model$matrix[, "age"]),
    c(31, 24, 44, 27, 22, 38))
  expect_false(identical(plan$scientific_plan_id, aged$scientific_plan_id))

  # And precision weighting is refused at plan construction.
  weighted <- catch_refusal(plan_population(fixture$subjects,
    fixture$transports, normalization = "precision_weighted"))
  expect_identical(weighted$capability, "precision_weighted_normalization")
  expect_true("per_subject_budget_variance_unavailable" %in% weighted$reasons)
})


# Section 2 --- the commutation ------------------------------------------------

test_that("the article's order reversal is a real one, and it commutes", {
  fixture <- vig_fixture()
  plan <- fixture$plan
  fit <- fixture$fit
  expect_identical(fit$basis, "query_bank")

  # The plumbing non-check: exact by construction, and the article says so.
  one <- estimate_population(plan, rbind(`face-house` = c(1, -1, 0)))
  plumbing <- max(vapply(names(plan$subjects), function(id) {
    native <- contrast_energy(plan$subjects[[id]], c(1, -1, 0))$total
    max(abs(transport_values(plan$transport[[id]], native) -
      one$values[, "face-house", id]))
  }, numeric(1)))
  expect_identical(plumbing, 0)

  # The real reversal: transport -> fit -> query, against query -> transport
  # -> fit. `population-form-v1` section 11 states 1e-12.
  form <- materialize_population(plan)
  expect_identical(form$basis, "complete_form")
  edges <- as.data.frame(rdm(form))
  edge <- edges[edges$left == "face" & edges$right == "house", ]
  node_ids <- as.character(fit$index$node)
  expect_identical(nrow(edge), length(node_ids))
  expect_identical(anyDuplicated(edge$node), 0L)
  expect_lt(max(abs(edge$estimate[match(node_ids, edge$node)] -
    fit$coefficients[node_ids, "face-house", "(Intercept)"])), 1e-12)

  # Not vacuous: the numbers being compared are real.
  expect_gt(max(abs(fit$coefficients)), 1)
})


# Section 3 --- conservation ---------------------------------------------------

test_that("the article's budget certificate holds, and the sink is load-bearing", {
  fit <- vig_fixture()$fit
  native_total <- fit$receipt$native_total
  carried <- apply(fit$values, c("subject", "query"), sum)
  l1 <- apply(abs(fit$values), c("subject", "query"), sum)

  expect_true(fit$receipt$budget$asserted)
  expect_identical(fit$receipt$budget$scale, "relative_to_ledger_l1_norm")
  expect_lt(fit$receipt$budget$max_relative_deviation, 1e-12)
  expect_lt(max(abs(carried - native_total) / l1), 1e-12)

  # Deleting the sink breaks it, by a real amount about a real coverage
  # failure -- roughly a third of the ledger's L1 norm on this fixture.
  group_only <- apply(fit$values[!fit$index$sink, , , drop = FALSE],
    c("subject", "query"), sum)
  expect_gt(max(abs(group_only - native_total) / l1), 1e-3)

  # The sink is fitted like any other row, and is in budget units.
  expect_identical(fit$index$units, rep("budget", 4L))
  expect_true(all(is.finite(fit$coefficients["<sink>", , ])))
  expect_gt(max(abs(fit$coefficients["<sink>", , "(Intercept)"])), 0)
  expect_identical(fit$ledger, "transported_total")
})


# Section 4 --- heterogeneity --------------------------------------------------

test_that("the article's planted odd participant is found by mode 1", {
  plan <- vig_fixture()$plan
  het <- heterogeneity(plan, estimator = "cross_fit")

  expect_identical(het$estimator, "cross_fit")
  expect_identical(het$space, "packed_form")
  expect_identical(het$residual_df, 5L)

  # Indefinite by construction, and reported as-is (section 6.2).
  expect_identical(sum(het$spectrum < -1e-8), 2L)
  expect_identical(het$latent$negative_modes, 2L)
  expect_gt(het$latent$moved_share, 0)

  # `s06` is the only participant carrying the planted extra direction, and
  # nothing told the estimator that.
  leading <- het$loadings[, 1L]
  others <- leading[names(leading) != "s06"]
  expect_identical(names(which.max(abs(leading))), "s06")
  expect_gt(max(abs(leading)), 4 * max(abs(others)))
  expect_gt(het$latent$cumulative[[1L]], 0.85)
})

test_that("the article's plug-in comparison keeps directions and refuses eigenvalues", {
  plan <- vig_fixture()$plan
  cross_fit <- heterogeneity(plan, estimator = "cross_fit")
  plug <- heterogeneity(plan, estimator = "plug_in")

  # Direction: the same participant, the same mode, to better than 0.99.
  expect_identical(names(which.max(abs(plug$loadings[, 1L]))), "s06")
  expect_gt(abs(sum(cross_fit$loadings[, 1L] * plug$loadings[, 1L])), 0.99)

  # Size: the plug-in books within-subject noise as heterogeneity. The +62.7%
  # magnitude is the contract's Monte Carlo number and stays there; one draw
  # of a trace estimates nothing, so only the direction is pinned here.
  expect_gt(sum(diag(plug$gram)), sum(diag(cross_fit$gram)))

  # And the package enforces section 6.3 rather than advising it.
  expect_null(plug$latent)
  refusal <- plug$receipt$latent_refusal
  expect_identical(refusal$capability, "plug_in_spectrum_functionals")
  expect_true(
    "plug_in_trace_inflated_62_7_percent_on_the_contract_fixture" %in%
      refusal$reasons)
})


# Section 5 --- the two error bars, and the count that is neither --------------

test_that("the article's between-subject layer is the group OLS's own", {
  fit <- vig_fixture()$fit
  uncertainty <- population_uncertainty(fit)
  between <- uncertainty$between

  expect_true(all(between$residual_df == 5L))
  hand_estimate <- apply(fit$values, c("node", "query"), mean)
  hand_se <- apply(fit$values, c("node", "query"),
    function(y) stats::sd(y) / sqrt(length(y)))
  expect_lt(max(abs(between$estimate[, , "(Intercept)"] - hand_estimate)),
    1e-12)
  expect_lt(max(abs(between$se[, , "(Intercept)"] - hand_se)), 1e-12)

  # The label the article explains at length.
  expect_identical(between$calibration, "uncalibrated")

  # The two layers are never pooled, and no field holds their sum.
  expect_identical(uncertainty$layers,
    c("between_subject", "within_subject"))
  expect_null(uncertainty$total)
  expect_identical(
    intersect(names(uncertainty), c("total", "pooled", "combined")),
    character(0))
  expect_identical(
    sort(unique(as.data.frame(uncertainty, layer = "between")$layer)),
    "between_subject")

  # The article's own six-participant transport admits no within layer: every
  # group node collects several native nodes.
  expect_null(uncertainty$within)
})

test_that("the article's HC3 sensitivity layer is fully identified", {
  fit <- vig_fixture()$fit
  classical <- population_uncertainty(fit)$between
  hc3 <- population_uncertainty(fit, estimator = "HC3")$between
  ratio <- hc3$se / classical$se

  expect_identical(hc3$estimator, "HC3")
  expect_true("heteroskedasticity_robust_sandwich" %in% hc3$assumptions)
  expect_identical(hc3$residual_df, fit$coverage$residual_df)
  expect_identical(dim(hc3$covariance), c(4L, 2L, 1L, 1L))
  expect_true(any(is.finite(ratio)))
})

test_that("the article's wild bootstrap preserves the subject unit", {
  fit <- vig_fixture()$fit
  wild <- population_wild_bootstrap(
    fit, "(Intercept)", null = 0, replicates = 199L, seed = 20260821L
  )
  table <- as.data.frame(wild)

  expect_identical(dim(wild$weights), c(6L, 199L))
  expect_identical(rownames(wild$weights), fit$coverage$planned_subjects)
  expect_true(all(wild$successful_replicates[!fit$index$sink, ] == 199L))
  expect_identical(wild$conditioning, fit$coverage$conditioning)
  expect_true(all(c("observed_t", "p_value", "monte_carlo_se",
    "successful_replicates", "status") %in% names(table)))
  expect_identical(sum(table$sink), nrow(fit$queries))
  expect_true(all(table$status[table$sink] %in% c("estimated", "refused")))
})

test_that("the article's within-subject layer is exact where admitted, absent elsewhere", {
  wu_subject <- function(id, features, seed) {
    set.seed(seed)
    domain <- abstract_domain(features,
      coordinates = cbind(x = seq_len(features) - 1),
      feature_ids = paste0("f", seq_len(features)), id = id)
    runs <- c("run-1", "run-2", "run-3")
    scans <- function(run) paste0(run, "-scan-", 1:6)
    design <- cbind(face = c(1, 0, 0, 1, 0, 0), house = c(0, 1, 0, 0, 1, 0),
      tool = c(0, 0, 1, 0, 0, 1))
    target <- diag(3)
    dimnames(target) <- list(c("face", "house", "tool"), colnames(design))
    fit <- estimate_relation(plan_relation(
      study(observations(
        stats::setNames(lapply(runs, function(run)
          matrix(stats::rnorm(6 * features), 6L, features)), runs),
        stats::setNames(lapply(runs, function(run)
          observation_index(scans(run), run)), runs),
        domain)),
      raw_design_model(stats::setNames(lapply(runs, function(run) {
        value <- design; rownames(value) <- scans(run); value
      }), runs)),
      raw_effect_map(target),
      observation_model("ols", sampling_unit = "scan")
    ))
    list(fit = fit, features = features,
      plan = plan_geometry(fit$relation, compile_frame(voxelwise(), domain),
        cross_partitions(fit$relation, independence = "independent")))
  }

  built <- list(u01 = wu_subject("u01", 4L, 41L),
    u02 = wu_subject("u02", 4L, 52L), u03 = wu_subject("u03", 5L, 63L))
  grid_carrier <- function(features) anatomical_transport(
    native_coords = cbind(seq_len(features) - 1), group_coords = cbind(0:3),
    semantics = "budget", native_index = paste0("f", seq_len(features)))

  plan <- plan_population(lapply(built, `[[`, "plan"),
    lapply(built, function(value) grid_carrier(value$features)))
  bank <- rbind(`face-house` = c(1, -1, 0), `face-tool` = c(1, 0, -1))
  fit <- estimate_population(plan, bank,
    uncertainty = lapply(built, function(value)
      rdm_sampling_covariance(value$plan, value$fit, target = "null",
        at = seq_len(value$features))))
  within <- population_uncertainty(fit)$within

  expect_identical(within$scope, "transported_single_source_column")
  expect_identical(within$assumption,
    "none: the cross-node terms carry weight zero")
  expect_identical(unname(within$admitted[, "u01"]), c(rep(TRUE, 4L), FALSE))
  expect_identical(unname(within$admitted[, "u02"]), c(rep(TRUE, 4L), FALSE))
  expect_identical(unname(within$admitted[, "u03"]),
    c(rep(TRUE, 3L), FALSE, FALSE))
  expect_identical(within$admitted_columns, 11L)
  expect_true(all(is.na(within$variance["<sink>", , ])))
  expect_true(all(is.na(within$variance["group4", , "u03"])))
  expect_true("transported_value_mixes_native_nodes" %in%
    within$refusal$reasons)
})

test_that("the article's prevalence stays descriptive and reports its coverage", {
  fixture <- vig_fixture()
  coverage_transports <- stats::setNames(lapply(names(vig_sizes), function(id) {
    anatomical_transport(
      native_coords = cbind(seq_len(vig_sizes[[id]]) - 1),
      group_coords = cbind(c(0, 5, 16)), semantics = "budget", radius = 2,
      native_index = paste0("f", seq_len(vig_sizes[[id]]))
    )
  }), names(vig_sizes))
  coverage_plan <- plan_population(fixture$subjects, coverage_transports,
    coverage_policy = "available_at_node")
  fit <- estimate_population(coverage_plan, vig_bank())
  before <- fit$uncertainty
  prevalence <- population_prevalence(fit, coverage_floor = 6L)

  expect_identical(prevalence$layer, "latent_descriptive")
  expect_identical(prevalence$reading,
    "latent descriptive layer; not for inference")
  expect_identical(prevalence$reference, 0.5)
  expect_identical(prevalence$receipt$prevalence$inference, "none_derivable")

  flatten <- function(x) {
    if (!is.list(x)) return(character())
    c(names(x), unlist(lapply(x, flatten), use.names = FALSE))
  }
  expect_identical(
    intersect(flatten(prevalence[c("sign", "alignment", "coverage")]),
      c("se", "t", "lower", "upper", "p_value", "p", "statistic", "level",
        "confidence")),
    character(0))

  # Constructing one does not reach into the inferential layer.
  expect_identical(fit$uncertainty, before)

  # The deliberate coverage example moves group3 to x = 16, so only the
  # larger native frames stand behind its available-at-node target.
  expect_identical(prevalence$coverage$floor, 6L)
  expect_identical(prevalence$coverage$below_floor, "group3")
  expect_lt(prevalence$coverage$minimum[["group3"]], 6L)
  expect_true(all(is.finite(prevalence$sign$fraction)))
  expect_true(all(prevalence$sign$fraction >= 0 &
    prevalence$sign$fraction <= 1))
})
