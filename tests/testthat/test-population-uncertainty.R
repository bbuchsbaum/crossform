# The group uncertainty layer (ticket E8).
#
# Three things are held here. The between-subject standard error must be the
# group OLS's own, checked against arithmetic done by hand and against base R's
# `lm()`, which shares no code with the population driver. The within-subject
# layer must be exact where it is admitted and absent where it is not, per
# participant and per group column, with no independence assumption anywhere.
# And the two must never be pooled: there is no field holding their sum and no
# table holding both, and this file fails if one appears.

# Fixtures ---------------------------------------------------------------------

# A participant with an error channel, so `rdm_sampling_covariance()` has
# something to read. The native node identifiers are the domain's feature ids,
# which is what the sampling covariance names its blocks by; a transport that
# wants the within layer has to agree with that naming.
pu_subject <- function(id, features, seed) {
  set.seed(seed)
  domain <- abstract_domain(features,
    coordinates = cbind(x = seq_len(features) - 1),
    feature_ids = paste0("f", seq_len(features)), id = id)
  runs <- c("run-1", "run-2", "run-3")
  responses <- stats::setNames(lapply(runs, function(run)
    matrix(stats::rnorm(6 * features), 6L, features)), runs)
  indices <- stats::setNames(lapply(runs, function(run)
    observation_index(paste0(run, "-scan-", 1:6), run)), runs)
  design <- cbind(face = c(1, 0, 0, 1, 0, 0), house = c(0, 1, 0, 0, 1, 0),
    tool = c(0, 0, 1, 0, 0, 1))
  designs <- stats::setNames(lapply(runs, function(run) {
    value <- design
    rownames(value) <- paste0(run, "-scan-", 1:6)
    value
  }), runs)
  target <- diag(3)
  dimnames(target) <- list(c("face", "house", "tool"), colnames(design))
  fit <- estimate_relation(plan_relation(
    study(observations(responses, indices, domain)),
    raw_design_model(designs), raw_effect_map(target),
    observation_model("ols", sampling_unit = "scan")
  ))
  list(
    fit = fit, features = features,
    plan = plan_geometry(fit$relation, compile_frame(voxelwise(), domain),
      cross_partitions(fit$relation, independence = "independent"))
  )
}

# A cheap participant with no error channel, for the between-subject arithmetic
# where no sampling covariance is needed.
pu_plain <- function(id, features, gain) {
  effects <- effect_space(c("face", "house", "tool"), basis_id = "pu:v1")
  domain <- abstract_domain(features,
    coordinates = cbind(x = seq_len(features) - 1),
    feature_ids = paste0("f", seq_len(features)), id = id)
  values <- function(divisor) matrix(
    gain * seq_len(3 * features) / (features * divisor), 3, features,
    dimnames = list(c("face", "house", "tool"), NULL)
  )
  relation <- relation(list(run1 = values(1), run2 = values(1.6)),
    effects = effects, domain = domain)
  plan_geometry(relation, compile_frame(voxelwise(), domain),
    cross_partitions(relation))
}

# Group nodes at 0, 4 and 9: an ordinary many-to-one anatomical assignment.
pu_carrier <- function(features, semantics = "budget") {
  anatomical_transport(
    native_coords = cbind(seq_len(features) - 1),
    group_coords = cbind(c(0, 4, 9)), semantics = semantics,
    native_index = paste0("f", seq_len(features))
  )
}

pu_bank <- function() {
  rbind(`face-house` = c(1, -1, 0), `face-tool` = c(1, 0, -1))
}

pu_population <- function(sizes, gains, model = ~ 1, data = NULL,
                          normalization = "none") {
  plan_population(
    stats::setNames(lapply(names(sizes), function(id)
      pu_plain(id, sizes[[id]], gains[[id]])), names(sizes)),
    lapply(sizes, pu_carrier),
    model = model, data = data, normalization = normalization
  )
}

# Every participant reaches every ordinary group node in the general
# inferential fixture. Variable-coverage inference is a separate estimand and
# is not smuggled into these full-design OLS oracles.
pu_sizes <- c(s01 = 10L, s02 = 11L, s03 = 12L, s04 = 13L, s05 = 14L,
  s06 = 15L)
pu_gains <- c(s01 = 1, s02 = 1.4, s03 = 0.7, s04 = 1.1, s05 = 0.9, s06 = 1.3)


# The between-subject layer ----------------------------------------------------

test_that("the between-subject standard error is the group OLS's own, exactly", {
  fit <- estimate_population(pu_population(pu_sizes, pu_gains), pu_bank())
  layer <- population_uncertainty(fit)

  # An intercept-only group model has a closed form nobody needs a package
  # for: the estimate is the mean of the participants' transported ledgers and
  # its standard error is their sample standard deviation over root N, on
  # `N - 1` degrees of freedom. Both are checked at every group node and query
  # rather than at one.
  expect_true(all(layer$between$residual_df == 6L - 1L))
  values <- fit$values
  hand_estimate <- apply(values, c(1L, 2L), mean)
  hand_se <- apply(values, c(1L, 2L), function(y) stats::sd(y) / sqrt(length(y)))
  expect_lt(max(abs(layer$between$estimate[, , "(Intercept)"] -
    hand_estimate)), 1e-14)
  expect_lt(max(abs(layer$between$se[, , "(Intercept)"] - hand_se)), 1e-14)
  # The sink row of this fixture is identically zero, so its residual scatter
  # is zero and its `t` is `NA` rather than `Inf`: a fit with no scatter is not
  # an infinitely precise estimate.
  estimable <- hand_se > 0
  expect_lt(max(abs((layer$between$t[, , "(Intercept)"] -
    hand_estimate / hand_se)[estimable])), 1e-12)
  expect_true(all(is.na(layer$between$t[, , "(Intercept)"][!estimable])))

  # The interval is the estimate plus or minus the t quantile times the same
  # standard error, and it is nominal: no correction is applied anywhere.
  # Where the statistic is refused the interval is refused with it, rather
  # than collapsing to a zero-width claim about the sink.
  half <- stats::qt(0.975, df = 5L) * hand_se
  expect_lt(max(abs((layer$between$lower[, , "(Intercept)"] -
    (hand_estimate - half))[estimable])), 1e-14)
  expect_lt(max(abs((layer$between$upper[, , "(Intercept)"] -
    (hand_estimate + half))[estimable])), 1e-14)
  expect_true(all(is.na(layer$between$lower[, , "(Intercept)"][!estimable])))
  expect_true(all(is.na(layer$between$upper[, , "(Intercept)"][!estimable])))
})

test_that("a covariate model agrees with base R's lm on every node and query", {
  covariates <- data.frame(
    age = c(21, 34, 27, 45, 31, 38),
    motion = c(0.4, 0.9, 0.2, 1.3, 0.7, 0.5),
    row.names = names(pu_sizes)
  )
  fit <- estimate_population(
    pu_population(pu_sizes, pu_gains, model = ~ age + motion,
      data = covariates),
    pu_bank()
  )
  layer <- population_uncertainty(fit)
  expect_true(all(layer$between$residual_df == 6L - 3L))
  expect_identical(layer$term, c("(Intercept)", "age", "motion"))

  # `lm()` is an independent court: it shares no code with the population
  # driver, and it is asked the same question at every cell rather than at a
  # convenient one.
  worst <- 0
  for (node in dimnames(fit$values)[[1L]]) {
    for (query in dimnames(fit$values)[[2L]]) {
      response <- fit$values[node, query, ]
      reference <- summary(stats::lm(
        response ~ covariates$age + covariates$motion
      ))$coefficients
      if (all(is.finite(reference[, "Std. Error"]))) {
        worst <- max(worst, abs(unname(reference[, "Std. Error"]) -
          unname(layer$between$se[node, query, ])))
      }
    }
  }
  expect_lt(worst, 1e-12)
})

test_that("a saturated group model refuses instead of reporting NA bars", {
  sizes <- pu_sizes[1:3]
  covariates <- data.frame(age = c(21, 34, 27), row.names = names(sizes))
  # Three participants, an intercept and two covariate columns: rank 3 on 3
  # rows leaves no residual degrees of freedom at all.
  covariates$motion <- c(0.4, 0.9, 0.2)
  fit <- estimate_population(
    pu_population(sizes, pu_gains[1:3], model = ~ age + motion,
      data = covariates),
    pu_bank()
  )
  expect_identical(fit$residual_df, 0L)
  refusal <- catch_refusal(population_uncertainty(fit))
  expect_identical(refusal$capability, "between_subject_residual_df")
  expect_identical(refusal$reasons, "saturated_group_model")
  expect_output(print(fit), "not estimable")
})

test_that("the streamed complete form refuses a between-subject layer", {
  plan <- pu_population(pu_sizes[1:3], pu_gains[1:3])
  form <- materialize_population(plan)
  refusal <- catch_refusal(population_uncertainty(form))
  expect_identical(refusal$capability, "population_between_subject_residuals")
  expect_identical(refusal$reasons, "streamed_route_retains_no_residuals")
  expect_match(refusal$remedies, "estimate_population", fixed = TRUE)
})


# The within-subject layer -----------------------------------------------------

test_that("the within layer is exact where a group column has one native row", {
  # Four group nodes on the native grid itself. The participant with four
  # native nodes maps one-to-one, so every group column is fed by exactly one
  # native row; the participant with five sends two native rows into the last
  # group node, which is the general case D8 refuses.
  built <- list(u01 = pu_subject("u01", 4L, 41L), u02 = pu_subject("u02", 5L, 52L))
  carrier <- function(features) anatomical_transport(
    native_coords = cbind(seq_len(features) - 1),
    group_coords = cbind(0:3), semantics = "budget",
    native_index = paste0("f", seq_len(features))
  )
  plan <- plan_population(
    lapply(built, `[[`, "plan"),
    lapply(built, function(value) carrier(value$features))
  )
  bank <- pu_bank()
  uncertainty <- lapply(built, function(value)
    rdm_sampling_covariance(value$plan, value$fit, target = "null",
      at = seq_len(value$features)))
  fit <- estimate_population(plan, bank, uncertainty = uncertainty)
  within <- fit$uncertainty$within

  expect_identical(within$scope, "transported_single_source_column")
  expect_identical(within$assumption,
    "none: the cross-node terms carry weight zero")
  # One-to-one for `u01` at every group node; the sink is fed by nothing and
  # `u02`'s last group node is fed by two native rows.
  expect_identical(unname(within$admitted[, "u01"]),
    c(TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_identical(unname(within$admitted[, "u02"]),
    c(TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_identical(within$admitted_columns, 7L)
  expect_identical(within$source_node[["group1", "u01"]], "f1")
  expect_identical(within$coefficient[["group1", "u01"]], 1)

  # Exact, not approximate: the transported value at a single-source column is
  # `w * z_x`, so its variance is `w^2 Var(z_x)` and `Var(z_x)` is the block
  # D8 supplies untransported.
  blocks <- sampling_covariance(uncertainty$u01, queries = bank)
  diagonal <- sampling_covariance(blocks[["f2"]], "diagonal")
  expect_identical(within$variance["group2", , "u01"],
    stats::setNames(as.numeric(diagonal[rownames(bank)]), rownames(bank)))

  # Where the column mixes native rows there is an absence, never a diagonal
  # approximation standing in for the missing cross terms.
  expect_true(all(is.na(within$variance["group4", , "u02"])))
  expect_true(all(is.na(within$variance["<sink>", , ])))
  expect_true("transported_value_mixes_native_nodes" %in%
    within$refusal$reasons)

  # The untransported blocks E4 ships are unchanged beside it.
  expect_identical(fit$uncertainty$scope, "native_node_marginal")
  expect_identical(fit$uncertainty$transported$capability,
    "transported_sampling_covariance")
})

test_that("the within layer refuses rather than assuming independent nodes", {
  built <- list(u01 = pu_subject("u01", 5L, 41L), u02 = pu_subject("u02", 6L, 52L))
  plan <- plan_population(
    lapply(built, `[[`, "plan"),
    lapply(built, function(value) pu_carrier(value$features))
  )
  bank <- pu_bank()
  uncertainty <- lapply(built, function(value)
    rdm_sampling_covariance(value$plan, value$fit, target = "null",
      at = seq_len(value$features)))
  fit <- estimate_population(plan, bank, uncertainty = uncertainty)
  within <- fit$uncertainty$within

  # An ordinary anatomical parcellation collects several native nodes into
  # each group node, so nothing is admitted. That is the refusal being
  # visible, not the layer being broken: a diagonal sum over overlapping
  # supports under spatially correlated noise underestimates the variance in a
  # known direction, which is worse than an absence.
  expect_false(any(within$admitted))
  expect_identical(within$admitted_columns, 0L)
  expect_true(all(is.na(within$variance)))
  expect_identical(within$status, "refused")
  expect_output(print(fit), "no group column is fed by a single native row")
})

test_that("the within layer is gated before it looks at any group column", {
  built <- list(u01 = pu_subject("u01", 4L, 41L), u02 = pu_subject("u02", 5L, 52L))
  carrier <- function(features) anatomical_transport(
    native_coords = cbind(seq_len(features) - 1),
    group_coords = cbind(0:3), semantics = "budget",
    native_index = paste0("f", seq_len(features))
  )
  bank <- pu_bank()
  uncertainty <- lapply(built, function(value)
    rdm_sampling_covariance(value$plan, value$fit, target = "null",
      at = seq_len(value$features)))

  # `unit_budget` divides by a total read from the same data. Section 4.3
  # records that the standard error of that divisor does not exist, so a
  # transported variance of the ratio would be inventing one.
  shares <- estimate_population(
    plan_population(lapply(built, `[[`, "plan"),
      lapply(built, function(value) carrier(value$features)),
      normalization = "unit_budget"),
    bank, uncertainty = uncertainty
  )
  expect_null(shares$uncertainty$within$admitted)
  expect_identical(shares$uncertainty$within$refusal$reasons,
    "same_data_ratio_normalization")

  # A transport whose native rows are named differently from the nodes the
  # sampling covariance read cannot be bound by position: that would attach
  # one node's error bar to another node's coefficient.
  unnamed <- function(features) anatomical_transport(
    native_coords = cbind(seq_len(features) - 1),
    group_coords = cbind(0:3), semantics = "budget"
  )
  misnamed <- estimate_population(
    plan_population(lapply(built, `[[`, "plan"),
      lapply(built, function(value) unnamed(value$features))),
    bank, uncertainty = uncertainty
  )
  expect_null(misnamed$uncertainty$within$admitted)
  expect_match(misnamed$uncertainty$within$refusal$reasons,
    "^native_node_labels_unaligned:")
  # The untransported blocks still ship: only the carve-out is unavailable.
  expect_identical(names(misnamed$uncertainty$native), c("u01", "u02"))
})


test_that("a density transport carries its data-free divisor into the carve-out", {
  # Under `"density"` the group rows are divided by transported row mass,
  # which is a property of the operator and the declared mass and not of the
  # data. The carried value is therefore still a fixed linear functional of
  # the native values, and the carve-out applies with the divisor folded into
  # the coefficient rather than being ignored.
  built <- list(u01 = pu_subject("u01", 4L, 61L), u02 = pu_subject("u02", 4L, 72L))
  operator <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L, 4L), j = c(1L, 2L, 2L, 3L),
    x = c(0.5, 1, 1, 0.5), dims = c(4L, 3L)
  )
  carrier <- location_transport(
    operator, native_index = paste0("f", 1:4),
    group_index = c("group1", "group2", "group3"), semantics = "density",
    provenance = list(method = "external", details = "density fixture"),
    row_mass = c(4, 1, 1, 0.5)
  )
  bank <- pu_bank()
  uncertainty <- lapply(built, function(value)
    rdm_sampling_covariance(value$plan, value$fit, target = "null",
      at = seq_len(value$features)))
  fit <- estimate_population(
    plan_population(lapply(built, `[[`, "plan"),
      list(u01 = carrier, u02 = carrier)),
    bank, uncertainty = uncertainty
  )
  within <- fit$uncertainty$within

  # `group1` is fed by `f1` alone with weight 0.5 against a transported mass
  # of `0.5 * 4 = 2`, so the coefficient is 0.25 and the variance carries its
  # square. `group3` is fed by `f4` alone with weight 0.5 against a mass of
  # `0.5 * 0.5 = 0.25`, so its coefficient is 2.
  expect_identical(unname(within$admitted[, "u01"]),
    c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(within$coefficient[["group1", "u01"]], 0.25)
  expect_equal(within$coefficient[["group3", "u01"]], 2)
  blocks <- sampling_covariance(uncertainty$u01, queries = bank)
  diagonal <- as.numeric(sampling_covariance(blocks[["f1"]], "diagonal"))
  expect_equal(unname(within$variance["group1", , "u01"]),
    0.25^2 * diagonal)

  # `group2` mixes two native rows and the sink mixes two others, so both are
  # absences rather than diagonal approximations.
  expect_true(all(is.na(within$variance["group2", , ])))
  expect_true(all(is.na(within$variance["<sink>", , ])))
})


# Separation -------------------------------------------------------------------

test_that("the two layers are reported separately and never pooled", {
  built <- list(u01 = pu_subject("u01", 4L, 41L), u02 = pu_subject("u02", 5L, 52L))
  carrier <- function(features) anatomical_transport(
    native_coords = cbind(seq_len(features) - 1),
    group_coords = cbind(0:3), semantics = "budget",
    native_index = paste0("f", seq_len(features))
  )
  bank <- pu_bank()
  uncertainty <- lapply(built, function(value)
    rdm_sampling_covariance(value$plan, value$fit, target = "null",
      at = seq_len(value$features)))
  fit <- estimate_population(plan_population(
    lapply(built, `[[`, "plan"),
    lapply(built, function(value) carrier(value$features))
  ), bank, uncertainty = uncertainty)
  layer <- population_uncertainty(fit)

  expect_identical(layer$layers, c("between_subject", "within_subject"))
  expect_match(layer$separation, "never pooled")
  expect_true(is.list(layer$between) && is.list(layer$within))

  # No field anywhere is the sum, and none is the root of the sum. The check is
  # numerical rather than a name search, because the failure this guards
  # against is a plausible-looking `total` column and not a badly named one.
  pooled <- layer$between$se[, , "(Intercept)"]^2 +
    layer$within$variance[, , "u01"]
  numeric_fields <- Filter(is.numeric, c(layer$between, layer$within))
  for (field in numeric_fields) {
    if (length(field) != length(pooled)) next
    expect_false(isTRUE(all.equal(as.numeric(field), as.numeric(pooled),
      tolerance = 1e-8)))
    expect_false(isTRUE(all.equal(as.numeric(field),
      sqrt(as.numeric(pooled)), tolerance = 1e-8)))
  }

  # There is no combined table either. The two layers have no shared index --
  # the within layer has a participant axis and no term axis -- so one frame
  # would have to invent a column to hold both, and a reader who found them in
  # one frame would add them.
  between_frame <- as.data.frame(layer, layer = "between")
  within_frame <- as.data.frame(layer, layer = "within")
  expect_identical(unique(between_frame$layer), "between_subject")
  expect_identical(unique(within_frame$layer), "within_subject")
  expect_true("term" %in% names(between_frame))
  expect_false("term" %in% names(within_frame))
  expect_true("subject" %in% names(within_frame))
  expect_false("subject" %in% names(between_frame))
  expect_identical(
    within_frame$se[within_frame$admitted],
    sqrt(within_frame$variance[within_frame$admitted])
  )

  # Asking for a within table that does not exist is a refusal naming why,
  # not an empty frame.
  bare <- population_uncertainty(
    estimate_population(pu_population(pu_sizes, pu_gains), pu_bank())
  )
  expect_null(bare$within)
  expect_identical(
    catch_refusal(as.data.frame(bare, layer = "within"))$reasons,
    "within_layer_absent"
  )
})


# Labels -----------------------------------------------------------------------

test_that("every printed and tabulated surface carries the uncalibrated label", {
  fit <- estimate_population(pu_population(pu_sizes, pu_gains), pu_bank())
  layer <- population_uncertainty(fit)

  expect_identical(layer$between$calibration, "uncalibrated")
  expect_true(all(as.data.frame(layer)$calibration == "uncalibrated"))
  expect_identical(fit$uncertainty$between$calibration, "uncalibrated")

  printed <- paste(utils::capture.output(print(layer)), collapse = " ")
  expect_match(printed, "uncalibrated")
  expect_match(printed, "UNCALIBRATED")
  expect_match(printed, "never pooled")
  expect_match(printed, "between-subject")
  expect_match(printed, "within-subject")
  # The measured coverage is quoted beside the label, so a reader is told both
  # that the arithmetic was checked and that the check does not transfer.
  expect_match(printed, "0.885")

  expect_match(format(layer), "^<effect_population_uncertainty:")
  expect_match(format(layer), "uncalibrated")
  expect_output(print(fit),
    "between-subject covariance classical/HC3, cell df 5 to 5")
  expect_output(print(fit), "population_uncertainty")
})

test_that("the reader verb validates its own arguments", {
  fit <- estimate_population(pu_population(pu_sizes, pu_gains), pu_bank())
  expect_error(population_uncertainty(fit, term = "absent"), "term")
  expect_error(population_uncertainty(fit, level = 1), "level")
  expect_error(population_uncertainty(fit, level = 0), "level")
  expect_identical(
    dim(population_uncertainty(fit, term = "(Intercept)")$between$se)[[3L]], 1L
  )
  expect_identical(population_uncertainty(fit, term = 1L)$term, "(Intercept)")
  expect_error(population_uncertainty(fit$receipt), "effect_population_result")

  # The level moves the interval and nothing else.
  wide <- population_uncertainty(fit, level = 0.99)
  narrow <- population_uncertainty(fit, level = 0.9)
  expect_identical(wide$between$se, narrow$between$se)
  finite <- is.finite(wide$between$upper) & is.finite(narrow$between$upper)
  expect_true(any(finite))
  expect_true(all(wide$between$upper[finite] >= narrow$between$upper[finite]))
  expect_false(identical(wide$scientific_plan_id, narrow$scientific_plan_id))
})


# The recorded null coverage ---------------------------------------------------

pu_coverage_summary <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "benchmark-results",
      "population-null-coverage-summary.csv"),
    tryCatch(system.file("extdata", "certification",
      "population-null-coverage-summary.csv", package = "crossform"),
      error = function(condition) ""),
    testthat::test_path("..", "..", "inst", "extdata", "certification",
      "population-null-coverage-summary.csv")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) {
    testthat::skip(paste0(
      "CERTIFICATION ABSENT: no recorded population null-coverage summary is ",
      "available here; run benchmarks/run-population-null-coverage.R."
    ))
  }
  utils::read.csv(candidates[[1L]], stringsAsFactors = FALSE)
}

test_that("the recorded null coverage still holds its declared bounds", {
  summary <- pu_coverage_summary()
  expect_setequal(summary$arm, c("gaussian", "heteroskedastic"))
  expect_setequal(summary$subjects, c(6L, 8L, 12L, 24L))
  expect_true(all(summary$replications >= 2000L))
  expect_identical(summary$residual_df, summary$subjects - 3L)

  # The correctly specified arm. Nominal coverage, and a null `t` that a KS
  # test cannot tell from `t_df`. This is the arm that says the arithmetic is
  # right.
  correct <- summary[summary$arm == "gaussian", ]
  expect_identical(nrow(correct), 4L)
  expect_true(all(correct$null_coverage > 0.93 & correct$null_coverage < 0.97))
  expect_true(all(correct$null_coverage_all_cells > 0.93 &
    correct$null_coverage_all_cells < 0.97))
  expect_true(all(correct$signal_coverage > 0.93 &
    correct$signal_coverage < 0.97))
  expect_true(all(correct$ks_statistic < 0.03))
  expect_true(all(correct$ks_p_value > 0.05))
  # Power rises with the participant count; the null rejection rate does not.
  ordered <- correct[order(correct$subjects), ]
  expect_true(all(diff(ordered$signal_power) > 0))

  # The misspecified arm, held as a floor on a *failure*. These bounds fail if
  # the demonstration the documentation cites stops demonstrating anything --
  # deleting the arm, or softening it until it looks nominal, is the move this
  # catches.
  misspecified <- summary[summary$arm == "heteroskedastic", ]
  expect_identical(nrow(misspecified), 4L)
  expect_true(all(misspecified$null_coverage < 0.94))
  worst <- misspecified[misspecified$subjects == 24L, ]
  expect_lt(worst$null_coverage, 0.90)
  expect_gt(worst$null_rejection, 0.09)
  expect_lt(worst$ks_p_value, 0.01)
  # It gets worse with more participants, which is the whole point: the bias
  # is in the standard error and not in the sample size.
  ordered <- misspecified[order(misspecified$subjects), ]
  expect_lt(ordered$null_coverage[[4L]], ordered$null_coverage[[1L]])

  # The numbers the Rd and the printed note quote are these numbers.
  quoted <- c(`6` = 0.9485, `8` = 0.9500, `12` = 0.9520, `24` = 0.9480)
  expect_lt(max(abs(correct$null_coverage[order(correct$subjects)] -
    quoted[order(as.integer(names(quoted)))])), 0.01)
  expect_lt(abs(worst$null_coverage - 0.8850), 0.01)
})

test_that("the recorded null-coverage artifact is bound to this source", {
  artifact <- certified_artifact(
    "population-null-coverage.rds", "run-population-null-coverage.R"
  )
  expect_identical(artifact$schema_version, 1L)
  expect_identical(artifact$verb, "population_uncertainty")
  expect_gte(artifact$replications, 2000L)
  expect_identical(artifact$level, 0.95)
  expect_identical(artifact$truth$null_coefficient, 0)

  # The runner's own wiring check: the simulated group layer must reproduce
  # `population_uncertainty()` exactly, or the coverage number is evidence
  # about the runner rather than about the shipped verb.
  expect_true(artifact$wiring$identical_se)
  expect_identical(artifact$wiring$se_max_absolute_difference, 0)
  expect_lt(artifact$wiring$lm_max_absolute_difference, 1e-10)

  # The committed CSV is the same table.
  recorded <- pu_coverage_summary()
  expect_identical(nrow(recorded), nrow(artifact$summary))
  expect_lt(max(abs(recorded$null_coverage - artifact$summary$null_coverage)),
    1e-12)
})
