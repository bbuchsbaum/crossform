# The group uncertainty layer --------------------------------------------------
#
# Two error bars can be put on a transported group ledger, they answer
# different questions, and `population-form-v1` section 7 requires them to be
# reported apart:
#
#   between-subject   how far the participants scatter around the group fit.
#                     Read off the group OLS at each group node and query, with
#                     that fit's own residual degrees of freedom. It needs no
#                     sampling law and no error channel: the participants are
#                     the replicates.
#
#   within-subject    the sampling variance of one participant's *transported*
#                     value, carried through the transport from D8's
#                     per-native-node blocks. It needs the covariance between
#                     native nodes, which D8 refuses, so it exists only for a
#                     group column fed by exactly one native row --- where the
#                     cross terms carry weight zero and the carve-out is exact.
#
# They are not summands. The between-subject residual already contains whatever
# measurement error survived into each participant's transported value, so
# adding the within-subject variance to it would double-count the part that is
# common to both and still miss the covariance structure a variance-components
# model would need. This file therefore reports two blocks and never a third,
# and `as.data.frame()` emits one layer at a time for the same reason.
#
# Neither layer is calibrated. The between-subject `t` is a ratio, its null
# distribution has been *measured* against `t_df` under a correctly specified
# group model (`benchmarks/run-population-null-coverage.R`), and that
# measurement says nothing about a real study, where the group model is
# misspecified to an unknown degree and the transport is heterogeneous across
# participants in ways section 7.5's diagnostics can show but not correct. The
# label stays `"uncalibrated"` in every field and on every printed line.

# The arithmetic, in one place.
#
# `benchmarks/run-population-null-coverage.R` calls this function directly on
# simulated group-layer responses, so the coverage number it records is a
# measurement of the shipped estimator rather than of a re-implementation of
# it. Everything here is a slice-wise scalar operation over the
# `node x query x term` array, so a non-finite response column --- an
# unresolved density node, an unadmitted `unit_budget` share --- carries `NA`
# through to its own standard error and to nothing else.
.population_between_statistics <- function(coefficients, residuals,
                                           residual_df, unscaled, level) {
  shape <- dim(coefficients)
  terms <- dimnames(coefficients)[[3L]]
  # One residual variance per (node, query): the same design served every
  # column of the fit, so the only thing that varies along those two axes is
  # the scatter of the participants about it.
  squares <- apply(residuals^2, c(1L, 2L), sum)
  residual_sd <- sqrt(squares / residual_df)
  leverage <- vapply(terms, function(term) unscaled[term, term], numeric(1))
  se <- array(rep(as.numeric(residual_sd), times = shape[[3L]]), shape) *
    rep(sqrt(leverage), each = shape[[1L]] * shape[[2L]])
  dimnames(se) <- dimnames(coefficients)
  statistic <- coefficients / se
  # A zero standard error is not an infinitely precise estimate; it is a fit
  # with no residual scatter at all --- an all-zero sink row on a transport
  # with full coverage is the ordinary way to get one. The ratio it produces
  # is not a statistic, and the zero-width interval it produces is not an
  # interval, so both are `NA` and the estimate itself is left alone: the
  # value is a measurement, and only its uncertainty is undefined.
  degenerate <- !is.finite(se) | se <= 0
  statistic[degenerate] <- NA_real_
  half <- stats::qt(1 - (1 - level) / 2, df = residual_df) * se
  half[degenerate] <- NA_real_
  list(
    estimate = coefficients,
    se = se,
    t = statistic,
    lower = coefficients - half,
    upper = coefficients + half,
    residual_sd = residual_sd,
    residual_df = as.integer(residual_df),
    level = as.numeric(level),
    calibration = "uncalibrated"
  )
}

.population_uncertainty_refuse_basis <- function() {
  .capability_refusal(paste0(
    "A between-subject standard error needs the participants' residuals ",
    "about the group fit, and the streamed complete-form route does not ",
    "retain them: its whole point is to bound the one array indexed by ",
    "participant, group node and packed coordinate. The coefficients it does ",
    "carry cannot reconstruct them."
  ),
    capability = "population_between_subject_residuals",
    namespace = "population_uncertainty",
    reasons = "streamed_route_retains_no_residuals",
    remedies = paste0(
      "Re-read the same population through `estimate_population(plan, ",
      "queries)` with the contrasts you want error bars for. That route ",
      "keeps `$residuals`, and its query bank is the axis the standard ",
      "errors are indexed by."
    ))
}

.population_uncertainty_refuse_df <- function(x) {
  .capability_refusal(sprintf(paste0(
    "The group model is saturated: %s and a design of rank %d leave %d ",
    "residual degrees of freedom, so the participants' scatter about the fit ",
    "is not estimable and no between-subject standard error exists. This is ",
    "a property of the estimand, not of the arithmetic."
  ), .msg_count(nrow(x$receipt$subjects), "participant"),
    nrow(x$receipt$subjects) - x$residual_df, x$residual_df),
    capability = "between_subject_residual_df",
    namespace = "population_uncertainty",
    reasons = "saturated_group_model",
    remedies = paste0(
      "Fit a smaller group model, or plan the population with more ",
      "participants. `residual_df` on the result is the number to watch."
    ))
}

.population_uncertainty_term <- function(x, term) {
  terms <- dimnames(x$coefficients)[[3L]]
  if (is.null(term)) return(terms)
  if (.is_number(term)) {
    position <- as.integer(term)
    if (position < 1L || position > length(terms)) {
      .input_error(sprintf(
        "`term` must select one of the %d group model columns; received %s.",
        length(terms), .msg_value(term)
      ), arg = "term", received = .msg_value(term),
        expected = paste0("one of ", .msg_names(terms)))
    }
    return(terms[[position]])
  }
  if (!.is_string(term) || !term %in% terms) {
    .input_error(sprintf(
      "`term` must name one group model column: %s.", .msg_names(terms)
    ), arg = "term", received = .msg_value(term),
      expected = paste0("one of ", .msg_names(terms)))
  }
  term
}

.population_uncertainty_id <- function(x, terms, level) {
  .sha256_signature(list(
    schema_version = 1L,
    contract_version = "population-form-v1",
    role = "population_uncertainty",
    parent = x$scientific_plan_id,
    terms = terms,
    level = as.numeric(level),
    layers = c("between_subject", "within_subject"),
    calibration = "uncalibrated"
  ), "population-sha256:")
}

.validate_population_uncertainty <- function(x) {
  expected <- c("basis", "layers", "separation", "between", "within",
    "index", "queries", "term", "ledger", "semantics", "normalization",
    "receipt", "scientific_plan_id")
  if (!.sealed_fields(x, "effect_population_uncertainty", expected) ||
      !is.list(x$between) || !is.data.frame(x$index) ||
      !.is_strings(x$term) || !.is_string(x$separation) ||
      !identical(x$between$calibration, "uncalibrated") ||
      !.strong_sha256(sub("^population-", "", x$scientific_plan_id))) {
    .input_error("Population-uncertainty fields are missing or noncanonical.")
  }
  if (!identical(dim(x$between$estimate), dim(x$between$se)) ||
      !identical(dim(x$between$estimate), dim(x$between$t))) {
    .contract_error(paste0(
      "A between-subject standard error must be indexed exactly as the ",
      "estimate it belongs to."
    ))
  }
  invisible(x)
}

#' Read the group uncertainty layers of an estimated population form
#'
#' `population_uncertainty()` reports two error bars on an
#' [estimate_population()] result and keeps them apart. The **between-subject**
#' layer is the scatter of the participants about the group fit: a standard
#' error per group node, query and model term from the group OLS residuals,
#' with that fit's residual degrees of freedom. The **within-subject** layer is
#' the sampling variance of one participant's transported value, and it exists
#' only where it is exact.
#'
#' @section The two layers are never pooled:
#' They answer different questions and are not summands. The between-subject
#' residual already contains whatever measurement error survived into each
#' participant's transported value, so adding the within-subject variance to it
#' would double-count the shared part and still miss the covariance a
#' variance-components model would need. The record carries `$between` and
#' `$within` as separate blocks, `as.data.frame()` emits one `layer` at a time,
#' and there is no field holding their sum.
#'
#' @section The `t` is uncalibrated, and stays uncalibrated:
#' `$between$t` is the estimate over its standard error. Its null distribution
#' has been *measured* against \eqn{t_{df}}{t_df} by
#' `benchmarks/run-population-null-coverage.R`, a 2,000-replication null
#' simulation of the group layer. Under a **correctly specified** group model
#' the nominal 95% interval covered the null term in **0.9485** of
#' replications at `N = 6`, **0.9500** at `N = 8`, **0.9520** at `N = 12` and
#' **0.9480** at `N = 24` (Monte Carlo standard error 0.005), and the
#' Kolmogorov-Smirnov distance between the null `t` and \eqn{t_{df}}{t_df} was
#' at most **0.0209** (`p >= 0.34`). The arithmetic is right, which was never
#' the part in doubt.
#'
#' The same simulation's second arm is why the label does not move. When each
#' participant's transported value carries a variance that depends on the group
#' covariates --- the transport-heterogeneity analogue, since a participant
#' whose warp is poor is noisier and warp quality is not independent of age,
#' motion or head size --- coverage of the same interval falls to **0.9230** at
#' `N = 6` and **0.8850** at `N = 24`. It gets **worse** with more
#' participants, because the bias is in the standard error and not in the
#' sample size: at `N = 24` the nominal 5% test rejects a true null **11.5%**
#' of the time, and the null `t` is distinguishable from \eqn{t_{df}}{t_df} at
#' `p = 1.2e-6`.
#'
#' A real population fit carries misspecification of unknown degree and a
#' transport whose displacement, entropy and subject coverage vary across
#' participants (`population-form-v1` section 7.5), which is exactly the
#' second arm's regime. Report the statistic; do not report a p-value derived
#' from it without an argument that those diagnostics are benign.
#'
#' @section What the within layer is admitted for:
#' Participant `i`'s transported value at group column `u` is a fixed linear
#' functional \eqn{w'z_i}{w'z_i} of that participant's native query values, so
#' its variance needs the covariance *between* native nodes. D8 refuses that
#' object (capability `cross_node_sampling_covariance`), so the layer is
#' admitted exactly where `w` has one nonzero entry: there the cross terms
#' carry weight zero and \eqn{\mathrm{Var} = w_x^2\,\mathrm{Var}(z_{ix})}{Var =
#' w_x^2 Var(z_ix)} is exact. **No independence assumption is made anywhere**,
#' and none would be defensible: overlapping searchlight supports under
#' spatially correlated noise are positively correlated, so a diagonal sum
#' would be an under-estimate in a known direction.
#'
#' Admission is per participant and per group column, and `$within$admitted`
#' is the matrix that records it. A hard anatomical parcellation usually
#' admits nothing, because every group node collects several native nodes;
#' that is the refusal being visible rather than the layer being broken.
#'
#' Two gates come before the per-column one, both recorded in
#' `$within$refusal$reasons` rather than raised:
#'
#' * `same_data_ratio_normalization` --- a `"unit_budget"` population divides
#'   each ledger by a total read from the same data, and section 4.3 records
#'   that the standard error of that divisor does not exist.
#' * `native_node_labels_unaligned` --- the covariance batch is named by the
#'   nodes `rdm_sampling_covariance(at = )` read and the transport by its own
#'   `native_index`. Binding two differently named node sets by position would
#'   attach one node's error bar to another's coefficient.
#'
#' @section Refusals:
#' Each is an `effect_capability_refusal` in namespace
#' `"population_uncertainty"` (see [catch_refusal()]).
#'
#' * `population_between_subject_residuals` --- a `"complete_form"` result.
#'   The streamed route does not retain the participant residuals, by design.
#' * `between_subject_residual_df` --- a saturated group model. With no
#'   residual degrees of freedom the participants' scatter is not estimable.
#'
#' @param x An `effect_population_result` from [estimate_population()].
#' @param term Which group model columns to report, by name or position.
#'   `NULL` (the default) reports every column.
#' @param level Nominal two-sided interval level for `$between$lower` and
#'   `$between$upper`. The default `0.95` is the level the recorded coverage
#'   simulation measured; the interval is labelled uncalibrated whatever level
#'   is asked for.
#' @return An `effect_population_uncertainty`: `$between` holding `$estimate`,
#'   `$se`, `$t`, `$lower`, `$upper` (each a `node`-by-`query`-by-`term`
#'   array), `$residual_sd`, `$residual_df`, `$level` and `$calibration`;
#'   `$within` holding the transported layer's `$admitted`, `$coefficient`,
#'   `$source_node`, `$variance` and `$refusal`, or `NULL` when
#'   [estimate_population()] was not given `uncertainty`; `$separation`
#'   stating that the two are never pooled; and the `$index`, `$queries`,
#'   `$ledger`, `$semantics`, `$normalization` and `$receipt` of the result it
#'   read. `as.data.frame(x, layer = )` returns one layer in long form.
#' @references `design/population-form-contract.md` (`population-form-v1`),
#'   section 7; `benchmarks/POPULATION-NULL-COVERAGE.md` for the recorded null
#'   simulation these numbers come from.
#' @family population transports
#' @seealso [estimate_population()] for the run this reads, and
#'   [population_views] for the group point estimates the standard errors
#'   belong to.
#' @examples
#' # Six participants on different native frames, one covariate at the group
#' # level, and a bank of two contrasts.
#' effects <- effect_space(c("face", "house"), basis_id = "popunc:v1")
#' subject <- function(id, n, gain) {
#'   domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
#'     feature_ids = paste0("f", seq_len(n)), id = id)
#'   values <- function(divisor) matrix(
#'     gain * seq_len(2 * n) / (n * divisor), 2, n,
#'     dimnames = list(c("face", "house"), NULL)
#'   )
#'   rel <- relation(list(run1 = values(1), run2 = values(1.7)),
#'     effects = effects, domain = domain)
#'   plan_geometry(rel, compile_frame(voxelwise(), domain),
#'     cross_partitions(rel))
#' }
#' carrier <- function(n) anatomical_transport(
#'   native_coords = cbind(seq_len(n) - 1), group_coords = cbind(c(0, 4)),
#'   semantics = "budget"
#' )
#' sizes <- c(s01 = 5L, s02 = 6L, s03 = 7L, s04 = 8L, s05 = 9L, s06 = 10L)
#' gains <- c(s01 = 1, s02 = 1.4, s03 = 0.7, s04 = 1.1, s05 = 0.9, s06 = 1.3)
#' subjects <- stats::setNames(lapply(names(sizes), function(id)
#'   subject(id, sizes[[id]], gains[[id]])), names(sizes))
#' covariates <- data.frame(age = c(21, 34, 27, 45, 31, 38))
#' rownames(covariates) <- names(sizes)
#' plan <- plan_population(subjects, lapply(sizes, carrier),
#'   model = ~ age, data = covariates)
#' fit <- estimate_population(plan, rbind(`face-house` = c(1, -1)))
#'
#' error_bars <- population_uncertainty(fit)
#' error_bars
#'
#' # The standard error of the group mean at each node, and its uncalibrated t.
#' error_bars$between$se[, , "(Intercept)"]
#' error_bars$between$t[, , "age"]
#'
#' # One layer at a time. There is no combined table, because the two layers
#' # are not summands.
#' head(as.data.frame(error_bars), 4)
#' @export
population_uncertainty <- function(x, term = NULL, level = 0.95) {
  .validate_population_result(x)
  if (!identical(x$basis, "query_bank")) {
    .population_uncertainty_refuse_basis()
  }
  level <- .check_number(level, "level",
    what = "one nominal two-sided interval level in (0, 1)")
  if (level <= 0 || level >= 1) {
    .input_error(
      "`level` must be one nominal two-sided interval level in (0, 1).",
      arg = "level", received = .msg_value(level),
      expected = "a number strictly between 0 and 1")
  }
  if (x$residual_df < 1L) {
    .population_uncertainty_refuse_df(x)
  }
  terms <- .population_uncertainty_term(x, term)
  between <- .population_between_statistics(
    x$coefficients[, , terms, drop = FALSE],
    x$residuals, x$residual_df, x$uncertainty$between$unscaled_covariance,
    level
  )
  between$estimator <- x$uncertainty$between$estimator
  between$calibration_evidence <- x$uncertainty$between$calibration_evidence
  value <- structure(list(
    basis = x$basis,
    layers = c("between_subject", "within_subject"),
    separation = x$uncertainty$separation,
    between = between,
    # The within-subject variance has no term axis at all --- it belongs to
    # one participant's transported value, not to a coefficient of the group
    # model --- which is the structural reason the two layers cannot be joined
    # into one table, and why `as.data.frame()` asks which one you want.
    within = x$uncertainty$within,
    index = x$index,
    queries = x$queries,
    term = terms,
    ledger = x$ledger,
    semantics = x$semantics,
    normalization = x$normalization,
    receipt = x$receipt,
    scientific_plan_id = .population_uncertainty_id(x, terms, level)
  ), class = "effect_population_uncertainty")
  .validate_population_uncertainty(value)
  value
}
