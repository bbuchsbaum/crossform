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
# Everything here is a cellwise operation over the `node x query x term`
# array. A non-finite response column --- an unresolved density node or an
# unadmitted `unit_budget` share --- carries `NA` through to its own standard
# error and to nothing else. The independent HC3 court deliberately writes out
# its own sandwich arithmetic; the classical coverage benchmark exercises the
# public reader and is maintained separately.
.population_bread <- function(factorization, columns) {
  rank <- factorization$rank
  triangle <- qr.R(factorization)[seq_len(rank), seq_len(rank), drop = FALSE]
  pivot <- factorization$pivot[seq_len(rank)]
  value <- matrix(NA_real_, columns, columns)
  value[pivot, pivot] <- chol2inv(triangle)
  value
}

# One covariance calculation per node-query cell. D1 made the contributing
# subject set a part of the result, so neither the design, its leverage nor its
# residual degrees of freedom may be borrowed from another cell. HC3 uses
#
#   (X'X)^-1 X' diag(e_i^2 / (1 - h_i)^2) X (X'X)^-1,
#
# and classical OLS uses s^2 (X'X)^-1. Both are returned as covariance
# matrices, not just diagonals, because coefficient recombination needs the
# off-diagonal entries just as much as a single coefficient needs its SE.
.population_between_statistics <- function(x, terms, level, estimator,
                                           leverage_tolerance) {
  all_terms <- dimnames(x$coefficients)[[3L]]
  term_position <- match(terms, all_terms)
  nodes <- dim(x$coefficients)[[1L]]
  queries <- dim(x$coefficients)[[2L]]
  subjects <- dim(x$values)[[3L]]
  design <- x$uncertainty$between$design
  labels <- list(
    node = dimnames(x$coefficients)[[1L]],
    query = dimnames(x$coefficients)[[2L]],
    term = terms
  )
  estimate <- x$coefficients[, , terms, drop = FALSE]
  se <- t_statistic <- lower <- upper <- array(NA_real_,
    c(nodes, queries, length(terms)), dimnames = labels)
  covariance <- array(NA_real_,
    c(nodes, queries, length(terms), length(terms)),
    dimnames = c(labels[1:2], list(term = terms, term2 = terms)))
  leverage <- array(NA_real_, c(nodes, queries, subjects),
    dimnames = list(
      node = labels$node, query = labels$query,
      subject = dimnames(x$values)[[3L]]
    ))
  adjusted_residual <- leverage
  residual_sd <- max_leverage <- matrix(NA_real_, nodes, queries,
    dimnames = labels[1:2])
  status <- matrix("refused", nodes, queries, dimnames = labels[1:2])
  reason <- matrix(NA_character_, nodes, queries, dimnames = labels[1:2])
  residual_df <- x$coverage$residual_df
  design_rank <- x$coverage$design_rank
  n <- x$coverage$n

  for (node in seq_len(nodes)) {
    for (query in seq_len(queries)) {
      available <- x$coverage$availability[node, query, ]
      if (!all(x$coverage$coefficient_estimable[node, query, ])) {
        reason[node, query] <- if (identical(
            x$coverage$status[node, query], "rank_deficient")) {
          "rank_deficient_design"
        } else {
          x$coverage$status[node, query]
        }
        next
      }
      if (residual_df[node, query] < 1L) {
        reason[node, query] <- "saturated_group_model"
        next
      }
      cell_design <- design[available, , drop = FALSE]
      cell_residual <- x$residuals[node, query, available]
      factorization <- qr(cell_design)
      if (factorization$rank < ncol(cell_design)) {
        reason[node, query] <- "rank_deficient_design"
        next
      }
      bread <- .population_bread(factorization, ncol(cell_design))
      cell_leverage <- rowSums((cell_design %*% bread) * cell_design)
      # Clamp only roundoff outside the mathematical [0, 1] interval. A value
      # genuinely close to one is retained and causes HC3 to refuse below.
      cell_leverage[cell_leverage < 0 &
        cell_leverage > -leverage_tolerance] <- 0
      cell_leverage[cell_leverage > 1 &
        cell_leverage < 1 + leverage_tolerance] <- 1
      leverage[node, query, available] <- cell_leverage
      max_leverage[node, query] <- max(cell_leverage)
      residual_sd[node, query] <- sqrt(
        sum(cell_residual^2) / residual_df[node, query]
      )
      if (identical(estimator, "HC3") &&
          any(1 - cell_leverage <= leverage_tolerance)) {
        reason[node, query] <- "leverage_near_one"
        next
      }
      adjusted_residual[node, query, available] <- if (
          identical(estimator, "HC3")) {
        cell_residual / (1 - cell_leverage)
      } else {
        cell_residual
      }
      cell_covariance <- if (identical(estimator, "HC3")) {
        adjusted_square <- adjusted_residual[node, query, available]^2
        meat <- crossprod(cell_design, cell_design * adjusted_square)
        bread %*% meat %*% bread
      } else {
        residual_sd[node, query]^2 * bread
      }
      cell_covariance <- (cell_covariance + t(cell_covariance)) / 2
      selected <- cell_covariance[term_position, term_position, drop = FALSE]
      cell_se <- sqrt(pmax(diag(selected), 0))
      covariance[node, query, , ] <- selected
      se[node, query, ] <- cell_se
      residual_scale_zero <- residual_sd[node, query] <=
        sqrt(.Machine$double.eps) * max(
          1, abs(x$values[node, query, available])
        )
      valid <- is.finite(cell_se) & cell_se > 0 & !residual_scale_zero
      if (any(valid)) {
        t_statistic[node, query, valid] <-
          estimate[node, query, valid] / cell_se[valid]
        half <- stats::qt(1 - (1 - level) / 2,
          df = residual_df[node, query]) * cell_se[valid]
        lower[node, query, valid] <- estimate[node, query, valid] - half
        upper[node, query, valid] <- estimate[node, query, valid] + half
      }
      if (all(valid)) {
        status[node, query] <- "estimated"
      } else {
        status[node, query] <- "degenerate_residual_scale"
        reason[node, query] <- "nonpositive_standard_error"
      }
    }
  }
  assumptions <- c(
    x$uncertainty$between$assumptions$common,
    x$uncertainty$between$assumptions[[estimator]]
  )
  list(
    estimator = if (identical(estimator, "HC3")) "HC3" else "classical_ols",
    assumptions = assumptions,
    estimate = estimate,
    covariance = covariance,
    se = se,
    t = t_statistic,
    lower = lower,
    upper = upper,
    residual_sd = residual_sd,
    n = n,
    design_rank = design_rank,
    residual_df = residual_df,
    leverage = leverage,
    adjusted_residual = adjusted_residual,
    max_leverage = max_leverage,
    leverage_tolerance = as.numeric(leverage_tolerance),
    status = status,
    reason = reason,
    level = as.numeric(level),
    calibration = "uncalibrated",
    conditioning = x$coverage$conditioning
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
  ranks <- x$coverage$design_rank
  dfs <- x$coverage$residual_df
  .capability_refusal(sprintf(paste0(
    "No population cell has positive residual degrees of freedom: the local ",
    "design ranks range from %d to %d and the local residual df from %d to ",
    "%d. The participants' scatter about a saturated or unidentified fit is ",
    "not estimable, so no between-subject covariance exists. This is a ",
    "property of the cellwise estimand, not of the arithmetic."
  ), min(ranks), max(ranks), min(dfs), max(dfs)),
    capability = "between_subject_residual_df",
    namespace = "population_uncertainty",
    reasons = "saturated_group_model",
    remedies = paste0(
      "Fit a smaller group model, or plan the population with more ",
      "participants. `$coverage$residual_df` on the result names the local ",
      "degrees of freedom to watch."
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

.population_uncertainty_id <- function(x, terms, level, estimator,
                                       leverage_tolerance) {
  .sha256_signature(list(
    schema_version = 1L,
    contract_version = "population-form-v1",
    role = "population_uncertainty",
    parent = x$scientific_plan_id,
    terms = terms,
    level = as.numeric(level),
    estimator = estimator,
    leverage_tolerance = if (identical(estimator, "HC3")) {
      as.numeric(leverage_tolerance)
    } else {
      NULL
    },
    layers = c("between_subject", "within_subject"),
    calibration = "uncalibrated"
  ), "population-sha256:")
}

.validate_population_uncertainty <- function(x) {
  expected <- c("basis", "layers", "separation", "estimator", "assumptions",
    "between", "within",
    "index", "queries", "term", "ledger", "semantics", "normalization",
    "receipt", "scientific_plan_id")
  if (!.sealed_fields(x, "effect_population_uncertainty", expected) ||
      !is.list(x$between) || !is.data.frame(x$index) ||
      !.is_strings(x$term) || !.is_string(x$separation) ||
      !.is_string(x$estimator) ||
      !x$estimator %in% c("classical_ols", "HC3") ||
      !.is_strings(x$assumptions, unique = TRUE) ||
      !identical(x$estimator, x$between$estimator) ||
      !identical(x$assumptions, x$between$assumptions) ||
      !identical(x$between$calibration, "uncalibrated") ||
      !.strong_sha256(sub("^population-", "", x$scientific_plan_id))) {
    .input_error("Population-uncertainty fields are missing or noncanonical.")
  }
  .validate_population_conditioning(
    x$between$conditioning, as.character(x$receipt$subjects$subject)
  )
  if (!identical(dim(x$between$estimate), dim(x$between$se)) ||
      !identical(dim(x$between$estimate), dim(x$between$t)) ||
      !identical(dim(x$between$estimate), dim(x$between$lower)) ||
      !identical(dim(x$between$estimate), dim(x$between$upper)) ||
      !identical(dim(x$between$covariance), c(
        dim(x$between$estimate)[1:2], rep(dim(x$between$estimate)[[3L]], 2L)
      )) ||
      !identical(dim(x$between$residual_df),
        dim(x$between$estimate)[1:2]) ||
      !identical(dim(x$between$status), dim(x$between$residual_df)) ||
      !identical(dim(x$between$reason), dim(x$between$residual_df)) ||
      !identical(dim(x$between$leverage), c(
        dim(x$between$estimate)[1:2], nrow(x$receipt$subjects)
      )) ||
      !identical(dim(x$between$adjusted_residual),
        dim(x$between$leverage))) {
    .contract_error(paste0(
      "A between-subject standard error must be indexed exactly as the ",
      "estimate it belongs to."
    ))
  }
  admitted <- x$between$status %in% c("estimated", "degenerate_residual_scale")
  refused_terms <- array(
    x$between$status == "refused", dim(x$between$se)
  )
  if (any(admitted & x$between$residual_df < 1L) ||
      any(x$between$status == "estimated" &
        !apply(is.finite(x$between$se), c(1L, 2L), all)) ||
      any(refused_terms & is.finite(x$between$se))) {
    .contract_error("Population uncertainty disagrees with its cell status.")
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
#' @section Classical OLS and HC3:
#' `estimator = "classical"` reproduces the ordinary OLS covariance
#' \eqn{s^2(X'X)^{-1}}{s^2 (X'X)^-1}. `estimator = "HC3"` returns the
#' heteroskedasticity-robust sandwich covariance
#' \eqn{(X'X)^{-1}X'\mathrm{diag}\{e_i^2/(1-h_i)^2\}X(X'X)^{-1}}{(X'X)^-1
#' X' diag(e_i^2 / (1-h_i)^2) X (X'X)^-1}. The latter is the recommended
#' sensitivity analysis when subject-level variance differs with covariates or
#' transport quality. It is not a distributional calibration theorem, so both
#' routes retain the `"uncalibrated"` label.
#'
#' Every node-query cell is computed from its exact contributing subject set.
#' `$between` therefore carries the full coefficient covariance, SEs, local
#' `n`, design rank and residual df, subject leverage, the residual adjustment
#' actually used in the sandwich, maximum leverage, conditioning statement,
#' and a cellwise status/reason. Rank-deficient and
#' saturated cells are refused. HC3 also refuses a cell when any
#' \eqn{1-h_i}{1-h_i} is no larger than `leverage_tolerance`; dividing by an
#' almost-zero leverage complement would turn a numerical singularity into a
#' confident-looking result.
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
#' for `estimator = "classical"` has been *measured* against
#' \eqn{t_{df}}{t_df} by
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
#' Those coverage figures do not transfer to HC3. Its randomized court proves
#' agreement with the defining sandwich covariance over full-rank designs and
#' adversarial leverage fixtures; it does not supply a finite-sample reference
#' distribution. HC3 therefore retains the same `"uncalibrated"` label.
#' The expanded paired benchmark in
#' `inst/extdata/certification/population-calibration-results.csv` records
#' classical, HC3, and wild-bootstrap behavior across eight synthetic regimes.
#' It supports only the regimes and conditional targets named there:
#' informative coverage remains ineligible for a marginal claim, and no arm is
#' marginal over transport estimation.
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
#' A future joint covariance may support transported precision or a separately
#' calibrated spatial procedure, but it is not required by the cellwise
#' between-subject covariance computed here. Dense/sparse representation, PSD
#' validation and scaling requirements are specified in
#' `design/cross-node-covariance-contract.md`. maxT and multiple-comparison
#' methods are optional later consumers, not part of this gate.
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
#' @section Transport conditioning:
#' Classical and HC3 intervals use the exact subject-level residuals and may
#' account for homoskedastic or heteroskedastic between-subject variation.
#' They remain conditional on each realized transport. In particular, HC3
#' does not include uncertainty from learning a functional transport or from
#' assigning its cross-fitting folds. The complete conditioning record is
#' retained in `$between$conditioning`.
#'
#' @param x An `effect_population_result` from [estimate_population()].
#' @param term Which group model columns to report, by name or position.
#'   `NULL` (the default) reports every column.
#' @param level Nominal two-sided interval level for `$between$lower` and
#'   `$between$upper`. The default `0.95` is the level the recorded coverage
#'   simulation measured; the interval is labelled uncalibrated whatever level
#'   is asked for.
#' @param estimator Between-subject covariance estimator. `"classical"` is
#'   ordinary homoskedastic OLS; `"HC3"` is the leverage-adjusted
#'   heteroskedasticity-robust sandwich estimator.
#' @param leverage_tolerance Positive tolerance below which HC3 treats
#'   \eqn{1-h_i}{1-h_i} as numerically zero and refuses that cell.
#' @return An `effect_population_uncertainty`: `$between` holding `$estimate`,
#'   `$covariance`, `$se`, `$t`, `$lower`, `$upper`, leverage, adjusted
#'   residuals and local design diagnostics, estimator identity, assumptions,
#'   status and refusal reasons;
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
population_uncertainty <- function(x, term = NULL, level = 0.95,
                                   estimator = c("classical", "HC3"),
                                   leverage_tolerance = 1e-8) {
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
  estimator <- match.arg(estimator)
  leverage_tolerance <- .check_number(
    leverage_tolerance, "leverage_tolerance",
    what = "one positive finite leverage-complement tolerance"
  )
  if (leverage_tolerance <= 0 || leverage_tolerance >= 1) {
    .input_error(
      "`leverage_tolerance` must lie strictly between zero and one.",
      arg = "leverage_tolerance", received = .msg_value(leverage_tolerance),
      expected = "one number strictly between zero and one"
    )
  }
  if (!any(x$coverage$residual_df > 0L &
      apply(x$coverage$coefficient_estimable, c(1L, 2L), all))) {
    .population_uncertainty_refuse_df(x)
  }
  terms <- .population_uncertainty_term(x, term)
  between <- .population_between_statistics(
    x, terms, level, estimator, leverage_tolerance
  )
  between$calibration_evidence <- x$uncertainty$between$calibration_evidence
  value <- structure(list(
    basis = x$basis,
    layers = c("between_subject", "within_subject"),
    separation = x$uncertainty$separation,
    estimator = between$estimator,
    assumptions = between$assumptions,
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
    scientific_plan_id = .population_uncertainty_id(
      x, terms, level, estimator, leverage_tolerance
    )
  ), class = "effect_population_uncertainty")
  .validate_population_uncertainty(value)
  value
}
