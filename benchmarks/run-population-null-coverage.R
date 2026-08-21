#!/usr/bin/env Rscript

# Null coverage of the between-subject uncertainty layer (ticket E8).
#
# WHAT THIS SIMULATION EXERCISES, AND WHAT IT DOES NOT
#
# It exercises the *group* layer of `estimate_population()`: the pivoted-QR
# OLS over the subject axis (`.population_ols()`), the unscaled covariance
# recorded at fit time (`.population_between_uncertainty()`), and the standard
# error, t statistic and interval `population_uncertainty()` reports
# (`.population_between_statistics()`). Those are the shipped functions, called
# here directly rather than re-implemented, and `wiring` below asserts that the
# fast path reproduces `population_uncertainty()` on a real fitted population
# bit for bit.
#
# It does **not** exercise the subject geometry kernels. The per-subject
# transported query values `z_{i,u,k}` are simulated directly, because the
# object under test is the sampling behaviour of a standard error computed from
# `N` numbers per group node and query, and running a geometry plan 2,000 times
# would measure the compiler rather than the estimator. The transport, the
# query and the fit commute (`population-form-v1` section 3), so the group
# layer sees exactly the array this script hands it.
#
# THE TWO ARMS
#
# `gaussian`        the group model is correct: homoskedastic normal subject
#                   effects, spatially correlated across group nodes. This is
#                   the case where nominal coverage is expected, and the number
#                   is recorded so the claim is a measurement and not a hope.
#
# `heteroskedastic` the group model is misspecified in the way a real
#                   population study is misspecified: each participant's
#                   transported value carries a variance that depends on a
#                   group covariate. Participants whose transport is poor are
#                   noisier, and transport quality is not independent of age,
#                   motion or head size. Coverage is recorded here too, and it
#                   is the reason the shipped label stays `uncalibrated`
#                   whatever the `gaussian` arm says.
#
# WHAT IS REPORTED
#
# Every distributional summary is computed from **one designated cell**
# (`group1`, the first query) per replication, so the `M` draws are
# independent, the Monte Carlo standard error is exact, and the KS statistic is
# valid. Coverage pooled over every node-query cell is reported beside it and
# is labelled dependent, because cells within a replication share subjects.

suppressPackageStartupMessages(devtools::load_all(quiet = TRUE))

arguments <- commandArgs(trailingOnly = TRUE)
replications <- as.integer(Sys.getenv("CROSSFORM_NULL_COVERAGE_REPS", "2000"))
if (length(arguments) > 0L && nzchar(arguments[[1L]])) {
  replications <- as.integer(arguments[[1L]])
}
if (is.na(replications) || replications < 500L) {
  stop("Null coverage requires at least 500 replications.", call. = FALSE)
}
results_dir <- if (length(arguments) >= 2L && nzchar(arguments[[2L]])) {
  arguments[[2L]]
} else {
  "benchmark-results"
}

repo <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(repo, "benchmarks", "provenance.R"), local = TRUE)
provenance <- crossform_benchmark_provenance(
  repo, "run-population-null-coverage.R"
)

level <- 0.95
subject_counts <- c(6L, 8L, 12L, 24L)
group_nodes <- c("group1", "group2", "group3", "<sink>")
query_labels <- c("face-house", "face-tool")
signal_beta <- 0.6
between_sd <- 0.8
within_sd <- 0.6
node_correlation <- 0.5

# The shipped group layer, called the way `estimate_population()` calls it.
# The reshape is the one the driver performs: `.population_ols()` returns
# `p`-by-`(m+1)K` and the result carries `node x query x term`.
group_layer <- function(model, stack, nodes, queries, subjects,
                        availability = is.finite(stack)) {
  fit <- crossform:::.population_ols(
    model, stack, availability, coverage_policy = "available"
  )
  shape <- function(values, third, name) {
    array(t(values), c(length(nodes), length(queries), length(third)),
      dimnames = stats::setNames(list(nodes, queries, third),
        c("node", "query", name)))
  }
  coverage_labels <- list(node = nodes, query = queries)
  result <- list(
    coefficients = shape(fit$coefficients, model$columns, "term"),
    values = shape(stack, subjects, "subject"),
    residuals = shape(fit$residuals, subjects, "subject"),
    uncertainty = list(between = crossform:::.population_between_uncertainty(
      model
    )),
    coverage = list(
      availability = shape(availability, subjects, "subject"),
      n = array(fit$n, c(length(nodes), length(queries)),
                dimnames = coverage_labels),
      design_rank = array(fit$rank, c(length(nodes), length(queries)),
                          dimnames = coverage_labels),
      residual_df = array(fit$residual_df,
                          c(length(nodes), length(queries)),
                          dimnames = coverage_labels),
      status = array(fit$status, c(length(nodes), length(queries)),
                     dimnames = coverage_labels),
      coefficient_estimable = shape(
        fit$estimable, model$columns, "term"
      )
    )
  )
  crossform:::.population_between_statistics(
    result, terms = model$columns, level = level, estimator = "classical",
    leverage_tolerance = sqrt(.Machine$double.eps)
  )
}

# One replication: draw the covariates, draw the transported values under the
# declared truth, and return the layer.
replicate_once <- function(n_subject, arm) {
  subjects <- sprintf("s%02d", seq_len(n_subject))
  covariates <- data.frame(
    null_x = stats::rnorm(n_subject),
    signal_x = stats::rnorm(n_subject),
    row.names = subjects
  )
  model <- crossform:::.population_model(
    ~ null_x + signal_x, covariates, subjects
  )
  nodes <- length(group_nodes)
  queries <- length(query_labels)
  # Spatially correlated subject error across group nodes: an AR(1) whitening
  # of independent draws. The marginal variance is unchanged, so the truth the
  # interval is checked against does not move; what changes is that the cells
  # of one replication are no longer independent, which is the realistic case.
  chol_factor <- chol(node_correlation^abs(outer(
    seq_len(nodes), seq_len(nodes), "-"
  )))
  draw <- function(sd) {
    noise <- array(stats::rnorm(nodes * queries * n_subject),
      c(nodes, queries, n_subject))
    for (subject in seq_len(n_subject)) {
      noise[, , subject] <- crossprod(chol_factor, noise[, , subject]) *
        sd[[subject]]
    }
    noise
  }
  scale <- if (identical(arm, "heteroskedastic")) {
    # Variance linked to *both* group covariates: the transport-heterogeneity
    # analogue. A participant whose warp is poor carries a noisier transported
    # value, and warp quality is not independent of age, motion or head size,
    # so the noise scale is a function of the same columns the group model
    # regresses on. Rescaled to keep the average variance equal to the
    # gaussian arm's, so the two arms differ in shape and not in magnitude.
    raw <- exp(0.7 * (covariates$null_x + covariates$signal_x) / sqrt(2))
    raw / sqrt(mean(raw^2))
  } else {
    rep(1, n_subject)
  }
  errors <- draw(between_sd * scale) + draw(within_sd * scale)
  # Truth: an intercept that varies over nodes and queries, no effect of
  # `null_x` anywhere, and a fixed effect of `signal_x` for power context.
  intercept <- outer(seq_len(nodes) / nodes, seq_len(queries))
  values <- array(0, c(nodes, queries, n_subject))
  for (subject in seq_len(n_subject)) {
    values[, , subject] <- intercept +
      signal_beta * covariates$signal_x[[subject]] + errors[, , subject]
  }
  stack <- matrix(aperm(values, c(3L, 1L, 2L)), n_subject,
    nodes * queries, dimnames = list(subjects, NULL))
  group_layer(model, stack, group_nodes, query_labels, subjects)
}

covers <- function(layer, term, truth) {
  layer$lower[, , term] <= truth & truth <= layer$upper[, , term]
}

run_cell <- function(n_subject, arm, seed) {
  set.seed(seed)
  residual_df <- n_subject - 3L
  null_t <- numeric(replications)
  null_cover <- logical(replications)
  signal_cover <- logical(replications)
  signal_reject <- logical(replications)
  pooled_null <- numeric(replications)
  pooled_signal <- numeric(replications)
  critical <- stats::qt(1 - (1 - level) / 2, df = residual_df)
  for (replication in seq_len(replications)) {
    layer <- replicate_once(n_subject, arm)
    null_t[[replication]] <- layer$t[1L, 1L, "null_x"]
    null_cover[[replication]] <- covers(layer, "null_x", 0)[1L, 1L]
    signal_cover[[replication]] <-
      covers(layer, "signal_x", signal_beta)[1L, 1L]
    signal_reject[[replication]] <-
      abs(layer$t[1L, 1L, "signal_x"]) > critical
    pooled_null[[replication]] <- mean(covers(layer, "null_x", 0))
    pooled_signal[[replication]] <-
      mean(covers(layer, "signal_x", signal_beta))
  }
  ks <- suppressWarnings(stats::ks.test(null_t, "pt", df = residual_df))
  coverage <- mean(null_cover)
  data.frame(
    arm = arm,
    subjects = n_subject,
    residual_df = residual_df,
    replications = replications,
    null_coverage = coverage,
    null_coverage_mcse = sqrt(coverage * (1 - coverage) / replications),
    null_rejection = 1 - coverage,
    null_coverage_all_cells = mean(pooled_null),
    signal_coverage = mean(signal_cover),
    signal_coverage_all_cells = mean(pooled_signal),
    signal_power = mean(signal_reject),
    t_mean = mean(null_t),
    t_sd = stats::sd(null_t),
    t_reference_sd = if (residual_df > 2L) {
      sqrt(residual_df / (residual_df - 2))
    } else {
      NA_real_
    },
    ks_statistic = unname(ks$statistic),
    ks_p_value = unname(ks$p.value),
    stringsAsFactors = FALSE
  )
}

# The wiring check. A real population is fitted through the public verbs, the
# fast path is run on the transported values that fit produced, and the two
# standard errors must agree exactly. Without this the coverage number would be
# evidence about this script rather than about `population_uncertainty()`.
wiring_check <- function() {
  effects <- effect_space(c("face", "house", "tool"), basis_id = "e8:null")
  subject <- function(id, n, gain) {
    domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
      feature_ids = paste0("f", seq_len(n)), id = id)
    values <- function(divisor) matrix(
      gain * seq_len(3 * n) / (n * divisor), 3, n,
      dimnames = list(c("face", "house", "tool"), NULL)
    )
    rel <- relation(list(run1 = values(1), run2 = values(1.6)),
      effects = effects, domain = domain)
    plan_geometry(rel, compile_frame(voxelwise(), domain),
      cross_partitions(rel))
  }
  sizes <- c(s01 = 5L, s02 = 6L, s03 = 7L, s04 = 8L, s05 = 9L, s06 = 10L)
  gains <- c(s01 = 1, s02 = 1.4, s03 = 0.7, s04 = 1.1, s05 = 0.9, s06 = 1.3)
  carrier <- function(n) anatomical_transport(
    native_coords = cbind(seq_len(n) - 1),
    group_coords = cbind(c(0, 4, 9)), semantics = "budget"
  )
  covariates <- data.frame(
    null_x = c(0.4, -1.1, 0.7, -0.3, 1.2, -0.9),
    signal_x = c(-0.8, 0.5, 1.3, -1.4, 0.2, 0.9),
    row.names = names(sizes)
  )
  plan <- plan_population(
    stats::setNames(lapply(names(sizes), function(id)
      subject(id, sizes[[id]], gains[[id]])), names(sizes)),
    lapply(sizes, carrier),
    model = ~ null_x + signal_x, data = covariates
  )
  fit <- estimate_population(plan,
    rbind(`face-house` = c(1, -1, 0), `face-tool` = c(1, 0, -1)))
  verb <- population_uncertainty(fit, level = level)
  stack <- matrix(aperm(fit$values, c(3L, 1L, 2L)), nrow(fit$receipt$subjects),
    dim(fit$values)[[1L]] * dim(fit$values)[[2L]],
    dimnames = list(dimnames(fit$values)[[3L]], NULL))
  availability <- matrix(
    aperm(fit$coverage$availability, c(3L, 1L, 2L)),
    nrow(fit$receipt$subjects),
    dim(fit$values)[[1L]] * dim(fit$values)[[2L]],
    dimnames = list(dimnames(fit$values)[[3L]], NULL)
  )
  fast <- group_layer(plan$model, stack, dimnames(fit$values)[[1L]],
    dimnames(fit$values)[[2L]], dimnames(fit$values)[[3L]], availability)
  se_difference <- abs(fast$se - verb$between$se)
  t_difference <- abs(fast$t - verb$between$t)
  fast_se_missing <- is.na(fast$se)
  verb_se_missing <- is.na(verb$between$se)
  list(
    se_max_absolute_difference = max(se_difference, na.rm = TRUE),
    t_max_absolute_difference = max(t_difference, na.rm = TRUE),
    residual_df = verb$between$residual_df,
    identical_se = identical(
      as.vector(fast_se_missing), as.vector(verb_se_missing)
    ) &&
      all(se_difference[is.finite(se_difference)] == 0),
    fast_missing_se = which(fast_se_missing, arr.ind = TRUE),
    verb_missing_se = which(verb_se_missing, arr.ind = TRUE),
    # An independent court for one cell: base R's own `lm()` on the same
    # response, which shares no code with the population driver.
    lm_max_absolute_difference = local({
      response <- fit$values["group1", "face-house", ]
      reference <- summary(stats::lm(
        response ~ covariates$null_x + covariates$signal_x
      ))$coefficients[, "Std. Error"]
      max(abs(unname(reference) - unname(verb$between$se["group1",
        "face-house", ])))
    })
  )
}

wiring <- wiring_check()
if (!isTRUE(wiring$identical_se) ||
    !isTRUE(wiring$se_max_absolute_difference < 1e-12) ||
    !isTRUE(wiring$lm_max_absolute_difference < 1e-10)) {
  print(wiring)
  stop("The simulated group layer does not reproduce population_uncertainty().",
    call. = FALSE)
}

cells <- list()
seed <- 81400L
for (arm in c("gaussian", "heteroskedastic")) {
  for (n_subject in subject_counts) {
    seed <- seed + 1L
    cells[[length(cells) + 1L]] <- run_cell(n_subject, arm, seed)
  }
}
summary <- do.call(rbind, cells)

record <- list(
  schema_version = 1L,
  provenance = provenance,
  role = "between_subject_uncertainty_null_coverage",
  verb = "population_uncertainty",
  exercises = paste0(
    "group layer only: .population_ols(), .population_between_uncertainty(), ",
    ".population_between_statistics(). Subject geometry kernels are not run; ",
    "transported query values are simulated directly."
  ),
  replications = replications,
  level = level,
  subject_counts = subject_counts,
  group_nodes = group_nodes,
  queries = query_labels,
  truth = list(
    null_term = "null_x",
    null_coefficient = 0,
    signal_term = "signal_x",
    signal_coefficient = signal_beta,
    between_subject_sd = between_sd,
    within_subject_sd = within_sd,
    node_correlation = node_correlation,
    heteroskedastic_link =
      "sd_i proportional to exp(0.7 * (null_x_i + signal_x_i) / sqrt(2))"
  ),
  wiring = wiring,
  summary = summary
)

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(record, file.path(results_dir, "population-null-coverage.rds"),
  version = 2L, compress = "xz")
utils::write.csv(summary,
  file.path(results_dir, "population-null-coverage-summary.csv"),
  row.names = FALSE)

print(summary, row.names = FALSE, digits = 4)
cat(sprintf(
  paste0("\nwiring: max |SE(fast) - SE(population_uncertainty)| = %.3g\n",
    "wiring: max |SE(population_uncertainty) - SE(lm)| = %.3g\n"),
  wiring$se_max_absolute_difference, wiring$lm_max_absolute_difference
))
