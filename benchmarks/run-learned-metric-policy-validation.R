#!/usr/bin/env Rscript

# Statistical validation of learned local precision policies.
#
# Run from the package root:
#
#   Rscript benchmarks/run-learned-metric-policy-validation.R 500 benchmark-results
#
# The script deliberately calls the exported product path for every estimate.
# It separates evaluator error relative to the realized learned metric from
# learner plus evaluator error relative to the population Mahalanobis target.
#
# Two generative regimes are used. well_specified is the Gaussian iid GLM
# under which fitted effects and residuals are independent. unmodelled_ar1
# generates stationary AR(1) observation noise but intentionally fits the
# identity observation metric. The latter is a scoped misspecification stress
# test, not a model of every possible analysis failure.

if (!requireNamespace("crossform", quietly = TRUE)) {
  stop("Install crossform before running policy validation.")
}

library(crossform)
repo <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(repo, "benchmarks", "provenance.R"), local = TRUE)
provenance <- crossform_benchmark_provenance(
  repo, "run-learned-metric-policy-validation.R"
)

arguments <- commandArgs(trailingOnly = TRUE)
replications <- if (length(arguments) >= 1L) {
  suppressWarnings(as.integer(arguments[[1L]]))
} else {
  500L
}
output_dir <- if (length(arguments) >= 2L) arguments[[2L]] else {
  "benchmark-results"
}
if (length(replications) != 1L || is.na(replications) ||
    replications < 8L) {
  stop("The replication count must be one integer of at least eight.")
}
if (!is.character(output_dir) || length(output_dir) != 1L ||
    is.na(output_dir) || !nzchar(output_dir)) {
  stop("The output directory must be one nonempty path.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 2026081307L
equivalence_margin <- 0.005
confidence_level <- 0.95
observations <- 48L
features <- 6L
partitions <- paste0("run", seq_len(4L))
selected_node <- 3L
shrinkage <- 0.2
ar1_phi <- 0.6

domain <- abstract_domain(
  features,
  coordinates = cbind(x = seq_len(features), y = 0, z = 0),
  feature_ids = paste0("feature", seq_len(features)),
  id = "learned-metric-policy-validation:v1",
  coordinate_units = "mm"
)
frame <- compile_frame(
  searchlights(radius = 1.01, normalization = "local"),
  domain
)
over <- cross_partitions(partitions, independence = "independent")
condition <- rep(c(-0.5, 0.5), each = observations / 2L)
null_regressor <- rep(c(-0.5, 0.5), times = observations / 2L)
design <- cbind(
  intercept = 1,
  null = null_regressor,
  signal = condition
)
effects <- rbind(
  baseline = c(1, 0, 0),
  null = c(0, 1, 0),
  signal = c(0, 0, 1)
)
contrasts <- list(
  null = c(baseline = 0, null = 1, signal = 0),
  signal = c(baseline = 0, null = 0, signal = 1)
)
neural_covariance <- toeplitz(0.55^(0:(features - 1L)))
neural_factor <- chol(neural_covariance)
signal_pattern <- c(0, 0.55, -0.35, 0.25, 0.45, 0)
recipe <- shrinkage_precision(
  shrinkage = shrinkage,
  relative_variance_floor = 1e-8,
  relative_spectral_floor = 1e-10
)
policies <- list(
  training_only = metric_training_policy("exclude_evaluation"),
  all_run_residuals = metric_training_policy(
    "all_partitions_residual_orthogonality",
    justification = paste(
      "Reuse is admitted only under the declared fitted-effect and GLM",
      "residual-orthogonality assumption; this validation includes a",
      "deliberately misspecified temporal-noise arm."
    )
  )
)

support <- crossform:::.support_index_support(
  frame$support_index, selected_node
)[[1L]]
frame_weights <- as.numeric(frame$weights[selected_node, support])
root_weights <- sqrt(frame_weights)
population_precision <- solve(
  neural_covariance[support, support, drop = FALSE]
)
weighted_signal <- signal_pattern[support] * root_weights
population_targets <- c(
  null = 0,
  signal = drop(
    weighted_signal %*% population_precision %*% weighted_signal
  )
)

generate_noise <- function(regime) {
  innovations <- matrix(
    stats::rnorm(observations * features), observations, features
  ) %*% neural_factor
  if (identical(regime, "well_specified")) return(innovations)
  value <- matrix(0, observations, features)
  value[1L, ] <- innovations[1L, ]
  innovation_scale <- sqrt(1 - ar1_phi^2)
  for (observation in 2:observations) {
    value[observation, ] <- ar1_phi * value[observation - 1L, ] +
      innovation_scale * innovations[observation, ]
  }
  value
}

run_replication <- function(regime, replication) {
  raw <- lapply(partitions, function(partition) {
    design %*% rbind(
      rep(0, features),
      rep(0, features),
      signal_pattern
    ) + generate_noise(regime)
  })
  names(raw) <- partitions
  fit <- lm_relation_fit(
    raw, design, effects, domain = domain,
    provenance = list(
      validation_regime = regime,
      replication = replication,
      observation_model = if (identical(regime, "well_specified")) {
        "gaussian_iid"
      } else {
        "gaussian_stationary_ar1_fitted_as_iid"
      }
    )
  )
  plans <- lapply(policies, function(policy) {
    plan_crossnobis(
      fit, frame, over,
      metric = recipe,
      training = policy,
      compute = compute_policy(workspace_bytes = 64 * 1024^2),
      residual_workspace_bytes = 64 * 1024^2
    )
  })
  if (length(unique(vapply(plans, `[[`, character(1), "lowering"))) != 1L ||
      length(unique(vapply(plans, `[[`, character(1),
        "kernel_version"))) != 1L ||
      length(unique(vapply(plans, function(plan) {
        plan$metric_schedule$statistics$signature
      }, character(1)))) != 1L ||
      length(unique(vapply(plans, function(plan) {
        plan$metric_schedule$signature
      }, character(1)))) != length(plans)) {
    stop("Policy plans do not share evidence/statistics execution with distinct identities.")
  }
  rows <- lapply(names(plans), function(policy_name) {
    plan <- plans[[policy_name]]
    evaluated <- lapply(contrasts, function(contrast) {
      crossnobis(plan, contrast)
    })
    if (any(vapply(evaluated, function(value) {
      !inherits(value, "effect_crossnobis_view") ||
        !identical(value$receipt$completion_status, "complete") ||
        !identical(value$receipt$kernel_version,
          "support-streamed-scheduled-metric-v1")
    }, logical(1)))) {
      stop("A policy estimate did not complete through the declared kernel.")
    }
    observed <- vapply(evaluated, function(value) {
      value$values[[selected_node]]
    }, numeric(1))
    conditional <- vapply(names(contrasts), function(estimand) {
      delta <- if (identical(estimand, "signal")) weighted_signal else {
        numeric(length(support))
      }
      sum(vapply(seq_len(nrow(over)), function(edge) {
        metric <- crossform:::materialize_metric(
          plan$metric_schedule, selected_node, edge
        )
        over$weight[[edge]] * drop(delta %*% metric$value %*% delta)
      }, numeric(1)))
    }, numeric(1))
    data.frame(
      regime = regime,
      replication = replication,
      policy = policy_name,
      estimand = names(observed),
      estimate = unname(observed),
      conditional_target = unname(conditional),
      population_target = unname(population_targets[names(observed)]),
      scientific_plan_id = plan$scientific_plan_id,
      metric_schedule_id = plan$metric_schedule$signature,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

set.seed(seed)
started <- proc.time()[["elapsed"]]
raw_rows <- vector("list", 2L * replications)
position <- 0L
for (regime in c("well_specified", "unmodelled_ar1")) {
  for (replication in seq_len(replications)) {
    position <- position + 1L
    raw_rows[[position]] <- run_replication(regime, replication)
    if (replication %% max(1L, floor(replications / 10L)) == 0L) {
      message(sprintf(
        "%s: %d/%d replications", regime, replication, replications
      ))
    }
  }
}
elapsed_seconds <- proc.time()[["elapsed"]] - started
raw_results <- do.call(rbind, raw_rows)
raw_results$evaluator_error <-
  raw_results$estimate - raw_results$conditional_target
raw_results$total_error <-
  raw_results$estimate - raw_results$population_target
raw_results$learner_target_error <-
  raw_results$conditional_target - raw_results$population_target

mean_se <- function(value) {
  c(mean = mean(value), se = stats::sd(value) / sqrt(length(value)))
}
summaries <- lapply(
  split(
    raw_results,
    interaction(
      raw_results$regime, raw_results$policy, raw_results$estimand,
      drop = TRUE
    )
  ),
  function(group) {
    estimate <- mean_se(group$estimate)
    evaluator <- mean_se(group$evaluator_error)
    total <- mean_se(group$total_error)
    learner <- mean_se(group$learner_target_error)
    data.frame(
      regime = group$regime[[1L]],
      policy = group$policy[[1L]],
      estimand = group$estimand[[1L]],
      replications = nrow(group),
      population_target = group$population_target[[1L]],
      mean_estimate = estimate[["mean"]],
      estimate_variance = stats::var(group$estimate),
      estimate_se = estimate[["se"]],
      conditional_target_bias = evaluator[["mean"]],
      conditional_target_bias_se = evaluator[["se"]],
      population_mahalanobis_bias = total[["mean"]],
      population_mahalanobis_bias_se = total[["se"]],
      learner_target_bias = learner[["mean"]],
      learner_target_bias_se = learner[["se"]],
      population_rmse = sqrt(mean(group$total_error^2)),
      stringsAsFactors = FALSE
    )
  }
)
summary_table <- do.call(rbind, summaries)
rownames(summary_table) <- NULL
summary_table <- summary_table[
  order(summary_table$regime, summary_table$estimand,
    summary_table$policy, method = "radix"),
  , drop = FALSE
]

paired_fields <- c(
  "regime", "replication", "estimand", "estimate",
  "conditional_target", "evaluator_error"
)
paired <- merge(
  raw_results[
    raw_results$policy == "training_only",
    paired_fields
  ],
  raw_results[
    raw_results$policy == "all_run_residuals",
    paired_fields
  ],
  by = c("regime", "replication", "estimand"),
  suffixes = c("_training_only", "_all_run")
)
paired_quantities <- do.call(rbind, lapply(
  c("estimate", "conditional_target", "evaluator_error"),
  function(quantity) {
    data.frame(
      regime = paired$regime,
      replication = paired$replication,
      estimand = paired$estimand,
      quantity = quantity,
      difference = paired[[paste0(quantity, "_all_run")]] -
        paired[[paste0(quantity, "_training_only")]],
      stringsAsFactors = FALSE
    )
  }
))
critical <- stats::qt(
  1 - (1 - confidence_level) / 2,
  df = replications - 1L
)
equivalence_rows <- lapply(
  split(
    paired_quantities,
    interaction(
      paired_quantities$regime, paired_quantities$estimand,
      paired_quantities$quantity, drop = TRUE
    )
  ),
  function(group) {
    difference <- mean_se(group$difference)
    lower <- difference[["mean"]] - critical * difference[["se"]]
    upper <- difference[["mean"]] + critical * difference[["se"]]
    data.frame(
      regime = group$regime[[1L]],
      estimand = group$estimand[[1L]],
      quantity = group$quantity[[1L]],
      replications = nrow(group),
      difference_direction = "all_run_minus_training_only",
      mean_difference = difference[["mean"]],
      paired_difference_se = difference[["se"]],
      confidence_level = confidence_level,
      ci_lower = lower,
      ci_upper = upper,
      equivalence_margin = equivalence_margin,
      equivalent_within_margin =
        lower > -equivalence_margin && upper < equivalence_margin,
      stringsAsFactors = FALSE
    )
  }
)
equivalence_table <- do.call(rbind, equivalence_rows)
rownames(equivalence_table) <- NULL
equivalence_table <- equivalence_table[
  order(equivalence_table$regime, equivalence_table$estimand,
    equivalence_table$quantity, method = "radix"),
  , drop = FALSE
]

scientific_contract <- list(
  validation_level = "statistical_recovery_and_policy_equivalence",
  seed = seed,
  replications = replications,
  regimes = list(
    well_specified = list(
      generator = "Gaussian iid observations under the fitted GLM",
      fitted_observation_metric = "identity",
      orthogonality_assumption = "satisfied"
    ),
    unmodelled_ar1 = list(
      generator = "Gaussian stationary AR(1) observations",
      phi = ar1_phi,
      fitted_observation_metric = "identity",
      orthogonality_assumption = "not guaranteed"
    )
  ),
  metric = list(
    estimator = recipe$hyperparameters$estimator,
    shrinkage = shrinkage,
    target = recipe$hyperparameters$target,
    target_estimated = recipe$hyperparameters$target_estimated,
    uncertainty_calibrated = FALSE
  ),
  equivalence = list(
    estimand = "all-run minus training-only learned-metric crossnobis",
    margin = equivalence_margin,
    confidence_level = confidence_level,
    rule = "two-sided confidence interval lies wholly within the margin"
  ),
  interpretation = list(
    conditional_target_bias =
      "effect-estimation and evaluator error for the realized learned metric",
    learner_target_bias =
      "learned-metric target minus population Mahalanobis target",
    population_mahalanobis_bias =
      "combined learned-metric and effect-estimation error"
  )
)
record <- list(
  provenance = provenance,
  contract = scientific_contract,
  fixture = list(
    observations_per_partition = observations,
    features = features,
    partitions = partitions,
    selected_node = selected_node,
    support = support,
    frame_weights = frame_weights,
    neural_covariance = neural_covariance,
    signal_pattern = signal_pattern,
    population_targets = population_targets,
    pairing = as.data.frame(over)
  ),
  summary = summary_table,
  equivalence = equivalence_table,
  raw = raw_results,
  elapsed_seconds = elapsed_seconds,
  session = utils::sessionInfo()
)

rds_path <- file.path(
  output_dir, "learned-metric-policy-validation.rds"
)
summary_path <- file.path(
  output_dir, "learned-metric-policy-summary.csv"
)
saveRDS(record, rds_path, compress = "xz")
export_columns <- union(names(summary_table), names(equivalence_table))
bind_export <- function(value, record_type) {
  missing <- setdiff(export_columns, names(value))
  value[missing] <- NA
  value <- value[export_columns]
  value$record_type <- record_type
  value
}
export_table <- rbind(
  bind_export(summary_table, "recovery"),
  bind_export(equivalence_table, "policy_equivalence")
)
utils::write.csv(export_table, summary_path, row.names = FALSE, na = "")

print(summary_table, row.names = FALSE, digits = 6)
cat("\nPolicy equivalence\n")
print(equivalence_table, row.names = FALSE, digits = 6)
cat(sprintf(
  "\nElapsed: %.1f seconds\nRDS: %s\nCSV: %s\n",
  elapsed_seconds, rds_path, summary_path
))
