# Frozen before the independent pilot and production calibration.
predictive_sampling_config <- list(
  contract = "predictive-geometry-v1", version = 1L,
  pilot_seed = 2026090413L, production_seed = 2026090403L,
  pilot_replicates = 2500L, minimum_replicates = 20000L,
  count_rounding = 1000L, pilot_variance_headroom = 1.5,
  confidence_multiplier = 5, relative_equivalence_margin = .02,
  pilot_target_fraction = .5, numerical_tolerance = 1e-10,
  arms = c("null", "aligned", "outside_model", "heterogeneous_metric"),
  predictors = c("model", "isotropic", "mismatched"),
  rank = 2L, penalty = .15,
  null_leakage_minimum = .005,
  # 12 identity checks and four generator-moment checks each use five MCSE.
  # Under the CLT the Bonferroni two-sided family bound is < 1e-5.
  familywise_rule = "Five MCSE for each of 12 identity and four generator-moment checks; no selected arm or seed.",
  selection_scope = "Fixed rank/penalty for every same-span comparator; zero predictor exact.")
