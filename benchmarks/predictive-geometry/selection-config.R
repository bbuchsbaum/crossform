# Frozen before pilot/production outcomes are available.
predictive_selection_config <- list(contract = "predictive-geometry-v1", version = 1L,
  arms = c("null", "aligned"), pilot_seed = 2026090414L, production_seed = 2026090404L,
  pilot_replicates = 1000L, minimum_replicates = 2000L, count_rounding = 1000L,
  confidence_multiplier = 5, relative_equivalence_margin = .02,
  pilot_target_fraction = .5, pilot_variance_headroom = 1.5, numerical_tolerance = 1e-10,
  ranks = 0:2, penalties = c(0, .15, .5), inner_weights = c(.2, .3, .5),
  scope = "global", measurement_weights = 1,
  leakage_minimum = .001,
  familywise_rule = "Two conditional-risk checks at five MCSE; no post-hoc arms or seeds.")
