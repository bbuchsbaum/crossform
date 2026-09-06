# Every committed certification artifact is read by at least one test.
#
# The scale gates with their own vertical slice live in the test file for that
# slice. This file covers the remaining recorded evidence: the brain-scale
# crossnobis gate, the executor admission record, and the two large
# statistical validation records that are deliberately kept out of the
# tarball. The large records are read when a local run has produced them and
# skipped with a message when it has not.

test_that("recorded brain-scale crossnobis evidence passes its gate", {
  artifact <- certified_artifact(
    "crossnobis-scale-gate.rds", "run-crossnobis-scale-gate.R"
  )
  expect_identical(artifact$schema_version, 1L)
  expect_identical(artifact$provenance$runner, "run-crossnobis-scale-gate.R")
  expect_identical(artifact$fixture$id, "crossnobis-scale-52k:v1")
  expect_gte(artifact$fixture$features, 50000L)
  expect_identical(artifact$execution$completion_status, "complete")
  # Ticket B3 (09e18c4) retired `effect_crossnobis_plan`, so crossnobis plans
  # now carry the `geometry-sha256:` identity every geometry plan carries.
  expect_match(artifact$execution$scientific_plan_id,
    "^geometry-sha256:[[:xdigit:]]{64}$")

  # Re-derive the recorded verdicts from the recorded measurements.
  expect_lte(
    artifact$memory$incremental_peak_rss_bytes,
    artifact$gate$maximum_incremental_rss_bytes
  )
  expect_lte(
    artifact$timing$analysis_seconds,
    artifact$gate$maximum_analysis_seconds
  )
  expect_true(artifact$output$exact_full_index_mapping)
  expect_true(artifact$output$finite_values)
  expect_identical(artifact$output$output_values, artifact$fixture$features)

  # The structural claim is the scientific one: no pair-atom field, no p-by-p
  # pair frame, and no node-by-edge factor table is ever retained.
  expect_false(artifact$work$pair_atoms_materialized)
  expect_false(artifact$work$pair_frame_materialized)
  expect_false(artifact$work$factor_table_retained)
  expect_identical(artifact$memory$persistent_factor_table_bytes, 0)

  expect_true(artifact$gate$memory_pass)
  expect_true(artifact$gate$runtime_pass)
  expect_true(artifact$gate$mapping_pass)
  expect_true(artifact$gate$structure_pass)
  expect_identical(
    artifact$gate$passed,
    all(unlist(artifact$gate[c("memory_pass", "runtime_pass", "mapping_pass",
      "structure_pass")], use.names = FALSE))
  )
  expect_true(artifact$gate$passed)
})

test_that("recorded measurement-form evidence retains the BLAS route", {
  artifact <- certified_artifact(
    "measurement-profile.rds", "run-measurement-profile.R"
  )
  expect_identical(artifact$schema_version, 1L)
  expect_identical(artifact$provenance$runner, "run-measurement-profile.R")
  expect_identical(artifact$dimensions$experimental_conditions, 16L)
  expect_identical(artifact$dimensions$domain_features, 128L)
  expect_gte(artifact$dimensions$scalar_nodes, 20L)
  expect_gte(artifact$dimensions$multivariate_nodes, 10L)
  expect_gte(artifact$dimensions$partitions, 3L)

  records <- artifact$records
  expect_gte(nrow(records), 10L)
  expect_setequal(records$frame, c("scalar", "requested_multivariate"))
  expect_true(all(is.finite(records$median_elapsed_seconds)))
  expect_true(all(records$median_elapsed_seconds > 0))
  expect_true(all(is.finite(records$peak_rss_bytes)))
  expect_true(all(records$planned_workspace_bytes > records$output_bytes))
  expect_lt(max(records$oracle_error), 1e-10)
  expect_true(artifact$checks$oracle_parity)
  expect_true(artifact$checks$plan_identity_stable)
  expect_true(artifact$checks$planned_workspace_exceeds_output)

  expect_true(is.finite(artifact$rprof$max_r_loop_share))
  expect_lt(artifact$rprof$max_r_loop_share, 0.15)
  expect_identical(artifact$decision, "no_rcpp_keep_blas")
  expect_match(artifact$decision_rule, "R-level loops are >= 15%")
  expect_match(artifact$decision_rule, "prototype projects >= 1.25x")
})

test_that("the fused pair kernel clears its cumulative-allocation court", {
  artifact <- certified_artifact(
    "native-pair-allocation.rds", "run-native-pair-allocation.R"
  )
  expect_identical(artifact$schema_version, 1L)
  expect_identical(
    artifact$provenance$runner, "run-native-pair-allocation.R"
  )
  expect_identical(artifact$dimensions$conditions, 48L)
  expect_identical(artifact$dimensions$features, 512L)
  expect_identical(artifact$dimensions$partitions, 4L)
  expect_identical(artifact$dimensions$ordered_edges, 12L)
  expect_identical(artifact$dimensions$pairs, 1128L)
  expect_gte(artifact$dimensions$repetitions, 5L)

  records <- artifact$records
  expect_identical(
    nrow(records), 2L * artifact$dimensions$repetitions
  )
  expect_setequal(records$route, c("native_fused", "r_two_pass_oracle"))
  expect_true(all(table(records$repetition, records$route) == 1L))
  expect_true(all(records$allocation_count > 0))
  expect_true(all(records$allocated_bytes >= records$output_bytes))
  expect_true(all(is.finite(records$elapsed_seconds)))
  expect_true(all(records$elapsed_seconds > 0))
  expect_equal(
    records$checksum,
    rep(artifact$numerical$checksum, nrow(records)),
    tolerance = 1e-12
  )
  expect_lte(
    artifact$numerical$max_abs_error,
    artifact$numerical$tolerance
  )

  median_for <- function(column, route) {
    stats::median(records[[column]][records$route == route])
  }
  native_allocated <- median_for("allocated_bytes", "native_fused")
  oracle_allocated <- median_for("allocated_bytes", "r_two_pass_oracle")
  allocation_ratio <- native_allocated / oracle_allocated
  expect_equal(
    artifact$summary$native_median_allocated_bytes, native_allocated,
    tolerance = 0
  )
  expect_equal(
    artifact$summary$oracle_median_allocated_bytes, oracle_allocated,
    tolerance = 0
  )
  expect_equal(
    artifact$summary$native_to_oracle_allocation_ratio, allocation_ratio,
    tolerance = 0
  )
  expect_identical(artifact$gate$maximum_allocation_ratio, 0.70)
  expect_true(artifact$gate$numerical_pass)
  expect_identical(
    artifact$gate$allocation_pass,
    allocation_ratio <= artifact$gate$maximum_allocation_ratio
  )
  expect_true(artifact$gate$allocation_pass)
  expect_true(artifact$gate$passed)
})

test_that("the recorded executor admission record refuses the adapter", {
  artifact <- certified_artifact(
    "shard-admission.rds", "run-shard-admission.R"
  )
  expect_identical(artifact$schema_version, 1L)
  expect_identical(artifact$provenance$runner, "run-shard-admission.R")

  # Parity and cleanup are the correctness half of admission and must hold.
  # The row guard keeps an empty summary from certifying vacuously through
  # `all(logical(0))`.
  expect_gte(nrow(artifact$summary), 5L)
  expect_setequal(artifact$summary$phase, c("cold", "warm"))
  expect_true(any(artifact$summary$parallel))
  expect_true(all(artifact$summary$total_parity))
  expect_true(all(artifact$summary$coherent_parity))
  expect_true(all(artifact$summary$pool_stopped))
  expect_true(all(artifact$summary$backing_removed))
  expect_true(artifact$cleanup_stress$caught)
  expect_true(artifact$cleanup_stress$pool_stopped)
  expect_true(artifact$cleanup_stress$backing_removed)
  expect_true(artifact$gate$numerical_parity)
  expect_true(artifact$gate$cleanup)
  expect_true(artifact$gate$cleanup_stress)

  # The runtime half did not clear its threshold, so the adapter is not
  # admitted. This test pins the refusal: an artifact recording admission
  # without a qualifying cold speedup would fail here.
  cold <- artifact$summary$parallel & artifact$summary$phase == "cold"
  expect_identical(
    artifact$gate$cold_speedup,
    any(artifact$summary$speedup_vs_sequential[cold] >
      artifact$gate$cold_speedup_threshold)
  )
  expect_identical(
    artifact$gate$admitted,
    artifact$gate$numerical_parity && artifact$gate$cleanup &&
      artifact$gate$cleanup_stress && artifact$gate$cold_speedup
  )
})

test_that("the recorded sampling-covariance validation is calibrated", {
  artifact <- certified_artifact(
    "sampling-covariance-validation.rds",
    "run-sampling-covariance-validation.R"
  )
  expect_identical(artifact$schema_version, 1L)
  expect_identical(artifact$equation, "Diedrichsen_et_al_2016_Eq13")
  expect_gte(artifact$repetitions, 1000L)

  # Bands are the ones the fast generative court declares for the same
  # estimator in test-evidence-sampling-generative.R.
  summary <- artifact$summary
  expect_setequal(summary$arm, c("null", "nonzero_signal"))
  expect_identical(nrow(summary), 6L)
  expect_true(all(summary$variance_ratio > 0.95 & summary$variance_ratio < 1.05))
  expect_true(all(summary$analytic_coverage > 0.93 &
    summary$analytic_coverage < 0.97))
  expect_lt(artifact$signal_full_covariance_relative_frobenius_error, 0.05)
  expect_lt(artifact$null_full_covariance_relative_frobenius_error, 0.05)
  expect_true(artifact$linear_variance_ratio > 0.95 &&
    artifact$linear_variance_ratio < 1.05)
  expect_true(artifact$linear_coverage > 0.93 && artifact$linear_coverage < 0.97)

  # The naive edge interval is materially anticonservative on the signal arm;
  # that gap is the reason the analytic law exists.
  signal <- summary[summary$arm == "nonzero_signal", ]
  expect_identical(nrow(signal), 3L)
  expect_true(all(signal$naive_edge_coverage < 0.82))
  expect_true(all(signal$analytic_coverage - signal$naive_edge_coverage > 0.13))
  expect_true(all(signal$plugin_coverage > 0.95 & signal$plugin_coverage < 0.995))
  expect_true(all(is.na(summary$plugin_coverage[summary$arm == "null"])))
})

test_that("the recorded learned-metric policy validation is equivalence-shaped", {
  artifact <- certified_artifact(
    "learned-metric-policy-validation.rds",
    "run-learned-metric-policy-validation.R"
  )
  expect_identical(artifact$contract$validation_level,
    "statistical_recovery_and_policy_equivalence")
  expect_identical(artifact$contract$replications, 500L)
  expect_identical(artifact$contract$equivalence$margin, 0.005)
  expect_false(artifact$contract$metric$uncertainty_calibrated)
  expect_setequal(names(artifact$contract$regimes),
    c("well_specified", "unmodelled_ar1"))

  # The declared rule: equivalent exactly when the interval lies wholly inside
  # the margin. Recorded verdicts are re-derived rather than believed.
  equivalence <- artifact$equivalence
  expect_identical(nrow(equivalence), 12L)
  expect_setequal(equivalence$quantity,
    c("conditional_target", "estimate", "evaluator_error"))
  expect_true(all(equivalence$equivalence_margin == 0.005))
  expect_true(all(equivalence$confidence_level == 0.95))
  expect_true(all(equivalence$replications == artifact$contract$replications))
  expect_identical(
    equivalence$equivalent_within_margin,
    equivalence$ci_lower > -equivalence$equivalence_margin &
      equivalence$ci_upper < equivalence$equivalence_margin
  )
  expect_true(all(equivalence$difference_direction ==
    "all_run_minus_training_only"))

  estimate <- equivalence[equivalence$quantity == "estimate", ]
  verdict <- stats::setNames(estimate$equivalent_within_margin,
    paste(estimate$regime, estimate$estimand, sep = ":"))
  expect_identical(verdict[["well_specified:null"]], TRUE)
  expect_identical(verdict[["well_specified:signal"]], TRUE)
  expect_identical(verdict[["unmodelled_ar1:null"]], TRUE)
  # Under unmodelled AR(1) the two residual policies are NOT equivalent on the
  # signal arm, and the difference is attributable to the learned metric
  # rather than to a second evaluation engine: evaluator error stays inside
  # the margin while the conditional learned target does not.
  expect_identical(verdict[["unmodelled_ar1:signal"]], FALSE)
  misspecified <- equivalence[equivalence$regime == "unmodelled_ar1" &
    equivalence$estimand == "signal", ]
  expect_identical(nrow(misspecified), 3L)
  expect_true(misspecified$equivalent_within_margin[
    misspecified$quantity == "evaluator_error"
  ])
  expect_false(misspecified$equivalent_within_margin[
    misspecified$quantity == "conditional_target"
  ])

  # Two regimes, two residual policies, and two estimands per replication.
  expect_identical(nrow(artifact$raw),
    8L * as.integer(artifact$contract$replications))
  # Ticket B3 (09e18c4) retired `effect_crossnobis_plan`; the identity scheme
  # is `geometry-sha256:` for these plans now.
  expect_true(all(grepl("^geometry-sha256:[[:xdigit:]]{64}$",
    artifact$raw$scientific_plan_id)))
})
