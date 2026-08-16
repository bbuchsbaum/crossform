test_that("benchmark scenarios span density and storage regimes", {
  scenarios <- crossform:::.memory_benchmark_scenarios()

  expect_true(all(c("memory", "block") %in% scenarios$storage))
  expect_true(all(c("cold", "warm") %in% scenarios$phase))
  expect_true(any(scenarios$density == 1))
  expect_true(any(scenarios$density < 0.1))
  expect_true(all(scenarios$features > scenarios$effects))
  expect_identical(anyDuplicated(scenarios$id), 0L)
})

test_that("allocation benchmark emits a complete evidence record", {
  scenario <- crossform:::.memory_benchmark_scenarios()[1, ]
  result <- crossform:::.run_memory_benchmark_case(scenario)

  expect_identical(result$schema_version, 3L)
  expect_gt(result$allocation$allocation_count, 0)
  expect_gt(result$allocation$allocated_bytes, 0)
  expect_gt(result$allocation$largest_allocation_bytes, 0)
  expect_s3_class(result$plan, "effect_memory_plan")
  expect_gt(result$total$relation_reads, 0)
  expect_gt(result$total$max_live_temporary_bytes, 0)
  expect_gt(result$total$durable_local_relation_bytes, 0)
  expect_gt(result$coherent$max_tile_bytes, 0)
  expect_gt(result$plan$categories[["frame"]], 0)
  expect_gt(result$plan$categories[["resident_source"]], 0)
  expect_true(is.finite(result$baseline_rss_bytes))
  expect_true(is.finite(result$sampled_incremental_peak_rss_bytes))
  expect_equal(result$worker_peak_rss_bytes, 0)
  expect_identical(result$rss_evidence,
    "sampled_only_until_child_process_runner_attaches_os_peak")
})

test_that("OS peak RSS parser recognizes GNU and macOS time formats", {
  expect_equal(
    crossform:::.parse_os_peak_rss(
      "Maximum resident set size (kbytes): 12345"
    ),
    12345 * 1024
  )
  expect_equal(
    crossform:::.parse_os_peak_rss("  987654 maximum resident set size"),
    987654
  )
  expect_true(is.na(crossform:::.parse_os_peak_rss("no measurement")))
})

test_that("recorded memory-benchmark evidence covers every declared scenario", {
  scenarios <- crossform:::.memory_benchmark_scenarios()
  for (index in seq_len(nrow(scenarios))) {
    id <- scenarios$id[[index]]
    artifact <- certified_artifact(
      paste0(id, ".rds"), "run-memory-benchmarks.R"
    )
    expect_identical(artifact$schema_version, 3L)
    expect_identical(artifact$provenance$runner, "run-memory-benchmarks.R")
    expect_identical(artifact$scenario$id, id)
    expect_identical(artifact$scenario$storage, scenarios$storage[[index]])
    expect_identical(artifact$scenario$phase, scenarios$phase[[index]])

    # The recorded RSS is a real OS measurement of an isolated child, taken
    # only after the fixture signalled ready.
    expect_identical(artifact$rss_evidence,
      "isolated_child_process_polled_os_peak")
    expect_gt(artifact$os_peak_rss_bytes, 0)
    expect_gte(artifact$incremental_peak_rss_bytes, 0)
    expect_true(is.finite(artifact$incremental_peak_rss_bytes))
    # Version 0.1 runs no workers, so an aggregated worker peak of zero is a
    # fact about this release rather than an unmeasured field.
    expect_identical(artifact$worker_peak_rss_bytes, 0)

    # The plan predicts crossform-owned workspace, not process RSS. The
    # prediction is checked against what it actually claims to bound.
    expect_s3_class(artifact$plan, "effect_memory_plan")
    expect_identical(artifact$plan$prediction_kind,
      "crossform_owned_workspace_upper_bound")
    expect_gte(artifact$plan$planned_workspace_bytes,
      artifact$plan$modeled_workspace_bytes)
    expect_gte(artifact$plan$planned_workspace_bytes,
      artifact$total$max_live_temporary_bytes +
        artifact$total$durable_local_relation_bytes)
    expect_gt(artifact$allocation$allocation_count, 0)
    expect_gt(artifact$allocation$allocated_bytes, 0)
  }
})
