test_that("benchmark scenarios span density and storage regimes", {
  scenarios <- effectagram:::.memory_benchmark_scenarios()

  expect_true(all(c("memory", "block") %in% scenarios$storage))
  expect_true(all(c("cold", "warm") %in% scenarios$phase))
  expect_true(any(scenarios$density == 1))
  expect_true(any(scenarios$density < 0.1))
  expect_true(all(scenarios$features > scenarios$effects))
  expect_identical(anyDuplicated(scenarios$id), 0L)
})

test_that("allocation benchmark emits a complete evidence record", {
  scenario <- effectagram:::.memory_benchmark_scenarios()[1, ]
  result <- effectagram:::.run_memory_benchmark_case(scenario)

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
    effectagram:::.parse_os_peak_rss(
      "Maximum resident set size (kbytes): 12345"
    ),
    12345 * 1024
  )
  expect_equal(
    effectagram:::.parse_os_peak_rss("  987654 maximum resident set size"),
    987654
  )
  expect_true(is.na(effectagram:::.parse_os_peak_rss("no measurement")))
})
