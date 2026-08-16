test_that("study binds clocks and reorders confounds without neural reads", {
  fixture <- bound_study_fixture(use_functions = TRUE)
  value <- study(
    fixture$observations,
    fixture$events,
    fixture$confounds,
    fixture$hierarchy
  )

  expect_s3_class(value, "effect_study")
  expect_identical(fixture$fixture$reads$count, 0L)
  expect_true(value$capabilities$aligned_observations)
  expect_true(value$capabilities$timing_resolved)
  expect_identical(value$lineage$`run-1`$observation_id, 1:7)
  expect_identical(value$lineage$`run-1`$retained,
    c(TRUE, FALSE, rep(TRUE, 5L)))
  expect_match(format(value), "2 partitions; timing resolved")
  expect_identical(fixture$fixture$reads$count, 0L)
})

test_that("censor realization is factual lineage, not positional luck", {
  fixture <- bound_study_fixture()
  value <- study(
    fixture$observations,
    fixture$events,
    fixture$confounds,
    fixture$hierarchy
  )
  source_rows <- value$lineage$`run-2`$source_row
  confound_rows <- value$lineage$`run-2`$confound_row

  expect_identical(source_rows, seq_len(fixture$fixture$counts[[2L]]))
  expect_false(identical(source_rows, confound_rows))
  expect_identical(
    fixture$confounds$data$observation_id[confound_rows],
    value$lineage$`run-2`$observation_id
  )
})

test_that("missing and extra confound rows refuse by partition with remedies", {
  fixture <- bound_study_fixture()
  table <- fixture$confounds$data
  selected <- !(table$partition == "run-1" & table$observation_id == 3L)
  table <- table[selected, ]
  table <- rbind(table, data.frame(
    partition = "run-1",
    observation_id = 999L,
    motion = 0,
    retained = TRUE
  ))
  confounds <- observation_confounds(table, censor = "retained")
  refusal <- catch_refusal(study(
    fixture$observations,
    fixture$events,
    confounds,
    fixture$hierarchy
  ))

  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "aligned_observations")
  expect_identical(refusal$namespace, "study")
  expect_match(conditionMessage(refusal), "run-1")
  expect_true(any(grepl("Missing confound rows: 3", refusal$reasons)))
  expect_true(any(grepl("Extra confound rows: 999", refusal$reasons)))
  expect_true(length(refusal$remedies) >= 1L)
})

test_that("clock unit mismatch and missing clocks refuse explicitly", {
  fixture <- bound_study_fixture()
  wrong_units <- fixture$events
  wrong_units$units <- "milliseconds"
  wrong_units$events_id <- crossform:::.validate_nonempty_id(
    wrong_units$events_id, "events_id"
  )
  # Rebuild rather than relying on a mutated identity.
  wrong_units <- events(wrong_units$data, units = "milliseconds")
  unit_refusal <- catch_refusal(study(
    fixture$observations,
    wrong_units,
    fixture$confounds,
    fixture$hierarchy
  ))
  expect_identical(unit_refusal$capability, "timing_resolved")
  expect_match(conditionMessage(unit_refusal), "different units")
  expect_length(unit_refusal$reasons, 2L)

  indexes <- fixture$fixture$indexes
  indexes$`run-2` <- observation_index(
    indexes$`run-2`$observation_id, "run-2"
  )
  untimed_observations <- observations(
    fixture$fixture$sources,
    indexes,
    fixture$fixture$domain
  )
  missing_refusal <- catch_refusal(study(
    untimed_observations,
    fixture$events,
    fixture$confounds,
    fixture$hierarchy
  ))
  expect_identical(missing_refusal$capability, "timing_resolved")
  expect_match(conditionMessage(missing_refusal), "run-2")
  expect_match(missing_refusal$remedies, "observation times")
})

test_that("event coverage failures name events and acquisition range", {
  fixture <- bound_study_fixture()
  table <- fixture$events$data
  selected <- table$partition == "run-2" & table$event_id == "run-2-event-3"
  table$onset[selected] <- 1000
  outside <- events(table)
  refusal <- catch_refusal(study(
    fixture$observations,
    outside,
    fixture$confounds,
    fixture$hierarchy
  ))

  expect_identical(refusal$capability, "timing_resolved")
  expect_match(conditionMessage(refusal), "run-2")
  expect_true(any(grepl("run-2-event-3", refusal$reasons)))
  expect_true(any(grepl("Observation coverage", refusal$reasons)))
})

test_that("hierarchy supplies explicit axis vocabulary only", {
  fixture <- bound_study_fixture()
  value <- study(
    fixture$observations,
    fixture$events,
    fixture$confounds,
    fixture$hierarchy
  )
  run <- study_axis(value, "run")
  session <- study_axis(value, "session")

  expect_s3_class(run, "effect_study_axis")
  expect_identical(run$parent, "session")
  expect_false(identical(run$signature, session$signature))
  expect_false(any(c("generalizes_over", "independence", "sampling_unit") %in%
    names(value)))
  expect_false(any(c("generalizes_over", "independence") %in% names(run)))
  expect_match(format(run), "study_axis<run; 2 levels>")

  refusal <- catch_refusal(study_axis(value, "site"))
  expect_identical(refusal$capability, "declared_generalization_axis")
  expect_match(refusal$reasons, "Available axes")
})

test_that("untimed reduced facts remain valid but do not earn timing", {
  fixture <- study_fact_fixture()
  indexes <- lapply(fixture$indexes, function(index) {
    observation_index(index$observation_id, index$partition)
  })
  observations <- observations(fixture$sources, indexes, fixture$domain)
  event_table <- data.frame(
    partition = fixture$partitions,
    event_id = c("block-a", "block-b"),
    condition = c("face", "place")
  )
  events <- events(event_table, onset = NULL, duration = NULL)
  value <- study(observations, events)

  expect_false(value$capabilities$timing_resolved)
  expect_true(value$capabilities$aligned_observations)
})

test_that("unknown partitions and hierarchy order refuse before source reads", {
  fixture <- bound_study_fixture(use_functions = TRUE)
  bad_events <- fixture$events$data
  bad_events$partition[1] <- "missing-run"
  refusal <- catch_refusal(study(
    fixture$observations,
    events(bad_events),
    fixture$confounds,
    fixture$hierarchy
  ))
  expect_identical(refusal$capability, "aligned_observations")
  expect_match(refusal$reasons, "missing-run")
  expect_identical(fixture$fixture$reads$count, 0L)

  reversed <- partition_hierarchy(fixture$hierarchy$data[2:1, ])
  hierarchy_refusal <- catch_refusal(study(
    fixture$observations,
    fixture$events,
    fixture$confounds,
    reversed
  ))
  expect_identical(hierarchy_refusal$capability, "aligned_observations")
  expect_match(conditionMessage(hierarchy_refusal), "hierarchy leaf")
  expect_identical(fixture$fixture$reads$count, 0L)
})

test_that("study identity includes factual binding and detects tampering", {
  fixture <- bound_study_fixture()
  value <- study(
    fixture$observations,
    fixture$events,
    fixture$confounds,
    fixture$hierarchy
  )
  changed <- study(
    fixture$observations,
    fixture$events,
    fixture$confounds,
    fixture$hierarchy,
    clock_tolerance = 0.25
  )
  expect_false(identical(value$study_id, changed$study_id))

  tampered <- value
  tampered$lineage$`run-1`$retained[1] <- FALSE
  expect_error(crossform:::.validate_study(tampered),
    "lineage or capabilities")
  expect_true(study_capabilities(value)$stable_source_revision)
})
