test_that("observation indexes bind ordered identifiers and optional clocks", {
  timed <- observation_index(letters[1:4], "run-1", c(0, 1.5, 3, 4.5), "s")
  untimed <- observation_index(1:3, "trial-block")

  expect_s3_class(timed, "effect_observation_index")
  expect_true(timed$timing)
  expect_false(untimed$timing)
  expect_false(identical(timed$signature, untimed$signature))
  expect_error(observation_index(1:3, "run", c(0, 2, 1), "s"),
    "strictly increasing", class = "effect_input_error")
  expect_error(observation_index(1:3, "run", units = "s"),
    "requires an observation", class = "effect_input_error")
})

test_that("observations support unequal partitions on one exact domain", {
  fixture <- study_fact_fixture()
  record <- observations(
    fixture$sources,
    fixture$indexes,
    fixture$domain
  )

  expect_s3_class(record, "effect_observations")
  expect_identical(record$partitions, fixture$partitions)
  expect_identical(unname(vapply(record$sources, function(value) {
    value$dim[[1L]]
  }, integer(1))), fixture$counts)
  expect_identical(record$domain, fixture$domain$reference)
  expect_true(all(vapply(record$capabilities, function(value) {
    crossform:::.strong_sha256(value$stable_revision)
  }, logical(1))))
})

test_that("metadata construction does not read lazy neural values", {
  fixture <- study_fact_fixture(use_functions = TRUE)
  record <- observations(
    fixture$sources,
    fixture$indexes,
    fixture$domain,
    source_dims = fixture$source_dims,
    capabilities = fixture$capabilities
  )

  expect_identical(fixture$reads$count, 0L)
  crossform:::.validate_observations(record)
  expect_identical(fixture$reads$count, 0L)
})

test_that("observation sources fail on row and neural-domain disagreement", {
  fixture <- study_fact_fixture()
  wrong_rows <- fixture$sources
  wrong_rows[["run-1"]] <- wrong_rows[["run-1"]][-1, , drop = FALSE]
  expect_error(
    observations(wrong_rows, fixture$indexes, fixture$domain),
    "source dimensions"
  , class = "effect_input_error")

  wrong_domain <- abstract_domain(4L, id = "wrong")
  expect_error(
    observations(fixture$sources, fixture$indexes, wrong_domain),
    "source dimensions"
  , class = "effect_input_error")
})

test_that("event tables carry schema but no permanent model roles", {
  table <- data.frame(
    partition = rep(c("run-1", "run-2"), each = 2L),
    event_id = paste0("event-", 1:4),
    onset = c(0, 4, 2, 8),
    duration = c(1, 1.5, 2, 1),
    condition = factor(c("face", "place", "face", "object")),
    motion_modulator = c(0.1, -0.2, 0.3, -0.1)
  )
  record <- observation_events(table)

  expect_s3_class(record, "effect_events")
  expect_true(record$timing)
  expect_identical(record$schema$condition$levels,
    levels(table$condition))
  expect_false(any(c("target", "nuisance", "role") %in% names(record)))

  untimed <- observation_events(
    table[, c("partition", "event_id", "condition")],
    onset = NULL, duration = NULL
  )
  expect_false(untimed$timing)
})

test_that("observation confounds remain on their own typed axis", {
  table <- data.frame(
    partition = rep(c("run-1", "run-2"), c(3, 2)),
    observation_id = c(1:3, 1:2),
    motion = c(0.1, 0.2, -0.1, 0.3, 0.4),
    retained = c(TRUE, FALSE, TRUE, TRUE, TRUE)
  )
  record <- observation_confounds(table, censor = "retained")

  expect_s3_class(record, "effect_observation_confounds")
  expect_identical(record$censor_column, "retained")
  expect_error(observation_confounds(rbind(table, table[1, ]),
    censor = "retained"), "unique", class = "effect_input_error")
})

test_that("partition hierarchies validate nesting and retain axis order", {
  hierarchy <- partition_hierarchy(data.frame(
    partition = c("run-1", "run-2", "run-3"),
    run = c("run-1", "run-2", "run-3"),
    session = c("session-1", "session-1", "session-2"),
    subject = "subject-1"
  ))

  expect_s3_class(hierarchy, "effect_partition_hierarchy")
  expect_identical(hierarchy$axes,
    c("partition", "run", "session", "subject"))
  expect_identical(unname(hierarchy$parent_maps$run),
    c("session-1", "session-1", "session-2"))

  bad <- data.frame(
    partition = c("p1", "p2"),
    run = c("run-1", "run-1"),
    session = c("s1", "s2")
  )
  expect_error(partition_hierarchy(bad), "does not nest uniquely",
    class = "effect_input_error")
})

test_that("fact identities detect mutation and semantic changes", {
  fixture <- study_fact_fixture()
  first <- observations(fixture$sources, fixture$indexes, fixture$domain)
  changed_index <- fixture$indexes
  changed_index$`run-1` <- observation_index(
    fixture$indexes$`run-1`$observation_id,
    "run-1",
    fixture$indexes$`run-1`$time + 0.25,
    "seconds"
  )
  second <- observations(fixture$sources, changed_index, fixture$domain)
  expect_false(identical(first$observations_id, second$observations_id))

  tampered <- first
  tampered$indexes$`run-1`$partition <- "changed"
  expect_error(crossform:::.validate_observations(tampered),
    "inconsistent|declares", class = "effect_contract_error")
})
