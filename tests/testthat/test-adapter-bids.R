.write_bids_tsv <- function(value, path) {
  utils::write.table(
    value, path, sep = "\t", row.names = FALSE, quote = FALSE,
    na = "n/a"
  )
  path
}

test_that("BIDS adapters preserve facts without assigning model roles", {
  bound <- bound_study_fixture()
  directory <- tempfile("crossform-bids-")
  dir.create(directory)
  event_files <- confound_files <- character(length(bound$fixture$partitions))
  names(event_files) <- names(confound_files) <- bound$fixture$partitions

  for (partition in bound$fixture$partitions) {
    event_rows <- bound$events$data$partition == partition
    event_table <- bound$events$data[event_rows,
      c("onset", "duration", "condition", "item")]
    event_table$trial_type <- event_table$condition
    event_files[[partition]] <- .write_bids_tsv(
      event_table, file.path(directory, paste0(partition, "_events.tsv"))
    )

    count <- bound$fixture$counts[[match(partition,
      bound$fixture$partitions)]]
    confound_table <- data.frame(
      trans_x = seq_len(count) / 100,
      non_steady_state_outlier00 = c(TRUE, rep(FALSE, count - 1L)),
      keep = seq_len(count) != 2L,
      stringsAsFactors = FALSE
    )
    confound_files[[partition]] <- .write_bids_tsv(
      confound_table,
      file.path(directory, paste0(partition, "_desc-confounds_timeseries.tsv"))
    )
  }

  event_record <- crossform:::bids_events(event_files, bound$fixture$partitions)
  expect_true(all(c("trial_type", "item") %in% names(event_record$data)))
  expect_false(any(c("role", "nuisance", "target") %in%
    names(event_record$data)))

  confound_record <- crossform:::bids_confounds(
    confound_files,
    bound$fixture$partitions,
    lapply(bound$fixture$indexes, `[[`, "observation_id"),
    censor = "keep"
  )
  expect_true("non_steady_state_outlier00" %in% names(confound_record$data))
  expect_identical(confound_record$censor_column, "keep")

  imported <- bids_study(
    bound$observations,
    event_files,
    confound_files,
    partitions = bound$fixture$partitions,
    censor = "keep",
    hierarchy = bound$hierarchy
  )
  expect_s3_class(imported, "effect_study")
  expect_identical(imported$partitions, bound$fixture$partitions)
  expect_equal(unname(vapply(imported$lineage,
    function(x) sum(x$retained), integer(1))),
  unname(bound$fixture$counts - 1L))
})

test_that("BIDS event timing and censoring ambiguities refuse by capability", {
  directory <- tempfile("crossform-bids-refusal-")
  dir.create(directory)
  invalid_events <- .write_bids_tsv(
    data.frame(onset = "n/a", duration = 1, trial_type = "face"),
    file.path(directory, "events.tsv")
  )
  timing <- catch_refusal(crossform:::bids_events(invalid_events, "run-1"))
  expect_s3_class(timing, "effect_capability_refusal")
  expect_identical(timing$capability, "timing_resolved")

  confounds <- .write_bids_tsv(
    data.frame(trans_x = c(0.1, 0.2), keep = c(1L, 0L)),
    file.path(directory, "confounds.tsv")
  )
  censoring <- catch_refusal(crossform:::bids_confounds(
    confounds, "run-1", list(`run-1` = 1:2), censor = "keep"
  ))
  expect_s3_class(censoring, "effect_capability_refusal")
  expect_identical(censoring$capability, "censoring_declared")
})

test_that("BIDS events bind runs with different columns and keep labels as text", {
  directory <- tempfile("crossform-bids-columns-")
  dir.create(directory)
  first <- file.path(directory, "run-1_events.tsv")
  second <- file.path(directory, "run-2_events.tsv")
  writeLines(c(
    "onset\tduration\ttrial_type\tresponse_time",
    "0\t0.5\t01\t0.8",
    "4\t0.5\t02\tn/a"
  ), first)
  writeLines(c(
    "onset\tduration\ttrial_type",
    "2\t1\t01"
  ), second)
  record <- crossform:::bids_events(
    c(`run-1` = first, `run-2` = second)
  )
  expect_identical(record$data$trial_type, c("01", "02", "01"))
  expect_identical(record$data$response_time, c("0.8", "n/a", "n/a"))
  expect_identical(record$data$onset, c(0, 4, 2))
  expect_identical(record$data$duration, c(0.5, 0.5, 1))
  expect_identical(record$data$.bids_partition, c("run-1", "run-1", "run-2"))
})

test_that("BIDS header-only events refuse with a clear reason", {
  path <- tempfile(fileext = ".tsv")
  writeLines("onset\tduration\ttrial_type", path)
  refusal <- catch_refusal(crossform:::bids_events(c(`run-1` = path)))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "timing_resolved")
  expect_match(conditionMessage(refusal), "no events")
})

test_that("BIDS confounds type complete columns and bind differing columns", {
  directory <- tempfile("crossform-bids-confounds-")
  dir.create(directory)
  first <- file.path(directory, "run-1_confounds.tsv")
  second <- file.path(directory, "run-2_confounds.tsv")
  writeLines(c(
    "framewise_displacement\ttrans_x\trot_x\tkeep",
    "n/a\t0.1\t1e-3\tTRUE",
    "0.2\t0.3\t-2\tFALSE"
  ), first)
  writeLines(c(
    "framewise_displacement\tkeep\trot_x\tmotion_outlier00",
    "n/a\tTRUE\t0.5\t1",
    "0.4\tTRUE\t0\t0"
  ), second)
  record <- crossform:::bids_confounds(
    c(`run-1` = first, `run-2` = second), censor = "keep"
  )
  # The BIDS missing marker is read as `NA` without changing the column's
  # type, including a column one run lacks altogether.
  expect_identical(record$data$framewise_displacement, c(NA, 0.2, NA, 0.4))
  expect_identical(record$data$rot_x, c(1e-3, -2, 0.5, 0))
  expect_identical(record$data$trans_x, c(0.1, 0.3, NA, NA))
  expect_identical(record$data$motion_outlier00, c(NA, NA, 1, 0))
  expect_identical(record$data$keep, c(TRUE, FALSE, TRUE, TRUE))

  # The same tables bind into a study; the missing values never reach the
  # censor policy, which reads only the complete retain column.
  value <- study(
    observations(
      list(`run-1` = matrix(0, 2L, 1L), `run-2` = matrix(0, 2L, 1L)),
      list(
        `run-1` = observation_index(1:2, "run-1"),
        `run-2` = observation_index(1:2, "run-2")
      ),
      abstract_domain(1L)
    ),
    confounds = record
  )
  expect_identical(value$confounds$data$framewise_displacement,
    c(NA, 0.2, NA, 0.4))
})

test_that("BIDS confounds refuse a censor column with missing values", {
  path <- tempfile("crossform-bids-censor-", fileext = ".tsv")
  writeLines(c(
    "framewise_displacement\tkeep",
    "n/a\tTRUE",
    "0.2\tn/a"
  ), path)
  refusal <- catch_refusal(
    crossform:::bids_confounds(c(`run-1` = path), censor = "keep")
  )
  expect_identical(refusal$capability, "censoring_declared")
  expect_match(paste(refusal$reasons, collapse = " "), "missing values",
    fixed = TRUE)
})

test_that("observation confounds allow missing values only in value columns", {
  confounds <- data.frame(
    partition = "run-1", observation_id = 1:3,
    framewise_displacement = c(NA, 0.2, 0.3),
    retained = c(TRUE, TRUE, FALSE)
  )
  record <- observation_confounds(confounds, censor = "retained")
  expect_identical(record$data$framewise_displacement, c(NA, 0.2, 0.3))
  expect_identical(
    crossform:::.validate_observation_confounds(record)$confounds_id,
    record$confounds_id
  )
  for (column in c("partition", "observation_id", "retained")) {
    broken <- confounds
    broken[[column]][[2L]] <- NA
    expect_error(observation_confounds(broken, censor = "retained"),
      column, class = "effect_input_error")
  }
  confounds$framewise_displacement[[2L]] <- Inf
  expect_error(observation_confounds(confounds, censor = "retained"),
    "framewise_displacement", class = "effect_input_error")
  confounds$framewise_displacement[[2L]] <- NaN
  expect_error(observation_confounds(confounds, censor = "retained"),
    "framewise_displacement", class = "effect_input_error")
})
