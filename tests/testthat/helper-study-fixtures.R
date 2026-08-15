study_fact_fixture <- function(use_functions = FALSE) {
  partitions <- c("run-1", "run-2")
  counts <- c(7L, 9L)
  domain <- abstract_domain(5L, id = "study-fixture-neural")
  indexes <- stats::setNames(lapply(seq_along(partitions), function(index) {
    observation_index(
      seq_len(counts[[index]]),
      partition = partitions[[index]],
      time = seq(0, by = 2, length.out = counts[[index]]),
      units = "seconds"
    )
  }), partitions)
  set.seed(2026081504)
  values <- lapply(counts, function(count) matrix(rnorm(count * 5L), count, 5L))
  names(values) <- partitions
  if (!use_functions) {
    sources <- values
    dimensions <- NULL
    capabilities <- NULL
    reads <- NULL
  } else {
    reads <- new.env(parent = emptyenv())
    reads$count <- 0L
    sources <- lapply(values, function(value) {
      force(value)
      function(features) {
        reads$count <- reads$count + 1L
        value[, features, drop = FALSE]
      }
    })
    names(sources) <- partitions
    dimensions <- lapply(values, dim)
    capabilities <- lapply(seq_along(partitions), function(index) {
      source_capabilities(
        block_read = TRUE,
        stable_revision = paste0(
          "sha256:", paste(rep(as.character(index), 64), collapse = "")
        )
      )
    })
    names(capabilities) <- partitions
  }
  list(
    partitions = partitions,
    counts = counts,
    domain = domain,
    indexes = indexes,
    sources = sources,
    source_dims = dimensions,
    capabilities = capabilities,
    reads = reads
  )
}

bound_study_fixture <- function(use_functions = FALSE) {
  fixture <- study_fact_fixture(use_functions)
  observation_record <- observations(
    fixture$sources,
    fixture$indexes,
    fixture$domain,
    source_dims = fixture$source_dims,
    capabilities = fixture$capabilities
  )
  event_table <- do.call(rbind, lapply(seq_along(fixture$partitions),
    function(index) {
      partition <- fixture$partitions[[index]]
      upper <- max(fixture$indexes[[partition]]$time)
      data.frame(
        partition = partition,
        event_id = paste0(partition, "-event-", 1:3),
        onset = c(0, upper / 3, upper * 2 / 3),
        duration = c(0.5, 1, 0.75),
        condition = c("face", "place", "object"),
        item = paste0("item-", 1:3),
        stringsAsFactors = FALSE
      )
    }))
  event_record <- events(event_table)

  confound_table <- do.call(rbind, lapply(rev(fixture$partitions),
    function(partition) {
      ids <- rev(fixture$indexes[[partition]]$observation_id)
      data.frame(
        partition = partition,
        observation_id = ids,
        motion = seq_along(ids) / 100,
        retained = ids != 2L,
        stringsAsFactors = FALSE
      )
    }))
  confound_record <- observation_confounds(
    confound_table, censor = "retained"
  )
  hierarchy <- partition_hierarchy(data.frame(
    partition = fixture$partitions,
    run = fixture$partitions,
    session = c("session-1", "session-2"),
    subject = "subject-1",
    stringsAsFactors = FALSE
  ))
  list(
    fixture = fixture,
    observations = observation_record,
    events = event_record,
    confounds = confound_record,
    hierarchy = hierarchy
  )
}
