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
