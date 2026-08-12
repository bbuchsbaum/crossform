args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: shard-admission-worker.R REPO MODE RESULT READY LIBRARY")
}
repo <- normalizePath(args[[1L]], mustWork = TRUE)
mode <- args[[2L]]
result_path <- args[[3L]]
ready_path <- args[[4L]]
library_path <- normalizePath(args[[5L]], mustWork = TRUE)
.libPaths(c(library_path, .libPaths()))

library(effectagram)
if (mode != "sequential-response") library(shard)

rss_self <- function() {
  lines <- system2("ps", c("-o", "rss=", "-p", as.character(Sys.getpid())),
    stdout = TRUE)
  as.double(trimws(lines[[1L]])) * 1024
}

set.seed(20260812)
scenario <- list(
  observations = 72L,
  effects = 8L,
  features = 49152L,
  partitions = 4L,
  measurements = 1536L,
  support = 64L,
  feature_block = 6144L,
  row_tile = 256L,
  coordinate_tile = 36L,
  workers = 2L
)
n <- scenario$observations
q <- scenario$effects
p <- scenario$features
r <- scenario$partitions
m <- scenario$measurements
partitions <- paste0("run", seq_len(r))
effects <- paste0("effect", seq_len(q))

responses <- matrix(rnorm(r * n * p), nrow = r * n, ncol = p)
extractors <- lapply(seq_len(r), function(index) {
  design <- matrix(rnorm(n * q), nrow = n, ncol = q)
  qr.solve(crossprod(design), t(design))
})
over <- cross_partitions(partitions)

rows <- rep(seq_len(m), each = scenario$support)
columns <- as.integer(unlist(lapply(seq_len(m), function(row) {
  start <- ((row - 1L) * 31L) %% p
  ((start + seq_len(scenario$support) - 1L) %% p) + 1L
}), use.names = FALSE))
weights <- Matrix::sparseMatrix(
  i = rows, j = columns, x = rep(1 / scenario$support, length(rows)),
  dims = c(m, p)
)
at <- additive_frame(weights)
feature_blocks <- split(seq_len(p),
  ceiling(seq_len(p) / scenario$feature_block))

materialize_relations <- function() {
  relations <- lapply(seq_len(r), function(index) {
    source_rows <- ((index - 1L) * n + 1L):(index * n)
    extractors[[index]] %*% responses[source_rows, , drop = FALSE]
  })
  do.call(rbind, relations)
}

task_from_response <- function(feature_ids, source) {
  relations <- stats::setNames(lapply(seq_len(r), function(index) {
    source_rows <- ((index - 1L) * n + 1L):(index * n)
    extractors[[index]] %*% source[source_rows, feature_ids, drop = FALSE]
  }), partitions)
  effectagram:::.crossgram_feature_task(
    relations, feature_ids, effects, partitions, over
  )
}

task_from_relation <- function(feature_ids, source) {
  relations <- stats::setNames(lapply(seq_len(r), function(index) {
    source_rows <- ((index - 1L) * q + 1L):(index * q)
    source[source_rows, feature_ids, drop = FALSE]
  }), partitions)
  effectagram:::.crossgram_feature_task(
    relations, feature_ids, effects, partitions, over
  )
}

reduce_tasks <- function(tasks) {
  reduced <- effectagram:::.reduce_crossgram_tasks(
    tasks, at, partitions, effects,
    output_width = q * (q + 1L) / 2L,
    row_tile = scenario$row_tile,
    coordinate_tile = scenario$coordinate_tile,
    retain_local_relations = TRUE
  )
  coherent <- effectagram:::.coherent_geometry_from_local(
    reduced$local_relations, over, Matrix::rowSums(weights),
    row_tile = scenario$row_tile
  )
  list(total = reduced$value, coherent = coherent$value,
    diagnostics = reduced$diagnostics)
}

shared <- NULL
pool <- NULL
cleanup <- list(pool_stopped = NA, source_closed = NA, backing_removed = NA)
on.exit({
  if (!is.null(pool)) {
    try(shard::pool_stop(pool), silent = TRUE)
  }
  if (!is.null(shared)) {
    try(close(shared), silent = TRUE)
  }
}, add = TRUE)

phase <- if (grepl("-warm$", mode)) "warm" else "cold"
staging <- if (grepl("response", mode)) "response" else "relation"
is_parallel <- !identical(mode, "sequential-response")
relation_matrix <- NULL
staging_seconds <- 0
pool_seconds <- 0
shared_bytes <- 0
shared_path <- NA_character_

if (is_parallel && phase == "warm") {
  stage_start <- proc.time()[["elapsed"]]
  stage_value <- if (identical(staging, "response")) {
    responses
  } else {
    relation_matrix <- materialize_relations()
    relation_matrix
  }
  shared <- shard::share(stage_value, backing = "mmap", readonly = TRUE)
  info <- shard::shared_info(shared)
  shared_bytes <- info$size
  shared_path <- info$path
  staging_seconds <- proc.time()[["elapsed"]] - stage_start
  pool_start <- proc.time()[["elapsed"]]
  pool <- shard::pool_create(scenario$workers, rss_limit = "4GB",
    rss_drift_threshold = Inf, packages = "effectagram")
  pool_seconds <- proc.time()[["elapsed"]] - pool_start
}

gc()
rss_before <- rss_self()
writeLines("ready", ready_path, useBytes = TRUE)
started <- proc.time()[["elapsed"]]

if (!is_parallel) {
  tasks <- lapply(feature_blocks, task_from_response, source = responses)
  dispatch_result <- NULL
} else {
  if (phase == "cold") {
    stage_start <- proc.time()[["elapsed"]]
    stage_value <- if (identical(staging, "response")) {
      responses
    } else {
      relation_matrix <- materialize_relations()
      relation_matrix
    }
    shared <- shard::share(stage_value, backing = "mmap", readonly = TRUE)
    info <- shard::shared_info(shared)
    shared_bytes <- info$size
    shared_path <- info$path
    staging_seconds <- proc.time()[["elapsed"]] - stage_start
    pool_start <- proc.time()[["elapsed"]]
    pool <- shard::pool_create(scenario$workers, rss_limit = "4GB",
      rss_drift_threshold = Inf, packages = "effectagram")
    pool_seconds <- proc.time()[["elapsed"]] - pool_start
  }
  blocks <- shard::shards_list(feature_blocks)
  if (identical(staging, "response")) {
    worker_fun <- function(shard, source, extractors, partitions, effects, over,
                           observations) {
      partition_ids <- unname(as.character(partitions))
      effect_ids <- unname(as.character(effects))
      relations <- stats::setNames(lapply(seq_along(partition_ids), function(index) {
        source_rows <- ((index - 1L) * observations + 1L):(index * observations)
        extractors[[index]] %*%
          source[source_rows, shard$idx, drop = FALSE]
      }), partition_ids)
      effectagram:::.crossgram_feature_task(
        relations, shard$idx, effect_ids, partition_ids, over
      )
    }
    dispatch_result <- shard::shard_map(
      blocks, worker_fun,
      borrow = list(
        source = shared, extractors = extractors, partitions = partitions,
        effects = effects, over = over, observations = n
      ),
      workers = scenario$workers, chunk_size = 1L, autotune = FALSE,
      diagnostics = TRUE, packages = "effectagram", recycle = FALSE
    )
  } else {
    worker_fun <- function(shard, source, partitions, effects, over,
                           effect_count) {
      partition_ids <- unname(as.character(partitions))
      effect_ids <- unname(as.character(effects))
      relations <- stats::setNames(lapply(seq_along(partition_ids), function(index) {
        source_rows <- ((index - 1L) * effect_count + 1L):(index * effect_count)
        source[source_rows, shard$idx, drop = FALSE]
      }), partition_ids)
      effectagram:::.crossgram_feature_task(
        relations, shard$idx, effect_ids, partition_ids, over
      )
    }
    dispatch_result <- shard::shard_map(
      blocks, worker_fun,
      borrow = list(
        source = shared, partitions = partitions, effects = effects,
        over = over, effect_count = q
      ),
      workers = scenario$workers, chunk_size = 1L, autotune = FALSE,
      diagnostics = TRUE, packages = "effectagram", recycle = FALSE
    )
  }
  if (!shard::succeeded(dispatch_result)) {
    stop("shard dispatch failed")
  }
  tasks <- shard::results(dispatch_result)
}

geometry <- reduce_tasks(tasks)
elapsed <- proc.time()[["elapsed"]] - started
rss_after <- rss_self()
geometry_digest <- digest::digest(
  list(total = geometry$total, coherent = geometry$coherent),
  algo = "sha256", serialize = TRUE
)
task_payload_bytes <- sum(vapply(
  tasks, function(task) as.double(utils::object.size(task)), numeric(1)
))
worker_memory <- if (is_parallel) shard::mem_report(pool) else NULL

if (!is.null(pool)) {
  shard::pool_stop(pool)
  pool <- NULL
  cleanup$pool_stopped <- is.null(shard::pool_get())
}
if (!is.null(shared)) {
  close(shared)
  shared <- NULL
  cleanup$source_closed <- TRUE
  cleanup$backing_removed <- !file.exists(shared_path)
}
gc()
cleanup$rss_after_cleanup <- rss_self()

saveRDS(list(
  schema_version = 1L,
  mode = mode,
  phase = phase,
  staging = staging,
  parallel = is_parallel,
  scenario = scenario,
  elapsed_seconds = elapsed,
  staging_seconds = staging_seconds,
  pool_seconds = pool_seconds,
  rss_before_bytes = rss_before,
  rss_after_bytes = rss_after,
  worker_total_rss_at_end_bytes = if (is.null(worker_memory)) 0 else
    worker_memory$total_rss,
  worker_peak_single_rss_at_end_bytes = if (is.null(worker_memory)) 0 else
    worker_memory$peak_rss,
  shared_bytes = shared_bytes,
  task_payload_bytes = task_payload_bytes,
  geometry_digest = geometry_digest,
  geometry = list(total = geometry$total, coherent = geometry$coherent),
  task_diagnostics = geometry$diagnostics,
  shard_diagnostics = if (is.null(dispatch_result)) NULL else
    dispatch_result$diagnostics,
  cleanup = cleanup
), result_path, version = 3)
