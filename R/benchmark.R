# Reproducible memory benchmark evidence ------------------------------------

.memory_benchmark_scenarios <- function() {
  data.frame(
    id = c(
      "small-dense-memory-cold", "medium-sparse-memory-cold",
      "medium-sparse-block-cold", "medium-sparse-memory-warm"
    ),
    features = c(128L, 512L, 512L, 512L),
    effects = c(8L, 12L, 12L, 12L),
    measurements = c(64L, 256L, 256L, 256L),
    partitions = c(3L, 3L, 3L, 3L),
    density = c(1, 0.05, 0.05, 0.05),
    row_tile = c(16L, 32L, 32L, 32L),
    coordinate_tile = c(12L, 24L, 24L, 24L),
    feature_tile = c(32L, 64L, 64L, 64L),
    storage = c("memory", "memory", "block", "memory"),
    phase = c("cold", "cold", "cold", "warm"),
    stringsAsFactors = FALSE
  )
}

.allocation_summary <- function(path) {
  lines <- readLines(path, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", lines)))
  bytes <- bytes[is.finite(bytes)]
  list(
    allocation_count = length(bytes),
    allocated_bytes = sum(bytes),
    largest_allocation_bytes = if (length(bytes)) max(bytes) else 0
  )
}

.current_rss_bytes <- function() {
  if (requireNamespace("ps", quietly = TRUE)) {
    value <- tryCatch(unname(ps::ps_memory_info()[["rss"]]),
      error = function(error) NA_real_)
    if (is.finite(value)) return(value)
  }
  lines <- tryCatch(
    system2("ps", c("-o", "rss=", "-p", as.character(Sys.getpid())),
      stdout = TRUE, stderr = FALSE),
    error = function(error) character()
  )
  value <- suppressWarnings(as.numeric(trimws((lines %||% NA_character_)[[1L]])))
  if (is.finite(value)) value * 1024 else NA_real_
}

`%||%` <- function(x, y) if (length(x) && !is.na(x[[1L]])) x else y

.run_memory_benchmark_case <- function(scenario, seed = 20260812L,
                                       ready_path = NULL) {
  if (is.data.frame(scenario)) scenario <- as.list(scenario[1L, , drop = FALSE])
  required <- c("id", "features", "effects", "measurements", "partitions",
    "density", "row_tile", "coordinate_tile", "feature_tile", "storage", "phase")
  if (!is.list(scenario) || !all(required %in% names(scenario))) {
    .input_error("Invalid memory benchmark scenario.")
  }
  set.seed(seed)
  p <- as.integer(scenario$features)
  q <- as.integer(scenario$effects)
  m <- as.integer(scenario$measurements)
  r <- as.integer(scenario$partitions)
  h <- q * (q + 1L) / 2L
  weights <- matrix(stats::runif(m * p), m, p)
  weights[matrix(stats::runif(m * p) > scenario$density, m, p)] <- 0
  weights[cbind(seq_len(m), ((seq_len(m) - 1L) %% p) + 1L)] <- 1
  weights <- Matrix::Matrix(weights, sparse = TRUE)
  frame <- additive_frame(weights)
  effects <- paste0("effect", seq_len(q))
  partitions <- paste0("partition", seq_len(r))
  relations <- stats::setNames(lapply(partitions, function(partition) {
    matrix(stats::rnorm(q * p), q, p)
  }), partitions)
  over <- cross_partitions(partitions, independence = "independent")
  reader <- function(partition, features) {
    relations[[partition]][, features, drop = FALSE]
  }
  allocation_log <- tempfile(fileext = ".Rprofmem")
  on.exit(unlink(allocation_log), add = TRUE)

  run_kernel <- function(total_store = NULL, coherent_store = NULL) {
    total_writer <- if (is.null(total_store)) NULL else
      function(rows, coordinates, increment) {
        existing <- .read_geometry_tile(total_store, rows, coordinates)
        .write_geometry_tile(total_store, rows, coordinates, existing + increment)
      }
    total <- .streamed_crossgram_contraction(
      frame, reader, partitions, effects, over,
      feature_block = scenario$feature_tile,
      row_tile = scenario$row_tile,
      coordinate_tile = scenario$coordinate_tile,
      accumulate_tile = total_writer,
      retain_local_relations = TRUE
    )
    coherent <- .coherent_geometry_from_local(
      total$local_relations, over, Matrix::rowSums(weights),
      row_tile = scenario$row_tile,
      write_tile = if (is.null(coherent_store)) NULL else
        function(rows, coordinates, value) {
          .write_geometry_tile(coherent_store, rows, coordinates, value)
        }
    )
    list(total = total, coherent = coherent)
  }

  if (identical(scenario$phase, "warm")) {
    warm <- run_kernel()
    rm(warm)
    gc()
  }

  total_path <- coherent_path <- NULL
  total_store <- coherent_store <- NULL
  if (identical(scenario$storage, "block")) {
    total_path <- tempfile(fileext = ".egm")
    coherent_path <- tempfile(fileext = ".egm")
    on.exit(unlink(c(total_path, coherent_path)), add = TRUE)
    total_store <- .file_geometry_store(total_path, c(m, h), create = TRUE)
    coherent_store <- .file_geometry_store(coherent_path, c(m, h), create = TRUE)
  }
  gc()
  rss_before <- .current_rss_bytes()
  if (!is.null(ready_path)) writeLines("ready", ready_path, useBytes = TRUE)
  utils::Rprofmem(allocation_log)
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(
    run_kernel(total_store, coherent_store),
    finally = utils::Rprofmem(NULL)
  )
  elapsed <- proc.time()[["elapsed"]] - started
  rss_after <- .current_rss_bytes()
  allocations <- .allocation_summary(allocation_log)
  total_diagnostics <- result$total$diagnostics
  coherent_diagnostics <- result$coherent$diagnostics
  output_bytes <- total_diagnostics$durable_output_bytes +
    total_diagnostics$durable_local_relation_bytes +
    coherent_diagnostics$durable_output_bytes
  total_contraction <- total_diagnostics$max_weight_slice_bytes +
    total_diagnostics$max_atom_slice_bytes +
    total_diagnostics$max_product_bytes
  local_contraction <- total_diagnostics$max_weight_slice_bytes +
    total_diagnostics$max_local_product_bytes
  coherent_contraction <- coherent_diagnostics$max_tile_bytes +
    coherent_diagnostics$max_work_bytes
  external_tile <- if (is.null(total_store)) 0 else as.double(utils::object.size(
    matrix(0, min(m, scenario$row_tile), min(h, scenario$coordinate_tile))
  ))
  plan <- memory_plan(
    frame_bytes = as.numeric(utils::object.size(weights)),
    resident_source_bytes = sum(vapply(relations,
      function(x) as.numeric(utils::object.size(x)), numeric(1))),
    relation_block_bytes = total_diagnostics$max_relation_bytes,
    atom_block_bytes = total_diagnostics$max_atom_bytes +
      total_diagnostics$max_atom_work_bytes,
    local_state_bytes = total_diagnostics$durable_local_relation_bytes,
    output_bytes = output_bytes - total_diagnostics$durable_local_relation_bytes,
    contraction_bytes = max(total_contraction, local_contraction,
      coherent_contraction),
    replacement_copy_bytes = max(
      total_diagnostics$max_existing_slice_bytes +
        total_diagnostics$max_replacement_bytes,
      total_diagnostics$max_local_existing_bytes +
        total_diagnostics$max_local_replacement_bytes,
      2 * external_tile
    ),
    safety_factor = 1.25
  )

  list(
    schema_version = 3L,
    scenario = scenario,
    elapsed_seconds = elapsed,
    rss_before_bytes = rss_before,
    rss_after_bytes = rss_after,
    baseline_rss_bytes = rss_before,
    sampled_absolute_peak_rss_bytes = suppressWarnings(max(c(rss_before, rss_after),
      na.rm = TRUE)),
    sampled_incremental_peak_rss_bytes = suppressWarnings(max(0,
      rss_after - rss_before, na.rm = TRUE)),
    os_peak_rss_bytes = NA_real_,
    worker_peak_rss_bytes = 0,
    aggregate_peak_rss_bytes = NA_real_,
    allocation = allocations,
    plan = plan,
    total = total_diagnostics,
    coherent = coherent_diagnostics,
    output_file_bytes = if (is.null(total_path)) 0 else
      file.info(total_path)$size + file.info(coherent_path)$size,
    rss_evidence = "sampled_only_until_child_process_runner_attaches_os_peak"
  )
}

.parse_os_peak_rss <- function(lines) {
  gnu <- grep("Maximum resident set size \\(kbytes\\):", lines, value = TRUE)
  if (length(gnu)) {
    value <- suppressWarnings(as.numeric(sub(".*: *", "", gnu[[length(gnu)]])))
    if (is.finite(value)) return(value * 1024)
  }
  mac <- grep("maximum resident set size", lines, value = TRUE, ignore.case = TRUE)
  if (length(mac)) {
    value <- suppressWarnings(as.numeric(sub("^ *([0-9]+).*$", "\\1", mac[[length(mac)]])))
    if (is.finite(value)) return(value)
  }
  NA_real_
}
