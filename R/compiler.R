# Public relation-to-geometry compiler -------------------------------------

.compiler_capabilities <- function(x) {
  if (is.null(x$capabilities)) {
    stop(paste0(
      "Opaque relation sources require explicit `source_capabilities()` ",
      "before execution."
    ), call. = FALSE)
  }
  capabilities <- lapply(x$capabilities, .validate_source_capabilities)
  if (!all(vapply(capabilities, function(value) isTRUE(value$block_read),
    logical(1)))) {
    stop("All relation sources must support bounded block reads.", call. = FALSE)
  }
  capabilities
}

.compiler_query <- function(query, effect_space) {
  effect_space <- .validate_effect_space(effect_space)
  q <- length(effect_space$coordinates)
  packed_width <- q * (q + 1L) / 2L
  if (inherits(query, "effect_query")) {
    .validate_query_for_compile(query)
    if (!identical(query$kind, "bilinear") || !isTRUE(query$fixed)) {
      stop("Direct execution requires a fixed bilinear query.", call. = FALSE)
    }
    if (nrow(query$operator) != q) {
      stop("The query operator dimension must equal the effect dimension.",
        call. = FALSE)
    }
    if (!is.null(query$effect_space) &&
        !.same_effect_space(query$effect_space, effect_space)) {
      stop("The query and relation effect spaces are incompatible.",
        call. = FALSE)
    }
    value <- matrix(.svec_symmetric(query$operator), ncol = 1L)
    colnames(value) <- "view1"
    return(value)
  }
  if (!is.matrix(query) || !is.numeric(query) || nrow(query) != packed_width ||
      ncol(query) < 1L || any(!is.finite(query))) {
    stop("`query` must be a finite packed-coordinate-by-view matrix.",
      call. = FALSE)
  }
  if (is.null(colnames(query))) colnames(query) <- paste0("view", seq_len(ncol(query)))
  query
}

.validate_compiler_inputs <- function(x, at, over) {
  .validate_relation(x)
  .validate_frame_for_compile(at)
  .validate_pairing(over)
  if (!identical(at$representation, "additive_diagonal")) {
    stop("effectagram 0.1 executes only additive diagonal frames.", call. = FALSE)
  }
  if (!.same_domain_reference(at$domain, x$domain)) {
    stop("Relation and frame exact neural-domain identities must agree.",
      call. = FALSE)
  }
  if (ncol(at$weights) != x$n_features) {
    stop("The frame feature dimension must equal the relation feature dimension.",
      call. = FALSE)
  }
  if (any(!c(over$left, over$right) %in% x$partitions)) {
    stop("Every pairing endpoint must identify a relation partition.",
      call. = FALSE)
  }
  invisible(TRUE)
}

.compiler_memory_plan <- function(x, at, compute, feature_block, row_tile,
                                  coordinate_tile, output_width, storage) {
  q <- length(x$effects)
  h <- q * (q + 1L) / 2L
  p <- x$n_features
  m <- nrow(at$weights)
  r <- length(x$partitions)
  f <- min(feature_block, p)
  rows <- min(row_tile, m)
  coordinates <- min(coordinate_tile, output_width)
  max_observations <- max(vapply(x$sources, function(source) source$dim[[1L]],
    integer(1)))
  source_shared <- sum(vapply(x$sources, function(source) {
    if (identical(source$kind, "matrix")) prod(as.double(source$dim)) * 8 else 0
  }, numeric(1)))
  source_block <- max_observations * f * 8
  relation_block <- r * q * f * 8
  atom_block <- f * (h + output_width) * 8
  local_bytes <- m * q * r * 8
  # Marginals are conservatively counted as two endpoint roles; undirected
  # pairings will use only half this allowance.
  marginal_bytes <- 2 * m * q * 8
  durable_geometry <- if (storage == "memory") 2 * m * output_width * 8 else 0
  output_bytes <- local_bytes + marginal_bytes + durable_geometry
  contraction <- (
    rows * f + f * coordinates + rows * coordinates + rows * q +
      rows * h
  ) * 8
  replacement <- (2 * rows * coordinates + 2 * rows * q) * 8

  memory_plan(
    shared_source_bytes = source_shared,
    source_block_bytes = source_block,
    relation_block_bytes = relation_block,
    atom_block_bytes = atom_block,
    output_bytes = output_bytes,
    contraction_bytes = contraction,
    replacement_copy_bytes = replacement,
    workers = 1L,
    n_active = 1L,
    budget_bytes = compute$memory_bytes
  )
}

.compiler_plan_id <- function(x, at, over, materialization, query, component) {
  semantic <- list(
    source_revisions = vapply(x$capabilities, `[[`, character(1),
      "stable_revision"),
    extractors = lapply(x$extractors, function(value) value$map),
    effect_space = x$effect_space,
    partitions = x$partitions,
    domain = x$domain,
    frame = list(weights = at$weights, normalization = at$normalization,
      domain = at$domain),
    pairing = list(edges = unclass(as.data.frame(over)),
      directed = attr(over, "directed"),
      self_pairs = attr(over, "self_pairs"),
      independence = attr(over, "independence")),
    materialization = materialization,
    query = query,
    component = component
  )
  paste0("sha256:", digest::digest(semantic, algo = "sha256", serialize = TRUE))
}

.compiler_index <- function(at) {
  index <- at$index
  if (is.null(index)) return(seq_len(nrow(at$weights)))
  if (is.data.frame(index) && "measurement" %in% names(index)) {
    return(index$measurement)
  }
  index
}

.compiler_blas <- function() {
  vendor <- tryCatch(unname(extSoftVersion()[["BLAS"]]), error = function(e) NULL)
  if (is.null(vendor) || !is.character(vendor) || length(vendor) != 1L ||
      is.na(vendor) || !nzchar(vendor)) vendor <- "unknown"
  list(vendor = vendor, threads = 1L)
}

.planned_compiler_receipt <- function(x, compute, memory, plan_id,
                                      feature_block, task_count) {
  execution_receipt(
    scientific_plan_id = plan_id,
    compute = compute,
    sources = x$capabilities,
    memory = memory,
    kernel_version = "additive-crossgram-v1",
    task_partition_id = sprintf("ascending-features-%d", feature_block),
    reduction_plan_id = "ascending-feature-row-coordinate-v1",
    numeric_contract = numerical_contract(),
    completion_status = "planned",
    task_count = task_count,
    completed_task_count = 0L,
    blas = .compiler_blas(),
    domain_signature = x$domain$signature
  )
}

.compiler_reporter <- function(reporter, final_receipt) {
  force(reporter)
  function(event) {
    if (event$type %in% c("complete", "failed", "interrupted")) {
      final_receipt$value <- event$receipt
    }
    if (!is.null(reporter)) reporter(event)
    invisible(NULL)
  }
}

.compiler_storage <- function(storage, storage_path, dim) {
  storage <- match.arg(storage, c("memory", "block"))
  if (storage == "memory") {
    if (!is.null(storage_path)) {
      stop("`storage_path` is only valid for block-backed geometry.",
        call. = FALSE)
    }
    return(list(kind = storage, total = NULL, coherent = NULL,
      created = character(), created_directory = FALSE))
  }
  if (!is.character(storage_path) || length(storage_path) != 1L ||
      is.na(storage_path) || !nzchar(storage_path)) {
    stop("Block-backed geometry requires one nonempty `storage_path`.",
      call. = FALSE)
  }
  created_directory <- FALSE
  if (!dir.exists(storage_path)) {
    parent <- dirname(storage_path)
    if (!dir.exists(parent)) {
      stop("The parent of `storage_path` must already exist.", call. = FALSE)
    }
    if (!dir.create(storage_path, recursive = FALSE, showWarnings = FALSE)) {
      stop("Could not create `storage_path`.", call. = FALSE)
    }
    created_directory <- TRUE
  }
  total_path <- file.path(storage_path, "total.egm")
  coherent_path <- file.path(storage_path, "coherent.egm")
  if (file.exists(total_path) || file.exists(coherent_path)) {
    stop("Refusing to overwrite an existing geometry component.", call. = FALSE)
  }
  total <- .file_geometry_store(total_path, dim, create = TRUE)
  coherent <- tryCatch(
    .file_geometry_store(coherent_path, dim, create = TRUE),
    error = function(error) {
      unlink(total_path)
      if (created_directory) unlink(storage_path, recursive = TRUE)
      stop(error)
    }
  )
  list(kind = storage, total = total, coherent = coherent,
    created = c(total_path, coherent_path),
    created_directory = created_directory)
}

.compiler_cleanup <- function(storage, completed) {
  force(storage)
  force(completed)
  function() {
    if (!isTRUE(completed$value) && length(storage$created)) {
      unlink(storage$created)
      if (storage$created_directory) unlink(dirname(storage$created[[1L]]),
        recursive = TRUE)
    }
    invisible(NULL)
  }
}

.run_geometry_compiler <- function(x, at, over, compute, storage,
                                   storage_path, query, component, reporter) {
  compute <- .validate_compute_policy(compute)
  .validate_compiler_inputs(x, at, over)
  if (!is.null(query)) query <- .compiler_query(query, x$effect_space)
  component <- if (is.null(query)) "full" else
    match.arg(component, c("total", "coherent", "configuration"))
  if (!is.null(query) && storage != "memory") {
    stop("Direct query execution returns an in-memory view and does not create geometry stores.",
      call. = FALSE)
  }
  capabilities <- .execution_preflight(compute, function() {
    .compiler_capabilities(x)
  })$source_capabilities
  x$capabilities <- capabilities

  q <- length(x$effect_space$coordinates)
  h <- q * (q + 1L) / 2L
  output_width <- if (is.null(query)) h else ncol(query)
  feature_block <- if (is.null(compute$block_features)) {
    min(1024L, x$n_features)
  } else {
    min(as.integer(compute$block_features), x$n_features)
  }
  row_tile <- min(256L, nrow(at$weights))
  coordinate_tile <- min(64L, output_width)
  memory <- .compiler_memory_plan(x, at, compute, feature_block, row_tile,
    coordinate_tile, output_width, storage)
  if (identical(memory$fits_budget, FALSE)) {
    stop(sprintf(
      "The conservative memory plan requires %.0f bytes, exceeding the %.0f-byte budget.",
      memory$conservative_peak_bytes, compute$memory_bytes
    ), call. = FALSE)
  }
  materialization <- if (is.null(query)) "full_geometry" else "direct_query"
  plan_id <- .compiler_plan_id(x, at, over, materialization, query, component)
  task_count <- ceiling(x$n_features / feature_block)
  planned_receipt <- .planned_compiler_receipt(x, compute, memory, plan_id,
    feature_block, task_count)
  geometry_storage <- .compiler_storage(storage, storage_path,
    c(nrow(at$weights), output_width))
  completed <- new.env(parent = emptyenv())
  completed$value <- FALSE
  final_receipt <- new.env(parent = emptyenv())
  final_receipt$value <- planned_receipt

  value <- .execute_guarded(
    compute = function() {
      total_accumulator <- NULL
      if (storage == "block") {
        total_accumulator <- function(rows, coordinates, increment) {
          existing <- .read_geometry_tile(geometry_storage$total, rows, coordinates)
          .write_geometry_tile(geometry_storage$total, rows, coordinates,
            existing + increment)
        }
      }
      streamed <- .streamed_crossgram_contraction(
        frame = at,
        read_relation = function(partition, features) {
          relation_block(x, partition, features)
        },
        partitions = x$partitions,
        effects = x$effect_space$coordinates,
        over = over,
        feature_block = feature_block,
        row_tile = row_tile,
        coordinate_tile = coordinate_tile,
        accumulate_tile = total_accumulator,
        retain_local_relations = TRUE,
        query = query
      )
      coherent_writer <- if (storage == "block") {
        function(rows, coordinates, value) {
          .write_geometry_tile(geometry_storage$coherent, rows, coordinates, value)
        }
      } else {
        NULL
      }
      coherent <- .coherent_geometry_from_local(
        streamed$local_relations, over, Matrix::rowSums(at$weights),
        row_tile = row_tile, write_tile = coherent_writer, query = query
      )
      total_value <- if (storage == "memory") streamed$value else
        geometry_storage$total
      coherent_value <- if (storage == "memory") coherent$value else
        geometry_storage$coherent
      marginals <- pairing_marginals(streamed$local_relations, over,
        mass = Matrix::rowSums(at$weights))
      metadata <- list(
        frame = list(representation = at$representation,
          normalization = at$normalization, domain = at$domain),
        pairing_estimate = attr(over, "estimate"),
        storage = storage,
        diagnostics = list(total = streamed$diagnostics,
          coherent = coherent$diagnostics),
        scientific_plan_id = plan_id
      )
      result <- if (is.null(query)) {
        effect_geometry(total_value, coherent_value, marginals,
          effects = x$effect_space, receipt = planned_receipt,
          index = .compiler_index(at), metadata = metadata)
      } else {
        total_matrix <- if (storage == "memory") total_value else
          .read_geometry_store(total_value)
        coherent_matrix <- if (storage == "memory") coherent_value else
          .read_geometry_store(coherent_value)
        values <- switch(component,
          total = total_matrix,
          coherent = coherent_matrix,
          configuration = total_matrix - coherent_matrix
        )
        colnames(values) <- colnames(query)
        effect_view(values, query, component, planned_receipt,
          effects = x$effect_space,
          index = .compiler_index(at), metadata = metadata)
      }
      completed$value <- TRUE
      result
    },
    receipt = planned_receipt,
    reporter = .compiler_reporter(reporter, final_receipt),
    cleanup = .compiler_cleanup(geometry_storage, completed)
  )
  value$receipt <- final_receipt$value
  value
}

#' Compile a complete cross-generalized geometry
#'
#' `geometry()` always materializes complete total and coherent packed geometry.
#' Use [evaluate_geometry()] when only one or a few fixed linear queries are
#' required.
#'
#' @param x An `effect_relation`.
#' @param at A compiled additive `effect_frame`.
#' @param over An `effect_pairing`.
#' @param storage Either `"memory"` or `"block"`.
#' @param storage_path Durable directory for block-backed geometry.
#' @param compute A sequential `compute_policy()`.
#' @param reporter Optional nonsemantic coordinator-side event reporter.
#' @return A complete `effect_geometry`.
#' @export
geometry <- function(x, at, over, storage = c("memory", "block"),
                     storage_path = NULL, compute = compute_policy(),
                     reporter = NULL) {
  storage <- match.arg(storage)
  .run_geometry_compiler(x, at, over, compute, storage, storage_path,
    query = NULL, component = NULL, reporter = reporter)
}

#' Evaluate a fixed query without materializing complete geometry
#'
#' @param x,at,over,compute,reporter As in [geometry()].
#' @param query A fixed `bilinear_query()` or packed-coordinate query matrix.
#' @param component One of `total`, `coherent`, or `configuration`.
#' @return A query-only `effect_view`.
#' @export
evaluate_geometry <- function(x, at, over, query,
                              component = c("total", "coherent", "configuration"),
                              compute = compute_policy(), reporter = NULL) {
  component <- match.arg(component)
  .run_geometry_compiler(x, at, over, compute, "memory", NULL,
    query = query, component = component, reporter = reporter)
}
