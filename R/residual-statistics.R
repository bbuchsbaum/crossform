# Canonical support-local residual sufficient statistics -------------------

.residual_statistics_support_index <- function(at, domain) {
  if (inherits(at, "effect_support_index")) {
    index <- .validate_support_index(at)
  } else {
    .validate_frame_for_compile(at)
    if (is.null(at$support_index)) {
      stop(paste0(
        "Residual pair statistics currently require a frame with an explicit ",
        "support index, such as a compiled searchlight frame."
      ), call. = FALSE)
    }
    index <- .validate_support_index(at$support_index)
  }
  if (!.same_domain_reference(index$domain, domain)) {
    stop("The residual fit and support index use different neural domains.",
      call. = FALSE)
  }
  index
}

.residual_pair_coordinates <- function(index) {
  index <- .validate_support_index(index)
  pattern <- methods::as(index$pair_pattern, "CsparseMatrix")
  columns <- rep.int(seq_len(ncol(pattern)), diff(pattern@p))
  rows <- as.integer(pattern@i + 1L)
  if (identical(pattern@uplo, "U")) {
    pair_i <- rows
    pair_j <- as.integer(columns)
  } else {
    pair_i <- as.integer(columns)
    pair_j <- rows
  }
  if (length(pair_i) < 1L) {
    stop("The support pair graph contains no neural feature pairs.",
      call. = FALSE)
  }
  ordered <- order(pair_j, pair_i, method = "radix")
  pair_i <- pair_i[ordered]
  pair_j <- pair_j[ordered]
  if (any(pair_i > pair_j) || any(pair_i < 1L) ||
      any(pair_j > index$domain$n_features)) {
    stop("The support pair graph has invalid canonical coordinates.",
      call. = FALSE)
  }
  list(i = pair_i, j = pair_j)
}

.residual_tile_width <- function(index) {
  index <- .validate_support_index(index)
  largest <- as.double(index$cost$support_size[["max"]])
  if (!is.finite(largest) || largest < 1) {
    stop("Support sizes cannot determine a canonical residual tile width.",
      call. = FALSE)
  }
  bounded <- min(128, largest)
  as.integer(2^ceiling(log2(bounded)))
}

.residual_pair_tile_plan <- function(pair_i, pair_j, features, width) {
  if (!is.integer(pair_i) || !is.integer(pair_j) ||
      length(pair_i) < 1L || !identical(length(pair_i), length(pair_j)) ||
      any(pair_i > pair_j) || any(pair_i < 1L) || any(pair_j > features)) {
    stop("Residual pair coordinates are invalid.", call. = FALSE)
  }
  width <- .validate_tile_size(width, "width")
  tiles <- as.integer(ceiling(features / width))
  left <- as.integer((pair_i - 1L) %/% width + 1L)
  right <- as.integer((pair_j - 1L) %/% width + 1L)
  key <- (as.double(right) - 1) * tiles + left
  ordered <- order(key, pair_j, pair_i, method = "radix")
  ordered_key <- key[ordered]
  starts <- c(1L, which(diff(ordered_key) != 0) + 1L)
  ends <- c(starts[-1L] - 1L, length(ordered))
  structure(list(
    width = width,
    tiles = tiles,
    pair_order = ordered,
    group_start = as.integer(starts),
    group_end = as.integer(ends),
    group_left = left[ordered[starts]],
    group_right = right[ordered[starts]],
    tile_pairs = as.integer(length(starts))
  ), class = "effect_residual_pair_tile_plan")
}

.residual_pair_memory_plan <- function(pair_count, partitions,
                                       observations, width, tiles,
                                       workspace_bytes) {
  values <- c(pair_count, partitions, observations, width, tiles)
  if (any(!is.finite(values)) || any(values < 1) || any(values %% 1 != 0)) {
    stop("Residual pair memory dimensions must be positive whole values.",
      call. = FALSE)
  }
  if (!is.numeric(workspace_bytes) || length(workspace_bytes) != 1L ||
      is.na(workspace_bytes) || !is.finite(workspace_bytes) ||
      workspace_bytes <= 0 || workspace_bytes > 2^53) {
    stop("`workspace_bytes` must be one positive finite byte budget.",
      call. = FALSE)
  }
  # Coordinates, values for every atomic partition, canonical grouping arrays,
  # and bounded extraction scratch. The estimate deliberately charges payload
  # rather than an R-version-specific object header.
  pair_state_bytes <- as.double(pair_count) * (28 + 8 * partitions)
  residual_tile_bytes <- 8 * max(observations) * width
  cross_product_bytes <- 8 * width^2
  resident_tile_floor <- min(2, tiles)
  # Eight tile payloads conservatively cover source read, pre-residualization
  # padding, whitening, projection, returned residuals, and two active inputs.
  minimum_workspace_bytes <- ceiling(
    pair_state_bytes + 2 * cross_product_bytes + 8 * residual_tile_bytes
  )
  if (minimum_workspace_bytes > workspace_bytes) {
    stop(sprintf(
      paste0(
        "Canonical residual pair accumulation requires at least %.0f bytes, ",
        "exceeding the %.0f-byte workspace budget."
      ), minimum_workspace_bytes, workspace_bytes
    ), call. = FALSE)
  }
  extra_tiles <- floor(
    (workspace_bytes - minimum_workspace_bytes) / residual_tile_bytes
  )
  cache_capacity <- min(
    tiles, 64L, as.integer(resident_tile_floor + extra_tiles)
  )
  cache_capacity <- max(resident_tile_floor, cache_capacity)
  planned_peak_bytes <- minimum_workspace_bytes +
    max(0, cache_capacity - resident_tile_floor) * residual_tile_bytes
  structure(list(
    workspace_bytes = as.double(workspace_bytes),
    pair_state_bytes = pair_state_bytes,
    residual_tile_bytes = residual_tile_bytes,
    cross_product_bytes = cross_product_bytes,
    minimum_workspace_bytes = minimum_workspace_bytes,
    cache_capacity = as.integer(cache_capacity),
    planned_peak_bytes = as.double(planned_peak_bytes),
    fixed_shape = c(as.integer(max(observations)), width)
  ), class = "effect_residual_pair_memory_plan")
}

.residual_tile_features <- function(tile, width, features) {
  first <- (as.integer(tile) - 1L) * width + 1L
  last <- min(first + width - 1L, features)
  seq.int(first, last)
}

.residual_pair_atomic_signature <- function(partition, values, residual_df,
                                             error_model, source_revision,
                                             residual_revision) {
  .sha256_signature(list(
    schema_version = 1L,
    partition = partition,
    cross_products = values,
    residual_df = residual_df,
    error_model = error_model,
    source_revision = source_revision,
    residual_revision = residual_revision
  ))
}

.accumulate_residual_pair_partition <- function(x, partition, pair_i, pair_j,
                                                tile_plan, memory_plan) {
  features <- x$relation$n_features
  width <- tile_plan$width
  capacity <- memory_plan$cache_capacity
  cache <- new.env(hash = TRUE, parent = emptyenv(), size = max(29L, capacity))
  lru <- integer()
  reads <- 0L
  hits <- 0L
  misses <- 0L

  get_tile <- function(tile) {
    key <- as.character(tile)
    if (exists(key, envir = cache, inherits = FALSE)) {
      hits <<- hits + 1L
      lru <<- c(lru[lru != tile], tile)
      return(get(key, envir = cache, inherits = FALSE))
    }
    misses <<- misses + 1L
    if (length(lru) >= capacity) {
      remove_key <- as.character(lru[[1L]])
      rm(list = remove_key, envir = cache)
      lru <<- lru[-1L]
    }
    block_features <- .residual_tile_features(tile, width, features)
    value <- .residual_padded_block(
      x, partition, block_features, width
    )
    reads <<- reads + 1L
    assign(key, value, envir = cache)
    lru <<- c(lru, tile)
    value
  }

  values <- numeric(length(pair_i))
  for (group in seq_len(tile_plan$tile_pairs)) {
    first <- tile_plan$group_start[[group]]
    last <- tile_plan$group_end[[group]]
    positions <- tile_plan$pair_order[seq.int(first, last)]
    left_tile <- tile_plan$group_left[[group]]
    right_tile <- tile_plan$group_right[[group]]
    left <- get_tile(left_tile)
    right <- if (left_tile == right_tile) left else get_tile(right_tile)
    product <- crossprod(left, right)
    left_local <- pair_i[positions] - (left_tile - 1L) * width
    right_local <- pair_j[positions] - (right_tile - 1L) * width
    values[positions] <- product[cbind(left_local, right_local)]
  }
  model <- x$error_models[[partition]]
  df <- model$residual_df
  structure(list(
    partition = partition,
    cross_products = values,
    residual_df = df,
    error_model = model$signature,
    source_revision = model$source_revision,
    residual_revision = model$residual_source$stable_revision,
    signature = .residual_pair_atomic_signature(
      partition, values, df, model$signature, model$source_revision,
      model$residual_source$stable_revision
    )
  ), class = "effect_atomic_residual_pair_statistics",
  execution = list(
    residual_reads = reads,
    cache_hits = hits,
    cache_misses = misses,
    cache_capacity = capacity,
    tile_pairs = tile_plan$tile_pairs
  ))
}

.residual_pair_statistics_signature <- function(x) {
  .sha256_signature(list(
    schema_version = 1L,
    relation_fit = x$relation_fit,
    domain = x$domain,
    support_index = x$support_index,
    pair_i = x$pair_i,
    pair_j = x$pair_j,
    partitions = x$partitions,
    atomic = vapply(x$atomic, `[[`, character(1), "signature"),
    numerical_contract = x$numerical_contract
  ))
}

.validate_atomic_residual_pair_statistics <- function(x, partition,
                                                      pair_count) {
  expected <- c("partition", "cross_products", "residual_df", "error_model",
    "source_revision", "residual_revision", "signature")
  if (!inherits(x, "effect_atomic_residual_pair_statistics") ||
      !is.list(x) || !identical(names(x), expected) ||
      !identical(x$partition, partition) ||
      !is.numeric(x$cross_products) || length(x$cross_products) != pair_count ||
      any(!is.finite(x$cross_products)) || !is.integer(x$residual_df) ||
      length(x$residual_df) != 1L || is.na(x$residual_df) ||
      x$residual_df < 1L || !.strong_sha256(x$error_model) ||
      !.strong_sha256(x$source_revision) ||
      !.strong_sha256(x$residual_revision) || !.strong_sha256(x$signature)) {
    stop("Atomic residual pair statistics are invalid.", call. = FALSE)
  }
  expected_signature <- .residual_pair_atomic_signature(
    x$partition, x$cross_products, x$residual_df, x$error_model,
    x$source_revision, x$residual_revision
  )
  if (!identical(x$signature, expected_signature)) {
    stop("Atomic residual pair identity is inconsistent.", call. = FALSE)
  }
  x
}

.validate_residual_pair_statistics <- function(x, deep = FALSE) {
  expected <- c("relation_fit", "domain", "support_index", "pair_i",
    "pair_j", "partitions", "atomic", "numerical_contract", "execution",
    "signature")
  if (!inherits(x, "effect_residual_pair_statistics") || !is.list(x) ||
      !identical(names(x), expected) || !.strong_sha256(x$relation_fit) ||
      !.strong_sha256(x$support_index) || !is.integer(x$pair_i) ||
      !is.integer(x$pair_j) || length(x$pair_i) < 1L ||
      !identical(length(x$pair_i), length(x$pair_j)) ||
      any(x$pair_i > x$pair_j) || !is.character(x$partitions) ||
      length(x$partitions) < 1L || anyNA(x$partitions) ||
      any(!nzchar(x$partitions)) || anyDuplicated(x$partitions) ||
      !is.list(x$atomic) || !identical(names(x$atomic), x$partitions) ||
      !is.list(x$numerical_contract) || !is.list(x$execution) ||
      !.strong_sha256(x$signature)) {
    stop("Residual-pair-statistics fields are missing or noncanonical.",
      call. = FALSE)
  }
  domain <- .validate_domain_reference(x$domain)
  if (any(x$pair_i < 1L) || any(x$pair_j > domain$n_features)) {
    stop("Residual pair coordinates do not belong to their neural domain.",
      call. = FALSE)
  }
  if (length(x$pair_i) > 1L) {
    canonical <- diff(x$pair_j) > 0L |
      (diff(x$pair_j) == 0L & diff(x$pair_i) > 0L)
    if (!all(canonical)) {
      stop("Residual pair coordinates must be unique and canonical.",
        call. = FALSE)
    }
  }
  lapply(seq_along(x$partitions), function(index) {
    .validate_atomic_residual_pair_statistics(
      x$atomic[[index]], x$partitions[[index]], length(x$pair_i)
    )
  })
  contract_names <- c("algorithm", "version", "tile_width",
    "tile_width_rule", "workspace_invariant_shape", "reduction")
  if (!identical(names(x$numerical_contract), contract_names) ||
      !identical(x$numerical_contract$algorithm,
        "fixed_shape_residual_pair_gemm") ||
      !identical(x$numerical_contract$version, 1L) ||
      !is.integer(x$numerical_contract$tile_width) ||
      length(x$numerical_contract$tile_width) != 1L ||
      x$numerical_contract$tile_width < 1L ||
      !identical(x$numerical_contract$tile_width_rule,
        "next_power_of_two_of_largest_support_capped_at_128") ||
      !identical(x$numerical_contract$workspace_invariant_shape, TRUE) ||
      !identical(x$numerical_contract$reduction,
        "one_crossproduct_per_canonical_tile_pair")) {
    stop("Residual pair numerical contract is invalid.", call. = FALSE)
  }
  if (isTRUE(deep) && !identical(x$signature,
      .residual_pair_statistics_signature(x))) {
    stop("Residual-pair-statistics identity is inconsistent.", call. = FALSE)
  }
  invisible(x)
}

#' Accumulate canonical residual pair sufficient statistics
#'
#' Computes residual cross-products only for neural feature pairs that coexist
#' in at least one requested support. The GEMM tile shape is derived from the
#' support topology and is independent of `workspace_bytes`; the budget changes
#' only how many already-computed residual tiles may be cached. Consequently,
#' changing the workspace plan cannot change the accumulated floating-point
#' values.
#'
#' @param x An `effect_relation_fit` with residual-block capability.
#' @param at A compiled frame carrying an explicit support index, such as a
#'   compiled searchlight frame. An internal `effect_support_index` is also
#'   accepted.
#' @param partitions Optional relation partitions. They are canonicalized to
#'   relation order.
#' @param workspace_bytes Positive crossform-owned workspace budget.
#' @return An `effect_residual_pair_statistics` object. `$pair_i`/`$pair_j`
#'   list the coexisting feature pairs, `$partitions` names the partitions,
#'   `$atomic` holds one `$cross_products` vector and `$residual_df` per
#'   partition, and `$numerical_contract` records the fixed tile shape.
#'   `$execution` diagnostics are excluded from the scientific `$signature`.
#' @seealso [plan_crossnobis()], which compiles these statistics into a
#'   learned metric schedule, and [rdm_sampling_covariance()], which can reuse
#'   them through `residual_strategy = "shared_pair_statistics"`.
#' @family sampling uncertainty
#' @examples
#' # Residual cross-products are accumulated only for feature pairs that
#' # coexist in some searchlight, so the cost tracks support topology rather
#' # than the square of the feature count.
#' example <- example_fmri_effects()
#' statistics <- residual_pair_statistics(example$fit, example$frame)
#' length(statistics$pair_i)
#' statistics$partitions
#'
#' # One atomic record per partition, each carrying its own residual df.
#' statistics$atomic[["run1"]]$residual_df
#'
#' # The workspace budget is a cache size, not part of the numerical shape,
#' # so shrinking it cannot change the accumulated values.
#' frugal <- residual_pair_statistics(
#'   example$fit, example$frame, workspace_bytes = 4 * 1024^2
#' )
#' identical(frugal$signature, statistics$signature)
#' @export
residual_pair_statistics <- function(
    x, at, partitions = NULL, workspace_bytes = 512 * 1024^2) {
  if (inherits(x, "effect_relation")) {
    .require_relation_fit_capability(x, "learned_metric_input")
  }
  .validate_relation_fit(x, deep = FALSE)
  index <- .residual_statistics_support_index(at, x$relation$domain)
  available <- x$relation$partitions
  if (is.null(partitions)) partitions <- available
  if (!is.character(partitions) || length(partitions) < 1L ||
      anyNA(partitions) || any(!nzchar(partitions)) ||
      anyDuplicated(partitions) || any(!partitions %in% available)) {
    stop("`partitions` must uniquely identify fitted relation partitions.",
      call. = FALSE)
  }
  partitions <- available[available %in% partitions]
  .require_relation_fit_capability(x, "learned_metric_input", partitions)
  pairs <- .residual_pair_coordinates(index)
  width <- .residual_tile_width(index)
  tile_plan <- .residual_pair_tile_plan(
    pairs$i, pairs$j, x$relation$n_features, width
  )
  observations <- vapply(partitions, function(partition) {
    x$error_models[[partition]]$residual_source$dim[[1L]]
  }, integer(1))
  memory <- .residual_pair_memory_plan(
    length(pairs$i), length(partitions), observations, width,
    tile_plan$tiles, workspace_bytes
  )
  atomic <- lapply(partitions, function(partition) {
    .accumulate_residual_pair_partition(
      x, partition, pairs$i, pairs$j, tile_plan, memory
    )
  })
  names(atomic) <- partitions
  execution_atomic <- lapply(atomic, attr, which = "execution", exact = TRUE)
  atomic <- lapply(atomic, function(value) {
    attr(value, "execution") <- NULL
    value
  })
  names(atomic) <- partitions
  numerical_contract <- list(
    algorithm = "fixed_shape_residual_pair_gemm",
    version = 1L,
    tile_width = width,
    tile_width_rule =
      "next_power_of_two_of_largest_support_capped_at_128",
    workspace_invariant_shape = TRUE,
    reduction = "one_crossproduct_per_canonical_tile_pair"
  )
  value <- structure(list(
    relation_fit = x$signature,
    domain = .domain_reference(x$relation$domain),
    support_index = index$signature,
    pair_i = pairs$i,
    pair_j = pairs$j,
    partitions = partitions,
    atomic = atomic,
    numerical_contract = numerical_contract,
    execution = list(
      memory = memory,
      atomic = execution_atomic
    ),
    signature = NA_character_
  ), class = "effect_residual_pair_statistics")
  value$signature <- .residual_pair_statistics_signature(value)
  .validate_residual_pair_statistics(value, deep = TRUE)
}

.residual_pair_scope <- function(x, partitions = x$partitions) {
  .validate_residual_pair_statistics(x)
  if (!is.character(partitions) || length(partitions) < 1L ||
      anyNA(partitions) || any(!nzchar(partitions)) ||
      anyDuplicated(partitions) || any(!partitions %in% x$partitions)) {
    stop("A residual training scope must select unique atomic partitions.",
      call. = FALSE)
  }
  partitions <- x$partitions[x$partitions %in% partitions]
  values <- lapply(x$atomic[partitions], `[[`, "cross_products")
  combined <- .canonical_reduce(values, partitions)
  df <- sum(vapply(x$atomic[partitions], `[[`, integer(1), "residual_df"))
  structure(list(
    partitions = partitions,
    cross_products = combined,
    residual_df = as.integer(df),
    covariance = combined / df,
    signature = .sha256_signature(list(
      schema_version = 1L,
      parent = x$signature,
      partitions = partitions,
      residual_df = df
    ))
  ), class = "effect_residual_pair_scope")
}

.residual_pair_scope_matrix <- function(x, partitions = x$partitions) {
  scope <- .residual_pair_scope(x, partitions)
  Matrix::sparseMatrix(
    i = x$pair_i,
    j = x$pair_j,
    x = scope$covariance,
    dims = rep(x$domain$n_features, 2L),
    symmetric = TRUE,
    giveCsparse = TRUE
  )
}
