# Complete geometry and query-only result contracts -------------------------

.memory_geometry_store <- function(value) {
  if (!is.matrix(value) || !is.numeric(value) || any(!is.finite(value))) {
    stop("Geometry components must be finite numeric matrices.", call. = FALSE)
  }
  structure(
    list(
      dim = dim(value),
      representation = "memory",
      manifest = list(schema_version = 1L, complete = TRUE, dim = dim(value),
        format = "packed-double-v1"),
      read = function(rows = NULL) {
        if (is.null(rows)) value else value[rows, , drop = FALSE]
      }
    ),
    class = "effect_geometry_store"
  )
}

.block_geometry_store <- function(dim, read) {
  if (!is.numeric(dim) || length(dim) != 2L || any(!is.finite(dim)) ||
      any(dim < 0) || any(dim %% 1 != 0)) {
    stop("A block store requires two nonnegative integer dimensions.", call. = FALSE)
  }
  if (!is.function(read)) {
    stop("A block store requires a `read` function.", call. = FALSE)
  }
  structure(
    list(
      dim = as.integer(dim),
      representation = "block_backed",
      manifest = list(schema_version = 1L, complete = TRUE,
        dim = as.integer(dim), format = "packed-double-v1"),
      read = read
    ),
    class = "effect_geometry_store"
  )
}

.as_geometry_store <- function(x) {
  if (inherits(x, "effect_geometry_store")) x else .memory_geometry_store(x)
}

.read_geometry_store <- function(store, rows = NULL) {
  if (!inherits(store, "effect_geometry_store")) {
    stop("Invalid geometry store.", call. = FALSE)
  }
  value <- store$read(rows)
  expected_rows <- if (is.null(rows)) store$dim[[1L]] else length(rows)
  if (!is.matrix(value) || !is.numeric(value) || any(!is.finite(value)) ||
      !identical(dim(value), c(as.integer(expected_rows), store$dim[[2L]]))) {
    stop("Geometry store returned a block with invalid shape or values.",
      call. = FALSE)
  }
  value
}

#' Construct a semantically complete effect geometry
#'
#' An `effect_geometry` contains complete packed total and coherent geometry.
#' Its storage may be in memory or block-backed, but storage never changes the
#' result's completeness. Configuration geometry is derived exactly as total
#' minus coherent.
#'
#' @param total,coherent Packed geometry matrices, with measurements in rows,
#'   or internal geometry stores having the same dimensions.
#' @param marginals Pairing-appropriate signed marginals. Undirected pairings
#'   contain `endpoint`; directed pairings contain `left` and `right`.
#' @param effects An `effect_space()` whose dimension must match the triangular
#'   packed-geometry width and marginal columns. Unique names are accepted as
#'   shorthand for an unspecified-basis space.
#' @param index Optional measurement index with one entry per geometry row.
#' @param receipt The `execution_receipt()` proving how the result was made.
#' @param metadata Optional compact semantic metadata.
#' @return A complete `effect_geometry`.
#' @export
effect_geometry <- function(total, coherent, marginals, effects, receipt, index = NULL,
                            metadata = list()) {
  total <- .as_geometry_store(total)
  coherent <- .as_geometry_store(coherent)
  if (is.null(index)) {
    index <- seq_len(total$dim[[1L]])
  }
  validated <- .validate_complete_geometry(total, coherent, marginals, effects,
    index, receipt, metadata)
  metadata$scientific_plan_id <- receipt$scientific_plan_id

  structure(
    list(
      total = total,
      coherent = coherent,
      marginals = marginals,
      effect_space = validated$effect_space,
      effects = validated$effect_space$coordinates,
      index = index,
      metadata = metadata,
      receipt = receipt,
      completeness = "full",
      storage = unique(c(total$representation, coherent$representation))
    ),
    class = "effect_geometry"
  )
}

.validate_geometry_store <- function(store, label) {
  if (!inherits(store, "effect_geometry_store") ||
      !is.numeric(store$dim) || length(store$dim) != 2L ||
      anyNA(store$dim) || any(!is.finite(store$dim)) || any(store$dim < 1L) ||
      any(store$dim %% 1 != 0) || !is.function(store$read)) {
    stop(sprintf("`%s` is not a valid complete geometry store.", label), call. = FALSE)
  }
  manifest <- store$manifest
  if (!is.list(manifest) || !identical(manifest$schema_version, 1L) ||
      !isTRUE(manifest$complete) || !identical(manifest$dim, store$dim) ||
      !identical(manifest$format, "packed-double-v1")) {
    stop(sprintf("`%s` has an incomplete or inconsistent store manifest.", label),
      call. = FALSE)
  }
  probes <- unique(c(1L, store$dim[[1L]]))
  tryCatch(
    .read_geometry_store(store, probes),
    error = function(error) stop(sprintf("`%s` reader cannot supply claimed geometry: %s",
      label, conditionMessage(error)), call. = FALSE)
  )
  invisible(store)
}

.validate_complete_geometry <- function(total, coherent, marginals, effects,
                                        index, receipt, metadata) {
  .validate_geometry_store(total, "total")
  .validate_geometry_store(coherent, "coherent")
  if (!identical(total$dim, coherent$dim)) {
    stop("`total` and `coherent` must have identical dimensions.", call. = FALSE)
  }
  packed_width <- total$dim[[2L]]
  effect_count <- (sqrt(8 * packed_width + 1) - 1) / 2
  if (!is.finite(effect_count) || effect_count %% 1 != 0) {
    stop("Packed geometry width is not triangular.", call. = FALSE)
  }
  effects <- .as_effect_space(effects, effect_count)
  effect_names <- effects$coordinates
  if (!inherits(marginals, "effect_marginals") || length(marginals) < 1L) {
    stop("`marginals` must be a nonempty pairing-appropriate marginal object.",
      call. = FALSE)
  }
  semantics <- attr(marginals, "semantics")
  expected_names <- switch(semantics,
    undirected_endpoint = "endpoint",
    directed_roles = c("left", "right"),
    stop("Marginal semantics are missing or invalid.", call. = FALSE)
  )
  if (!identical(names(marginals), expected_names)) {
    stop("Marginal members do not match their declared endpoint semantics.",
      call. = FALSE)
  }
  if (!all(vapply(marginals, function(x) {
    is.matrix(x) && is.numeric(x) && identical(dim(x),
      c(total$dim[[1L]], as.integer(effect_count))) && all(is.finite(x)) &&
      identical(colnames(x), effect_names)
  }, logical(1)))) {
    stop("Every marginal must match measurement rows and named effect columns.",
      call. = FALSE)
  }
  if (length(index) != total$dim[[1L]] || anyNA(index) || anyDuplicated(index)) {
    stop("`index` must uniquely identify every measurement row.", call. = FALSE)
  }
  if (!is.list(metadata)) stop("`metadata` must be a list.", call. = FALSE)
  .validate_execution_receipt(receipt)
  if (!is.null(metadata$scientific_plan_id) &&
      !identical(metadata$scientific_plan_id, receipt$scientific_plan_id)) {
    stop("Result metadata and execution receipt identify different scientific plans.",
      call. = FALSE)
  }
  list(effect_space = effects)
}

#' Construct an explicit query-only effect view
#'
#' An `effect_view` contains derived values but not the geometry from which they
#' arose. It is never returned where a complete `effect_geometry` is promised.
#'
#' @param values A finite numeric measurement-by-view matrix.
#' @param query The compiled query matrix or a descriptive query value.
#' @param component Geometry component used to create the view.
#' @param index Optional measurement index.
#' @param receipt The `execution_receipt()` proving how the view was made.
#' @param metadata Optional compact semantic metadata.
#' @param effects The `effect_space()` to which the view coordinates refer.
#' @return A query-only `effect_view`.
#' @export
effect_view <- function(values, query, component, receipt, index = NULL,
                        metadata = list(), effects = NULL) {
  if (!is.matrix(values) || !is.numeric(values) || any(!is.finite(values))) {
    stop("`values` must be a finite numeric matrix.", call. = FALSE)
  }
  component <- match.arg(component, c("total", "coherent", "configuration"))
  if (is.null(index)) index <- seq_len(nrow(values))
  if (length(index) != nrow(values)) {
    stop("`index` must have one entry per measurement.", call. = FALSE)
  }
  if (!is.list(metadata)) {
    stop("`metadata` must be a list.", call. = FALSE)
  }
  .validate_execution_receipt(receipt)
  if (is.null(effects) && inherits(query, "effect_query") &&
      !is.null(query$effect_space)) {
    effects <- query$effect_space
  }
  if (is.null(effects)) {
    packed_width <- if (inherits(query, "effect_query")) {
      nrow(query$operator) * (nrow(query$operator) + 1L) / 2L
    } else if (is.matrix(query)) {
      nrow(query)
    } else {
      NA_real_
    }
    effect_count <- (sqrt(8 * packed_width + 1) - 1) / 2
    if (!is.finite(effect_count) || effect_count < 1L || effect_count %% 1 != 0) {
      stop("`effects` is required when the query does not identify a packed effect space.",
        call. = FALSE)
    }
    effects <- effect_space(paste0("effect", seq_len(effect_count)))
  } else {
    effects <- .validate_effect_space(effects)
  }

  structure(
    list(
      values = values,
      query = query,
      effect_space = effects,
      component = component,
      index = index,
      metadata = metadata,
      receipt = receipt,
      completeness = "query_only"
    ),
    class = "effect_view"
  )
}

#' Read one component of a complete geometry
#'
#' @param x An `effect_geometry`.
#' @param component One of `total`, `coherent`, or `configuration`.
#' @param rows Optional measurement rows to read.
#' @return A packed numeric geometry matrix.
#' @export
geometry_component <- function(x, component = "total", rows = NULL) {
  if (!inherits(x, "effect_geometry")) {
    stop("`x` must be a complete effect_geometry, not a query-only view.",
      call. = FALSE)
  }
  component <- match.arg(component, c("total", "coherent", "configuration"))
  if (!is.null(rows) &&
      (!is.numeric(rows) || anyNA(rows) || any(rows %% 1 != 0) ||
       any(rows < 1) || any(rows > x$total$dim[[1L]]))) {
    stop("`rows` contains invalid measurement indices.", call. = FALSE)
  }
  total <- if (component != "coherent") .read_geometry_store(x$total, rows) else NULL
  coherent <- if (component != "total") .read_geometry_store(x$coherent, rows) else NULL

  switch(component,
    total = total,
    coherent = coherent,
    configuration = total - coherent
  )
}

#' Apply a linear query to a complete geometry
#'
#' @param x An `effect_geometry`.
#' @param query A finite packed-coordinate-by-view numeric matrix, or a fixed
#'   symmetric `bilinear_query()` whose operator matches the effect space.
#' @param component Geometry component to query.
#' @param row_block Positive number of measurement rows read at once. This
#'   bounds packed-geometry memory for block-backed stores.
#' @return An `effect_view`, not another geometry.
#' @export
query_geometry <- function(x, query, component = "total", row_block = 1024L) {
  validated <- .validate_geometry_query(x, query, component, row_block)
  query <- validated$query
  component <- validated$component
  row_block <- validated$row_block

  values <- matrix(0, x$total$dim[[1L]], ncol(query))
  colnames(values) <- colnames(query)
  for (start in .tile_starts(x$total$dim[[1L]], row_block)) {
    rows <- start:min(start + row_block - 1L, x$total$dim[[1L]])
    values[rows, ] <- geometry_component(x, component, rows = rows) %*% query
  }
  effect_view(
    values = values,
    query = query,
    component = component,
    receipt = x$receipt,
    index = x$index,
    metadata = x$metadata
    , effects = x$effect_space
  )
}

.validate_geometry_query <- function(x, query, component, row_block) {
  if (!inherits(x, "effect_geometry") || !identical(x$completeness, "full") ||
      !inherits(x$total, "effect_geometry_store") ||
      !is.numeric(x$total$dim) || length(x$total$dim) != 2L ||
      any(!is.finite(x$total$dim)) || any(x$total$dim < 1L)) {
    stop("`x` must be a structurally complete effect_geometry.", call. = FALSE)
  }
  effect_space <- .validate_effect_space(x$effect_space)
  if (!identical(x$effects, effect_space$coordinates)) {
    stop("Geometry coordinate labels are inconsistent with its effect space.",
      call. = FALSE)
  }
  component <- match.arg(component, c("total", "coherent", "configuration"))
  row_block <- .validate_tile_size(row_block, "row_block")

  if (inherits(query, "effect_query")) {
    .validate_query_for_compile(query)
    if (!identical(query$kind, "bilinear") || !isTRUE(query$fixed)) {
      stop("Geometry operator queries must be fixed bilinear queries.", call. = FALSE)
    }
    if (nrow(query$operator) != length(x$effect_space$coordinates)) {
      stop("The query operator dimension must equal the experimental dimension.",
        call. = FALSE)
    }
    if (!is.null(query$effect_space) &&
        !.same_effect_space(query$effect_space, x$effect_space)) {
      stop("The query and geometry effect spaces are incompatible.",
        call. = FALSE)
    }
    query <- matrix(.svec_symmetric(query$operator), ncol = 1L,
      dimnames = list(NULL, "view1"))
  }
  if (!is.matrix(query) || !is.numeric(query) ||
      nrow(query) < 1L || ncol(query) < 1L || any(!is.finite(query))) {
    stop("`query` must be a finite, nonempty packed-coordinate matrix.",
      call. = FALSE)
  }
  if (nrow(query) != x$total$dim[[2L]]) {
    stop("The query input dimension must equal the packed geometry width.",
      call. = FALSE)
  }
  if (is.null(colnames(query))) colnames(query) <- paste0("view", seq_len(ncol(query)))
  list(query = query, component = component, row_block = row_block)
}

.svec_symmetric <- function(x) {
  q <- nrow(x)
  out <- numeric(q * (q + 1L) / 2L)
  k <- 0L
  for (column in seq_len(q)) {
    for (row in column:q) {
      k <- k + 1L
      out[[k]] <- x[row, column] * if (row == column) 1 else sqrt(2)
    }
  }
  out
}
