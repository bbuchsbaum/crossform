# Explicit operation-stage identity ----------------------------------------

.new_edge_normalizer <- function(kind, zero_policy = "error") {
  structure(
    list(schema_version = 1L, kind = kind, zero_policy = zero_policy),
    class = "effect_edge_normalizer"
  )
}

.validate_edge_normalizer <- function(x) {
  if (!inherits(x, "effect_edge_normalizer") || !is.list(x) ||
      !identical(names(x), c("schema_version", "kind", "zero_policy")) ||
      !identical(x$schema_version, 1L) ||
      !x$kind %in% c("inner_product", "covariance", "cosine", "correlation") ||
      !x$zero_policy %in% c("error", "zero")) {
    stop("Invalid edge normalizer.", call. = FALSE)
  }
  x
}

#' Declare uncentered inner-product edge geometry
#'
#' @return An immutable edge-normalization specification.
#' @export
inner_product <- function() .new_edge_normalizer("inner_product")

# These constructors remain private vocabulary for validating stored research
# plans. No exported workflow accepts them, so they are deliberately absent
# from the public API until an executable, scale-qualified path exists.
covariance <- function() .new_edge_normalizer("covariance")

cosine <- function(zero_norm = c("error", "zero")) {
  .new_edge_normalizer("cosine", match.arg(zero_norm))
}

correlation <- function(zero_variance = c("error", "zero")) {
  .new_edge_normalizer("correlation", match.arg(zero_variance))
}

.identity_edge_transform <- function() {
  structure(
    list(schema_version = 1L, kind = "identity", boundary = NULL,
      delta = NULL, ties = NULL),
    class = "effect_edge_transform"
  )
}

.new_edge_transform <- function(kind, boundary = NULL, delta = NULL,
                                ties = NULL) {
  structure(
    list(schema_version = 1L, kind = kind, boundary = boundary,
      delta = delta, ties = ties),
    class = "effect_edge_transform"
  )
}

.validate_edge_transform <- function(x) {
  expected <- c("schema_version", "kind", "boundary", "delta", "ties")
  if (!inherits(x, "effect_edge_transform") || !is.list(x) ||
      !identical(names(x), expected) || !identical(x$schema_version, 1L) ||
      !x$kind %in% c("identity", "fisher_z", "rank_edges")) {
    stop("Invalid edge transform.", call. = FALSE)
  }
  if (x$kind == "identity" &&
      (!is.null(x$boundary) || !is.null(x$delta) || !is.null(x$ties))) {
    stop("Identity edge transforms cannot carry policy parameters.",
      call. = FALSE)
  }
  if (x$kind == "fisher_z") {
    if (!x$boundary %in% c("error", "clip") || !is.null(x$ties) ||
        (x$boundary == "error" && !is.null(x$delta)) ||
        (x$boundary == "clip" &&
          (!is.numeric(x$delta) || length(x$delta) != 1L || is.na(x$delta) ||
           !is.finite(x$delta) || x$delta <= 0 || x$delta >= 1))) {
      stop("Invalid Fisher boundary policy.", call. = FALSE)
    }
  }
  allowed_ties <- c("average", "first", "last", "min", "max")
  if (x$kind == "rank_edges" &&
      (!is.null(x$boundary) || !is.null(x$delta) ||
       !is.character(x$ties) || length(x$ties) != 1L || is.na(x$ties) ||
       !x$ties %in% allowed_ties)) {
    stop("Invalid edge-ranking tie policy.", call. = FALSE)
  }
  x
}

# Private constructors used only by internal research-stage views.
fisher_z <- function(boundary = c("error", "clip"), delta = NULL) {
  boundary <- match.arg(boundary)
  value <- .new_edge_transform("fisher_z", boundary = boundary, delta = delta)
  .validate_edge_transform(value)
}

rank_edges <- function(ties = c("average", "first", "last", "min", "max")) {
  value <- .new_edge_transform("rank_edges", ties = match.arg(ties))
  .validate_edge_transform(value)
}

.new_partition_reducer <- function(order) {
  structure(
    list(kind = "weighted_sum", weight_convention = "normalized_unit_mass",
      order = order),
    class = "effect_partition_reducer"
  )
}

#' Reduce normalized and transformed partition edges
#'
#' This is the default estimand: normalize each declared edge, transform it,
#' and only then apply the partition weights.
#'
#' @return An immutable partition-reducer specification.
#' @export
reduce_partitions <- function() .new_partition_reducer("edge_first")

#' Aggregate edge sufficient statistics before normalization
#'
#' This names a distinct estimand. Raw sufficient statistics are combined by
#' partition weights before normalization and transformation.
#'
#' @return An immutable partition-reducer specification.
#' @export
aggregate_first <- function() .new_partition_reducer("aggregate_first")

.edge_operation_plan <- function(normalizer = inner_product(),
                                 transform = NULL,
                                 reducer = reduce_partitions()) {
  normalizer <- .validate_edge_normalizer(normalizer)
  if (is.null(transform)) transform <- .identity_edge_transform()
  transform <- .validate_edge_transform(transform)
  reducer <- .validate_partition_reducer(reducer)
  if (transform$kind == "fisher_z" && normalizer$kind != "correlation") {
    stop("Fisher transformation requires correlation-valued edge input.",
      call. = FALSE)
  }
  lowering <- if (transform$kind != "identity") {
    "required_edge_materialization"
  } else {
    switch(normalizer$kind,
      inner_product = "bilinear",
      covariance = "centered_bilinear",
      cosine = "normalized_bilinear",
      correlation = "normalized_bilinear"
    )
  }
  semantic <- list(
    schema_version = 1L,
    operation_order = c(
      "ordered_edge_product", "spatial_normalization", "edge_transform",
      "partition_reduction", "pair_query"
    ),
    normalizer = unclass(normalizer),
    transform = unclass(transform),
    reducer = unclass(reducer),
    lowering = lowering
  )
  structure(c(semantic, list(
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE
    ))
  )), class = "effect_edge_operation_plan")
}

.apply_edge_transform <- function(value, transform) {
  if (transform$kind == "identity") return(value)
  if (transform$kind == "rank_edges") {
    return(matrix(rank(as.vector(value), ties.method = transform$ties),
      nrow(value), ncol(value)))
  }
  if (transform$boundary == "error" && any(abs(value) >= 1)) {
    stop("Fisher transformation encountered an absolute-correlation boundary.",
      call. = FALSE)
  }
  if (transform$boundary == "clip") {
    value <- pmax(-1 + transform$delta, pmin(1 - transform$delta, value))
  }
  atanh(value)
}
