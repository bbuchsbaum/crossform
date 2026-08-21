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
    .input_error("Invalid edge normalizer.")
  }
  x
}

#' Declare uncentered inner-product edge geometry
#'
#' `inner_product()` names the default edge normalization: raw uncentered
#' products, with no centering, scaling, or correlation normalization applied
#' to an edge before it is reduced. Use it wherever an operation stage must
#' record that the geometry is the plain bilinear form.
#'
#' @return An `effect_edge_normalizer` recording `$kind` and its
#'   `$zero_policy`. It is a declaration, not a computation.
#' @seealso [reduce_partitions()] and [aggregate_first()] for the partition
#'   stage that follows normalization; [connectivity()] for the normalized
#'   views that are gated behind explicit capabilities.
#' @examples
#' # The normalizer only names the stage; nothing is computed here.
#' normalizer <- crossform:::inner_product()
#' normalizer$kind
#'
#' # Because no denominator is involved, there is no zero-norm case to
#' # resolve at execution time.
#' normalizer$zero_policy
#' @keywords internal
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
  if (!.sealed_fields(x, "effect_edge_transform", expected) ||
      !identical(x$schema_version, 1L) ||
      !x$kind %in% c("identity", "fisher_z", "rank_edges")) {
    .input_error("Invalid edge transform.")
  }
  if (x$kind == "identity" &&
      (!is.null(x$boundary) || !is.null(x$delta) || !is.null(x$ties))) {
    .input_error("Identity edge transforms cannot carry policy parameters.")
  }
  if (x$kind == "fisher_z") {
    if (!x$boundary %in% c("error", "clip") || !is.null(x$ties) ||
        (x$boundary == "error" && !is.null(x$delta)) ||
        (x$boundary == "clip" &&
          (!is.numeric(x$delta) || length(x$delta) != 1L || is.na(x$delta) ||
           !is.finite(x$delta) || x$delta <= 0 || x$delta >= 1))) {
      .input_error("Invalid Fisher boundary policy.")
    }
  }
  allowed_ties <- c("average", "first", "last", "min", "max")
  if (x$kind == "rank_edges" &&
      (!is.null(x$boundary) || !is.null(x$delta) ||
       !is.character(x$ties) || length(x$ties) != 1L || is.na(x$ties) ||
       !x$ties %in% allowed_ties)) {
    .input_error("Invalid edge-ranking tie policy.")
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

# The reducer's validator, beside the reducer. It was written in
# `R/pairing.R`, one file above, so the record and the only statement of what
# makes it canonical sat on opposite sides of a cycle: operations called up to
# pairing to check a value pairing had called down to operations to build.
.validate_partition_reducer <- function(reducer) {
  expected <- c("kind", "weight_convention", "order")
  if (!.sealed_fields(reducer, "effect_partition_reducer", expected) ||
      !identical(reducer$kind, "weighted_sum") ||
      !identical(reducer$weight_convention, "normalized_unit_mass") ||
      !reducer$order %in% c("edge_first", "aggregate_first")) {
    .input_error("Partition reducer fields are missing or noncanonical.")
  }
  invisible(reducer)
}

#' Reduce normalized and transformed partition edges
#'
#' This is the default estimand: normalize each declared edge, transform it,
#' and only then apply the partition weights. Choose it when the quantity you
#' mean is the average of per-partition normalized edges.
#'
#' @return An `effect_partition_reducer` recording `$kind`,
#'   `$weight_convention`, and the stage `$order` (`"edge_first"`).
#' @seealso [aggregate_first()], the other stage order, and [pairing()],
#'   which supplies the partition weights this reducer applies.
#' @family generalization pairings
#' @examples
#' # Edge-first: each partition pair is normalized and transformed before the
#' # pairing weights are applied.
#' reducer <- reduce_partitions()
#' reducer$order
#'
#' # The two reducers name different estimands, so they are not
#' # interchangeable defaults.
#' c(edge_first = reduce_partitions()$order,
#'   aggregate_first = aggregate_first()$order)
#' @export
reduce_partitions <- function() .new_partition_reducer("edge_first")

#' Aggregate edge sufficient statistics before normalization
#'
#' This names a distinct estimand. Raw sufficient statistics are combined by
#' partition weights before normalization and transformation. It is the
#' default for [measurement_form()] and [coupling()], where averaging raw
#' partition products and normalizing once is the intended estimator.
#'
#' @return An `effect_partition_reducer` recording `$kind`,
#'   `$weight_convention`, and the stage `$order` (`"aggregate_first"`).
#' @seealso [reduce_partitions()] for the edge-first order, and
#'   [measurement_form()], which records this choice in its plan identity.
#' @family geometry plans and views
#' @examples
#' # Aggregate-first: raw partition products are pooled, then normalized once.
#' reducer <- aggregate_first()
#' reducer$order
#' reducer$weight_convention
#'
#' # This is the default measurement-pipeline order, so it is what
#' # `measurement_form()` and `coupling()` record when `reducer` is omitted.
#' identical(reducer, eval(formals(measurement_form)$reducer))
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
    .input_error(
      "Fisher transformation requires correlation-valued edge input."
    )
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
    signature = .sha256_signature(semantic)
  )), class = "effect_edge_operation_plan")
}

.apply_edge_transform <- function(value, transform) {
  if (transform$kind == "identity") return(value)
  if (transform$kind == "rank_edges") {
    return(matrix(rank(as.vector(value), ties.method = transform$ties),
      nrow(value), ncol(value)))
  }
  if (transform$boundary == "error" && any(abs(value) >= 1)) {
    .input_error(
      "Fisher transformation encountered an absolute-correlation boundary."
    )
  }
  if (transform$boundary == "clip") {
    value <- pmax(-1 + transform$delta, pmin(1 - transform$delta, value))
  }
  atanh(value)
}
