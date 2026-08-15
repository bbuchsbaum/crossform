# Additive spatial frame specifications -----------------------------------

.frame_spec <- function(kind, normalization, ...) {
  normalization <- match.arg(normalization, c("none", "local", "conservative"))
  structure(
    c(list(kind = kind, normalization = normalization), list(...)),
    class = "effect_frame_spec"
  )
}

#' Specify point measurements
#'
#' @param normalization Explicit frame normalization.
#' @return An additive frame specification.
#' @export
voxels <- function(normalization = "conservative") {
  .frame_spec("voxels", normalization)
}

#' Specify Euclidean searchlights
#'
#' The default `"local"` normalization matches the conventional
#' center-assigned searchlight map, in which overlapping neighborhoods
#' repeatedly count shared features. Local values under this default cannot
#' be summed into a whole-brain quantity. Request
#' `normalization = "conservative"` for exact local-to-global accounting of
#' the `total` component, and check any frame with [frame_conservation()].
#'
#' @param radius Positive radius in domain coordinate units.
#' @param normalization Explicit frame normalization.
#' @return An additive frame specification.
#' @export
searchlights <- function(radius, normalization = "local") {
  if (!is.numeric(radius) || length(radius) != 1L || is.na(radius) ||
      !is.finite(radius) || radius <= 0) {
    stop("`radius` must be one positive finite number.", call. = FALSE)
  }
  .frame_spec("searchlights", normalization, radius = radius)
}

#' Specify region measurements
#'
#' @param labels One region label per neural feature. Missing labels are
#'   excluded unless conservative normalization is requested.
#' @param normalization Explicit frame normalization.
#' @return An additive frame specification.
#' @export
regions <- function(labels, normalization = "local") {
  if (length(labels) < 1L || !(is.atomic(labels) || is.factor(labels))) {
    stop("`labels` must provide one atomic region label per feature.",
      call. = FALSE)
  }
  .frame_spec("regions", normalization, labels = labels)
}

#' Specify a whole-brain additive measurement
#'
#' @param normalization Explicit frame normalization.
#' @return An additive frame specification.
#' @export
whole_brain <- function(normalization = "local") {
  .frame_spec("whole_brain", normalization)
}

#' Compile a spatial specification against a neural domain
#'
#' Voxel, searchlight, region, and whole-brain scopes all compile to the same
#' sparse additive-frame representation.
#'
#' @param specification An additive frame specification.
#' @param domain An `effect_domain`.
#' @return An `effect_frame` whose weights are a sparse measurement-by-feature
#'   matrix.
#' @export
compile_frame <- function(specification, domain) {
  .validate_domain(domain)
  specification <- .validate_frame_specification(specification)
  n <- domain$n_features
  support_index <- NULL
  kind <- specification$kind
  if (kind == "voxels") {
    weights <- Matrix::Diagonal(n, x = 1)
    index <- data.frame(measurement = domain$feature_ids,
      stringsAsFactors = FALSE)
  } else if (kind == "whole_brain") {
    weights <- Matrix::sparseMatrix(i = rep(1L, n), j = seq_len(n), x = 1,
      dims = c(1L, n))
    index <- data.frame(measurement = "whole_brain", stringsAsFactors = FALSE)
  } else if (kind == "regions") {
    labels <- specification$labels
    if (length(labels) != n) {
      stop("Region labels must have one entry per domain feature.", call. = FALSE)
    }
    present <- !is.na(labels) & nzchar(as.character(labels))
    if (!any(present)) stop("At least one region label is required.", call. = FALSE)
    region_ids <- unique(as.character(labels[present]))
    region_index <- match(as.character(labels[present]), region_ids)
    weights <- Matrix::sparseMatrix(
      i = region_index, j = which(present), x = 1,
      dims = c(length(region_ids), n)
    )
    index <- data.frame(measurement = region_ids, stringsAsFactors = FALSE)
  } else if (kind == "searchlights") {
    radius <- specification$radius
    if (!is.numeric(radius) || length(radius) != 1L || is.na(radius) ||
        !is.finite(radius) || radius <= 0) {
      stop("Searchlight radius is invalid.", call. = FALSE)
    }
    support_index <- .euclidean_support_index(domain, radius)
    weights <- .support_index_membership(support_index)
    index <- data.frame(measurement = domain$feature_ids,
      stringsAsFactors = FALSE)
  } else {
    stop("Unknown additive frame specification.", call. = FALSE)
  }

  weights <- .normalize_frame(weights, specification$normalization)
  result <- additive_frame(weights, normalization = specification$normalization,
    domain = domain)
  result$index <- index
  result$domain_kind <- domain$kind
  result$specification <- specification
  if (!is.null(support_index)) {
    result$support_index <- support_index
  }
  result
}

.validate_frame_specification <- function(specification) {
  if (!inherits(specification, "effect_frame_spec") ||
      !is.list(specification) || !is.character(specification$kind) ||
      length(specification$kind) != 1L || is.na(specification$kind) ||
      !is.character(specification$normalization) ||
      length(specification$normalization) != 1L ||
      is.na(specification$normalization)) {
    stop("`specification` must be a valid additive frame specification.",
      call. = FALSE)
  }
  rebuilt <- switch(specification$kind,
    voxels = voxels(specification$normalization),
    searchlights = {
      if (!identical(names(specification), c("kind", "normalization", "radius"))) {
        stop("Frame specification fields are missing or noncanonical.", call. = FALSE)
      }
      searchlights(specification$radius, specification$normalization)
    },
    regions = {
      if (!identical(names(specification), c("kind", "normalization", "labels"))) {
        stop("Frame specification fields are missing or noncanonical.", call. = FALSE)
      }
      regions(specification$labels, specification$normalization)
    },
    whole_brain = whole_brain(specification$normalization),
    stop("Unknown additive frame specification.", call. = FALSE)
  )
  if (!identical(specification, rebuilt)) {
    stop("Frame specification fields are missing or noncanonical.", call. = FALSE)
  }
  rebuilt
}

#' Diagnose local-to-global conservation of a compiled frame
#'
#' Under a conservative frame every domain feature carries total weight mass
#' one, so local `total` geometries sum exactly to the global geometry:
#' \deqn{\sum_x G_x^{\mathrm{total}} = G_\Omega^{\mathrm{total}}.}
#' Overlapping neighborhoods under `normalization = "local"` double-count
#' shared features, so their local values are not contributions to a
#' whole-brain quantity. Two preconditions matter when checking the law:
#' the global comparator must be the unnormalized operator
#' (`whole_brain("none")`), because default local normalization divides by
#' the feature count; and the law covers the `total` component only —
#' `coherent` is defined by each measurement's own weighted common mode, and
#' local coherent values do not sum to the global coherent component.
#'
#' @param x A compiled `effect_frame`.
#' @param tolerance Nonnegative absolute per-feature mass tolerance.
#' @return An `effect_frame_conservation` list reporting `conserved`, the
#'   covered `component` (`"total"`), the frame `normalization`, the maximum
#'   per-feature mass deviation from one, and the per-feature mass vector.
#' @examples
#' domain <- abstract_domain(4, id = "conservation-example")
#' conservative <- compile_frame(voxels(), domain)
#' frame_conservation(conservative)$conserved
#' @export
frame_conservation <- function(x, tolerance = 1e-10) {
  .validate_frame_for_compile(x)
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || !is.finite(tolerance) || tolerance < 0) {
    stop("`tolerance` must be one nonnegative finite number.", call. = FALSE)
  }
  mass <- as.numeric(Matrix::colSums(x$weights))
  max_deviation <- max(abs(mass - 1))
  structure(list(
    conserved = max_deviation <= tolerance,
    component = "total",
    normalization = x$normalization,
    max_deviation = max_deviation,
    feature_mass = mass,
    tolerance = tolerance
  ), class = "effect_frame_conservation")
}

.normalize_frame <- function(weights, normalization) {
  row_mass <- Matrix::rowSums(weights)
  if (any(!is.finite(row_mass)) || any(row_mass <= 0)) {
    stop("Every frame measurement must contain at least one feature.",
      call. = FALSE)
  }
  if (normalization == "local") {
    weights <- Matrix::Diagonal(x = 1 / row_mass) %*% weights
  } else if (normalization == "conservative") {
    coverage <- Matrix::colSums(weights)
    if (any(!is.finite(coverage)) || any(coverage <= 0)) {
      stop("Conservative frames must cover every domain feature.", call. = FALSE)
    }
    weights <- weights %*% Matrix::Diagonal(x = 1 / coverage)
  }
  methods::as(methods::as(weights, "generalMatrix"), "CsparseMatrix")
}
