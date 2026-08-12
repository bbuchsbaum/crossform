# Neural feature domains ---------------------------------------------------

#' Construct an abstract neural feature domain
#'
#' @param n_features Positive neural feature count.
#' @param coordinates Optional finite feature-by-coordinate matrix used by
#'   spatial frame builders.
#' @param feature_ids Optional unique feature identifiers.
#' @param id Stable nonempty domain identity.
#' @return An `effect_domain`.
#' @export
abstract_domain <- function(n_features, coordinates = NULL,
                            feature_ids = NULL, id = "abstract") {
  n_features <- .domain_count(n_features)
  if (is.null(feature_ids)) feature_ids <- seq_len(n_features)
  if (length(feature_ids) != n_features || anyNA(feature_ids) ||
      anyDuplicated(feature_ids)) {
    stop("`feature_ids` must uniquely identify every neural feature.",
      call. = FALSE)
  }
  if (!is.null(coordinates) &&
      (!is.matrix(coordinates) || !is.numeric(coordinates) ||
       nrow(coordinates) != n_features || ncol(coordinates) < 1L ||
       any(!is.finite(coordinates)))) {
    stop("`coordinates` must be a finite feature-by-coordinate matrix.",
      call. = FALSE)
  }
  .domain_id(id)
  structure(
    list(
      id = id,
      kind = "abstract",
      n_features = n_features,
      feature_ids = feature_ids,
      coordinates = coordinates,
      metadata = list()
    ),
    class = "effect_domain"
  )
}

#' Construct a native volumetric neural feature domain
#'
#' @param mask A three-dimensional logical or numeric mask. Finite nonzero
#'   entries are included.
#' @param spacing Three positive finite voxel spacings.
#' @param id Stable domain identity.
#' @return An `effect_domain` with native voxel and physical coordinates.
#' @export
volume_domain <- function(mask, spacing = c(1, 1, 1), id = "native-volume") {
  if (!is.array(mask) || length(dim(mask)) != 3L || any(dim(mask) < 1L) ||
      !(is.logical(mask) || is.numeric(mask)) ||
      (is.numeric(mask) && any(!is.finite(mask)))) {
    stop("`mask` must be a finite three-dimensional logical or numeric array.",
      call. = FALSE)
  }
  included <- if (is.logical(mask)) !is.na(mask) & mask else mask != 0
  if (!any(included)) stop("`mask` must include at least one feature.", call. = FALSE)
  if (!is.numeric(spacing) || length(spacing) != 3L ||
      any(!is.finite(spacing)) || any(spacing <= 0)) {
    stop("`spacing` must contain three positive finite values.", call. = FALSE)
  }
  .domain_id(id)
  voxel <- arrayInd(which(included), dim(mask), useNames = FALSE)
  physical <- sweep(voxel - 1, 2L, spacing, `*`)
  structure(
    list(
      id = id,
      kind = "volume",
      n_features = as.integer(nrow(voxel)),
      feature_ids = which(included),
      coordinates = physical,
      metadata = list(dim = as.integer(dim(mask)), spacing = spacing,
        voxel = voxel, mask = included)
    ),
    class = "effect_domain"
  )
}

.domain_count <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1L || x %% 1 != 0 || x > .Machine$integer.max) {
    stop("`n_features` must be one positive integer.", call. = FALSE)
  }
  as.integer(x)
}

.domain_id <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("Domain `id` must be one nonempty identifier.", call. = FALSE)
  }
  invisible(x)
}

.validate_domain <- function(x) {
  expected <- c("id", "kind", "n_features", "feature_ids", "coordinates",
    "metadata")
  if (!inherits(x, "effect_domain") || !is.list(x) ||
      !identical(names(x), expected)) {
    stop("Domain fields are missing or noncanonical.", call. = FALSE)
  }
  .domain_id(x$id)
  n <- .domain_count(x$n_features)
  if (!is.character(x$kind) || length(x$kind) != 1L ||
      !x$kind %in% c("abstract", "volume")) {
    stop("Domain kind is invalid.", call. = FALSE)
  }
  if (length(x$feature_ids) != n || anyNA(x$feature_ids) ||
      anyDuplicated(x$feature_ids)) {
    stop("Domain feature identities are invalid.", call. = FALSE)
  }
  if (!is.null(x$coordinates) &&
      (!is.matrix(x$coordinates) || !is.numeric(x$coordinates) ||
       nrow(x$coordinates) != n || ncol(x$coordinates) < 1L ||
       any(!is.finite(x$coordinates)))) {
    stop("Domain coordinates are invalid.", call. = FALSE)
  }
  if (!is.list(x$metadata)) stop("Domain metadata must be a list.", call. = FALSE)
  invisible(x)
}
