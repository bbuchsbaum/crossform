# Neural feature domains ---------------------------------------------------

#' Construct an abstract neural feature domain
#'
#' @param n_features Positive neural feature count.
#' @param coordinates Optional finite feature-by-coordinate matrix used by
#'   spatial frame builders.
#' @param feature_ids Optional unique feature identifiers.
#' @param id Stable nonempty domain identity.
#' @param coordinate_units One coordinate unit or one per coordinate axis.
#' @return An `effect_domain`.
#' @export
abstract_domain <- function(n_features, coordinates = NULL,
                            feature_ids = NULL, id = "abstract",
                            coordinate_units = "arbitrary") {
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
  coordinate_units <- .domain_coordinate_units(coordinate_units, coordinates)
  .new_domain(id, "abstract", feature_ids, coordinates, coordinate_units,
    metadata = list())
}

#' Construct a native volumetric neural feature domain
#'
#' @param mask A three-dimensional logical or numeric mask. Finite nonzero
#'   entries are included.
#' @param spacing Three positive finite voxel spacings.
#' @param id Stable domain identity.
#' @param coordinate_units One physical coordinate unit or one per axis.
#' @return An `effect_domain` with native voxel and physical coordinates.
#' @export
volume_domain <- function(mask, spacing = c(1, 1, 1), id = "native-volume",
                          coordinate_units = "mm") {
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
  coordinate_units <- .domain_coordinate_units(coordinate_units, physical)
  .new_domain(id, "volume", which(included), physical, coordinate_units,
    metadata = list(dim = as.integer(dim(mask)), spacing = spacing,
      voxel = voxel, mask = included))
}

.new_domain <- function(id, kind, feature_ids, coordinates, coordinate_units,
                        metadata) {
  geometry_signature <- .domain_geometry_signature(kind, coordinates,
    coordinate_units, metadata)
  reference <- .new_domain_reference(id, feature_ids, coordinate_units,
    geometry_signature)
  structure(
    list(
      id = id,
      kind = kind,
      n_features = as.integer(length(feature_ids)),
      feature_ids = feature_ids,
      coordinates = coordinates,
      coordinate_units = coordinate_units,
      geometry_signature = geometry_signature,
      reference = reference,
      metadata = metadata
    ),
    class = "effect_domain"
  )
}

.domain_coordinate_units <- function(x, coordinates) {
  axes <- if (is.null(coordinates)) 1L else ncol(coordinates)
  if (!is.character(x) || !length(x) %in% c(1L, axes) || anyNA(x) ||
      any(!nzchar(x))) {
    stop("`coordinate_units` must provide one nonempty unit or one per axis.",
      call. = FALSE)
  }
  if (length(x) == 1L) x <- rep(x, axes)
  unname(x)
}

.domain_geometry_signature <- function(kind, coordinates, coordinate_units,
                                       metadata) {
  semantic <- list(kind = kind, coordinates = coordinates,
    coordinate_units = coordinate_units, metadata = metadata)
  paste0("sha256:", digest::digest(semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L))
}

.new_domain_reference <- function(id, feature_ids, coordinate_units,
                                  geometry_signature) {
  .domain_id(id)
  if (length(feature_ids) < 1L || anyNA(feature_ids) || anyDuplicated(feature_ids)) {
    stop("Domain feature identities are invalid.", call. = FALSE)
  }
  if (!is.character(coordinate_units) || length(coordinate_units) < 1L ||
      anyNA(coordinate_units) || any(!nzchar(coordinate_units))) {
    stop("Domain coordinate units are invalid.", call. = FALSE)
  }
  if (!is.character(geometry_signature) || length(geometry_signature) != 1L ||
      is.na(geometry_signature) ||
      !grepl("^sha256:[[:xdigit:]]{64}$", geometry_signature)) {
    stop("Domain geometry signature is invalid.", call. = FALSE)
  }
  semantic <- list(
    id = id,
    n_features = as.integer(length(feature_ids)),
    feature_ids = feature_ids,
    coordinate_units = coordinate_units,
    geometry_signature = geometry_signature
  )
  signature <- paste0("sha256:", digest::digest(semantic, algo = "sha256",
    serialize = TRUE, serializeVersion = 2L))
  structure(c(semantic, list(signature = signature)),
    class = "effect_domain_reference")
}

.positional_domain_reference <- function(n_features, id = "abstract") {
  n_features <- .domain_count(n_features)
  coordinate_units <- "arbitrary"
  geometry_signature <- .domain_geometry_signature("abstract", NULL,
    coordinate_units, list())
  .new_domain_reference(id, seq_len(n_features), coordinate_units,
    geometry_signature)
}

.domain_reference <- function(x) {
  if (inherits(x, "effect_domain_reference")) {
    return(.validate_domain_reference(x))
  }
  .validate_domain(x)
  x$reference
}

.validate_domain_reference <- function(x) {
  expected <- c("id", "n_features", "feature_ids", "coordinate_units",
    "geometry_signature", "signature")
  if (!inherits(x, "effect_domain_reference") || !is.list(x) ||
      !identical(names(x), expected)) {
    stop("Domain-reference fields are missing or noncanonical.", call. = FALSE)
  }
  rebuilt <- .new_domain_reference(x$id, x$feature_ids, x$coordinate_units,
    x$geometry_signature)
  if (!identical(x, rebuilt)) {
    stop("Domain-reference metadata or signature is inconsistent.",
      call. = FALSE)
  }
  rebuilt
}

.same_domain_reference <- function(x, y) {
  x <- .validate_domain_reference(x)
  y <- .validate_domain_reference(y)
  identical(x$signature, y$signature) && identical(x, y)
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
    "coordinate_units", "geometry_signature", "reference", "metadata")
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
  units <- .domain_coordinate_units(x$coordinate_units, x$coordinates)
  geometry_signature <- .domain_geometry_signature(x$kind, x$coordinates,
    units, x$metadata)
  if (!identical(x$geometry_signature, geometry_signature)) {
    stop("Domain geometry signature is inconsistent.", call. = FALSE)
  }
  reference <- .new_domain_reference(x$id, x$feature_ids, units,
    geometry_signature)
  if (!identical(x$reference, reference)) {
    stop("Domain reference is inconsistent with the full domain.", call. = FALSE)
  }
  invisible(x)
}
