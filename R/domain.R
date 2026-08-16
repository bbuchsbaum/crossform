# Neural feature domains ---------------------------------------------------

#' Construct an abstract neural feature domain
#'
#' Use this when the neural features are a plain ordered set — surface nodes,
#' parcels, electrodes, or columns of a beta matrix — and any geometry you have
#' is supplied directly as a coordinate matrix rather than read from a volume.
#'
#' @param n_features Positive neural feature count.
#' @param coordinates Optional finite feature-by-coordinate matrix used by
#'   spatial frame builders.
#' @param feature_ids Optional unique feature identifiers.
#' @param id Stable nonempty domain identity.
#' @param coordinate_units One coordinate unit or one per coordinate axis.
#' @return An `effect_domain`: a list with `$id`, `$kind` (`"abstract"`),
#'   `$n_features`, `$feature_ids`, the optional `$coordinates` matrix,
#'   `$coordinate_units`, a `$geometry_signature`, a compact `$reference` that
#'   other objects store instead of the full domain, and `$metadata`.
#' @family neural domains and frames
#' @seealso [volume_domain()] for a mask-derived domain and [compile_frame()]
#'   to turn a domain into measurements.
#' @examples
#' # Twelve features with no geometry: enough for region or whole-brain scopes.
#' domain <- abstract_domain(12L, id = "from-observations:v1")
#' domain$n_features
#' domain$kind
#'
#' # Supply coordinates when a spatial scope such as searchlights() will be
#' # compiled; the radius is then read in these coordinate units.
#' grid <- as.matrix(expand.grid(x = 1:3, y = 1:4))
#' placed <- abstract_domain(nrow(grid), coordinates = grid, id = "planar-grid")
#' dim(compile_frame(searchlights(1.5), placed)$weights)
#' @export
abstract_domain <- function(n_features, coordinates = NULL,
                            feature_ids = NULL, id = "abstract",
                            coordinate_units = "arbitrary") {
  if (missing(n_features)) {
    stop(paste0(
      "`n_features` is required: give the number of neural features, for ",
      "example `abstract_domain(ncol(betas))`."
    ), call. = FALSE)
  }
  n_features <- .domain_count(n_features)
  if (is.null(feature_ids)) feature_ids <- seq_len(n_features)
  if (length(feature_ids) != n_features || anyNA(feature_ids) ||
      anyDuplicated(feature_ids)) {
    stop(sprintf(paste0(
      "`feature_ids` must uniquely identify every neural feature: expected ",
      "%s, received %s%s."
    ), .msg_count(n_features, "identifier"),
      .msg_count(length(feature_ids), "value"),
      if (anyDuplicated(feature_ids)) " with repeats" else ""),
      call. = FALSE)
  }
  if (!is.null(coordinates) &&
      (!is.matrix(coordinates) || !is.numeric(coordinates) ||
       ncol(coordinates) < 1L)) {
    stop(sprintf(paste0(
      "`coordinates` must be a numeric feature-by-axis matrix with one row ",
      "per feature and at least one coordinate axis; received %s."
    ), .msg_value(coordinates)), call. = FALSE)
  }
  if (!is.null(coordinates) && nrow(coordinates) != n_features) {
    stop(sprintf(paste0(
      "`coordinates` has %s but the domain declares %s; supply one ",
      "coordinate row per feature, in feature order."
    ), .msg_count(nrow(coordinates), "row"),
      .msg_count(n_features, "feature")), call. = FALSE)
  }
  if (!is.null(coordinates) && any(!is.finite(coordinates))) {
    stop(sprintf(paste0(
      "`coordinates` must be finite; %d of %d values are NA, NaN, or Inf."
    ), sum(!is.finite(coordinates)), length(coordinates)), call. = FALSE)
  }
  .domain_id(id)
  coordinate_units <- .domain_coordinate_units(coordinate_units, coordinates)
  .new_domain(id, "abstract", feature_ids, coordinates, coordinate_units,
    metadata = list())
}

#' Construct a native volumetric neural feature domain
#'
#' Use this when features are the in-mask voxels of a 3D volume. The compact
#' feature order follows the full-volume index order of the mask, so results
#' can always be written back to the original array positions.
#'
#' @param mask A three-dimensional logical or numeric mask. Finite nonzero
#'   entries are included.
#' @param spacing Three positive finite voxel spacings.
#' @param id Stable domain identity.
#' @param coordinate_units One physical coordinate unit or one per axis.
#' @return An `effect_domain` with `$kind` `"volume"`, `$feature_ids` giving
#'   the stable full-volume indices of the included voxels, `$coordinates`
#'   holding their physical positions, `$metadata` carrying `dim`, `spacing`,
#'   `voxel` indices and the logical `mask`, plus the usual
#'   `$geometry_signature` and `$reference`.
#' @family neural domains and frames
#' @seealso [abstract_domain()] for non-volumetric features,
#'   [neuroim2_volume_domain()] to build the same object from a `NeuroVol`, and
#'   [compile_frame()] to place searchlights on it.
#' @examples
#' # A 5 x 5 x 3 volume with one voxel excluded from the mask.
#' mask <- array(TRUE, c(5L, 5L, 3L))
#' mask[1, 1, 1] <- FALSE
#' domain <- volume_domain(mask, spacing = c(3, 3, 3), id = "example-volume")
#' domain$n_features
#'
#' # Coordinates are millimeters, so a searchlight radius is physical too.
#' utils::head(domain$coordinates, 3)
#' nrow(compile_frame(searchlights(4), domain)$weights)
#'
#' # feature_ids are full-volume indices, which is how compact results are
#' # written back into the original array.
#' utils::head(domain$feature_ids, 3)
#' @export
volume_domain <- function(mask, spacing = c(1, 1, 1), id = "native-volume",
                          coordinate_units = "mm") {
  if (missing(mask)) {
    stop(paste0(
      "`mask` is required: pass a three-dimensional logical or numeric array ",
      "whose nonzero entries are the in-mask voxels."
    ), call. = FALSE)
  }
  if (!is.array(mask) || length(dim(mask)) != 3L || any(dim(mask) < 1L) ||
      !(is.logical(mask) || is.numeric(mask)) ||
      (is.numeric(mask) && any(!is.finite(mask)))) {
    stop(sprintf(paste0(
      "`mask` must be a finite three-dimensional logical or numeric array; ",
      "received %s."
    ), .msg_value(mask)), call. = FALSE)
  }
  included <- if (is.logical(mask)) !is.na(mask) & mask else mask != 0
  if (!any(included)) {
    stop(paste0(
      "`mask` selects no voxels: every entry is FALSE, zero, or missing, so ",
      "the domain would have no features."
    ), call. = FALSE)
  }
  if (!is.numeric(spacing) || length(spacing) != 3L ||
      any(!is.finite(spacing)) || any(spacing <= 0)) {
    stop(sprintf(paste0(
      "`spacing` must contain three positive finite voxel sizes, one per ",
      "array axis; received %s."
    ), .msg_value(spacing)), call. = FALSE)
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
  .sha256_signature(semantic)
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
  if (!.strong_sha256(geometry_signature)) {
    stop("Domain geometry signature is invalid.", call. = FALSE)
  }
  semantic <- list(
    id = id,
    n_features = as.integer(length(feature_ids)),
    feature_ids = feature_ids,
    coordinate_units = coordinate_units,
    geometry_signature = geometry_signature
  )
  signature <- .sha256_signature(semantic)
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
    stop(sprintf(paste0(
      "`n_features` must be one positive whole number giving the count of ",
      "neural features; received %s."
    ), if (is.numeric(x) && length(x) == 1L && !is.na(x)) {
      paste0("`", format(x), "`")
    } else {
      .msg_value(x)
    }), call. = FALSE)
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
  if (!inherits(x, "effect_domain")) {
    stop(sprintf(paste0(
      "Expected an `effect_domain` (see `abstract_domain()`, ",
      "`volume_domain()`, or `neuroim2_volume_domain()`); received %s."
    ), .msg_value(x)), call. = FALSE)
  }
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
