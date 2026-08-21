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
    .input_error(paste0(
      "`n_features` is required: give the number of neural features, for ",
      "example `abstract_domain(ncol(betas))`."
    ),
      arg = "n_features", received = "no argument",
      expected = "the number of neural features")
  }
  n_features <- .domain_count(n_features)
  if (is.null(feature_ids)) feature_ids <- seq_len(n_features)
  if (length(feature_ids) != n_features || anyNA(feature_ids) ||
      anyDuplicated(feature_ids)) {
    .input_error(sprintf(paste0(
      "`feature_ids` must uniquely identify every neural feature: expected ",
      "%s, received %s%s."
    ), .msg_count(n_features, "identifier"),
      .msg_count(length(feature_ids), "value"),
      if (anyDuplicated(feature_ids)) " with repeats" else ""),
      arg = "feature_ids",
      received = .msg_count(length(feature_ids), "value"),
      expected = .msg_count(n_features, "unique identifier"))
  }
  if (!is.null(coordinates) &&
      (!is.matrix(coordinates) || !is.numeric(coordinates) ||
       ncol(coordinates) < 1L)) {
    .input_error(sprintf(paste0(
      "`coordinates` must be a numeric feature-by-axis matrix with one row ",
      "per feature and at least one coordinate axis; received %s."
    ), .msg_value(coordinates)),
      arg = "coordinates", received = .msg_value(coordinates),
      expected = "a numeric feature-by-axis matrix")
  }
  if (!is.null(coordinates) && nrow(coordinates) != n_features) {
    .input_error(sprintf(paste0(
      "`coordinates` has %s but the domain declares %s; supply one ",
      "coordinate row per feature, in feature order."
    ), .msg_count(nrow(coordinates), "row"),
      .msg_count(n_features, "feature")),
      arg = "coordinates",
      received = .msg_count(nrow(coordinates), "row"),
      expected = .msg_count(n_features, "row"))
  }
  if (!is.null(coordinates) && any(!is.finite(coordinates))) {
    .input_error(sprintf(paste0(
      "`coordinates` must be finite; %d of %d values are NA, NaN, or Inf."
    ), sum(!is.finite(coordinates)), length(coordinates)),
      arg = "coordinates",
      received = sprintf("%d non-finite of %d", sum(!is.finite(coordinates)),
        length(coordinates)),
      expected = "all coordinates finite")
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
#' @param metadata Optional uniquely named list of extra facts about where the
#'   geometry came from, recorded alongside the grid facts the constructor
#'   derives. See *Provider metadata*.
#' @return An `effect_domain` with `$kind` `"volume"`, `$feature_ids` giving
#'   the stable full-volume indices of the included voxels, `$coordinates`
#'   holding their physical positions, `$metadata` carrying `dim`, `spacing`,
#'   `voxel` indices and the logical `mask` followed by anything `metadata`
#'   added, plus the usual `$geometry_signature` and `$reference`.
#' @section Provider metadata:
#' A domain built by an adapter usually knows something about the geometry
#' that the array itself does not carry --- the native header it came from,
#' the file, the transform it was resampled under. `metadata` is where that
#' goes, and it is not decoration: the domain's `$geometry_signature` covers
#' it, so two domains that agree on every voxel but disagree about their
#' provenance are correctly *different* domains, and anything holding a
#' compact result vector against one of them refuses the other.
#'
#' That is the point. `neuroim2_volume_domain()` records the hash of the full
#' `neuroim2` space this way, which is what makes writing a result back to
#' voxels safe. The grid facts the constructor derives itself --- `dim`,
#' `spacing`, `voxel`, `mask` --- cannot be overridden.
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
#'
#' # A provider records where the geometry came from, and the record is part
#' # of what the domain is: the same voxels under a different provenance are a
#' # different domain, not the same one.
#' stamped <- volume_domain(mask, spacing = c(3, 3, 3), id = "example-volume",
#'   metadata = list(source_file = "sub-01_mask.nii.gz"))
#' stamped$metadata$source_file
#' identical(stamped$reference, domain$reference)
#' @export
volume_domain <- function(mask, spacing = c(1, 1, 1), id = "native-volume",
                          coordinate_units = "mm", metadata = list()) {
  if (missing(mask)) {
    .input_error(paste0(
      "`mask` is required: pass a three-dimensional logical or numeric array ",
      "whose nonzero entries are the in-mask voxels."
    ),
      arg = "mask", received = "no argument",
      expected = "a three-dimensional logical or numeric array")
  }
  if (!is.array(mask) || length(dim(mask)) != 3L || any(dim(mask) < 1L) ||
      !(is.logical(mask) || is.numeric(mask)) ||
      (is.numeric(mask) && any(!is.finite(mask)))) {
    .input_error(sprintf(paste0(
      "`mask` must be a finite three-dimensional logical or numeric array; ",
      "received %s."
    ), .msg_value(mask)),
      arg = "mask", received = .msg_value(mask),
      expected = "a finite three-dimensional logical or numeric array")
  }
  included <- if (is.logical(mask)) !is.na(mask) & mask else mask != 0
  if (!any(included)) {
    .input_error(paste0(
      "`mask` selects no voxels: every entry is FALSE, zero, or missing, so ",
      "the domain would have no features."
    ),
      arg = "mask", received = "no nonzero entries",
      expected = "at least one in-mask voxel")
  }
  if (!.is_finite_numeric(spacing) || length(spacing) != 3L ||
      any(spacing <= 0)) {
    .input_error(sprintf(paste0(
      "`spacing` must contain three positive finite voxel sizes, one per ",
      "array axis; received %s."
    ), .msg_value(spacing)),
      arg = "spacing", received = .msg_value(spacing),
      expected = "three positive finite voxel sizes")
  }
  grid_facts <- c("dim", "spacing", "voxel", "mask")
  if (!is.list(metadata) || (length(metadata) &&
      (is.null(names(metadata)) || anyNA(names(metadata)) ||
       any(!nzchar(names(metadata))) || anyDuplicated(names(metadata))))) {
    .input_error(sprintf(paste0(
      "`metadata` must be a uniquely named list of extra facts about where ",
      "this geometry came from; received %s."
    ), .msg_value(metadata)),
      arg = "metadata", received = .msg_value(metadata),
      expected = "a uniquely named list")
  }
  overridden <- intersect(names(metadata), grid_facts)
  if (length(overridden)) {
    .input_error(sprintf(paste0(
      "`metadata` may not restate the grid facts this constructor derives ",
      "from `mask` and `spacing` (%s); received %s."
    ), paste(grid_facts, collapse = ", "), .msg_names(overridden)),
      arg = "metadata", received = .msg_names(overridden),
      expected = "names outside the derived grid facts")
  }
  .domain_id(id)
  voxel <- arrayInd(which(included), dim(mask), useNames = FALSE)
  physical <- sweep(voxel - 1, 2L, spacing, `*`)
  coordinate_units <- .domain_coordinate_units(coordinate_units, physical)
  .new_domain(id, "volume", which(included), physical, coordinate_units,
    metadata = c(list(dim = as.integer(dim(mask)), spacing = spacing,
      voxel = voxel, mask = included), metadata))
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
  if (!.is_strings(x) || !length(x) %in% c(1L, axes)) {
    .input_error(
      "`coordinate_units` must provide one nonempty unit or one per axis."
    )
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
    .input_error("Domain feature identities are invalid.")
  }
  if (!.is_strings(coordinate_units) || length(coordinate_units) < 1L) {
    .input_error("Domain coordinate units are invalid.")
  }
  if (!.strong_sha256(geometry_signature)) {
    .input_error("Domain geometry signature is invalid.")
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
  if (!.sealed_fields(x, "effect_domain_reference", expected)) {
    .input_error("Domain-reference fields are missing or noncanonical.")
  }
  rebuilt <- .new_domain_reference(x$id, x$feature_ids, x$coordinate_units,
    x$geometry_signature)
  if (!identical(x, rebuilt)) {
    .contract_error("Domain-reference metadata or signature is inconsistent.")
  }
  rebuilt
}

.same_domain_reference <- function(x, y) {
  x <- .validate_domain_reference(x)
  y <- .validate_domain_reference(y)
  identical(x$signature, y$signature) && identical(x, y)
}

.domain_count <- function(x) {
  if (!.is_number(x) || x < 1L || x %% 1 != 0 || x > .Machine$integer.max) {
    .input_error(sprintf(paste0(
      "`n_features` must be one positive whole number giving the count of ",
      "neural features; received %s."
    ), if (is.numeric(x) && length(x) == 1L && !is.na(x)) {
      paste0("`", format(x), "`")
    } else {
      .msg_value(x)
    }),
      arg = "n_features", received = .msg_value(x),
      expected = "one positive whole number")
  }
  as.integer(x)
}

.domain_id <- function(x) {
  if (!.is_string(x)) {
    .input_error("Domain `id` must be one nonempty identifier.")
  }
  invisible(x)
}

.validate_domain <- function(x) {
  if (!inherits(x, "effect_domain")) {
    .input_error(sprintf(paste0(
      "Expected an `effect_domain` (see `abstract_domain()`, ",
      "`volume_domain()`, or `neuroim2_volume_domain()`); received %s."
    ), .msg_value(x)),
      arg = "x", received = .msg_value(x),
      expected = "an `effect_domain`")
  }
  expected <- c("id", "kind", "n_features", "feature_ids", "coordinates",
    "coordinate_units", "geometry_signature", "reference", "metadata")
  if (!.sealed_fields(x, "effect_domain", expected)) {
    .input_error("Domain fields are missing or noncanonical.")
  }
  .domain_id(x$id)
  n <- .domain_count(x$n_features)
  if (!is.character(x$kind) || length(x$kind) != 1L ||
      !x$kind %in% c("abstract", "volume")) {
    .input_error("Domain kind is invalid.")
  }
  if (length(x$feature_ids) != n || anyNA(x$feature_ids) ||
      anyDuplicated(x$feature_ids)) {
    .input_error("Domain feature identities are invalid.")
  }
  if (!is.null(x$coordinates) &&
      (!is.matrix(x$coordinates) || !is.numeric(x$coordinates) ||
       nrow(x$coordinates) != n || ncol(x$coordinates) < 1L ||
       any(!is.finite(x$coordinates)))) {
    .input_error("Domain coordinates are invalid.")
  }
  if (!is.list(x$metadata)) .input_error("Domain metadata must be a list.")
  units <- .domain_coordinate_units(x$coordinate_units, x$coordinates)
  geometry_signature <- .domain_geometry_signature(x$kind, x$coordinates,
    units, x$metadata)
  .check_signature(
    x$geometry_signature, geometry_signature,
    "Domain geometry signature is inconsistent."
  )
  reference <- .new_domain_reference(x$id, x$feature_ids, units,
    geometry_signature)
  if (!identical(x$reference, reference)) {
    .contract_error("Domain reference is inconsistent with the full domain.")
  }
  invisible(x)
}
