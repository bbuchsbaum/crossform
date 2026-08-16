# Conditional neuroim2 volume geometry adapter ------------------------------

.require_neuroim2_searchlight_indices <- function() {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("The neuroim2 adapter requires the suggested package `neuroim2`.",
      call. = FALSE)
  }
  if (!exists("searchlight_indices", envir = asNamespace("neuroim2"),
      inherits = FALSE)) {
    stop(paste0(
      "The installed neuroim2 does not provide `searchlight_indices()`; ",
      "install neuroim2 0.19.0 or later (upstream commit 77b1ddb)."
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Construct a crossform domain from a neuroim2 volume mask
#'
#' This conditional adapter records the mask's stable full-volume indices,
#' physical coordinates, spacing, and full neuroim2 spatial metadata. It does
#' not extract response data or import neuroim2 result types.
#'
#' @param mask A three-dimensional neuroim2 `NeuroVol` mask.
#' @param id Stable domain identity.
#' @return An `effect_domain` with `$kind` `"volume"`, `$feature_ids` giving
#'   the mask's full-volume indices, `$coordinates` in millimeters, and
#'   `$metadata` carrying `dim`, `spacing`, `voxel` indices, the logical `mask`,
#'   and a `neuroim2_space_sha256` hash of the full neuroim2 space.
#' @family neural domains and frames
#' @seealso [volume_domain()] for the same object from a plain array,
#'   [neuroim2_searchlights()] for matching neighborhoods, and [as_neurovol()]
#'   to write compact results back out.
#' @examples
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     utils::packageVersion("neuroim2") >= "0.19.0") {
#'   values <- array(FALSE, c(5L, 5L, 4L))
#'   values[2:4, 2:4, 2:3] <- TRUE
#'   mask <- neuroim2::LogicalNeuroVol(
#'     values, neuroim2::NeuroSpace(c(5L, 5L, 4L), spacing = c(3, 3, 3))
#'   )
#'   domain <- neuroim2_volume_domain(mask, id = "subject-mask")
#'
#'   # Compact features, addressed by their stable full-volume indices, so a
#'   # result vector can always be placed back in the original array.
#'   print(domain$n_features)
#'   print(identical(domain$feature_ids, which(values)))
#'
#'   # Physical spacing and a hash of the full neuroim2 space are recorded;
#'   # any frame built later must agree with that geometry.
#'   print(domain$metadata$spacing)
#'   print(substr(domain$metadata$neuroim2_space_sha256, 1, 24))
#' }
#' @export
neuroim2_volume_domain <- function(mask, id = "neuroim2-volume") {
  .require_neuroim2_searchlight_indices()
  if (!inherits(mask, "NeuroVol") || length(dim(mask)) != 3L) {
    stop("`mask` must be a three-dimensional neuroim2 NeuroVol.", call. = FALSE)
  }
  values <- as.array(mask)
  included <- is.finite(values) & values != 0
  if (!any(included)) stop("`mask` must include at least one feature.", call. = FALSE)
  feature_ids <- which(included)
  grid <- neuroim2::index_to_grid(mask, feature_ids)
  spacing <- as.numeric(neuroim2::spacing(mask))[1:3]
  physical <- sweep(grid - 1, 2L, spacing, `*`)
  spatial_metadata <- serialize(neuroim2::space(mask), NULL, version = 3)
  domain <- abstract_domain(
    length(feature_ids), coordinates = physical, feature_ids = feature_ids,
    id = id, coordinate_units = "mm"
  )
  domain$kind <- "volume"
  domain$metadata <- list(
    dim = as.integer(dim(mask)),
    spacing = spacing,
    voxel = unname(grid),
    mask = included,
    neuroim2_space_sha256 = .sha256_string(spatial_metadata)
  )
  .new_domain(domain$id, domain$kind, domain$feature_ids, domain$coordinates,
    domain$coordinate_units, domain$metadata)
}

#' Compile neuroim2 searchlight indices into a crossform frame
#'
#' The function calls only `neuroim2::searchlight_indices()` and maps its stable
#' full-volume indices to the ordered compact feature columns of `domain`.
#'
#' @param mask A three-dimensional neuroim2 `NeuroVol` mask.
#' @param radius Positive spherical radius in millimeters.
#' @param domain An exact domain from [neuroim2_volume_domain()]. When omitted,
#'   it is constructed from `mask`.
#' @param normalization One of `none`, `local`, or `conservative`.
#' @param nonzero Passed to `neuroim2::searchlight_indices()`; version 0.1
#'   requires `TRUE` so every member belongs to the compact domain.
#' @return An `effect_frame` whose `$weights` are the sparse
#'   measurement-by-feature operator, with `$index$measurement` holding the
#'   full-volume center indices, `$normalization`, a `$specification`
#'   recording the radius and the pinned `upstream_commit`, and a
#'   `$support_index`.
#' @family neural domains and frames
#' @seealso [searchlights()] plus [compile_frame()] for the built-in
#'   neighborhood builder, [neuroim2_volume_domain()] for the domain, and
#'   [frame_conservation()] to check the normalization you chose.
#' @examples
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     utils::packageVersion("neuroim2") >= "0.19.0") {
#'   values <- array(FALSE, c(5L, 5L, 4L))
#'   values[2:4, 2:4, 2:3] <- TRUE
#'   mask <- neuroim2::LogicalNeuroVol(
#'     values, neuroim2::NeuroSpace(c(5L, 5L, 4L), spacing = c(3, 3, 3))
#'   )
#'   domain <- neuroim2_volume_domain(mask)
#'
#'   # Neighborhoods come from neuroim2; crossform only maps their full-volume
#'   # member indices onto the ordered compact feature columns.
#'   frame <- neuroim2_searchlights(mask, radius = 4, domain = domain)
#'   print(dim(frame$weights))
#'   print(identical(frame$index$measurement, domain$feature_ids))
#'
#'   # The pinned upstream geometry is part of the frame's specification.
#'   print(frame$specification$upstream_commit)
#'
#'   # A mask whose geometry differs from the declared domain is refused.
#'   moved <- neuroim2::LogicalNeuroVol(
#'     values, neuroim2::NeuroSpace(c(5L, 5L, 4L), spacing = c(2, 2, 4))
#'   )
#'   print(try(neuroim2_searchlights(moved, radius = 4, domain = domain)))
#' }
#' @export
neuroim2_searchlights <- function(mask, radius, domain = NULL,
                                  normalization = "local", nonzero = TRUE) {
  if (missing(mask) || missing(radius)) {
    stop(paste0(
      "`mask` and `radius` are both required: pass the `NeuroVol` mask and a ",
      "positive spherical radius in millimetres, for example ",
      "`neuroim2_searchlights(mask, radius = 8)`."
    ), call. = FALSE)
  }
  .require_neuroim2_searchlight_indices()
  if (!is.logical(nonzero) || length(nonzero) != 1L || is.na(nonzero) ||
      !nonzero) {
    stop("crossform neuroim2 searchlights require `nonzero = TRUE`.",
      call. = FALSE)
  }
  if (is.null(domain)) domain <- neuroim2_volume_domain(mask)
  .validate_domain(domain)
  if (!identical(domain$kind, "volume")) {
    stop(sprintf(paste0(
      "`domain` must be a volume domain from `neuroim2_volume_domain()` or ",
      "`volume_domain()`; received a `%s` domain."
    ), domain$kind), call. = FALSE)
  }
  mask_domain <- neuroim2_volume_domain(mask, id = domain$id)
  if (!.same_domain_reference(domain$reference, mask_domain$reference)) {
    stop(sprintf(paste0(
      "`mask` and `domain` have different volume geometry, so a compact ",
      "index in one does not name the same voxel in the other. The mask is ",
      "%s with %s and %s; domain `%s` is %s with %s and %s. Pass the mask ",
      "the domain was built from."
    ),
      paste(mask_domain$metadata$dim, collapse = " x "),
      paste0("spacing ", paste(format(mask_domain$metadata$spacing),
        collapse = " x ")),
      .msg_count(mask_domain$n_features, "in-mask voxel"),
      domain$id,
      paste(domain$metadata$dim, collapse = " x "),
      paste0("spacing ", paste(format(domain$metadata$spacing),
        collapse = " x ")),
      .msg_count(domain$n_features, "feature")), call. = FALSE)
  }
  neighborhoods <- neuroim2::searchlight_indices(mask, radius,
    nonzero = TRUE)
  centers <- attr(neighborhoods, "center_indices", exact = TRUE)
  if (!identical(centers, domain$feature_ids)) {
    stop("neuroim2 searchlight centers do not match the ordered domain features.",
      call. = FALSE)
  }
  members <- lapply(neighborhoods, function(indices) {
    mapped <- match(indices, domain$feature_ids)
    if (anyNA(mapped)) {
      stop("A neuroim2 searchlight member lies outside the declared domain.",
        call. = FALSE)
    }
    as.integer(mapped)
  })
  counts <- lengths(members)
  if (length(members) < 1L || any(counts < 1L)) {
    stop("Every neuroim2 searchlight must contain at least one domain feature.",
      call. = FALSE)
  }
  support_index <- .support_index_from_members(
    members, domain, centers,
    construction = list(
      kind = "euclidean_ball",
      provider = "neuroim2_searchlight_indices",
      radius = as.numeric(radius),
      coordinate_units = domain$coordinate_units,
      upstream_commit = "77b1ddb"
    )
  )
  weights <- .support_index_membership(support_index)
  weights <- .normalize_frame(weights, normalization)
  result <- additive_frame(weights, normalization = normalization,
    domain = domain)
  result$index <- data.frame(measurement = centers, stringsAsFactors = FALSE)
  result$domain_kind <- domain$kind
  result$specification <- list(kind = "neuroim2_searchlights",
    radius = as.numeric(radius), units = "mm", nonzero = TRUE,
    upstream_commit = "77b1ddb")
  result$support_index <- support_index
  result
}

#' Map a compact result vector back to a neuroim2 volume
#'
#' The compact values are inserted at the exact full-volume indices carried by
#' a crossform volume domain. Features outside the domain receive `fill`.
#' This is an output adapter only; it performs no interpolation, smoothing, or
#' coordinate reinterpretation.
#'
#' @param values One finite numeric value per compact domain *feature* (voxel),
#'   in `domain$feature_ids` order. crossform result views carry one value per
#'   *measurement* instead, which coincides with the features only for a
#'   voxelwise or searchlight frame. For a coarser frame, expand first with the
#'   frame's membership pattern —
#'   `as.numeric(Matrix::crossprod(frame$weights != 0, values))` — as shown in
#'   the "Measurements are not features" section of
#'   `vignette("neuroim2-data")`.
#' @param mask The three-dimensional neuroim2 `NeuroVol` whose geometry defined
#'   `domain`.
#' @param domain The exact domain from [neuroim2_volume_domain()].
#' @param fill Finite value written outside the compact domain.
#' @param label Optional result-volume label.
#' @return A neuroim2 `NeuroVol` on the same space as `mask`, carrying `values`
#'   at `domain$feature_ids` and `fill` everywhere else.
#' @family neural domains and frames
#' @seealso [neuroim2_volume_domain()] for the domain whose `feature_ids` fix
#'   the output positions, and [geometry_component()] for one source of the
#'   compact vector.
#' @examples
#' if (requireNamespace("neuroim2", quietly = TRUE) &&
#'     utils::packageVersion("neuroim2") >= "0.19.0") {
#'   values <- array(FALSE, c(5L, 5L, 4L))
#'   values[2:4, 2:4, 2:3] <- TRUE
#'   mask <- neuroim2::LogicalNeuroVol(
#'     values, neuroim2::NeuroSpace(c(5L, 5L, 4L), spacing = c(3, 3, 3))
#'   )
#'   domain <- neuroim2_volume_domain(mask)
#'
#'   # One number per compact feature, in domain feature order.
#'   statistic <- seq_len(domain$n_features) / domain$n_features
#'   volume <- as_neurovol(statistic, mask, domain, label = "example statistic")
#'
#'   # Values land at exactly the mask indices; everything else stays `fill`.
#'   print(dim(volume))
#'   print(identical(as.numeric(volume[domain$feature_ids]), statistic))
#'   print(all(is.na(as.array(volume)[!values])))
#' }
#' @export
as_neurovol <- function(values, mask, domain = NULL, fill = NA_real_,
                        label = "crossform result") {
  if (missing(values) || missing(mask)) {
    stop(paste0(
      "`values` and `mask` are both required: `as_neurovol()` writes one ",
      "number per compact domain feature onto the space of the `NeuroVol` ",
      "mask the domain was built from."
    ), call. = FALSE)
  }
  .require_neuroim2_searchlight_indices()
  if (is.null(domain)) domain <- neuroim2_volume_domain(mask)
  .validate_domain(domain)
  if (!identical(domain$kind, "volume")) {
    stop(sprintf(paste0(
      "`domain` must be a volume domain from `neuroim2_volume_domain()` or ",
      "`volume_domain()`; received a `%s` domain."
    ), domain$kind), call. = FALSE)
  }
  mask_domain <- neuroim2_volume_domain(mask, id = domain$id)
  if (!.same_domain_reference(domain$reference, mask_domain$reference)) {
    stop(sprintf(paste0(
      "`mask` and `domain` have different volume geometry, so a compact ",
      "index in one does not name the same voxel in the other. The mask is ",
      "%s with %s and %s; domain `%s` is %s with %s and %s. Pass the mask ",
      "the domain was built from."
    ),
      paste(mask_domain$metadata$dim, collapse = " x "),
      paste0("spacing ", paste(format(mask_domain$metadata$spacing),
        collapse = " x ")),
      .msg_count(mask_domain$n_features, "in-mask voxel"),
      domain$id,
      paste(domain$metadata$dim, collapse = " x "),
      paste0("spacing ", paste(format(domain$metadata$spacing),
        collapse = " x ")),
      .msg_count(domain$n_features, "feature")), call. = FALSE)
  }
  if (!is.numeric(values) || is.matrix(values)) {
    stop(sprintf(paste0(
      "`values` must be a numeric vector with one value per compact domain ",
      "feature; received %s. A view's `$values` matrix has one column per ",
      "query, so select the column you want to write out."
    ), .msg_value(values)), call. = FALSE)
  }
  if (length(values) != domain$n_features) {
    stop(sprintf(paste0(
      "`values` has %s but domain `%s` has %s. crossform result views are ",
      "one value per *measurement* (a searchlight, a region, the whole ",
      "brain), while `as_neurovol()` writes one value per *feature* (a ",
      "voxel). They coincide only for a searchlight or voxelwise frame. To ",
      "write a coarser frame out, expand its measurements to voxels first ",
      "with the frame's own membership pattern: ",
      "`as.numeric(Matrix::crossprod(frame$weights != 0, values))`. Use the ",
      "membership pattern rather than `frame$weights`, whose local ",
      "normalization would rescale the numbers. See the \"Measurements are ",
      "not features\" section of `vignette(\"neuroim2-data\")`."
    ), .msg_count(length(values), "value"), domain$id,
      .msg_count(domain$n_features, "feature")), call. = FALSE)
  }
  if (any(!is.finite(values))) {
    stop(sprintf(paste0(
      "`values` must be finite; %d of %d are NA, NaN, or Inf. `as_neurovol()` ",
      "will not guess what a missing measurement means -- replace them ",
      "deliberately, or build the domain from the coverage you actually have."
    ), sum(!is.finite(values)), length(values)), call. = FALSE)
  }
  if (!is.numeric(fill) || length(fill) != 1L || is.nan(fill) ||
      is.infinite(fill)) {
    stop(sprintf(paste0(
      "`fill` must be one numeric value (`NA_real_` is allowed) written ",
      "outside the compact domain; received %s."
    ), .msg_value(fill)), call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label)) {
    stop(sprintf("`label` must be one character string; received %s.",
      .msg_value(label)), call. = FALSE)
  }
  payload <- array(as.double(fill), dim = dim(mask))
  payload[domain$feature_ids] <- as.double(values)
  neuroim2::NeuroVol(payload, neuroim2::space(mask), label = label)
}
