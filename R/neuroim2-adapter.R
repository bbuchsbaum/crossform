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

#' Construct an effectagram domain from a neuroim2 volume mask
#'
#' This conditional adapter records the mask's stable full-volume indices,
#' physical coordinates, spacing, and full neuroim2 spatial metadata. It does
#' not extract response data or import neuroim2 result types.
#'
#' @param mask A three-dimensional neuroim2 `NeuroVol` mask.
#' @param id Stable domain identity.
#' @return An `effect_domain` whose `feature_ids` are full-volume indices.
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
    neuroim2_space_sha256 = paste0("sha256:", digest::digest(spatial_metadata,
      algo = "sha256", serialize = FALSE))
  )
  .new_domain(domain$id, domain$kind, domain$feature_ids, domain$coordinates,
    domain$coordinate_units, domain$metadata)
}

#' Compile neuroim2 searchlight indices into an effectagram frame
#'
#' The function calls only `neuroim2::searchlight_indices()` and maps its stable
#' full-volume indices to the ordered compact feature columns of `domain`.
#'
#' @param mask A three-dimensional neuroim2 `NeuroVol` mask.
#' @param radius Positive spherical radius in millimetres.
#' @param domain An exact domain from [neuroim2_volume_domain()]. When omitted,
#'   it is constructed from `mask`.
#' @param normalization One of `none`, `local`, or `conservative`.
#' @param nonzero Passed to `neuroim2::searchlight_indices()`; version 0.1
#'   requires `TRUE` so every member belongs to the compact domain.
#' @return A sparse additive `effect_frame`.
#' @export
neuroim2_searchlights <- function(mask, radius, domain = NULL,
                                  normalization = "local", nonzero = TRUE) {
  .require_neuroim2_searchlight_indices()
  if (!is.logical(nonzero) || length(nonzero) != 1L || is.na(nonzero) ||
      !nonzero) {
    stop("effectagram neuroim2 searchlights require `nonzero = TRUE`.",
      call. = FALSE)
  }
  if (is.null(domain)) domain <- neuroim2_volume_domain(mask)
  .validate_domain(domain)
  if (!identical(domain$kind, "volume")) {
    stop("`domain` must be a volume domain.", call. = FALSE)
  }
  mask_domain <- neuroim2_volume_domain(mask, id = domain$id)
  if (!.same_domain_reference(domain$reference, mask_domain$reference)) {
    stop("`mask` and `domain` have incompatible exact volume geometry.",
      call. = FALSE)
  }
  neighborhoods <- neuroim2::searchlight_indices(mask, radius,
    nonzero = TRUE)
  centers <- attr(neighborhoods, "center_indices", exact = TRUE)
  if (!identical(centers, domain$feature_ids)) {
    stop("neuroim2 searchlight centres do not match the ordered domain features.",
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
#' an effectagram volume domain. Features outside the domain receive `fill`.
#' This is an output adapter only; it performs no interpolation, smoothing, or
#' coordinate reinterpretation.
#'
#' @param values One finite numeric value per compact domain feature.
#' @param mask The three-dimensional neuroim2 `NeuroVol` whose geometry defined
#'   `domain`.
#' @param domain The exact domain from [neuroim2_volume_domain()].
#' @param fill Finite value written outside the compact domain.
#' @param label Optional result-volume label.
#' @return A neuroim2 `NeuroVol` with values at `domain$feature_ids`.
#' @export
as_neurovol <- function(values, mask, domain = NULL, fill = NA_real_,
                        label = "effectagram result") {
  .require_neuroim2_searchlight_indices()
  if (is.null(domain)) domain <- neuroim2_volume_domain(mask)
  .validate_domain(domain)
  if (!identical(domain$kind, "volume")) {
    stop("`domain` must be a volume domain.", call. = FALSE)
  }
  mask_domain <- neuroim2_volume_domain(mask, id = domain$id)
  if (!.same_domain_reference(domain$reference, mask_domain$reference)) {
    stop("`mask` and `domain` have incompatible exact volume geometry.",
      call. = FALSE)
  }
  if (!is.numeric(values) || is.matrix(values) ||
      length(values) != domain$n_features || any(!is.finite(values))) {
    stop("`values` must provide one finite number per compact domain feature.",
      call. = FALSE)
  }
  if (!is.numeric(fill) || length(fill) != 1L || is.nan(fill) ||
      is.infinite(fill)) {
    stop("`fill` must be one finite or missing numeric value.", call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label)) {
    stop("`label` must be one character string.", call. = FALSE)
  }
  payload <- array(as.double(fill), dim = dim(mask))
  payload[domain$feature_ids] <- as.double(values)
  neuroim2::NeuroVol(payload, neuroim2::space(mask), label = label)
}
