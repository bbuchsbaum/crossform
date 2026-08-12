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
  weights <- Matrix::sparseMatrix(
    i = rep(seq_along(members), counts),
    j = unlist(members, use.names = FALSE),
    x = 1,
    dims = c(length(members), domain$n_features)
  )
  weights <- .normalize_frame(weights, normalization)
  result <- additive_frame(weights, normalization = normalization,
    domain = domain)
  result$index <- data.frame(measurement = centers, stringsAsFactors = FALSE)
  result$domain_kind <- domain$kind
  result$specification <- list(kind = "neuroim2_searchlights",
    radius = as.numeric(radius), units = "mm", nonzero = TRUE,
    upstream_commit = "77b1ddb")
  result
}
