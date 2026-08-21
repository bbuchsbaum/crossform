# Conditional neuroim2 volume geometry adapter ------------------------------

.require_neuroim2_searchlight_indices <- function() {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    .input_error(
      "The neuroim2 adapter requires the suggested package `neuroim2`."
    )
  }
  if (!exists("searchlight_indices", envir = asNamespace("neuroim2"),
      inherits = FALSE)) {
    .input_error(paste0(
      "The installed neuroim2 does not provide `searchlight_indices()`; ",
      "install neuroim2 0.19.0 or later (upstream commit 77b1ddb)."
    ))
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
    .input_error("`mask` must be a three-dimensional neuroim2 NeuroVol.")
  }
  # Below the two questions only neuroim2 can answer -- what is in the mask,
  # and what is its physical spacing -- this is the ordinary volume domain any
  # provider builds, so the ordinary public constructor builds it. What the
  # adapter adds is the one fact `volume_domain()` cannot derive from an
  # array: the identity of the full neuroim2 space these voxels are addressed
  # in, which is what makes writing a result back to them safe later.
  volume_domain(
    array(as.vector(mask), dim = dim(mask)),
    spacing = as.numeric(neuroim2::spacing(mask))[1:3],
    id = id, coordinate_units = "mm",
    metadata = list(
      neuroim2_space_sha256 = .sha256_string(
        serialize(neuroim2::space(mask), NULL, version = 3)
      )
    )
  )
}

# Both entry points in this file take a `NeuroVol` mask and an optional
# domain, and both must establish the same thing before they can map a compact
# index onto a voxel: that the domain is a volume domain and that it names the
# same volume geometry as the mask. Getting that wrong writes numbers into the
# wrong voxels silently, which is why the refusal reports both geometries
# side by side rather than saying they differ.
.neuroim2_domain_for_mask <- function(mask, domain) {
  if (is.null(domain)) domain <- neuroim2_volume_domain(mask)
  .check_class(domain, "effect_domain", "domain",
    from = "neuroim2_volume_domain()")
  if (!identical(domain$kind, "volume")) {
    .input_error(sprintf(paste0(
      "`domain` must be a volume domain from `neuroim2_volume_domain()`; ",
      "received a `%s` domain."
    ), domain$kind))
  }
  # The comparison is against a whole domain rebuilt from the mask, not
  # against its reference, and that is deliberate: a domain identical to one
  # this constructor would have built is a *valid* domain, so the agreement
  # test and the validity test are the same test, and both are reachable from
  # outside the package. A domain whose recorded identity disagreed with its
  # own fields fails here for the same reason a domain from another mask does.
  mask_domain <- neuroim2_volume_domain(mask, id = domain$id)
  if (!identical(domain, mask_domain)) {
    .contract_error(sprintf(paste0(
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
      .msg_count(domain$n_features, "feature")))
  }
  domain
}

#' Compile neuroim2 searchlight indices into a crossform frame
#'
#' The function calls only `neuroim2::searchlight_indices()` and maps its stable
#' full-volume indices to the ordered compact feature columns of `domain`.
#'
#' Several radii request one conservative frame per radius, stacked into a
#' multiscale [frame_family()]; see *Multiscale families*.
#'
#' @param mask A three-dimensional neuroim2 `NeuroVol` mask.
#' @param radius Positive spherical radius in millimeters. Several radii
#'   request a multiscale family, one member frame per radius.
#' @param domain An exact domain from [neuroim2_volume_domain()]. When omitted,
#'   it is constructed from `mask`.
#' @param normalization One of `none`, `local`, or `conservative`. Several
#'   radii admit only `conservative`.
#' @param nonzero Passed to `neuroim2::searchlight_indices()`; version 0.1
#'   requires `TRUE` so every member belongs to the compact domain.
#' @param weights Family weights for a multiscale request: one positive weight
#'   per radius, summing to one, matched to the radii in order or by the
#'   `"radius-<r>"` names. `NULL` (the default) weights the radii equally. A
#'   single radius is one frame with no budget to divide, so `weights` is
#'   refused there rather than ignored.
#' @return An `effect_frame` whose `$weights` are the sparse
#'   measurement-by-feature operator, with `$index$measurement` holding the
#'   full-volume center indices, `$normalization`, a `$specification`
#'   recording the radius and the pinned `upstream_commit`, and a
#'   `$support_index`. Several radii return a [frame_family()] instead: its
#'   `$index` carries one row per measurement with that row's `family`,
#'   `node`, `scale`, `center`, and `alpha`.
#' @section Multiscale families:
#' `neuroim2_searchlights(mask, c(4, 8), normalization = "conservative")`
#' builds one conservative frame per radius from the same neuroim2
#' neighborhoods, names them `"radius-4"` and `"radius-8"`, and stacks them
#' with [frame_family()] under family weights `weights` (equal by default).
#' Only conservative normalization is admitted, for the reason [searchlights()]
#' gives: locally normalized values are not contributions to any total, so a
#' family of them has no budget for `weights` to divide.
#'
#' Per-scale totals of such a family are \eqn{\alpha_s G_\Omega}{alpha_s *
#' G_Omega} by construction, so a total-energy-by-scale panel reports the
#' analyst's own `weights`, not the data. Only the coherent share of each
#' block's fixed budget varies informatively with scale
#' (`design/conservative-geometry-contract.md` sections 3.1 and 3.2).
#' @family neural domains and frames
#' @seealso [searchlights()] plus [compile_frame()] for the built-in
#'   neighborhood builder, [frame_family()] for the family several radii
#'   compile to, [neuroim2_volume_domain()] for the domain, and
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
#'   # Several radii stack into a conservative family, one member per radius.
#'   family <- neuroim2_searchlights(mask, c(4, 6), domain = domain,
#'     normalization = "conservative", weights = c(0.4, 0.6))
#'   print(unique(family$index[, c("family", "scale", "alpha")]))
#'
#'   # Each block carries exactly its weight, so per-scale energy is the
#'   # `weights` vector and only the coherent share is a finding.
#'   print(frame_conservation(family)$members)
#'
#'   # A mask whose geometry differs from the declared domain is refused.
#'   moved <- neuroim2::LogicalNeuroVol(
#'     values, neuroim2::NeuroSpace(c(5L, 5L, 4L), spacing = c(2, 2, 4))
#'   )
#'   print(try(neuroim2_searchlights(moved, radius = 4, domain = domain)))
#' }
#' @export
neuroim2_searchlights <- function(mask, radius, domain = NULL,
                                  normalization = "local", nonzero = TRUE,
                                  weights = NULL) {
  if (missing(mask) || missing(radius)) {
    .input_error(paste0(
      "`mask` and `radius` are both required: pass the `NeuroVol` mask and a ",
      "positive spherical radius in millimetres, for example ",
      "`neuroim2_searchlights(mask, radius = 8)`."
    ))
  }
  .require_neuroim2_searchlight_indices()
  if (!.is_flag(nonzero) || !nonzero) {
    .input_error("crossform neuroim2 searchlights require `nonzero = TRUE`.")
  }
  # A multiscale request is validated by the constructor that owns the rule --
  # `searchlights()` -- so both spatial providers refuse the same things in the
  # same words, and the returned specification carries the scale names and the
  # applied weights. Only the neighborhoods differ between the two providers.
  multiscale <- length(radius) != 1L || !is.null(weights)
  request <- if (multiscale) {
    searchlights(radius, normalization, weights)
  } else {
    NULL
  }
  domain <- .neuroim2_domain_for_mask(mask, domain)
  if (multiscale) {
    members <- lapply(request$radius, function(one) {
      .neuroim2_searchlight_frame(mask, one, domain, "conservative")
    })
    names(members) <- names(request$weights)
    family <- do.call(frame_family,
      c(members, list(alpha = request$weights)))
    family$specification$request <- list(
      kind = "neuroim2_searchlight_family",
      radius = request$radius, weights = request$weights, units = "mm",
      nonzero = TRUE, upstream_commit = "77b1ddb"
    )
    return(family)
  }
  .neuroim2_searchlight_frame(mask, radius, domain, normalization)
}

# One neuroim2 neighborhood frame at one radius, over a domain that has
# already been checked against the mask. Split out of `neuroim2_searchlights()`
# so a multiscale request builds its members without re-resolving the domain
# or re-running the adapter's argument checks once per radius.
.neuroim2_searchlight_frame <- function(mask, radius, domain, normalization) {
  neighborhoods <- neuroim2::searchlight_indices(mask, radius,
    nonzero = TRUE)
  centers <- attr(neighborhoods, "center_indices", exact = TRUE)
  if (!identical(centers, domain$feature_ids)) {
    .contract_error(
      "neuroim2 searchlight centers do not match the ordered domain features."
    )
  }
  members <- lapply(neighborhoods, function(indices) {
    mapped <- match(indices, domain$feature_ids)
    if (anyNA(mapped)) {
      .input_error(
        "A neuroim2 searchlight member lies outside the declared domain."
      )
    }
    as.integer(mapped)
  })
  counts <- lengths(members)
  if (length(members) < 1L || any(counts < 1L)) {
    .input_error(
      "Every neuroim2 searchlight must contain at least one domain feature."
    )
  }
  # The provider's whole contribution is `members`: which compact features
  # each neighborhood covers, and what produced them. The normalization law,
  # the membership operator, and the support bookkeeping belong to the frame
  # constructor, so they are asked for rather than reimplemented here.
  additive_frame(
    members = members, measurements = centers,
    normalization = normalization, domain = domain,
    construction = list(
      kind = "euclidean_ball",
      provider = "neuroim2_searchlight_indices",
      radius = as.numeric(radius),
      coordinate_units = domain$coordinate_units,
      upstream_commit = "77b1ddb"
    ),
    specification = list(kind = "neuroim2_searchlights",
      radius = as.numeric(radius), units = "mm", nonzero = TRUE,
      upstream_commit = "77b1ddb")
  )
}

# `as_neurovol()` is the package's only output adapter, and it is a generic so
# that a package holding its own result type can teach crossform to write that
# type out without crossform importing it. Dispatch is the only thing the
# generic does: a method written for another class may legitimately need no
# mask, so the argument checks belong to crossform's own methods rather than to
# the generic.
.as_neurovol_required_arguments <- function() {
  .input_error(paste0(
    "`values` and `mask` are both required: `as_neurovol()` writes one ",
    "number per compact domain feature onto the space of the `NeuroVol` ",
    "mask the domain was built from."
  ))
}

# The single body behind both shipped methods. `as_neurovol.default()` is not a
# refusal stub: a bare numeric vector is exactly what this function accepted
# before it became generic, and a classed numeric vector reached the same code,
# so the default keeps writing anything numeric and refuses everything else
# with the message it has always raised -- in the same order, so a mask that
# disagrees with its domain is still reported before a type complaint.
.as_neurovol_compact <- function(values, mask, domain, fill, label, ...) {
  dots <- list(...)
  if (length(dots)) {
    .input_error(sprintf(paste0(
      "`as_neurovol()` writes a compact numeric vector and takes only ",
      "`values`, `mask`, `domain`, `fill`, and `label`; received %s. A ",
      "method registered for another class may take more; this one does not."
    ), .msg_count(length(dots), "further argument")))
  }
  .require_neuroim2_searchlight_indices()
  domain <- .neuroim2_domain_for_mask(mask, domain)
  if (!is.numeric(values) || is.matrix(values)) {
    .input_error(sprintf(paste0(
      "`values` must be a numeric vector with one value per compact domain ",
      "feature; received %s. A view's `$values` matrix has one column per ",
      "query, so select the column you want to write out."
    ), .msg_value(values)))
  }
  if (length(values) != domain$n_features) {
    .input_error(sprintf(paste0(
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
      .msg_count(domain$n_features, "feature")))
  }
  if (any(!is.finite(values))) {
    .input_error(sprintf(paste0(
      "`values` must be finite; %d of %d are NA, NaN, or Inf. `as_neurovol()` ",
      "will not guess what a missing measurement means -- replace them ",
      "deliberately, or build the domain from the coverage you actually have."
    ), sum(!is.finite(values)), length(values)))
  }
  if (!is.numeric(fill) || length(fill) != 1L || is.nan(fill) ||
      is.infinite(fill)) {
    .input_error(sprintf(paste0(
      "`fill` must be one numeric value (`NA_real_` is allowed) written ",
      "outside the compact domain; received %s."
    ), .msg_value(fill)))
  }
  if (!.is_string(label, allow_empty = TRUE)) {
    .input_error(sprintf("`label` must be one character string; received %s.",
      .msg_value(label)))
  }
  payload <- array(as.double(fill), dim = dim(mask))
  payload[domain$feature_ids] <- as.double(values)
  neuroim2::NeuroVol(payload, neuroim2::space(mask), label = label)
}

#' Map a compact result vector back to a neuroim2 volume
#'
#' The compact values are inserted at the exact full-volume indices carried by
#' a crossform volume domain. Features outside the domain receive `fill`.
#' This is an output adapter only; it performs no interpolation, smoothing, or
#' coordinate reinterpretation.
#'
#' @details
#' `as_neurovol()` is an S3 generic dispatching on `values`, so a package that
#' owns its own result type can write that type out without crossform having to
#' know about it. Register a method the ordinary way --- `S3method(as_neurovol,
#' my_result)` in your NAMESPACE --- and it receives `mask`, `domain`, `fill`,
#' and `label` unchanged; it is expected to return a `NeuroVol` on the space of
#' `mask`. The generic validates nothing itself, so a method is free to require
#' different arguments, or none beyond the object.
#'
#' crossform ships the numeric-vector method described here. The default method
#' behaves identically, so any numeric vector still writes out whether or not it
#' carries a class, and anything that is not a numeric vector is refused.
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
#' @param ... Arguments passed on to methods. The methods crossform ships take
#'   no further arguments and refuse any.
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
#'
#'   # The generic is the extension point: a package with its own result type
#'   # registers a method for it and delegates the writing back here.
#'   as_neurovol.crossform_example_map <- function(values, mask, ...) {
#'     as_neurovol(values$statistic, mask, ...)
#'   }
#'   boxed <- structure(list(statistic = statistic),
#'     class = "crossform_example_map")
#'   print(identical(
#'     as.numeric(as_neurovol(boxed, mask, domain)[domain$feature_ids]),
#'     statistic
#'   ))
#' }
#' @export
as_neurovol <- function(values, ...) {
  if (missing(values)) .as_neurovol_required_arguments()
  UseMethod("as_neurovol")
}

#' @rdname as_neurovol
#' @export
as_neurovol.numeric <- function(values, mask, domain = NULL, fill = NA_real_,
                                label = "crossform result", ...) {
  if (missing(mask)) .as_neurovol_required_arguments()
  .as_neurovol_compact(values, mask, domain, fill, label, ...)
}

#' @rdname as_neurovol
#' @export
as_neurovol.default <- function(values, mask, domain = NULL, fill = NA_real_,
                                label = "crossform result", ...) {
  if (missing(mask)) .as_neurovol_required_arguments()
  .as_neurovol_compact(values, mask, domain, fill, label, ...)
}
