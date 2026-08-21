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
#' Declares one measurement per neural feature. This is the finest scope and,
#' under the default conservative normalization, the only spatial scope whose
#' local values sum exactly to the whole-brain value.
#'
#' @param normalization Explicit frame normalization.
#' @return An `effect_frame_spec` with `$kind` `"voxels"` and the requested
#'   `$normalization`. Pass it to [compile_frame()] to obtain measurements.
#' @family neural domains and frames
#' @seealso [compile_frame()], and [searchlights()], [regions()], or
#'   [whole_brain()] for coarser scopes.
#' @examples
#' domain <- abstract_domain(4L, id = "voxelwise-example")
#'
#' # One measurement per feature: the operator is the identity.
#' frame <- compile_frame(voxelwise(), domain)
#' dim(frame$weights)
#'
#' # Conservative normalization is the default here, so local totals add up to
#' # the whole-brain total exactly.
#' frame_conservation(frame)$conserved
#' @export
voxelwise <- function(normalization = "conservative") {
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
#' Several radii request one conservative frame per radius, stacked into a
#' multiscale [frame_family()]; see *Multiscale families*.
#'
#' @param radius Positive radius in domain coordinate units. Several radii
#'   request a multiscale family, one member frame per radius.
#' @param normalization Explicit frame normalization. Several radii admit only
#'   `"conservative"`.
#' @param weights Family weights for a multiscale request: one positive weight
#'   per radius, summing to one, matched to the radii in order or by the
#'   `"radius-<r>"` names. `NULL` (the default) weights the radii equally. A
#'   single radius is one frame with no budget to divide, so `weights` is
#'   refused there rather than ignored.
#' @return An `effect_frame_spec`. With one radius its `$kind` is
#'   `"searchlights"` and it carries the requested `$radius` and
#'   `$normalization`, exactly as before. With several its `$kind` is
#'   `"searchlight_family"` and it carries the `$radius` vector and the
#'   `$weights` that will be applied. Pass either to [compile_frame()].
#' @section Multiscale families:
#' `compile_frame(searchlights(c(4, 8, 12), "conservative"), domain)` returns
#' a [frame_family()]: one conservative member frame per radius, named
#' `"radius-4"`, `"radius-8"`, `"radius-12"`, stacked with family weights
#' `weights` (equal by default). Every row of the compiled frame's `$index`
#' carries its own `family`, `scale` (the radius it came from), `center`, and
#' `alpha`, so a result can be grouped by scale after the fact.
#'
#' Only `normalization = "conservative"` admits several radii. Locally
#' normalized values are not contributions to any total, so a family of them
#' has no budget for `weights` to divide and the per-scale law below is
#' undefined (`design/conservative-geometry-contract.md` section 3.1).
#'
#' **What a multiscale family can and cannot show.** Because every member is
#' column normalized and the weights sum to one, the family conserves block by
#' block: the `total` component summed over the rows of scale `s` is exactly
#' \eqn{\alpha_s G_\Omega}{alpha_s * G_Omega}, the scale's weight times the
#' whole-brain total, whatever the data say. A panel of *total energy by
#' scale* is therefore a plot of the analyst's own `weights` vector and is not
#' a finding. What does vary informatively with scale is the split of each
#' block's fixed budget into coherent and configuration parts: the coherent
#' share is invariant to `alpha` and is the scale-resolved quantity a
#' multiscale family exists to report
#' (`design/conservative-geometry-contract.md` sections 3.1 and 3.2).
#' @family neural domains and frames
#' @seealso [compile_frame()], [frame_conservation()] to check the
#'   normalization you chose, [frame_family()] for the family a multiscale
#'   request compiles to, and [neuroim2_searchlights()] for neighborhoods
#'   built by neuroim2 instead.
#' @examples
#' grid <- as.matrix(expand.grid(x = 1:4, y = 1:4))
#' domain <- abstract_domain(
#'   nrow(grid), coordinates = grid, id = "searchlight-example"
#' )
#'
#' # One center-assigned neighborhood per feature; the radius is read in the
#' # domain's own coordinate units.
#' local <- compile_frame(searchlights(1.5), domain)
#' dim(local$weights)
#'
#' # The default "local" normalization double-counts features shared between
#' # overlapping neighborhoods, so its local values are not contributions to a
#' # whole-brain total. Ask for "conservative" when they must be.
#' frame_conservation(local)$conserved
#' frame_conservation(
#'   compile_frame(searchlights(1.5, normalization = "conservative"), domain)
#' )$conserved
#'
#' # Several radii request a multiscale family instead: one conservative
#' # member frame per radius, each row labelled with the scale it came from.
#' family <- compile_frame(
#'   searchlights(c(1.5, 2.5), "conservative", weights = c(0.25, 0.75)), domain
#' )
#' unique(family$index[, c("family", "scale", "alpha")])
#'
#' # Each block carries exactly its weight, whatever the data: per-scale
#' # energy is the `weights` vector, so only the coherent share is a finding.
#' frame_conservation(family)$members
#' @export
searchlights <- function(radius, normalization = "local", weights = NULL) {
  if (missing(radius)) {
    .input_error(paste0(
      "`radius` is required: pass one positive radius in the domain's own ",
      "coordinate units, for example `searchlights(8)` for 8 mm on a volume ",
      "domain."
    ),
      arg = "radius", received = "no argument",
      expected = "one positive radius in the domain's coordinate units")
  }
  if (length(radius) == 1L) {
    if (!.is_number(radius) || radius <= 0) {
      .input_error(sprintf(paste0(
        "`radius` must be one positive finite number in the domain's ",
        "coordinate units; received %s."
      ), if (is.numeric(radius) && length(radius) == 1L) {
        paste0("`", format(radius), "`")
      } else {
        .msg_value(radius)
      }),
        arg = "radius", received = .msg_value(radius),
        expected = "one positive finite number")
    }
    .reject_single_scale_weights(weights)
    return(.frame_spec("searchlights", normalization, radius = radius))
  }
  if (!is.numeric(radius) || length(radius) < 1L ||
      any(!is.finite(radius)) || any(radius <= 0)) {
    .input_error(sprintf(paste0(
      "`radius` must be one positive finite number in the domain's ",
      "coordinate units, or several of them for a multiscale family; ",
      "received %s."
    ), .msg_value(radius)),
      arg = "radius", received = .msg_value(radius),
      expected = "one or more positive finite numbers")
  }
  .require_multiscale_normalization(normalization, length(radius))
  scales <- .searchlight_scale_names(radius)
  .reject_indistinct_scales(scales, radius)
  weights <- .frame_family_alpha(weights, scales, 1e-12, arg = "weights",
    unit = "radius", plural = "radii", subject = "radius")
  .frame_spec("searchlight_family", "conservative",
    radius = as.numeric(radius), weights = weights)
}

# The multiscale rules, shared by `searchlights()` and, through it,
# `neuroim2_searchlights()`. They are stated without naming a constructor
# because both entry points raise them and the argument names are the same.

# A family name per radius. It is what `$index$family` reports and what named
# `weights` are matched against, so it must be stable and readable: 4 prints
# as "radius-4", not "radius-4.000000".
.searchlight_scale_names <- function(radius) {
  paste0("radius-", vapply(as.numeric(radius), function(value) {
    format(value, trim = TRUE, digits = 15)
  }, character(1)))
}

.reject_single_scale_weights <- function(weights) {
  if (is.null(weights)) return(invisible(NULL))
  .input_error(sprintf(paste0(
    "`weights` divides a multiscale family's conserved budget between its ",
    "radii, so it applies only when `radius` names more than one scale; ",
    "received one radius and %s. One radius is one frame, which already ",
    "carries the whole budget, so the weight is refused rather than ignored: ",
    "drop `weights`, or ask for several radii."
  ), .msg_count(length(weights), "weight")),
    arg = "weights", received = "weights alongside a single radius",
    expected = "NULL, or several radii")
}

.require_multiscale_normalization <- function(normalization, count) {
  if (.is_string(normalization) && identical(normalization, "conservative")) {
    return(invisible(NULL))
  }
  received <- if (.is_string(normalization)) {
    paste0("\"", normalization, "\"")
  } else {
    .msg_value(normalization)
  }
  .input_error(sprintf(paste0(
    "A multiscale request needs `normalization = \"conservative\"`; %s ",
    "asked for %s. The family law sum over the rows of scale s of G equals ",
    "alpha_s times G_Omega is defined only for column-normalized members, so ",
    "locally normalized scales have no conserved budget for `weights` to ",
    "divide (`design/conservative-geometry-contract.md` section 3.1). Pass ",
    "`normalization = \"conservative\"`, or one radius."
  ), .msg_count(count, "radius", "radii"), received),
    arg = "normalization", received = received,
    expected = "\"conservative\"")
}

.reject_indistinct_scales <- function(scales, radius) {
  if (!anyDuplicated(scales)) return(invisible(NULL))
  duplicated_scale <- scales[[anyDuplicated(scales)]]
  .input_error(sprintf(paste0(
    "Every radius must name a distinct scale, but `%s` names more than one ",
    "of the %s supplied. The scale name is each member's `family` identity, ",
    "so radii that print alike cannot be told apart in `$index`."
  ), duplicated_scale, .msg_count(length(radius), "radius", "radii")),
    arg = "radius", received = sprintf("duplicated scale `%s`",
      duplicated_scale),
    expected = "distinct radii")
}

#' Specify region measurements
#'
#' Declares one measurement per distinct label, so an atlas or region-label
#' vector becomes a set of non-overlapping regional measurements.
#'
#' @param labels One region label per neural feature. Missing labels are
#'   excluded unless conservative normalization is requested.
#' @param normalization Explicit frame normalization.
#' @return An `effect_frame_spec` with `$kind` `"regions"`, the supplied
#'   `$labels`, and `$normalization`. Pass it to [compile_frame()].
#' @family neural domains and frames
#' @seealso [compile_frame()], and [whole_brain()] for the single-region case.
#' @examples
#' domain <- abstract_domain(12L, id = "region-example")
#' labels <- rep(paste0("roi-", 1:3), each = 4L)
#'
#' # Measurement order follows first appearance of each label.
#' frame <- compile_frame(regions(labels), domain)
#' frame$index$measurement
#'
#' # Unlabeled features are dropped: with "none" the row sums are the member
#' # counts, so the effect of an NA label is visible.
#' as.numeric(Matrix::rowSums(
#'   compile_frame(regions(labels, normalization = "none"), domain)$weights
#' ))
#' labels[[1L]] <- NA
#' as.numeric(Matrix::rowSums(
#'   compile_frame(regions(labels, normalization = "none"), domain)$weights
#' ))
#' @export
regions <- function(labels, normalization = "local") {
  if (missing(labels)) {
    .input_error(paste0(
      "`labels` is required: pass one region label per neural feature, in ",
      "domain feature order, for example `regions(atlas_labels)`."
    ),
      arg = "labels", received = "no argument",
      expected = "one region label per neural feature")
  }
  if (length(labels) < 1L || !(is.atomic(labels) || is.factor(labels))) {
    .input_error(sprintf(paste0(
      "`labels` must be a nonempty atomic or factor vector with one region ",
      "label per neural feature; received %s."
    ), .msg_value(labels)),
      arg = "labels", received = .msg_value(labels),
      expected = "a nonempty atomic or factor vector, one label per feature")
  }
  .frame_spec("regions", normalization, labels = labels)
}

#' Specify a whole-brain additive measurement
#'
#' Declares a single measurement covering every neural feature in the domain.
#' Use `normalization = "none"` when this is the global comparator against
#' which local measurements are checked, because the default averages instead
#' of summing.
#'
#' @param normalization Explicit frame normalization.
#' @return An `effect_frame_spec` with `$kind` `"whole_brain"` and the
#'   requested `$normalization`. Pass it to [compile_frame()].
#' @family neural domains and frames
#' @seealso [compile_frame()], [frame_conservation()], and [regions()] for
#'   several measurements instead of one.
#' @examples
#' domain <- abstract_domain(6L, id = "whole-brain-example")
#'
#' # A single measurement row spanning the whole domain.
#' frame <- compile_frame(whole_brain(), domain)
#' dim(frame$weights)
#'
#' # The default "local" normalization averages over features (mass 1); the
#' # unnormalized operator is the one that sums, and is the correct global
#' # comparator when checking local-to-global conservation.
#' as.numeric(Matrix::rowSums(frame$weights))
#' as.numeric(Matrix::rowSums(
#'   compile_frame(whole_brain("none"), domain)$weights
#' ))
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
#' @return An `effect_frame`: a list whose `$weights` are the sparse
#'   measurement-by-feature operator, with `$normalization`, an `$index` data
#'   frame naming each measurement, `$domain` (the domain reference),
#'   `$domain_kind`, the originating `$specification`, and, for neighborhood
#'   scopes, a `$support_index`.
#' @section Structure:
#' A compiled frame is the spatial operator together with the record of how
#' it was built.
#'
#' - `$weights`: the sparse measurement-by-feature operator. Row `m` holds
#'   the weight each domain feature contributes to measurement `m`, in domain
#'   feature order, after normalization was applied.
#' - `$index`: one row per measurement, in `$weights` row order. Its
#'   `measurement` column names each measurement: domain feature identifiers
#'   for [voxelwise()] and [searchlights()], the distinct labels in first
#'   appearance order for [regions()], and `"whole_brain"` for
#'   [whole_brain()]. Views carry these identifiers through as their `$index`.
#' - `$normalization`: the normalization that was applied, which is what
#'   [frame_conservation()] reports against.
#' - `$specification`: the `effect_frame_spec` the frame was compiled from,
#'   so the scope and its arguments travel with the operator.
#' - `$domain`: the neural domain the columns are bound to, carrying its
#'   `id`, `n_features`, and `feature_ids`.
#'
#' Any element not listed here, including `$support_index` and the
#' representation flags, is internal and may change.
#'
#' One specification is not one operator: a multiscale [searchlights()]
#' request compiles to a [frame_family()], whose `$index` carries the extra
#' per-row `family`, `node`, `scale`, `center`, and `alpha` columns and whose
#' `$specification` records every member rather than one scope.
#' @family neural domains and frames
#' @seealso [voxelwise()], [searchlights()], [regions()], [whole_brain()] for
#'   the specifications, [frame_conservation()] to check normalization, and
#'   [plan_geometry()], which consumes the compiled frame.
#' @examples
#' domain <- abstract_domain(9L, id = "compile-frame-example")
#'
#' # Every scope compiles to one sparse measurement-by-feature operator; only
#' # the number of measurement rows differs.
#' scopes <- list(
#'   voxelwise(), regions(rep(c("a", "b", "c"), each = 3L)), whole_brain()
#' )
#' vapply(scopes, function(scope) nrow(compile_frame(scope, domain)$weights),
#'   integer(1))
#'
#' # The compiled frame remembers how it was built, which is what downstream
#' # receipts record.
#' frame <- compile_frame(voxelwise(), domain)
#' frame$normalization
#' frame$index$measurement
#' @export
compile_frame <- function(specification, domain) {
  if (missing(domain) || !inherits(domain, "effect_domain")) {
    .input_error(sprintf(paste0(
      "`domain` must be an `effect_domain` (see `abstract_domain()`, ",
      "`volume_domain()`, or `neuroim2_volume_domain()`); received %s."
    ), if (missing(domain)) "no argument" else .msg_value(domain)),
      arg = "domain",
      received = if (missing(domain)) "no argument" else .msg_value(domain),
      expected = "an `effect_domain`")
  }
  .validate_domain(domain)
  if (missing(specification) || !inherits(specification, "effect_frame_spec")) {
    .input_error(sprintf(paste0(
      "`specification` must be a frame specification from `voxelwise()`, ",
      "`searchlights()`, `regions()`, or `whole_brain()`; received %s."
    ), if (missing(specification)) "no argument" else
      .msg_value(specification)),
      arg = "specification",
      received = if (missing(specification)) "no argument" else
        .msg_value(specification),
      expected = paste("a frame specification from `voxelwise()`,",
        "`searchlights()`, `regions()`, or `whole_brain()`"))
  }
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
      .input_error(sprintf(paste0(
        "`regions()` supplied %s but the domain has %s. Labels are read in ",
        "domain feature order, one per feature."
      ), .msg_count(length(labels), "label"),
        .msg_count(n, "feature")),
        arg = "specification",
        received = .msg_count(length(labels), "label"),
        expected = .msg_count(n, "label"))
    }
    present <- !is.na(labels) & nzchar(as.character(labels))
    if (!any(present)) {
      .input_error(paste0(
        "Every region label is missing or empty, so the frame would have no ",
        "measurements. At least one feature must carry a label."
      ),
        arg = "specification", received = "every label missing or empty",
        expected = "at least one feature carrying a label")
    }
    region_ids <- unique(as.character(labels[present]))
    region_index <- match(as.character(labels[present]), region_ids)
    weights <- Matrix::sparseMatrix(
      i = region_index, j = which(present), x = 1,
      dims = c(length(region_ids), n)
    )
    index <- data.frame(measurement = region_ids, stringsAsFactors = FALSE)
  } else if (kind == "searchlights") {
    radius <- specification$radius
    if (!.is_number(radius) || radius <= 0) {
      .input_error("Searchlight radius is invalid.")
    }
    support_index <- .euclidean_support_index(domain, radius)
    weights <- .support_index_membership(support_index)
    index <- data.frame(measurement = domain$feature_ids,
      stringsAsFactors = FALSE)
  } else if (kind == "searchlight_family") {
    # A multiscale request is not one operator: it compiles to one member
    # frame per radius, stacked by `frame_family()`, so it leaves here rather
    # than falling through to the single-operator normalization below.
    return(.compile_searchlight_family(specification, domain))
  } else {
    .input_error("Unknown additive frame specification.")
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
      !is.list(specification) ||
      !.is_string(specification$kind, allow_empty = TRUE) ||
      !.is_string(specification$normalization, allow_empty = TRUE)) {
    .input_error(
      "`specification` must be a valid additive frame specification."
    )
  }
  rebuilt <- switch(specification$kind,
    voxels = voxelwise(specification$normalization),
    searchlights = {
      if (!identical(names(specification), c("kind", "normalization", "radius"))) {
        .input_error("Frame specification fields are missing or noncanonical.")
      }
      searchlights(specification$radius, specification$normalization)
    },
    searchlight_family = {
      if (!identical(names(specification),
          c("kind", "normalization", "radius", "weights"))) {
        .input_error("Frame specification fields are missing or noncanonical.")
      }
      searchlights(specification$radius, specification$normalization,
        specification$weights)
    },
    regions = {
      if (!identical(names(specification), c("kind", "normalization", "labels"))) {
        .input_error("Frame specification fields are missing or noncanonical.")
      }
      regions(specification$labels, specification$normalization)
    },
    whole_brain = whole_brain(specification$normalization),
    .input_error("Unknown additive frame specification.")
  )
  if (!identical(specification, rebuilt)) {
    .input_error("Frame specification fields are missing or noncanonical.")
  }
  rebuilt
}

# One conservative member frame per radius, stacked with the requested family
# weights. The members are ordinary single-radius frames, so a multiscale
# family is exactly what an analyst could have assembled by hand -- what the
# constructor adds is that the weights were checked before any frame was
# built, and that the originating request survives on `$specification$request`
# where the stacked specification would otherwise remember only the members.
.compile_searchlight_family <- function(specification, domain) {
  scales <- names(specification$weights)
  members <- lapply(specification$radius, function(radius) {
    compile_frame(searchlights(radius, "conservative"), domain)
  })
  names(members) <- scales
  family <- do.call(frame_family,
    c(members, list(alpha = specification$weights)))
  family$specification$request <- specification
  family
}

# Frame families ------------------------------------------------------------

#' Combine conservative frames into one alpha-weighted family
#'
#' Stacks several conservative frames over one domain into a single compiled
#' frame whose rows are the members' rows scaled by family weights `alpha`.
#' Because each member is column-normalized on its own and the weights sum to
#' one, the stacked columns still sum to one, so the family is conservative and
#' the `total` component conserves both overall and block by block:
#' \deqn{\sum_{x \in s} G_{s,x} = \alpha_s\,G_\Omega,\qquad
#'       \sum_x G_x = G_\Omega.}{sum_{x in s} G_{s,x} = alpha_s G_Omega, and sum_x G_x = G_Omega.}
#'
#' The consequence is that per-scale *energy* is fixed by `alpha` alone and is
#' never a finding. What varies with the data is the split of each block's
#' fixed budget into coherent and configuration parts. See
#' `design/conservative-geometry-contract.md` section 3.
#'
#' @param ... Two or more compiled `effect_frame`s over one domain, ideally
#'   named. A member's name becomes its `family` identity; unnamed members are
#'   named `frame1`, `frame2`, and so on by position. Pass a list of frames
#'   with `do.call(frame_family, c(frames, list(alpha = alpha)))`.
#' @param alpha One positive family weight per member, summing to one. Named
#'   weights are matched to member names; `NULL` (the default) weights the
#'   members equally. Weights are never renormalized for you, because the
#'   per-block law reads against the weight actually applied.
#' @param normalization Frame normalization of the family. Only
#'   `"conservative"` is defined: the per-block law needs column-normalized
#'   members, and a family of locally normalized frames has no conserved
#'   budget to weight.
#' @param tolerance Nonnegative absolute tolerance for the two construction
#'   checks: that `alpha` sums to one, and that each member's columns sum to
#'   one on its own.
#' @return An `effect_frame` usable anywhere a compiled frame is, carrying
#'   `$weights` (the alpha-scaled row-bind), `$index` (one row per
#'   measurement, see *Per-row metadata*), and a `$specification` recording
#'   every member specification together with its applied weight.
#' @section Per-row metadata:
#' A family row must be self-describing, because its own scale and provenance
#' can no longer be read off a frame-wide specification. `$index` therefore
#' carries one row per measurement, in `$weights` row order:
#'
#' - `measurement`: the row's identity, `"<family>::<node>"`. It is unique
#'   across the family, which the same node label appearing at several scales
#'   is not, and it is what reaches a result's `$index`.
#' - `family`: the member the row came from.
#' - `node`: the row's label inside its own member, exactly as that member's
#'   `$index$measurement` had it.
#' - `scale`: the member's scale parameter -- the radius for a searchlight
#'   member, `NA` for a member that has no scale.
#' - `center`: the anchor feature identifier, for members whose rows are
#'   anchored at a feature (points and searchlights); `NA` otherwise.
#' - `alpha`: the family weight applied to the row.
#'
#' Join a result back to this table by `measurement` to group values by scale,
#' by center, or by member.
#' @family neural domains and frames
#' @seealso [compile_frame()] for the members, [frame_conservation()], which
#'   certifies a family both overall and block by block, and
#'   [searchlights()], whose multiscale form is the shorthand for a family of
#'   conservative searchlight frames at several radii.
#' @examples
#' domain <- abstract_domain(
#'   9L, coordinates = cbind(seq_len(9L) - 1, 0), id = "frame-family-example"
#' )
#'
#' # Two conservative scales over one domain, weighted one quarter and three
#' # quarters.
#' family <- frame_family(
#'   point = compile_frame(voxelwise("conservative"), domain),
#'   narrow = compile_frame(searchlights(1.01, "conservative"), domain),
#'   alpha = c(point = 0.25, narrow = 0.75)
#' )
#' head(family$index, 3L)
#'
#' # The stacked family conserves, and each block carries exactly its alpha.
#' report <- frame_conservation(family)
#' report$conserved
#' report$members
#'
#' # Weights that do not sum to one are refused rather than renormalized: the
#' # per-block law reads against the weight actually applied.
#' refused <- try(
#'   frame_family(
#'     point = compile_frame(voxelwise("conservative"), domain),
#'     narrow = compile_frame(searchlights(1.01, "conservative"), domain),
#'     alpha = c(1, 1)
#'   ),
#'   silent = TRUE
#' )
#' conditionMessage(attr(refused, "condition"))
#' @export
frame_family <- function(..., alpha = NULL, normalization = "conservative",
                         tolerance = 1e-12) {
  members <- list(...)
  if (!length(members)) {
    .input_error(paste0(
      "`frame_family()` needs at least one compiled frame: pass the members ",
      "as named arguments, for example `frame_family(point = ..., ",
      "narrow = ...)`."
    ),
      arg = "...", received = "no frames",
      expected = "one or more compiled `effect_frame`s")
  }
  if (!.is_string(normalization) ||
      !identical(normalization, "conservative")) {
    received <- if (.is_string(normalization)) {
      paste0("\"", normalization, "\"")
    } else {
      .msg_value(normalization)
    }
    .input_error(sprintf(paste0(
      "`normalization` must be \"conservative\"; received %s. The family law ",
      "sum over the rows of scale s of G equals alpha_s times G_Omega is ",
      "defined only for column-normalized members, so a family of locally ",
      "normalized frames has no conserved budget for `alpha` to divide."
    ), received),
      arg = "normalization", received = received,
      expected = "\"conservative\"")
  }
  .check_number(tolerance, "tolerance", nonnegative = TRUE)
  names(members) <- .frame_family_names(names(members), length(members))
  for (position in seq_along(members)) {
    .validate_frame_family_member(members[[position]], names(members)[[position]],
      tolerance)
  }
  domain <- .frame_family_domain(members)
  alpha <- .frame_family_alpha(alpha, names(members), tolerance)

  blocks <- Map(function(member, weight) {
    weight * .frame_family_block(member$weights)
  }, members, alpha)
  weights <- if (length(blocks) == 1L) blocks[[1L]] else do.call(rbind, blocks)
  column_mass <- as.numeric(Matrix::colSums(weights))
  deviation <- max(abs(column_mass - 1))
  if (!all(is.finite(column_mass)) || deviation > 1e-12) {
    .input_error(sprintf(paste0(
      "The stacked family columns deviate from unit mass by %s, so the ",
      "family is not conservative. A conservative frame's columns must sum to ",
      "one within 1e-12 whatever `tolerance` the per-member checks were given: ",
      "every member must be column normalized on its own, and `alpha` must sum ",
      "to one."
    ), format(deviation, digits = 3)),
      arg = "...", received = sprintf("column mass off by %s",
        format(deviation, digits = 3)),
      expected = "columns summing to one")
  }

  index <- do.call(rbind, Map(.frame_family_member_index, members,
    names(members), alpha))
  rownames(index) <- NULL
  if (anyDuplicated(index$measurement)) {
    duplicated_label <- index$measurement[[anyDuplicated(index$measurement)]]
    .input_error(sprintf(paste0(
      "Family measurement identifiers must be unique, but `%s` appears more ",
      "than once. Identifiers are `\"<family>::<node>\"`, so give the members ",
      "distinct names, or the offending member distinct node labels."
    ), duplicated_label),
      arg = "...", received = sprintf("duplicated identifier `%s`",
        duplicated_label),
      expected = "one identifier per measurement")
  }

  frame <- additive_frame(weights, normalization = "conservative",
    domain = domain)
  frame$index <- index
  domain_kind <- .frame_family_domain_kind(members)
  if (!is.null(domain_kind)) frame$domain_kind <- domain_kind
  frame$specification <- .frame_family_specification(members, alpha)
  .validate_frame_for_compile(frame)
  frame
}

.frame_family_names <- function(supplied, count) {
  names <- if (is.null(supplied)) rep("", count) else supplied
  names[is.na(names)] <- ""
  blank <- !nzchar(names)
  names[blank] <- paste0("frame", seq_len(count))[blank]
  if (anyDuplicated(names)) {
    duplicated_name <- names[[anyDuplicated(names)]]
    .input_error(sprintf(paste0(
      "Family member names must be unique, but `%s` names more than one ",
      "member. The name is the row's `family` identity, so it has to ",
      "distinguish the members it labels."
    ), duplicated_name),
      arg = "...", received = sprintf("duplicated name `%s`", duplicated_name),
      expected = "one name per member")
  }
  names
}

.validate_frame_family_member <- function(member, name, tolerance) {
  if (!inherits(member, "effect_frame")) {
    .input_error(sprintf(paste0(
      "Family member `%s` must be a compiled `effect_frame` from ",
      "`compile_frame()`; received %s. Frame weights and `alpha` are ",
      "separate arguments: pass the frames in `...` and the weights in ",
      "`alpha`."
    ), name, .msg_value(member)),
      arg = name, received = .msg_value(member),
      expected = "a compiled `effect_frame`")
  }
  .validate_frame_for_compile(member)
  if (!identical(member$representation, "additive_diagonal")) {
    .input_error(sprintf(paste0(
      "Family member `%s` uses the `%s` representation; only additive ",
      "diagonal frames stack into a family."
    ), name, member$representation),
      arg = name, received = member$representation,
      expected = "an additive diagonal frame")
  }
  if (!is.null(member$metric_folded)) {
    .input_error(sprintf(paste0(
      "Family member `%s` has a diagonal metric folded into its weights, so ",
      "its columns carry the metric diagonal rather than unit mass and the ",
      "per-block law has no fixed budget to divide. Build the family from ",
      "declared frames and supply the metric to `plan_geometry()` instead."
    ), name),
      arg = name, received = "a metric-folded frame",
      expected = "a frame whose columns sum to one")
  }
  # Gap G2 of the contract: the stacked family can conserve while no single
  # block does, so column normalization is checked per member and never on the
  # stack. A member that leaves a feature uncovered cannot be column
  # normalized at all, and is refused here rather than absorbed.
  column_mass <- as.numeric(Matrix::colSums(member$weights))
  deviation <- max(abs(column_mass - 1))
  if (!all(is.finite(column_mass)) || deviation > tolerance) {
    uncovered <- sum(column_mass <= 0)
    .input_error(sprintf(paste0(
      "Family member `%s` is not column normalized on its own: its per-",
      "feature mass is off by %s%s. The per-block law needs every member to ",
      "partition the whole domain separately, which a conserving stack does ",
      "not imply. Compile it with `normalization = \"conservative\"`."
    ), name, format(deviation, digits = 3),
      if (uncovered > 0L) {
        sprintf(", and %s carry no mass at all",
          .msg_count(uncovered, "feature"))
      } else {
        ""
      }),
      arg = name,
      received = sprintf("per-feature mass off by %s",
        format(deviation, digits = 3)),
      expected = "a member whose columns sum to one")
  }
  invisible(member)
}

.frame_family_domain <- function(members) {
  domain <- .domain_reference(members[[1L]]$domain)
  for (position in seq_along(members)) {
    other <- .domain_reference(members[[position]]$domain)
    if (!.same_domain_reference(other, domain)) {
      .contract_error(sprintf(paste0(
        "Family member `%s` is bound to neural domain `%s`, but member `%s` ",
        "is bound to `%s`. A family stacks rows over one domain, so every ",
        "member must be compiled against the same one."
      ), names(members)[[position]], other$id, names(members)[[1L]],
        domain$id),
        arg = names(members)[[position]], received = other$id,
        expected = domain$id)
    }
  }
  domain
}

# One weight per member, checked before anything is stacked. The argument is
# `alpha` on `frame_family()` and `weights` on a multiscale `searchlights()`
# request, and the members are family members in one case and radii in the
# other, so the four naming arguments exist to keep both refusals literal
# about what the caller actually typed. The rules themselves are identical,
# and G1 of the contract requires that they be: one applied weight per member,
# positive, summing to one, never renormalized.
.frame_family_alpha <- function(alpha, names, tolerance, arg = "alpha",
                                unit = "family member", plural = "members",
                                subject = "member") {
  count <- length(names)
  if (is.null(alpha)) alpha <- rep(1 / count, count)
  if (!is.numeric(alpha) || length(alpha) != count ||
      any(!is.finite(alpha))) {
    .input_error(sprintf(paste0(
      "`%s` must be %s, one per %s, and every weight finite; received %s."
    ), arg, .msg_count(count, "finite number"), unit, .msg_value(alpha)),
      arg = arg, received = .msg_value(alpha),
      expected = .msg_count(count, "finite weight"))
  }
  if (!is.null(names(alpha))) {
    if (!setequal(names(alpha), names) || anyDuplicated(names(alpha))) {
      .input_error(sprintf(paste0(
        "Named `%s` must name every %s exactly once. The %s are %s; `%s` ",
        "names %s."
      ), arg, unit, plural, .frame_family_name_list(names), arg,
        .frame_family_name_list(names(alpha))),
        arg = arg, received = .frame_family_name_list(names(alpha)),
        expected = .frame_family_name_list(names))
    }
    alpha <- alpha[names]
  }
  if (any(alpha <= 0)) {
    .input_error(sprintf(paste0(
      "Every family weight must be positive; `%s` holds %s. A nonpositive ",
      "weight does not remove a %s: it contributes rows of zero or ",
      "negative mass, which is not a share of the budget. Drop the %s ",
      "from the family instead."
    ), arg, format(min(alpha)), subject, subject),
      arg = arg, received = sprintf("a weight of %s", format(min(alpha))),
      expected = "positive weights")
  }
  total <- sum(alpha)
  if (abs(total - 1) > tolerance) {
    .input_error(sprintf(paste0(
      "Family weights must sum to one; `%s` sums to %s, off by %s. The ",
      "per-block law reads the weight that was actually applied, so `%s` ",
      "is never renormalized for you: pass weights summing to one, for ",
      "example `%s / sum(%s)`."
    ), arg, format(total, digits = 12), format(abs(total - 1), digits = 3),
      arg, arg, arg),
      arg = arg, received = sprintf("weights summing to %s",
        format(total, digits = 12)),
      expected = "weights summing to one")
  }
  stats::setNames(as.numeric(alpha), names)
}

.frame_family_name_list <- function(names) {
  paste0("`", names, "`", collapse = ", ")
}

.frame_family_block <- function(weights) {
  if (!inherits(weights, "Matrix")) {
    weights <- Matrix::Matrix(weights, sparse = TRUE)
  }
  weights <- methods::as(weights, "dMatrix")
  methods::as(methods::as(weights, "generalMatrix"), "CsparseMatrix")
}

# One row of per-row metadata per measurement of one member. A row's scale and
# center used to be readable only from the frame-wide `$specification`, which
# stacking destroys; carrying them here is what lets a transported or grouped
# result say which instrument produced a number.
.frame_family_member_index <- function(member, family, alpha) {
  count <- nrow(member$weights)
  specification <- member$specification
  node <- .frame_family_nodes(member, count)
  scale <- if (.is_number(specification$radius)) {
    as.numeric(specification$radius)
  } else {
    NA_real_
  }
  # A center exists exactly when the row is anchored at a feature: a
  # neighborhood records its anchor in `$support_index$node_ids`, and a point
  # frame's row *is* one feature. A region or whole-brain row has none.
  center <- if (!is.null(member$support_index)) {
    as.character(member$support_index$node_ids)
  } else if (identical(specification$kind, "voxels")) {
    node
  } else {
    rep(NA_character_, count)
  }
  data.frame(
    measurement = paste0(family, "::", node),
    family = rep(family, count),
    node = node,
    scale = rep(scale, count),
    center = center,
    alpha = rep(as.numeric(alpha), count),
    stringsAsFactors = FALSE
  )
}

.frame_family_nodes <- function(member, count) {
  index <- member$index
  if (is.data.frame(index) && "measurement" %in% names(index) &&
      nrow(index) == count) {
    return(as.character(index$measurement))
  }
  as.character(seq_len(count))
}

.frame_family_domain_kind <- function(members) {
  kinds <- unique(unlist(lapply(members, `[[`, "domain_kind")))
  if (length(kinds) != 1L) NULL else kinds
}

.frame_family_specification <- function(members, alpha) {
  records <- Map(function(member, family, weight) {
    list(
      family = family,
      alpha = as.numeric(weight),
      scale = if (.is_number(member$specification$radius)) {
        as.numeric(member$specification$radius)
      } else {
        NA_real_
      },
      measurements = nrow(member$weights),
      declared_normalization = member$normalization,
      specification = member$specification
    )
  }, members, names(members), alpha)
  structure(list(
    kind = "frame_family",
    normalization = "conservative",
    members = records,
    alpha = alpha
  ), class = "effect_frame_family_spec")
}

# The per-block certificate of contract claim 3b, read off the weights alone:
# a family conserves block by block exactly when every member's columns sum to
# its own alpha. It costs one column sum per member and needs no geometry.
.frame_family_conservation <- function(frame, tolerance) {
  specification <- frame$specification
  index <- frame$index
  if (!identical(specification$kind, "frame_family") ||
      !is.data.frame(index) || !"family" %in% names(index)) {
    return(NULL)
  }
  records <- specification$members
  rows <- split(seq_len(nrow(frame$weights)),
    factor(index$family, levels = names(records)))
  blocks <- lapply(names(records), function(family) {
    positions <- rows[[family]]
    weight <- records[[family]]$alpha
    mass <- as.numeric(Matrix::colSums(
      frame$weights[positions, , drop = FALSE]
    ))
    deviation <- max(abs(mass - weight))
    data.frame(family = family, alpha = weight,
      measurements = length(positions), max_deviation = deviation,
      conserved = deviation <= tolerance, stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, blocks)
  rownames(result) <- NULL
  result
}

#' Diagnose local-to-global conservation of a compiled frame
#'
#' Under a conservative frame every domain feature carries total weight mass
#' one, so local `total` geometries sum exactly to the global geometry:
#' \deqn{\sum_x G_x^{\mathrm{total}} = G_\Omega^{\mathrm{total}}.}{sum_x G_x^total = G_Omega^total.}
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
#'   covered `component` (`"total"`), the frame `normalization`, its
#'   `declared_normalization` together with `metric_folded` (whether a
#'   diagonal metric has been folded into the weights) and the `composition`
#'   law in force (`"none"` for a declared frame, `"diagonal_metric_fold"` for
#'   a folded one), the maximum
#'   per-feature deviation from the conserving `reference_mass`, and the
#'   per-feature mass vector. `reference_mass` is one for a declared frame and
#'   the folded metric diagonal for a metric-folded one, because that frame's
#'   global comparator is read under the same metric. A [frame_family()]
#'   additionally reports `members`, one row per family member giving its
#'   `alpha`, its measurement count, and the deviation of its own per-feature
#'   mass from that `alpha` -- the block-by-block half of the law, which the
#'   whole family conserving does not imply.
#' @family neural domains and frames
#' @seealso [compile_frame()] and the normalization argument of [voxelwise()],
#'   [searchlights()], [regions()], and [whole_brain()].
#' @examples
#' domain <- abstract_domain(4, id = "conservation-example")
#'
#' # A conservative frame gives every feature total weight mass one, so local
#' # `total` geometries sum exactly to the global one.
#' conservative <- compile_frame(voxelwise(), domain)
#' frame_conservation(conservative)$conserved
#'
#' # Under "local" each measurement is rescaled instead, and the per-feature
#' # mass tells you by how much the accounting is off.
#' report <- frame_conservation(compile_frame(whole_brain(), domain))
#' report$conserved
#' report$max_deviation
#' @export
frame_conservation <- function(x, tolerance = 1e-10) {
  .validate_frame_for_compile(x)
  .check_number(tolerance, "tolerance", nonnegative = TRUE)
  mass <- as.numeric(Matrix::colSums(x$weights))
  # The target mass is one for a declared frame. A frame that has had a
  # diagonal metric folded into it carries weights `w_xv d_v`, and its global
  # comparator is read under the same metric, so the conserving target is the
  # metric diagonal rather than one. Comparing such a frame against one would
  # report a conservative frame as unconserved purely because of the metric.
  reference <- .frame_conservation_reference(x)
  fold <- x$metric_folded
  max_deviation <- max(abs(mass - reference))
  members <- .frame_family_conservation(x, tolerance)
  report <- structure(list(
    conserved = max_deviation <= tolerance,
    component = "total",
    normalization = x$normalization,
    declared_normalization = if (is.null(fold)) {
      x$normalization
    } else {
      fold$declared_normalization
    },
    metric_folded = !is.null(fold),
    # The composition law these weights were built under, so a reader can tell
    # what `reference_mass` is a mass *of*. `"none"` is a declared frame whose
    # columns carry unit mass. A plan's metric composition is a separate
    # statement and lives on `$metric_schedule`: `composition = "whitened"`
    # transforms the effect coordinates, never the frame, so a frame used by a
    # whitened plan correctly reports `"none"` here.
    composition = if (is.null(fold)) "none" else fold$composition,
    max_deviation = max_deviation,
    feature_mass = mass,
    reference_mass = reference,
    tolerance = tolerance
  ), class = "effect_frame_conservation")
  if (!is.null(members)) report$members <- members
  report
}
