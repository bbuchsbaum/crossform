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
#' @param radius Positive radius in domain coordinate units.
#' @param normalization Explicit frame normalization.
#' @return An `effect_frame_spec` with `$kind` `"searchlights"`, the requested
#'   `$radius`, and `$normalization`. Pass it to [compile_frame()].
#' @family neural domains and frames
#' @seealso [compile_frame()], [frame_conservation()] to check the
#'   normalization you chose, and [neuroim2_searchlights()] for neighborhoods
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
#' @export
searchlights <- function(radius, normalization = "local") {
  if (missing(radius)) {
    stop(paste0(
      "`radius` is required: pass one positive radius in the domain's own ",
      "coordinate units, for example `searchlights(8)` for 8 mm on a volume ",
      "domain."
    ), call. = FALSE)
  }
  if (!is.numeric(radius) || length(radius) != 1L || is.na(radius) ||
      !is.finite(radius) || radius <= 0) {
    stop(sprintf(paste0(
      "`radius` must be one positive finite number in the domain's ",
      "coordinate units; received %s."
    ), if (is.numeric(radius) && length(radius) == 1L) {
      paste0("`", format(radius), "`")
    } else {
      .msg_value(radius)
    }), call. = FALSE)
  }
  .frame_spec("searchlights", normalization, radius = radius)
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
    stop(paste0(
      "`labels` is required: pass one region label per neural feature, in ",
      "domain feature order, for example `regions(atlas_labels)`."
    ), call. = FALSE)
  }
  if (length(labels) < 1L || !(is.atomic(labels) || is.factor(labels))) {
    stop(sprintf(paste0(
      "`labels` must be a nonempty atomic or factor vector with one region ",
      "label per neural feature; received %s."
    ), .msg_value(labels)), call. = FALSE)
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
    stop(sprintf(paste0(
      "`domain` must be an `effect_domain` (see `abstract_domain()`, ",
      "`volume_domain()`, or `neuroim2_volume_domain()`); received %s."
    ), if (missing(domain)) "no argument" else .msg_value(domain)),
      call. = FALSE)
  }
  .validate_domain(domain)
  if (missing(specification) || !inherits(specification, "effect_frame_spec")) {
    stop(sprintf(paste0(
      "`specification` must be a frame specification from `voxelwise()`, ",
      "`searchlights()`, `regions()`, or `whole_brain()`; received %s."
    ), if (missing(specification)) "no argument" else
      .msg_value(specification)), call. = FALSE)
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
      stop(sprintf(paste0(
        "`regions()` supplied %s but the domain has %s. Labels are read in ",
        "domain feature order, one per feature."
      ), .msg_count(length(labels), "label"),
        .msg_count(n, "feature")), call. = FALSE)
    }
    present <- !is.na(labels) & nzchar(as.character(labels))
    if (!any(present)) {
      stop(paste0(
        "Every region label is missing or empty, so the frame would have no ",
        "measurements. At least one feature must carry a label."
      ), call. = FALSE)
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
    voxels = voxelwise(specification$normalization),
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
  # A pattern (n*) matrix has no numeric slot; coerce to numeric so unnormalized
  # membership frames carry explicit unit weights.
  weights <- methods::as(weights, "dMatrix")
  methods::as(methods::as(weights, "generalMatrix"), "CsparseMatrix")
}
