# Declarative compiler capabilities -----------------------------------------

#' Describe an additive diagonal spatial frame
#'
#' An additive frame is the only frame representation whose locations collapse
#' to rows of one spatial contraction in the version 0.1 compiler.
#'
#' @param weights A finite, nonnegative measurement-by-feature base or sparse
#'   `Matrix` matrix. Supply this or `members`, never both.
#' @param normalization One of `none`, `local` (row sums equal one), or
#'   `conservative` (column sums equal one).
#' @param domain_id Stable identity of the neural feature domain.
#' @param domain Optional exact `effect_domain` or internal domain reference.
#'   The `members` route requires an exact domain.
#' @param members A nonempty list of integer vectors, one per measurement,
#'   giving the domain feature *positions* that measurement is supported on.
#'   This is the neighborhood route described under *Two routes*.
#' @param measurements Identifiers for the rows of the `members` route, one per
#'   measurement, recorded as `$index$measurement`. Defaults to the measurement
#'   positions.
#' @param construction Named list recording how a `members` neighborhood was
#'   produced --- the rule, the provider, its parameters --- kept with the
#'   support pattern and folded into its identity.
#' @param specification Optional named list recording what generated a
#'   `members` frame, kept as `$specification` the way a [compile_frame()]
#'   result keeps the scope it was compiled from.
#' @return An `effect_frame` with `representation = "additive_diagonal"`,
#'   carrying the `$weights` matrix, its `$normalization`, and the `$domain`
#'   reference the weights are bound to.
#' @section Structure:
#' A declared frame carries the operator and the domain it claims, and nothing
#' about how it was chosen.
#'
#' - `$weights`: the measurement-by-feature operator exactly as supplied. Row
#'   `m` holds the weight each domain feature contributes to measurement `m`,
#'   in domain feature order.
#' - `$normalization`: the normalization asserted about those rows. It is
#'   checked, not applied.
#' - `$domain`: the neural domain reference the columns are bound to, and
#'   `$domain_id` its identity.
#'
#' Unlike a [compile_frame()] result, a frame declared from `weights` has no
#' `$index` and no `$specification`: nothing generated it, so views index its
#' measurements by position. Any other element is internal and may change.
#' @section Two routes:
#' A frame can be declared from an operator or built from neighborhoods, and
#' the difference is who applies the normalization.
#'
#' - **`weights`** is the declaration: the operator is taken exactly as
#'   supplied and `normalization` is *checked* against it, so unnormalized
#'   rows are refused rather than silently rescaled.
#' - **`members`** is the neighborhood route, for a spatial provider that
#'   computes its own supports --- searchlights from an external package,
#'   parcels from an atlas, any rule crossform does not implement. Each element
#'   of `members` lists the domain feature positions one measurement covers;
#'   `normalization` is then *applied* to that membership pattern, exactly as
#'   [compile_frame()] applies it to the neighborhoods it computes itself. What
#'   comes back is a frame of the same shape a compiled one has, carrying
#'   `$index`, the `$specification` you recorded, and the support pattern that
#'   makes a locality-aware plan possible.
#'
#' The neighborhood route is the reason an external searchlight provider does
#' not have to reproduce crossform's normalization law or its support
#' bookkeeping to hand back a frame the rest of the package accepts;
#' `neuroim2_searchlights()` is written this way and is the worked example.
#' @seealso [compile_frame()] with [searchlights()], [regions()], or
#'   [voxelwise()], which build additive frames from a neural domain;
#'   [neuroim2_searchlights()] for the `members` route in use; and
#'   [measurement_frame()], which adapts one into oriented measurements.
#' @family neural domains and frames
#' @examples
#' # Two overlapping local averages over four features, declared directly.
#' # `"local"` asserts that each measurement's weights already sum to one, so
#' # a frame row is a weighted mean rather than a weighted sum.
#' domain <- abstract_domain(4, id = "demo:native:v1")
#' frame <- additive_frame(
#'   matrix(c(
#'     1 / 2, 1 / 2, 0, 0,
#'     0, 1 / 3, 1 / 3, 1 / 3
#'   ), 2, 4, byrow = TRUE),
#'   normalization = "local", domain = domain
#' )
#' dim(frame$weights)
#' rowSums(as.matrix(frame$weights))
#'
#' # The assertion is checked, not applied: unnormalized rows are refused
#' # rather than silently rescaled.
#' unnormalized <- try(
#'   additive_frame(matrix(1, 1, 4), normalization = "local", domain = domain),
#'   silent = TRUE
#' )
#' conditionMessage(attr(unnormalized, "condition"))
#'
#' # The declared width must match the domain it claims.
#' wrong <- try(additive_frame(matrix(1, 1, 3), domain = domain), silent = TRUE)
#' conditionMessage(attr(wrong, "condition"))
#'
#' # The neighborhood route: hand over the supports a provider computed and
#' # the declared normalization is applied to them, not asserted about them.
#' neighborhoods <- additive_frame(
#'   members = list(1:2, 2:4), measurements = c("left", "right"),
#'   normalization = "local", domain = domain,
#'   construction = list(kind = "declared_neighborhoods", provider = "example")
#' )
#' rowSums(as.matrix(neighborhoods$weights))
#'
#' # It comes back shaped like a compiled frame: the measurements are named,
#' # not merely positional.
#' neighborhoods$index$measurement
#' @export
additive_frame <- function(weights, normalization = "none",
                           domain_id = "abstract", domain = NULL,
                           members = NULL, measurements = NULL,
                           construction = list(), specification = NULL) {
  support_index <- NULL
  domain_kind <- NULL
  if (!is.null(members)) {
    # The neighborhood route. Everything below the declaration -- the
    # membership pattern, the normalization law, the support bookkeeping -- is
    # crossform's own, so a provider that computes supports supplies exactly
    # the supports and nothing else.
    if (!missing(weights)) {
      .input_error(paste0(
        "Supply `weights` or `members`, not both: `weights` declares an ",
        "operator already normalized, `members` asks for one to be built ",
        "from neighborhoods."
      ), arg = "members", received = "both `weights` and `members`",
        expected = "exactly one of them")
    }
    if (!inherits(domain, "effect_domain")) {
      .input_error(sprintf(paste0(
        "The `members` route needs the exact `effect_domain` its supports ",
        "index into (see `abstract_domain()`, `volume_domain()`, or ",
        "`neuroim2_volume_domain()`); received %s."
      ), .msg_value(domain)),
        arg = "domain", received = .msg_value(domain),
        expected = "an `effect_domain`")
    }
    if (!.is_string(normalization) ||
        !normalization %in% c("none", "local", "conservative")) {
      .input_error(sprintf(paste0(
        "`normalization` must be one of \"none\", \"local\", or ",
        "\"conservative\"; received %s."
      ), .msg_value(normalization)),
        arg = "normalization", received = .msg_value(normalization),
        expected = "one of \"none\", \"local\", or \"conservative\"")
    }
    if (is.null(measurements)) measurements <- seq_along(members)
    if (!is.null(specification) &&
        (!is.list(specification) || !length(specification) ||
         is.null(names(specification)) || anyNA(names(specification)) ||
         any(!nzchar(names(specification))) ||
         anyDuplicated(names(specification)))) {
      .input_error(paste0(
        "`specification` must be a uniquely named list recording what ",
        "generated these measurements."
      ), arg = "specification", received = .msg_value(specification),
        expected = "a uniquely named list")
    }
    domain_kind <- domain$kind
    support_index <- .support_index_from_members(members, domain, measurements,
      construction)
    weights <- .normalize_frame(.support_index_membership(support_index),
      normalization)
  }
  width <- if (length(dim(weights)) == 2L) ncol(weights) else NA_integer_
  domain <- if (is.null(domain)) {
    .positional_domain_reference(width, domain_id)
  } else {
    reference <- .domain_reference(domain)
    if (!identical(domain_id, "abstract") && !identical(domain_id, reference$id)) {
      .contract_error(
        "`domain` and `domain_id` identify different neural domains."
      )
    }
    if (!identical(as.integer(width), reference$n_features)) {
      .input_error("The frame width must match its exact neural domain.")
    }
    reference
  }
  frame <- structure(
    list(
      representation = "additive_diagonal",
      fixed = TRUE,
      locally_estimated = FALSE,
      weights = weights,
      normalization = normalization,
      domain = domain,
      domain_id = domain$id
    ),
    class = "effect_frame"
  )
  .validate_frame_for_compile(frame)
  if (!is.null(support_index)) {
    # The same fields, in the same order, that `compile_frame()` attaches to a
    # frame it built itself, so a family may stack members from either route.
    frame$index <- data.frame(measurement = measurements,
      stringsAsFactors = FALSE)
    frame$domain_kind <- domain_kind
    if (!is.null(specification)) frame$specification <- specification
    frame$support_index <- support_index
  }
  frame
}

# The normalization law itself: what `additive_frame(members = )` applies and
# what `compile_frame()` applies to the neighborhoods it computes. It lives
# beside the frame value rather than beside the frame *specifications*,
# because it is a statement about a weights matrix and about nothing else.
.normalize_frame <- function(weights, normalization) {
  row_mass <- Matrix::rowSums(weights)
  if (any(!is.finite(row_mass)) || any(row_mass <= 0)) {
    .input_error("Every frame measurement must contain at least one feature.")
  }
  if (normalization == "local") {
    weights <- Matrix::Diagonal(x = 1 / row_mass) %*% weights
  } else if (normalization == "conservative") {
    coverage <- Matrix::colSums(weights)
    if (any(!is.finite(coverage)) || any(coverage <= 0)) {
      .input_error("Conservative frames must cover every domain feature.")
    }
    weights <- weights %*% Matrix::Diagonal(x = 1 / coverage)
  }
  # A pattern (n*) matrix has no numeric slot; coerce to numeric so unnormalized
  # membership frames carry explicit unit weights.
  weights <- methods::as(weights, "dMatrix")
  methods::as(methods::as(weights, "generalMatrix"), "CsparseMatrix")
}

#' Describe a fixed or locally estimated factor frame
#'
#' Factor frames use a separate contraction lowering. A locally estimated
#' factor retains location-dependent fitting and is never treated as an
#' additive-frame collapse.
#'
#' @param factors A nonempty list of numeric factor matrices.
#' @param locally_estimated Whether the factors are estimated independently at
#'   each location.
#' @param domain_id Stable identity of the neural feature domain.
#' @param domain Optional exact `effect_domain` or internal domain reference.
#' @return A declarative frame value.
#' @keywords internal
factor_frame <- function(factors, locally_estimated = FALSE,
                         domain_id = "abstract", domain = NULL) {
  if (!is.list(factors) || length(factors) < 1L ||
      !all(vapply(factors, function(x) is.matrix(x) && is.numeric(x), logical(1)))) {
    .input_error("`factors` must be a nonempty list of numeric matrices.")
  }
  .check_flag(locally_estimated, "locally_estimated")
  domain <- if (is.null(domain)) {
    .positional_domain_reference(ncol(factors[[1L]]), domain_id)
  } else {
    reference <- .domain_reference(domain)
    if (!identical(domain_id, "abstract") && !identical(domain_id, reference$id)) {
      .contract_error(
        "`domain` and `domain_id` identify different neural domains."
      )
    }
    reference
  }

  structure(
    list(
      representation = "factor",
      fixed = !locally_estimated,
      locally_estimated = locally_estimated,
      factors = factors,
      domain = domain,
      domain_id = domain$id
    ),
    class = "effect_frame"
  )
}

#' Describe an axis-bound pair query
#'
#' A pair query is the universal rectangular readout. Its rows are bound to one
#' ordered left effect space and its columns to one ordered right effect space;
#' equal dimensions or labels do not substitute for identity.
#'
#' @param H A finite nonempty left-by-right numeric base or `Matrix` matrix.
#' @param left_space,right_space The ordered `effect_space()` identities bound
#'   to the rows and columns of `H`. Unique character coordinates are accepted
#'   as shorthand for an unspecified-basis space.
#' @param metadata Optional compact semantic metadata, used by higher-level
#'   pair-design constructors for balance and design diagnostics.
#' @return An `effect_pair_query` carrying the `$operator`, its bound
#'   `$left_space` and `$right_space`, and any `$metadata` a higher-level
#'   constructor attached.
#' @seealso [pair_lm_query()] and [coupling_contrast()], which compile
#'   designed pair operators; [bilinear_query()] for the symmetric
#'   single-space case; [evaluate_geometry()] for reading a rectangular plan
#'   with one.
#' @family coupling and connectivity views
#' @examples
#' # An encoding-retrieval readout: 2 encoding conditions by 3 retrieval
#' # conditions, so the query is rectangular and cannot be a square RDM.
#' encoding <- effect_space(c("encode_a", "encode_b"), basis_id = "demo:enc")
#' retrieval <- effect_space(
#'   c("retrieve_a", "retrieve_b", "lure"), basis_id = "demo:ret"
#' )
#' query <- pair_query(
#'   matrix(c(1, -0.5, 0, 0, 0.5, -1), 2, 3, byrow = TRUE),
#'   encoding, retrieval
#' )
#' dim(query$operator)
#' query$left_space$coordinates
#'
#' # Axis identity is checked, not inferred: equal dimensions do not make two
#' # different experimental spaces interchangeable.
#' swapped <- try(
#'   pair_query(matrix(1, 3, 2), encoding, retrieval), silent = TRUE
#' )
#' conditionMessage(attr(swapped, "condition"))
#' @export
pair_query <- function(H, left_space, right_space, metadata = list()) {
  matrix_like <- (is.matrix(H) && is.numeric(H)) || inherits(H, "Matrix")
  if (!matrix_like || any(dim(H) < 1L) || any(!is.finite(H))) {
    .input_error("`H` must be a finite, nonempty numeric matrix.")
  }
  left_space <- .as_effect_space(left_space, nrow(H))
  right_space <- .as_effect_space(right_space, ncol(H))
  if ((!is.null(rownames(H)) &&
      !identical(rownames(H), left_space$coordinates)) ||
      (!is.null(colnames(H)) &&
      !identical(colnames(H), right_space$coordinates))) {
    .input_error(
      "Named `H` axes must exactly match their ordered effect spaces."
    )
  }
  if (!is.list(metadata)) {
    .input_error("Pair-query `metadata` must be a list.")
  }
  structure(
    list(
      kind = "pair",
      fixed = TRUE,
      operator = H,
      left_space = left_space,
      right_space = right_space,
      metadata = metadata
    ),
    class = c("effect_pair_query", "effect_query")
  )
}

#' Describe a bilinear geometry query
#'
#' This compatibility constructor describes a symmetric query on one effect
#' space. Use `pair_query()` for distinct or unequal axes.
#'
#' @param operator A finite square symmetric numeric matrix.
#' @param fixed Whether the query is fixed before local data are inspected.
#' @param effects Optional `effect_space()` binding for the operator axes.
#' @return An `effect_query` with `kind = "bilinear"`, carrying the symmetric
#'   `$operator`, the `$fixed` flag, and an optional `$effect_space` binding.
#' @seealso [evaluate_geometry()], which reads a plan with this query, and
#'   [pair_query()] for distinct or unequal axes.
#' @family geometry plans and views
#' @examples
#' # A contrast read as a rank-one bilinear query: `t(c) G c` for c = a - b.
#' contrast <- c(1, -1)
#' query <- bilinear_query(tcrossprod(contrast))
#' query$operator
#'
#' # Evaluate it against a plan without materializing complete geometry.
#' domain <- abstract_domain(3, id = "bilinear-example")
#' relation <- relation(
#'   list(run1 = rbind(a = c(1, 0, 2), b = c(0, 1, 1)),
#'        run2 = rbind(a = c(1.1, 0.1, 1.9), b = c(0.1, 0.9, 1.2))),
#'   domain = domain
#' )
#' plan <- plan_geometry(
#'   relation, compile_frame(whole_brain(), domain),
#'   cross_partitions(relation, independence = "independent")
#' )
#' evaluate_geometry(plan, query = query)$values
#'
#' # A bilinear query is symmetric by construction; asymmetry is refused.
#' asymmetric <- try(bilinear_query(matrix(c(1, 2, 3, 4), 2)), silent = TRUE)
#' conditionMessage(attr(asymmetric, "condition"))
#' @export
bilinear_query <- function(operator, fixed = TRUE, effects = NULL) {
  if (!.is_finite_matrix(operator) || nrow(operator) != ncol(operator) ||
      nrow(operator) < 1L) {
    .input_error("`operator` must be a finite, nonempty square numeric matrix.")
  }
  if (max(abs(operator - t(operator))) > 1e-12) {
    .input_error("`operator` must be symmetric.")
  }
  .check_flag(fixed, "fixed")

  if (!is.null(effects)) effects <- .as_effect_space(effects, nrow(operator))
  structure(
    list(kind = "bilinear", fixed = fixed, operator = operator,
      effect_space = effects),
    class = "effect_query"
  )
}

#' Describe a nonlinear geometry readout
#'
#' @param fun A function applied after local geometry is constructed.
#' @return A declarative query value.
#' @keywords internal
nonlinear_query <- function(fun) {
  if (!is.function(fun)) {
    .input_error("`fun` must be a function.")
  }
  structure(
    list(kind = "nonlinear", fixed = TRUE, fun = fun),
    class = "effect_query"
  )
}

#' Compile a frame-query pair to its algebraic lowering
#'
#' Only an additive diagonal frame paired with a fixed bilinear query admits
#' the searchlight-collapse lowering. Factor frames, locally estimated
#' transforms, adaptive queries, and nonlinear readouts remain distinct work.
#'
#' @param frame An `effect_frame`.
#' @param query An `effect_query`.
#' @return A small compiler-decision value.
#' @keywords internal
compile_lowering <- function(frame, query) {
  .check_class(frame, "effect_frame", "frame", what = "an effect frame")
  .check_class(query, "effect_query", "query", what = "an effect query")
  .validate_frame_for_compile(frame)
  .validate_query_for_compile(query)

  if (isTRUE(frame$locally_estimated)) {
    return(.new_lowering("location_dependent_fit", FALSE,
      "the spatial transform is estimated separately at each location"))
  }
  if (!isTRUE(query$fixed)) {
    return(.new_lowering("adaptive_query", FALSE,
      "the query depends on local data"))
  }
  if (identical(query$kind, "pair")) {
    return(.new_lowering("two_sided_pair_form", FALSE,
      "ordered two-sided relation tasks are not compiled by this lowering yet"))
  }
  if (!identical(query$kind, "bilinear")) {
    return(.new_lowering("nonlinear_readout", FALSE,
      "the readout is not bilinear in the local geometry"))
  }
  if (!identical(frame$representation, "additive_diagonal")) {
    return(.new_lowering("factor_contraction", FALSE,
      "a generic factor frame is not an additive diagonal frame"))
  }

  .new_lowering("additive_contraction", TRUE,
    "fixed bilinear readout over an additive diagonal frame")
}

.frame_values <- function(weights) {
  if (is.matrix(weights) && is.numeric(weights)) return(as.numeric(weights))
  if (inherits(weights, "Matrix")) {
    if (inherits(weights, "sparseMatrix") && "x" %in% methods::slotNames(weights)) {
      return(methods::slot(weights, "x"))
    }
    return(as.numeric(weights))
  }
  .input_error(
    "Additive weights must be a numeric base matrix or Matrix object."
  )
}

.validate_domain_id <- function(domain_id) {
  if (!.is_string(domain_id)) {
    .input_error("Frame `domain_id` must be one nonempty identifier.")
  }
}

.validate_frame_for_compile <- function(frame) {
  if (!inherits(frame, "effect_frame")) {
    .input_error(sprintf(paste0(
      "Expected a compiled `effect_frame` from `compile_frame()` (or ",
      "`neuroim2_searchlights()`); received %s."
    ), .msg_value(frame)))
  }
  if (!is.list(frame)) {
    .input_error("Frame fields are missing or noncanonical.")
  }
  expected_names <- if (identical(frame$representation, "additive_diagonal")) {
    c("representation", "fixed", "locally_estimated", "weights",
      "normalization", "domain", "domain_id")
  } else if (identical(frame$representation, "factor")) {
    c("representation", "fixed", "locally_estimated", "factors", "domain",
      "domain_id")
  } else {
    .input_error("Unknown frame representation.")
  }
  # Optional fields keep one canonical order, but they are not all produced
  # together: `compile_frame()` attaches the leading three (plus
  # `support_index` for neighborhood scopes), while `.metric_additive_frame()`
  # attaches `metric_folded` to a declared frame that has none of them. The
  # admitted trailing block is therefore an ordered *subsequence* of
  # `optional_names`, not a prefix of it. Unknown names, duplicates, and
  # reordering are still refused.
  optional_names <- c("index", "domain_kind", "specification", "support_index",
    "metric_folded")
  trailing_names <- names(frame)[length(expected_names) +
    seq_len(max(0L, length(frame) - length(expected_names)))]
  positions <- match(trailing_names, optional_names)
  valid_trailing <- !anyNA(positions) &&
    (length(positions) < 2L || all(diff(positions) > 0L))
  if (!identical(names(frame)[seq_along(expected_names)], expected_names) ||
      !valid_trailing || anyDuplicated(names(frame))) {
    .input_error("Frame fields are missing or noncanonical.")
  }
  .validate_domain_id(frame$domain_id)
  domain <- .validate_domain_reference(frame$domain)
  if (!identical(frame$domain_id, domain$id)) {
    .contract_error(
      "Frame domain label is inconsistent with its exact domain reference."
    )
  }
  if (identical(frame$representation, "additive_diagonal")) {
    if (!isTRUE(frame$fixed) || isTRUE(frame$locally_estimated)) {
      .input_error(
        "An additive collapse frame must be fixed and not locally estimated."
      )
    }
    weights <- frame$weights
    values <- .frame_values(weights)
    if (length(dim(weights)) != 2L || any(dim(weights) < 1L) ||
        any(!is.finite(values)) || any(values < 0)) {
      .input_error(paste0(
        "Additive weights must have positive dimensions and finite ",
        "nonnegative values."
      ))
    }
    if (!identical(as.integer(ncol(weights)), domain$n_features)) {
      .contract_error(
        "Frame width is inconsistent with its exact neural domain."
      )
    }
    if (!is.null(frame$support_index)) {
      support_index <- .validate_support_index(frame$support_index)
      if (!.same_domain_reference(support_index$domain, domain) ||
          length(support_index$node_ids) != nrow(weights) ||
          (!is.null(frame$index) &&
           !identical(support_index$node_ids, frame$index$measurement))) {
        .contract_error(
          "Frame support topology is inconsistent with its spatial frame."
        )
      }
    }
    normalization <- frame$normalization
    if (!is.character(normalization) || length(normalization) != 1L ||
        !normalization %in% c("none", "local", "conservative")) {
      .input_error("Frame normalization must be none, local, or conservative.")
    }
    # A frame index is optional, but when present it is the row labelling the
    # rest of the package reads: the executor takes a result's measurement
    # identifiers from it, and a frame family's per-row metadata is carried in
    # it. An index that does not line up with the weights, or that names one
    # row twice, would mislabel every value downstream rather than fail, so it
    # is refused here.
    if (!is.null(frame$index)) {
      measurement <- if (is.data.frame(frame$index)) {
        frame$index$measurement
      } else {
        NULL
      }
      if (is.null(measurement) || nrow(frame$index) != nrow(weights) ||
          anyNA(measurement) || anyDuplicated(measurement)) {
        .input_error(
          "A frame index must name every measurement row exactly once."
        )
      }
    }
    row_mass <- if (inherits(weights, "Matrix")) Matrix::rowSums(weights) else rowSums(weights)
    if (any(!is.finite(row_mass)) || any(row_mass <= 0)) {
      .input_error("Every additive frame row must have finite positive mass.")
    }
    tolerance <- 1e-12
    if (normalization == "local" && any(abs(row_mass - 1) > tolerance)) {
      .input_error("Locally normalized frame rows must sum to one.")
    }
    if (normalization == "conservative") {
      column_mass <- if (inherits(weights, "Matrix"))
        Matrix::colSums(weights) else colSums(weights)
      if (any(!is.finite(column_mass)) || any(abs(column_mass - 1) > tolerance)) {
        .input_error("Conservative frame columns must sum to one.")
      }
    }
  } else if (identical(frame$representation, "factor")) {
    if (!is.list(frame$factors) || length(frame$factors) < 1L ||
        !all(vapply(frame$factors, function(x) {
          is.matrix(x) && is.numeric(x) && all(dim(x) > 0) && all(is.finite(x))
        }, logical(1)))) {
      .input_error("Factor frames require finite nonempty numeric matrices.")
    }
    widths <- vapply(frame$factors, ncol, integer(1))
    if (length(unique(widths)) != 1L) {
      .input_error(
        "All factor-frame elements must share one feature dimension."
      )
    }
    if (!identical(widths[[1L]], domain$n_features)) {
      .contract_error(
        "Factor-frame width is inconsistent with its exact neural domain."
      )
    }
  }
  if (!is.null(frame$metric_folded)) .validate_frame_metric_fold(frame)
  invisible(frame)
}

# A metric-folded frame is a composed object: its weights are a declared
# frame's weights times a diagonal metric. Folding rescales every column by
# that diagonal, so folded weights cannot satisfy the "conservative"
# column-mass law (nor, generally, the "local" row-mass law) that
# `.validate_frame_for_compile()` checks. The composed frame therefore
# declares `normalization = "none"`, which is true of the weights it actually
# carries, and keeps the frame's declared normalization here as provenance
# rather than asserting it of the composed weights. `reference_mass` is the
# per-feature mass a conserving fold must reach -- the metric diagonal -- so
# `frame_conservation()` can certify the fold instead of comparing against a
# unit mass the folded weights were never supposed to have.
.validate_frame_metric_fold <- function(frame) {
  fold <- frame$metric_folded
  expected <- c("folded", "declared_normalization", "metric_kind",
    "metric_signature", "schedule_kind", "composition", "reference_mass")
  if (!identical(frame$representation, "additive_diagonal")) {
    .input_error("Only an additive frame can carry a folded metric.")
  }
  if (!is.list(fold) || !identical(names(fold), expected) ||
      !isTRUE(fold$folded) ||
      !.is_string(fold$declared_normalization) ||
      !fold$declared_normalization %in% c("none", "local", "conservative") ||
      !.is_string(fold$metric_kind) || !.is_string(fold$metric_signature) ||
      !.is_string(fold$schedule_kind) || !.is_string(fold$composition) ||
      !identical(frame$normalization, "none")) {
    .input_error("Frame metric-fold provenance is missing or noncanonical.")
  }
  reference <- fold$reference_mass
  if (!is.numeric(reference) || length(reference) != ncol(frame$weights) ||
      any(!is.finite(reference)) || any(reference < 0)) {
    .input_error("Frame metric-fold provenance is missing or noncanonical.")
  }
  invisible(frame)
}

# The per-feature mass a frame's columns must carry for local `total`
# geometries to sum to the whole-support geometry. It is one for a declared
# frame and the folded metric diagonal for a metric-folded frame, because the
# global comparator is then read under the same metric.
.frame_conservation_reference <- function(frame) {
  fold <- frame$metric_folded
  if (is.null(fold)) rep(1, ncol(frame$weights)) else fold$reference_mass
}

.validate_query_for_compile <- function(query) {
  if (!inherits(query, "effect_query") || !is.list(query) ||
      !.is_string(query$kind, allow_empty = TRUE)) {
    .input_error("Query fields are missing or noncanonical.")
  }
  if (identical(query$kind, "pair")) {
    expected <- c("kind", "fixed", "operator", "left_space", "right_space",
      "metadata")
    if (!inherits(query, "effect_pair_query") ||
        !identical(names(query), expected) || !identical(query$fixed, TRUE)) {
      .input_error("Pair-query fields are missing or noncanonical.")
    }
    operator <- query$operator
    matrix_like <- (is.matrix(operator) && is.numeric(operator)) ||
      inherits(operator, "Matrix")
    if (!matrix_like || any(dim(operator) < 1L) || any(!is.finite(operator))) {
      .input_error("Pair-query operators must be finite nonempty matrices.")
    }
    left_space <- .validate_effect_space(query$left_space)
    right_space <- .validate_effect_space(query$right_space)
    if (length(left_space$coordinates) != nrow(operator) ||
        length(right_space$coordinates) != ncol(operator)) {
      .contract_error(
        "Pair-query dimensions must match their bound effect spaces."
      )
    }
    if (!is.list(query$metadata)) {
      .input_error("Pair-query metadata must be a list.")
    }
    rebuilt <- pair_query(operator, left_space, right_space, query$metadata)
    if (!identical(query, rebuilt)) {
      .input_error("Pair-query fields are missing or noncanonical.")
    }
    return(invisible(rebuilt))
  }
  if (!identical(query$kind, "bilinear")) {
    if (!identical(query$kind, "nonlinear") ||
        !identical(names(query), c("kind", "fixed", "fun")) ||
        !identical(query$fixed, TRUE) || !is.function(query$fun)) {
      .input_error("Query fields are missing or noncanonical.")
    }
    return(invisible(query))
  }
  if (!identical(names(query), c("kind", "fixed", "operator", "effect_space"))) {
    .input_error("Query fields are missing or noncanonical.")
  }
  operator <- query$operator
  if (!.is_finite_matrix(operator) || nrow(operator) < 1L ||
      nrow(operator) != ncol(operator) ||
      max(abs(operator - t(operator))) > 1e-12) {
    .input_error(
      "Bilinear query operators must be finite, square, and symmetric."
    )
  }
  if (!.is_flag(query$fixed)) {
    .input_error("Query fixedness must be one logical value.")
  }
  if (!is.null(query$effect_space)) {
    effect_space <- .validate_effect_space(query$effect_space)
    if (length(query$effect_space$coordinates) != nrow(operator)) {
      .contract_error("Query effect-space dimension must match its operator.")
    }
  } else {
    effect_space <- NULL
  }
  rebuilt <- bilinear_query(operator, query$fixed, effect_space)
  if (!identical(query, rebuilt)) {
    .input_error("Query fields are missing or noncanonical.")
  }
  invisible(rebuilt)
}

.new_lowering <- function(kind, collapsed, reason) {
  structure(
    list(kind = kind, collapsed = collapsed, reason = reason),
    class = "effect_lowering"
  )
}
