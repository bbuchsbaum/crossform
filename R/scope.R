# Declarative compiler capabilities -----------------------------------------

#' Describe an additive diagonal spatial frame
#'
#' An additive frame is the only frame representation whose locations collapse
#' to rows of one spatial contraction in the version 0.1 compiler.
#'
#' @param weights A finite, nonnegative measurement-by-feature base or sparse
#'   `Matrix` matrix.
#' @param normalization One of `none`, `local` (row sums equal one), or
#'   `conservative` (column sums equal one).
#' @param domain_id Stable identity of the neural feature domain.
#' @param domain Optional exact `effect_domain` or internal domain reference.
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
#' Unlike a [compile_frame()] result, a declared frame has no `$index` and no
#' `$specification`: nothing generated it, so views index its measurements by
#' position. Any other element is internal and may change.
#' @seealso [compile_frame()] with [searchlights()], [regions()], or
#'   [voxelwise()], which build additive frames from a neural domain;
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
#' @export
additive_frame <- function(weights, normalization = "none",
                           domain_id = "abstract", domain = NULL) {
  width <- if (length(dim(weights)) == 2L) ncol(weights) else NA_integer_
  domain <- if (is.null(domain)) {
    .positional_domain_reference(width, domain_id)
  } else {
    reference <- .domain_reference(domain)
    if (!identical(domain_id, "abstract") && !identical(domain_id, reference$id)) {
      stop("`domain` and `domain_id` identify different neural domains.",
        call. = FALSE)
    }
    if (!identical(as.integer(width), reference$n_features)) {
      stop("The frame width must match its exact neural domain.", call. = FALSE)
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
  frame
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
    stop("`factors` must be a nonempty list of numeric matrices.", call. = FALSE)
  }
  if (!is.logical(locally_estimated) || length(locally_estimated) != 1L ||
      is.na(locally_estimated)) {
    stop("`locally_estimated` must be TRUE or FALSE.", call. = FALSE)
  }
  domain <- if (is.null(domain)) {
    .positional_domain_reference(ncol(factors[[1L]]), domain_id)
  } else {
    reference <- .domain_reference(domain)
    if (!identical(domain_id, "abstract") && !identical(domain_id, reference$id)) {
      stop("`domain` and `domain_id` identify different neural domains.",
        call. = FALSE)
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
    stop("`H` must be a finite, nonempty numeric matrix.", call. = FALSE)
  }
  left_space <- .as_effect_space(left_space, nrow(H))
  right_space <- .as_effect_space(right_space, ncol(H))
  if ((!is.null(rownames(H)) &&
      !identical(rownames(H), left_space$coordinates)) ||
      (!is.null(colnames(H)) &&
      !identical(colnames(H), right_space$coordinates))) {
    stop("Named `H` axes must exactly match their ordered effect spaces.",
      call. = FALSE)
  }
  if (!is.list(metadata)) {
    stop("Pair-query `metadata` must be a list.", call. = FALSE)
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
  if (!is.matrix(operator) || !is.numeric(operator) ||
      nrow(operator) != ncol(operator) || nrow(operator) < 1L ||
      any(!is.finite(operator))) {
    stop("`operator` must be a finite, nonempty square numeric matrix.", call. = FALSE)
  }
  if (max(abs(operator - t(operator))) > 1e-12) {
    stop("`operator` must be symmetric.", call. = FALSE)
  }
  if (!is.logical(fixed) || length(fixed) != 1L || is.na(fixed)) {
    stop("`fixed` must be TRUE or FALSE.", call. = FALSE)
  }

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
    stop("`fun` must be a function.", call. = FALSE)
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
  if (!inherits(frame, "effect_frame")) {
    stop("`frame` must be an effect frame.", call. = FALSE)
  }
  if (!inherits(query, "effect_query")) {
    stop("`query` must be an effect query.", call. = FALSE)
  }
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
  stop("Additive weights must be a numeric base matrix or Matrix object.",
    call. = FALSE)
}

.validate_domain_id <- function(domain_id) {
  if (!is.character(domain_id) || length(domain_id) != 1L ||
      is.na(domain_id) || !nzchar(domain_id)) {
    stop("Frame `domain_id` must be one nonempty identifier.", call. = FALSE)
  }
}

.validate_frame_for_compile <- function(frame) {
  if (!inherits(frame, "effect_frame")) {
    stop(sprintf(paste0(
      "Expected a compiled `effect_frame` from `compile_frame()` (or ",
      "`neuroim2_searchlights()`); received %s."
    ), .msg_value(frame)), call. = FALSE)
  }
  if (!is.list(frame)) {
    stop("Frame fields are missing or noncanonical.", call. = FALSE)
  }
  expected_names <- if (identical(frame$representation, "additive_diagonal")) {
    c("representation", "fixed", "locally_estimated", "weights",
      "normalization", "domain", "domain_id")
  } else if (identical(frame$representation, "factor")) {
    c("representation", "fixed", "locally_estimated", "factors", "domain",
      "domain_id")
  } else {
    stop("Unknown frame representation.", call. = FALSE)
  }
  optional_names <- c("index", "domain_kind", "specification", "support_index")
  trailing_names <- names(frame)[length(expected_names) +
    seq_len(max(0L, length(frame) - length(expected_names)))]
  valid_trailing <- identical(trailing_names, character()) ||
    identical(trailing_names, utils::head(optional_names,
      length(trailing_names)))
  if (!identical(names(frame)[seq_along(expected_names)], expected_names) ||
      !valid_trailing || anyDuplicated(names(frame))) {
    stop("Frame fields are missing or noncanonical.", call. = FALSE)
  }
  .validate_domain_id(frame$domain_id)
  domain <- .validate_domain_reference(frame$domain)
  if (!identical(frame$domain_id, domain$id)) {
    stop("Frame domain label is inconsistent with its exact domain reference.",
      call. = FALSE)
  }
  if (identical(frame$representation, "additive_diagonal")) {
    if (!isTRUE(frame$fixed) || isTRUE(frame$locally_estimated)) {
      stop("An additive collapse frame must be fixed and not locally estimated.",
        call. = FALSE)
    }
    weights <- frame$weights
    values <- .frame_values(weights)
    if (length(dim(weights)) != 2L || any(dim(weights) < 1L) ||
        any(!is.finite(values)) || any(values < 0)) {
      stop("Additive weights must have positive dimensions and finite nonnegative values.",
        call. = FALSE)
    }
    if (!identical(as.integer(ncol(weights)), domain$n_features)) {
      stop("Frame width is inconsistent with its exact neural domain.",
        call. = FALSE)
    }
    if (!is.null(frame$support_index)) {
      support_index <- .validate_support_index(frame$support_index)
      if (!.same_domain_reference(support_index$domain, domain) ||
          length(support_index$node_ids) != nrow(weights) ||
          (!is.null(frame$index) &&
           !identical(support_index$node_ids, frame$index$measurement))) {
        stop("Frame support topology is inconsistent with its spatial frame.",
          call. = FALSE)
      }
    }
    normalization <- frame$normalization
    if (!is.character(normalization) || length(normalization) != 1L ||
        !normalization %in% c("none", "local", "conservative")) {
      stop("Frame normalization must be none, local, or conservative.",
        call. = FALSE)
    }
    row_mass <- if (inherits(weights, "Matrix")) Matrix::rowSums(weights) else rowSums(weights)
    if (any(!is.finite(row_mass)) || any(row_mass <= 0)) {
      stop("Every additive frame row must have finite positive mass.", call. = FALSE)
    }
    tolerance <- 1e-12
    if (normalization == "local" && any(abs(row_mass - 1) > tolerance)) {
      stop("Locally normalized frame rows must sum to one.", call. = FALSE)
    }
    if (normalization == "conservative") {
      column_mass <- if (inherits(weights, "Matrix"))
        Matrix::colSums(weights) else colSums(weights)
      if (any(!is.finite(column_mass)) || any(abs(column_mass - 1) > tolerance)) {
        stop("Conservative frame columns must sum to one.", call. = FALSE)
      }
    }
  } else if (identical(frame$representation, "factor")) {
    if (!is.list(frame$factors) || length(frame$factors) < 1L ||
        !all(vapply(frame$factors, function(x) {
          is.matrix(x) && is.numeric(x) && all(dim(x) > 0) && all(is.finite(x))
        }, logical(1)))) {
      stop("Factor frames require finite nonempty numeric matrices.", call. = FALSE)
    }
    widths <- vapply(frame$factors, ncol, integer(1))
    if (length(unique(widths)) != 1L) {
      stop("All factor-frame elements must share one feature dimension.",
        call. = FALSE)
    }
    if (!identical(widths[[1L]], domain$n_features)) {
      stop("Factor-frame width is inconsistent with its exact neural domain.",
        call. = FALSE)
    }
  }
  invisible(frame)
}

.validate_query_for_compile <- function(query) {
  if (!inherits(query, "effect_query") || !is.list(query) ||
      !is.character(query$kind) || length(query$kind) != 1L || is.na(query$kind)) {
    stop("Query fields are missing or noncanonical.", call. = FALSE)
  }
  if (identical(query$kind, "pair")) {
    expected <- c("kind", "fixed", "operator", "left_space", "right_space",
      "metadata")
    if (!inherits(query, "effect_pair_query") ||
        !identical(names(query), expected) || !identical(query$fixed, TRUE)) {
      stop("Pair-query fields are missing or noncanonical.", call. = FALSE)
    }
    operator <- query$operator
    matrix_like <- (is.matrix(operator) && is.numeric(operator)) ||
      inherits(operator, "Matrix")
    if (!matrix_like || any(dim(operator) < 1L) || any(!is.finite(operator))) {
      stop("Pair-query operators must be finite nonempty matrices.", call. = FALSE)
    }
    left_space <- .validate_effect_space(query$left_space)
    right_space <- .validate_effect_space(query$right_space)
    if (length(left_space$coordinates) != nrow(operator) ||
        length(right_space$coordinates) != ncol(operator)) {
      stop("Pair-query dimensions must match their bound effect spaces.",
        call. = FALSE)
    }
    if (!is.list(query$metadata)) {
      stop("Pair-query metadata must be a list.", call. = FALSE)
    }
    rebuilt <- pair_query(operator, left_space, right_space, query$metadata)
    if (!identical(query, rebuilt)) {
      stop("Pair-query fields are missing or noncanonical.", call. = FALSE)
    }
    return(invisible(rebuilt))
  }
  if (!identical(query$kind, "bilinear")) {
    if (!identical(query$kind, "nonlinear") ||
        !identical(names(query), c("kind", "fixed", "fun")) ||
        !identical(query$fixed, TRUE) || !is.function(query$fun)) {
      stop("Query fields are missing or noncanonical.", call. = FALSE)
    }
    return(invisible(query))
  }
  if (!identical(names(query), c("kind", "fixed", "operator", "effect_space"))) {
    stop("Query fields are missing or noncanonical.", call. = FALSE)
  }
  operator <- query$operator
  if (!is.matrix(operator) || !is.numeric(operator) || nrow(operator) < 1L ||
      nrow(operator) != ncol(operator) || any(!is.finite(operator)) ||
      max(abs(operator - t(operator))) > 1e-12) {
    stop("Bilinear query operators must be finite, square, and symmetric.",
      call. = FALSE)
  }
  if (!is.logical(query$fixed) || length(query$fixed) != 1L || is.na(query$fixed)) {
    stop("Query fixedness must be one logical value.", call. = FALSE)
  }
  if (!is.null(query$effect_space)) {
    effect_space <- .validate_effect_space(query$effect_space)
    if (length(query$effect_space$coordinates) != nrow(operator)) {
      stop("Query effect-space dimension must match its operator.", call. = FALSE)
    }
  } else {
    effect_space <- NULL
  }
  rebuilt <- bilinear_query(operator, query$fixed, effect_space)
  if (!identical(query, rebuilt)) {
    stop("Query fields are missing or noncanonical.", call. = FALSE)
  }
  invisible(rebuilt)
}

.new_lowering <- function(kind, collapsed, reason) {
  structure(
    list(kind = kind, collapsed = collapsed, reason = reason),
    class = "effect_lowering"
  )
}
