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
#' @return A declarative frame value.
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

#' Describe a bilinear geometry query
#'
#' @param operator A finite square numeric matrix.
#' @param fixed Whether the query is fixed before local data are inspected.
#' @param effects Optional `effect_space()` binding for the operator axes.
#' @return A declarative query value.
#' @export
bilinear_query <- function(operator, fixed = TRUE, effects = NULL) {
  if (!is.matrix(operator) || !is.numeric(operator) ||
      nrow(operator) != ncol(operator) || nrow(operator) < 1L ||
      any(!is.finite(operator))) {
    stop("`operator` must be a finite, nonempty square numeric matrix.", call. = FALSE)
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
  if (!inherits(frame, "effect_frame") || !is.list(frame)) {
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
  optional_names <- c("index", "domain_kind", "specification")
  trailing_names <- names(frame)[length(expected_names) +
    seq_len(max(0L, length(frame) - length(expected_names)))]
  if (!identical(names(frame)[seq_along(expected_names)], expected_names) ||
      !(identical(trailing_names, character()) ||
        identical(trailing_names, optional_names)) || anyDuplicated(names(frame))) {
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
