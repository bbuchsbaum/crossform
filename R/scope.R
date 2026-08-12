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
#' @return A declarative frame value.
#' @export
additive_frame <- function(weights, normalization = "none",
                           domain_id = "abstract") {
  frame <- structure(
    list(
      representation = "additive_diagonal",
      fixed = TRUE,
      locally_estimated = FALSE,
      weights = weights,
      normalization = normalization,
      domain_id = domain_id
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
#' @return A declarative frame value.
#' @export
factor_frame <- function(factors, locally_estimated = FALSE,
                         domain_id = "abstract") {
  if (!is.list(factors) || length(factors) < 1L ||
      !all(vapply(factors, function(x) is.matrix(x) && is.numeric(x), logical(1)))) {
    stop("`factors` must be a nonempty list of numeric matrices.", call. = FALSE)
  }
  if (!is.logical(locally_estimated) || length(locally_estimated) != 1L ||
      is.na(locally_estimated)) {
    stop("`locally_estimated` must be TRUE or FALSE.", call. = FALSE)
  }

  structure(
    list(
      representation = "factor",
      fixed = !locally_estimated,
      locally_estimated = locally_estimated,
      factors = factors,
      domain_id = domain_id
    ),
    class = "effect_frame"
  )
}

#' Describe a bilinear geometry query
#'
#' @param operator A finite square numeric matrix.
#' @param fixed Whether the query is fixed before local data are inspected.
#' @return A declarative query value.
#' @export
bilinear_query <- function(operator, fixed = TRUE) {
  if (!is.matrix(operator) || !is.numeric(operator) ||
      nrow(operator) != ncol(operator) || nrow(operator) < 1L ||
      any(!is.finite(operator))) {
    stop("`operator` must be a finite, nonempty square numeric matrix.", call. = FALSE)
  }
  if (!is.logical(fixed) || length(fixed) != 1L || is.na(fixed)) {
    stop("`fixed` must be TRUE or FALSE.", call. = FALSE)
  }

  structure(
    list(kind = "bilinear", fixed = fixed, operator = operator),
    class = "effect_query"
  )
}

#' Describe a nonlinear geometry readout
#'
#' @param fun A function applied after local geometry is constructed.
#' @return A declarative query value.
#' @export
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
#' @export
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
  .validate_domain_id(frame$domain_id)
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
  } else {
    stop("Unknown frame representation.", call. = FALSE)
  }
  invisible(frame)
}

.validate_query_for_compile <- function(query) {
  if (!identical(query$kind, "bilinear")) return(invisible(query))
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
  invisible(query)
}

.new_lowering <- function(kind, collapsed, reason) {
  structure(
    list(kind = kind, collapsed = collapsed, reason = reason),
    class = "effect_lowering"
  )
}
