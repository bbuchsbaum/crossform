# Scientific views of complete geometry -----------------------------------

.align_contrast <- function(value, effects) {
  if (!is.numeric(value) || length(value) != length(effects) ||
      any(!is.finite(value))) {
    stop("A contrast must contain one finite numeric weight per effect.",
      call. = FALSE)
  }
  if (!is.null(names(value))) {
    if (anyNA(names(value)) || any(!nzchar(names(value))) ||
        anyDuplicated(names(value)) || !setequal(names(value), effects)) {
      stop("Named contrast weights must identify every effect exactly once.",
        call. = FALSE)
    }
    value <- value[effects]
  }
  stats::setNames(as.numeric(value), effects)
}

#' Read a contrast as signed, coherent, configuration, and total evidence
#'
#' This view never removes a mean pattern or reruns an analysis. It reads the
#' exact coherent/configuration decomposition already contained in `x`.
#'
#' @param x An `effect_geometry_plan` or complete `effect_geometry`.
#' @param weights One finite contrast weight per named experimental effect.
#' @return An `effect_contrast_view` containing signed marginals and the three
#'   energy components. `coherence_fraction` is reported only where the raw
#'   cross-generalized components form a nonnegative partition.
#' @export
contrast <- function(x, weights) {
  if (inherits(x, "effect_geometry_plan")) {
    .validate_geometry_plan(x)
    weights <- .align_contrast(
      weights, x$task$left_relation$effect_space$coordinates
    )
    return(.run_geometry_compiler(
      x,
      query = bilinear_query(tcrossprod(weights)),
      component = "contrast",
      signed_query = weights
    ))
  }
  if (!inherits(x, "effect_geometry")) {
    stop("`x` must be a geometry plan or complete effect_geometry.",
      call. = FALSE)
  }
  .validate_effect_geometry(x, probe = FALSE)
  weights <- .align_contrast(weights, x$effects)
  query <- bilinear_query(tcrossprod(weights))
  total <- drop(query_geometry(x, query, "total")$values)
  coherent <- drop(query_geometry(x, query, "coherent")$values)
  .new_effect_contrast_view(
    total, coherent, x$marginals, weights, x$index, x$receipt
  )
}

.new_effect_contrast_view <- function(total, coherent, marginals, weights,
                                      index, receipt) {
  total <- drop(total)
  coherent <- drop(coherent)
  configuration <- total - coherent
  signed <- lapply(marginals, function(value) drop(value %*% weights))
  if (identical(attr(marginals, "semantics"), "undirected_endpoint")) {
    signed <- signed$endpoint
  }
  fraction_valid <- is.finite(total) & total > 0 & coherent >= 0 &
    configuration >= 0
  coherence_fraction <- rep(NA_real_, length(total))
  coherence_fraction[fraction_valid] <- coherent[fraction_valid] /
    total[fraction_valid]

  structure(
    list(
      signed = signed,
      coherent = coherent,
      configuration = configuration,
      total = total,
      coherence_fraction = coherence_fraction,
      coherence_fraction_valid = fraction_valid,
      weights = weights,
      index = index,
      receipt = receipt
    ),
    class = "effect_contrast_view"
  )
}

.self_geometry_source <- function(x, operation, complete = FALSE) {
  if (inherits(x, "effect_geometry_plan")) {
    if (isTRUE(complete)) {
      stop(sprintf("%s requires a complete effect form.", operation),
        call. = FALSE)
    }
    .validate_geometry_plan(x)
    space <- x$task$left_relation$effect_space
    return(list(
      kind = "plan",
      effects = space$coordinates,
      codec = "symmetric_packed",
      index = .compiler_index(x$frame)
    ))
  }
  if (!inherits(x, "effect_form")) {
    stop(sprintf(
      "%s requires a geometry plan or complete effect form.", operation
    ), call. = FALSE)
  }
  .validate_effect_form(x, probe = FALSE)
  if (!isTRUE(x$capabilities$self_form) ||
      !isTRUE(x$capabilities$symmetric)) {
    stop(sprintf("%s requires a symmetric self form.", operation), call. = FALSE)
  }
  list(
    kind = "form",
    effects = x$effects,
    codec = x$codec,
    index = x$index
  )
}

.require_self_form_view <- function(x, operation) {
  .self_geometry_source(x, operation, complete = TRUE)
  invisible(x)
}

.rdm_query <- function(effects, codec = "symmetric_packed") {
  q <- length(effects)
  if (q < 2L) stop("An RDM requires at least two experimental effects.",
    call. = FALSE)
  pairs <- utils::combn(seq_len(q), 2L)
  width <- if (identical(codec, "symmetric_packed")) {
    q * (q + 1L) / 2L
  } else {
    q * q
  }
  query <- matrix(0, width, ncol(pairs))
  labels <- character(ncol(pairs))
  for (pair in seq_len(ncol(pairs))) {
    difference <- numeric(q)
    difference[pairs[1L, pair]] <- 1
    difference[pairs[2L, pair]] <- -1
    operator <- tcrossprod(difference)
    query[, pair] <- if (identical(codec, "symmetric_packed")) {
      .svec_symmetric(operator)
    } else {
      as.vector(operator)
    }
    labels[[pair]] <- paste(effects[pairs[, pair]], collapse = " - ")
  }
  colnames(query) <- labels
  list(query = query, pairs = data.frame(
    left = effects[pairs[1L, ]], right = effects[pairs[2L, ]],
    stringsAsFactors = FALSE
  ))
}

#' Read squared experimental distances from geometry
#'
#' The returned values are point estimates. For equal-weight all-partition-
#' pairs geometry with a common fixed metric, [rdm_sampling_covariance()] can
#' construct a separate analytic within-measurement covariance law when the
#' plan was built from `lm_relation_fit()` and its residual error channel is
#' still available. Pair rows that share a partition are not independent
#' replicates for a spread-across-pairs standard error.
#'
#' @param x An `effect_geometry_plan` or a complete effect form carrying the
#'   symmetric self-form capability.
#' @param component One of `total`, `coherent`, or `configuration`.
#' @return An `effect_rdm_view`; rows are spatial measurements and columns are
#'   unique experimental pairs. Cross-generalized distances may be negative.
#' @export
rdm <- function(x, component = c("total", "coherent", "configuration")) {
  source <- .self_geometry_source(x, "RDM")
  component <- match.arg(component)
  compiled <- .rdm_query(source$effects, source$codec)
  view <- if (source$kind == "plan") {
    evaluate_geometry(x, query = compiled$query, component = component)
  } else {
    query_geometry(x, compiled$query, component)
  }
  structure(
    list(
      values = view$values,
      pairs = compiled$pairs,
      component = component,
      index = view$index,
      receipt = view$receipt
    ),
    class = "effect_rdm_view"
  )
}

.validate_rdm_models <- function(models, effects, label) {
  q <- length(effects)
  if (is.null(models)) return(list())
  if (is.matrix(models)) models <- list(model = models)
  if (!is.list(models) || length(models) < 1L) {
    stop(sprintf("`%s` must be a matrix or nonempty named list of matrices.",
      label), call. = FALSE)
  }
  if (is.null(names(models)) || anyNA(names(models)) ||
      any(!nzchar(names(models))) || anyDuplicated(names(models))) {
    stop(sprintf("`%s` must have unique nonempty names.", label), call. = FALSE)
  }
  lapply(models, function(value) {
    if (!is.matrix(value) || !is.numeric(value) ||
        !identical(dim(value), c(q, q)) || any(!is.finite(value)) ||
        max(abs(value - t(value))) > 1e-12 ||
        max(abs(diag(value))) > 1e-12) {
      stop(sprintf("Every `%s` RDM must be finite, symmetric, q-by-q, and zero-diagonal.",
        label), call. = FALSE)
    }
    row_ids <- rownames(value)
    column_ids <- colnames(value)
    if (!is.null(row_ids) || !is.null(column_ids)) {
      if (is.null(row_ids) || is.null(column_ids) ||
          anyNA(row_ids) || anyNA(column_ids) ||
          anyDuplicated(row_ids) || anyDuplicated(column_ids) ||
          !setequal(row_ids, effects) || !setequal(column_ids, effects)) {
        stop(sprintf("Named `%s` RDM axes must identify every effect exactly once.",
          label), call. = FALSE)
      }
      value <- value[effects, effects, drop = FALSE]
    }
    value
  })
}

.rdm_vector <- function(value) {
  pairs <- utils::combn(seq_len(nrow(value)), 2L)
  value[cbind(pairs[1L, ], pairs[2L, ])]
}

#' Fit multiple-regression RSA as one compiled geometry query
#'
#' @param x An `effect_geometry_plan` or a complete effect form carrying the
#'   symmetric self-form capability.
#' @param models Named model RDMs, each matching the experimental dimension.
#' @param nuisance Optional named nuisance RDMs.
#' @param intercept Whether to include an intercept in RDM space.
#' @param component Geometry component to read.
#' @return An `effect_rsa_view` with one coefficient per requested model,
#'   nuisance model, and optional intercept.
#' @export
rsa <- function(x, models, nuisance = NULL, intercept = TRUE,
                component = c("total", "coherent", "configuration")) {
  source <- .self_geometry_source(x, "RSA")
  if (!is.logical(intercept) || length(intercept) != 1L || is.na(intercept)) {
    stop("`intercept` must be TRUE or FALSE.", call. = FALSE)
  }
  component <- match.arg(component)
  q <- length(source$effects)
  models <- .validate_rdm_models(models, source$effects, "models")
  nuisance <- .validate_rdm_models(nuisance, source$effects, "nuisance")
  if (any(names(models) %in% names(nuisance))) {
    stop("Model and nuisance names must be distinct.", call. = FALSE)
  }
  predictors <- c(models, nuisance)
  design <- do.call(cbind, lapply(predictors, .rdm_vector))
  colnames(design) <- names(predictors)
  roles <- c(rep("model", length(models)), rep("nuisance", length(nuisance)))
  if (intercept) {
    design <- cbind(`(Intercept)` = 1, design)
    roles <- c("intercept", roles)
  }
  qr_design <- qr(design, LAPACK = FALSE)
  if (qr_design$rank != ncol(design)) {
    stop("The RSA design is rank deficient; remove redundant RDMs.",
      call. = FALSE)
  }
  rdm_compiled <- .rdm_query(source$effects, source$codec)
  coefficient_map <- qr.coef(qr_design, diag(nrow(design)))
  query <- rdm_compiled$query %*% t(coefficient_map)
  colnames(query) <- colnames(design)
  view <- if (source$kind == "plan") {
    evaluate_geometry(x, query = query, component = component)
  } else {
    query_geometry(x, query, component)
  }

  structure(
    list(
      coefficients = view$values,
      terms = data.frame(term = colnames(design), role = roles,
        stringsAsFactors = FALSE),
      component = component,
      index = view$index,
      receipt = view$receipt,
      query = query
    ),
    class = "effect_rsa_view"
  )
}

.unsvec_symmetric <- function(value, q) {
  if (!is.numeric(value) || length(value) != q * (q + 1L) / 2L ||
      any(!is.finite(value))) {
    stop("Packed geometry has the wrong width or non-finite values.",
      call. = FALSE)
  }
  out <- matrix(0, q, q)
  coordinate <- 0L
  for (column in seq_len(q)) {
    for (row in column:q) {
      coordinate <- coordinate + 1L
      entry <- value[[coordinate]] / if (row == column) 1 else sqrt(2)
      out[row, column] <- entry
      out[column, row] <- entry
    }
  }
  out
}

#' Read the signed eigenvalue spectrum of cross-generalized geometry
#'
#' @param x A complete effect form carrying the symmetric self-form capability.
#' @param component Geometry component to decompose.
#' @param row_block Positive number of measurement rows read per block.
#' @return An `effect_spectrum_view`. Eigenvalues are ordered from largest to
#'   smallest and are never truncated at zero.
#' @export
geometry_spectrum <- function(x,
                              component = c("total", "coherent", "configuration"),
                              row_block = 1024L) {
  .require_self_form_view(x, "A geometry spectrum")
  component <- match.arg(component)
  row_block <- .validate_tile_size(row_block, "row_block")
  .validate_effect_form(x)
  .require_effect_form_component(x, component)
  q <- length(x$effects)
  values <- matrix(0, x$total$dim[[1L]], q,
    dimnames = list(NULL, paste0("root", seq_len(q))))
  for (start in .tile_starts(nrow(values), row_block)) {
    rows <- start:min(start + row_block - 1L, nrow(values))
    packed <- .geometry_component_validated(x, component, rows)
    for (position in seq_along(rows)) {
      form <- if (identical(x$codec, "symmetric_packed")) {
        .unsvec_symmetric(packed[position, ], q)
      } else {
        matrix(packed[position, ], q, q)
      }
      values[rows[[position]], ] <- eigen(
        form, symmetric = TRUE,
        only.values = TRUE
      )$values
    }
  }
  structure(
    list(values = values, component = component, index = x$index,
      indefinite_estimates_preserved = TRUE, receipt = x$receipt),
    class = "effect_spectrum_view"
  )
}
