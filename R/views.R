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
#' @param x A complete `effect_geometry`.
#' @param weights One finite contrast weight per named experimental effect.
#' @return An `effect_contrast_view` containing signed marginals and the three
#'   energy components. `coherence_fraction` is reported only where the raw
#'   cross-generalized components form a nonnegative partition.
#' @export
contrast <- function(x, weights) {
  if (!inherits(x, "effect_geometry")) {
    stop("`x` must be a complete effect_geometry.", call. = FALSE)
  }
  .validate_effect_geometry(x, probe = FALSE)
  weights <- .align_contrast(weights, x$effects)
  query <- bilinear_query(tcrossprod(weights))
  total <- drop(query_geometry(x, query, "total")$values)
  coherent <- drop(query_geometry(x, query, "coherent")$values)
  configuration <- total - coherent
  signed <- lapply(x$marginals, function(value) drop(value %*% weights))
  if (identical(attr(x$marginals, "semantics"), "undirected_endpoint")) {
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
      index = x$index,
      receipt = x$receipt
    ),
    class = "effect_contrast_view"
  )
}

.rdm_query <- function(effects) {
  q <- length(effects)
  if (q < 2L) stop("An RDM requires at least two experimental effects.",
    call. = FALSE)
  pairs <- utils::combn(seq_len(q), 2L)
  query <- matrix(0, q * (q + 1L) / 2L, ncol(pairs))
  labels <- character(ncol(pairs))
  for (pair in seq_len(ncol(pairs))) {
    difference <- numeric(q)
    difference[pairs[1L, pair]] <- 1
    difference[pairs[2L, pair]] <- -1
    query[, pair] <- .svec_symmetric(tcrossprod(difference))
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
#' @param x A complete `effect_geometry`.
#' @param component One of `total`, `coherent`, or `configuration`.
#' @return An `effect_rdm_view`; rows are spatial measurements and columns are
#'   unique experimental pairs. Cross-generalized distances may be negative.
#' @export
rdm <- function(x, component = c("total", "coherent", "configuration")) {
  if (!inherits(x, "effect_geometry")) {
    stop("`x` must be a complete effect_geometry.", call. = FALSE)
  }
  .validate_effect_geometry(x, probe = FALSE)
  component <- match.arg(component)
  compiled <- .rdm_query(x$effects)
  view <- query_geometry(x, compiled$query, component)
  structure(
    list(
      values = view$values,
      pairs = compiled$pairs,
      component = component,
      index = x$index,
      receipt = x$receipt
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
#' @param x A complete `effect_geometry`.
#' @param models Named model RDMs, each matching the experimental dimension.
#' @param nuisance Optional named nuisance RDMs.
#' @param intercept Whether to include an intercept in RDM space.
#' @param component Geometry component to read.
#' @return An `effect_rsa_view` with one coefficient per requested model,
#'   nuisance model, and optional intercept.
#' @export
rsa <- function(x, models, nuisance = NULL, intercept = TRUE,
                component = c("total", "coherent", "configuration")) {
  if (!inherits(x, "effect_geometry")) {
    stop("`x` must be a complete effect_geometry.", call. = FALSE)
  }
  .validate_effect_geometry(x, probe = FALSE)
  if (!is.logical(intercept) || length(intercept) != 1L || is.na(intercept)) {
    stop("`intercept` must be TRUE or FALSE.", call. = FALSE)
  }
  component <- match.arg(component)
  q <- length(x$effects)
  models <- .validate_rdm_models(models, x$effects, "models")
  nuisance <- .validate_rdm_models(nuisance, x$effects, "nuisance")
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
  rdm_compiled <- .rdm_query(x$effects)
  coefficient_map <- qr.coef(qr_design, diag(nrow(design)))
  query <- rdm_compiled$query %*% t(coefficient_map)
  colnames(query) <- colnames(design)
  view <- query_geometry(x, query, component)

  structure(
    list(
      coefficients = view$values,
      terms = data.frame(term = colnames(design), role = roles,
        stringsAsFactors = FALSE),
      component = component,
      index = x$index,
      receipt = x$receipt,
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
#' @param x A complete `effect_geometry`.
#' @param component Geometry component to decompose.
#' @param row_block Positive number of measurement rows read per block.
#' @return An `effect_spectrum_view`. Eigenvalues are ordered from largest to
#'   smallest and are never truncated at zero.
#' @export
geometry_spectrum <- function(x,
                              component = c("total", "coherent", "configuration"),
                              row_block = 1024L) {
  if (!inherits(x, "effect_geometry")) {
    stop("`x` must be a complete effect_geometry.", call. = FALSE)
  }
  component <- match.arg(component)
  row_block <- .validate_tile_size(row_block, "row_block")
  .validate_effect_geometry(x)
  q <- length(x$effects)
  values <- matrix(0, x$total$dim[[1L]], q,
    dimnames = list(NULL, paste0("root", seq_len(q))))
  for (start in .tile_starts(nrow(values), row_block)) {
    rows <- start:min(start + row_block - 1L, nrow(values))
    packed <- .geometry_component_validated(x, component, rows)
    for (position in seq_along(rows)) {
      values[rows[[position]], ] <- eigen(
        .unsvec_symmetric(packed[position, ], q), symmetric = TRUE,
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
