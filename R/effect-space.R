# Canonical experimental-coordinate identity --------------------------------

#' Define an experimental coordinate space
#'
#' An effect space gives experimental coordinates a semantic identity beyond
#' their display labels. Its stable signature covers coordinate order, basis,
#' units, scaling, and provenance.
#'
#' @param coordinates Unique ordered coordinate identifiers.
#' @param basis_id One nonempty semantic basis identifier.
#' @param units One unit shared by all coordinates, or one per coordinate.
#' @param scale One positive finite scale shared by all coordinates, or one per
#'   coordinate.
#' @param provenance Optional design or contrast provenance as a list.
#' @return An immutable-by-convention `effect_space` value.
#' @export
effect_space <- function(coordinates, basis_id = "unspecified",
                         units = "arbitrary", scale = 1,
                         provenance = list()) {
  if (!is.character(coordinates) || length(coordinates) < 1L) {
    stop("`coordinates` must contain at least one effect coordinate.",
      call. = FALSE)
  }
  coordinates <- .validate_effect_names(coordinates, length(coordinates))
  if (!is.character(basis_id) || length(basis_id) != 1L || is.na(basis_id) ||
      !nzchar(basis_id)) {
    stop("`basis_id` must be one nonempty identifier.", call. = FALSE)
  }
  if (!is.character(units) || !length(units) %in% c(1L, length(coordinates)) ||
      anyNA(units) || any(!nzchar(units))) {
    stop("`units` must provide one nonempty unit or one per coordinate.",
      call. = FALSE)
  }
  if (length(units) == 1L) units <- rep(units, length(coordinates))
  if (!is.numeric(scale) || !length(scale) %in% c(1L, length(coordinates)) ||
      anyNA(scale) || any(!is.finite(scale)) || any(scale <= 0)) {
    stop("`scale` must provide one positive finite value or one per coordinate.",
      call. = FALSE)
  }
  if (length(scale) == 1L) scale <- rep(as.numeric(scale), length(coordinates))
  .validate_effect_provenance(provenance)

  units <- stats::setNames(unname(units), coordinates)
  scale <- stats::setNames(as.numeric(scale), coordinates)
  semantic <- list(
    schema_version = 1L,
    coordinates = unname(coordinates),
    basis_id = unname(basis_id),
    units = units,
    scale = scale,
    provenance = provenance
  )
  signature <- paste0(
    "sha256:", digest::digest(semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L)
  )
  structure(
    c(semantic[-1L], list(signature = signature)),
    class = "effect_space"
  )
}

.validate_effect_provenance <- function(x, path = "provenance") {
  if (!is.list(x)) {
    stop("`provenance` must be a list.", call. = FALSE)
  }
  inspect <- function(value, location) {
    if (is.environment(value) || is.function(value) ||
        typeof(value) %in% c("externalptr", "weakref")) {
      stop(sprintf("`%s` contains a nonportable value.", location), call. = FALSE)
    }
    if (is.list(value)) {
      for (index in seq_along(value)) {
        label <- names(value)[index]
        if (is.null(label) || is.na(label) || !nzchar(label)) label <- index
        inspect(value[[index]], paste0(location, "$", label))
      }
    }
    invisible(NULL)
  }
  inspect(x, path)
  invisible(x)
}

.validate_effect_space <- function(x) {
  expected <- c("coordinates", "basis_id", "units", "scale", "provenance",
    "signature")
  if (!inherits(x, "effect_space") || !is.list(x) ||
      !identical(names(x), expected)) {
    stop("Effect-space fields are missing or noncanonical.", call. = FALSE)
  }
  rebuilt <- effect_space(x$coordinates, x$basis_id, x$units, x$scale,
    x$provenance)
  if (!identical(x, rebuilt)) {
    stop("Effect-space metadata or signature is inconsistent.", call. = FALSE)
  }
  rebuilt
}

.as_effect_space <- function(x, expected = NULL) {
  if (is.null(x) && !is.null(expected)) {
    x <- paste0("effect", seq_len(expected))
  }
  value <- if (inherits(x, "effect_space")) {
    .validate_effect_space(x)
  } else {
    effect_space(x)
  }
  if (!is.null(expected) && length(value$coordinates) != expected) {
    stop("The effect-space dimension does not match the declared data.",
      call. = FALSE)
  }
  value
}

.same_effect_space <- function(x, y) {
  x <- .validate_effect_space(x)
  y <- .validate_effect_space(y)
  identical(x$signature, y$signature) && identical(x, y)
}

.effect_names <- function(x) .validate_effect_space(x)$coordinates
