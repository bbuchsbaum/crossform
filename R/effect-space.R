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
#' @return An `effect_space`: a list with `$coordinates` (the ordered names),
#'   `$basis_id`, `$units` and `$scale` (both named by coordinate),
#'   `$provenance`, and a `$signature` that changes whenever any of those does.
#'   Treat the value as immutable.
#' @family studies and effect maps
#' @seealso [condition_space()] for the mean-model vocabulary an
#'   [effect_map()] is declared against, and [effect_extractor()] for the map
#'   that produces coordinates in this space.
#' @examples
#' # Name the coordinates a fit will produce, and the basis they live in.
#' space <- effect_space(
#'   c("face", "house"),
#'   basis_id = "condition-means:v1", units = "arbitrary-BOLD"
#' )
#' space$coordinates
#' space$units[["face"]]
#'
#' # The signature covers meaning, not display labels: the same two names in
#' # the default unspecified basis are a different space, so objects built
#' # against one cannot be silently combined with the other.
#' identical(space$signature, effect_space(c("face", "house"))$signature)
#' @export
effect_space <- function(coordinates, basis_id = "unspecified",
                         units = "arbitrary", scale = 1,
                         provenance = list()) {
  if (missing(coordinates)) {
    stop(paste0(
      "`coordinates` is required: name the experimental effects this space ",
      "holds, for example `effect_space(c(\"face\", \"house\"))`."
    ), call. = FALSE)
  }
  if (!is.character(coordinates) || length(coordinates) < 1L) {
    stop(sprintf(paste0(
      "`coordinates` must be a nonempty character vector naming the ",
      "experimental effects; received %s."
    ), .msg_value(coordinates)), call. = FALSE)
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
