# Condition-vocabulary effect maps ------------------------------------------

.validate_positive_scale <- function(value, size, name = "scale") {
  if (!.is_finite_numeric(value) || !length(value) %in% c(1L, size) ||
      anyNA(value) || any(value <= 0)) {
    .input_error(sprintf(
      "`%s` must provide one positive finite value or one per coordinate.",
      name
    ))
  }
  if (length(value) == 1L) value <- rep(as.numeric(value), size)
  as.numeric(value)
}

.validate_units <- function(value, size, name = "units") {
  if (!.is_strings(value) || !length(value) %in% c(1L, size)) {
    .input_error(sprintf(
      "`%s` must provide one nonempty unit or one per coordinate.", name
    ))
  }
  if (length(value) == 1L) value <- rep(value, size)
  unname(value)
}

#' Define a semantic condition space
#'
#' A condition space identifies the named mean-model coordinates against which
#' an [effect_map()] is declared. It is independent of design-matrix coding and
#' column order.
#'
#' @param coordinates Unique ordered semantic coordinates, such as condition
#'   names or condition-by-basis names.
#' @param basis_id One identifier for the targeted mean-model basis.
#' @param units One unit shared by all coordinates, or one per coordinate.
#' @param scale One positive scale shared by all coordinates, or one per
#'   coordinate.
#' @param provenance Portable semantic provenance.
#' @return An `effect_condition_space`: a list with `$coordinates`,
#'   `$basis_id`, `$units` and `$scale` (named by coordinate), `$provenance`,
#'   and a `$signature`. Two objects bind the same condition space only when
#'   their signatures agree.
#' @family studies and effect maps
#' @seealso [effect_map()] to declare effects in this vocabulary and
#'   [coefficient_parameterization()] to bind it to a compiled design.
#' @examples
#' # The scientific vocabulary: three condition means in BOLD units.
#' conditions <- condition_space(
#'   c("face", "body", "tool"),
#'   basis_id = "scan-level-condition-mean:v1", units = "arbitrary-BOLD"
#' )
#' conditions
#' conditions$coordinates
#'
#' # These coordinates are semantic, not design columns: the same space can be
#' # reached by cell-means or treatment coding, and the shared signature is
#' # what lets an effect map and a design model be checked against each other.
#' substr(conditions$signature, 1, 24)
#' @export
condition_space <- function(coordinates, basis_id = "condition-means",
                            units = "arbitrary", scale = 1,
                            provenance = list()) {
  if (!is.character(coordinates) || length(coordinates) < 1L) {
    .input_error("`coordinates` must contain at least one semantic coordinate.")
  }
  coordinates <- .validate_effect_names(coordinates, length(coordinates))
  basis_id <- .validate_nonempty_id(basis_id, "basis_id")
  units <- .validate_units(units, length(coordinates))
  scale <- .validate_positive_scale(scale, length(coordinates))
  .validate_effect_provenance(provenance)

  units <- stats::setNames(units, coordinates)
  scale <- stats::setNames(scale, coordinates)
  semantic <- list(
    schema_version = 1L,
    coordinates = unname(coordinates),
    basis_id = basis_id,
    units = units,
    scale = scale,
    provenance = provenance
  )
  structure(c(semantic[-1L], list(
    signature = .sha256_signature(semantic, "condition-space-sha256:")
  )), class = "effect_condition_space")
}

.validate_condition_space <- function(value) {
  expected <- c(
    "coordinates", "basis_id", "units", "scale", "provenance", "signature"
  )
  if (!.sealed_fields(value, "effect_condition_space", expected)) {
    .input_error("Condition-space fields are missing or noncanonical.")
  }
  rebuilt <- condition_space(
    value$coordinates, value$basis_id, value$units, value$scale,
    value$provenance
  )
  if (!identical(value, rebuilt)) {
    .contract_error("Condition-space metadata or signature is inconsistent.")
  }
  rebuilt
}

.as_condition_space <- function(value, expected = NULL) {
  result <- if (inherits(value, "effect_condition_space")) {
    .validate_condition_space(value)
  } else {
    condition_space(value)
  }
  if (!is.null(expected) && length(result$coordinates) != expected) {
    .contract_error(
      "The condition-space dimension does not match the declared map."
    )
  }
  result
}

#' Declare effects in a semantic condition vocabulary
#'
#' The rows of `weights` are named scientific effects and the columns are bound
#' semantic condition coordinates. The returned identity never contains a
#' compiled design matrix or coefficient coding.
#'
#' @param weights A finite effect-by-condition matrix.
#' @param conditions A [condition_space()] or condition names. Defaults to the
#'   column names of `weights`.
#' @param effects An [effect_space()] or output effect names. Defaults to the
#'   row names of `weights`.
#' @param units Output units. By default the shared condition-space unit.
#' @param scale Output scale.
#' @param component One targeted temporal-basis component.
#' @param provenance Portable provenance for the declared functional.
#' @return An `effect_condition_map`: a list with the effect-by-condition
#'   `$weights`, the bound `$condition_space` and derived `$effect_space`, the
#'   targeted `$component`, and an `$effect_map_id` covering all of them.
#' @family studies and effect maps
#' @seealso [condition_space()] for the input vocabulary,
#'   [coefficient_parameterization()] plus [lower_effect_map()] to reach a
#'   coefficient axis, and [plan_relation()] to request these effects from a
#'   study.
#' @examples
#' conditions <- condition_space(c("face", "body", "tool"))
#'
#' # Request the three condition means and one named contrast at once, so
#' # later views can form the contrast without refitting the run models.
#' weights <- rbind(
#'   face = c(1, 0, 0), body = c(0, 1, 0), tool = c(0, 0, 1),
#'   `face-body` = c(1, -1, 0)
#' )
#' colnames(weights) <- conditions$coordinates
#' effects <- effect_map(weights, conditions)
#' effects
#' effects$weights["face-body", ]
#'
#' # The derived effect space inherits the condition units, and its identity
#' # records the functional rather than any design coding.
#' effects$effect_space$basis_id
#' @export
effect_map <- function(weights, conditions = colnames(weights),
                       effects = rownames(weights), units = NULL, scale = 1,
                       component = "canonical-amplitude",
                       provenance = list()) {
  if (!.is_finite_matrix(weights) || any(dim(weights) < 1L)) {
    .input_error(
      "`weights` must be a finite nonempty effect-by-condition matrix."
    )
  }
  if (is.null(conditions)) {
    .input_error("Supply a condition space or name every column of `weights`.")
  }
  conditions <- .as_condition_space(conditions, ncol(weights))
  if (!is.null(colnames(weights)) &&
      !identical(colnames(weights), conditions$coordinates)) {
    .input_error(
      "`weights` columns must follow the condition-space coordinates."
    )
  }
  colnames(weights) <- conditions$coordinates
  component <- .validate_nonempty_id(component, "component")
  .validate_effect_provenance(provenance)

  if (inherits(effects, "effect_space")) {
    if (!is.null(units) || !identical(scale, 1)) {
      .input_error("Do not supply `units` or `scale` with an `effect_space`.")
    }
    effects <- .as_effect_space(effects, nrow(weights))
  } else {
    if (is.null(effects)) {
      .input_error("Supply output effect names or name every row of `weights`.")
    }
    effects <- .validate_effect_names(effects, nrow(weights))
    if (is.null(units)) {
      shared_units <- unique(unname(conditions$units))
      if (length(shared_units) != 1L) {
        .input_error(
          "Supply `units` when condition coordinates use different units."
        )
      }
      units <- shared_units
    }
    functional_revision <- .sha256_signature(list(
      rows = effects,
      columns = conditions$coordinates,
      weights = unname(weights)
    ))
    effects <- effect_space(
      effects,
      basis_id = paste0("condition-functional:", component),
      units = units,
      scale = scale,
      provenance = list(
        condition_space_signature = conditions$signature,
        functional_revision = functional_revision,
        declaration = provenance
      )
    )
  }
  rownames(weights) <- effects$coordinates

  semantic <- list(
    schema_version = 1L,
    condition_space = conditions,
    effect_space = effects,
    weights = weights,
    dimnames = dimnames(weights),
    component = component,
    provenance = provenance
  )
  structure(c(semantic[-1L], list(
    effect_map_id = .sha256_signature(semantic, "effect-map-sha256:")
  )), class = "effect_condition_map")
}

.validate_effect_map <- function(value) {
  expected <- c(
    "condition_space", "effect_space", "weights", "dimnames", "component",
    "provenance", "effect_map_id"
  )
  if (!.sealed_fields(value, "effect_condition_map", expected)) {
    .input_error("Effect-map fields are missing or noncanonical.")
  }
  conditions <- .validate_condition_space(value$condition_space)
  effects <- .validate_effect_space(value$effect_space)
  weights <- value$weights
  if (!is.matrix(weights)) {
    weights <- matrix(weights, nrow = length(effects$coordinates))
  }
  dimnames(weights) <- value$dimnames
  rebuilt <- effect_map(
    weights,
    conditions = conditions,
    effects = effects,
    component = value$component,
    provenance = value$provenance
  )
  if (!identical(value, rebuilt)) {
    .contract_error("Effect-map metadata or identity is inconsistent.")
  }
  rebuilt
}

#' Identify one concrete coefficient parameterization
#'
#' `map` declares `mu = map %*% beta`, where `mu` follows the semantic
#' condition coordinates and `beta` follows the compiled coefficient axis.
#' The object belongs to a design receipt, not scientific plan identity.
#'
#' @param map A finite condition-by-coefficient matrix.
#' @param conditions The bound [condition_space()].
#' @param coefficients Unique ordered coefficient names. Defaults to `map`
#'   column names.
#' @param coding_id One coding identifier recorded in the receipt.
#' @param provenance Portable compiler provenance.
#' @param tolerance Positive rank tolerance.
#' @return An `effect_coefficient_parameterization`: a list with the
#'   condition-by-coefficient `$map`, the bound `$condition_space`, the
#'   `$coefficients` axis, the `$coding_id`, the achieved `$semantic_rank`, and
#'   a `$parameterization_id`.
#' @family studies and effect maps
#' @seealso [lower_effect_map()] to apply it to an [effect_map()], and
#'   [design_model()] which carries one parameterization per partition.
#' @examples
#' conditions <- condition_space(c("face", "body", "tool"))
#'
#' # The compiled design adds a drift column. This declares that the three
#' # condition means are the first three coefficients and ignore drift, so the
#' # effect request survives a change of coding.
#' map <- cbind(diag(3), drift = 0)
#' dimnames(map) <- list(
#'   conditions$coordinates, c(conditions$coordinates, "drift")
#' )
#' coding <- coefficient_parameterization(
#'   map, conditions, coding_id = "cell-means-plus-drift"
#' )
#' coding
#' coding$semantic_rank
#'
#' # A coding that cannot identify every condition is refused, not rounded.
#' collapsed <- map
#' collapsed["body", ] <- collapsed["face", ]
#' refusal <- catch_refusal(coefficient_parameterization(
#'   collapsed, conditions, coding_id = "face-and-body-collapsed"
#' ))
#' refusal$capability
#' @export
coefficient_parameterization <- function(
    map, conditions, coefficients = colnames(map), coding_id,
    provenance = list(), tolerance = sqrt(.Machine$double.eps)) {
  if (!.is_finite_matrix(map) || any(dim(map) < 1L)) {
    .input_error(
      "`map` must be a finite nonempty condition-by-coefficient matrix."
    )
  }
  conditions <- .as_condition_space(conditions, nrow(map))
  if (!is.null(rownames(map)) &&
      !identical(rownames(map), conditions$coordinates)) {
    .input_error("`map` rows must follow the condition-space coordinates.")
  }
  if (is.null(coefficients)) {
    .input_error(
      "Supply unique coefficient names or name every column of `map`."
    )
  }
  coefficients <- .validate_effect_names(coefficients, ncol(map))
  if (!is.null(colnames(map)) && !identical(colnames(map), coefficients)) {
    .input_error("`map` columns must follow `coefficients`.")
  }
  coding_id <- .validate_nonempty_id(coding_id, "coding_id")
  .check_number(tolerance, "tolerance", positive = TRUE)
  .validate_effect_provenance(provenance)
  rownames(map) <- conditions$coordinates
  colnames(map) <- coefficients
  rank <- qr(map, tol = tolerance, LAPACK = FALSE)$rank
  if (rank != nrow(map)) {
    .capability_refusal(
      paste0(
        "The coefficient parameterization identifies only ", rank, " of ",
        nrow(map), " semantic coordinates."
      ),
      capability = "valid_effect_lowering",
      namespace = "relation_plan",
      reasons = sprintf(
        "Coding `%s` has semantic rank %d for %d coordinates.",
        coding_id, rank, nrow(map)
      ),
      remedies = paste0(
        "Supply a full-row-rank semantic-to-coefficient map or narrow the ",
        "declared condition space."
      )
    )
  }

  semantic <- list(
    schema_version = 1L,
    condition_space = conditions,
    coefficients = coefficients,
    map = map,
    dimnames = dimnames(map),
    coding_id = coding_id,
    provenance = provenance,
    tolerance = as.numeric(tolerance),
    semantic_rank = as.integer(rank)
  )
  structure(c(semantic[-1L], list(
    parameterization_id = .sha256_signature(semantic, "coefficient-parameterization-sha256:")
  )), class = "effect_coefficient_parameterization")
}

.validate_coefficient_parameterization <- function(value) {
  expected <- c(
    "condition_space", "coefficients", "map", "dimnames", "coding_id",
    "provenance", "tolerance", "semantic_rank", "parameterization_id"
  )
  if (!.sealed_fields(value, "effect_coefficient_parameterization", expected)) {
    .input_error(
      "Coefficient-parameterization fields are missing or noncanonical."
    )
  }
  conditions <- .validate_condition_space(value$condition_space)
  map <- value$map
  if (!is.matrix(map)) {
    map <- matrix(map, nrow = length(conditions$coordinates))
  }
  dimnames(map) <- value$dimnames
  rebuilt <- coefficient_parameterization(
    map,
    conditions = conditions,
    coefficients = value$coefficients,
    coding_id = value$coding_id,
    provenance = value$provenance,
    tolerance = value$tolerance
  )
  if (!identical(value, rebuilt)) {
    .contract_error(
      "Coefficient-parameterization metadata or identity is inconsistent."
    )
  }
  rebuilt
}

#' Lower semantic effects into a concrete coefficient basis
#'
#' @param effects An [effect_map()].
#' @param parameterization A [coefficient_parameterization()] for the same
#'   semantic condition space.
#' @return An `effect_lowered_map`: a list with the effect-by-coefficient
#'   `$target`, the `$effect_space` and `$condition_space` it came from, the
#'   `$effect_map_id`, `$parameterization_id`, and `$coding_id` receipt fields,
#'   `$capabilities` (all `TRUE` on this route), and a `$lowering_id`.
#' @family studies and effect maps
#' @seealso [raw_effect_map()] for the degenerate route with no condition
#'   space, and [plan_relation()] which lowers effects for every partition.
#' @examples
#' conditions <- condition_space(c("face", "body"))
#' weights <- rbind(`face-body` = c(1, -1))
#' colnames(weights) <- conditions$coordinates
#' effects <- effect_map(weights, conditions)
#'
#' # A design with a drift nuisance column, and the coding that names which
#' # coefficients carry the condition means.
#' map <- cbind(diag(2), drift = 0)
#' dimnames(map) <- list(
#'   conditions$coordinates, c(conditions$coordinates, "drift")
#' )
#' coding <- coefficient_parameterization(
#'   map, conditions, coding_id = "cell-means-plus-drift"
#' )
#'
#' # Lowering moves the request from condition names onto the design's
#' # coefficient axis without changing what was asked for.
#' lowered <- lower_effect_map(effects, coding)
#' lowered$target
#' @export
lower_effect_map <- function(effects, parameterization) {
  effects <- .validate_effect_map(effects)
  parameterization <- .validate_coefficient_parameterization(parameterization)
  if (!identical(
      effects$condition_space$signature,
      parameterization$condition_space$signature
    )) {
    .capability_refusal(
      "The effect and coefficient parameterization bind different condition spaces.",
      capability = "valid_effect_lowering",
      namespace = "relation_plan",
      reasons = c(
        paste0("Effect map: ", effects$condition_space$signature, "."),
        paste0(
          "Parameterization: ",
          parameterization$condition_space$signature, "."
        )
      ),
      remedies = "Compile both objects against the same condition space."
    )
  }
  target <- effects$weights %*% parameterization$map
  dimnames(target) <- list(
    effects$effect_space$coordinates,
    parameterization$coefficients
  )
  semantic <- list(
    schema_version = 1L,
    target = target,
    dimnames = dimnames(target),
    effect_space = effects$effect_space,
    condition_space = effects$condition_space,
    effect_map_id = effects$effect_map_id,
    parameterization_id = parameterization$parameterization_id,
    coding_id = parameterization$coding_id,
    capabilities = list(
      symbolic_effects = TRUE,
      valid_effect_lowering = TRUE,
      coding_invariant = TRUE
    )
  )
  structure(c(semantic[-1L], list(
    lowering_id = .sha256_signature(semantic, "effect-lowering-sha256:")
  )), class = "effect_lowered_map")
}

.validate_lowered_effect_map <- function(value) {
  common <- c(
    "target", "dimnames", "effect_space", "condition_space",
    "effect_map_id", "parameterization_id", "coding_id", "capabilities"
  )
  expected <- if (inherits(value, "effect_raw_map")) {
    c(common, "provenance", "lowering_id")
  } else {
    c(common, "lowering_id")
  }
  if (!.sealed_fields(value, "effect_lowered_map", expected)) {
    .input_error("Lowered-effect-map fields are missing or noncanonical.")
  }
  effects <- .validate_effect_space(value$effect_space)
  target <- value$target
  if (!.is_finite_matrix(target) ||
      nrow(target) != length(effects$coordinates) ||
      !identical(dimnames(target), value$dimnames) ||
      !identical(rownames(target), effects$coordinates)) {
    .input_error("Lowered target values or axes are inconsistent.")
  }
  .validate_nonempty_id(value$effect_map_id, "effect_map_id")
  .validate_nonempty_id(value$coding_id, "coding_id")
  if (!is.list(value$capabilities) ||
      !identical(names(value$capabilities), c(
        "symbolic_effects", "valid_effect_lowering", "coding_invariant"
      )) || !all(vapply(value$capabilities, is.logical, logical(1))) ||
      !all(lengths(value$capabilities) == 1L)) {
    .input_error("Lowered-effect capabilities are missing or noncanonical.")
  }
  raw <- inherits(value, "effect_raw_map")
  if (raw) {
    if (!is.null(value$condition_space) ||
        !is.null(value$parameterization_id) ||
        any(unlist(value$capabilities, use.names = FALSE))) {
      .input_error("Raw targets cannot claim symbolic lowering capabilities.")
    }
    .validate_effect_provenance(value$provenance)
    expected_map_id <- .sha256_signature(list(
      target = unname(target),
      dimnames = dimnames(target),
      effect_space = effects,
      provenance = value$provenance
    ), "raw-effect-map-sha256:")
    if (!identical(value$effect_map_id, expected_map_id)) {
      .contract_error("Raw effect-map identity is inconsistent.")
    }
  } else {
    .validate_condition_space(value$condition_space)
    .validate_nonempty_id(value$parameterization_id, "parameterization_id")
    if (!all(unlist(value$capabilities, use.names = FALSE))) {
      .input_error("A semantic lowering is missing a construction guarantee.")
    }
  }
  semantic <- c(list(schema_version = 1L), value[setdiff(
    names(value), "lowering_id"
  )])
  prefix <- if (raw) {
    "raw-effect-lowering-sha256:"
  } else {
    "effect-lowering-sha256:"
  }
  if (!identical(value$lowering_id, .sha256_signature(semantic, prefix))) {
    .contract_error("Lowered-effect-map identity is inconsistent.")
  }
  value
}

#' Bind a raw coefficient-space target
#'
#' This is the honest degenerate route for a supplied `T` matrix. Its identity
#' includes the coefficient order and target values, and it cannot claim
#' symbolic or coding-invariant effects.
#'
#' @param target A finite effect-by-coefficient target matrix.
#' @param effects An [effect_space()] or output effect names.
#' @param coefficients Unique ordered coefficient names.
#' @param units,scale Output units and scale when `effects` is a character
#'   vector.
#' @param provenance Portable target provenance.
#' @return An `effect_raw_map` (also an `effect_lowered_map`): a list with the
#'   effect-by-coefficient `$target`, its `$effect_space`, a `$coding_id` of
#'   `"raw-X-T"`, `$capabilities` with `symbolic_effects`,
#'   `valid_effect_lowering`, and `coding_invariant` all `FALSE`, and a
#'   `$lowering_id` that includes the target values themselves.
#' @family studies and effect maps
#' @seealso [lower_effect_map()] for the semantic route, and
#'   [raw_design_model()], the design model this map must be paired with.
#' @examples
#' # Use this when you have a target matrix but no condition-space meaning for
#' # the coefficient columns it multiplies.
#' target <- rbind(`face-body` = c(1, -1, 0))
#' colnames(target) <- c("beta_face", "beta_body", "beta_drift")
#' raw <- raw_effect_map(target, units = "arbitrary-BOLD")
#' raw
#'
#' # The honest cost: the numeric values are part of the identity, and no
#' # coding-invariance claim is made, unlike lower_effect_map().
#' unlist(raw$capabilities)
#' @export
raw_effect_map <- function(target, effects = rownames(target),
                           coefficients = colnames(target),
                           units = "arbitrary", scale = 1,
                           provenance = list()) {
  if (!.is_finite_matrix(target) || any(dim(target) < 1L)) {
    .input_error(
      "`target` must be a finite nonempty effect-by-coefficient matrix."
    )
  }
  if (is.null(effects)) {
    .input_error("Supply output effect names or name every row of `target`.")
  }
  effects <- if (inherits(effects, "effect_space")) {
    .as_effect_space(effects, nrow(target))
  } else {
    effect_space(effects, basis_id = "raw-coefficient-target",
      units = units, scale = scale, provenance = provenance)
  }
  if (is.null(coefficients)) {
    .input_error("Supply coefficient names or name every column of `target`.")
  }
  coefficients <- .validate_effect_names(coefficients, ncol(target))
  .validate_effect_provenance(provenance)
  dimnames(target) <- list(effects$coordinates, coefficients)
  semantic <- list(
    schema_version = 1L,
    target = target,
    dimnames = dimnames(target),
    effect_space = effects,
    condition_space = NULL,
    effect_map_id = .sha256_signature(list(
      target = unname(target),
      dimnames = dimnames(target),
      effect_space = effects,
      provenance = provenance
    ), "raw-effect-map-sha256:"),
    parameterization_id = NULL,
    coding_id = "raw-X-T",
    capabilities = list(
      symbolic_effects = FALSE,
      valid_effect_lowering = FALSE,
      coding_invariant = FALSE
    ),
    provenance = provenance
  )
  lowering_id <- .sha256_signature(semantic, "raw-effect-lowering-sha256:")
  structure(c(semantic[-1L], list(lowering_id = lowering_id)),
    class = c("effect_raw_map", "effect_lowered_map"))
}

#' @export
format.effect_condition_space <- function(x, ...) {
  x <- .validate_condition_space(x)
  sprintf("condition_space<%d; %s>", length(x$coordinates), x$basis_id)
}

#' @export
print.effect_condition_space <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
format.effect_condition_map <- function(x, ...) {
  x <- .validate_effect_map(x)
  sprintf(
    "effect_map<%d effects x %d conditions; %s>",
    nrow(x$weights), ncol(x$weights), x$component
  )
}

#' @export
print.effect_condition_map <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
format.effect_lowered_map <- function(x, ...) {
  x <- .validate_lowered_effect_map(x)
  sprintf(
    "%s<%d effects x %d coefficients; %s>",
    if (inherits(x, "effect_raw_map")) "raw_effect_map" else "lowered_effect_map",
    nrow(x$target), ncol(x$target), x$coding_id
  )
}

#' @export
print.effect_lowered_map <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
format.effect_coefficient_parameterization <- function(x, ...) {
  x <- .validate_coefficient_parameterization(x)
  sprintf(
    "coefficient_parameterization<%d conditions x %d coefficients; %s>",
    nrow(x$map), ncol(x$map), x$coding_id
  )
}

#' @export
print.effect_coefficient_parameterization <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}
