# Pair-space query constructors ---------------------------------------------

.pair_axis_index <- function(value, space, label) {
  if (!.is_strings(value)) {
    .input_error(sprintf("`%s` identifiers must be nonempty strings.", label))
  }
  index <- match(value, space$coordinates)
  if (anyNA(index)) {
    .input_error(
      sprintf("Every `%s` identifier must belong to its effect space.", label)
    )
  }
  index
}

.pair_eligibility <- function(eligible, left_space, right_space) {
  shape <- c(length(left_space$coordinates), length(right_space$coordinates))
  if (is.null(eligible)) return(matrix(TRUE, shape[[1L]], shape[[2L]],
    dimnames = list(left_space$coordinates, right_space$coordinates)))
  if (is.matrix(eligible) && is.logical(eligible) &&
      identical(dim(eligible), shape) && !anyNA(eligible)) {
    dimnames(eligible) <- list(left_space$coordinates, right_space$coordinates)
    return(eligible)
  }
  if (!is.data.frame(eligible) ||
      !identical(names(eligible), c("left", "right")) || nrow(eligible) < 1L) {
    .input_error(
      "`eligible` must be NULL, a logical pair matrix, or a left/right table."
    )
  }
  value <- matrix(FALSE, shape[[1L]], shape[[2L]],
    dimnames = list(left_space$coordinates, right_space$coordinates))
  left <- .pair_axis_index(eligible$left, left_space, "left")
  right <- .pair_axis_index(eligible$right, right_space, "right")
  value[cbind(left, right)] <- TRUE
  value
}

.new_pair_coupling <- function(kind, value, eligible, left_space, right_space) {
  structure(list(
    schema_version = 1L,
    kind = kind,
    value = value,
    eligible = eligible,
    left_space = left_space,
    right_space = right_space
  ), class = "effect_pair_coupling")
}

.validate_pair_coupling <- function(x) {
  expected <- c("schema_version", "kind", "value", "eligible", "left_space",
    "right_space")
  if (!.sealed_fields(x, "effect_pair_coupling", expected) ||
      !identical(x$schema_version, 1L) || !x$kind %in% c("match", "control") ||
      !.is_finite_matrix(x$value) || any(x$value < 0) ||
      !is.matrix(x$eligible) || !is.logical(x$eligible) || anyNA(x$eligible)) {
    .input_error("Invalid pair coupling.")
  }
  left_space <- .validate_effect_space(x$left_space)
  right_space <- .validate_effect_space(x$right_space)
  shape <- c(length(left_space$coordinates), length(right_space$coordinates))
  if (!identical(dim(x$value), shape) || !identical(dim(x$eligible), shape) ||
      any(x$value[!x$eligible] != 0)) {
    .input_error("Pair coupling shape or eligibility is inconsistent.")
  }
  x
}

#' Mark matched encoding-retrieval effect pairs
#'
#' Duplicate rows are retained as multiplicity, so repeated retrievals are
#' explicit rather than silently deduplicated.
#'
#' @param left,right Matched coordinate identifiers of equal positive length.
#' @param left_space,right_space Ordered effect-space identities.
#' @param eligible Optional restricted eligible-pair set.
#' @return An `effect_pair_coupling` with `kind = "match"`, a nonnegative
#'   `$value` matrix counting the multiplicity of each matched cell, the
#'   logical `$eligible` mask, and the bound `$left_space`/`$right_space`.
#' @seealso [control_coupling()] for the complementary cells,
#'   [coupling_contrast()] and [match_control()] to turn the pair into a
#'   query.
#' @family coupling and connectivity views
#' @examples
#' # Three studied items and their matched retrieval probes.
#' encoding <- effect_space(
#'   c("item1", "item2", "item3"), basis_id = "demo:encoding:v1"
#' )
#' retrieval <- effect_space(
#'   c("probe1", "probe2", "probe3"), basis_id = "demo:retrieval:v1"
#' )
#' matches <- match_coupling(
#'   c("item1", "item2", "item3"), c("probe1", "probe2", "probe3"),
#'   encoding, retrieval
#' )
#' matches$value
#'
#' # Repeated retrievals of one item are multiplicity, not duplicates, so the
#' # cell count rises rather than being silently collapsed.
#' match_coupling(
#'   c("item1", "item1", "item2"), c("probe1", "probe1", "probe2"),
#'   encoding, retrieval
#' )$value
#'
#' # Identifiers must belong to their declared axis.
#' refused <- try(
#'   match_coupling("item1", "item2", encoding, retrieval), silent = TRUE
#' )
#' conditionMessage(attr(refused, "condition"))
#' @export
match_coupling <- function(left, right, left_space, right_space,
                           eligible = NULL) {
  if (length(left) != length(right) || length(left) < 1L) {
    .input_error(
      "Matched `left` and `right` identifiers must have equal positive length."
    )
  }
  left_space <- .as_effect_space(left_space)
  right_space <- .as_effect_space(right_space)
  eligible <- .pair_eligibility(eligible, left_space, right_space)
  left_index <- .pair_axis_index(left, left_space, "left")
  right_index <- .pair_axis_index(right, right_space, "right")
  if (any(!eligible[cbind(left_index, right_index)])) {
    .input_error("Every matched pair must be eligible.")
  }
  value <- matrix(0, nrow(eligible), ncol(eligible), dimnames = dimnames(eligible))
  for (row in seq_along(left_index)) {
    value[left_index[[row]], right_index[[row]]] <-
      value[left_index[[row]], right_index[[row]]] + 1
  }
  .new_pair_coupling("match", value, eligible, left_space, right_space)
}

#' Mark eligible control pairs
#'
#' `control_coupling()` names the comparison set for a matched-pair analysis:
#' the eligible cells that are not matches. Pair it with the matches through
#' [coupling_contrast()] to obtain the matched-versus-control query.
#'
#' @param matches A `match_coupling()` value.
#' @param include_matches Whether matched cells are also controls. The default
#'   excludes them.
#' @return An `effect_pair_coupling` with `kind = "control"`, a 0/1 `$value`
#'   indicator over the same axes and `$eligible` set as `matches`.
#' @seealso [match_coupling()] and [coupling_contrast()].
#' @family coupling and connectivity views
#' @examples
#' encoding <- effect_space(
#'   c("item1", "item2", "item3"), basis_id = "demo:encoding:v1"
#' )
#' retrieval <- effect_space(
#'   c("probe1", "probe2", "probe3"), basis_id = "demo:retrieval:v1"
#' )
#' matches <- match_coupling(
#'   c("item1", "item2", "item3"), c("probe1", "probe2", "probe3"),
#'   encoding, retrieval
#' )
#'
#' # The default controls are the six mismatched cells.
#' controls <- control_coupling(matches)
#' controls$value
#'
#' # Including the matches gives every eligible cell, which is the marginal
#' # baseline rather than a contrast partner.
#' control_coupling(matches, include_matches = TRUE)$value
#'
#' # With no eligible non-matched cell left, there is nothing to contrast.
#' saturated <- match_coupling(
#'   rep(c("item1", "item2", "item3"), each = 3),
#'   rep(c("probe1", "probe2", "probe3"), 3), encoding, retrieval
#' )
#' refused <- try(control_coupling(saturated), silent = TRUE)
#' conditionMessage(attr(refused, "condition"))
#' @export
control_coupling <- function(matches, include_matches = FALSE) {
  matches <- .validate_pair_coupling(matches)
  if (matches$kind != "match" || !.is_flag(include_matches)) {
    .input_error("Controls require a match coupling and one inclusion flag.")
  }
  selected <- matches$eligible & (include_matches | matches$value == 0)
  if (!any(selected)) .input_error("No eligible control pairs remain.")
  value <- matrix(as.numeric(selected), nrow(selected), ncol(selected),
    dimnames = dimnames(selected))
  .new_pair_coupling("control", value, matches$eligible,
    matches$left_space, matches$right_space)
}

.pair_balance_diagnostics <- function(H, tolerance = 1e-12) {
  row <- rowSums(H)
  column <- colSums(H)
  list(
    row_marginals = row,
    column_marginals = column,
    maximum_absolute_row_marginal = max(abs(row)),
    maximum_absolute_column_marginal = max(abs(column)),
    zero_row_marginals = all(abs(row) <= tolerance),
    zero_column_marginals = all(abs(column) <= tolerance),
    additive_baseline_invariant = all(abs(row) <= tolerance) &&
      all(abs(column) <= tolerance),
    tolerance = tolerance
  )
}

#' Contrast matched and control pair couplings
#'
#' `coupling_contrast()` subtracts the control coupling from the matched one
#' to give the matched-versus-control readout as a single fixed pair query.
#' Normalizing first makes the two sides comparable when they contain
#' different numbers of cells.
#'
#' @param matches,controls Compatible match and control couplings.
#' @param normalize Normalize each coupling to unit total mass before
#'   subtraction.
#' @return An axis-bound [pair_query()] whose `$operator` is the contrast, and
#'   whose `$metadata$balance` reports the row and column marginals and
#'   whether the operator is `additive_baseline_invariant`.
#' @seealso [match_coupling()], [control_coupling()], and [match_control()],
#'   which instead compiles the same comparison as a regression coefficient
#'   with item nuisance effects.
#' @family coupling and connectivity views
#' @examples
#' encoding <- effect_space(
#'   c("item1", "item2", "item3"), basis_id = "demo:encoding:v1"
#' )
#' retrieval <- effect_space(
#'   c("probe1", "probe2", "probe3"), basis_id = "demo:retrieval:v1"
#' )
#' matches <- match_coupling(
#'   c("item1", "item2", "item3"), c("probe1", "probe2", "probe3"),
#'   encoding, retrieval
#' )
#' contrast <- coupling_contrast(matches, control_coupling(matches))
#' round(as.matrix(contrast$operator), 3)
#'
#' # Both marginals are zero here, so the readout is unchanged by adding a
#' # constant to any item or probe effect.
#' contrast$metadata$balance$additive_baseline_invariant
#'
#' # That balance is an observed property of this operator, not a promise
#' # about anything applied downstream.
#' contrast$metadata$claim
#' @export
coupling_contrast <- function(matches, controls, normalize = TRUE) {
  matches <- .validate_pair_coupling(matches)
  controls <- .validate_pair_coupling(controls)
  if (matches$kind != "match" || controls$kind != "control" ||
      !identical(matches$left_space, controls$left_space) ||
      !identical(matches$right_space, controls$right_space) ||
      !identical(matches$eligible, controls$eligible) || !.is_flag(normalize)) {
    .contract_error(
      "Pair contrast couplings or normalization are incompatible."
    )
  }
  matched <- matches$value
  control <- controls$value
  if (normalize) {
    if (sum(matched) <= 0 || sum(control) <= 0) {
      .input_error(
        "Normalized pair contrasts require positive mass on both sides."
      )
    }
    matched <- matched / sum(matched)
    control <- control / sum(control)
  }
  H <- matched - control
  pair_query(H, matches$left_space, matches$right_space, metadata = list(
    constructor = "coupling_contrast",
    balance = .pair_balance_diagnostics(H),
    claim = "observed_operator_balance_only"
  ))
}

.pair_design_matrix <- function(design, left_space, right_space,
                                encoding_nuisance, retrieval_nuisance) {
  if (!is.data.frame(design) || nrow(design) < 1L ||
      !all(c("left", "right") %in% names(design))) {
    .input_error(
      "`design` must be a nonempty data frame with `left` and `right`."
    )
  }
  left_index <- .pair_axis_index(as.character(design$left), left_space, "left")
  right_index <- .pair_axis_index(as.character(design$right), right_space, "right")
  reserved <- c("left", "right", "weight")
  predictors <- setdiff(names(design), reserved)
  if (length(predictors) < 1L ||
      any(!vapply(design[predictors], is.numeric, logical(1))) ||
      any(!vapply(design[predictors], function(value) all(is.finite(value)),
        logical(1)))) {
    .input_error("Pair designs require at least one finite numeric predictor.")
  }
  X <- cbind(`(Intercept)` = 1, as.matrix(design[predictors]))
  if (encoding_nuisance) {
    if (length(unique(left_index)) < 2L) {
      .input_error(
        "Encoding nuisance effects are infeasible with fewer than two items."
      )
    }
    nuisance <- stats::model.matrix(~ factor(left_index))[, -1L, drop = FALSE]
    colnames(nuisance) <- paste0("encoding:", colnames(nuisance))
    X <- cbind(X, nuisance)
  }
  if (retrieval_nuisance) {
    if (length(unique(right_index)) < 2L) {
      .input_error(
        "Retrieval nuisance effects are infeasible with fewer than two items."
      )
    }
    nuisance <- stats::model.matrix(~ factor(right_index))[, -1L, drop = FALSE]
    colnames(nuisance) <- paste0("retrieval:", colnames(nuisance))
    X <- cbind(X, nuisance)
  }
  list(X = X, left_index = left_index, right_index = right_index,
    predictors = predictors)
}

#' Compile a weighted pair-space linear-model coefficient to an effect query
#'
#' The returned matrix is the exact linear map from eligible pair values to a
#' requested weighted least-squares coefficient. Encoding and retrieval
#' nuisance effects are ordinary design columns, not a new analysis class.
#'
#' @param design A data frame with `left`, `right`, and one or more finite
#'   numeric predictor columns. Duplicate pair rows are allowed.
#' @param coefficient One coefficient name, or a named numeric contrast over
#'   the compiled design columns.
#' @param left_space,right_space Ordered effect-space identities.
#' @param weights Optional positive row weights, or the name of a design weight
#'   column. A `weight` column is used automatically when present.
#' @param encoding_nuisance,retrieval_nuisance Include fixed-effect nuisance
#'   columns for the respective item axis.
#' @param sparse Return `H` as a sparse `Matrix` object.
#' @return An axis-bound [pair_query()] whose `$operator` maps eligible pair
#'   values to the requested coefficient, with
#'   `$metadata$coefficient` (the contrast over compiled columns) and
#'   `$metadata$diagnostics` reporting `rank`, `columns`, `observations`,
#'   `unique_pairs`, and operator `balance`.
#' @seealso [match_control()] for the matched-versus-control special case,
#'   and [pair_query()] for a hand-written operator.
#' @family coupling and connectivity views
#' @examples
#' # A pair-space regression: how does encoding-retrieval similarity change
#' # with study-test lag, adjusting for a match indicator?
#' encoding <- effect_space(
#'   c("item1", "item2", "item3"), basis_id = "demo:encoding:v1"
#' )
#' retrieval <- effect_space(
#'   c("probe1", "probe2", "probe3"), basis_id = "demo:retrieval:v1"
#' )
#' design <- expand.grid(
#'   left = encoding$coordinates, right = retrieval$coordinates,
#'   stringsAsFactors = FALSE
#' )
#' design$lag <- abs(
#'   match(design$left, encoding$coordinates) -
#'     match(design$right, retrieval$coordinates)
#' )
#' design$match <- as.numeric(design$lag == 0)
#'
#' # The result is the exact linear map from pair values to the `lag`
#' # coefficient, compiled once and reusable as a fixed query.
#' query <- pair_lm_query(design, "lag", encoding, retrieval)
#' round(as.matrix(query$operator), 3)
#' query$metadata$diagnostics[c("rank", "columns", "observations")]
#'
#' # A collinear predictor makes the coefficient undefined, and the design is
#' # rejected before any geometry is read.
#' design$lag_copy <- design$lag
#' refused <- try(
#'   pair_lm_query(design, "lag", encoding, retrieval), silent = TRUE
#' )
#' conditionMessage(attr(refused, "condition"))
#' @export
pair_lm_query <- function(design, coefficient, left_space, right_space,
                          weights = NULL, encoding_nuisance = FALSE,
                          retrieval_nuisance = FALSE, sparse = FALSE) {
  left_space <- .as_effect_space(left_space)
  right_space <- .as_effect_space(right_space)
  flags <- c(encoding_nuisance, retrieval_nuisance)
  if (!is.logical(flags) || length(flags) != 2L || anyNA(flags)) {
    .input_error("Nuisance flags must be TRUE or FALSE.")
  }
  .check_flag(sparse, "sparse")
  compiled <- .pair_design_matrix(
    design, left_space, right_space,
    encoding_nuisance, retrieval_nuisance
  )
  X <- compiled$X
  if (is.character(weights) && length(weights) == 1L && !is.na(weights)) {
    if (!weights %in% names(design)) .input_error("Unknown design weight column.")
    weights <- design[[weights]]
  } else if (is.null(weights) && "weight" %in% names(design)) {
    weights <- design$weight
  } else if (is.null(weights)) {
    weights <- rep(1, nrow(design))
  }
  if (!.is_finite_numeric(weights) || length(weights) != nrow(design) ||
      anyNA(weights) || any(weights <= 0)) {
    .input_error("Pair-design weights must be positive and finite.")
  }
  if (is.character(coefficient) && length(coefficient) == 1L &&
      !is.na(coefficient)) {
    if (!coefficient %in% colnames(X)) .input_error("Unknown pair-design coefficient.")
    contrast <- stats::setNames(rep(0, ncol(X)), colnames(X))
    contrast[[coefficient]] <- 1
  } else {
    if (!.is_finite_numeric(coefficient) || is.null(names(coefficient)) ||
        anyNA(coefficient) || anyDuplicated(names(coefficient)) ||
        !setequal(names(coefficient), colnames(X))) {
      .input_error(
        "A coefficient contrast must name every compiled design column."
      )
    }
    contrast <- coefficient[colnames(X)]
  }
  root_weight <- sqrt(weights)
  weighted_X <- X * root_weight
  qr_X <- qr(weighted_X, tol = 1e-10)
  if (qr_X$rank != ncol(X)) {
    .input_error(
      "Pair-space design is rank deficient; revise predictors or nuisances."
    )
  }
  information <- crossprod(X, weights * X)
  weighted_transpose <- sweep(t(X), 2L, weights, `*`)
  influence <- drop(contrast %*% solve(information, weighted_transpose))
  H <- matrix(0, length(left_space$coordinates),
    length(right_space$coordinates),
    dimnames = list(left_space$coordinates, right_space$coordinates))
  for (row in seq_len(nrow(design))) {
    H[compiled$left_index[[row]], compiled$right_index[[row]]] <-
      H[compiled$left_index[[row]], compiled$right_index[[row]]] +
      influence[[row]]
  }
  diagnostics <- list(
    rank = qr_X$rank,
    columns = colnames(X),
    observations = nrow(X),
    unique_pairs = sum(H != 0),
    encoding_nuisance = encoding_nuisance,
    retrieval_nuisance = retrieval_nuisance,
    balance = .pair_balance_diagnostics(H)
  )
  if (sparse) H <- Matrix::Matrix(H, sparse = TRUE)
  pair_query(H, left_space, right_space, metadata = list(
    constructor = "pair_lm_query",
    coefficient = contrast,
    diagnostics = diagnostics,
    claim = "balance_reported_after_final_weighting"
  ))
}

#' Compile a matched-versus-control pair-space coefficient
#'
#' This convenience uses the eligible pair set carried by `matches`, includes
#' both item nuisance families by default, and reports the balance of the final
#' operator. It does not claim invariance after any later external weighting.
#'
#' @param matches A `match_coupling()` value.
#' @param weights Optional positive eligible-pair weights.
#' @param encoding_nuisance,retrieval_nuisance Nuisance-effect flags.
#' @return An axis-bound [pair_query()] for the matched-versus-control
#'   coefficient, with `$metadata$constructor` set to `"match_control"`,
#'   `$metadata$diagnostics` (design `rank`, `columns`, operator `balance`),
#'   and a `$metadata$claim` naming exactly what the reported balance covers.
#' @seealso [match_coupling()] for the input, [coupling_contrast()] for the
#'   plain difference without nuisance effects, and [pair_lm_query()] for the
#'   general designed coefficient.
#' @family coupling and connectivity views
#' @examples
#' encoding <- effect_space(
#'   c("item1", "item2", "item3"), basis_id = "demo:encoding:v1"
#' )
#' retrieval <- effect_space(
#'   c("probe1", "probe2", "probe3"), basis_id = "demo:retrieval:v1"
#' )
#' matches <- match_coupling(
#'   c("item1", "item2", "item3"), c("probe1", "probe2", "probe3"),
#'   encoding, retrieval
#' )
#'
#' # Both item nuisance families are included by default, so item-specific
#' # offsets cannot masquerade as a matching effect.
#' query <- match_control(matches)
#' query$metadata$diagnostics$columns
#' round(as.matrix(query$operator), 3)
#'
#' # The result records what its balance diagnostics actually claim.
#' query$metadata$claim
#'
#' # Dropping both nuisance families gives the simpler unadjusted contrast.
#' unadjusted <- match_control(
#'   matches, encoding_nuisance = FALSE, retrieval_nuisance = FALSE
#' )
#' unadjusted$metadata$diagnostics$columns
#' @export
match_control <- function(matches, weights = NULL,
                          encoding_nuisance = TRUE,
                          retrieval_nuisance = TRUE) {
  matches <- .validate_pair_coupling(matches)
  cells <- which(matches$eligible, arr.ind = TRUE)
  design <- data.frame(
    left = matches$left_space$coordinates[cells[, 1L]],
    right = matches$right_space$coordinates[cells[, 2L]],
    match = as.numeric(matches$value[cells] > 0),
    stringsAsFactors = FALSE
  )
  if (!is.null(weights)) design$weight <- weights
  query <- pair_lm_query(
    design, "match", matches$left_space, matches$right_space,
    encoding_nuisance = encoding_nuisance,
    retrieval_nuisance = retrieval_nuisance
  )
  query$metadata$constructor <- "match_control"
  query$metadata$claim <- "observed_operator_balance_not_downstream_invariance"
  query
}
