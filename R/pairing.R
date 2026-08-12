# Partition pairings and their signed marginals -----------------------------

#' Construct a partition pairing
#'
#' @param left,right Equal-length vectors naming partition endpoints.
#' @param weight Optional finite nonnegative edge weights. They are normalized
#'   to sum to one.
#' @param directed Whether endpoint roles have scientific meaning.
#' @param self_pairs Whether diagonal self-products are forbidden or explicitly
#'   admitted as noise-biased estimates.
#' @param independence Whether distinct endpoint estimates are declared
#'   independent. Self-products can only be marked `not_independent`.
#' @return An edge table with explicit endpoint semantics.
#' @export
pairing <- function(left, right, weight = NULL, directed = FALSE,
                    self_pairs = c("forbid", "allow_biased"),
                    independence = c("independent", "not_independent")) {
  self_pairs <- match.arg(self_pairs)
  independence <- match.arg(independence)
  if (length(left) != length(right) || length(left) < 1L) {
    stop("`left` and `right` must have the same positive length.", call. = FALSE)
  }
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("`directed` must be TRUE or FALSE.", call. = FALSE)
  }

  left <- as.character(left)
  right <- as.character(right)
  if (anyNA(left) || anyNA(right) || any(left == "") || any(right == "")) {
    stop("Pairing endpoints must be non-missing, nonempty identifiers.", call. = FALSE)
  }

  if (is.null(weight)) {
    weight <- rep(1, length(left))
  }
  if (!is.numeric(weight) || length(weight) != length(left) ||
      any(!is.finite(weight)) || any(weight < 0) || max(weight) <= 0) {
    stop("`weight` must be finite, nonnegative, and have positive total mass.",
      call. = FALSE)
  }

  is_self <- left == right
  if (any(is_self) && any(!is_self)) {
    stop("A pairing cannot mix self-products with cross-partition products.",
      call. = FALSE)
  }
  if (any(is_self) && self_pairs != "allow_biased") {
    stop("Self-products require `self_pairs = \"allow_biased\"`.", call. = FALSE)
  }
  if (any(is_self) && independence != "not_independent") {
    stop("Self-products must declare `independence = \"not_independent\"`.",
      call. = FALSE)
  }

  key <- if (directed) {
    paste(left, right, sep = "\r")
  } else {
    paste(pmin(left, right), pmax(left, right), sep = "\r")
  }
  if (anyDuplicated(key)) {
    stop("The pairing contains duplicate edges.", call. = FALSE)
  }

  scaled_weight <- weight / max(weight)
  normalized_weight <- scaled_weight / sum(scaled_weight)
  if (any(!is.finite(normalized_weight)) || sum(normalized_weight) <= 0 ||
      abs(sum(normalized_weight) - 1) > 1e-12) {
    stop("Normalized pairing weights must be finite and have positive unit mass.",
      call. = FALSE)
  }

  ans <- data.frame(
    left = left,
    right = right,
    weight = normalized_weight,
    stringsAsFactors = FALSE
  )
  attr(ans, "directed") <- directed
  attr(ans, "self_pairs") <- self_pairs
  attr(ans, "independence") <- independence
  attr(ans, "estimate") <- if (any(is_self)) {
    "self_product_biased"
  } else if (independence == "independent") {
    "cross_generalized"
  } else {
    "nonindependent_cross_product"
  }
  class(ans) <- c("effect_pairing", "data.frame")
  .validate_pairing(ans)
  ans
}

#' Pair every distinct partition once
#'
#' @param partitions Partition identifiers or an `effect_relation`.
#' @return An undirected pairing containing one row per unordered pair.
#' @export
cross_partitions <- function(partitions) {
  if (inherits(partitions, "effect_relation")) {
    .validate_relation(partitions)
    partitions <- partitions$partitions
  }
  partitions <- unique(as.character(partitions))
  if (length(partitions) < 2L || anyNA(partitions) || any(partitions == "")) {
    stop("At least two distinct, non-missing partitions are required.", call. = FALSE)
  }
  edges <- utils::combn(partitions, 2L)
  pairing(edges[1L, ], edges[2L, ], directed = FALSE)
}

#' Compute pairing-appropriate signed relation marginals
#'
#' @param local A measurement-by-effect-by-partition numeric array containing
#'   spatially weighted relation sums.
#' @param over An `effect_pairing`.
#' @param mass Positive frame mass for each measurement.
#' @return For undirected pairings, a list containing only `endpoint`; for
#'   directed pairings, a list containing `left` and `right`.
#' @export
pairing_marginals <- function(local, over, mass = 1) {
  if (!is.array(local) || length(dim(local)) != 3L || !is.numeric(local) ||
      any(!is.finite(local))) {
    stop("`local` must be a finite numeric measurement-by-effect-by-partition array.",
      call. = FALSE)
  }
  if (!inherits(over, "effect_pairing")) {
    stop("`over` must be an effect pairing.", call. = FALSE)
  }
  .validate_pairing(over)

  partition_ids <- dimnames(local)[[3L]]
  if (is.null(partition_ids) || anyNA(partition_ids) || any(partition_ids == "")) {
    stop("The partition dimension of `local` must have unique names.", call. = FALSE)
  }
  if (anyDuplicated(partition_ids)) {
    stop("The partition dimension of `local` must have unique names.", call. = FALSE)
  }
  left_index <- match(over$left, partition_ids)
  right_index <- match(over$right, partition_ids)
  if (anyNA(left_index) || anyNA(right_index)) {
    stop("Every pairing endpoint must identify a partition in `local`.", call. = FALSE)
  }

  m <- dim(local)[1L]
  q <- dim(local)[2L]
  if (!is.numeric(mass) || length(mass) == 0L ||
      !(length(mass) %in% c(1L, m)) || any(!is.finite(mass)) || any(mass <= 0)) {
    stop("`mass` must be one positive finite value or one per measurement.",
      call. = FALSE)
  }
  mass <- rep_len(mass, m)

  normalized_partition <- function(k) {
    matrix(local[, , k, drop = FALSE], nrow = m, ncol = q) / mass
  }
  accumulate_role <- function(indices) {
    out <- matrix(0, nrow = m, ncol = q)
    for (edge in seq_len(nrow(over))) {
      out <- out + over$weight[[edge]] * normalized_partition(indices[[edge]])
    }
    dimnames(out) <- dimnames(local)[1:2]
    out
  }

  if (isTRUE(attr(over, "directed"))) {
    ans <- list(
      left = accumulate_role(left_index),
      right = accumulate_role(right_index)
    )
    attr(ans, "semantics") <- "directed_roles"
  } else {
    endpoint <- matrix(0, nrow = m, ncol = q)
    for (edge in seq_len(nrow(over))) {
      endpoint <- endpoint + over$weight[[edge]] * 0.5 *
        (normalized_partition(left_index[[edge]]) +
          normalized_partition(right_index[[edge]]))
    }
    dimnames(endpoint) <- dimnames(local)[1:2]
    ans <- list(endpoint = endpoint)
    attr(ans, "semantics") <- "undirected_endpoint"
  }

  class(ans) <- c("effect_marginals", "list")
  ans
}

.validate_pairing <- function(over) {
  if (!inherits(over, "effect_pairing") || !is.data.frame(over) ||
      !identical(names(over), c("left", "right", "weight")) || nrow(over) < 1L) {
    stop("Pairing objects must be nonempty left/right/weight edge tables.",
      call. = FALSE)
  }
  directed <- attr(over, "directed", exact = TRUE)
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("Pairing directedness must be one logical value.", call. = FALSE)
  }
  self_policy <- attr(over, "self_pairs", exact = TRUE)
  independence <- attr(over, "independence", exact = TRUE)
  estimate <- attr(over, "estimate", exact = TRUE)
  if (!is.character(self_policy) || length(self_policy) != 1L ||
      !self_policy %in% c("forbid", "allow_biased")) {
    stop("Pairing self-product policy is missing or invalid.", call. = FALSE)
  }
  if (!is.character(independence) || length(independence) != 1L ||
      !independence %in% c("independent", "not_independent")) {
    stop("Pairing independence semantics are missing or invalid.", call. = FALSE)
  }
  if (!is.character(over$left) || !is.character(over$right) ||
      anyNA(over$left) || anyNA(over$right) ||
      any(over$left == "") || any(over$right == "")) {
    stop("Pairing endpoints must be non-missing, nonempty identifiers.", call. = FALSE)
  }
  if (!is.numeric(over$weight) || length(over$weight) != nrow(over) ||
      any(!is.finite(over$weight)) || any(over$weight < 0) ||
      sum(over$weight) <= 0 || abs(sum(over$weight) - 1) > 1e-12) {
    stop("Pairing weights must be finite, nonnegative, and normalized to unit mass.",
      call. = FALSE)
  }

  is_self <- over$left == over$right
  if (any(is_self) && any(!is_self)) {
    stop("A pairing cannot mix self-products with cross-partition products.",
      call. = FALSE)
  }
  if (any(is_self) && self_policy != "allow_biased") {
    stop("Self-products require explicit biased-estimate admission.", call. = FALSE)
  }
  if (any(is_self) && independence != "not_independent") {
    stop("Self-products cannot be declared independent.", call. = FALSE)
  }
  expected_estimate <- if (any(is_self)) {
    "self_product_biased"
  } else if (independence == "independent") {
    "cross_generalized"
  } else {
    "nonindependent_cross_product"
  }
  if (!identical(estimate, expected_estimate)) {
    stop("Pairing estimate semantics are inconsistent with its edges.", call. = FALSE)
  }

  key <- if (directed) {
    paste(over$left, over$right, sep = "\r")
  } else {
    paste(pmin(over$left, over$right), pmax(over$left, over$right), sep = "\r")
  }
  if (anyDuplicated(key)) {
    stop("The pairing contains duplicate edges.", call. = FALSE)
  }
  invisible(over)
}
