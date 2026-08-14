# Independent numerical oracles for the normative effect-form laws. These
# helpers deliberately do not call production effectagram functions.

oracle_effect_form <- function(left, right, weight = NULL) {
  if (!is.matrix(left) || !is.matrix(right) || ncol(left) != ncol(right)) {
    stop("Oracle relations must be matrices over one common measurement axis.")
  }
  if (is.null(weight)) weight <- rep(1, ncol(left))
  left %*% diag(weight, nrow = length(weight)) %*% t(right)
}

oracle_query <- function(form, query) sum(form * query)

oracle_bind_query <- function(query, left_space, right_space) {
  if (!is.matrix(query) || !is.numeric(query) || any(!is.finite(query))) {
    stop("query must be one finite matrix")
  }
  if (!is.list(left_space) || !is.list(right_space) ||
      !identical(names(left_space), c("id", "coordinates")) ||
      !identical(names(right_space), c("id", "coordinates")) ||
      !is.character(left_space$id) || length(left_space$id) != 1L ||
      !is.character(right_space$id) || length(right_space$id) != 1L ||
      !identical(dim(query), c(length(left_space$coordinates),
        length(right_space$coordinates)))) {
    stop("query axes do not match their bound spaces")
  }
  list(query = query, left_space = left_space, right_space = right_space)
}

oracle_apply_bound_query <- function(form, bound_query, left_space, right_space) {
  if (!identical(bound_query$left_space, left_space) ||
      !identical(bound_query$right_space, right_space)) {
    stop("query and form axis identities differ")
  }
  oracle_query(form, bound_query$query)
}

oracle_capabilities <- function(self_form, symmetric, guaranteed_psd) {
  flags <- c(self_form, symmetric, guaranteed_psd)
  if (!is.logical(flags) || length(flags) != 3L || anyNA(flags)) {
    stop("capabilities must be logical guarantees")
  }
  if (symmetric && !self_form) stop("symmetric implies self_form")
  if (guaranteed_psd && !symmetric) stop("guaranteed_psd implies symmetric")
  stats::setNames(flags, c("self_form", "symmetric", "guaranteed_psd"))
}

oracle_result_contract <- function(capabilities, codec, completeness) {
  codec <- match.arg(codec, c("rectangular", "symmetric_packed"))
  completeness <- match.arg(completeness, c("complete_form", "query_only"))
  if (codec == "symmetric_packed" && !isTRUE(capabilities[["symmetric"]])) {
    stop("packed codec requires a symmetric guarantee")
  }
  list(capabilities = capabilities, codec = codec, completeness = completeness)
}

oracle_vec <- function(form) as.vector(form)

oracle_svec <- function(form) {
  stopifnot(is.matrix(form), nrow(form) == ncol(form))
  out <- numeric(nrow(form) * (nrow(form) + 1L) / 2L)
  coordinate <- 0L
  for (column in seq_len(ncol(form))) {
    for (row in column:nrow(form)) {
      coordinate <- coordinate + 1L
      out[[coordinate]] <- form[row, column] *
        if (row == column) 1 else sqrt(2)
    }
  }
  out
}

oracle_centered_statistics <- function(left, right, weight) {
  mass <- sum(weight)
  total <- oracle_effect_form(left, right, weight)
  left_first <- as.vector(left %*% weight)
  right_first <- as.vector(right %*% weight)
  coherent <- outer(left_first, right_first) / mass
  list(
    mass = mass,
    total = total,
    left_first = left_first,
    right_first = right_first,
    coherent = coherent,
    centered = total - coherent,
    covariance = (total - coherent) / mass
  )
}

oracle_explicit_center <- function(left, right, weight) {
  mass <- sum(weight)
  left_mean <- as.vector(left %*% weight) / mass
  right_mean <- as.vector(right %*% weight) / mass
  centered_left <- sweep(left, 1L, left_mean, "-")
  centered_right <- sweep(right, 1L, right_mean, "-")
  oracle_effect_form(centered_left, centered_right, weight)
}

oracle_weighted_correlation <- function(left, right, weight,
                                        zero_variance = "error") {
  zero_variance <- match.arg(zero_variance, c("error", "zero"))
  cross <- oracle_centered_statistics(left, right, weight)$centered
  left_ss <- diag(oracle_centered_statistics(left, left, weight)$centered)
  right_ss <- diag(oracle_centered_statistics(right, right, weight)$centered)
  zero_left <- left_ss == 0
  zero_right <- right_ss == 0
  if (zero_variance == "error" && (any(zero_left) || any(zero_right))) {
    stop("zero weighted variance")
  }
  denominator <- sqrt(outer(left_ss, right_ss))
  out <- cross / denominator
  out[outer(zero_left, zero_right, `|`)] <- 0
  out
}

oracle_fisher <- function(value, boundary = "error", delta = NULL) {
  boundary <- match.arg(boundary, c("error", "clip"))
  if (boundary == "error") {
    if (any(abs(value) >= 1)) stop("Fisher boundary")
    return(atanh(value))
  }
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) ||
      delta <= 0 || delta >= 1) {
    stop("invalid Fisher clipping delta")
  }
  atanh(pmin(1 - delta, pmax(-1 + delta, value)))
}

oracle_bridge_form <- function(left, right, left_leg, right_leg) {
  measured_left <- left %*% t(left_leg)
  measured_right <- right %*% t(right_leg)
  oracle_effect_form(measured_left, measured_right)
}

oracle_direct_bridge_form <- function(left, right, left_leg, right_leg) {
  bridge <- t(left_leg) %*% right_leg
  left %*% bridge %*% t(right)
}

oracle_require_bridge <- function(left_neural_id, right_neural_id, bridge = NULL) {
  if (identical(left_neural_id, right_neural_id)) return(invisible(TRUE))
  if (is.null(bridge)) stop("distinct neural spaces require a bridge")
  invisible(TRUE)
}
