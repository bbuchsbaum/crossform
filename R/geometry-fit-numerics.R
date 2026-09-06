# Regularized prediction of a signed form (layer 5) -------------------------
#
# predictive-geometry-v1: the fixed model pool supplies U D U'. First
# restrict S to its positive support, then shift by lambda D^-1. The latent
# primitive supplies one symmetric eigendecomposition. Rank selection only
# truncates this frozen spectrum; it neither rotates nor reads more data.
# The original signed source and the shifted source have separate records.

.geometry_fit_symmetric <- function(S) {
  if (!.is_finite_matrix(S) || nrow(S) < 1L || nrow(S) != ncol(S)) {
    .input_error("`S` must be a nonempty finite square symmetric signed form.")
  }
  if (max(abs(S - t(S))) > 1e-12 * max(abs(S))) {
    .input_error("`S` must be symmetric; an asymmetric input is not a declared form.")
  }
  S / 2 + t(S) / 2
}

.geometry_fit_path <- function(S, pool, penalty) {
  S <- .geometry_fit_symmetric(S)
  .check_number(penalty, "penalty", nonnegative = TRUE)
  U <- pool$vectors
  d <- pool$values
  if (!.is_finite_matrix(U) || nrow(U) != nrow(S) ||
      !.is_finite_numeric(d) || length(d) != ncol(U) || !length(d) ||
      any(d <= 0) || max(abs(crossprod(U) - diag(length(d)))) > 1e-8) {
    .input_error("The model pool must supply orthonormal positive support matching `S`.")
  }
  restricted <- crossprod(U, S %*% U)
  shifted <- restricted / 2 + t(restricted) / 2 - diag(penalty / d, length(d))
  if (!.is_finite_matrix(shifted)) {
    .input_error("The model penalty shift overflows; rescale the form, kernel or penalty.")
  }
  spectral <- .latent_rank_psd_form(shifted)
  signed <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(c(spectral$values, signed)))) {
    .input_error("The form spectrum overflows; rescale the input.")
  }
  list(source = S, signed_spectrum = signed,
    source_negative_mass = -sum(signed[signed < 0]),
    shifted = shifted, shifted_spectrum = spectral$values,
    shifted_vectors = spectral$vectors,
    shifted_negative_mass = spectral$clipped_negative_mass,
    pool_vectors = U, pool_values = d, pool_signature = pool$signature,
    support_dimension = length(d), union_dimension = nrow(S),
    penalty = unname(penalty), tolerance = pool$tolerance)
}

# A positive tied group has invariant aggregate evidence. Its individual
# basis vectors do not. Group membership follows the declared spectral gap
# tolerance; a rank cut through such a group is refused before this step.
.geometry_fit_mode_groups <- function(amplitudes, spectral_scale, tolerance) {
  if (!length(amplitudes)) return(integer())
  as.integer(cumsum(c(TRUE, abs(diff(amplitudes)) > tolerance * spectral_scale)))
}

.geometry_fit_rank <- function(path, rank) {
  rank <- .check_count(rank, "rank", min = 0L, max = .Machine$integer.max)
  budget <- min(rank, path$support_dimension)
  .model_basis_rank_cut(path$shifted_spectrum, budget, path$tolerance,
    "the response form")
  projection <- .latent_psd_projection(matrix(path$shifted_spectrum, nrow = 1L),
    budget)
  amplitudes <- as.numeric(projection$spectrum[1L, ])
  active <- which(amplitudes > 0)
  amplitudes <- amplitudes[active]
  modes <- path$pool_vectors %*% path$shifted_vectors[, active, drop = FALSE]
  factor <- sweep(modes, 2L, sqrt(amplitudes), "*")
  A <- tcrossprod(factor)
  dimnames(A) <- dimnames(path$source)
  norm_sq <- sum(amplitudes^2)
  penalty_trace <- sum(rowSums(crossprod(path$pool_vectors, factor)^2) / path$pool_values)
  residual_sq <- sum((path$source - A)^2)
  objective <- 0.5 * residual_sq + path$penalty * penalty_trace
  if (any(!is.finite(c(A, norm_sq, penalty_trace, residual_sq, objective)))) {
    .input_error("The fitted form or loss overflows; rescale the input before fitting.")
  }
  list(form = A, factor = factor, modes = modes, amplitudes = amplitudes,
    groups = .geometry_fit_mode_groups(amplitudes,
      max(abs(path$shifted_spectrum)), path$tolerance),
    rank_requested = rank, rank_budget = as.integer(budget),
    rank_clamped = rank > budget, rank_effective = as.integer(length(active)),
    prediction_norm_sq = norm_sq, penalty_trace = penalty_trace,
    compressed_residual_sq = residual_sq, objective = objective,
    source_negative_mass = path$source_negative_mass,
    shifted_negative_mass = path$shifted_negative_mass,
    shifted_truncated_positive_mass = projection$truncated_positive_mass[[1L]])
}
