# Exact structured sampling-covariance forms -------------------------------
#
# One class carries every structured sampling covariance. What differs between
# an RDM covariance and the general evidence form is not its fields but the
# basis its coordinates live in, so `basis` is a field of the form rather than
# a subclass: `"rdm"` for crossvalidated squared distances, `"evidence"` for
# the general evidence coordinates the kernel builds, `"query_bank"` for a
# declared bank of contrast-energy queries reduced out of the distance basis.
# Readers of the object (print, format, plot, the query entry point) branch on
# that field, and the basis enters the content signature because two forms over
# the same numbers in different bases are different artifacts.

.sampling_covariance_bases <- c("evidence", "rdm", "query_bank")

.sampling_covariance_signature <- function(fields) {
  .sha256_signature(list(
    schema_version = 1L,
    contract = "evidence-sampling-v1",
    plan = fields$plan$scientific_plan_id,
    signal_factor = fields$signal_factor,
    xi_factor = fields$xi_factor,
    noise_trace = fields$noise_trace,
    partitions = fields$partitions,
    labels = fields$labels,
    basis = fields$basis,
    source = fields$source
  ))
}

.sampling_covariance_form <- function(plan, signal_factor, xi_factor,
                                      noise_trace,
                                      partitions, labels = NULL,
                                      source = list(),
                                      basis = "evidence") {
  .validate_evidence_sampling_plan(plan, deep = FALSE)
  .require_sampling_covariance(plan)
  if (!.is_finite_matrix(signal_factor) || nrow(signal_factor) < 1L ||
      ncol(signal_factor) < 1L || !.is_finite_matrix(xi_factor) ||
      nrow(xi_factor) != nrow(signal_factor) || ncol(xi_factor) < 1L) {
    .input_error(paste0(
      "Sampling covariance requires finite row-factor matrices with one ",
      "shared evidence axis."
    ))
  }
  if (!.is_number(noise_trace) || noise_trace < 0) {
    .input_error(
      "The residual noise trace must be one finite nonnegative value."
    )
  }
  if (!.is_number(partitions) || partitions %% 1 != 0 || partitions < 2L ||
      as.integer(partitions) != plan$partition$count) {
    .contract_error(
      "The covariance partition count must match the sampling plan."
    )
  }
  dimension <- nrow(signal_factor)
  if (is.null(labels)) labels <- paste0("evidence", seq_len(dimension))
  if (!.is_strings(labels, unique = TRUE) || length(labels) != dimension) {
    .input_error(paste0(
      "Sampling covariance labels must uniquely identify every evidence ",
      "coordinate."
    ))
  }
  if (!is.list(source)) {
    .input_error("Sampling covariance source metadata must be a list.")
  }
  if (!.is_string(basis) || !basis %in% .sampling_covariance_bases) {
    .input_error(sprintf(
      "Sampling covariance basis must be one of %s.",
      paste(sQuote(.sampling_covariance_bases, q = FALSE), collapse = ", ")
    ))
  }
  rownames(signal_factor) <- labels
  rownames(xi_factor) <- labels
  fields <- list(
    plan = plan,
    signal_factor = signal_factor,
    xi_factor = xi_factor,
    noise_trace = as.double(noise_trace),
    partitions = as.integer(partitions),
    labels = labels,
    dimension = as.integer(dimension),
    basis = basis,
    source = source,
    signature = NA_character_
  )
  fields$signature <- .sampling_covariance_signature(fields)
  value <- structure(fields, class = "effect_sampling_covariance")
  .validate_sampling_covariance(value, deep = TRUE)
  value
}

# A rank-revealing PSD square root: `value == root %*% t(root)` with one
# column per retained direction. The spectral factorization costs O(n^3) in
# the size of `value`, which for the residual covariance is the local support
# a caller has already materialized densely, not the whole domain.
.sampling_psd_root <- function(value, what, empty = c("refuse", "zero")) {
  empty <- match.arg(empty)
  spectrum <- eigen(value, symmetric = TRUE)
  scale <- max(1, max(abs(spectrum$values)))
  if (min(spectrum$values) < -1e-10 * scale) {
    .input_error(sprintf("%s must be positive semidefinite.", what))
  }
  retained <- spectrum$values > 1e-12 * scale
  if (!any(retained)) {
    if (identical(empty, "refuse")) {
      .input_error(sprintf("%s has no positive sampling direction.", what))
    }
    # A degenerate covariance contributes nothing; the rank-one zero column
    # keeps the factorized form well formed instead of losing its axis.
    return(matrix(0, nrow(value), 1L))
  }
  spectrum$vectors[, retained, drop = FALSE] %*%
    diag(sqrt(spectrum$values[retained]), sum(retained))
}

# The signal-independent term of the sampling law is QUADRATIC in the
# whitened residual covariance, through tr(Sigma_w^2). A caller that supplies
# the true Sigma_w may use that trace directly, but a caller that supplies the
# pooled sample residual covariance S_w must not: for S_w ~ W_P(nu, Sigma_w)/nu,
#
#   E tr(S_w^2) = ((nu + 1) / nu) tr(Sigma_w^2) + tr(Sigma_w)^2 / nu,
#
# which overstates tr(Sigma_w^2) by a factor that grows with the ratio of the
# support size to the residual degrees of freedom. On the standard-error scale
# the inflation is sqrt(1 + (1 + P_eff) / nu) with
# P_eff = tr(Sigma_w)^2 / tr(Sigma_w^2), so a 50-voxel searchlight at nu = 168
# reports standard errors 14% too large and an 800-voxel one more than twice
# too large. The Wishart-unbiased estimator of the quadratic functional is
#
#   trhat(Sigma_w^2) = nu^2 / ((nu - 1)(nu + 2)) (tr(S_w^2) - tr(S_w)^2 / nu),
#
# and `residual_df` is how a caller declares that its residual covariance is a
# plug-in with nu degrees of freedom. `NULL` means "this is the covariance
# itself", which is what an oracle or a known-Sigma caller supplies. The
# signal term is LINEAR in Sigma_w and needs no correction.
.sampling_unbiased_noise_trace <- function(residual_covariance,
                                           residual_df = NULL) {
  raw <- sum(residual_covariance * residual_covariance)
  total <- sum(diag(residual_covariance))
  # tr(Sigma_w)^2 / tr(Sigma_w^2) is the participation ratio: the number of
  # residual directions the support actually spends its variance on. It, not
  # the raw support size, is what the sufficiency check compares against nu,
  # and it is read off whichever estimator of tr(Sigma_w^2) is in force.
  if (is.null(residual_df)) {
    return(list(
      value = raw, corrected = FALSE, raw = raw,
      residual_df = NA_integer_,
      effective_dimension = if (raw > 0) total^2 / raw else 0
    ))
  }
  if (!.is_number(residual_df) || residual_df %% 1 != 0 || residual_df < 1L) {
    .input_error(paste0(
      "Residual degrees of freedom must be one positive integer when the ",
      "residual covariance is a plug-in estimate."
    ))
  }
  nu <- as.double(residual_df)
  # At nu = 1 the (nu - 1) denominator vanishes and the estimator is
  # undefined; `.require_sufficient_residual_df()` refuses rather than divide.
  value <- if (nu >= 2) {
    (nu^2 / ((nu - 1) * (nu + 2))) * (raw - total^2 / nu)
  } else {
    NA_real_
  }
  list(
    value = if (is.na(value)) NA_real_ else max(value, 0),
    corrected = TRUE, raw = raw,
    residual_df = as.integer(residual_df),
    effective_dimension = if (!is.na(value) && value > 0) {
      total^2 / value
    } else {
      Inf
    }
  )
}

# What the corrected estimator needs is fixed by the estimator itself, not by
# a comparison of nu with the size of the support:
#
#   * nu >= 2, because the (nu - 1) denominator vanishes at nu = 1; and
#   * a positive estimate. Since rank(S_w) <= nu, tr(S_w^2) >= tr(S_w)^2 / nu
#     always, so the bracket is never negative. It is zero only when S_w is
#     zero or exactly isotropic on nu directions, where the plug-in says
#     nothing about tr(Sigma_w^2) and a zero noise term would turn into a
#     confidently small standard error.
#
# There is deliberately NO requirement that P_eff <= nu. The estimator is
# ratio-consistent as P / nu -> infinity (Srivastava 2005), and its relative
# sampling error FALLS as P_eff grows. Simulated over Wishart draws with a
# mildly anisotropic Sigma at nu = 10, the coefficient of variation of the
# estimate of tr(Sigma_w^2) is about 0.94 at P_eff = 1, 0.55 at P_eff = 3.7,
# 0.26 at P_eff = 41 and 0.19 at P_eff = 327; at P = 800, nu = 168 it is
# about 0.014. The worst case is the voxelwise one (P_eff = 1), whose error is
# that of an ordinary variance-of-a-variance at nu df and which has always
# been admitted, so the old P_eff <= nu rule refused precisely the supports
# the estimator handles best. The floor is the one the formula sets, nu >= 2;
# no stricter policy is imposed on one support shape and not another.
.require_sufficient_residual_df <- function(quadratic) {
  if (!isTRUE(quadratic$corrected)) return(invisible(NULL))
  nu <- quadratic$residual_df
  if (nu < 2L) {
    .capability_refusal(sprintf(paste0(
      "This measurement's residual covariance is estimated from %d residual ",
      "degree of freedom, and the unbiased estimator of the noise term ",
      "tr(Sigma^2) of the sampling law needs at least two; it is undefined ",
      "here."
    ), nu),
      capability = "sufficient_residual_df",
      namespace = "evidence_sampling",
      reasons = "residual_df_below_two",
      remedies = paste0(
        "Fit more partitions or observations, which raises the residual df."
      )
    )
  }
  if (is.finite(quadratic$value) && quadratic$value > 0) {
    return(invisible(NULL))
  }
  .capability_refusal(sprintf(paste0(
    "This measurement's plug-in residual covariance (%d residual degrees of ",
    "freedom) is zero or exactly isotropic on its residual directions, so the ",
    "unbiased estimate of the noise term tr(Sigma^2) of the sampling law is ",
    "zero and carries no information. Reporting it would give a confidently ",
    "small standard error."
  ), nu),
    capability = "sufficient_residual_df",
    namespace = "evidence_sampling",
    reasons = "noise_trace_estimate_nonpositive",
    remedies = c(
      "Check that the residuals are not identically zero (a saturated fit).",
      "Fit more partitions or observations, which raises the residual df."
    )
  )
}

# Exact sampling law of the equal-weight all-partition-pairs crossvalidated
# distance estimator, in the coordinates crossform already whitens into.
# Write mu_r = c_r U L' for the whitened contrast pattern of distance r
# (`contrasts %*% signal_patterns`), Sigma_w for the whitened residual
# covariance, Xi = C Sigma_K C' for the effect-coordinate cross-products,
# M for the partition count, and nu for the distance normalization. Then
#
#   Cov(d_r, d_s) = (4 / M) Xi_rs (mu_r Sigma_w mu_s') / nu^2
#                 + 2 / (M (M - 1)) Xi_rs^2 tr(Sigma_w Sigma_w) / nu^2.
#
# `normalization` is nu: the divisor already applied to every distance. Both
# terms carry 1 / nu^2 because both are quadratic in the estimate. The
# product path folds its frame weights into the metric and passes nu = 1.
#
# `residual_df` declares that `residual_covariance` is a plug-in estimate with
# that many degrees of freedom, which the quadratic noise term corrects for;
# see `.sampling_unbiased_noise_trace()`.
.sampling_covariance_from_components <- function(
    plan, contrasts, signal_patterns, effect_covariance,
    residual_covariance, normalization = ncol(signal_patterns),
    residual_df = NULL, labels = NULL, source = list(),
    xi_factor = NULL, basis = "evidence",
    signal_target = c("fixed", "partition_mean_plugin")) {
  .validate_evidence_sampling_plan(plan, deep = FALSE)
  .require_sampling_covariance(plan)
  plan_coordinates <-
    plan$evidence_plan$task$left_relation$effect_space$coordinates
  plan_effects <- length(plan_coordinates)
  if (is.matrix(contrasts) && ncol(contrasts) != plan_effects) {
    .contract_error(sprintf(paste0(
      "Sampling contrasts declare %d experimental effects but the bound ",
      "evidence plan has %d; a covariance artifact must carry the identity ",
      "of the plan it claims."
    ), ncol(contrasts), plan_effects))
  }
  if (is.matrix(contrasts) && !is.null(colnames(contrasts)) &&
      !identical(colnames(contrasts), plan_coordinates)) {
    .contract_error(paste0(
      "Sampling contrast effect names do not match the bound evidence ",
      "plan's effect space."
    ))
  }
  if (!.is_finite_matrix(contrasts) || nrow(contrasts) < 1L ||
      ncol(contrasts) < 2L || !.is_finite_matrix(signal_patterns) ||
      nrow(signal_patterns) != ncol(contrasts) || ncol(signal_patterns) < 1L ||
      !.is_finite_matrix(effect_covariance) ||
      !identical(dim(effect_covariance), as.integer(rep(ncol(contrasts), 2L))) ||
      !.sampling_symmetric_within(effect_covariance) ||
      !.is_finite_matrix(residual_covariance) ||
      nrow(residual_covariance) != ncol(signal_patterns) ||
      ncol(residual_covariance) != ncol(signal_patterns) ||
      !.sampling_symmetric_within(residual_covariance)) {
    .input_error("Sampling components have incompatible or non-finite axes.")
  }
  # Both arrive as products (`R S R'`, `C Sigma_K C'`) whose rounding leaves a
  # last-bit asymmetry; the tolerance above is relative to their scale, and
  # the exact symmetrization here makes everything downstream see one matrix.
  effect_covariance <- 0.5 * (effect_covariance + t(effect_covariance))
  residual_covariance <- 0.5 * (residual_covariance + t(residual_covariance))
  signal_target <- match.arg(signal_target)
  if (!.is_number(normalization) || normalization <= 0) {
    .input_error(
      "Sampling covariance normalization must be one positive finite value."
    )
  }
  if (is.null(xi_factor)) {
    effect_root <- .sampling_psd_root(
      effect_covariance, "Effect covariance", empty = "refuse"
    )
    xi_factor <- contrasts %*% effect_root
  } else if (!.is_finite_matrix(xi_factor) ||
      nrow(xi_factor) != nrow(contrasts)) {
    .contract_error(
      "A hoisted sampling xi factor must share the contrast evidence axis."
    )
  }
  # The signal term of the exact law is
  #
  #   (4 / M) Xi_rs (mu_r Sigma_w mu_s') / nu^2,
  #
  # so the residual covariance enters the signal term as a genuine metric on
  # the whitened contrast patterns, not as an isotropic scalar. Factoring
  # Sigma_w = L L' keeps the rank-preserving row-factor form the query
  # machinery requires while carrying the full anisotropy of Sigma_w.
  residual_root <- .sampling_psd_root(
    residual_covariance, "Residual covariance", empty = "zero"
  )
  features <- ncol(signal_patterns)
  quadratic <- .sampling_unbiased_noise_trace(
    residual_covariance, residual_df
  )
  .require_sufficient_residual_df(quadratic)
  noise_trace <- quadratic$value / normalization^2
  if (identical(signal_target, "partition_mean_plugin")) {
    debiased <- .sampling_debiased_plugin_signal(
      signal_patterns, residual_root, effect_covariance, quadratic$value,
      plan$partition$count
    )
    signal_factor <- (contrasts %*% debiased$root) / normalization
  } else {
    debiased <- NULL
    signal_factor <- (contrasts %*% signal_patterns %*% residual_root) /
      normalization
  }
  .sampling_covariance_form(
    plan,
    signal_factor = signal_factor,
    xi_factor = xi_factor,
    noise_trace = noise_trace,
    partitions = plan$partition$count,
    labels = labels,
    basis = basis,
    source = c(source, list(
      construction = "diedrichsen_eq13_components_general_metric",
      signal_rank = ncol(signal_factor),
      effect_covariance_rank = ncol(xi_factor),
      residual_dimension = as.integer(features),
      normalization = as.double(normalization),
      noise_trace_estimator = if (quadratic$corrected) {
        "wishart_unbiased_quadratic"
      } else {
        "known_residual_covariance"
      },
      residual_effective_dimension = quadratic$effective_dimension
    ), if (!is.null(debiased)) list(
      signal_target_correction = "partition_mean_bias_removed",
      signal_clamped_directions = debiased$clamped
    ))
  )
}

# Relative symmetry: a covariance built as `R S R'` from data on any scale
# carries a rounding asymmetry proportional to its own magnitude, so an
# absolute tolerance refuses well-formed input merely for being large.
.sampling_symmetric_within <- function(value, tolerance = 1e-10) {
  scale <- max(abs(value))
  max(abs(value - t(value))) <= tolerance * scale
}

# The `partition_mean_plugin` target substitutes the partition mean of the
# whitened estimates, B_bar, for the unknown signal. Each partition estimate
# is B + E_m with Cov(vec E_m) = Sigma_K (x) Sigma_w, so B_bar carries noise
# Sigma_K (x) Sigma_w / M and, with S_w independent of B_bar and unbiased,
#
#   E[B_bar S_w B_bar'] = B Sigma_w B' + Sigma_K tr(Sigma_w^2) / M.
#
# Taken at face value the signal term (4 / M) Xi_rs (C H C')_rs therefore
# carries an extra 4 Xi_rs^2 tr(Sigma_w^2) / M^2. Against the noise term
# 2 Xi_rs^2 tr(Sigma_w^2) / (M (M - 1)) that is a ratio of 2 (M - 1) / M, so
# the plug-in roughly TRIPLES the noise contribution under weak signal; it is
# not a mild O(M^-2) effect. The unbiased signal Gram subtracts the bias with
# the unbiased estimate of tr(Sigma_w^2) already used by the noise term:
#
#   H_hat = B_bar S_w B_bar' - Sigma_K trhat(Sigma_w^2) / M.
#
# H_hat need not be positive semidefinite (it rarely is when the true signal
# is weak), and an indefinite H_hat can make whole variances negative, which
# the factorized form cannot carry. It is projected onto the PSD cone in the
# Sigma_K-whitened coordinates: with Z = Sigma_K^{-1/2} H Sigma_K^{-1/2},
# the eigenvalues of Z are replaced by max(lambda - trhat / M, 0). That is the
# positive-part estimator of the signal-to-noise spectrum. It is exactly
# unbiased whenever no eigenvalue is clamped, and in the Loewner order it
# always lies between H_hat and the raw plug-in H, so it is never more biased
# upward than the plug-in it replaces. Where Sigma_K is singular the
# projection falls back to the Euclidean one on H_hat.
.sampling_debiased_plugin_signal <- function(signal_patterns, residual_root,
                                             effect_covariance, noise_trace,
                                             partitions) {
  pattern <- signal_patterns %*% residual_root
  gram <- tcrossprod(pattern)
  gram <- 0.5 * (gram + t(gram))
  shift <- noise_trace / partitions
  spectrum <- eigen(effect_covariance, symmetric = TRUE)
  scale <- max(abs(spectrum$values))
  if (min(spectrum$values) > 1e-10 * scale) {
    root <- spectrum$vectors %*% diag(sqrt(spectrum$values),
      length(spectrum$values)) %*% t(spectrum$vectors)
    inverse_root <- spectrum$vectors %*% diag(1 / sqrt(spectrum$values),
      length(spectrum$values)) %*% t(spectrum$vectors)
    whitened <- inverse_root %*% gram %*% inverse_root
    whitened <- 0.5 * (whitened + t(whitened))
    inner <- eigen(whitened, symmetric = TRUE)
    values <- inner$values - shift
    vectors <- root %*% inner$vectors
  } else {
    inner <- eigen(gram - effect_covariance * shift, symmetric = TRUE)
    values <- inner$values
    vectors <- inner$vectors
  }
  retained <- values > 0
  clamped <- sum(!retained)
  if (!any(retained)) {
    return(list(root = matrix(0, nrow(gram), 1L), clamped = clamped))
  }
  list(
    root = vectors[, retained, drop = FALSE] %*%
      diag(sqrt(values[retained]), sum(retained)),
    clamped = as.integer(clamped)
  )
}

.validate_sampling_covariance <- function(x, deep = TRUE) {
  if (.validated_before(x, "sampling_covariance", deep)) return(invisible(x))
  expected <- c("plan", "signal_factor", "xi_factor", "noise_trace", "partitions",
    "labels", "dimension", "basis", "source", "signature")
  if (!.sealed_fields(x, "effect_sampling_covariance", expected) ||
      !.is_string(x$basis) || !x$basis %in% .sampling_covariance_bases ||
      !is.matrix(x$signal_factor) || !is.matrix(x$xi_factor) ||
      nrow(x$signal_factor) != nrow(x$xi_factor) ||
      !.is_number(x$noise_trace) || x$noise_trace < 0 ||
      !is.integer(x$partitions) || length(x$partitions) != 1L ||
      x$partitions < 2L || !.is_strings(x$labels, unique = TRUE) ||
      length(x$labels) != nrow(x$signal_factor) || !is.integer(x$dimension) ||
      length(x$dimension) != 1L || x$dimension != nrow(x$signal_factor) ||
      !is.list(x$source) || !.strong_sha256(x$signature)) {
    .input_error("Sampling-covariance fields are missing or noncanonical.")
  }
  .validate_evidence_sampling_plan(x$plan, deep = FALSE)
  .require_sampling_covariance(x$plan)
  if (x$partitions != x$plan$partition$count ||
      !identical(rownames(x$signal_factor), x$labels) ||
      !identical(rownames(x$xi_factor), x$labels) ||
      any(!is.finite(x$signal_factor)) || any(!is.finite(x$xi_factor))) {
    .contract_error(
      "Sampling-covariance axes, plan, or row factors are inconsistent."
    )
  }
  if (isTRUE(deep) && !identical(x$signature,
      .sampling_covariance_signature(x))) {
    .contract_error("Sampling-covariance identity is inconsistent.")
  }
  .record_validated(x, "sampling_covariance", deep)
  invisible(x)
}

.sampling_covariance_coefficients <- function(x) {
  .validate_sampling_covariance(x, deep = FALSE)
  list(
    signal = 4 / x$partitions,
    noise = 2 * x$noise_trace /
      (x$partitions * (x$partitions - 1L))
  )
}

.sampling_row_inner_products <- function(factor, row, column) {
  unname(rowSums(
    factor[row, , drop = FALSE] * factor[column, , drop = FALSE]
  ))
}

.sampling_covariance_entries <- function(x, row, column) {
  .validate_sampling_covariance(x, deep = FALSE)
  if (!.is_finite_numeric(row) || !.is_finite_numeric(column) ||
      length(row) < 1L || !identical(length(row), length(column)) ||
      anyNA(row) || anyNA(column) || any(row %% 1 != 0) || any(column %% 1 != 0) ||
      any(row < 1L) || any(column < 1L) || any(row > x$dimension) ||
      any(column > x$dimension)) {
    .input_error("Sampling-covariance entry indices are invalid.")
  }
  row <- as.integer(row)
  column <- as.integer(column)
  coefficient <- .sampling_covariance_coefficients(x)
  signal <- .sampling_row_inner_products(x$signal_factor, row, column)
  xi <- .sampling_row_inner_products(x$xi_factor, row, column)
    coefficient$signal *
    signal * xi + coefficient$noise * xi^2
}

.sampling_align_covariance_value <- function(value, labels,
                                             axis = c("rows", "columns"),
                                             what = "sampling query") {
  axis <- match.arg(axis)
  identifiers <- if (is.null(dim(value))) {
    names(value)
  } else if (identical(axis, "rows")) {
    rownames(value)
  } else {
    colnames(value)
  }
  if (is.null(identifiers)) return(value)
  if (!.is_strings(identifiers, unique = TRUE) ||
      !setequal(identifiers, labels)) {
    .input_error(sprintf(
      "Named %s axes must identify every sampling-covariance coordinate exactly once.",
      what
    ))
  }
  if (is.null(dim(value))) {
    value[labels]
  } else if (identical(axis, "rows")) {
    value[labels, , drop = FALSE]
  } else {
    value[, labels, drop = FALSE]
  }
}

.sampling_hadamard_gram_apply <- function(left, right, value,
                                          workspace_bytes = 32 * 1024^2) {
  dimension <- nrow(left)
  rank_left <- ncol(left)
  rank_right <- ncol(right)
  if (!.is_number(workspace_bytes) || workspace_bytes < 8) {
    .input_error(
      "The Hadamard-Gram workspace budget must be positive and finite."
    )
  }
  # For one right-hand side v,
  #
  #   ((L L') o (R R')) v
  #     = row_i { l_i' [L' diag(v) R] r_i }.
  #
  # This is the same rowwise-tensor contraction as the literal formula, but
  # it uses only the compact rL-by-rR middle product.  In particular, neither
  # the D-by-D covariance nor the D-by-(rL*rR) tensor factor is materialized.
  result <- matrix(0, dimension, ncol(value))
  middle_bytes <- 8 * as.double(rank_left) * rank_right
  if (!is.finite(middle_bytes) || middle_bytes > workspace_bytes) {
    .input_error(sprintf(
      paste0("Exact Hadamard-Gram application requires at least %.0f bytes ",
        "for its compact interaction workspace, exceeding the %.0f-byte budget."),
      middle_bytes, workspace_bytes
    ))
  }
  temporary_rank <- min(rank_left, rank_right)
  available <- max(8, workspace_bytes - middle_bytes)
  rows_per_block <- max(1L, min(dimension,
    as.integer(floor(available / max(16 * temporary_rank, 8)))))
  starts <- seq.int(1L, dimension, by = rows_per_block)
  for (column in seq_len(ncol(value))) {
    weights <- value[, column]
    if (rank_left <= rank_right) {
      middle <- crossprod(left * weights, right)
    } else {
      middle <- crossprod(left, right * weights)
    }
    for (first in starts) {
      rows <- first:min(dimension, first + rows_per_block - 1L)
      if (rank_right <= rank_left) {
        temporary <- left[rows, , drop = FALSE] %*% middle
        result[rows, column] <- rowSums(
          temporary * right[rows, , drop = FALSE]
        )
      } else {
        temporary <- right[rows, , drop = FALSE] %*% t(middle)
        result[rows, column] <- rowSums(
          temporary * left[rows, , drop = FALSE]
        )
      }
    }
  }
  result
}

.sampling_hadamard_gram_quadratic <- function(left, right, value) {
  # v'[(L L') o (R R')]v = ||L' diag(v) R||_F^2.
  if (ncol(left) <= ncol(right)) {
    middle <- crossprod(left * value, right)
  } else {
    middle <- crossprod(left, right * value)
  }
  sum(middle * middle)
}

.sampling_hadamard_gram_transport <- function(left, right, map,
                                              workspace_bytes = 64 * 1024^2) {
  output_dimension <- nrow(map)
  pair_rank <- as.double(ncol(left)) * ncol(right)
  required <- 8 * pair_rank * output_dimension
  if (!is.finite(required) || required > workspace_bytes) {
    .input_error(sprintf(
      paste0("Exact Hadamard-Gram transport requires %.0f bytes for its ",
        "factorized output workspace, exceeding the %.0f-byte budget."),
      required, workspace_bytes
    ))
  }
  # If a_j is one output row, its transported row-tensor factor is
  # vec(L' diag(a_j) R). Their Gram matrix is exactly
  # A[(LL') o (RR')]A', without constructing the evidence covariance.
  factor <- matrix(0, output_dimension, as.integer(pair_rank))
  for (output in seq_len(output_dimension)) {
    weights <- map[output, ]
    middle <- if (ncol(left) <= ncol(right)) {
      crossprod(left * weights, right)
    } else {
      crossprod(left, right * weights)
    }
    factor[output, ] <- as.vector(middle)
  }
  tcrossprod(factor)
}

.sampling_covariance_apply <- function(x, value) {
  .validate_sampling_covariance(x, deep = FALSE)
  vector_input <- is.atomic(value) && is.null(dim(value))
  value <- .sampling_align_covariance_value(
    value, x$labels, axis = "rows", what = "covariance-action"
  )
  if (vector_input) value <- matrix(value, ncol = 1L)
  if (!.is_finite_matrix(value) || nrow(value) != x$dimension ||
      ncol(value) < 1L) {
    .input_error(
      "Sampling covariance can act only on finite evidence-coordinate vectors."
    )
  }
  coefficient <- .sampling_covariance_coefficients(x)
  result <- coefficient$signal * .sampling_hadamard_gram_apply(
    x$signal_factor, x$xi_factor, value
  ) + coefficient$noise * .sampling_hadamard_gram_apply(
    x$xi_factor, x$xi_factor, value
  )
  rownames(result) <- x$labels
  if (vector_input) {
    stats::setNames(drop(result), x$labels)
  } else {
    result
  }
}

.sampling_covariance_diagonal <- function(x) {
  .validate_sampling_covariance(x, deep = FALSE)
  values <- .sampling_covariance_entries(
    x, seq_len(x$dimension), seq_len(x$dimension)
  )
  stats::setNames(values, x$labels)
}

.sampling_covariance_quadratic <- function(x, value) {
  .validate_sampling_covariance(x, deep = FALSE)
  value <- .sampling_align_covariance_value(
    value, x$labels, axis = "rows", what = "quadratic-form"
  )
  if (!.is_finite_numeric(value) || is.matrix(value) ||
      length(value) != x$dimension) {
    .input_error(paste0(
      "A covariance quadratic form requires one finite evidence-coordinate ",
      "vector."
    ))
  }
  coefficient <- .sampling_covariance_coefficients(x)
  coefficient$signal * .sampling_hadamard_gram_quadratic(
    x$signal_factor, x$xi_factor, value
  ) + coefficient$noise * .sampling_hadamard_gram_quadratic(
    x$xi_factor, x$xi_factor, value
  )
}

.sampling_covariance_transport <- function(x, map) {
  .validate_sampling_covariance(x, deep = FALSE)
  map <- .sampling_align_covariance_value(
    map, x$labels, axis = "columns", what = "transport-input"
  )
  if (!.is_finite_matrix(map) || ncol(map) != x$dimension || nrow(map) < 1L) {
    .input_error(
      "A covariance transport must be a finite output-by-evidence matrix."
    )
  }
  coefficient <- .sampling_covariance_coefficients(x)
  value <- coefficient$signal * .sampling_hadamard_gram_transport(
    x$signal_factor, x$xi_factor, map
  ) + coefficient$noise * .sampling_hadamard_gram_transport(
    x$xi_factor, x$xi_factor, map
  )
  value <- 0.5 * (value + t(value))
  if (!is.null(rownames(map))) {
    dimnames(value) <- list(rownames(map), rownames(map))
  }
  value
}

.sampling_covariance_materialize <- function(
    x, max_bytes = 512 * 1024^2) {
  .validate_sampling_covariance(x, deep = FALSE)
  if (!.is_number(max_bytes) || max_bytes < 8) {
    .input_error(
      "`max_bytes` must be one finite positive materialization budget."
    )
  }
  output_bytes <- 8 * as.double(x$dimension)^2
  # One exact row block needs the signal Gram, Xi Gram, and result workspace.
  # Reserve conservatively for R's intermediate arithmetic as well.
  minimum_workspace <- 64 * as.double(x$dimension)
  required <- output_bytes + minimum_workspace
  if (!is.finite(required) || required > max_bytes) {
    .input_error(sprintf(
      paste0("Dense sampling covariance requires at least %.0f bytes including ",
        "one exact working row, exceeding the %.0f-byte materialization budget."),
      required, max_bytes
    ))
  }
  coefficient <- .sampling_covariance_coefficients(x)
  available <- max_bytes - output_bytes
  rows_per_block <- max(1L, min(x$dimension,
    as.integer(floor(available / (64 * x$dimension)))))
  value <- matrix(0, x$dimension, x$dimension)
  starts <- seq.int(1L, x$dimension, by = rows_per_block)
  for (first in starts) {
    rows <- first:min(x$dimension, first + rows_per_block - 1L)
    signal <- tcrossprod(
      x$signal_factor[rows, , drop = FALSE], x$signal_factor
    )
    xi <- tcrossprod(x$xi_factor[rows, , drop = FALSE], x$xi_factor)
    value[rows, ] <- coefficient$signal * signal * xi +
      coefficient$noise * xi^2
  }
  value <- 0.5 * (value + t(value))
  dimnames(value) <- list(x$labels, x$labels)
  value
}

.execute_evidence_sampling_plan <- function(plan, covariance,
                                            max_bytes = 512 * 1024^2) {
  .validate_evidence_sampling_plan(plan)
  .validate_sampling_covariance(covariance, deep = FALSE)
  if (!identical(plan$evidence$signature,
      covariance$plan$evidence$signature) ||
      !identical(plan$error_channel$signature,
        covariance$plan$error_channel$signature) ||
      !identical(plan$metric$signature, covariance$plan$metric$signature) ||
      !identical(plan$partition$signature,
        covariance$plan$partition$signature) ||
      !identical(plan$target$signature, covariance$plan$target$signature) ||
      !identical(plan$sampling_axis, covariance$plan$sampling_axis) ||
      !identical(plan$spatial_scope, covariance$plan$spatial_scope)) {
    .contract_error(paste0(
      "The sampling operation and covariance source have different scientific ",
      "identities."
    ))
  }
  kind <- plan$operation$operation
  argument <- plan$operation$argument
  switch(kind,
    diagonal = .sampling_covariance_diagonal(covariance),
    selected_entries = {
      if (!is.matrix(argument) || ncol(argument) != 2L) {
        .input_error("Selected entries require a two-column index matrix.")
      }
      .sampling_covariance_entries(covariance, argument[, 1L], argument[, 2L])
    },
    apply = .sampling_covariance_apply(covariance, argument),
    quadratic_form = .sampling_covariance_quadratic(covariance, argument),
    transport = .sampling_covariance_transport(covariance, argument),
    materialize = .sampling_covariance_materialize(covariance, max_bytes),
    .input_error("Unknown sampling-covariance operation.")
  )
}
