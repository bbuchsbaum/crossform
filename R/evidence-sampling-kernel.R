# Exact structured sampling-covariance forms -------------------------------

.sampling_covariance_signature <- function(fields) {
  paste0("sha256:", digest::digest(list(
    schema_version = 1L,
    contract = "evidence-sampling-v1",
    plan = fields$plan$scientific_plan_id,
    signal_factor = fields$signal_factor,
    xi_factor = fields$xi_factor,
    noise_trace = fields$noise_trace,
    partitions = fields$partitions,
    labels = fields$labels,
    source = fields$source
  ), algo = "sha256", serialize = TRUE, serializeVersion = 2L))
}

.sampling_covariance_form <- function(plan, signal_factor, xi_factor,
                                      noise_trace,
                                      partitions, labels = NULL,
                                      source = list()) {
  .validate_evidence_sampling_plan(plan, deep = FALSE)
  .require_sampling_covariance(plan)
  if (!is.matrix(signal_factor) || !is.numeric(signal_factor) ||
      nrow(signal_factor) < 1L || ncol(signal_factor) < 1L ||
      any(!is.finite(signal_factor)) || !is.matrix(xi_factor) ||
      !is.numeric(xi_factor) || nrow(xi_factor) != nrow(signal_factor) ||
      ncol(xi_factor) < 1L || any(!is.finite(xi_factor))) {
    stop("Sampling covariance requires finite row-factor matrices with one shared evidence axis.",
      call. = FALSE)
  }
  if (!is.numeric(noise_trace) || length(noise_trace) != 1L ||
      is.na(noise_trace) || !is.finite(noise_trace) ||
      noise_trace < 0) {
    stop("The residual noise trace must be one finite nonnegative value.",
      call. = FALSE)
  }
  if (!is.numeric(partitions) || length(partitions) != 1L ||
      is.na(partitions) || !is.finite(partitions) ||
      partitions %% 1 != 0 || partitions < 2L ||
      as.integer(partitions) != plan$partition$count) {
    stop("The covariance partition count must match the sampling plan.",
      call. = FALSE)
  }
  dimension <- nrow(signal_factor)
  if (is.null(labels)) labels <- paste0("evidence", seq_len(dimension))
  if (!is.character(labels) || length(labels) != dimension || anyNA(labels) ||
      any(!nzchar(labels)) || anyDuplicated(labels)) {
    stop("Sampling covariance labels must uniquely identify every evidence coordinate.",
      call. = FALSE)
  }
  if (!is.list(source)) {
    stop("Sampling covariance source metadata must be a list.", call. = FALSE)
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
    stop(sprintf("%s must be positive semidefinite.", what), call. = FALSE)
  }
  retained <- spectrum$values > 1e-12 * scale
  if (!any(retained)) {
    if (identical(empty, "refuse")) {
      stop(sprintf("%s has no positive sampling direction.", what),
        call. = FALSE)
    }
    # A degenerate covariance contributes nothing; the rank-one zero column
    # keeps the factorized form well formed instead of losing its axis.
    return(matrix(0, nrow(value), 1L))
  }
  spectrum$vectors[, retained, drop = FALSE] %*%
    diag(sqrt(spectrum$values[retained]), sum(retained))
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
.sampling_covariance_from_components <- function(
    plan, contrasts, signal_patterns, effect_covariance,
    residual_covariance, normalization = ncol(signal_patterns),
    labels = NULL, source = list()) {
  .validate_evidence_sampling_plan(plan, deep = FALSE)
  .require_sampling_covariance(plan)
  plan_coordinates <-
    plan$evidence_plan$task$left_relation$effect_space$coordinates
  plan_effects <- length(plan_coordinates)
  if (is.matrix(contrasts) && ncol(contrasts) != plan_effects) {
    stop(sprintf(paste0(
      "Sampling contrasts declare %d experimental effects but the bound ",
      "evidence plan has %d; a covariance artifact must carry the identity ",
      "of the plan it claims."
    ), ncol(contrasts), plan_effects), call. = FALSE)
  }
  if (is.matrix(contrasts) && !is.null(colnames(contrasts)) &&
      !identical(colnames(contrasts), plan_coordinates)) {
    stop(paste0(
      "Sampling contrast effect names do not match the bound evidence ",
      "plan's effect space."
    ), call. = FALSE)
  }
  if (!is.matrix(contrasts) || !is.numeric(contrasts) ||
      nrow(contrasts) < 1L || ncol(contrasts) < 2L ||
      any(!is.finite(contrasts)) ||
      !is.matrix(signal_patterns) || !is.numeric(signal_patterns) ||
      nrow(signal_patterns) != ncol(contrasts) ||
      ncol(signal_patterns) < 1L || any(!is.finite(signal_patterns)) ||
      !is.matrix(effect_covariance) || !is.numeric(effect_covariance) ||
      !identical(dim(effect_covariance),
        as.integer(rep(ncol(contrasts), 2L))) ||
      any(!is.finite(effect_covariance)) ||
      max(abs(effect_covariance - t(effect_covariance))) > 1e-10 ||
      !is.matrix(residual_covariance) ||
      !is.numeric(residual_covariance) ||
      nrow(residual_covariance) != ncol(signal_patterns) ||
      ncol(residual_covariance) != ncol(signal_patterns) ||
      any(!is.finite(residual_covariance)) ||
      max(abs(residual_covariance - t(residual_covariance))) > 1e-10) {
    stop("Sampling components have incompatible or non-finite axes.",
      call. = FALSE)
  }
  if (!is.numeric(normalization) || length(normalization) != 1L ||
      is.na(normalization) || !is.finite(normalization) ||
      normalization <= 0) {
    stop("Sampling covariance normalization must be one positive finite value.",
      call. = FALSE)
  }
  effect_root <- .sampling_psd_root(
    effect_covariance, "Effect covariance", empty = "refuse"
  )
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
  noise_trace <- sum(residual_covariance * residual_covariance) /
    normalization^2
  signal_factor <- (contrasts %*% signal_patterns %*% residual_root) /
    normalization
  xi_factor <- contrasts %*% effect_root
  .sampling_covariance_form(
    plan,
    signal_factor = signal_factor,
    xi_factor = xi_factor,
    noise_trace = noise_trace,
    partitions = plan$partition$count,
    labels = labels,
    source = c(source, list(
      construction = "diedrichsen_eq13_components_general_metric",
      signal_rank = ncol(signal_factor),
      effect_covariance_rank = ncol(xi_factor),
      residual_dimension = as.integer(features),
      normalization = as.double(normalization)
    ))
  )
}

.validate_sampling_covariance <- function(x, deep = TRUE) {
  if (.validated_before(x, "sampling_covariance", deep)) return(invisible(x))
  expected <- c("plan", "signal_factor", "xi_factor", "noise_trace", "partitions",
    "labels", "dimension", "source", "signature")
  if (!inherits(x, "effect_sampling_covariance") || !is.list(x) ||
      !identical(names(x), expected) || !is.matrix(x$signal_factor) ||
      !is.matrix(x$xi_factor) ||
      nrow(x$signal_factor) != nrow(x$xi_factor) ||
      !is.numeric(x$noise_trace) || length(x$noise_trace) != 1L ||
      is.na(x$noise_trace) || !is.finite(x$noise_trace) ||
      x$noise_trace < 0 || !is.integer(x$partitions) ||
      length(x$partitions) != 1L || x$partitions < 2L ||
      !is.character(x$labels) || length(x$labels) != nrow(x$signal_factor) ||
      anyNA(x$labels) || any(!nzchar(x$labels)) || anyDuplicated(x$labels) ||
      !is.integer(x$dimension) || length(x$dimension) != 1L ||
      x$dimension != nrow(x$signal_factor) || !is.list(x$source) ||
      !.strong_sha256(x$signature)) {
    stop("Sampling-covariance fields are missing or noncanonical.",
      call. = FALSE)
  }
  .validate_evidence_sampling_plan(x$plan, deep = FALSE)
  .require_sampling_covariance(x$plan)
  if (x$partitions != x$plan$partition$count ||
      !identical(rownames(x$signal_factor), x$labels) ||
      !identical(rownames(x$xi_factor), x$labels) ||
      any(!is.finite(x$signal_factor)) || any(!is.finite(x$xi_factor))) {
    stop("Sampling-covariance axes, plan, or row factors are inconsistent.",
      call. = FALSE)
  }
  if (isTRUE(deep) && !identical(x$signature,
      .sampling_covariance_signature(x))) {
    stop("Sampling-covariance identity is inconsistent.", call. = FALSE)
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
  if (!is.numeric(row) || !is.numeric(column) || length(row) < 1L ||
      !identical(length(row), length(column)) || anyNA(row) ||
      anyNA(column) || any(!is.finite(row)) || any(!is.finite(column)) ||
      any(row %% 1 != 0) || any(column %% 1 != 0) || any(row < 1L) ||
      any(column < 1L) || any(row > x$dimension) ||
      any(column > x$dimension)) {
    stop("Sampling-covariance entry indices are invalid.", call. = FALSE)
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
  if (!is.character(identifiers) || anyNA(identifiers) ||
      any(!nzchar(identifiers)) || anyDuplicated(identifiers) ||
      !setequal(identifiers, labels)) {
    stop(sprintf(
      "Named %s axes must identify every sampling-covariance coordinate exactly once.",
      what
    ), call. = FALSE)
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
  if (!is.numeric(workspace_bytes) || length(workspace_bytes) != 1L ||
      is.na(workspace_bytes) || !is.finite(workspace_bytes) ||
      workspace_bytes < 8) {
    stop("The Hadamard-Gram workspace budget must be positive and finite.",
      call. = FALSE)
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
    stop(sprintf(
      paste0("Exact Hadamard-Gram application requires at least %.0f bytes ",
        "for its compact interaction workspace, exceeding the %.0f-byte budget."),
      middle_bytes, workspace_bytes
    ), call. = FALSE)
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
    stop(sprintf(
      paste0("Exact Hadamard-Gram transport requires %.0f bytes for its ",
        "factorized output workspace, exceeding the %.0f-byte budget."),
      required, workspace_bytes
    ), call. = FALSE)
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
  if (!is.matrix(value) || !is.numeric(value) ||
      nrow(value) != x$dimension || ncol(value) < 1L ||
      any(!is.finite(value))) {
    stop("Sampling covariance can act only on finite evidence-coordinate vectors.",
      call. = FALSE)
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
  if (!is.numeric(value) || is.matrix(value) ||
      length(value) != x$dimension || any(!is.finite(value))) {
    stop("A covariance quadratic form requires one finite evidence-coordinate vector.",
      call. = FALSE)
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
  if (!is.matrix(map) || !is.numeric(map) || ncol(map) != x$dimension ||
      nrow(map) < 1L || any(!is.finite(map))) {
    stop("A covariance transport must be a finite output-by-evidence matrix.",
      call. = FALSE)
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
  if (!is.numeric(max_bytes) || length(max_bytes) != 1L ||
      is.na(max_bytes) || !is.finite(max_bytes) || max_bytes < 8) {
    stop("`max_bytes` must be one finite positive materialization budget.",
      call. = FALSE)
  }
  output_bytes <- 8 * as.double(x$dimension)^2
  # One exact row block needs the signal Gram, Xi Gram, and result workspace.
  # Reserve conservatively for R's intermediate arithmetic as well.
  minimum_workspace <- 64 * as.double(x$dimension)
  required <- output_bytes + minimum_workspace
  if (!is.finite(required) || required > max_bytes) {
    stop(sprintf(
      paste0("Dense sampling covariance requires at least %.0f bytes including ",
        "one exact working row, exceeding the %.0f-byte materialization budget."),
      required, max_bytes
    ), call. = FALSE)
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

.sampling_covariance_receipt <- function(x) {
  .validate_sampling_covariance(x, deep = FALSE)
  structure(list(
    contract = x$plan$contract,
    sampling_plan_id = x$plan$scientific_plan_id,
    evidence_plan_id = x$plan$evidence$plan_id,
    relation_fit = x$plan$error_channel$channel_identity,
    metric = x$plan$metric$identity,
    metric_status = x$plan$metric$status,
    partition_model = x$plan$partition$model,
    partitions = x$partitions,
    sampling_axis = x$plan$sampling_axis,
    target = x$plan$target$target,
    target_policy = x$plan$target$policy,
    spatial_scope = x$plan$spatial_scope,
    source = x$source,
    covariance = x$signature
  ), class = "effect_sampling_receipt")
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
    stop("The sampling operation and covariance source have different scientific identities.",
      call. = FALSE)
  }
  kind <- plan$operation$operation
  argument <- plan$operation$argument
  switch(kind,
    diagonal = .sampling_covariance_diagonal(covariance),
    selected_entries = {
      if (!is.matrix(argument) || ncol(argument) != 2L) {
        stop("Selected entries require a two-column index matrix.",
          call. = FALSE)
      }
      .sampling_covariance_entries(covariance, argument[, 1L], argument[, 2L])
    },
    apply = .sampling_covariance_apply(covariance, argument),
    quadratic_form = .sampling_covariance_quadratic(covariance, argument),
    transport = .sampling_covariance_transport(covariance, argument),
    materialize = .sampling_covariance_materialize(covariance, max_bytes),
    stop("Unknown sampling-covariance operation.", call. = FALSE)
  )
}
