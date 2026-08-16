# Independent numerical laws for evidence-sampling-v1. These helpers use base
# R only and deliberately call no crossform function.

sampling_oracle_matrix <- function(value, name, symmetric = FALSE) {
  if (!is.matrix(value) || !is.numeric(value) || any(dim(value) < 1L) ||
      any(!is.finite(value))) {
    stop(name, " must be one finite nonempty numeric matrix")
  }
  if (isTRUE(symmetric) &&
      (nrow(value) != ncol(value) ||
       max(abs(value - t(value))) > 1e-10)) {
    stop(name, " must be symmetric")
  }
  value
}

sampling_oracle_condition_contrasts <- function(conditions) {
  if (!is.numeric(conditions) || length(conditions) != 1L ||
      is.na(conditions) || !is.finite(conditions) ||
      conditions %% 1 != 0 || conditions < 2L) {
    stop("conditions must be one integer of at least two")
  }
  conditions <- as.integer(conditions)
  pairs <- utils::combn(conditions, 2L)
  contrasts <- matrix(0, ncol(pairs), conditions)
  for (distance in seq_len(ncol(pairs))) {
    contrasts[distance, pairs[1L, distance]] <- 1
    contrasts[distance, pairs[2L, distance]] <- -1
  }
  rownames(contrasts) <- paste0(pairs[1L, ], "-", pairs[2L, ])
  colnames(contrasts) <- paste0("condition", seq_len(conditions))
  contrasts
}

sampling_oracle_all_pairs <- function(partition_differences) {
  partition_differences <- sampling_oracle_matrix(
    partition_differences, "partition differences"
  )
  partitions <- nrow(partition_differences)
  features <- ncol(partition_differences)
  if (partitions < 2L) stop("at least two partitions are required")
  edges <- utils::combn(partitions, 2L)
  values <- vapply(seq_len(ncol(edges)), function(edge) {
    sum(partition_differences[edges[1L, edge], ] *
      partition_differences[edges[2L, edge], ]) / features
  }, numeric(1))
  mean(values)
}

sampling_oracle_leave_one_out <- function(partition_differences) {
  partition_differences <- sampling_oracle_matrix(
    partition_differences, "partition differences"
  )
  partitions <- nrow(partition_differences)
  features <- ncol(partition_differences)
  if (partitions < 2L) stop("at least two partitions are required")
  values <- vapply(seq_len(partitions), function(partition) {
    other_mean <- colMeans(partition_differences[-partition, , drop = FALSE])
    sum(partition_differences[partition, ] * other_mean) / features
  }, numeric(1))
  mean(values)
}

sampling_oracle_partitions <- function(partitions, minimum = 2L) {
  if (!is.numeric(partitions) || length(partitions) != 1L ||
      is.na(partitions) || !is.finite(partitions) ||
      partitions %% 1 != 0 || partitions < minimum) {
    stop("partitions must be one integer of at least ", minimum)
  }
  as.integer(partitions)
}

sampling_oracle_normalization <- function(normalization) {
  if (!is.numeric(normalization) || length(normalization) != 1L ||
      is.na(normalization) || !is.finite(normalization) ||
      normalization <= 0) {
    stop("normalization must be one positive finite value")
  }
  as.double(normalization)
}

# The exact sampling law of the equal-weight all-unordered-pairs
# crossvalidated distance estimator, given
#
#   mu_r    row r of `differences`, the whitened contrast pattern;
#   Sigma_R `sigma_r`, the whitened residual covariance;
#   Xi      `xi`, the effect-coordinate contrast cross-products;
#   M       `partitions`; and
#   nu      `normalization`, the divisor already applied to each distance:
#
#   Cov(d_r, d_s) = (4 / M) Xi_rs (mu_r Sigma_R mu_s') / nu^2
#                 + 2 / (M (M - 1)) Xi_rs^2 tr(Sigma_R Sigma_R) / nu^2.
#
# The signal term uses Sigma_R as a metric on the whitened patterns. Only
# when Sigma_R = (tr(Sigma_R^2) / nu) I does it collapse to a multiple of
# the plain Gram mu_r mu_s'.
sampling_oracle_eq13_terms <- function(signal_gram, xi, noise_trace,
                                       partitions) {
  signal_gram <- sampling_oracle_matrix(
    signal_gram, "signal Gram", symmetric = TRUE
  )
  xi <- sampling_oracle_matrix(xi, "Xi", symmetric = TRUE)
  if (!identical(dim(signal_gram), dim(xi))) {
    stop("signal Gram and Xi dimensions must agree")
  }
  if (!is.numeric(noise_trace) || length(noise_trace) != 1L ||
      is.na(noise_trace) || !is.finite(noise_trace) || noise_trace < 0) {
    stop("noise trace must be one finite nonnegative value")
  }
  partitions <- sampling_oracle_partitions(partitions)
  signal <- 4 * (signal_gram * xi) / partitions
  noise <- 2 * (xi * xi) * noise_trace /
    (partitions * (partitions - 1L))
  list(signal = signal, noise = noise, covariance = signal + noise)
}

sampling_oracle_eq13 <- function(differences, xi, sigma_r, partitions,
                                 normalization = ncol(differences)) {
  differences <- sampling_oracle_matrix(
    differences, "whitened contrast differences"
  )
  sigma_r <- sampling_oracle_matrix(
    sigma_r, "residual covariance", symmetric = TRUE
  )
  if (ncol(differences) != nrow(sigma_r)) {
    stop("difference and residual-covariance feature axes must agree")
  }
  normalization <- sampling_oracle_normalization(normalization)
  sampling_oracle_eq13_terms(
    signal_gram = differences %*% sigma_r %*% t(differences) /
      normalization^2,
    xi = xi,
    noise_trace = sum(sigma_r * sigma_r) / normalization^2,
    partitions = partitions
  )
}

sampling_oracle_endpoint_enumeration <- function(differences, xi, sigma_r,
                                                  partitions,
                                                  normalization =
                                                    ncol(differences)) {
  differences <- sampling_oracle_matrix(
    differences, "whitened contrast differences"
  )
  xi <- sampling_oracle_matrix(xi, "Xi", symmetric = TRUE)
  sigma_r <- sampling_oracle_matrix(
    sigma_r, "residual covariance", symmetric = TRUE
  )
  if (ncol(differences) != nrow(sigma_r)) {
    stop("difference and residual-covariance feature axes must agree")
  }
  if (!identical(dim(xi), as.integer(rep(nrow(differences), 2L)))) {
    stop("Xi must be one distance-by-distance matrix")
  }
  partitions <- sampling_oracle_partitions(partitions)
  normalization <- sampling_oracle_normalization(normalization)
  edges <- t(utils::combn(partitions, 2L))
  edge_count <- nrow(edges)
  shared_endpoints <- 0
  identical_edges <- 0
  for (first in seq_len(edge_count)) {
    for (second in seq_len(edge_count)) {
      shared_endpoints <- shared_endpoints +
        length(intersect(edges[first, ], edges[second, ]))
      identical_edges <- identical_edges + as.integer(first == second)
    }
  }
  signal_gram <- differences %*% sigma_r %*% t(differences) /
    normalization^2
  noise_trace <- sum(sigma_r * sigma_r) / normalization^2
  signal_coefficient <- shared_endpoints / edge_count^2
  noise_coefficient <- identical_edges / edge_count^2
  list(
    shared_endpoints = shared_endpoints,
    identical_edges = identical_edges,
    signal_coefficient = signal_coefficient,
    noise_coefficient = noise_coefficient,
    covariance = signal_coefficient * signal_gram * xi +
      noise_coefficient * xi * xi * noise_trace
  )
}

sampling_oracle_components <- function(patterns, sigma_k, sigma_r) {
  patterns <- sampling_oracle_matrix(patterns, "true patterns")
  sigma_k <- sampling_oracle_matrix(
    sigma_k, "condition covariance", symmetric = TRUE
  )
  sigma_r <- sampling_oracle_matrix(
    sigma_r, "residual covariance", symmetric = TRUE
  )
  conditions <- nrow(patterns)
  features <- ncol(patterns)
  if (!identical(dim(sigma_k), c(conditions, conditions)) ||
      !identical(dim(sigma_r), c(features, features))) {
    stop("pattern and covariance dimensions do not agree")
  }
  contrasts <- sampling_oracle_condition_contrasts(conditions)
  differences <- contrasts %*% patterns
  list(
    contrasts = contrasts,
    differences = differences,
    delta = tcrossprod(differences) / features,
    xi = contrasts %*% sigma_k %*% t(contrasts)
  )
}

# The same law written entry by entry with explicit scalar loops, sharing no
# code path with `sampling_oracle_eq13()`. Slow on purpose: it exists so an
# exact algebraic comparison has somewhere independent to stand.
sampling_oracle_scalar_law <- function(differences, xi, sigma_r, partitions,
                                       normalization) {
  differences <- sampling_oracle_matrix(
    differences, "whitened contrast differences"
  )
  xi <- sampling_oracle_matrix(xi, "Xi", symmetric = TRUE)
  sigma_r <- sampling_oracle_matrix(
    sigma_r, "residual covariance", symmetric = TRUE
  )
  if (ncol(differences) != nrow(sigma_r) ||
      !identical(dim(xi), as.integer(rep(nrow(differences), 2L)))) {
    stop("scalar law axes do not agree")
  }
  partitions <- sampling_oracle_partitions(partitions)
  normalization <- sampling_oracle_normalization(normalization)
  features <- ncol(differences)
  dimension <- nrow(differences)
  noise_trace <- 0
  for (i in seq_len(features)) {
    for (j in seq_len(features)) {
      noise_trace <- noise_trace + sigma_r[i, j] * sigma_r[j, i]
    }
  }
  value <- matrix(0, dimension, dimension)
  for (row in seq_len(dimension)) {
    for (column in seq_len(dimension)) {
      quadratic <- 0
      for (i in seq_len(features)) {
        for (j in seq_len(features)) {
          quadratic <- quadratic +
            differences[row, i] * sigma_r[i, j] * differences[column, j]
        }
      }
      value[row, column] <-
        4 * xi[row, column] * quadratic /
          (partitions * normalization^2) +
        2 * xi[row, column]^2 * noise_trace /
          (partitions * (partitions - 1) * normalization^2)
    }
  }
  value
}

sampling_oracle_psd_root <- function(value, name) {
  value <- sampling_oracle_matrix(value, name, symmetric = TRUE)
  spectrum <- eigen(value, symmetric = TRUE)
  scale <- max(1, max(abs(spectrum$values)))
  if (min(spectrum$values) < -1e-10 * scale) {
    stop(name, " must be positive semidefinite")
  }
  retained <- spectrum$values > 1e-12 * scale
  if (!any(retained)) return(matrix(0, nrow(value), 1L))
  spectrum$vectors[, retained, drop = FALSE] %*%
    diag(sqrt(spectrum$values[retained]), sum(retained))
}

# Draw `replications` independent realizations of the equal-weight
# all-unordered-pairs distance estimator directly in distance coordinates.
# Contrast rows of a matrix-normal partition estimate are themselves
# matrix-normal: C E ~ MN(0, C Sigma_K C', Sigma_R) = MN(0, Xi, Sigma_R), so
# no condition-space or design-space simulation is needed and no crossform
# code is involved.
sampling_oracle_distance_draws <- function(differences, xi, sigma_r,
                                           partitions, normalization,
                                           replications, seed) {
  differences <- sampling_oracle_matrix(
    differences, "whitened contrast differences"
  )
  sigma_r <- sampling_oracle_matrix(
    sigma_r, "residual covariance", symmetric = TRUE
  )
  if (ncol(differences) != nrow(sigma_r)) {
    stop("difference and residual-covariance feature axes must agree")
  }
  dimension <- nrow(differences)
  features <- ncol(differences)
  if (!identical(dim(xi), as.integer(rep(dimension, 2L)))) {
    stop("Xi must be one distance-by-distance matrix")
  }
  partitions <- sampling_oracle_partitions(partitions)
  normalization <- sampling_oracle_normalization(normalization)
  if (!is.numeric(replications) || length(replications) != 1L ||
      is.na(replications) || !is.finite(replications) ||
      replications %% 1 != 0 || replications < 2L) {
    stop("replications must be one integer of at least two")
  }
  replications <- as.integer(replications)
  root_xi <- sampling_oracle_psd_root(xi, "Xi")
  root_r <- sampling_oracle_psd_root(sigma_r, "residual covariance")
  rank_xi <- ncol(root_xi)
  rank_r <- ncol(root_r)
  set.seed(seed)
  total <- vector("list", features)
  for (feature in seq_len(features)) {
    total[[feature]] <- matrix(0, replications, dimension)
  }
  squares <- matrix(0, replications, dimension)
  for (partition in seq_len(partitions)) {
    noise <- array(
      matrix(
        rnorm(replications * rank_xi * rank_r),
        replications * rank_xi, rank_r
      ) %*% t(root_r),
      c(replications, rank_xi, features)
    )
    for (feature in seq_len(features)) {
      draw <- matrix(noise[, , feature], replications, rank_xi) %*%
        t(root_xi) + rep(differences[, feature], each = replications)
      total[[feature]] <- total[[feature]] + draw
      squares <- squares + draw * draw
    }
  }
  # sum_{a<b} <g_a, g_b> = (||sum_m g_m||^2 - sum_m ||g_m||^2) / 2.
  crossed <- matrix(0, replications, dimension)
  for (feature in seq_len(features)) {
    crossed <- crossed + total[[feature]] * total[[feature]]
  }
  edges <- choose(partitions, 2L)
  values <- (crossed - squares) / (2 * edges * normalization)
  colnames(values) <- rownames(differences)
  values
}

sampling_oracle_naive_signal_factors <- function(partitions) {
  if (!is.numeric(partitions) || length(partitions) != 1L ||
      is.na(partitions) || !is.finite(partitions) ||
      partitions %% 1 != 0 || partitions < 3L) {
    stop("partitions must be one integer of at least three")
  }
  partitions <- as.integer(partitions)
  edges <- t(utils::combn(partitions, 2L))
  edge_count <- nrow(edges)
  shared_distinct <- 0L
  for (first in seq_len(edge_count)) {
    for (second in seq_len(edge_count)) {
      if (first != second) {
        shared_distinct <- shared_distinct +
          length(intersect(edges[first, ], edges[second, ]))
      }
    }
  }
  marginal_signal_variance <- 2
  mean_off_diagonal_covariance <- shared_distinct /
    (edge_count * (edge_count - 1L))
  expected_sample_variance <- marginal_signal_variance -
    mean_off_diagonal_covariance
  true_mean_variance <- 4 / partitions
  sample_sd_naive_variance <- expected_sample_variance / edge_count
  marginal_edge_naive_variance <- marginal_signal_variance / edge_count
  c(
    sample_sd = sqrt(true_mean_variance / sample_sd_naive_variance),
    marginal_edge = sqrt(true_mean_variance /
      marginal_edge_naive_variance)
  )
}

sampling_oracle_matrix_normal_draw <- function(mean, sigma_k, sigma_r) {
  mean <- sampling_oracle_matrix(mean, "matrix-normal mean")
  sigma_k <- sampling_oracle_matrix(
    sigma_k, "effect-coordinate covariance", symmetric = TRUE
  )
  sigma_r <- sampling_oracle_matrix(
    sigma_r, "neural covariance", symmetric = TRUE
  )
  if (!identical(dim(sigma_k), c(nrow(mean), nrow(mean))) ||
      !identical(dim(sigma_r), c(ncol(mean), ncol(mean)))) {
    stop("matrix-normal covariance dimensions do not match the mean")
  }
  left <- t(chol(sigma_k))
  right <- chol(sigma_r)
  mean + left %*% matrix(rnorm(length(mean)), nrow(mean), ncol(mean)) %*%
    right
}

sampling_oracle_partition_edge_distances <- function(patterns, contrasts,
                                                       normalization) {
  if (!is.list(patterns) || length(patterns) < 2L) {
    stop("at least two partition pattern matrices are required")
  }
  reference <- sampling_oracle_matrix(patterns[[1L]], "partition patterns")
  if (!all(vapply(patterns, function(value) {
    is.matrix(value) && is.numeric(value) && identical(dim(value), dim(reference)) &&
      all(is.finite(value))
  }, logical(1)))) {
    stop("partition pattern matrices must share one finite shape")
  }
  contrasts <- sampling_oracle_matrix(contrasts, "distance contrasts")
  if (ncol(contrasts) != nrow(reference) || !is.numeric(normalization) ||
      length(normalization) != 1L || is.na(normalization) ||
      !is.finite(normalization) || normalization <= 0) {
    stop("distance contrasts or normalization are invalid")
  }
  differences <- lapply(patterns, function(value) contrasts %*% value)
  edges <- t(utils::combn(length(patterns), 2L))
  value <- vapply(seq_len(nrow(edges)), function(edge) {
    rowSums(
      differences[[edges[edge, 1L]]] * differences[[edges[edge, 2L]]]
    ) / normalization
  }, numeric(nrow(contrasts)))
  if (nrow(contrasts) == 1L) value <- matrix(value, nrow = 1L)
  dimnames(value) <- list(rownames(contrasts), NULL)
  value
}

sampling_oracle_calibration_experiment <- function(
    true_patterns, sigma_k, sigma_r, partitions, repetitions,
    seed, plugin = TRUE) {
  true_patterns <- sampling_oracle_matrix(true_patterns, "true patterns")
  sigma_k <- sampling_oracle_matrix(
    sigma_k, "effect-coordinate covariance", symmetric = TRUE
  )
  sigma_r <- sampling_oracle_matrix(
    sigma_r, "neural covariance", symmetric = TRUE
  )
  if (!is.numeric(partitions) || length(partitions) != 1L ||
      is.na(partitions) || !is.finite(partitions) ||
      partitions %% 1 != 0 || partitions < 2L ||
      !is.numeric(repetitions) || length(repetitions) != 1L ||
      is.na(repetitions) || !is.finite(repetitions) ||
      repetitions %% 1 != 0 || repetitions < 2L ||
      !is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed)) {
    stop("experiment counts and seed are invalid")
  }
  partitions <- as.integer(partitions)
  repetitions <- as.integer(repetitions)
  contrasts <- sampling_oracle_condition_contrasts(nrow(true_patterns))
  dimension <- nrow(contrasts)
  edge_count <- choose(partitions, 2L)
  normalization <- ncol(true_patterns)
  estimates <- matrix(NA_real_, repetitions, dimension)
  naive_se <- matrix(NA_real_, repetitions, dimension)
  plugin_se <- if (isTRUE(plugin)) {
    matrix(NA_real_, repetitions, dimension)
  } else {
    NULL
  }
  set.seed(seed)
  for (replicate in seq_len(repetitions)) {
    draws <- lapply(seq_len(partitions), function(index) {
      sampling_oracle_matrix_normal_draw(true_patterns, sigma_k, sigma_r)
    })
    edge_values <- sampling_oracle_partition_edge_distances(
      draws, contrasts, normalization
    )
    estimates[replicate, ] <- rowMeans(edge_values)
    naive_se[replicate, ] <- apply(edge_values, 1L, stats::sd) /
      sqrt(edge_count)
    if (isTRUE(plugin)) {
      mean_patterns <- Reduce(`+`, draws) / partitions
      components <- sampling_oracle_components(
        mean_patterns, sigma_k, sigma_r
      )
      plugin_se[replicate, ] <- sqrt(diag(sampling_oracle_eq13(
        components$differences, components$xi, sigma_r, partitions,
        normalization
      )$covariance))
    }
  }
  components <- sampling_oracle_components(
    true_patterns, sigma_k, sigma_r
  )
  covariance <- sampling_oracle_eq13(
    components$differences, components$xi, sigma_r, partitions,
    normalization
  )$covariance
  list(
    estimates = estimates,
    truth = diag(components$delta),
    covariance = covariance,
    plugin_se = plugin_se,
    naive_se = naive_se,
    contrasts = contrasts,
    partitions = partitions,
    repetitions = repetitions
  )
}
