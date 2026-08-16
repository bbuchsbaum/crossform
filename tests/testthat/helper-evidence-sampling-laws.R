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

sampling_oracle_eq13 <- function(delta, xi, sigma_r, partitions) {
  delta <- sampling_oracle_matrix(delta, "Delta", symmetric = TRUE)
  xi <- sampling_oracle_matrix(xi, "Xi", symmetric = TRUE)
  sigma_r <- sampling_oracle_matrix(
    sigma_r, "residual covariance", symmetric = TRUE
  )
  if (!identical(dim(delta), dim(xi))) {
    stop("Delta and Xi dimensions must agree")
  }
  if (!is.numeric(partitions) || length(partitions) != 1L ||
      is.na(partitions) || !is.finite(partitions) ||
      partitions %% 1 != 0 || partitions < 2L) {
    stop("partitions must be one integer of at least two")
  }
  partitions <- as.integer(partitions)
  features <- nrow(sigma_r)
  residual_factor <- sum(sigma_r * sigma_r) / features^2
  signal <- 4 * (delta * xi) / partitions * residual_factor
  noise <- 2 * (xi * xi) /
    (partitions * (partitions - 1L)) * residual_factor
  list(signal = signal, noise = noise, covariance = signal + noise)
}

sampling_oracle_endpoint_enumeration <- function(delta, xi, sigma_r,
                                                  partitions) {
  delta <- sampling_oracle_matrix(delta, "Delta", symmetric = TRUE)
  xi <- sampling_oracle_matrix(xi, "Xi", symmetric = TRUE)
  sigma_r <- sampling_oracle_matrix(
    sigma_r, "residual covariance", symmetric = TRUE
  )
  if (!identical(dim(delta), dim(xi))) {
    stop("Delta and Xi dimensions must agree")
  }
  if (!is.numeric(partitions) || length(partitions) != 1L ||
      is.na(partitions) || !is.finite(partitions) ||
      partitions %% 1 != 0 || partitions < 2L) {
    stop("partitions must be one integer of at least two")
  }
  partitions <- as.integer(partitions)
  features <- nrow(sigma_r)
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
  residual_factor <- sum(sigma_r * sigma_r) / features^2
  signal_coefficient <- shared_endpoints / edge_count^2
  noise_coefficient <- identical_edges / edge_count^2
  list(
    shared_endpoints = shared_endpoints,
    identical_edges = identical_edges,
    signal_coefficient = signal_coefficient,
    noise_coefficient = noise_coefficient,
    covariance = (
      signal_coefficient * delta * xi +
      noise_coefficient * xi * xi
    ) * residual_factor
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
        components$delta, components$xi, sigma_r, partitions
      )$covariance))
    }
  }
  components <- sampling_oracle_components(
    true_patterns, sigma_k, sigma_r
  )
  covariance <- sampling_oracle_eq13(
    components$delta, components$xi, sigma_r, partitions
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
