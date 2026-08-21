# Tiny dense oracle for the *future* cross-node covariance contract.
#
# This file defines algebra and validation fixtures only. It is not sourced by
# the package, exports no inference function, and must not be cited as evidence
# that cross-node covariance or precision weighting is implemented.

future_cross_node_dense <- function(covariance, tolerance = sqrt(.Machine$double.eps),
                                    max_dense_bytes = 1024^2) {
  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      nrow(covariance) != ncol(covariance) || any(!is.finite(covariance))) {
    stop("future covariance must be one finite square numeric matrix")
  }
  bytes <- 8 * length(covariance)
  if (bytes > max_dense_bytes) {
    stop("future dense covariance exceeds its declared compute budget")
  }
  scale <- max(1, max(abs(covariance)))
  if (max(abs(covariance - t(covariance))) > tolerance * scale) {
    stop("future covariance is not symmetric within tolerance")
  }
  covariance <- (covariance + t(covariance)) / 2
  minimum <- min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values)
  if (minimum < -tolerance * scale) {
    stop("future covariance is not positive semidefinite within tolerance")
  }
  covariance
}

future_cross_node_transport <- function(operator, covariance) {
  if (!is.matrix(operator) || !is.numeric(operator) ||
      nrow(operator) != nrow(covariance) ||
      max(abs(rowSums(operator) - 1)) > 1e-12) {
    stop("future transport must be numeric, aligned, and row stochastic")
  }
  crossprod(operator, covariance %*% operator)
}

cross_node_covariance_oracle <- function() {
  covariance <- matrix(c(
    1.00, 0.35, 0.20,
    0.35, 1.50, 0.40,
    0.20, 0.40, 0.80
  ), 3L, 3L, byrow = TRUE)
  covariance <- future_cross_node_dense(covariance)
  operator <- rbind(
    c(1.0, 0.0, 0.0),
    c(0.4, 0.6, 0.0),
    c(0.0, 0.7, 0.3)
  )
  transported <- future_cross_node_transport(operator, covariance)
  full_budget_variance <- as.numeric(crossprod(rep(1, 3L),
    covariance %*% rep(1, 3L)))
  transported_budget_variance <- as.numeric(crossprod(rep(1, 3L),
    transported %*% rep(1, 3L)))
  diagonal_shortcut <- sum(diag(covariance))

  sparse <- Matrix::Matrix(covariance, sparse = TRUE)
  sparse_roundtrip <- as.matrix(sparse)
  list(
    covariance = covariance,
    operator = operator,
    transported = transported,
    sparse = sparse,
    sparse_roundtrip = sparse_roundtrip,
    full_budget_variance = full_budget_variance,
    transported_budget_variance = transported_budget_variance,
    diagonal_shortcut = diagonal_shortcut,
    psd_minimum = min(eigen(covariance, symmetric = TRUE,
      only.values = TRUE)$values)
  )
}

if (sys.nframe() == 0L) {
  oracle <- cross_node_covariance_oracle()
  stopifnot(
    oracle$psd_minimum >= 0,
    identical(oracle$sparse_roundtrip, oracle$covariance),
    abs(oracle$transported_budget_variance -
      oracle$full_budget_variance) < 1e-12,
    abs(oracle$diagonal_shortcut - oracle$full_budget_variance) > 0.5
  )
  cat("cross-node-covariance-future-oracle PASS\n")
}
