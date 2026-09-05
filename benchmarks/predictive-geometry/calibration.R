# Shared scientific generator and independent dense risk oracle. The fitting
# calls use production numerical routines; dense targets/moments do not.

pgc_design <- function(arm) {
  n <- 6L; p <- 5L
  U <- qr.Q(qr(cbind(rep(1, n), diag(n))))[, -1L, drop = FALSE]
  Q <- U[, 1:2, drop = FALSE]
  dimnames(Q) <- list(paste0("condition", seq_len(n)), NULL)
  selected <- if (arm == "outside_model") c(3L, 4L) else 1:2
  B <- U[, selected, drop = FALSE] %*% diag(c(1.5, .8)) %*% cbind(diag(2), matrix(0, 2, p - 2))
  if (arm == "null") B[] <- 0
  rownames(B) <- rownames(Q)
  metric <- if (arm == "heterogeneous_metric") diag(seq(.7, 1.5, length.out = p)) + matrix(.1, p, p) else diag(p)
  V <- if (arm == "heterogeneous_metric") diag(c(.25, .5, 1, 1.5, 2)) else diag(p)
  sigma <- if (arm == "heterogeneous_metric") c(.4, .8, 1.3, 1.1, .5, 1.4) else rep(1, 6)
  train_edges <- data.frame(left = c(1L, 1L, 2L), right = c(2L, 3L, 3L), weight = c(.2, .3, .5))
  test_edges <- data.frame(left = c(4L, 4L, 5L), right = c(5L, 6L, 6L), weight = c(.5, .2, .3))
  truth <- B %*% metric %*% t(B)
  variance <- pgc_noise_moment(B, metric, V, sigma, test_edges)
  list(arm = arm, n = n, p = p, Q = Q, B = B, metric = metric,
    V = V, noise_factor = t(chol(V)), sigma = sigma,
    train_edges = train_edges, test_edges = test_edges,
    truth = truth, noise_moment = variance, reference_scale = sum(truth^2) + variance)
}

pgc_dense <- function(parts, metric, edges) {
  result <- matrix(0, nrow(parts[[1]]), nrow(parts[[1]]))
  for (i in seq_len(nrow(edges))) {
    product <- parts[[edges$left[i]]] %*% metric %*% t(parts[[edges$right[i]]])
    result <- result + edges$weight[i] * (product + t(product)) / 2
  }
  result
}

pgc_noise_moment <- function(B, metric, V, sigma, edges) {
  # E_a has independent condition rows, each with covariance sigma_a^2 V.
  # Distinct quadratic edges have zero covariance, even when sharing one
  # endpoint. Linear terms collect incident weights before squaring them.
  degrees <- numeric(length(sigma))
  for (i in seq_len(nrow(edges))) {
    degrees[edges$left[i]] <- degrees[edges$left[i]] + edges$weight[i]
    degrees[edges$right[i]] <- degrees[edges$right[i]] + edges$weight[i]
  }
  n <- nrow(B)
  linear <- (n + 1) / 2 * sum(diag(B %*% metric %*% V %*% metric %*% t(B))) * sum(degrees^2 * sigma^2)
  quadratic <- n * (n + 1) / 2 * sum(diag(metric %*% V %*% metric %*% V)) *
    sum(edges$weight^2 * sigma[edges$left]^2 * sigma[edges$right]^2)
  linear + quadratic
}

pgc_parts <- function(design, stream) {
  assign(".Random.seed", stream, envir = .GlobalEnv)
  lapply(seq_along(design$sigma), function(i) {
    design$B + design$sigma[i] * matrix(rnorm(design$n * design$p), design$n) %*% t(design$noise_factor)
  })
}

pgc_models <- function(design) {
  eig <- list(model = c(.8, .2), isotropic = c(.5, .5), mismatched = c(.2, .8))
  lapply(eig, function(d) {
    K <- design$Q %*% diag(d) %*% t(design$Q)
    dimnames(K) <- list(rownames(design$Q), rownames(design$Q))
    basis <- crossform::model_basis(kernels = list(model = K), normalize = "trace")
    list(basis = basis, pool = crossform:::.model_basis_pool(basis, NULL))
  })
}

pgc_one <- function(design, models, parts, rank, penalty) {
  train <- pgc_dense(parts, design$metric, design$train_edges)
  test <- pgc_dense(parts, design$metric, design$test_edges)
  lapply(models, function(model) {
    Q <- model$basis$Q
    path <- crossform:::.geometry_fit_path(crossprod(Q, train %*% Q), model$pool, penalty)
    fit <- crossform:::.geometry_fit_rank(path, rank)
    prediction <- list(factor = fit$factor, amplitudes = fit$amplitudes, shifted_spectrum = path$shifted_spectrum)
    score <- crossform:::.geometry_score_row(prediction, crossprod(Q, test %*% Q), model$basis$tolerance)
    F <- Q %*% fit$form %*% t(Q)
    risk_gain <- sum(design$truth^2) - sum((design$truth - F)^2)
    leakage <- 2 * sum(F * train) - sum(F^2) - risk_gain
    c(error = score$gain - risk_gain, gain = score$gain,
      norm_sq = sum(F^2), target_gain = risk_gain, projector_energy = sum(score$evidence),
      noise_moment = sum((test - design$truth)^2), leakage_error = leakage)
  })
}

pgc_run_arm <- function(arm, count, seed, config, progress = TRUE) {
  design <- pgc_design(arm); models <- pgc_models(design)
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed); stream <- .Random.seed
  results <- array(0, c(count, 7L, length(models)), dimnames = list(NULL,
    c("error", "gain", "norm_sq", "target_gain", "projector_energy", "noise_moment", "leakage_error"), names(models)))
  for (i in seq_len(count)) {
    parts <- pgc_parts(design, stream)
    result <- pgc_one(design, models, parts, config$rank, config$penalty)
    for (j in seq_along(models)) results[i, , j] <- result[[j]]
    stream <- parallel::nextRNGStream(stream)
    if (progress && i %% 5000L == 0L) { cat(arm, i, "/", count, "datasets\n"); flush.console() }
  }
  list(design = design, count = count, seed = seed, values = results)
}

pgc_summary <- function(result, config) {
  b <- result$design$reference_scale
  do.call(rbind, lapply(dimnames(result$values)[[3L]], function(predictor) {
    values <- result$values[, , predictor]
    error <- values[, "error"]; estimate <- mean(error); mcse <- sd(error) / sqrt(nrow(values))
    delta <- config$relative_equivalence_margin * b
    data.frame(arm = result$design$arm, predictor = predictor, n = nrow(values),
      reference_scale = b, error = estimate, normalized_error = estimate / b, mcse = mcse,
      margin = delta, bias_pass = abs(estimate) <= config$confidence_multiplier * mcse + config$numerical_tolerance,
      precision_pass = abs(estimate) + config$confidence_multiplier * mcse <= delta,
      mean_gain = mean(values[, "gain"]), mean_prediction_norm_sq = mean(values[, "norm_sq"]),
      mean_projector_energy = mean(values[, "projector_energy"]),
      mean_leakage_error = mean(values[, "leakage_error"]),
      leakage_mcse = sd(values[, "leakage_error"]) / sqrt(nrow(values)))
  }))
}
