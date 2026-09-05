# Source calibration.R first. This repeats the full adaptive procedure on
# small matrices and is checked against .geometry_nested's public fit/score
# route. All choices use production candidate IDs and the production selector.

pgs_design <- function(arm, config) {
  design <- pgc_design(arm)
  design$sigma <- rep(1, 8)
  pairs <- t(combn(1:6, 2))
  design$train_edges <- data.frame(left = pairs[, 1], right = pairs[, 2], weight = (1:15)/sum(1:15))
  design$test_edges <- data.frame(left = 7L, right = 8L, weight = 1)
  splits <- list(list(train = 1:4, validate = 5:6), list(train = 3:6, validate = 1:2),
    list(train = c(1L, 2L, 5L, 6L), validate = 3:4))
  subset_edges <- function(parts) {
    edges <- design$train_edges[design$train_edges$left %in% parts & design$train_edges$right %in% parts, ]
    edges$weight <- edges$weight / sum(edges$weight)
    edges
  }
  design$inner <- lapply(seq_along(splits), function(i) list(
    train = subset_edges(splits[[i]]$train), validate = subset_edges(splits[[i]]$validate),
    weight = config$inner_weights[i]))
  model <- pgc_models(design)$model
  grid <- expand.grid(rank = config$ranks, penalty = config$penalties)
  design$candidates <- crossform:::.geometry_candidates(lapply(seq_len(nrow(grid)), function(i)
    list(basis = model$basis, rank = grid$rank[i], penalty = grid$penalty[i])))
  design$pool <- model$pool
  design$noise_moment <- pgc_noise_moment(design$B, design$metric, design$V, design$sigma, design$test_edges)
  design$reference_scale <- sum(design$truth^2) + design$noise_moment
  design
}

pgs_one <- function(design, parts) {
  candidates <- design$candidates; nc <- length(candidates)
  Q <- candidates[[1]]$basis$Q
  fit_candidates <- function(S) {
    penalties <- unique(vapply(candidates, `[[`, 0, "penalty"))
    paths <- lapply(penalties, function(p) crossform:::.geometry_fit_path(S, design$pool, p))
    lapply(candidates, function(candidate) {
      path <- paths[[match(candidate$penalty, penalties)]]
      crossform:::.geometry_fit_rank(path, candidate$rank)
    })
  }
  gains <- matrix(0, 1L, nc)
  for (fold in design$inner) {
    S <- crossprod(Q, pgc_dense(parts, design$metric, fold$train) %*% Q)
    test <- crossprod(Q, pgc_dense(parts, design$metric, fold$validate) %*% Q)
    fits <- fit_candidates(S)
    for (j in seq_len(nc)) {
      A <- fits[[j]]$form
      gains[1, j] <- gains[1, j] + fold$weight * (2 * sum(A * test) - sum(A^2))
    }
  }
  chosen <- crossform:::.geometry_select_indices(gains, candidates, "global", 1)
  S <- crossprod(Q, pgc_dense(parts, design$metric, design$train_edges) %*% Q)
  test <- pgc_dense(parts, design$metric, design$test_edges)
  fits <- fit_candidates(S)
  scores <- vapply(fits, function(fit) {
    F <- Q %*% fit$form %*% t(Q)
    2 * sum(F * test) - sum(F^2)
  }, 0)
  risk_gains <- vapply(fits, function(fit) {
    F <- Q %*% fit$form %*% t(Q)
    sum(design$truth^2) - sum((design$truth - F)^2)
  }, 0)
  test_selected <- which.max(scores) # explicitly invalid negative control
  c(error = unname(scores[chosen] - risk_gains[chosen]), gain = unname(scores[chosen]),
    norm_sq = fits[[chosen]]$prediction_norm_sq, target_gain = unname(risk_gains[chosen]),
    leakage_error = unname(scores[test_selected] - risk_gains[test_selected]),
    selected = chosen, selected_rank = candidates[[chosen]]$rank, selected_penalty = candidates[[chosen]]$penalty)
}

pgs_run_arm <- function(arm, count, seed, config) {
  design <- pgs_design(arm, config)
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed); stream <- .Random.seed
  result <- matrix(0, count, 8L, dimnames = list(NULL,
    c("error", "gain", "norm_sq", "target_gain", "leakage_error", "selected", "selected_rank", "selected_penalty")))
  for (i in seq_len(count)) {
    result[i, ] <- pgs_one(design, pgc_parts(design, stream))
    stream <- parallel::nextRNGStream(stream)
    if (i %% 1000L == 0L) { cat(arm, i, "/", count, "nested datasets\n"); flush.console() }
  }
  list(design = design, count = count, seed = seed, values = result)
}

pgs_summary <- function(result, config) {
  estimate <- mean(result$values[, "error"]); mcse <- sd(result$values[, "error"])/sqrt(result$count)
  b <- result$design$reference_scale; delta <- config$relative_equivalence_margin * b
  data.frame(arm = result$design$arm, n = result$count, error = estimate, normalized_error = estimate/b,
    mcse = mcse, reference_scale = b, margin = delta,
    bias_pass = abs(estimate) <= config$confidence_multiplier * mcse + config$numerical_tolerance,
    precision_pass = abs(estimate) + config$confidence_multiplier * mcse <= delta,
    leakage_bias = mean(result$values[, "leakage_error"]),
    leakage_mcse = sd(result$values[, "leakage_error"])/sqrt(result$count),
    mean_gain = mean(result$values[, "gain"]))
}
