residual_statistics_fixture <- function() {
  set.seed(8301)
  observations <- 37L
  features <- 19L
  coordinates <- cbind(x = seq_len(features), y = 0, z = 0)
  domain <- abstract_domain(
    features, coordinates = coordinates, id = "residual-pair-fixture",
    coordinate_units = "mm"
  )
  frame <- compile_frame(
    searchlights(2.01, normalization = "local"), domain
  )
  design <- cbind(
    intercept = 1,
    condition = rep(c(-0.5, 0.5), length.out = observations),
    drift = seq(-1, 1, length.out = observations)
  )
  effects <- rbind(
    condition = c(0, 1, 0),
    drift = c(0, 0, 1)
  )
  latent <- matrix(rnorm(observations * 4L), observations, 4L)
  loading <- matrix(rnorm(4L * features), 4L, features)
  response <- lapply(seq_len(3L), function(partition) {
    latent %*% loading +
      matrix(rnorm(observations * features, sd = 1e-3),
        observations, features) + partition / 10
  })
  names(response) <- paste0("run", seq_along(response))
  reads <- stats::setNames(integer(length(response)), names(response))
  sources <- lapply(names(response), function(partition) {
    force(partition)
    function(features) {
      reads[[partition]] <<- reads[[partition]] + 1L
      response[[partition]][, features, drop = FALSE]
    }
  })
  names(sources) <- names(response)
  revisions <- lapply(seq_along(response), function(index) {
    paste0("sha256:", paste(rep(letters[[index]], 64L), collapse = ""))
  })
  capabilities <- lapply(revisions, function(revision) {
    source_capabilities(block_read = TRUE, stable_revision = revision)
  })
  fit <- lm_relation_fit(
    sources, design, effects,
    source_dims = rep(list(c(observations, features)), length(response)),
    capabilities = capabilities, domain = domain
  )
  reads[] <- 0L
  list(
    fit = fit,
    frame = frame,
    domain = domain,
    response = response,
    reads = function() reads,
    reset_reads = function() reads[] <<- 0L
  )
}

residual_statistics_budgets <- function(fixture) {
  index <- fixture$frame$support_index
  pairs <- effectagram:::.residual_pair_coordinates(index)
  width <- effectagram:::.residual_tile_width(index)
  observations <- vapply(
    fixture$fit$error_models, function(model) model$residual_source$dim[[1L]],
    integer(1)
  )
  generous <- effectagram:::.residual_pair_memory_plan(
    length(pairs$i), length(observations), observations, width,
    ceiling(fixture$fit$relation$n_features / width), 2^40
  )
  list(
    minimum = generous$minimum_workspace_bytes,
    wider = generous$minimum_workspace_bytes +
      generous$residual_tile_bytes * 8
  )
}

oracle_local_residual_covariance <- function(fit, partitions, support) {
  cross_products <- lapply(partitions, function(partition) {
    crossprod(residual_block(fit, partition, support))
  })
  df <- sum(vapply(partitions, function(partition) {
    residual_df(fit, partition)
  }, integer(1)))
  Reduce(`+`, cross_products) / df
}

oracle_shrinkage_covariance <- function(raw, recipe) {
  variance <- diag(raw)
  positive <- variance[variance > 0]
  floor <- max(
    recipe$hyperparameters$absolute_variance_floor,
    recipe$hyperparameters$relative_variance_floor * mean(positive)
  )
  target_variance <- pmax(variance, floor)
  alpha <- recipe$hyperparameters$shrinkage
  value <- (1 - alpha) * raw + alpha * diag(target_variance, nrow(raw))
  spectrum <- eigen(value, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(abs(spectrum), max(abs(diag(value))), .Machine$double.xmin)
  ridge <- max(0,
    recipe$hyperparameters$relative_spectral_floor * scale - min(spectrum))
  if (ridge > 0) value <- value + diag(ridge, nrow(value))
  value
}

metric_learning_setup <- function(workspace = c("wide", "narrow")) {
  workspace <- match.arg(workspace)
  fixture <- residual_statistics_fixture()
  budgets <- residual_statistics_budgets(fixture)
  statistics <- residual_pair_statistics(
    fixture$fit, fixture$frame,
    workspace_bytes = budgets[[if (workspace == "wide") "wider" else
      "minimum"]]
  )
  fixture$reset_reads()
  over <- pairing("run1", "run2", independence = "independent")
  list(fixture = fixture, statistics = statistics, over = over,
    budgets = budgets)
}
