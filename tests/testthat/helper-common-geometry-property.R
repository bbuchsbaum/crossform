cg_oracle <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    cached <<- new.env(parent = globalenv())
    sys.source(test_path("..", "..", "design", "oracles",
      "common-geometry-equivalence.R"), envir = cached)
    cached
  }
})

cg_plan <- function(blocks, metric, id) {
  features <- ncol(blocks[[1L]])
  domain <- abstract_domain(features, id = id)
  relation <- relation(blocks, domain = domain)
  plan_geometry(
    relation,
    compile_frame(whole_brain(), domain),
    cross_partitions(relation, independence = "independent"),
    metric = noise_precision(metric, domain,
      provenance = list(source = "seeded common-geometry court"))
  )
}

cg_rdm_matrix <- function(values, effects) {
  pairs <- utils::combn(seq_along(effects), 2L)
  matrix_value <- matrix(0, length(effects), length(effects),
    dimnames = list(effects, effects))
  matrix_value[cbind(pairs[1L, ], pairs[2L, ])] <- values
  matrix_value <- matrix_value + t(matrix_value)
  matrix_value
}

cg_case <- function(seed, rank_deficient = FALSE) {
  set.seed(seed)
  q <- sample(3:5, 1L)
  d <- sample(3:6, 1L)
  partitions <- sample(2:4, 1L)
  effects <- paste0("e", seq_len(q))
  runs <- paste0("run", seq_len(partitions))
  blocks <- stats::setNames(lapply(runs, function(run) {
    matrix(stats::rnorm(q * d), q, d, dimnames = list(effects, NULL))
  }), runs)
  metric_rank <- if (rank_deficient) max(1L, d - 1L) else d
  factor <- matrix(stats::rnorm(d * metric_rank), d, metric_rank)
  metric <- tcrossprod(factor)
  if (!rank_deficient) metric <- metric + diag(0.2, d)
  gamma <- matrix(1, partitions, partitions,
    dimnames = list(runs, runs))
  diag(gamma) <- 0
  gamma <- gamma / sum(gamma)
  contrast <- stats::setNames(stats::rnorm(q), effects)
  contrast <- contrast - mean(contrast)
  contrast2 <- stats::setNames(stats::rnorm(q), effects)
  contrast2 <- contrast2 - mean(contrast2)

  pair_count <- q * (q - 1L) / 2L
  model_vectors <- cbind(
    model1 = stats::rnorm(pair_count),
    model2 = stats::rnorm(pair_count)
  )
  design <- cbind(`(Intercept)` = 1, model_vectors)
  if (qr(design)$rank != ncol(design)) {
    model_vectors[, 2L] <- seq_len(pair_count)
    design <- cbind(`(Intercept)` = 1, model_vectors)
  }
  models <- lapply(seq_len(ncol(model_vectors)), function(column) {
    cg_rdm_matrix(model_vectors[, column], effects)
  })
  names(models) <- colnames(model_vectors)

  oracle <- cg_oracle()
  effective_metric <- metric / d
  geometry <- oracle$oracle_common_geometry(blocks, gamma, effective_metric)
  rdm <- oracle$oracle_rdm(geometry)
  coefficient_map <- solve(crossprod(design), t(design))
  rsa <- as.numeric(coefficient_map %*% rdm)
  names(rsa) <- colnames(design)
  condition <- oracle$oracle_effective_condition(effective_metric)
  list(
    seed = seed, q = q, d = d, partitions = partitions,
    rank = metric_rank, rank_deficient = rank_deficient,
    effects = effects, blocks = blocks, metric = metric, gamma = gamma,
    contrast = contrast, contrast2 = contrast2, models = models,
    design = design, geometry = geometry, rdm = rdm, rsa = rsa,
    plan = cg_plan(blocks, metric, paste0("cg-property-", seed)),
    condition = condition,
    replay = sprintf(
      "seed=%d q=%d d=%d partitions=%d metric_rank=%d kappa=%.6g",
      seed, q, d, partitions, metric_rank, condition
    )
  )
}

cg_property_seeds <- function() {
  if (identical(Sys.getenv("CROSSFORM_DEEP_EQUIVALENCE"), "true")) {
    1201:1250
  } else {
    1201:1206
  }
}
