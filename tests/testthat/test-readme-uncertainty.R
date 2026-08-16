test_that("the README fixed-metric uncertainty journey runs end to end", {
  set.seed(81370)
  conditions <- c("face", "house", "object")
  features <- 6L
  condition <- factor(rep(conditions, each = 8L), levels = conditions)
  X <- stats::model.matrix(~ 0 + condition)
  colnames(X) <- conditions
  C <- diag(length(conditions))
  rownames(C) <- conditions
  truth <- rbind(
    face = c(0.6, 0.3, 0.1, -0.1, 0, 0.2),
    house = c(0, 0.2, 0.6, 0.3, 0.1, 0),
    object = c(0.1, 0.1, 0.3, 0.5, 0.3, 0.1)
  )
  domain <- abstract_domain(features, id = "readme-uncertainty")
  response_runs <- stats::setNames(lapply(seq_len(4L), function(run) {
    X %*% truth + matrix(rnorm(nrow(X) * features, sd = 0.25),
      nrow(X), features)
  }), paste0("run", seq_len(4L)))
  fit <- lm_relation_fit(
    response_runs, design = X, effects = C,
    effect_names = conditions, sampling_unit = "trial", domain = domain
  )
  known_precision <- noise_precision(
    diag(features), domain, covariance = diag(features),
    provenance = list(source = "fixed-before-evaluation")
  )
  fixed_plan <- plan_geometry(
    fit$relation,
    at = compile_frame(whole_brain(), domain),
    over = cross_partitions(fit$relation, independence = "independent"),
    metric = known_precision
  )
  distance_estimate <- rdm(fixed_plan)
  distance_covariance <- rdm_sampling_covariance(
    fixed_plan, fit, target = "plugin", at = 1L
  )
  distance_variance <- sampling_covariance(distance_covariance)

  expect_s3_class(distance_estimate, "effect_rdm_view")
  expect_s3_class(distance_covariance, "effect_rdm_sampling_covariance")
  expect_identical(names(distance_variance),
    colnames(distance_estimate$values))
  expect_true(all(is.finite(distance_variance)))
  expect_true(all(distance_variance > 0))
})
