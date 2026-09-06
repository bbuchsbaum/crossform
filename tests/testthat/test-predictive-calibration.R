pgc_test_helpers <- function() {
  path <- testthat::test_path("..", "..", "benchmarks", "predictive-geometry", "calibration.R")
  skip_if_not(file.exists(path), "Source-tree calibration harness is not installed with test files")
  env <- new.env(parent = globalenv()); sys.source(path, env)
  env
}

test_that("[T38 T45 T46] generator noise moment has an exact independent finite-noise witness", {
  env <- pgc_test_helpers()
  B <- matrix(c(1, -.3), 2)
  metric <- matrix(1.3); V <- matrix(.4); sigma <- c(.6, .9, 1.2)
  edges <- data.frame(left = c(1L, 1L, 2L), right = c(2L, 3L, 3L), weight = c(.2, .3, .5))
  # Six independent signs enumerate all 64 datasets, including dependent
  # cross-product edges. The second-moment oracle needs no Gaussian fourth
  # moment because all products are between independent partitions.
  signs <- expand.grid(rep(list(c(-1, 1)), 6))
  truth <- B %*% metric %*% t(B)
  squared_error <- numeric(nrow(signs))
  for (i in seq_len(nrow(signs))) {
    errors <- matrix(as.numeric(signs[i, ]), 2, 3)
    parts <- lapply(1:3, function(j) B + sigma[j] * sqrt(V[1]) * errors[, j, drop = FALSE])
    actual <- env$pgc_dense(parts, metric, edges)
    squared_error[i] <- sum((actual - truth)^2)
  }
  expect_equal(mean(squared_error), env$pgc_noise_moment(B, metric, V, sigma, edges), tolerance = 1e-12)
})

test_that("[T24 T45 T46 T54] calibration uses the same estimator as the full public route", {
  env <- pgc_test_helpers()
  withr::local_seed(84023)
  RNGkind("L'Ecuyer-CMRG"); set.seed(2026090403)
  for (arm in c("null", "aligned", "outside_model", "heterogeneous_metric")) {
    design <- env$pgc_design(arm); models <- env$pgc_models(design)
    parts <- env$pgc_parts(design, .Random.seed)
    expected <- env$pgc_one(design, models, parts, rank = 2, penalty = .15)
    names(parts) <- paste0("run", 1:6)
    domain <- abstract_domain(design$p, id = "calibration-public")
    rel <- relation(parts, domain = domain, provenance = list(observation_origins = list(
      id = "calibration-six-runs", partitions = as.list(stats::setNames(names(parts), names(parts))),
      independence = "independent", assumption = "Independent Gaussian partition errors with fixed design.")))
    at <- compile_frame(whole_brain(normalization = "none"), domain)
    metric <- neural_metric(design$metric, domain)
    make_plan <- function(edges) plan_geometry(rel, at,
      pairing(names(parts)[edges$left], names(parts)[edges$right], weight = edges$weight,
        independence = "independent", generalizes_over = "run"), metric = metric)
    train <- make_plan(design$train_edges); test <- make_plan(design$test_edges)
    for (name in names(models)) {
      fit <- fit_geometry(train, models[[name]]$basis, rank = 2, penalty = .15)
      score <- score_geometry(fit, test)
      expect_equal(score$table$gain, unname(expected[[name]]["gain"]), tolerance = 1e-10)
      expect_equal(fit$diagnostics$prediction_norm_sq, unname(expected[[name]]["norm_sq"]), tolerance = 1e-10)
      F <- pg_forms(fit)[[1]]
      expect_equal(score$table$gain - (sum(design$truth^2) - sum((design$truth - F)^2)),
        unname(expected[[name]]["error"]), tolerance = 1e-10)
    }
  }
})

test_that("[T47] known correlated origins refuse independently of the numerical score", {
  f <- pg_fixture()
  # A shared preprocessing dependency is enough to violate the required
  # train/test separation, even though the displayed run labels are distinct.
  manifest <- list(id = "shared-preprocessing", partitions = list(a = "a", b = "b", c = "c", d = "d"),
    dependencies = list(a = "calibration-a", b = character(), c = "calibration-a", d = character()),
    independence = "independent", assumption = "Nominally independent runs.")
  rel <- relation(list(a = f$B, b = f$B, c = f$B, d = f$B),
    domain = f$domain, provenance = list(observation_origins = manifest))
  fit <- fit_geometry(pg_plan(f, "a", "b", rel = rel), f$basis, rank = 2, penalty = .2)
  expect_identical(catch_refusal(score_geometry(fit, pg_plan(f, rel = rel)))$reasons[[1]],
    "overlapping_training_evaluation_origins")
})
