test_that("[T06 T08 T23] model units and duplicate split weights preserve public predictions and scores", {
  f <- pg_fixture(n = 6L)
  V <- qr.Q(qr(cbind(1, diag(6))))[, 2:5, drop = FALSE]
  A <- V[, 1:2] %*% diag(c(2, .7))
  rownames(A) <- rownames(f$V)
  for (kind in c("independent", "overlap", "near_overlap")) {
    B <- switch(kind, independent = V[, 3:4], overlap = V[, 2:3],
      near_overlap = cbind(V[, 1] + 1e-4 * V[, 3], V[, 2] + 1e-4 * V[, 4]))
    B <- B %*% diag(c(1.3, .4)); rownames(B) <- rownames(A)
    declare <- function(models) model_basis(features = models,
      conditions = f$rel$effect_space, normalize = "trace", tolerance = 1e-12)
    b <- declare(list(a = A, b = B))
    fit <- fit_geometry(f$plan, b, weights = c(a = .4, b = .6), rank = 1, penalty = .02)
    gain <- as.data.frame(score_geometry(fit, pg_plan(f)))$gain
    duplicated <- declare(list(a = A, b = B, duplicate = A))
    duplicate_fit <- fit_geometry(f$plan, duplicated,
      weights = c(a = .1, b = .6, duplicate = .3), rank = 1, penalty = .02)
    expect_lt(nrow(duplicated$R), ncol(duplicated$R))
    expect_equal(pg_forms(duplicate_fit), pg_forms(fit), tolerance = 1e-7, info = kind)
    expect_equal(as.data.frame(score_geometry(duplicate_fit, pg_plan(f)))$gain,
      gain, tolerance = 1e-7, info = kind)
    for (scale in c(1e-8, 1e8)) {
      scaled <- declare(list(a = A * sqrt(scale), b = B / sqrt(scale)))
      scaled_fit <- fit_geometry(f$plan, scaled,
        weights = c(a = .4, b = .6), rank = 1, penalty = .02)
      expect_equal(pg_forms(scaled_fit), pg_forms(fit), tolerance = 1e-7,
        info = paste(kind, scale))
      expect_equal(as.data.frame(score_geometry(scaled_fit, pg_plan(f)))$gain,
        gain, tolerance = 1e-7, info = paste(kind, scale))
    }
  }
})

test_that("[T10] centered public fits and signed compression exclude a neural baseline", {
  f <- pg_fixture()
  fit <- fit_geometry(f$plan, f$basis, rank = 2, penalty = .2, retain_source = TRUE)
  baseline <- c(100, -250)
  parts <- list(a = f$B, b = f$B, c = f$B * .8, d = f$B * 1.1)
  shifted <- relation(lapply(parts, function(B) sweep(B, 2, baseline, "+")),
    domain = f$domain, provenance = f$rel$provenance)
  shifted_fit <- fit_geometry(pg_plan(f, "a", "b", rel = shifted), f$basis,
    rank = 2, penalty = .2, retain_source = TRUE)
  expect_equal(pg_forms(shifted_fit), pg_forms(fit), tolerance = 1e-10)
  expect_equal(geometry_component(shifted_fit$source, "total"),
    geometry_component(fit$source, "total"), tolerance = 1e-10)
  expect_equal(as.data.frame(score_geometry(shifted_fit, pg_plan(f, rel = shifted)))$gain,
    as.data.frame(score_geometry(fit, pg_plan(f)))$gain, tolerance = 1e-10)
})

test_that("[T42 T55] fitted amplitudes cannot inherit raw spectrum or covariance inference", {
  f <- pg_fixture()
  fit <- fit_geometry(f$plan, f$basis, rank = 1, penalty = .2)
  expect_error(geometry_spectrum(fit), class = "effect_input_error")
  expect_error(sampling_covariance(fit), "must come from `rdm_sampling_covariance", fixed = TRUE)
  expect_s3_class(score_geometry(fit, pg_plan(f)), "effect_geometry_score")
  B <- f$B; rownames(B) <- paste0("new-condition-", seq_len(nrow(B)))
  changed <- relation(list(a = B, b = B, c = B, d = B),
    domain = f$domain, provenance = f$rel$provenance)
  testthat::local_mocked_bindings(.geometry_fit_materialize = function(...) stop("neural read"),
    .package = "crossform")
  expect_identical(catch_refusal(score_geometry(fit, pg_plan(f, rel = changed)))$reasons[[1]],
    "prediction_effect_space_mismatch")
})

test_that("[T23] centered public fits retain the dense objective across dimensions, gaps and scales", {
  withr::local_seed(62203)
  for (n in c(6L, 9L)) for (r in c(2L, 4L)) for (condition in c(10, 1e7)) {
    V <- qr.Q(qr(cbind(1, diag(n))))[, seq.int(2L, r + 1L), drop = FALSE]
    d <- exp(seq(0, -log(condition), length.out = r)); d <- d / sum(d)
    K <- V %*% diag(d) %*% t(V)
    O <- qr.Q(qr(matrix(rnorm(r*r), r)))
    S <- O %*% diag(seq(5, -1, length.out = r)) %*% t(O)
    ids <- paste0("effect", seq_len(n)); dimnames(K) <- list(ids, ids)
    basis <- model_basis(kernels = K, normalize = "trace", tolerance = 1e-12)
    for (scale in c(1e-8, 1, 1e8)) {
      G <- V %*% (scale * S) %*% t(V)
      B1 <- diag(n); B2 <- G; rownames(B1) <- rownames(B2) <- ids
      domain <- abstract_domain(n, id = "stress-objective")
      rel <- relation(list(a = B1, b = B2), domain = domain)
      plan <- plan_geometry(rel, compile_frame(whole_brain(normalization = "none"), domain),
        pairing("a", "b", independence = "independent", generalizes_over = "run"))
      penalty <- .003 * scale
      fit <- fit_geometry(plan, basis, rank = 1, penalty = penalty)
      actual <- unname(pg_forms(fit)[[1]])
      e <- eigen(scale * S - penalty * diag(1/d), symmetric = TRUE)
      reference <- max(e$values[1], 0) * tcrossprod(V %*% e$vectors[, 1])
      # Compare the original full-space objective using the independently
      # declared kernel eigenvectors, not the compiler's compressed basis.
      objective <- function(F) .5 * sum((G - F)^2) +
        penalty * sum(diag(crossprod(V, F %*% V)) / d)
      expect_equal(actual / scale, reference / scale, tolerance = 1e-7)
      expect_equal(objective(actual) / scale^2, objective(reference) / scale^2, tolerance = 1e-7)
      expect_lte(objective(actual), objective(matrix(0, n, n)) + 1e-7 * scale^2)
      expect_true(all(is.finite(actual)))
      expect_lt(max(abs(actual - t(actual))) / scale, 1e-10)
    }
  }
})
