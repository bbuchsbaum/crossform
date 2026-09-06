test_that("[T10 T16 T39] a fit freezes reduced predictions, signed spectra and its exact scope", {
  f <- pg_fixture()
  fit <- fit_geometry(f$plan, f$basis, rank = 2, penalty = .2, retain_source = TRUE, row_block = 1)
  expect_s3_class(fit, "effect_geometry_fit")
  expect_false(inherits(fit, "effect_form"))
  expect_false(fit$interpretation$unbiased)
  expect_identical(fit$interpretation$preference, "kernel_regularized")
  expect_identical(fit$training$observations, c("raw-a", "raw-b"))
  expect_identical(fit$target$measurement_ids, c("first", "second"))
  expect_identical(fit$parameters$rank, 2L)
  expect_identical(fit$parameters$penalty, .2)
  expect_true(fit$complete)
  expected <- f$V %*% diag(c(6 - .2/.6, 3 - .2/.4)) %*% t(f$V)
  expect_equal(unname(pg_forms(fit)[[1]]), unname(expected), tolerance = 1e-11)
  expect_identical(fit$values$dim, c(2L, 10L)) # 2*2 factors, 2 amplitudes, 2+2 spectra
  expect_identical(fit$source$effect_space, f$basis$effect_space)
  expect_equal(crossform:::.geometry_prediction_rows(fit, 1)[[1]]$signed_spectrum, c(6, 3), tolerance = 1e-11)
  expect_identical(as.data.frame(fit), fit$diagnostics)
  expect_output(print(fit), "fitting is biased")
  expect_error(geometry_spectrum(fit), class = "effect_input_error")
  expect_error(as.data.frame(fit, rank = 1), class = "effect_input_error")
})

test_that("[T16 T21] zero ranks, clamping and full-span interpretation are explicit", {
  f <- pg_fixture()
  zero <- fit_geometry(f$plan, f$basis, rank = 0, penalty = .2)
  expect_identical(zero$layout$budget, 0L)
  expect_identical(zero$diagnostics$rank, c(0L, 0L))
  expect_identical(zero$diagnostics$prediction_norm_sq, c(0, 0))
  expect_true(all(vapply(pg_forms(zero), function(A) all(A == 0), FALSE)))
  oversized <- fit_geometry(f$plan, f$basis, rank = 99, penalty = .2)
  expect_identical(oversized$parameters$rank, 99L)
  expect_identical(oversized$layout$budget, 2L)
  full <- pg_fixture(n = 3L)
  unregularized <- fit_geometry(full$plan, full$basis, rank = 2, penalty = 0)
  regularized <- fit_geometry(full$plan, full$basis, rank = 2, penalty = .2)
  expect_true(unregularized$interpretation$generic_full_span_baseline)
  expect_false(regularized$interpretation$generic_full_span_baseline)
  expect_identical(regularized$interpretation$preference, "kernel_regularized")
})

test_that("[T29 T41] serialization and block reopening require no training reads", {
  f <- pg_fixture()
  memory <- fit_geometry(f$plan, f$basis, rank = 2, penalty = .2, row_block = 1)
  path <- tempfile("prediction-fit-")
  withr::defer(unlink(path, recursive = TRUE))
  blocked <- fit_geometry(f$plan, f$basis, rank = 2, penalty = .2,
    storage = "block", storage_path = path, row_block = 1)
  expect_null(blocked$source)
  expect_false(dir.exists(file.path(path, "training")))
  reopened <- readRDS(file.path(path, "fit.rds"))
  expect_identical(reopened$prediction_id, memory$prediction_id)
  expect_equal(pg_forms(reopened), pg_forms(memory), tolerance = 1e-12)
  expect_identical(unserialize(serialize(memory, NULL))$prediction_id, memory$prediction_id)
  for (p in f$rel$partitions) f$rel$sources[[p]]$read <- function(...) stop("training reread")
  expect_s3_class(crossform:::.validate_geometry_fit(reopened), "effect_geometry_fit")
  expect_equal(pg_forms(reopened), pg_forms(memory), tolerance = 1e-12)
  expect_error(fit_geometry(f$plan, f$basis, rank = 2, penalty = .2,
    storage = "block", storage_path = path), "overwrite")
})

test_that("[T29 T42] identity-bearing and derived fields cannot silently change", {
  f <- pg_fixture()
  fit <- fit_geometry(f$plan, f$basis, rank = 2, penalty = .2)
  for (edit in list(
    function(x) { x$parameters$penalty <- .1; x },
    function(x) { x$diagnostics$prediction_norm_sq[1] <- 1; x },
    function(x) { x$target$measurement_ids <- rev(x$target$measurement_ids); x },
    function(x) { x$interpretation$unbiased <- TRUE; x },
    function(x) { x$complete <- FALSE; x },
    function(x) { x$training$observations <- "raw-c"; x }
  )) expect_error(crossform:::.validate_geometry_fit(edit(fit)), class = "effect_contract_error")
  original_read <- fit$values$read
  changed <- fit
  changed$values$read <- function(rows = NULL) {
    value <- original_read(rows); value[1, 1] <- value[1, 1] + .1; value
  }
  expect_error(crossform:::.geometry_prediction_rows(changed, 1), "frozen row identities")
  # A changed factor with freshly recomputed byte digests still has to obey
  # its amplitude/support identities; byte integrity is not the sole check.
  payload <- original_read(1:2); payload[1, 1] <- payload[1, 1] + .2
  changed <- fit; changed$values <- crossform:::.memory_geometry_store(payload, "rectangular")
  changed$row_signatures <- vapply(1:2, function(i) crossform:::.sha256_signature(as.numeric(payload[i, ])), "")
  changed$signature <- crossform:::.sha256_signature(crossform:::.geometry_fit_semantic(changed), "geometry-fit-sha256:")
  expect_error(crossform:::.geometry_prediction_rows(changed, 1), "derived fit diagnostics")
})

test_that("[T30 T55] fitting without ancestry stays a latent fit and unsupported losses refuse", {
  f <- pg_fixture(origins = FALSE)
  fit <- fit_geometry(f$plan, f$basis, rank = 1, penalty = .2)
  expect_false(fit$training$known)
  expect_output(print(fit), "independent scoring is unavailable")
  expect_identical(catch_refusal(fit_geometry(f$plan, f$basis, rank = 1, penalty = .2,
    loss = "covariance_weighted_rdm"))$reasons[[1]], "predictive_loss_not_implemented")
  expect_error(fit_geometry(f$plan, f$basis, rank = 1), "both")
  expect_error(fit_geometry(f$plan, f$basis, penalty = .2), "both")
})
