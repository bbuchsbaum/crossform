pg_prediction <- function(A) {
  eig <- eigen(A, symmetric = TRUE)
  amp <- pmax(eig$values, 0)
  list(factor = sweep(eig$vectors, 2, sqrt(amp), "*"), amplitudes = amp,
    shifted_spectrum = eig$values)
}
pg_score_row <- function(A, S) crossform:::.geometry_score_row(pg_prediction(A), S, 1e-10)

test_that("[T33 T34] gain agrees with independent dense and isometric packed arithmetic", {
  A <- diag(c(2, 1)); S <- diag(c(1, -1))
  result <- pg_score_row(A, S)
  expect_equal(result$inner_product, 1)
  expect_equal(result$prediction_norm_sq, 5)
  expect_equal(result$gain, -3)
  A <- matrix(c(1, .5, .5, 1), 2); S <- matrix(c(0, 1, 1, 0), 2)
  result <- pg_score_row(A, S)
  expect_equal(result$inner_product, 1)
  expect_equal(result$prediction_norm_sq, 2.5)
  expect_equal(result$gain, -.5)
  pack <- function(M) c(M[1, 1], sqrt(2) * M[2, 1], M[2, 2])
  expect_equal(result$inner_product, sum(pack(A) * pack(S)))
})

test_that("[T35 T36 T38] negative modes survive and the conditional risk identity is exact", {
  A <- diag(c(2, 1)); S <- diag(c(2, 0))
  result <- pg_score_row(A, S)
  expect_equal(result$mode_gain, c(4, -1))
  expect_equal(result$gain, 3)
  expect_identical(result$rank, 2L)
  over <- pg_score_row(matrix(3), matrix(1))
  expect_equal(over$evidence, 1)
  expect_equal(over$gain, -3)
  negative <- pg_score_row(A, diag(c(-1, -2)))
  expect_true(all(negative$evidence < 0))
  expect_true(all(negative$mode_gain < 0))
  noise <- matrix(c(1, 2, 2, -3), 2)
  for (truth in list(matrix(0, 2, 2), matrix(c(2, .5, .5, 1), 2))) {
    mean_gain <- mean(vapply(list(noise, -noise), function(N) pg_score_row(A, truth + N)$gain, 0))
    expect_equal(mean_gain, sum(truth^2) - sum((truth - A)^2))
  }
})

test_that("[T34] a frozen full operator and reduced packed scoring give the same off-diagonal evidence", {
  f <- pg_fixture()
  A <- matrix(c(1, .5, .5, 1), 2); S <- matrix(c(0, 1, 1, 0), 2)
  Z <- pg_prediction(A)$factor
  rel <- relation(list(a = f$V %*% Z, b = f$V %*% Z, c = f$V, d = f$V %*% S),
    domain = f$domain, provenance = f$rel$provenance)
  fit <- fit_geometry(pg_plan(f, "a", "b", rel = rel), f$basis, rank = 2, penalty = 0)
  plan <- pg_plan(f, rel = rel)
  score <- score_geometry(fit, plan)
  F <- pg_forms(fit)[[1]]
  query <- bilinear_query(F / 2 + t(F) / 2, effects = rel$effect_space)
  fixed <- evaluate_geometry(plan, query = query)
  expect_equal(as.numeric(fixed$values[1, ]), 1, tolerance = 1e-11)
  expect_equal(score$table$inner_product[1], as.numeric(fixed$values[1, ]), tolerance = 1e-11)
  expect_equal(score$table$gain[1], -.5, tolerance = 1e-11)
})

test_that("[T24 T33 T37 T40] the public score matches dense raw geometry and frozen component evidence", {
  f <- pg_fixture()
  fit <- fit_geometry(f$plan, f$basis, rank = 2, penalty = .2)
  plan <- pg_plan(f)
  score <- score_geometry(fit, plan, modes = TRUE, components = TRUE, row_block = 1)
  expect_s3_class(score, "effect_geometry_score")
  expect_false(inherits(score, "effect_form"))
  for (i in 1:2) {
    F <- pg_forms(fit)[[i]]
    w <- as.numeric(f$at$weights[i, ])
    G <- .88 * f$B %*% diag(w, 2) %*% t(f$B)
    expect_equal(score$table$inner_product[i], sum(F * G), tolerance = 1e-11)
    expect_equal(score$table$prediction_norm_sq[i], sum(F^2), tolerance = 1e-11)
    expect_equal(score$table$gain[i], 2 * sum(F * G) - sum(F^2), tolerance = 1e-11)
    u <- w / sqrt(sum(w))
    coherent <- .88 * tcrossprod(f$B %*% u)
    expect_equal(score$components$coherent_inner_product[i], sum(F * coherent), tolerance = 1e-11)
    expect_equal(score$components$configuration_inner_product[i], sum(F * (G - coherent)), tolerance = 1e-11)
  }
  expect_equal(2 * rowSums(score$components[-1]) - score$table$prediction_norm_sq, score$table$gain)
  details <- as.data.frame(score, view = "modes")
  expect_equal(as.numeric(tapply(details$gain, details$measurement, sum)), score$table$gain)
  expect_identical(as.data.frame(score), score$table)
  expect_output(print(score), "squared geometry units")
  reversed <- additive_frame(members = list(1L, 1:2), domain = f$domain, measurements = c("second", "first"))
  reordered <- score_geometry(fit, pg_plan(f, at = reversed), modes = TRUE, components = TRUE)
  expect_equal(as.data.frame(reordered), as.data.frame(score))
  expect_equal(as.data.frame(reordered, view = "modes"), details)
  expect_error(geometry_spectrum(score), class = "effect_input_error")
})

test_that("[T16 T22 T41] zero prediction and retained tied groups have invariant readouts", {
  f <- pg_fixture()
  zero <- score_geometry(fit_geometry(f$plan, f$basis, rank = 0, penalty = .2), pg_plan(f), modes = TRUE)
  expect_identical(zero$table$gain, c(0, 0))
  expect_identical(zero$table$prediction_norm_sq, c(0, 0))
  expect_equal(nrow(as.data.frame(zero, view = "modes")), 0L)
  # Equal positive amplitudes allow an arbitrary within-group rotation.
  original <- pg_prediction(diag(c(2, 2)))
  O <- matrix(c(cos(.7), sin(.7), -sin(.7), cos(.7)), 2)
  rotated <- original; rotated$factor <- original$factor %*% O
  S <- matrix(c(1, .3, .3, -2), 2)
  a <- crossform:::.geometry_score_row(original, S, 1e-10)
  b <- crossform:::.geometry_score_row(rotated, S, 1e-10)
  expect_identical(a$groups, c(1L, 1L))
  expect_equal(a$gain, b$gain)
  expect_equal(sum(a$evidence), sum(b$evidence))
  expect_false(isTRUE(all.equal(a$evidence, b$evidence)))
  B <- f$V %*% diag(sqrt(c(2, 2)))
  rel <- relation(list(a = B, b = B, c = f$B * .8, d = f$B * 1.1),
    domain = f$domain, provenance = f$rel$provenance)
  tied_fit <- fit_geometry(pg_plan(f, "a", "b", rel = rel), f$basis, rank = 2, penalty = 0)
  tied <- score_geometry(tied_fit, pg_plan(f, rel = rel), modes = TRUE)
  grouped <- as.data.frame(tied, view = "groups")
  expect_equal(grouped$size[grouped$measurement == "first"], 2L)
  expect_equal(grouped$gain[grouped$measurement == "first"], tied$table$gain[1])
})

test_that("[T29 T32 T41 T42] scoring is frozen, reopens, and detects inconsistent evidence", {
  f <- pg_fixture(); fit <- fit_geometry(f$plan, f$basis, rank = 2, penalty = .2)
  # R may JIT-compile read closures during evaluation. Check scientific
  # metadata and numeric content, not incidental interpreter bytecode.
  fit_before <- serialize(crossform:::.geometry_fit_semantic(fit), NULL)
  payload_before <- serialize(fit$values$read(), NULL)
  reference <- score_geometry(fit, pg_plan(f), modes = TRUE)
  testthat::local_mocked_bindings(fit_geometry = function(...) stop("refitted"),
    .geometry_fit_path = function(...) stop("refitted eigensystem"), .package = "crossform")
  path <- tempfile("prediction-score-"); withr::defer(unlink(path, recursive = TRUE))
  blocked <- score_geometry(fit, pg_plan(f), modes = TRUE, row_block = 1,
    storage = "block", storage_path = path)
  reopened <- readRDS(file.path(path, "score.rds"))
  expect_equal(as.data.frame(reopened), as.data.frame(reference))
  expect_equal(as.data.frame(reopened, view = "modes"), as.data.frame(reference, view = "modes"))
  expect_identical(serialize(crossform:::.geometry_fit_semantic(fit), NULL), fit_before)
  expect_identical(serialize(fit$values$read(), NULL), payload_before)
  expect_false(dir.exists(file.path(path, "evaluation")))
  expect_error(score_geometry(fit, pg_plan(f), storage = "block", storage_path = path), "overwrite")
  for (edit in list(
    function(x) { x$table$gain[1] <- 10; x },
    function(x) { x$table$prediction_norm_sq[1] <- 0; x },
    function(x) { x$table$rank[1] <- 999L; x },
    function(x) { x$validation$mapping <- 2:1; x },
    function(x) { x$parameters$rank <- 0; x },
    function(x) { x$complete <- FALSE; x }
  )) expect_error(as.data.frame(edit(reference)), class = "effect_contract_error")
  raw_read <- reference$mode_values$read
  broken <- reference; broken$mode_values$read <- function(rows) { v <- raw_read(rows); v[1, 1] <- v[1, 1] + 1; v }
  expect_error(as.data.frame(broken, view = "modes"), "frozen row identities")
  parts <- list(a = f$B, b = f$B, c = -f$B, d = f$B)
  altered <- relation(parts, domain = f$domain, provenance = f$rel$provenance)
  changed <- score_geometry(fit, pg_plan(f, rel = altered), modes = TRUE)
  expect_identical(changed$prediction_id, reference$prediction_id)
  expect_true(all(changed$table$gain < 0))
  expect_true(all(as.data.frame(changed, view = "modes")$evidence < 0))
})
test_that("[T35 T41 T60] score views preserve requested occurrences and stable empty schemas", {
  f <- pg_fixture()
  fit <- fit_geometry(f$plan, f$basis, rank = 2, penalty = .1)
  score <- score_geometry(fit, pg_plan(f), modes = TRUE)
  schemas <- list(summary = c("measurement", "gain", "rank", "inner_product", "prediction_norm_sq"),
    modes = c("measurement", "mode", "group", "amplitude", "evidence", "prediction_cost", "gain"),
    groups = c("measurement", "group", "size", "amplitude_sum", "evidence_sum", "prediction_cost", "gain"))
  for (view in names(schemas)) {
    empty <- as.data.frame(score, view = view, rows = integer())
    expect_s3_class(empty, "data.frame")
    expect_identical(names(empty), schemas[[view]])
    expect_identical(nrow(empty), 0L)
    selected <- as.data.frame(score, view = view, rows = c(2, 1, 1))
    reference <- do.call(rbind, lapply(c(2, 1, 1), function(row) as.data.frame(score, view = view, rows = row)))
    rownames(reference) <- NULL
    expect_equal(selected, reference)
  }
  zero <- score_geometry(fit_geometry(f$plan, f$basis, rank = 0, penalty = .1), pg_plan(f))
  for (view in c("modes", "groups")) {
    expect_identical(names(as.data.frame(zero, view = view)), schemas[[view]])
    expect_identical(nrow(as.data.frame(zero, view = view)), 0L)
  }
  omitted <- score_geometry(fit, pg_plan(f))
  expect_error(as.data.frame(omitted, view = "modes"), "not retained")
})
