test_that("first-journey objects print compact summaries", {
  geometry <- result_fixture()
  query <- query_geometry(geometry, matrix(c(1, 0, -1), ncol = 1))
  contrast_view <- contrast(geometry, c(a = 1, b = -1))
  rdm_view <- rdm(geometry)
  rsa_view <- rsa(geometry, models = list(separation = matrix(
    c(0, 1, 1, 0), 2, 2
  )), intercept = FALSE)
  spectrum_view <- geometry_spectrum(geometry)

  expect_snapshot(print(geometry))
  expect_snapshot(print(query))
  expect_snapshot(print(contrast_view))
  expect_snapshot(print(rdm_view))
  expect_snapshot(print(rsa_view))
  expect_snapshot(print(spectrum_view))

  expect_match(format(geometry), "2 effects")
  expect_match(format(query), "2 measurements")
  expect_match(format(contrast_view), "energy decomposition")
  expect_match(format(rdm_view), "1 distances")
  expect_match(format(rsa_view), "1 coefficients")
  expect_match(format(spectrum_view), "2 roots")
})

test_that("materialized views coerce losslessly to row-wise data frames", {
  geometry <- result_fixture()
  query <- query_geometry(geometry,
    matrix(c(1, 0, -1, 2, 0.5, 1), nrow = 3,
      dimnames = list(NULL, c("first", "second"))))
  contrast_view <- contrast(geometry, c(a = 1, b = -1))
  rdm_view <- rdm(geometry)
  rsa_view <- rsa(geometry, models = list(separation = matrix(
    c(0, 1, 1, 0), 2, 2
  )), intercept = FALSE)
  spectrum_view <- geometry_spectrum(geometry)

  query_data <- as.data.frame(query)
  expect_identical(query_data$measurement, geometry$index)
  expect_equal(as.matrix(query_data[c("first", "second")]), query$values)

  contrast_data <- as.data.frame(contrast_view)
  expect_equal(contrast_data$total, contrast_view$total)
  expect_equal(contrast_data$configuration, contrast_view$configuration)

  rdm_data <- as.data.frame(rdm_view)
  expect_equal(as.matrix(rdm_data[-1L]), rdm_view$values,
    ignore_attr = TRUE)

  rsa_data <- as.data.frame(rsa_view)
  expect_equal(as.matrix(rsa_data[-1L]), rsa_view$coefficients,
    ignore_attr = TRUE)

  spectrum_data <- as.data.frame(spectrum_view)
  expect_equal(as.matrix(spectrum_data[-1L]), spectrum_view$values,
    ignore_attr = TRUE)
})

test_that("crossnobis views print and coerce without losing signs", {
  domain <- abstract_domain(3, id = "print-crossnobis")
  run1 <- rbind(a = c(1, 0, 0), b = c(0, 1, 0))
  run2 <- rbind(a = c(1.1, 0, 0), b = c(0, 0.9, 0))
  relation <- relation(list(run1 = run1, run2 = run2), domain = domain)
  plan <- plan_geometry(
    relation,
    compile_frame(voxels(), domain),
    cross_partitions(relation, independence = "independent"),
    metric = noise_precision(diag(3), domain)
  )
  view <- crossnobis(plan, c(a = 1, b = -1))

  data <- as.data.frame(view)
  expect_identical(data$measurement, view$index)
  expect_identical(data$crossnobis, view$values)
  expect_snapshot(print(view))
  expect_match(format(plan), "effect_geometry_plan")
  expect_match(format(view), "signed estimates")
})

test_that("relation fits print capabilities rather than nested sources", {
  design <- cbind(intercept = 1, condition = rep(c(-0.5, 0.5), 4))
  response <- list(run1 = matrix(seq_len(24) / 10, 8, 3))
  fit <- lm_relation_fit(
    response, design,
    rbind(condition = c(0, 1)),
    domain = abstract_domain(3, id = "print-fit")
  )

  expect_snapshot(print(fit))
  expect_match(format(fit), "3 features")
})
