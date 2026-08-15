route_fixture <- function() {
  set.seed(71203)
  design <- cbind(1, condition = rep(c(-0.5, 0.5), 8))
  effects <- rbind(level = c(1, 0), condition = c(0, 1), extra = c(1, -1))
  sources <- stats::setNames(lapply(seq_len(3L), function(index) {
    matrix(rnorm(16 * 6), 16, 6)
  }), paste0("run", seq_len(3L)))
  domain <- abstract_domain(6L, id = "route-identity-domain")
  fit <- lm_relation_fit(
    sources, design, effects, sampling_unit = "trial", domain = domain
  )
  plan <- plan_geometry(
    fit$relation, compile_frame(voxels(), domain),
    cross_partitions(fit$relation)
  )
  list(plan = plan, form = geometry(plan))
}

test_that("rdm carries one estimand identity on both execution routes", {
  fixture <- route_fixture()
  fused <- rdm(fixture$plan)
  projected <- rdm(fixture$form)
  expect_identical(
    fused$receipt$scientific_plan_id,
    projected$receipt$scientific_plan_id
  )
  expect_equal(fused$values, projected$values, tolerance = 1e-12)
  # The receipts still record different executions of that one estimand.
  expect_false(identical(
    fused$receipt$task_partition_id,
    projected$receipt$task_partition_id
  ))
})

test_that("contrast carries one estimand identity on both execution routes", {
  fixture <- route_fixture()
  weights <- c(level = 0, condition = 1, extra = 0)
  fused <- contrast(fixture$plan, weights)
  projected <- contrast(fixture$form, weights)
  expect_identical(
    fused$receipt$scientific_plan_id,
    projected$receipt$scientific_plan_id
  )
  expect_equal(fused$total, projected$total, tolerance = 1e-12)
  expect_equal(fused$coherent, projected$coherent, tolerance = 1e-12)
})

test_that("rsa carries one estimand identity on both execution routes", {
  fixture <- route_fixture()
  model <- matrix(0, 3, 3,
    dimnames = list(
      c("level", "condition", "extra"), c("level", "condition", "extra")
    )
  )
  model["level", "condition"] <- model["condition", "level"] <- 1
  model["condition", "extra"] <- model["extra", "condition"] <- 2
  model["level", "extra"] <- model["extra", "level"] <- 3
  fused <- rsa(fixture$plan, models = list(category = model))
  projected <- rsa(fixture$form, models = list(category = model))
  expect_identical(
    fused$receipt$scientific_plan_id,
    projected$receipt$scientific_plan_id
  )
  expect_equal(
    fused$coefficients, projected$coefficients, tolerance = 1e-12
  )
})

test_that("distinct view estimands stay distinct from each other", {
  fixture <- route_fixture()
  ids <- c(
    plan = fixture$plan$scientific_plan_id,
    full = fixture$form$receipt$scientific_plan_id,
    rdm_total = rdm(fixture$plan)$receipt$scientific_plan_id,
    rdm_configuration = rdm(
      fixture$plan, component = "configuration"
    )$receipt$scientific_plan_id,
    contrast = contrast(
      fixture$plan, c(level = 0, condition = 1, extra = 0)
    )$receipt$scientific_plan_id
  )
  # The plan and its full materialization are one estimand; every distinct
  # view is its own estimand derived from that parent.
  expect_identical(ids[["plan"]], ids[["full"]])
  expect_length(unique(ids[-1L]), 4L)
})
