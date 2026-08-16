axis_fixture <- function(domain_id = "generalization-axis-domain") {
  set.seed(60411)
  design <- cbind(1, condition = rep(c(-0.5, 0.5), 8))
  effects <- rbind(level = c(1, 0), condition = c(0, 1))
  sources <- stats::setNames(lapply(seq_len(3L), function(index) {
    matrix(rnorm(16 * 5), 16, 5)
  }), paste0("p", seq_len(3L)))
  domain <- abstract_domain(5L, id = domain_id)
  fit <- lm_relation_fit(
    sources, design, effects, sampling_unit = "trial", domain = domain
  )
  list(fit = fit, frame = compile_frame(whole_brain(), domain))
}

test_that("the generalization axis is declared, validated, and retained", {
  over <- cross_partitions(paste0("p", 1:3), generalizes_over = "run")
  expect_identical(attr(over, "generalizes_over", exact = TRUE), "run")
  expect_silent(crossform:::.validate_pairing(over))
  expect_error(
    pairing("a", "b", generalizes_over = ""),
    "one nonempty axis name"
  )
})

test_that("cross-run and cross-session estimands differ under generic labels", {
  fixture <- axis_fixture()
  run_plan <- plan_geometry(
    fixture$fit$relation, fixture$frame,
    cross_partitions(fixture$fit$relation, generalizes_over = "run")
  )
  session_plan <- plan_geometry(
    fixture$fit$relation, fixture$frame,
    cross_partitions(fixture$fit$relation, generalizes_over = "session")
  )
  undeclared_plan <- plan_geometry(
    fixture$fit$relation, fixture$frame,
    cross_partitions(fixture$fit$relation)
  )
  identities <- c(
    run = run_plan$scientific_plan_id,
    session = session_plan$scientific_plan_id,
    undeclared = undeclared_plan$scientific_plan_id
  )
  expect_length(unique(identities), 3L)

  # The axis names the estimand, not the numbers: values agree exactly.
  weights <- c(level = 0, condition = 1)
  expect_identical(
    contrast_energy(run_plan, weights)$total,
    contrast_energy(session_plan, weights)$total
  )
})

test_that("the declared axis reaches metric-pairing identity", {
  run_pairing <- cross_partitions(paste0("p", 1:3), generalizes_over = "run")
  session_pairing <- cross_partitions(
    paste0("p", 1:3), generalizes_over = "session"
  )
  undeclared_pairing <- cross_partitions(paste0("p", 1:3))
  run_identity <- crossform:::.metric_pairing_identity(run_pairing)
  session_identity <- crossform:::.metric_pairing_identity(session_pairing)
  undeclared_identity <- crossform:::.metric_pairing_identity(
    undeclared_pairing
  )
  expect_identical(run_identity$generalizes_over, "run")
  expect_false(identical(run_identity, session_identity))
  expect_false("generalizes_over" %in% names(undeclared_identity))
})
