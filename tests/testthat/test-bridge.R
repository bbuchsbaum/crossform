bridge_task_relation <- function(partitions, effects, domain, reads = NULL,
                                 revision_character = "a") {
  effect_count <- length(effects$coordinates)
  if (is.null(reads)) {
    sources <- stats::setNames(lapply(seq_along(partitions), function(index) {
      matrix(index, effect_count, domain$reference$n_features)
    }), partitions)
    return(relation(sources, effects = effects, domain = domain))
  }
  source <- function(features) {
    reads$count <- reads$count + 1L
    matrix(1, effect_count, length(features))
  }
  revision <- paste0(
    "sha256:", paste(rep(revision_character, 64), collapse = "")
  )
  relation(
    stats::setNames(rep(list(source), length(partitions)), partitions),
    source_dims = rep(list(c(effect_count, domain$reference$n_features)),
      length(partitions)),
    effects = effects,
    domain = domain,
    capabilities = source_capabilities(TRUE, stable_revision = revision)
  )
}

bridge_fixture <- function(lazy = FALSE) {
  left_domain <- abstract_domain(3, id = "left-neural:v1")
  right_domain <- abstract_domain(4, id = "right-neural:v1")
  common <- crossform:::measurement_space(2, "common-measurement:v1",
    provenance = list(method = "fixed-reference"))
  left_leg <- matrix(c(1, 0, 2, 0, 1, -1), 2, 3, byrow = TRUE)
  right_leg <- matrix(c(1, 2, 0, -1, 0, 1, 3, 2), 2, 4, byrow = TRUE)
  bridge <- crossform:::measurement_bridge(
    left_leg, right_leg, left_domain, right_domain, common,
    provenance = list(revision = "bridge:v1", learned = FALSE)
  )
  reads <- new.env(parent = emptyenv())
  reads$count <- 0L
  left_space <- effect_space(c("l1", "l2"), basis_id = "left-effects:v1")
  right_space <- effect_space(c("r1", "r2", "r3"),
    basis_id = "right-effects:v1")
  if (lazy) {
    left <- bridge_task_relation(
      c("e1", "e2"), left_space, left_domain, reads, "d"
    )
    right <- bridge_task_relation(
      c("r1", "r2"), right_space, right_domain, reads, "e"
    )
  } else {
    left <- bridge_task_relation(c("e1", "e2"), left_space, left_domain)
    right <- bridge_task_relation(c("r1", "r2"), right_space, right_domain)
  }
  list(
    left_domain = left_domain,
    right_domain = right_domain,
    common = common,
    left_leg = left_leg,
    right_leg = right_leg,
    bridge = bridge,
    left = left,
    right = right,
    left_space = left_space,
    right_space = right_space,
    reads = reads,
    over = pairing(c("e1", "e2"), c("r1", "r2"),
      c(1, 2), directed = TRUE)
  )
}

test_that("factorized common-coordinate evaluation equals direct K", {
  fixture <- bridge_fixture()
  set.seed(98)
  left <- matrix(rnorm(2 * 3), 2, 3)
  right <- matrix(rnorm(3 * 4), 3, 4)
  common <- left %*% t(fixture$left_leg) %*%
    t(right %*% t(fixture$right_leg))
  K <- t(fixture$left_leg) %*% fixture$right_leg
  direct <- left %*% K %*% t(right)

  expect_equal(common, direct, tolerance = 1e-14)
  applied <- crossform:::.apply_measurement_bridge(
    list(left = left), list(right = right), fixture$bridge
  )
  expect_equal(applied$left$left %*% t(applied$right$right),
    direct, tolerance = 1e-14)
  expect_identical(applied$bridge_signature, fixture$bridge$signature)
})

test_that("joint measured block is positive semidefinite", {
  fixture <- bridge_fixture()
  set.seed(23)
  left <- matrix(rnorm(2 * 3), 2, 3)
  right <- matrix(rnorm(3 * 4), 3, 4)
  joint <- crossform:::.bridge_joint_gram(left, right, fixture$bridge)

  expect_equal(joint, t(joint), tolerance = 1e-15)
  expect_gte(min(eigen(joint, symmetric = TRUE, only.values = TRUE)$values),
    -1e-12)
})

test_that("bridge reversal transposes the induced form", {
  fixture <- bridge_fixture()
  set.seed(11)
  left <- matrix(rnorm(2 * 3), 2, 3)
  right <- matrix(rnorm(3 * 4), 3, 4)
  reversed <- crossform:::reverse_bridge(fixture$bridge)
  forward <- left %*% t(fixture$left_leg) %*%
    t(right %*% t(fixture$right_leg))
  backward <- right %*% t(reversed$left_leg) %*%
    t(left %*% t(reversed$right_leg))

  expect_equal(backward, t(forward), tolerance = 1e-14)
  expect_identical(reversed$left_domain, fixture$bridge$right_domain)
  expect_identical(reversed$right_domain, fixture$bridge$left_domain)
  expect_identical(reversed$common_space, fixture$common)
  expect_identical(crossform:::reverse_bridge(reversed), fixture$bridge)
})

test_that("bridge legs have canonical numerical identity", {
  fixture <- bridge_fixture()
  named_integer_left <- matrix(as.integer(fixture$left_leg), 2, 3,
    dimnames = list(c("m1", "m2"), c("x", "y", "z")))
  canonical <- crossform:::measurement_bridge(
    named_integer_left, fixture$right_leg,
    fixture$left_domain, fixture$right_domain, fixture$common,
    provenance = list(revision = "bridge:v1", learned = FALSE)
  )

  expect_identical(canonical, fixture$bridge)
  expect_error(crossform:::measurement_space(1.5, "bad"), "n_measurements",
    class = "effect_input_error")
})

test_that("distinct neural spaces require an explicit bridge before reads", {
  fixture <- bridge_fixture(lazy = TRUE)
  expect_error(crossform:::.compile_effect_task(
    fixture$left, fixture$over, fixture$right
  ), "require an explicit", class = "effect_contract_error")
  expect_identical(fixture$reads$count, 0L)

  task <- crossform:::.compile_effect_task(
    fixture$left, fixture$over, fixture$right,
    bridge = fixture$bridge
  )
  expect_identical(task$bridge, fixture$bridge)
  expect_identical(fixture$reads$count, 0L)
  expect_silent(crossform:::.validate_compiled_effect_task(task))
})

test_that("bridge incompatibilities fail before lazy source reads", {
  fixture <- bridge_fixture(lazy = TRUE)
  wrong_left <- abstract_domain(3, feature_ids = 3:1, id = "left-neural:v1")
  wrong <- crossform:::measurement_bridge(
    fixture$left_leg, fixture$right_leg,
    wrong_left, fixture$right_domain, fixture$common
  )
  expect_error(crossform:::.compile_effect_task(
    fixture$left, fixture$over, fixture$right, bridge = wrong
  ), "left neural-space", class = "effect_contract_error")
  expect_identical(fixture$reads$count, 0L)

  expect_error(crossform:::measurement_bridge(
    fixture$left_leg[-1, , drop = FALSE], fixture$right_leg,
    fixture$left_domain, fixture$right_domain, fixture$common
  ), "common measurement dimension", class = "effect_input_error")
  forged <- fixture$bridge
  forged$common_space <- crossform:::measurement_space(2, "other-common:v1")
  expect_error(crossform:::.validate_measurement_bridge(
    forged, fixture$left, fixture$right
  ), "identity is inconsistent", class = "effect_contract_error")
  expect_identical(fixture$reads$count, 0L)
})

test_that("bridge identity participates in compiled task identity", {
  fixture <- bridge_fixture()
  first <- crossform:::.compile_effect_task(
    fixture$left, fixture$over, fixture$right, bridge = fixture$bridge
  )
  changed_bridge <- crossform:::measurement_bridge(
    fixture$left_leg * 2, fixture$right_leg,
    fixture$left_domain, fixture$right_domain, fixture$common,
    provenance = list(revision = "bridge:v2")
  )
  changed <- crossform:::.compile_effect_task(
    fixture$left, fixture$over, fixture$right, bridge = changed_bridge
  )

  expect_false(identical(first$task_id, changed$task_id))
  expect_false(identical(first$bridge$signature, changed$bridge$signature))
})
