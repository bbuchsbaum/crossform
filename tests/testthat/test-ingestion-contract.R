# Independent algebra and identity oracles for first-moment-relation-v1.
# These fixtures deliberately do not call production effectagram functions.

ingestion_oracle_digest <- function(value, prefix) {
  paste0(prefix, digest::digest(value, algo = "sha256", serialize = TRUE))
}

ingestion_oracle_inverse <- function(value, tolerance = 1e-12) {
  decomposition <- svd(value, nu = min(dim(value)), nv = ncol(value))
  keep <- decomposition$d > tolerance * max(decomposition$d)
  decomposition$v[, keep, drop = FALSE] %*%
    (t(decomposition$u[, keep, drop = FALSE]) / decomposition$d[keep])
}

ingestion_oracle_extractor <- function(design, target, whitener,
                                       route = c("qr", "svd")) {
  route <- match.arg(route)
  whitened <- whitener %*% design
  inverse <- if (route == "qr") {
    solve(crossprod(whitened), t(whitened))
  } else {
    ingestion_oracle_inverse(whitened)
  }
  target %*% inverse %*% whitener
}

ingestion_oracle_fixture <- function(observations = 13L) {
  stopifnot(observations >= 9L)
  condition <- rep(c("face", "place", "object"), length.out = observations)
  semantic <- stats::model.matrix(~ condition - 1)
  colnames(semantic) <- c("face", "object", "place")
  semantic <- semantic[, c("face", "place", "object"), drop = FALSE]
  drift <- seq(-0.83, 1.17, length.out = observations)
  cell_design <- cbind(semantic, drift = drift)

  # beta_cell = A %*% beta_treatment, hence X_treatment = X_cell %*% A.
  coding_map <- rbind(
    c(1, 0, 0, 0),
    c(1, 1, 0, 0),
    c(1, 0, 1, 0),
    c(0, 0, 0, 1)
  )
  treatment_design <- cell_design %*% coding_map
  colnames(treatment_design) <- c(
    "intercept", "place_minus_face", "object_minus_face", "drift"
  )

  functional <- rbind(
    place_minus_object = c(0, 1, -1),
    face_vs_others = c(1, -0.5, -0.5)
  )
  cell_target <- cbind(functional, drift = 0)
  treatment_target <- cell_target %*% coding_map
  whitener <- diag(seq(0.81, 1.19, length.out = observations))

  set.seed(2026081502)
  response <- matrix(rnorm(observations * 7L), observations, 7L)
  list(
    semantic_coordinates = colnames(semantic),
    functional = functional,
    cell_design = cell_design,
    treatment_design = treatment_design,
    cell_target = cell_target,
    treatment_target = treatment_target,
    whitener = whitener,
    response = response
  )
}

ingestion_oracle_semantic_plan <- function(
    functional,
    sampling_unit = "scan",
    observation_kind = "fixed_gls") {
  list(
    contract = "first-moment-relation-v1",
    study = list(
      observation_axis = "scan-within-run",
      partitions = c("run-1", "run-2"),
      neural_domain = "voxel-grid-v1"
    ),
    model = list(
      semantic_mean = "canonical-hrf-amplitude-by-condition",
      conditions = c("face", "place", "object"),
      nuisance = "linear-drift"
    ),
    effects = list(
      coordinates = rownames(functional),
      functional = unname(functional),
      units = "percent-signal-change"
    ),
    observation_model = list(
      kind = observation_kind,
      sampling_unit = sampling_unit,
      whitener_status = if (identical(observation_kind, "fixed_gls")) {
        "fixed"
      } else {
        "learned"
      }
    )
  )
}

ingestion_oracle_plan_id <- function(plan) {
  ingestion_oracle_digest(plan, "relation-plan-sha256:")
}

ingestion_oracle_receipt_id <- function(plan_id, design, target,
                                        coding, solver) {
  ingestion_oracle_digest(list(
    plan_id = plan_id,
    coding = coding,
    solver = solver,
    design = unname(design),
    lowering = unname(target)
  ), "design-receipt-sha256:")
}

ingestion_oracle_fit_id <- function(plan_id, receipt_id, source_revision) {
  ingestion_oracle_digest(list(
    plan_id = plan_id,
    receipt_id = receipt_id,
    source_revision = source_revision
  ), "relation-fit-sha256:")
}

ingestion_oracle_raw_plan_id <- function(design, target) {
  ingestion_oracle_digest(list(
    route = "raw-X-T",
    design = unname(design),
    target = unname(target)
  ), "raw-relation-plan-sha256:")
}

ingestion_oracle_estimability <- function(design, target,
                                          tolerance = 1e-12) {
  decomposition <- svd(design, nu = min(dim(design)), nv = ncol(design))
  keep <- decomposition$d > tolerance * max(decomposition$d)
  basis <- decomposition$v[, keep, drop = FALSE]
  projected <- target %*% basis %*% t(basis)
  scale <- pmax(1, sqrt(rowSums(target^2)))
  sqrt(rowSums((target - projected)^2)) / scale <= tolerance * 10
}

ingestion_oracle_generalization_request <- function(hierarchy, axis) {
  if (missing(axis)) {
    stop("generalizes_over must be requested explicitly", call. = FALSE)
  }
  if (!axis %in% names(hierarchy)) {
    stop("the requested generalization axis is absent", call. = FALSE)
  }
  list(generalizes_over = axis, levels = hierarchy[[axis]])
}

ingestion_oracle_capabilities <- function(observation_kind) {
  fixed <- identical(observation_kind, "fixed_gls")
  list(
    identified_relation = TRUE,
    residual_channel = TRUE,
    fixed_observation_model = fixed,
    analytic_effect_covariance = fixed,
    learned_observation_model = !fixed
  )
}

test_that("condition-space effects are invariant to supported coding routes", {
  fixture <- ingestion_oracle_fixture()

  cell <- ingestion_oracle_extractor(
    fixture$cell_design,
    fixture$cell_target,
    fixture$whitener,
    "qr"
  )
  treatment <- ingestion_oracle_extractor(
    fixture$treatment_design,
    fixture$treatment_target,
    fixture$whitener,
    "qr"
  )

  expect_equal(cell, treatment, tolerance = 1e-12)
  expect_equal(
    cell %*% fixture$response,
    treatment %*% fixture$response,
    tolerance = 1e-12
  )
})

test_that("route-stable plan identity is separate from receipts and fits", {
  fixture <- ingestion_oracle_fixture()
  semantic_plan <- ingestion_oracle_semantic_plan(fixture$functional)
  plan_id <- ingestion_oracle_plan_id(semantic_plan)

  cell_receipt <- ingestion_oracle_receipt_id(
    plan_id, fixture$cell_design, fixture$cell_target, "cell-means", "qr"
  )
  treatment_receipt <- ingestion_oracle_receipt_id(
    plan_id,
    fixture$treatment_design,
    fixture$treatment_target,
    "treatment",
    "svd"
  )

  expect_identical(plan_id, ingestion_oracle_plan_id(semantic_plan))
  expect_false(identical(cell_receipt, treatment_receipt))
  expect_false(identical(
    ingestion_oracle_fit_id(plan_id, cell_receipt, "source-a"),
    ingestion_oracle_fit_id(plan_id, treatment_receipt, "source-a")
  ))
  expect_false(identical(
    ingestion_oracle_fit_id(plan_id, cell_receipt, "source-a"),
    ingestion_oracle_fit_id(plan_id, cell_receipt, "source-b")
  ))
})

test_that("raw matrix routes honestly bind parameterization into identity", {
  fixture <- ingestion_oracle_fixture()
  expect_false(identical(
    ingestion_oracle_raw_plan_id(
      fixture$cell_design, fixture$cell_target
    ),
    ingestion_oracle_raw_plan_id(
      fixture$treatment_design, fixture$treatment_target
    )
  ))
})

test_that("semantic requests and assumptions change plan identity", {
  fixture <- ingestion_oracle_fixture()
  base <- ingestion_oracle_semantic_plan(fixture$functional)
  changed_functional <- fixture$functional
  changed_functional[1, ] <- c(-1, 1, 0)

  ids <- c(
    ingestion_oracle_plan_id(base),
    ingestion_oracle_plan_id(ingestion_oracle_semantic_plan(
      changed_functional
    )),
    ingestion_oracle_plan_id(ingestion_oracle_semantic_plan(
      fixture$functional, sampling_unit = "trial"
    )),
    ingestion_oracle_plan_id(ingestion_oracle_semantic_plan(
      fixture$functional, observation_kind = "learned_ar"
    ))
  )
  expect_length(unique(ids), 4L)
})

test_that("QR and SVD routes agree on a floating-point relation", {
  fixture <- ingestion_oracle_fixture()
  qr_map <- ingestion_oracle_extractor(
    fixture$cell_design,
    fixture$cell_target,
    fixture$whitener,
    "qr"
  )
  svd_map <- ingestion_oracle_extractor(
    fixture$cell_design,
    fixture$cell_target,
    fixture$whitener,
    "svd"
  )
  expect_equal(qr_map, svd_map, tolerance = 1e-12)
  expect_equal(
    qr_map %*% fixture$response,
    svd_map %*% fixture$response,
    tolerance = 1e-12
  )
})

test_that("unequal partitions retain one semantic effect axis", {
  short <- ingestion_oracle_fixture(11L)
  long <- ingestion_oracle_fixture(17L)
  short_fit <- ingestion_oracle_extractor(
    short$cell_design, short$cell_target, short$whitener, "qr"
  ) %*% short$response
  long_fit <- ingestion_oracle_extractor(
    long$treatment_design, long$treatment_target, long$whitener, "svd"
  ) %*% long$response

  expect_identical(rownames(short_fit), rownames(long_fit))
  expect_identical(rownames(short_fit), rownames(short$functional))
  expect_identical(ncol(short_fit), ncol(long_fit))
})

test_that("estimability is a property of named functionals per partition", {
  condition <- rep(c("face", "place", "object"), each = 3L)
  indicators <- stats::model.matrix(~ condition - 1)
  design <- cbind(intercept = 1, indicators)
  colnames(design) <- c("intercept", "face", "object", "place")

  individual <- rbind(face_coefficient = c(0, 1, 0, 0))
  difference <- rbind(face_minus_place = c(0, 1, 0, -1))
  expect_false(ingestion_oracle_estimability(design, individual))
  expect_true(ingestion_oracle_estimability(design, difference))
})

test_that("hierarchy supplies vocabulary but never a generalization default", {
  hierarchy <- list(
    run = c("run-1", "run-2", "run-3"),
    session = c("session-1", "session-2")
  )
  expect_error(
    ingestion_oracle_generalization_request(hierarchy),
    "requested explicitly"
  )
  run_request <- ingestion_oracle_generalization_request(hierarchy, "run")
  session_request <- ingestion_oracle_generalization_request(
    hierarchy, "session"
  )
  expect_false(identical(run_request, session_request))
})

test_that("fixed and learned observation models earn different guarantees", {
  fixed <- ingestion_oracle_capabilities("fixed_gls")
  learned <- ingestion_oracle_capabilities("learned_ar")

  expect_true(fixed$analytic_effect_covariance)
  expect_true(fixed$fixed_observation_model)
  expect_false(learned$analytic_effect_covariance)
  expect_true(learned$learned_observation_model)
})
