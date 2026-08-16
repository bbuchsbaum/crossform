# A validation message is only useful if the reader can act on it, so the
# public entry points must state what was received, what was expected, and the
# call that fixes it. These tests pin that contract on the misuses a first-hour
# user actually commits.

dx_fixture <- function(effects = c("face", "house", "body", "tool"),
                       id = "dx-message-domain") {
  set.seed(4021)
  domain <- abstract_domain(6L, id = id)
  betas <- function() {
    value <- matrix(rnorm(length(effects) * 6L), length(effects), 6L)
    rownames(value) <- effects
    value
  }
  rel <- relation(list(run1 = betas(), run2 = betas()), domain = domain)
  frame <- compile_frame(voxelwise(), domain)
  plan <- plan_geometry(
    rel, frame,
    cross_partitions(rel, independence = "independent",
      generalizes_over = "run")
  )
  list(domain = domain, relation = rel, frame = frame, plan = plan,
    effects = effects)
}

test_that("contrast weight errors name the count and the declared effects", {
  fixture <- dx_fixture()

  expect_error(contrast_energy(fixture$plan, c(1, -1, 0)),
    "has 3 values but the relation declares 4 effects")
  expect_error(contrast_energy(fixture$plan, c(1, -1, 0)),
    "`face`, `house`, `body`, `tool`")

  expect_error(
    contrast_energy(fixture$plan,
      c(face = 1, houses = -1, body = 0, tool = 0)),
    "`houses` is not a declared effect; `house` has no weight")
  expect_error(
    contrast_energy(fixture$plan, c(face = 1, face = -1, body = 0, tool = 0)),
    "names `face` more than once")
  expect_error(
    contrast_energy(fixture$plan,
      c(face = 1, house = NA, body = 0, tool = 0)),
    "the weight for `house` is NA, NaN, or Inf")
  expect_error(contrast_energy(fixture$plan),
    "`weights` is required")
  expect_error(contrast_energy(fixture$relation, c(1, -1, 0, 0)),
    "received an object of class `effect_relation`")
})

test_that("rdm and rsa name the unknown effect and the offending model", {
  fixture <- dx_fixture()

  expect_error(rdm(fixture$plan, pairs = cbind("face", "hous")),
    "`pairs` names `hous`, which is not a declared experimental effect")

  wrong_size <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3L, 3L)
  expect_error(rsa(fixture$plan, models = list(category = wrong_size)),
    "`models` RDM `category` is 3 x 3; the relation declares 4 effects")

  similarity <- matrix(1, 4L, 4L)
  expect_error(rsa(fixture$plan, models = list(category = similarity)),
    "nonzero diagonal")

  mislabelled <- matrix(0, 4L, 4L)
  mislabelled[upper.tri(mislabelled)] <- 1
  mislabelled <- mislabelled + t(mislabelled)
  wrong_names <- c("face", "house", "body", "toolz")
  dimnames(mislabelled) <- list(wrong_names, wrong_names)
  expect_error(rsa(fixture$plan, models = list(category = mislabelled)),
    "`toolz` is not a declared effect")

  expect_error(rsa(fixture$plan), "`models` is required")
})

test_that("rank deficiency names the intercept and the collinear columns", {
  fixture <- dx_fixture()

  constant <- matrix(1, 4L, 4L)
  diag(constant) <- 0
  message <- tryCatch(
    rsa(fixture$plan, models = list(flat = constant)),
    error = conditionMessage
  )
  expect_match(message, "rank 1 of 2 columns")
  expect_match(message, "`flat` is an exact linear combination of `\\(Intercept\\)`")
  expect_match(message, "intercept = FALSE")
  # And the suggested remedy is real: the same models fit without it.
  expect_s3_class(
    rsa(fixture$plan, models = list(flat = constant), intercept = FALSE),
    "effect_rsa_view"
  )

  duplicate <- matrix(0, 4L, 4L)
  duplicate[1L, 2L] <- duplicate[2L, 1L] <- 1
  pair_message <- tryCatch(
    rsa(fixture$plan, models = list(first = duplicate, second = duplicate)),
    error = conditionMessage
  )
  expect_match(pair_message,
    "`second` is an exact linear combination of `first`")
  expect_no_match(pair_message, "intercept = FALSE")
})

test_that("frame and domain errors report the value that was supplied", {
  fixture <- dx_fixture()

  expect_error(compile_frame(voxelwise(), NULL),
    "`domain` must be an `effect_domain`.*received NULL")
  expect_error(compile_frame(voxelwise(), data.frame(a = 1)),
    "received a data frame with 1 row and 1 column")
  expect_error(compile_frame(voxelwise()),
    "received no argument")
  expect_error(compile_frame("voxelwise", fixture$domain),
    "must be a frame specification from `voxelwise\\(\\)`")

  expect_error(searchlights(-2), "received `-2`")
  expect_error(searchlights(c(3, 4)), "received a numeric vector of length 2")
  expect_error(compile_frame(regions(c("a", "b")), fixture$domain),
    "supplied 2 labels but the domain has 6 features")
})

test_that("required plan arguments name the call that supplies them", {
  fixture <- dx_fixture()

  expect_error(plan_geometry(fixture$relation),
    "`at` is required: pass a compiled frame from `compile_frame\\(\\)`")
  expect_error(plan_geometry(fixture$relation, fixture$frame),
    "`over` is required: pass a pairing from `cross_partitions\\(\\)`")
  expect_error(plan_geometry("betas", fixture$frame, fixture$plan$pairing),
    "Expected an `effect_relation`")
  expect_error(plan_geometry(fixture$relation, "voxels", fixture$plan$pairing),
    "Expected a compiled `effect_frame`")
  expect_error(plan_geometry(fixture$relation, fixture$frame, "runs"),
    "Expected an `effect_pairing`")
})

test_that("relation and pairing errors quote the partitions involved", {
  fixture <- dx_fixture()

  single <- relation(
    list(run1 = matrix(1, 4L, 6L,
      dimnames = list(fixture$effects, NULL))),
    domain = fixture$domain
  )
  expect_error(cross_partitions(single, independence = "independent"),
    "this relation has 1 partition \\(`run1`\\)")

  contaminated <- matrix(rnorm(24), 4L, 6L,
    dimnames = list(fixture$effects, NULL))
  contaminated[2L, 3L] <- NA
  expect_error(
    relation(list(run1 = contaminated), domain = fixture$domain),
    "Partition `run1`: 1 of 24 values is non-finite"
  )
  expect_error(
    relation(list(run1 = contaminated), domain = fixture$domain),
    "in row `house`"
  )
})

test_that("as_neurovol explains the measurement-versus-feature distinction", {
  skip_if_not_installed("neuroim2")
  skip_if_not(utils::packageVersion("neuroim2") >= "0.19.0")
  values <- array(FALSE, c(5L, 5L, 4L))
  values[2:4, 2:4, 2:3] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(
    values, neuroim2::NeuroSpace(c(5L, 5L, 4L), spacing = c(3, 3, 3))
  )
  domain <- neuroim2_volume_domain(mask, id = "dx-message-volume")

  message <- tryCatch(as_neurovol(c(1, 2, 3), mask, domain),
    error = conditionMessage)
  expect_match(message, "has 3 values but domain `dx-message-volume` has 18")
  expect_match(message, "one value per \\*measurement\\*")
  expect_match(message, "one value per \\*feature\\*")
  expect_match(message, "Matrix::crossprod\\(frame\\$weights != 0, values\\)")
  expect_match(message, "neuroim2-data")

  # The documented remedy is the one the vignette teaches, and it works.
  regional <- compile_frame(
    regions(rep(c("a", "b", "c"), length.out = domain$n_features)), domain
  )
  regional_values <- c(1, 2, 3)
  expanded <- as.numeric(
    Matrix::crossprod(regional$weights != 0, regional_values)
  )
  expect_length(expanded, domain$n_features)
  expect_s4_class(as_neurovol(expanded, mask, domain), "NeuroVol")
})

test_that("aligned contrast weights carry the relation's effect names", {
  fixture <- dx_fixture()
  positional <- contrast_energy(fixture$plan, c(1, -1, 0, 0))
  expect_identical(names(positional$weights), fixture$effects)

  named <- contrast_energy(fixture$plan,
    c(tool = 0, house = -1, face = 1, body = 0))
  expect_identical(named$weights, positional$weights)
  expect_equal(named$total, positional$total, tolerance = 0)
})

test_that("an out-of-range measurement index names the argument and range", {
  example <- example_fmri_effects()
  plan <- plan_geometry(
    example$fit$relation, example$frame,
    cross_partitions(example$fit$relation, independence = "independent")
  )
  measurements <- plan$measurements

  expect_error(
    rdm_sampling_covariance(plan, example$fit, target = "null", at = 99999),
    sprintf("`at` = 99999 is outside the plan's 1\\.\\.%d measurements",
      measurements)
  )
  expect_error(
    rdm_sampling_covariance(plan, example$fit, target = "null", at = 0),
    sprintf("`at` = 0 is outside the plan's 1\\.\\.%d measurements",
      measurements)
  )
  expect_error(
    rdm_sampling_covariance(plan, example$fit, target = "null", at = "peak"),
    sprintf(
      "`at` must be one measurement index in 1\\.\\.%d; received a character",
      measurements)
  )

  # The compiled accessor reports the same way when it is reached directly.
  read_node <- crossform:::.frame_metric_node_accessor(example$frame)
  expect_error(read_node(1e6),
    sprintf("`at` = 1000000 is outside the frame's 1\\.\\.%d measurements",
      measurements))
})
