# A view prints its frame's `$index$measurement` under `measurement`. On a
# masked volume those are full-volume voxel indices, not row positions, and a
# bare number copied from a print used to select the row at that position --
# a different voxel, returned without complaint. These tests pin the fix: a
# bare number is refused wherever it could mean either, and
# `measurement_id()` / `measurement_position()` (or a character id) say which.

masked_selector_fixture <- function(holes = TRUE) {
  set.seed(4417)
  mask <- array(TRUE, c(4L, 4L, 2L))
  if (holes) {
    mask[1:2, , 1] <- FALSE
    mask[4, 4, 2] <- FALSE
  }
  domain <- volume_domain(mask, spacing = c(3, 3, 3),
    id = if (holes) "selector-masked" else "selector-full")
  n <- domain$n_features
  observations <- 32L
  condition <- factor(rep(c("a", "b", "c"), length.out = observations))
  design <- stats::model.matrix(~ 0 + condition)
  colnames(design) <- c("a", "b", "c")
  effects <- diag(3)
  rownames(effects) <- colnames(design)
  sources <- stats::setNames(lapply(seq_len(3L), function(run) {
    design %*% matrix(stats::rnorm(3L * n, sd = 0.3), 3L, n) +
      matrix(stats::rnorm(observations * n), observations, n)
  }), paste0("run", seq_len(3L)))
  fit <- lm_relation_fit(sources, design, effects, sampling_unit = "trial",
    domain = domain)
  frame <- compile_frame(searchlights(3.1), domain)
  plan <- plan_geometry(fit$relation, frame,
    cross_partitions(fit$relation, independence = "independent"),
    metric = noise_precision(diag(n), domain, covariance = diag(n),
      provenance = list(source = "selector-fixture")))
  list(domain = domain, fit = fit, frame = frame, plan = plan)
}

test_that("a masked frame's printed identifiers are not its row positions", {
  fixture <- masked_selector_fixture()
  ids <- fixture$frame$index$measurement
  expect_identical(ids, fixture$domain$feature_ids)
  expect_false(identical(as.integer(ids), seq_along(ids)))
  view <- rdm(fixture$plan)
  # The view carries the identifiers as a plain vector, the frame's column.
  expect_identical(view$index, ids)
  expect_output(print(view), "\\b12\\b")
})

test_that("rdm_sampling_covariance() refuses a bare number on a masked frame", {
  fixture <- masked_selector_fixture()
  ids <- fixture$frame$index$measurement
  id <- 12L
  position <- match(id, ids)
  # Both readings are valid measurements, and they are different ones: the
  # case that used to return the wrong voxel silently.
  expect_true(id <= length(ids) && ids[[id]] != id)

  condition <- tryCatch(
    rdm_sampling_covariance(fixture$plan, fixture$fit, target = "null",
      at = id),
    error = identity)
  expect_s3_class(condition, "effect_input_error")
  expect_match(conditionMessage(condition), "ambiguous", fixed = TRUE)
  expect_match(conditionMessage(condition), "measurement_id(12)", fixed = TRUE)
  expect_match(conditionMessage(condition), "measurement_position(12)",
    fixed = TRUE)

  by_id <- rdm_sampling_covariance(fixture$plan, fixture$fit, target = "null",
    at = measurement_id(id))
  by_string <- rdm_sampling_covariance(fixture$plan, fixture$fit,
    target = "null", at = "12")
  by_position <- rdm_sampling_covariance(fixture$plan, fixture$fit,
    target = "null", at = measurement_position(position))
  expect_identical(by_id$source$node, "12")
  expect_identical(by_string$source$node, "12")
  expect_identical(by_position$source$node, "12")
  expect_equal(sampling_covariance(by_id), sampling_covariance(by_position),
    tolerance = 1e-12)

  # The other reading names a different voxel; it is reachable, but only by
  # saying so.
  other <- rdm_sampling_covariance(fixture$plan, fixture$fit, target = "null",
    at = measurement_position(id))
  expect_identical(other$source$node, as.character(ids[[id]]))
  expect_false(identical(other$source$node, "12"))

  batch <- rdm_sampling_covariance(fixture$plan, fixture$fit, target = "null",
    at = measurement_id(c(12, 17)))
  expect_identical(vapply(batch, function(x) x$source$node, character(1)),
    c("12", "17"))
  expect_error(
    rdm_sampling_covariance(fixture$plan, fixture$fit, target = "null",
      at = c(12L, 17L)),
    "ambiguous", class = "effect_input_error")
  expect_error(
    rdm_sampling_covariance(fixture$plan, fixture$fit, target = "null",
      at = measurement_id(c(12, 12))),
    "`12` repeated", class = "effect_input_error")
  expect_error(
    rdm_sampling_covariance(fixture$plan, fixture$fit, target = "null",
      at = measurement_id(999)),
    "does not carry: `999`", class = "effect_input_error")
  expect_error(
    rdm_sampling_covariance(fixture$plan, fixture$fit, target = "null",
      at = measurement_position(length(ids) + 1L)),
    "outside the plan's", class = "effect_input_error")
})

test_that("plot() reads measurement and highlight the same way", {
  fixture <- masked_selector_fixture()
  ids <- fixture$frame$index$measurement
  view <- rdm(fixture$plan)
  energy <- contrast_energy(fixture$plan, c(a = 1, b = -1, c = 0))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_error(plot(view, measurement = 12), "ambiguous",
    class = "effect_input_error")
  expect_error(plot(energy, highlight = c(12, 17)), "ambiguous",
    class = "effect_input_error")
  expect_silent(plot(view, measurement = measurement_id(12)))
  expect_silent(plot(view, measurement = "12"))
  expect_silent(plot(view, measurement = measurement_position(1L)))
  # The default measurement is a position the package chose itself.
  expect_silent(plot(view))
  expect_silent(plot(energy, highlight = measurement_id(c(12, 17))))
  expect_silent(plot(energy, highlight = ids == 12L))
  expect_error(plot(view, measurement = measurement_id(999)),
    "does not identify a measurement", class = "effect_input_error")
})

test_that("bare numbers stay positions where identifiers equal positions", {
  fixture <- masked_selector_fixture(holes = FALSE)
  ids <- fixture$frame$index$measurement
  expect_identical(as.integer(ids), seq_along(ids))
  bare <- rdm_sampling_covariance(fixture$plan, fixture$fit, target = "null",
    at = 5L)
  wrapped <- rdm_sampling_covariance(fixture$plan, fixture$fit,
    target = "null", at = measurement_id(5))
  expect_identical(bare$source$node, "5")
  expect_identical(wrapped$source$node, "5")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(rdm(fixture$plan), measurement = 5L))
})

test_that("non-numeric identifiers never shadow positions", {
  expect_false(crossform:::.msg_ids_shadow_positions(c("V1", "V2"), 2L))
  expect_false(crossform:::.msg_ids_shadow_positions(1:3, 3L))
  expect_true(crossform:::.msg_ids_shadow_positions(c(3L, 1L, 2L), 3L))
  expect_true(crossform:::.msg_ids_shadow_positions(c(1e5, 2e5), 2L))
  expect_identical(crossform:::.msg_id_key(c(1e5, 7)), c("100000", "7"))
})

test_that("the selector constructors validate their input", {
  expect_s3_class(measurement_id(c(3, 4)), "crossform_measurement_id")
  expect_s3_class(measurement_id("V1"), "crossform_measurement_id")
  expect_s3_class(measurement_position(2L), "crossform_measurement_position")
  expect_output(print(measurement_id(117)), "<measurement_id> 117")
  expect_output(print(measurement_position(5)), "<measurement_position> 5")
  expect_error(measurement_id(NA), class = "effect_input_error")
  expect_error(measurement_id(list(1)), class = "effect_input_error")
  expect_error(measurement_id(""), class = "effect_input_error")
  expect_error(measurement_position(1.5), class = "effect_input_error")
  expect_error(measurement_position("1"), class = "effect_input_error")
})

# Item: a named contrast that leaves effects out is refused with the complete
# vector to paste back, never defaulted to zero.
test_that("a partially named contrast names the missing effects", {
  effects <- c("face", "house", "body", "tool")
  condition <- tryCatch(
    crossform:::.align_contrast(c(face = 1, house = -1), effects),
    error = identity)
  expect_s3_class(condition, "effect_input_error")
  message <- conditionMessage(condition)
  expect_match(message, "`body`, `tool` have no weight", fixed = TRUE)
  expect_match(message, "c(face = 1, house = -1, body = 0, tool = 0)",
    fixed = TRUE)
  expect_false(grepl("name the weights to have them aligned", message))
  # Still a refusal, not a silent zero.
  expect_error(crossform:::.align_contrast(c(face = 1, house = -1), effects),
    class = "effect_input_error")
  expect_equal(
    crossform:::.align_contrast(
      c(tool = 0, body = 0, house = -1, face = 1), effects),
    c(face = 1, house = -1, body = 0, tool = 0))
  expect_error(crossform:::.align_contrast(c(1, -1), effects),
    "c(face = 1, house = -1, body = 0, tool = 0)", fixed = TRUE,
    class = "effect_input_error")
})

test_that("contrast_energy() reports the missing effects end to end", {
  fixture <- masked_selector_fixture()
  expect_error(contrast_energy(fixture$plan, c(a = 1, b = -1)),
    "`c` has no weight.*c\\(a = 1, b = -1, c = 0\\)",
    class = "effect_input_error")
})

# Item: a data-frame source used to report "has  rows" with a blank count.
test_that("relation() names a data-frame source instead of a blank count", {
  table <- data.frame(v1 = c(1, 2), v2 = c(3, 4),
    row.names = c("face", "house"))
  for (call in list(
    quote(relation(list(r1 = table, r2 = table),
      effects = effect_space(c("face", "house")))),
    quote(relation(list(r1 = table, r2 = table)))
  )) {
    condition <- tryCatch(eval(call), error = identity)
    expect_s3_class(condition, "effect_input_error")
    expect_false(grepl("has  ", conditionMessage(condition), fixed = TRUE))
    expect_match(conditionMessage(condition), "data frame")
    expect_match(conditionMessage(condition), "as.matrix()", fixed = TRUE)
  }
  converted <- relation(list(r1 = as.matrix(table), r2 = as.matrix(table)))
  expect_s3_class(converted, "effect_relation")
  expect_error(
    relation(list(r1 = function(features) NULL)),
    "source_dims", class = "effect_input_error")
})
