make_neuroim2_mask <- function(dims = c(5L, 5L, 4L),
                               spacing = c(1, 2, 4), active = NULL) {
  testthat::skip_if_not_installed("neuroim2", minimum_version = "0.19.0")
  values <- array(FALSE, dims)
  if (is.null(active)) values[2:4, 2:4, 2:3] <- TRUE else values[active] <- TRUE
  neuroim2::LogicalNeuroVol(values,
    neuroim2::NeuroSpace(dims, spacing = spacing))
}

test_that("neuroim2 volume domains preserve full indices and spacing", {
  mask <- make_neuroim2_mask(active = c(1L, 2L, 6L, 55L))
  domain <- neuroim2_volume_domain(mask, id = "subject-mask")

  expect_identical(domain$feature_ids, which(mask != 0))
  expect_identical(domain$metadata$spacing, c(1, 2, 4))
  expect_equal(domain$coordinates,
    sweep(neuroim2::index_to_grid(mask, domain$feature_ids) - 1, 2,
      c(1, 2, 4), `*`), tolerance = 0)
  expect_match(domain$metadata$neuroim2_space_sha256,
    "^sha256:[[:xdigit:]]{64}$")
})

test_that("neuroim2 neighborhoods map exactly into compact frame columns", {
  mask <- make_neuroim2_mask()
  domain <- neuroim2_volume_domain(mask)
  indices <- neuroim2::searchlight_indices(mask, radius = 4, nonzero = TRUE)
  compiled <- neuroim2_searchlights(mask, radius = 4, domain = domain,
    normalization = "none")

  expect_identical(compiled$domain, domain$reference)
  expect_identical(compiled$index$measurement,
    attr(indices, "center_indices"))
  for (center in seq_along(indices)) {
    expect_setequal(which(as.matrix(compiled$weights)[center, ] != 0),
      match(indices[[center]], domain$feature_ids))
  }
})

test_that("adapter preserves boundaries, anisotropy, and future plan", {
  testthat::skip_if_not_installed("future")
  mask <- make_neuroim2_mask(dims = c(4L, 5L, 3L),
    active = c(1L, 2L, 5L, 6L, 21L))
  before <- future::plan()
  frame <- neuroim2_searchlights(mask, radius = 4, normalization = "none")
  after <- future::plan()

  expect_identical(after, before)
  expect_true(all(frame$index$measurement %in% which(mask != 0)))
  expect_true(all(Matrix::rowSums(frame$weights) >= 1))
  expect_identical(frame$specification$upstream_commit, "77b1ddb")
})

test_that("adapter rejects incompatible domains and unrestricted members", {
  mask <- make_neuroim2_mask()
  changed <- make_neuroim2_mask(spacing = c(2, 2, 4))
  domain <- neuroim2_volume_domain(mask)

  expect_error(neuroim2_searchlights(changed, 4, domain = domain),
    "different volume geometry", class = "effect_contract_error")
  expect_error(neuroim2_searchlights(mask, 4, domain = domain,
    nonzero = FALSE), "require.*TRUE", class = "effect_input_error")
})

test_that("compact result values map to exact NeuroVol indices", {
  mask <- make_neuroim2_mask(active = c(1L, 2L, 6L, 55L))
  domain <- neuroim2_volume_domain(mask, id = "result-mask")
  values <- c(-1.5, 0, 2.25, 7)
  result <- as_neurovol(values, mask, domain, fill = NA_real_,
    label = "signed crossnobis")
  array_value <- as.array(result)

  expect_s4_class(result, "NeuroVol")
  expect_identical(dim(result), dim(mask))
  expect_identical(neuroim2::space(result), neuroim2::space(mask))
  expect_equal(array_value[domain$feature_ids], values, tolerance = 0)
  expect_true(all(is.na(array_value[-domain$feature_ids])))

  expect_error(as_neurovol(values[-1L], mask, domain),
    "has 3 values but domain `result-mask` has 4 features",
    class = "effect_input_error")
  expect_error(as_neurovol(values[-1L], mask, domain),
    "one value per \\*measurement\\*", class = "effect_input_error")
  expect_error(as_neurovol(replace(values, 1L, NA_real_), mask, domain),
    "1 of 4 are NA, NaN, or Inf", class = "effect_input_error")
  changed <- make_neuroim2_mask(spacing = c(2, 2, 4),
    active = c(1L, 2L, 6L, 55L))
  expect_error(as_neurovol(values, changed, domain),
    "different volume geometry", class = "effect_contract_error")
})

test_that("missing upstream index API fails with an actionable message", {
  testthat::local_mocked_bindings(
    .require_neuroim2_searchlight_indices = function() {
      stop(paste0("The installed neuroim2 does not provide ",
        "`searchlight_indices()`; install neuroim2 0.19.0 or later ",
        "(upstream commit 77b1ddb)."), call. = FALSE)
    },
    .package = "crossform"
  )
  expect_error(neuroim2_volume_domain(structure(list(), class = "NeuroVol")),
    "0.19.0.*77b1ddb")
})

# Registration into a loaded namespace is what an extension package's
# `S3method()` directive does at load time, so the test does it the same way
# and then undoes it. `on.exit()` lives in a real function frame here because
# it would not fire reliably from inside the `test_that()` expression.
with_registered_as_neurovol <- function(class, method, code) {
  namespace <- asNamespace("crossform")
  registerS3method("as_neurovol", class, method, envir = namespace)
  table <- get(".__S3MethodsTable__.", envir = namespace)
  on.exit(rm(list = paste0("as_neurovol.", class), envir = table), add = TRUE)
  force(code)
}

test_that("as_neurovol dispatches to a method registered by another package", {
  # No neuroim2 skip: the generic validates nothing, so an extension type that
  # brings its own writer never touches the suggested package through us.
  method <- function(values, mask, domain = NULL, fill = NA_real_,
                     label = "crossform result", ...) {
    list(dispatched = TRUE, mask = mask, label = label, extra = list(...))
  }
  thing <- structure(list(), class = "crossform_test_thing")
  result <- with_registered_as_neurovol("crossform_test_thing", method,
    as_neurovol(thing, mask = "extension mask", label = "third party",
      units = "z"))

  expect_true(result$dispatched)
  expect_identical(result$mask, "extension mask")
  expect_identical(result$label, "third party")
  expect_identical(result$extra, list(units = "z"))

  # The registration is the only reason dispatch worked, and it is gone again.
  table <- get(".__S3MethodsTable__.", envir = asNamespace("crossform"))
  expect_false(exists("as_neurovol.crossform_test_thing", envir = table,
    inherits = FALSE))
  # Which refusal comes back depends on whether neuroim2 is installed, so only
  # the fact that the extension method no longer answers is asserted here.
  expect_error(as_neurovol(thing, mask = "extension mask"),
    class = "effect_input_error")
})

test_that("as_neurovol finds a method defined where the generic is called", {
  as_neurovol.crossform_local_thing <- function(values, ...) "dispatched"
  local_thing <- structure(list(), class = "crossform_local_thing")

  expect_identical(as_neurovol(local_thing), "dispatched")
})

test_that("as_neurovol refuses input that no method knows how to write", {
  mask <- make_neuroim2_mask(active = c(1L, 2L, 6L, 55L))
  domain <- neuroim2_volume_domain(mask, id = "result-mask")

  expect_error(as_neurovol(letters[1:4], mask, domain),
    "`values` must be a numeric vector", class = "effect_input_error")
  expect_error(
    as_neurovol(structure(list(1, 2, 3, 4), class = "crossform_unregistered"),
      mask, domain),
    "`values` must be a numeric vector", class = "effect_input_error")
  expect_error(as_neurovol(matrix(1, 4L, 2L), mask, domain),
    "select the column you want to write out", class = "effect_input_error")
  expect_error(as_neurovol(c(-1.5, 0, 2.25, 7), mask, domain, bogus = 1),
    "takes only", class = "effect_input_error")

  # A classed numeric vector was writable before the function became generic
  # and still is: the default method accepts anything numeric.
  labelled <- structure(c(-1.5, 0, 2.25, 7), class = "crossform_labelled")
  expect_s4_class(as_neurovol(labelled, mask, domain), "NeuroVol")
})

test_that("as_neurovol still names both required arguments", {
  expect_error(as_neurovol(), "`values` and `mask` are both required",
    class = "effect_input_error")
  expect_error(as_neurovol(c(1, 2, 3)),
    "`values` and `mask` are both required", class = "effect_input_error")
  expect_error(as_neurovol(structure(list(), class = "crossform_unregistered")),
    "`values` and `mask` are both required", class = "effect_input_error")
})
