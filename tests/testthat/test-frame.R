grid_domain <- function() {
  coordinates <- as.matrix(expand.grid(x = 0:2, y = 0:2))
  abstract_domain(nrow(coordinates), coordinates = coordinates,
    feature_ids = paste0("v", seq_len(nrow(coordinates))), id = "grid3x3")
}

test_that("all spatial scopes compile to one sparse additive representation", {
  domain <- grid_domain()
  frames <- list(
    point = compile_frame(voxelwise(), domain),
    searchlight = compile_frame(searchlights(1.01), domain),
    region = compile_frame(regions(rep(c("left", "right", "middle"), each = 3)), domain),
    global = compile_frame(whole_brain(), domain)
  )

  for (candidate in frames) {
    expect_s3_class(candidate, "effect_frame")
    expect_true(inherits(candidate$weights, "sparseMatrix"))
    expect_identical(candidate$representation, "additive_diagonal")
    expect_identical(candidate$domain_id, domain$id)
    expect_identical(candidate$domain, domain$reference)
  }
  expect_equal(dim(frames$point$weights), c(9, 9))
  expect_equal(dim(frames$searchlight$weights), c(9, 9))
  expect_equal(dim(frames$region$weights), c(3, 9))
  expect_equal(dim(frames$global$weights), c(1, 9))
})

test_that("local normalization produces unit measurement mass", {
  domain <- grid_domain()
  for (specification in list(searchlights(1.5),
    regions(rep(c("a", "b", "c"), each = 3)), whole_brain())) {
    compiled <- compile_frame(specification, domain)
    expect_equal(Matrix::rowSums(compiled$weights),
      rep(1, nrow(compiled$weights)), tolerance = 1e-14)
  }
})

test_that("conservative frames partition global feature mass", {
  domain <- grid_domain()
  specs <- list(
    voxelwise(),
    searchlights(1.5, normalization = "conservative"),
    regions(rep(c("a", "b", "c"), each = 3), normalization = "conservative"),
    whole_brain(normalization = "conservative")
  )
  for (specification in specs) {
    compiled <- compile_frame(specification, domain)
    expect_equal(Matrix::colSums(compiled$weights), rep(1, 9),
      tolerance = 1e-14)
  }
})

test_that("searchlight neighborhoods respect domain geometry", {
  domain <- grid_domain()
  compiled <- compile_frame(searchlights(1.01, normalization = "none"), domain)
  sizes <- Matrix::rowSums(compiled$weights)

  expect_equal(sizes[c(1, 3, 7, 9)], rep(3, 4))
  expect_equal(sizes[5], 5)
})

test_that("frame compilation rejects mismatches and uncovered conservation", {
  domain <- grid_domain()
  expect_error(compile_frame(regions(c("a", "b")), domain), "one entry")
  expect_error(compile_frame(
    regions(c(rep("a", 8), NA), normalization = "conservative"), domain
  ), "cover every")
  no_coordinates <- abstract_domain(3)
  expect_error(compile_frame(searchlights(1), no_coordinates), "require domain coordinates")

  forged <- searchlights(1)
  forged$invented <- TRUE
  expect_error(compile_frame(forged, domain), "noncanonical")
})

test_that("compiled frames remain fail closed after mutation", {
  compiled <- compile_frame(searchlights(1.5), grid_domain())
  compiled$weights@x[[1]] <- -1
  expect_error(compile_lowering(compiled, bilinear_query(diag(2))), "nonnegative")
})
