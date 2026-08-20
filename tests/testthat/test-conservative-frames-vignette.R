# `vignettes/conservative-frames.Rmd` -- ticket D10.
#
# The vignette is `design/conservative-geometry-contract.md` made runnable, and
# its assertions are its content: every identity is asserted in *visible* code,
# so the article stops knitting the day one of them stops holding. That much is
# enforced by `R CMD build`. What a build cannot notice is a chunk quietly
# going away -- an article that no longer demonstrates the smoothed-ledger
# caveat still knits perfectly. This file pins the demonstrations themselves,
# and it recomputes nothing.

conservative_vignette_path <- function(name = "conservative-frames.Rmd") {
  # `R CMD build` keeps the vignette source in `inst/doc`, so the installed
  # package carries it under check. A source checkout has neither, and reads
  # `vignettes/` directly.
  installed <- tryCatch(
    system.file("doc", name, package = "crossform"),
    error = function(condition) ""
  )
  candidates <- c(
    if (nzchar(installed)) installed,
    testthat::test_path("..", "..", "vignettes", name),
    testthat::test_path("..", "..", "inst", "doc", name)
  )
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(candidate)
    }
  }
  NA_character_
}

conservative_vignette_text <- function() {
  path <- conservative_vignette_path()
  if (is.na(path)) {
    skip(paste(
      "The conservative-frames vignette source is not available here:",
      "neither system.file('doc') nor the checkout's vignettes/ has it."
    ))
  }
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

test_that("the conservative-frames guide contrasts both instruments", {
  text <- conservative_vignette_text()

  # Section 1. Detection and attribution built from the same fit, and the
  # comparator taken under the *unnormalized* whole-domain operator -- the
  # default `"local"` divides by the feature count and manufactures a spurious
  # failure (contract section 2, precondition 1).
  expect_match(text, 'searchlights(4, "local")', fixed = TRUE)
  expect_match(text, 'searchlights(4, "conservative")', fixed = TRUE)
  expect_match(text, 'whole_brain("none")', fixed = TRUE)
  expect_match(text, "frame_conservation(", fixed = TRUE)

  # The refusal, exercised rather than described: a detection map has no
  # budget, so `contribution()` refuses it by capability.
  expect_match(text, "conservative_frame", fixed = TRUE)
})

test_that("the guide demonstrates the smoothed ledger and the alpha law", {
  text <- conservative_vignette_text()

  # Section 2, claim 3: the conservative total is the frame contraction of a
  # per-voxel ledger, which `voxelwise()` supplies as a point frame. Dropping
  # this chunk would leave the "no cross-voxel information" caveat asserted
  # only in prose.
  expect_match(text, "voxelwise()", fixed = TRUE)
  expect_match(text, "attribution$weights %*% ledger$total", fixed = TRUE)

  # Section 3.1: per-scale totals are alpha_s times the whole-domain total, so
  # a multiscale energy panel is a picture of `weights`. The multiscale route
  # is a radius vector, and the prohibition is enforced, not documented --
  # `plot(spectrum, which = "profile")` must be *exercised* here.
  expect_match(text, 'searchlights(radii, "conservative")', fixed = TRUE)
  expect_match(text, 'plot(spectrum, which = "profile")', fixed = TRUE)
  expect_match(text, "scale_energy_panel", fixed = TRUE)
})

test_that("the guide measures alpha-invariance rather than citing it", {
  text <- conservative_vignette_text()

  # Section 3.2 is the reason the spectrum is well posed, and the vignette has
  # to *rerun* the family under a lopsided weighting to show the energy moving
  # while the share does not.
  expect_match(text, "weights = c(0.6, 0.3, 0.1)", fixed = TRUE)
  expect_match(text, "coherence_spectrum(", fixed = TRUE)
  expect_match(text, "by_location = TRUE", fixed = TRUE)
  expect_match(text, "alpha_invariant", fixed = TRUE)
  expect_match(text, "alpha_fixed", fixed = TRUE)
})

test_that("the guide keeps the ledger and latent-layer disciplines", {
  text <- conservative_vignette_text()

  # Section 4: budget exactness, and the two things a territory ledger refuses
  # to give -- a summed density and a frame-independent coherent budget.
  expect_match(text, "contribution(attribution_map, by = territory)",
    fixed = TRUE)
  expect_match(text, "budget_exact", fixed = TRUE)
  expect_match(text, "frame_relative", fixed = TRUE)
  expect_match(text, "masked", fixed = TRUE)

  # Section 5: the signed estimation layer versus the declared projection.
  expect_match(text, "latent_geometry(", fixed = TRUE)
  expect_match(text, "psd_projection", fixed = TRUE)
  expect_match(text, "moved_share", fixed = TRUE)
  expect_match(text, "n_eff", fixed = TRUE)
  expect_match(text, "latent_projection_source", fixed = TRUE)
  expect_match(text, "not for inference", fixed = TRUE)
})

test_that("the guide keeps the metric composition demonstration", {
  text <- conservative_vignette_text()

  # Section 6: diagonal metrics fold and conserve, dense ones break native
  # conservation, and `"whitened"` restores it as a *different estimand* whose
  # difference has to enter plan identity.
  expect_match(text, "metric_capabilities(", fixed = TRUE)
  expect_match(text, "feature_additive", fixed = TRUE)
  expect_match(text, 'composition = "whitened"', fixed = TRUE)
  expect_match(text, "scientific_plan_id", fixed = TRUE)

  # The dense failure is asserted as "not small", never as an equality: its
  # size and sign are fixture specific (contract section 5.1 and section 9).
  expect_match(text, "> 0.01", fixed = TRUE)
})

test_that("the guide asserts its identities in visible code", {
  text <- conservative_vignette_text()

  # The article's premise is that a broken identity breaks the knit. That only
  # holds if the assertions are in evaluated chunks: an `include = FALSE`
  # check still runs, but a chunk marked `eval = FALSE` does not, and a reader
  # cannot see either. Both halves are pinned -- the assertions exist, and no
  # chunk in this article opts out of evaluation.
  expect_gte(lengths(regmatches(text, gregexpr("stopifnot(", text,
    fixed = TRUE)))[[1L]], 12L)
  expect_false(grepl("eval = FALSE", text, fixed = TRUE))
  expect_false(grepl("eval=FALSE", text, fixed = TRUE))

  # Contract section 9's tolerance for the conservation laws. A vignette that
  # relaxed these to `1e-6` would still knit and would still look like a
  # demonstration.
  expect_match(text, "1e-12", fixed = TRUE)

  # Tier discipline: this is an advanced-tier article, so it may reach past the
  # core views -- but not into the package's internals.
  expect_false(grepl("crossform:::", text, fixed = TRUE))
})
