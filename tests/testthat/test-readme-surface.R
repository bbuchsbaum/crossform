# The README is the first hour. Ticket A9's rule is that the first hour costs
# the reader exactly one tier: every crossform function the README *executes*
# must be core or core-ingestion, so a newcomer who works through it has met
# the entry price and nothing beyond it. Advanced, experimental-coupling, and
# adapter material may be named in prose (`?pairing`, `measurement_form()`,
# `neuroim2_volume_domain()`), and prose is not what this file reads: it parses
# the fenced ```r chunks and looks only at names in call position.
#
# This test never *runs* README code. Executing it needs a fitted relation, a
# volume, and several seconds; the vignettes and exemplars carry that burden.
# Here the README is treated as a document whose R is parsed for its API
# surface.
#
# The two tier vectors below are transcribed from `design/api-tiers.md` (the
# tier ledger, sections "Tier: core (21)" and "Tier: core-ingestion (20)").
# They are embedded rather than parsed out of the ledger so that this file
# states the promise it enforces. If this test fails, a README example has
# reached outside core; widening the vectors here is therefore a *deliberate
# edit* — it re-tiers a function in the ledger's terms, and the ledger must be
# edited first. Do not add a name to make a failure go away.

readme_core_tier <- c(
  "abstract_domain", "catch_refusal", "compile_frame", "compute_policy",
  "contrast_energy", "cross_partitions", "effect_space",
  "example_fmri_effects", "lm_relation_fit", "plan_geometry", "rdm",
  "rdm_sampling_covariance", "regions", "relation", "rsa",
  "sampling_capabilities", "sampling_covariance", "searchlights",
  "volume_domain", "voxelwise", "whole_brain"
)

readme_core_ingestion_tier <- c(
  "coefficient_parameterization", "compiler_conformance", "condition_space",
  "design_model", "effect_map", "estimate_relation", "lower_effect_map",
  "observation_confounds", "observation_events", "observation_index",
  "observation_model", "observations", "partition_hierarchy",
  "plan_relation", "raw_design_model", "raw_effect_map",
  "relation_plan_receipts", "study", "study_axis", "study_capabilities"
)

# Fenced R chunks only: ```r ... ``` (and ```{r ...} if the README ever moves
# to knitr). A ```sh or ```yaml fence is skipped, as is anything outside a
# fence.
readme_r_chunks <- function(lines) {
  opens <- grepl("^\\s*```\\s*[{]?r[},[:space:]]*$", lines)
  fences <- grepl("^\\s*```", lines)
  chunks <- list()
  i <- 1L
  n <- length(lines)
  while (i <= n) {
    if (opens[[i]]) {
      j <- i + 1L
      while (j <= n && !fences[[j]]) j <- j + 1L
      if (j > i + 1L) chunks[[length(chunks) + 1L]] <- lines[(i + 1L):(j - 1L)]
      i <- j + 1L
    } else {
      i <- i + 1L
    }
  }
  chunks
}

# Names in call position, recursively. Comments (`#>` output blocks included)
# never reach here because `parse()` has already dropped them, and a
# namespace-qualified call such as `Matrix::rowSums()` is reported as
# `rowSums`, not as a bare crossform name.
called_names <- function(expr) {
  if (is.call(expr)) {
    head <- expr[[1L]]
    here <- if (is.symbol(head)) as.character(head) else character(0)
    if (is.call(head) && is.symbol(head[[1L]]) &&
        as.character(head[[1L]]) %in% c("::", ":::")) {
      here <- character(0)
    }
    rest <- unlist(lapply(as.list(expr)[-1L], called_names), use.names = FALSE)
    c(here, if (is.call(head)) called_names(head), rest)
  } else if (is.pairlist(expr) || is.expression(expr) || is.list(expr)) {
    unlist(lapply(as.list(expr), called_names), use.names = FALSE)
  } else {
    character(0)
  }
}

readme_called_names <- function(path) {
  chunks <- readme_r_chunks(readLines(path, warn = FALSE))
  out <- character(0)
  for (k in seq_along(chunks)) {
    parsed <- tryCatch(
      parse(text = chunks[[k]], keep.source = FALSE),
      error = function(e) {
        testthat::fail(sprintf(
          "README R chunk %d does not parse: %s", k, conditionMessage(e)
        ))
        NULL
      }
    )
    if (!is.null(parsed)) {
      out <- c(out, unlist(lapply(as.list(parsed), called_names),
                           use.names = FALSE))
    }
  }
  list(chunks = chunks, names = sort(unique(out)))
}

readme_path <- testthat::test_path("..", "..", "README.md")

test_that("the README's R chunks parse and are actually found", {
  skip_if_not(file.exists(readme_path), "README.md not available")
  found <- readme_called_names(readme_path)

  # A regex that silently stops matching would turn this whole file into a
  # no-op that passes, so pin the extraction itself.
  expect_gt(length(found$chunks), 5L)
  expect_true("plan_geometry" %in% found$names)
  expect_true("contrast_energy" %in% found$names)

  # The shell fence must not be swept up with the R ones.
  expect_false("R" %in% found$names)
})

test_that("every crossform function the README executes is core-tier", {
  skip_if_not(file.exists(readme_path), "README.md not available")
  found <- readme_called_names(readme_path)

  exports <- getNamespaceExports("crossform")
  executed <- sort(intersect(found$names, exports))
  first_hour <- c(readme_core_tier, readme_core_ingestion_tier)

  outside <- setdiff(executed, first_hour)
  expect_identical(
    outside, character(0),
    info = paste0(
      "README chunks call non-core crossform exports: ",
      paste(outside, collapse = ", "),
      ". Move the example behind a vignette, or re-tier the function in ",
      "design/api-tiers.md and update this test deliberately."
    )
  )

  # The README is supposed to *demonstrate* the core spine, not merely avoid
  # leaving it, so require that it exercises a real slice of it.
  expect_gt(length(executed), 10L)
})

test_that("the embedded tier vectors still match the ledger's shape", {
  # Transcription guard. The ledger states 21 core and 20 core-ingestion
  # exports; both sets must be exported, and they must not overlap.
  expect_identical(length(readme_core_tier), 21L)
  expect_identical(length(readme_core_ingestion_tier), 20L)
  expect_identical(anyDuplicated(readme_core_tier), 0L)
  expect_identical(anyDuplicated(readme_core_ingestion_tier), 0L)
  expect_identical(
    intersect(readme_core_tier, readme_core_ingestion_tier), character(0)
  )

  exports <- getNamespaceExports("crossform")
  expect_identical(setdiff(readme_core_tier, exports), character(0))
  expect_identical(setdiff(readme_core_ingestion_tier, exports), character(0))
})
