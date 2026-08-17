# The package is a layered compiler: primitives -> values -> plans ->
# compiler/execution -> results/views -> adapters. A file may call downward or
# sideways, never upward. The rule is enforced here rather than by convention,
# because 91% of cross-file calls in this package are to dot-internals and a
# single upward call is enough to turn the whole tree into one strongly
# connected component again.
#
# `allowed_upward` below is the debt register: every entry is a real violation
# that exists today, named with the follow-up that removes it. Adding an entry
# is a deliberate act; the test fails on any upward edge that is not listed,
# and also fails on a listed edge that no longer exists, so the register
# cannot rot.
#
# See design/architecture.md for the prose version.

layer_of <- c(
  ## 1. primitives — pure leaves, no crossform dependency outside this layer.
  ## `RcppExports.R` is Rcpp-generated `.Call()` glue: it depends on nothing in
  ## the package, so it is a leaf by construction.
  "primitives.R" = 1L, "message-helpers.R" = 1L, "conditions.R" = 1L,
  "check.R" = 1L, "RcppExports.R" = 1L,

  ## 2. values — domain, frame, pairing, relation, effect space, metric values
  "domain.R" = 2L, "frame.R" = 2L, "pairing.R" = 2L, "relation.R" = 2L,
  "effect-space.R" = 2L, "effect-map.R" = 2L, "metric.R" = 2L,
  "metric-learning.R" = 2L, "source.R" = 2L, "study.R" = 2L,
  "study-facts.R" = 2L, "design-model.R" = 2L, "observation-model.R" = 2L,
  "extractor.R" = 2L, "scope.R" = 2L, "support-index.R" = 2L,
  "numerics.R" = 2L, "query-structured.R" = 2L, "pair-query.R" = 2L,
  "operations.R" = 2L, "receipt.R" = 2L, "reliability.R" = 2L,
  "compute-policy.R" = 2L, "memory-plan.R" = 2L, "capabilities.R" = 2L,
  "validation-memo.R" = 2L, "measurement.R" = 2L,
  "measurement-storage.R" = 2L, "relation-fit.R" = 2L,
  "residual-statistics.R" = 2L, "bridge.R" = 2L, "crossform-package.R" = 2L,

  ## 3. plans — scientific estimands and their identities
  "geometry-plan.R" = 3L, "relation-plan.R" = 3L, "crossnobis.R" = 3L,
  "evidence-task.R" = 3L, "evidence-sampling.R" = 3L,
  "evidence-sampling-kernel.R" = 3L, "evidence-sampling-product.R" = 3L,
  "compiler-conformance.R" = 3L,

  ## 4. compiler / execution, and the result records execution produces.
  ## `result.R` holds the sealed `effect_form` / `effect_geometry` /
  ## `effect_view` / `effect_contrast_view` records and their validators. It
  ## is the executor's output contract, not a presentation of one -- it calls
  ## nothing above layer 3, and every reader of a result (views, formatting,
  ## printing, plotting) sits above it.
  "compiler.R" = 4L, "execution-driver.R" = 4L, "kernel.R" = 4L,
  "task.R" = 4L, "storage.R" = 4L, "measurement-kernel.R" = 4L,
  "result.R" = 4L, "crossnobis-driver.R" = 4L,

  ## 5. results / views
  "views.R" = 5L, "geometry-entry.R" = 5L,
  "coupling-views.R" = 5L,
  "tomography.R" = 5L, "measurement-result.R" = 5L,
  "measurement-decomposition.R" = 5L, "format-results.R" = 5L,
  "print-methods.R" = 5L, "plot-methods.R" = 5L,

  ## 6. adapters and the public facade
  "adapter-bids.R" = 6L, "adapter-fmridesign.R" = 6L,
  "adapter-fmrireg.R" = 6L, "neuroim2-adapter.R" = 6L,
  "benchmark.R" = 6L, "example-data.R" = 6L, "evidence-api.R" = 6L
)

# Known upward edges, as "caller.R -> callee.R", each with the follow-up that
# would remove it. Sorted by layer of the caller.
# The register is empty. Every upward edge the 2026-08-15 audit found has been
# removed rather than tolerated, so any new entry here is a new violation and
# has to be argued for on its own terms.
allowed_upward <- character()

describe_edges <- function(rows) {
  if (!nrow(rows)) return(character())
  sort(paste(rows$from, "->", rows$to, ":", rows$symbol))
}

find_source_dir <- function() {
  candidates <- c(
    "../../R",                       # testthat::test_local / devtools::test
    "../../../R",
    "../../00_pkg_src/crossform/R"   # R CMD check on a built tarball
  )
  hit <- candidates[dir.exists(candidates)]
  if (!length(hit)) return(NULL)
  normalizePath(hit[[1L]])
}

.graph_cache <- new.env(parent = emptyenv())

internal_call_graph <- function(dir) {
  if (!is.null(.graph_cache$edges)) return(.graph_cache$edges)
  files <- list.files(dir, pattern = "[.]R$", full.names = TRUE)
  defs <- new.env(parent = emptyenv())
  parsed <- lapply(files, parse, keep.source = FALSE)
  names(parsed) <- basename(files)
  for (name in names(parsed)) {
    for (e in parsed[[name]]) {
      if (is.call(e) && (identical(e[[1L]], quote(`<-`)) ||
          identical(e[[1L]], quote(`=`))) && is.name(e[[2L]]) &&
          is.call(e[[3L]]) && identical(e[[3L]][[1L]], quote(`function`))) {
        assign(as.character(e[[2L]]), name, envir = defs)
      }
    }
  }
  collect <- function(x, out) {
    if (is.name(x)) {
      n <- as.character(x)
      if (nzchar(n) && exists(n, envir = defs, inherits = FALSE)) {
        assign(n, TRUE, envir = out)
      }
      return(invisible())
    }
    if (is.call(x) || is.pairlist(x)) {
      for (i in seq_along(x)) {
        el <- x[[i]]
        if (missing(el) || is.null(el)) next
        collect(el, out)
      }
    }
    invisible()
  }
  rows <- list()
  for (name in names(parsed)) {
    out <- new.env(parent = emptyenv())
    for (e in parsed[[name]]) collect(e, out)
    syms <- ls(out, all.names = TRUE)
    if (!length(syms)) next
    to <- vapply(syms, get, character(1), envir = defs)
    keep <- to != name
    if (!any(keep)) next
    rows[[length(rows) + 1L]] <- data.frame(from = name, to = unname(to[keep]),
      symbol = syms[keep], stringsAsFactors = FALSE)
  }
  .graph_cache$edges <- do.call(rbind, rows)
  .graph_cache$edges
}

test_that("every source file has a declared layer", {
  dir <- find_source_dir()
  skip_if(is.null(dir), "package sources are not available under this runner")
  files <- basename(list.files(dir, pattern = "[.]R$"))

  expect_setequal(files, names(layer_of))
})

test_that("no file calls upward through the layering, outside the register", {
  dir <- find_source_dir()
  skip_if(is.null(dir), "package sources are not available under this runner")
  edges <- internal_call_graph(dir)
  edges <- edges[edges$from %in% names(layer_of) &
    edges$to %in% names(layer_of), ]
  upward <- edges[layer_of[edges$to] > layer_of[edges$from], ]
  # `paste()` recycles its length-one separator against two zero-length
  # columns and returns `" -> "`, so the clean case has to be spelled out:
  # no upward edges means no observed edges, not one edge with empty ends.
  observed <- if (nrow(upward)) {
    sort(unique(paste(upward$from, "->", upward$to)))
  } else {
    character()
  }

  # No unregistered upward edge.
  expect_identical(setdiff(observed, allowed_upward), character())
  # No stale register entry: the debt list must shrink honestly.
  expect_identical(setdiff(allowed_upward, observed), character())
})

test_that("the compute core never calls into results, views, or printing", {
  # This is the edge this refactor removed: kernel.R used to reach up into
  # views.R (`.unsvec_symmetric()`, `.align_contrast()`) and result.R
  # (`.svec_symmetric()`). Those primitives now live in R/primitives.R.
  #
  # `result.R` stays on this list even though it is layer 4 now: the layering
  # would permit a sideways call, but a kernel that builds a result record is
  # a kernel that has stopped being a numerical primitive. Only the executor
  # constructs results.
  dir <- find_source_dir()
  skip_if(is.null(dir), "package sources are not available under this runner")
  edges <- internal_call_graph(dir)
  presentation <- c("result.R", "views.R", "coupling-views.R", "tomography.R",
    "measurement-result.R", "measurement-decomposition.R",
    "format-results.R", "print-methods.R", "plot-methods.R")
  offending <- edges[edges$from %in% c("kernel.R", "task.R") &
    edges$to %in% presentation, ]

  expect_identical(describe_edges(offending), character())
})

test_that("primitives depend on nothing above the primitive layer", {
  dir <- find_source_dir()
  skip_if(is.null(dir), "package sources are not available under this runner")
  edges <- internal_call_graph(dir)
  primitive_files <- names(layer_of)[layer_of == 1L]
  offending <- edges[edges$from %in% primitive_files &
    !edges$to %in% primitive_files, ]

  expect_identical(describe_edges(offending), character())
})

test_that("values never reach into the compiler or the executor", {
  # Registered exceptions are declaration constructors and capability records
  # rather than compilation; see `allowed_upward`.
  dir <- find_source_dir()
  skip_if(is.null(dir), "package sources are not available under this runner")
  edges <- internal_call_graph(dir)
  engine <- names(layer_of)[layer_of == 4L]
  values <- names(layer_of)[layer_of == 2L]
  offending <- edges[edges$from %in% values & edges$to %in% engine, ]
  pairs <- if (nrow(offending)) {
    sort(unique(paste(offending$from, "->", offending$to)))
  } else {
    character()
  }

  expect_identical(setdiff(pairs, allowed_upward), character())
})
