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
  ## Layer 2 has an internal order of its own, recorded in
  ## design/architecture.md under "The value order": spaces -> maps and
  ## extractors -> sources -> relations -> fits -> queries, pairings and
  ## frames. It is not enforced numerically here, because a value file may
  ## call sideways; what the order buys is that layer 2 now contains no
  ## cycle at all, which the receipt test below relies on.
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
  "relation-session.R" = 2L,
  "residual-statistics.R" = 2L, "bridge.R" = 2L,
  "transport.R" = 2L, "crossform-package.R" = 2L,

  ## 3. plans — scientific estimands and their identities
  "geometry-plan.R" = 3L, "relation-plan.R" = 3L, "population-plan.R" = 3L,
  "crossnobis.R" = 3L,
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

# The largest strongly connected component the file call graph may contain.
#
# The value is 1: every component is a single file, i.e. `R/` is acyclic and
# every file can be read, changed, and reasoned about without holding another
# one open. Layer 2 was made acyclic in commit 9c8be27 ("Record the layer-2
# value order and untangle the nine-file SCC"), which was the last cycle in
# the package; the whole graph has been acyclic since.
#
# This is a ratchet, not a threshold. Lowering it is impossible (1 is the
# floor -- a component is at least the file itself). Raising it means
# accepting a cycle, and the only honest way to do that is in a commit that
# says which files are in it and what would remove it, the way
# `allowed_upward` above carries its debt. Prefer the alternative: cycles in
# this package have every time been one edge that pointed the wrong way, and
# `Rscript benchmarks/call-graph-scc.R` prints the component so you can find
# it. The layering test above catches upward edges *across* layers; this one
# catches the sideways loop inside a layer that the layering permits.
scc_ratchet <- 1L

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
    ## `x$field` and `x@slot` are not calls to `field` or `slot`: the
    ## right-hand symbol names a component of `x`, so descend into the
    ## left-hand side only. Without this rule a record field that shares a
    ## name with a function defined in another file -- and in a package whose
    ## values are named after the things they hold, `x$relation`,
    ## `plan$pairing`, `x$task$left_relation$effect_space` all are -- records
    ## a call the file never makes. There were 44 such phantom edges, and
    ## they were not harmless: they inflated the file graph by 21 file-pairs
    ## and made files look mutually dependent when only their vocabulary
    ## overlapped.
    if (is.call(x) && length(x) == 3L && is.name(x[[3L]]) &&
        (identical(x[[1L]], quote(`$`)) || identical(x[[1L]], quote(`@`)))) {
      collect(x[[2L]], out)
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

# Tarjan's strongly connected components, the same implementation
# benchmarks/call-graph-scc.R runs, so a failure here and the benchmark's
# output name the same component. Recursion depth is bounded by the file
# count, which is under a hundred.
tarjan_scc <- function(nodes, edges) {
  adj <- split(edges$to, factor(edges$from, levels = nodes))
  index <- setNames(rep(NA_integer_, length(nodes)), nodes)
  low <- index
  on_stack <- setNames(rep(FALSE, length(nodes)), nodes)
  stack <- character()
  counter <- 0L
  comps <- list()
  strongconnect <- function(v) {
    counter <<- counter + 1L
    index[[v]] <<- counter
    low[[v]] <<- counter
    stack <<- c(stack, v)
    on_stack[[v]] <<- TRUE
    for (w in unique(adj[[v]])) {
      if (is.na(index[[w]])) {
        strongconnect(w)
        low[[v]] <<- min(low[[v]], low[[w]])
      } else if (on_stack[[w]]) {
        low[[v]] <<- min(low[[v]], index[[w]])
      }
    }
    if (low[[v]] == index[[v]]) {
      pos <- match(v, stack)
      comp <- stack[pos:length(stack)]
      stack <<- stack[seq_len(pos - 1L)]
      on_stack[comp] <<- FALSE
      comps[[length(comps) + 1L]] <<- sort(comp)
    }
    invisible()
  }
  for (v in nodes) if (is.na(index[[v]])) strongconnect(v)
  comps
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

test_that("the file call graph has no component larger than the ratchet", {
  # The layering test is a statement about direction; this one is a statement
  # about shape, and it is the stronger of the two. A file may call sideways
  # within its layer without violating the layering, and two files that call
  # each other sideways are one unit of code wearing two filenames: neither
  # can be read alone, neither can be changed alone, and the pair will grow
  # until someone breaks it. The specific cycles this package has had --
  # layer 2's nine-file tangle, the receipt's one edge back into it, the
  # evidence-sampling triple -- each have a test of their own above, because
  # each has a story worth telling. This test is the general case: it does
  # not know which cycle to expect, only that there is to be none.
  dir <- find_source_dir()
  skip_if(is.null(dir), "package sources are not available under this runner")
  edges <- internal_call_graph(dir)
  # Isolated files are nodes too -- a file nothing calls and that calls
  # nothing must still appear, or the graph the components are computed over
  # is not the package.
  nodes <- sort(unique(c(edges$from, edges$to,
    basename(list.files(dir, pattern = "[.]R$")))))
  comps <- tarjan_scc(nodes, edges)

  # Name the members on failure. A bare size tells you a cycle exists;
  # the members tell you where to look, and `SCC_FOCUS=<file> Rscript
  # benchmarks/call-graph-scc.R` then prints the edges that hold it shut.
  offending <- comps[lengths(comps) > scc_ratchet]
  described <- if (length(offending)) {
    sort(vapply(offending, function(cm) {
      paste0("SCC (", length(cm), "): ", paste(cm, collapse = " "))
    }, character(1)))
  } else {
    character()
  }

  expect_identical(described, character())
  expect_lte(max(lengths(comps)), scc_ratchet)
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

test_that("a receipt records canonical values and depends on none of them", {
  # A receipt is a record of what an execution used, so nothing it names may
  # be something it also helps define. Until 2026-08-17 `receipt.R` had one
  # outgoing edge into the layer-2 value tangle -- it called `memory_plan()`
  # to rebuild a plan from its own categories and check that the derived
  # totals followed -- and that single edge held five other files inside the
  # cycle with it: the tangle was fifteen files, and dropped to nine when the
  # rebuild moved to `.workspace_plan_record()` in the primitive layer.
  #
  # The list below is that fifteen-file component as design/architecture.md
  # recorded it, minus `receipt.R`. Files may still call *into* the receipt --
  # `relation.R` and `study-facts.R` did, until the value they were calling
  # for moved out of it -- and that is the direction the design wants. What
  # must stay empty is the other one.
  #
  # `capabilities.R` is not on the list, and its absence is the point of the
  # second assertion. `source_capabilities()` and its validator were defined
  # in `receipt.R`, because a receipt is where a capability is finally
  # written down; but `relation.R`, `relation-fit.R`, and `study-facts.R` all
  # have to *state* a capability long before any receipt exists, so the
  # record moved to `capabilities.R`, the file named for it, and the receipt
  # now calls down to the value it records. That is a legitimate edge and a
  # list check cannot tell it from an illegitimate one. What distinguishes
  # them is whether it closes a loop, so that is what is asserted.
  dir <- find_source_dir()
  skip_if(is.null(dir), "package sources are not available under this runner")
  edges <- internal_call_graph(dir)
  value_cycle <- c("effect-map.R", "effect-space.R",
    "extractor.R", "measurement.R", "memory-plan.R", "metric.R",
    "operations.R", "pairing.R", "relation-fit.R", "relation.R", "scope.R",
    "source.R", "study-facts.R")
  offending <- edges[edges$from == "receipt.R" & edges$to %in% value_cycle, ]

  expect_identical(describe_edges(offending), character())

  # Nothing the receipt calls may reach the receipt again, by any route.
  reachable <- character()
  frontier <- unique(edges$to[edges$from == "receipt.R"])
  while (length(frontier)) {
    reachable <- union(reachable, frontier)
    frontier <- setdiff(unique(edges$to[edges$from %in% frontier]), reachable)
  }

  expect_false("receipt.R" %in% reachable)
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

test_that("the evidence-sampling triple is a DAG", {
  # The layering test above only orders files across layers; these three sit
  # together in layer 3 and so are free to call sideways. Their intended
  # order is the data flow -- plan (evidence-sampling.R) -> covariance kernel
  # -> relation-fit product -- and nothing enforces it but this test. A single
  # sideways call back down the flow turns the triple into one component and
  # makes the three files unreadable in isolation, which is what ticket C3
  # removed; see design/architecture.md, "The evidence-sampling triple".
  dir <- find_source_dir()
  skip_if(is.null(dir), "package sources are not available under this runner")
  triple <- c("evidence-sampling.R", "evidence-sampling-kernel.R",
    "evidence-sampling-product.R")
  edges <- internal_call_graph(dir)
  induced <- edges[edges$from %in% triple & edges$to %in% triple, ]

  # Transitive closure of the induced subgraph (Floyd-Warshall over three
  # nodes). A file that reaches itself sits on a cycle, whether that cycle is
  # a mutual pair or the full three-step loop.
  reach <- matrix(FALSE, length(triple), length(triple),
    dimnames = list(triple, triple))
  if (nrow(induced)) {
    reach[cbind(induced$from, induced$to)] <- TRUE
  }
  for (k in triple) {
    for (i in triple) {
      if (reach[i, k]) reach[i, ] <- reach[i, ] | reach[k, ]
    }
  }
  on_cycle <- triple[as.logical(diag(reach))]

  # Report the offending edges, not just the file names, so a failure names
  # the call to move rather than sending the reader back to the call graph.
  expect_identical(
    if (length(on_cycle)) describe_edges(induced) else character(),
    character()
  )
  expect_identical(on_cycle, character())
})
