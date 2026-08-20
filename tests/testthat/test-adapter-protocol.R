# The developer protocol, enforced ------------------------------------------
#
# `design/api-tiers.md` ("Tier: developer") claims that an extension package
# needs five sanctioned entry points -- `file_matrix_source()`,
# `source_capabilities()`, `effect_extractor()`, `relation_fit()`,
# `relation_block()` -- plus the ordinary exported core, and nothing else.
# `vignette("crossform-extending")` publishes that claim as a promise.
#
# The four in-tree adapters are the only evidence the package has for it. This
# file holds them to it: every crossform function an adapter calls must be
#
#   (a) a sanctioned developer entry point, or
#   (b) an exported function (any tier of the ledger -- every ledger row is an
#       export, so `getNamespaceExports()` is the ledger's own universe), or
#   (c) a dot-internal defined in layer 1 (`primitives.R`, `check.R`,
#       `conditions.R`, `message-helpers.R`) -- the shared leaf layer named in
#       design/architecture.md, which an adapter may use for argument checks
#       and refusal raising because it is what every layer uses.
#
# Anything else is an adapter reaching into machinery an external package
# cannot reach, which would make the published protocol false.
#
# `unsanctioned_internal_calls` below is the debt register for the calls that
# are outside the protocol today. It works like `allowed_upward` in
# test-architecture.R: the test fails on an unregistered violation *and* on a
# registered one that no longer exists, so the list is green but ratcheting and
# cannot rot. Every entry is a real gap in the published protocol -- an
# external package cannot write these adapters without `:::`.
#
# The call graph is parsed here rather than imported from test-architecture.R
# because testthat gives each test file its own environment; only helpers are
# shared, and this file should not add one for two small functions.

sanctioned_developer_entry_points <- c(
  "file_matrix_source", "source_capabilities", "effect_extractor",
  "relation_fit", "relation_block"
)

# design/architecture.md layer 1: pure leaves with no crossform dependency
# outside the layer. `RcppExports.R` is Rcpp-generated `.Call()` glue.
layer_one_files <- c(
  "primitives.R", "message-helpers.R", "conditions.R", "check.R",
  "RcppExports.R"
)

adapter_files <- c(
  "adapter-bids.R", "adapter-fmridesign.R", "adapter-fmrireg.R",
  "neuroim2-adapter.R"
)

# "<adapter file> -> <internal symbol>", with the reason it is not reachable
# from outside the package today. Sorted by file, then symbol.
unsanctioned_internal_calls <- c(
  # Fact validators. An external adapter that builds an `observations()` or a
  # `study()` re-validates its own input before handing it on; the protocol
  # offers no public re-validation, so the check is skipped or duplicated.
  "adapter-bids.R -> .validate_observations",
  "adapter-bids.R -> .validate_partition_names",
  "adapter-fmridesign.R -> .validate_study",
  # Adapter version pinning. `.require_adapter_version()` is what turns "the
  # installed fmridesign is not the certified one" into an
  # `installed_compiler_adapter` / `supported_compiler_version` refusal. An
  # external adapter must hand-roll the equivalent refusal.
  "adapter-fmridesign.R -> .require_adapter_version",
  "adapter-fmrireg.R -> .require_adapter_version",
  # Executing a `plan_relation()` from outside the package. These two are the
  # whole of the "external executor" seam and neither is public: an external
  # engine cannot read a plan's planned sources without `:::`.
  "adapter-fmrireg.R -> .planned_observation_sources",
  "adapter-fmrireg.R -> .validate_relation_plan",
  # Domain and frame construction below the public constructors. An external
  # spatial provider can reach `abstract_domain()` and `additive_frame()`, but
  # not the support-index machinery that makes a searchlight frame carry its
  # membership pattern, nor the domain-identity comparison that makes writing
  # results back to voxels safe.
  "neuroim2-adapter.R -> .new_domain",
  "neuroim2-adapter.R -> .normalize_frame",
  "neuroim2-adapter.R -> .same_domain_reference",
  "neuroim2-adapter.R -> .support_index_from_members",
  "neuroim2-adapter.R -> .support_index_membership",
  "neuroim2-adapter.R -> .validate_domain"
)

# The exported crossform functions each adapter calls, as `design/api-tiers.md`
# records them in its "Verdict on the developer set". Kept here so the ledger's
# published count stays true of the tree.
adapter_exported_calls <- list(
  # `bids_confounds()` and `bids_events()` left this list when the subtraction
  # release demoted them (ledger decision 1); `adapter-bids.R` still calls
  # both, but they are now its own internals rather than exports.
  "adapter-bids.R" = c(
    "observation_confounds", "observation_events", "study"
  ),
  "adapter-fmridesign.R" = c(
    "coefficient_parameterization", "condition_space", "design_model"
  ),
  "adapter-fmrireg.R" = c("effect_extractor", "relation", "relation_fit"),
  "neuroim2-adapter.R" = c(
    "abstract_domain", "additive_frame", "neuroim2_volume_domain"
  )
)

protocol_source_dir <- function() {
  candidates <- c(
    "../../R",                       # testthat::test_local / devtools::test
    "../../../R",
    "../../00_pkg_src/crossform/R"   # R CMD check on a built tarball
  )
  hit <- candidates[dir.exists(candidates)]
  if (!length(hit)) return(NULL)
  normalizePath(hit[[1L]])
}

# Every top-level `name <- function(...)` in R/, mapped to its defining file.
protocol_definitions <- function(parsed) {
  defs <- new.env(parent = emptyenv())
  for (file in names(parsed)) {
    for (expression in parsed[[file]]) {
      if (is.call(expression) &&
          (identical(expression[[1L]], quote(`<-`)) ||
           identical(expression[[1L]], quote(`=`))) &&
          is.name(expression[[2L]]) && is.call(expression[[3L]]) &&
          identical(expression[[3L]][[1L]], quote(`function`))) {
        assign(as.character(expression[[2L]]), file, envir = defs)
      }
    }
  }
  defs
}

# Two collectors, because the two questions want different precision.
#
# `heads` takes only what is actually called -- `f(x)`, `crossform:::f(x)` --
# which is what "the functions an adapter calls" means when the answer is
# published as a count.
#
# `names_used` takes every symbol, so a function passed as a value
# (`vapply(x, .validate_source_capabilities, ...)`) still counts as a use. It
# over-collects: an argument named `study` is indistinguishable from a call to
# `study()`. That is the safe direction for the purity check, and harmless,
# since a local variable never begins with a dot.
protocol_symbols <- function(expression, mode) {
  found <- new.env(parent = emptyenv())
  walk <- function(x) {
    if (is.name(x)) {
      if (identical(mode, "names")) assign(as.character(x), TRUE, envir = found)
      return(invisible())
    }
    if (is.call(x)) {
      if (identical(mode, "heads")) {
        head <- x[[1L]]
        if (is.call(head) && is.name(head[[1L]]) &&
            as.character(head[[1L]]) %in% c("::", ":::")) {
          head <- head[[3L]]
        }
        if (is.name(head)) assign(as.character(head), TRUE, envir = found)
      }
    }
    if (is.call(x) || is.pairlist(x)) {
      for (index in seq_along(x)) {
        element <- x[[index]]
        if (missing(element) || is.null(element)) next
        walk(element)
      }
    }
    invisible()
  }
  walk(expression)
  sort(ls(found, all.names = TRUE))
}

# Crossform symbols an adapter file uses, with the file that defines each.
protocol_uses <- function(parsed, defs, file, mode) {
  used <- character()
  for (expression in parsed[[file]]) {
    used <- union(used, protocol_symbols(expression, mode))
  }
  used <- sort(used[vapply(used, exists, logical(1),
    envir = defs, inherits = FALSE)])
  stats::setNames(vapply(used, get, character(1), envir = defs), used)
}

protocol_parsed_sources <- function(dir) {
  files <- list.files(dir, pattern = "[.]R$", full.names = TRUE)
  parsed <- lapply(files, parse, keep.source = FALSE)
  names(parsed) <- basename(files)
  parsed
}

test_that("the five sanctioned developer entry points are exported", {
  exported <- getNamespaceExports("crossform")

  expect_true(all(sanctioned_developer_entry_points %in% exported))
  # `as_neurovol()` is the sixth extension point and a different kind: an S3
  # generic a third party registers a method on (commit f38dcef), so what has
  # to hold is that it still dispatches.
  expect_true("as_neurovol" %in% exported)
  expect_true(any(grepl("UseMethod", deparse(body(crossform::as_neurovol)),
    fixed = TRUE)))
})

test_that("adapters call no crossform internal outside the protocol", {
  dir <- protocol_source_dir()
  skip_if(is.null(dir), "package sources are not available under this runner")
  parsed <- protocol_parsed_sources(dir)
  expect_true(all(adapter_files %in% names(parsed)))
  defs <- protocol_definitions(parsed)
  exported <- getNamespaceExports("crossform")

  observed <- character()
  for (file in adapter_files) {
    uses <- protocol_uses(parsed, defs, file, mode = "names")
    # Same-file helpers are the adapter's own; layer 1 and the exported
    # surface are the protocol. Everything else is a reach.
    outside <- uses[uses != file & !uses %in% layer_one_files &
      !names(uses) %in% exported]
    if (length(outside)) {
      observed <- c(observed, paste(file, "->", names(outside)))
    }
  }
  observed <- sort(observed)

  # No adapter reaches into an internal that is not already registered.
  expect_identical(setdiff(observed, unsanctioned_internal_calls), character())
  # No stale register entry: the debt list must shrink honestly.
  expect_identical(setdiff(unsanctioned_internal_calls, observed), character())
})

test_that("adapters use only layer-1 internals plus the exported surface", {
  # The same rule stated positively, so a reader of the failure sees what the
  # protocol allows rather than only what it forbids.
  dir <- protocol_source_dir()
  skip_if(is.null(dir), "package sources are not available under this runner")
  parsed <- protocol_parsed_sources(dir)
  defs <- protocol_definitions(parsed)
  exported <- getNamespaceExports("crossform")
  registered <- sub("^.* -> ", "", unsanctioned_internal_calls)

  for (file in adapter_files) {
    uses <- protocol_uses(parsed, defs, file, mode = "names")
    internal <- uses[uses != file & !names(uses) %in% exported]
    unexplained <- names(internal)[
      !uses[names(internal)] %in% layer_one_files &
        !names(internal) %in% registered
    ]

    expect_identical(unexplained, character())
  }
})

test_that("the ledger's record of what an adapter calls is still true", {
  # design/api-tiers.md publishes "about fourteen crossform functions between
  # them (twelve after decision 1), of which [five] are developer-tier". That
  # number is a claim about this tree, so it is checked against this tree.
  dir <- protocol_source_dir()
  skip_if(is.null(dir), "package sources are not available under this runner")
  parsed <- protocol_parsed_sources(dir)
  defs <- protocol_definitions(parsed)
  exported <- getNamespaceExports("crossform")

  observed <- lapply(adapter_files, function(file) {
    uses <- protocol_uses(parsed, defs, file, mode = "heads")
    sort(intersect(names(uses), exported))
  })
  names(observed) <- adapter_files

  expect_identical(observed, adapter_exported_calls)
  expect_identical(length(unique(unlist(observed, use.names = FALSE))), 12L)

  # Of those, exactly two are developer-tier: the design-in and
  # error-channel-in seams. The other three sanctioned entry points are
  # unexercised because no in-tree adapter brings its own out-of-memory
  # source -- see the ledger's "Verdict on the developer set".
  developer_used <- intersect(
    unique(unlist(observed, use.names = FALSE)),
    sanctioned_developer_entry_points
  )

  expect_setequal(developer_used, c("effect_extractor", "relation_fit"))
})
