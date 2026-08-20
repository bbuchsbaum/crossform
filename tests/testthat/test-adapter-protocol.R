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
# **The register is empty.** It held twelve entries in four groups, and the
# follow-up ticket to this file closed all twelve, so the four adapters are now
# written entirely within the protocol they document. What that cost, group by
# group, is worth recording here because it is the evidence for the claim:
#
#   fact re-validation      -- deleted. Each adapter had been re-checking a
#     record it was about to hand to a constructor that checks it on intake
#     (`study()`), or re-checking its own character vector. `bids_study()` now
#     gates on the class and lets `study()` raise the refusal it always
#     raised; `fmridesign_design_model()` asks `study_capabilities()`, the
#     public verb that runs the same validation.
#   adapter version pinning -- exported, as `adapter_version_certificate()`.
#     `vignette("crossform-extending")` *obliges* an external adapter to
#     certify against installed versions; a protocol that requires a behaviour
#     and hides its only implementation is not a protocol.
#   external plan execution -- rewritten. `relation_plan_receipts()` validates
#     the plan and yields the receipts, and the planned sources are built in
#     the adapter from the plan's documented fields plus
#     `source_capabilities()`, which is what an external executor must do.
#   spatial provider internals -- moved into the public constructors.
#     `volume_domain(metadata = )` records what a provider knows and the array
#     does not; `additive_frame(members = )` compiles a provider's own
#     neighborhoods the way `compile_frame()` compiles crossform's.
#
# Keep the register. A new gap belongs in it, with its justification, rather
# than being argued away -- and the ratchet means it can only shrink again.
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
# from outside the package today. Sorted by file, then symbol. Empty since the
# twelve entries above were closed; the two tests below still hold it to both
# halves of the ratchet, so an entry added without a closure fails as loudly as
# one removed without a fix.
unsanctioned_internal_calls <- character()

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
  # `study_capabilities()` joined when the adapter stopped calling the study
  # validator directly: it is the public verb that runs that validation, and
  # asking a study to prove itself through it is what an external adapter does.
  "adapter-fmridesign.R" = c(
    "adapter_version_certificate", "coefficient_parameterization",
    "condition_space", "design_model", "study_capabilities"
  ),
  # The two that closed the external-executor gap. `relation_plan_receipts()`
  # is how a plan is validated from outside, and `source_capabilities()` is how
  # the derived, row-restricted sources an executor hands to `relation()`
  # declare what they can honestly do.
  "adapter-fmrireg.R" = c(
    "adapter_version_certificate", "effect_extractor", "relation",
    "relation_fit", "relation_plan_receipts", "source_capabilities"
  ),
  # `searchlights()` joined this list with ticket D3: a multiscale
  # `neuroim2_searchlights()` request delegates its refusals to the
  # constructor that owns the multiscale rules, so both spatial providers
  # refuse the same things in the same words. That route also calls
  # `frame_family()`, which the two purity checks above do see and this
  # heads-mode count does not, because `do.call()` makes it an argument rather
  # than a call head. Both are exported, so the protocol holds either way;
  # only the published count reads the narrower of the two.
  #
  # `abstract_domain()` left when the adapter stopped assembling a volume
  # domain out of an abstract one and asked `volume_domain()` for the volume
  # domain it actually wanted.
  "neuroim2-adapter.R" = c(
    "additive_frame", "neuroim2_volume_domain", "searchlights",
    "volume_domain"
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
  # them (twelve after decision 1, thirteen after ticket D3, seventeen once the
  # debt register was closed), of which [five] are developer-tier". That number
  # is a claim about this tree, so it is checked against this tree. It went up
  # by four when the register closed, which is the shape of the whole exercise:
  # an adapter that reaches for an internal makes no exported call, and an
  # adapter written within the protocol makes one.
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
  expect_identical(length(unique(unlist(observed, use.names = FALSE))), 17L)

  # Of those, three are developer-tier: the design-in and error-channel-in
  # seams, plus `source_capabilities()`, which `fmrireg_relation()` began
  # calling when it started declaring its planned sources itself instead of
  # borrowing the plan's internal builder. The remaining two entry points are
  # unexercised because no in-tree adapter brings its own out-of-memory
  # source -- see the ledger's "Verdict on the developer set".
  developer_used <- intersect(
    unique(unlist(observed, use.names = FALSE)),
    sanctioned_developer_entry_points
  )

  expect_setequal(developer_used,
    c("effect_extractor", "relation_fit", "source_capabilities"))
})
