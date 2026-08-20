#!/usr/bin/env Rscript
#
# Build the crossform external review bundle.
#
#   Rscript exemplars/review-bundle/build-bundle.R [options]
#
# Assembles dist/crossform-review-bundle-<date>.tar.gz from the CURRENT working
# tree: contract documents, standalone oracles, exemplar evidence without the
# fetched data, the shipped certification artifacts and the helper that binds
# them to a source digest, and instructions for re-running the fast test suite.
#
# The bundle is evidence about one source state, so the manifest records the
# state it was cut from -- git commit, whether the tree was dirty, and the
# aggregate SHA-256 of R/*.R that certification artifacts bind against. A
# bundle built from a dirty tree is still built; it says so.
#
# This script only reads the tree. It writes nothing outside dist/.
#
# Options:
#   --output-dir=DIR    Where to write (default: <repo>/dist)
#   --date=YYYY-MM-DD   Bundle date stamp (default: today, UTC)
#   --max-file-mb=N     Per-file inclusion cap in MB (default: 4). Files over
#                       the cap are excluded from the payload and recorded in
#                       MANIFEST.txt with their sha256 and where to get them.
#   --keep-staging      Do not delete the staging directory after tarring.
#   --force             Overwrite an existing tarball of the same name.
#   --verify=PATH       Do not build. Extract PATH and re-check its MANIFEST.
#
# Exit status is non-zero on any failure, including a required file that is
# missing from the tree or a manifest hash that does not verify.

options(warn = 1L)

# ---------------------------------------------------------------- arguments --

parse_arguments <- function(arguments) {
  values <- list(
    output_dir = NULL,
    date = format(Sys.time(), "%Y-%m-%d", tz = "UTC"),
    max_file_mb = 4,
    keep_staging = FALSE,
    force = FALSE,
    verify = NULL
  )
  for (argument in arguments) {
    if (identical(argument, "--keep-staging")) {
      values$keep_staging <- TRUE
    } else if (identical(argument, "--force")) {
      values$force <- TRUE
    } else if (grepl("^--output-dir=", argument)) {
      values$output_dir <- sub("^--output-dir=", "", argument)
    } else if (grepl("^--date=", argument)) {
      values$date <- sub("^--date=", "", argument)
    } else if (grepl("^--max-file-mb=", argument)) {
      values$max_file_mb <- as.numeric(sub("^--max-file-mb=", "", argument))
    } else if (grepl("^--verify=", argument)) {
      values$verify <- sub("^--verify=", "", argument)
    } else if (argument %in% c("-h", "--help")) {
      cat(readLines(script_path())[2:34], sep = "\n")
      quit(status = 0L)
    } else {
      stop("unknown argument: ", argument, call. = FALSE)
    }
  }
  if (!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", values$date)) {
    stop("--date must be YYYY-MM-DD, got: ", values$date, call. = FALSE)
  }
  if (!is.finite(values$max_file_mb) || values$max_file_mb <= 0) {
    stop("--max-file-mb must be a positive number", call. = FALSE)
  }
  values
}

script_path <- function() {
  arguments <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", arguments, value = TRUE)
  if (length(hit)) {
    return(normalizePath(sub("^--file=", "", hit[[1L]]), mustWork = FALSE))
  }
  normalizePath("exemplars/review-bundle/build-bundle.R", mustWork = FALSE)
}

# The repository root is two levels above this script, whatever the caller's
# working directory is. Nothing here depends on being run from the root.
repository_root <- function() {
  root <- normalizePath(file.path(dirname(script_path()), "..", ".."),
    mustWork = FALSE)
  if (!file.exists(file.path(root, "DESCRIPTION"))) {
    stop("cannot locate the crossform repository root (looked in ", root, ")",
      call. = FALSE)
  }
  root
}

# ------------------------------------------------------------------ helpers --

say <- function(...) cat(..., "\n", sep = "")

file_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256")
}

# `git` may be absent, or the bundle may be cut from an exported tree with no
# history. Neither is fatal: the manifest records what it could learn.
git_field <- function(root, arguments) {
  value <- tryCatch(
    suppressWarnings(system2("git", c("-C", shQuote(root), arguments),
      stdout = TRUE, stderr = FALSE)),
    error = function(condition) character()
  )
  status <- attr(value, "status")
  if (!is.null(status) && !identical(as.integer(status), 0L)) {
    return(NA_character_)
  }
  if (!length(value)) {
    return("")
  }
  paste(value, collapse = "\n")
}

# ------------------------------------------------------------- the contents --

# Every entry is (source relative path, destination relative path, required).
# `required = TRUE` aborts the build when the file is absent, so a rename in
# the tree fails loudly instead of quietly shipping an incomplete bundle.
entry <- function(source, destination = source, required = TRUE) {
  list(source = source, destination = destination, required = required)
}

# A bounded, explicit expansion: one directory, one pattern, no recursion.
# Never a blind glob -- exemplar `data/` directories hold hundreds of MB and
# are reachable from none of these.
expand <- function(root, directory, pattern, destination_directory,
                   minimum = 1L) {
  absolute <- file.path(root, directory)
  files <- character()
  if (dir.exists(absolute)) {
    files <- sort(list.files(absolute, pattern = pattern, all.files = FALSE),
      method = "radix")
    files <- files[!grepl("^\\.", files)]
    files <- files[!dir.exists(file.path(absolute, files))]
  }
  if (length(files) < minimum) {
    stop(sprintf(
      "expected at least %d file(s) matching %s in %s, found %d",
      minimum, pattern, directory, length(files)
    ), call. = FALSE)
  }
  lapply(files, function(file) {
    entry(file.path(directory, file), file.path(destination_directory, file))
  })
}

bundle_contents <- function(root) {
  contents <- list()
  push <- function(...) contents <<- c(contents, list(...))
  push_all <- function(items) contents <<- c(contents, items)

  # (a) The contract documents. These are what a reviewer reads first; the
  # numbers everywhere else are claims against them.
  for (document in c(
    "conservative-geometry-contract.md",
    "population-form-contract.md",
    "crossform-execution-design.md",
    "architecture.md",
    "api-tiers.md",
    "certification-report.md"
  )) {
    push(entry(file.path("design", document)))
  }

  # (b) The oracles. Standalone-runnable; the three population-* ones do not
  # load crossform at all, which is the point of shipping them.
  push_all(expand(root, "design/oracles", "\\.R$", "design/oracles",
    minimum = 6L))

  # (c) Exemplar evidence, without the fetched data. `exemplars/*/data/` is
  # git-ignored in the source tree and is deliberately unreachable from here.
  for (file in c("README.md", "requirements.txt", "run-all.sh")) {
    push(entry(file.path("exemplars/rsatoolbox-parity", file)))
  }
  push_all(expand(root, "exemplars/rsatoolbox-parity", "\\.(R|py)$",
    "exemplars/rsatoolbox-parity", minimum = 5L))
  push_all(expand(root, "exemplars/rsatoolbox-parity/results", "\\.csv$",
    "exemplars/rsatoolbox-parity/results", minimum = 10L))

  push(entry("exemplars/er-rsa/README.md"))
  push_all(expand(root, "exemplars/er-rsa", "\\.R$", "exemplars/er-rsa",
    minimum = 4L))
  push_all(expand(root, "exemplars/er-rsa/results", "\\.csv$",
    "exemplars/er-rsa/results", minimum = 10L))

  # Haxby: the README, the scripts that produced the receipts, and the
  # committed smoke receipts themselves. Enumerated, not globbed over the
  # directory tree -- `exemplars/haxby2001/data/` is ~604 MB locally.
  push(entry("exemplars/haxby2001/README.md"))
  push_all(expand(root, "exemplars/haxby2001", "\\.R$", "exemplars/haxby2001",
    minimum = 7L))
  for (receipt in c(
    "smoke-report.md",
    "smoke-comparison.rds",
    "smoke-crossform.rds",
    "smoke-effectagram.rds",
    "smoke-rmvpa.rds",
    "smoke-uncertainty.rds",
    "coherent-configuration.rds",
    "coherent-configuration.png"
  )) {
    # Receipt names have moved once already (`effectagram` -> `crossform`), so
    # each is optional individually; the count is checked after resolution.
    push(entry(file.path("exemplars/haxby2001/results", receipt),
      required = FALSE))
  }

  for (file in c("DECISION.md", "README.md", "manifest.csv", "fetch.sh")) {
    push(entry(file.path("exemplars/population-slice2", file)))
  }

  # (d) The shipped certification artifacts: the gate records themselves plus
  # the per-gate summary CSVs.
  push_all(expand(root, "inst/extdata/certification", ".",
    "certification/artifacts", minimum = 15L))

  # (e) The binding helper, the digest function it uses, and how to re-certify.
  push(entry("tests/testthat/helper-certification.R",
    "certification/helper-certification.R"))
  push(entry("benchmarks/provenance.R", "certification/provenance.R"))
  push(entry("benchmarks/promote-artifacts.R",
    "certification/promote-artifacts.R"))
  push(entry("benchmarks/RECERTIFY.md", "certification/RECERTIFY.md"))
  push(entry("benchmarks/README.md", "certification/benchmarks-README.md"))
  push(entry("benchmarks/admission-coverage.R",
    "certification/admission-coverage.R"))

  # The bundle's own front matter, and the script that cut it. Shipping the
  # builder makes the bundle self-describing and lets a reviewer re-verify the
  # manifest with `Rscript build-bundle.R --verify=<tarball>`; the verify path
  # never looks for a repository root, so it works from anywhere.
  push(entry("exemplars/review-bundle/README.md", "README.md"))
  push(entry("exemplars/review-bundle/build-bundle.R", "build-bundle.R"))

  contents
}

# ---------------------------------------------------- provenance of the cut --

source_tree_digest <- function(root) {
  helper <- file.path(root, "benchmarks", "provenance.R")
  if (!file.exists(helper)) {
    return(NA_character_)
  }
  environment <- new.env(parent = globalenv())
  sys.source(helper, envir = environment)
  environment$.crossform_source_tree_digest(root)
}

description_version <- function(root) {
  fields <- tryCatch(read.dcf(file.path(root, "DESCRIPTION"),
    fields = "Version"), error = function(condition) NULL)
  if (is.null(fields) || is.na(fields[[1L, 1L]])) {
    return(NA_character_)
  }
  as.character(fields[[1L, 1L]])
}

short_digest <- function(digest) {
  if (!is.character(digest) || !length(digest) || is.na(digest[[1L]])) {
    return("<none>")
  }
  substr(sub("^sha256:", "", digest[[1L]]), 1L, 12L)
}

# Read the provenance block out of every shipped gate artifact and say, per
# artifact, whether it still binds to the current R/ digest. This is the one
# table a reviewer needs in order to know what the recorded gates are evidence
# for; see README.md, claim 5.
certification_binding <- function(root, current_digest) {
  directory <- file.path(root, "inst", "extdata", "certification")
  files <- sort(list.files(directory, pattern = "\\.rds$"), method = "radix")
  rows <- lapply(files, function(file) {
    artifact <- tryCatch(readRDS(file.path(directory, file)),
      error = function(condition) NULL)
    provenance <- if (is.list(artifact)) artifact$provenance else NULL
    recorded <- if (is.list(provenance)) provenance$source_digest else NULL
    recorded <- if (is.character(recorded) && length(recorded) == 1L &&
      !is.na(recorded) && nzchar(recorded)) recorded else NA_character_
    state <- if (is.na(recorded)) {
      "UNBOUND"
    } else if (identical(recorded, current_digest)) {
      "BOUND"
    } else {
      "STALE"
    }
    data.frame(
      artifact = file,
      recorded_digest = short_digest(recorded),
      state = state,
      package_version = if (is.list(provenance) &&
        !is.null(provenance$package_version)) {
        as.character(provenance$package_version)[[1L]]
      } else {
        NA_character_
      },
      recorded_at = if (is.list(provenance) &&
        !is.null(provenance$recorded_at)) {
        as.character(provenance$recorded_at)[[1L]]
      } else {
        NA_character_
      },
      runner = if (is.list(provenance) && !is.null(provenance$runner)) {
        as.character(provenance$runner)[[1L]]
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  })
  if (!length(rows)) {
    return(NULL)
  }
  do.call(rbind, rows)
}

# ------------------------------------------------------------ the fast-suite --

running_tests_page <- function(root, version) {
  c(
    "# Re-running the fast test suite from a source checkout",
    "",
    "One page. This is the layer a reviewer can run without fetching any",
    "data, without a Python environment, and without the ~12-minute",
    "re-certification sequence.",
    "",
    "## What you need",
    "",
    paste0("- R >= 4.3 (the recorded runs are R 4.5.1). crossform ", version,
      " is a source package; it has compiled code, so you need the usual"),
    "  toolchain (Xcode command line tools on macOS, `r-base-dev` on Debian).",
    "- The crossform source checkout. **This bundle is not the checkout.** It",
    "  is evidence cut from one; the manifest names the commit it was cut",
    "  from.",
    "- R packages: `devtools`, `testthat` (>= 3), `digest`, `Matrix`, `Rcpp`,",
    "  and the package's own Imports. `remotes::install_deps(dependencies =",
    "  TRUE)` from the checkout root installs them.",
    "",
    "## The run",
    "",
    "```sh",
    "cd /path/to/crossform",
    "Rscript -e 'devtools::test()'",
    "```",
    "",
    "That is the fast suite. It compiles `src/`, loads the source tree, and",
    "runs every file in `tests/testthat/`. On an idle Apple silicon laptop it",
    "takes a few minutes.",
    "",
    "Two tiers are held out of it deliberately:",
    "",
    "- **Scale gates** are opt-in. They are skipped unless",
    "  `CROSSFORM_RUN_SCALE_TESTS=true` is set, because they allocate at",
    "  brain scale and measure wall clock. Run them alone on an idle machine",
    "  or the timing gates report contention as regression.",
    "",
    "  ```sh",
    "  CROSSFORM_RUN_SCALE_TESTS=true Rscript -e 'devtools::test()'",
    "  ```",
    "",
    "- **Certification blocks** read the recorded gate artifacts in",
    "  `inst/extdata/certification/`. They do not re-run the gates; they read",
    "  a record, and they refuse to read it unless it still binds to the code",
    "  under test. See the next section.",
    "",
    "## Reading the skips, which is the point",
    "",
    "`certification/helper-certification.R` in this bundle is the mechanism.",
    "A recorded artifact is evidence for exactly one source tree, identified",
    "by the aggregate SHA-256 over the sorted per-file SHA-256 of `R/*.R`",
    "(`certification/provenance.R`). Three outcomes, all loud:",
    "",
    "| message | meaning |",
    "| --- | --- |",
    "| `CERTIFICATION STALE` | the artifact records a different `R/` digest than the checkout has. Any edit to any file under `R/`, including a comment, produces this. |",
    "| `CERTIFICATION UNBOUND` | the artifact carries no source digest at all, or only the weaker package-version tier was available (under `R CMD check` the sources are not shipped). |",
    "| `CERTIFICATION ABSENT` | no recorded artifact is available here. |",
    "",
    "A skip is not a pass. Grep for it:",
    "",
    "```sh",
    "Rscript -e 'devtools::test()' 2>&1 | grep -i 'CERTIFICATION'",
    "```",
    "",
    "`certification/BINDING.txt` in this bundle records, per artifact, which",
    "digest it binds to and whether that was the current one when the bundle",
    "was cut. Re-certifying is `certification/RECERTIFY.md`: about twelve",
    "minutes of runners on a frozen `R/`, then",
    "`certification/promote-artifacts.R`.",
    "",
    "## `R CMD check`",
    "",
    "```sh",
    "R CMD build /path/to/crossform && R CMD check --no-manual crossform_*.tar.gz",
    "```",
    "",
    "Under check the source tree is not shipped, so certification blocks fall",
    "back to the package-version tier and say so with a",
    "`CERTIFICATION TIER package_version` message. A version-level binding is",
    "not a source-level one; the helper prints the distinction rather than",
    "letting a pass be read as stronger than it is.",
    "",
    "## The other layers, and where they are documented",
    "",
    "| layer | how to run | where |",
    "| --- | --- | --- |",
    "| oracles | `Rscript design/oracles/<name>.R` | `README.md` claim 3 |",
    "| rsatoolbox parity | pinned Python env, then four scripts | `exemplars/rsatoolbox-parity/README.md` |",
    "| ER-RSA | three R scripts, no data fetch | `exemplars/er-rsa/README.md` |",
    "| Haxby | ~300 MB fetch, then five scripts | `exemplars/haxby2001/README.md` |",
    "| certification gates | the twelve-minute sequence | `certification/RECERTIFY.md` |",
    ""
  )
}

# ------------------------------------------------------------------- staging --

stage_bundle <- function(root, contents, staging, max_bytes) {
  included <- list()
  oversized <- list()
  missing <- character()

  for (item in contents) {
    source <- file.path(root, item$source)
    if (!file.exists(source) || dir.exists(source)) {
      if (isTRUE(item$required)) {
        stop("required file missing from the tree: ", item$source, call. = FALSE)
      }
      missing <- c(missing, item$source)
      next
    }
    size <- file.info(source)$size
    if (size > max_bytes) {
      # Do not silently drop it and do not silently blow up the tarball:
      # record what it is and how to get it.
      oversized[[length(oversized) + 1L]] <- data.frame(
        path = item$source, bytes = size, sha256 = file_sha256(source),
        stringsAsFactors = FALSE
      )
      next
    }
    destination <- file.path(staging, item$destination)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(source, destination, overwrite = TRUE, copy.date = TRUE)) {
      stop("failed to copy ", item$source, call. = FALSE)
    }
    included[[length(included) + 1L]] <- item$destination
  }

  list(
    included = unlist(included, use.names = FALSE),
    oversized = if (length(oversized)) do.call(rbind, oversized) else NULL,
    missing = missing
  )
}

write_manifest <- function(staging, provenance, oversized) {
  # Hash everything actually present in the staging directory, including the
  # files this script generated. MANIFEST.txt and SHA256SUMS are the only
  # exclusions; a file cannot carry its own hash.
  files <- sort(list.files(staging, recursive = TRUE, all.files = TRUE,
    no.. = TRUE), method = "radix")
  files <- setdiff(files, c("MANIFEST.txt", "SHA256SUMS"))
  files <- files[!dir.exists(file.path(staging, files))]

  information <- file.info(file.path(staging, files))
  hashes <- vapply(file.path(staging, files), file_sha256, character(1),
    USE.NAMES = FALSE)

  lines <- c(
    "# crossform external review bundle -- MANIFEST",
    "#",
    "# Every payload file below, with its SHA-256 and byte count. The bundle",
    "# is evidence about one source state; that state is recorded here.",
    "#",
    paste0("# bundle                 ", provenance$bundle_name),
    paste0("# built at (UTC)         ", provenance$built_at),
    paste0("# package                crossform ", provenance$version),
    paste0("# source git commit      ", provenance$git_commit),
    paste0("# source git branch      ", provenance$git_branch),
    paste0("# working tree at build  ", provenance$git_state),
    paste0("# R/ source digest       ", provenance$source_digest),
    paste0("# built with             ", provenance$r_version, " / ",
      provenance$platform),
    paste0("# payload files          ", length(files)),
    paste0("# payload bytes          ",
      format(sum(information$size), big.mark = ",", scientific = FALSE)),
    "#"
  )

  if (!identical(provenance$git_state, "clean")) {
    lines <- c(lines,
      "# NOTE: the working tree was not clean when this bundle was cut, so the",
      "#       commit hash above does not fully identify its contents. The",
      "#       modified paths were:",
      paste0("#         ", provenance$git_dirty_paths),
      "#")
  }

  if (!is.null(oversized)) {
    lines <- c(lines,
      paste0("# EXCLUDED for size (over ", provenance$max_file_mb,
        " MB per file). Hash and origin recorded;"),
      "# fetch these from the source checkout at the commit above.",
      "#",
      paste0("#   ", oversized$sha256, "  ",
        format(oversized$bytes, big.mark = ",", scientific = FALSE), "  ",
        oversized$path),
      "#")
  }

  lines <- c(lines,
    "# NOT included at any size: exemplars/*/data/. The Haxby subject-1",
    "# tarball (~300 MB) and the OpenNeuro ds003745 subset (7.21 GB) are",
    "# fetched, never versioned. Their checksums and fetch scripts are in the",
    "# exemplar READMEs and in exemplars/population-slice2/manifest.csv.",
    "#",
    "# sha256                                                            bytes  path",
    paste0(hashes, "  ",
      formatC(information$size, width = 9L, format = "d"), "  ", files))

  writeLines(lines, file.path(staging, "MANIFEST.txt"))

  # A plain `shasum -a 256 -c` input, so a reviewer without R can verify.
  writeLines(paste0(hashes, "  ", files), file.path(staging, "SHA256SUMS"))

  list(files = files, bytes = sum(information$size), hashes = hashes)
}

write_binding <- function(staging, binding, current_digest) {
  lines <- c(
    "# Certification artifacts: what source each one binds to",
    "#",
    "# tests/testthat/helper-certification.R (shipped here as",
    "# certification/helper-certification.R) refuses to read a recorded gate",
    "# unless the R/ digest it recorded still equals the current one. This is",
    "# that comparison, evaluated when the bundle was cut.",
    "#",
    paste0("# current R/ digest at build: ", current_digest),
    "#",
    "#   BOUND   -- recorded digest equals the current one; the test reads it.",
    "#   STALE   -- recorded against a different R/; the test SKIPS, loudly.",
    "#   UNBOUND -- no source digest recorded; the test SKIPS, loudly.",
    "#"
  )
  if (is.null(binding)) {
    lines <- c(lines, "# (no .rds artifacts found)")
  } else {
    width <- max(nchar(binding$artifact))
    lines <- c(lines,
      sprintf("%-*s  %-13s  %-8s  %-12s  %s", width, "artifact",
        "recorded R/", "state", "version", "recorded at"),
      sprintf("%-*s  %-13s  %-8s  %-12s  %s", width, binding$artifact,
        binding$recorded_digest, binding$state,
        ifelse(is.na(binding$package_version), "-", binding$package_version),
        ifelse(is.na(binding$recorded_at), "-", binding$recorded_at)))
    stale <- sum(binding$state != "BOUND")
    lines <- c(lines, "",
      sprintf("%d of %d artifacts do not bind to the R/ tree this bundle was cut from.",
        stale, nrow(binding)),
      if (stale > 0L) {
        paste0("Those gates are recorded results, not live ones. See README.md ",
          "claim 5 and certification/RECERTIFY.md.")
      })
  }
  writeLines(lines, file.path(staging, "certification", "BINDING.txt"))
  invisible(binding)
}

# -------------------------------------------------------------------- verify --

verify_bundle <- function(path) {
  if (!file.exists(path)) {
    stop("no such bundle: ", path, call. = FALSE)
  }
  scratch <- tempfile("crossform-review-verify-")
  dir.create(scratch, recursive = TRUE)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  status <- system2("tar", c("-xzf", shQuote(normalizePath(path)),
    "-C", shQuote(scratch)))
  if (!identical(as.integer(status), 0L)) {
    stop("extraction failed", call. = FALSE)
  }
  roots <- list.files(scratch)
  if (length(roots) != 1L) {
    stop("expected exactly one top-level directory, found ", length(roots),
      call. = FALSE)
  }
  extracted <- file.path(scratch, roots[[1L]])
  sums <- file.path(extracted, "SHA256SUMS")
  if (!file.exists(sums)) {
    stop("bundle has no SHA256SUMS", call. = FALSE)
  }
  records <- readLines(sums)
  records <- records[nzchar(records)]
  expected <- sub("^([0-9a-f]{64})  (.*)$", "\\1", records)
  paths <- sub("^([0-9a-f]{64})  (.*)$", "\\2", records)

  present <- sort(list.files(extracted, recursive = TRUE, all.files = TRUE,
    no.. = TRUE), method = "radix")
  present <- setdiff(present, c("MANIFEST.txt", "SHA256SUMS"))

  bad <- character()
  absent <- character()
  for (index in seq_along(paths)) {
    target <- file.path(extracted, paths[[index]])
    if (!file.exists(target)) {
      absent <- c(absent, paths[[index]])
      next
    }
    if (!identical(file_sha256(target), expected[[index]])) {
      bad <- c(bad, paths[[index]])
    }
  }
  extra <- setdiff(present, paths)

  say("verifying ", basename(path))
  say("  manifest records : ", length(paths), " files")
  say("  hashes matched   : ", length(paths) - length(bad) - length(absent))
  say("  hash mismatches  : ", length(bad))
  say("  missing from tar : ", length(absent))
  say("  unlisted in tar  : ", length(extra))
  if (length(bad)) say("  MISMATCH: ", paste(bad, collapse = ", "))
  if (length(absent)) say("  MISSING: ", paste(absent, collapse = ", "))
  if (length(extra)) say("  UNLISTED: ", paste(extra, collapse = ", "))
  ok <- !length(bad) && !length(absent) && !length(extra)
  say(if (ok) "  RESULT: OK" else "  RESULT: FAILED")
  if (!ok) quit(status = 1L)
  invisible(TRUE)
}

# ---------------------------------------------------------------------- main --

main <- function() {
  arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("the 'digest' package is required (crossform imports it)",
      call. = FALSE)
  }
  if (!is.null(arguments$verify)) {
    verify_bundle(arguments$verify)
    return(invisible(TRUE))
  }

  root <- repository_root()
  output_dir <- if (is.null(arguments$output_dir)) {
    file.path(root, "dist")
  } else {
    arguments$output_dir
  }
  name <- paste0("crossform-review-bundle-", arguments$date)
  tarball <- file.path(output_dir, paste0(name, ".tar.gz"))
  staging <- file.path(output_dir, name)

  if (file.exists(tarball) && !arguments$force) {
    stop("refusing to overwrite ", tarball, " (pass --force)", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  # dist/ holds only build output. Ignore all of it, including this file, so
  # a staging directory never lands in `git status` and never gets committed.
  writeLines(c(
    "# Build output only: the review bundle is cut from the tree on demand.",
    "# Nothing under dist/ is versioned, including this file.",
    "*"
  ), file.path(output_dir, ".gitignore"))

  unlink(staging, recursive = TRUE)
  dir.create(staging, recursive = TRUE)

  say("crossform review bundle")
  say("  root    : ", root)
  say("  output  : ", tarball)

  current_digest <- source_tree_digest(root)
  version <- description_version(root)
  dirty <- git_field(root, c("status", "--porcelain"))
  dirty_paths <- if (is.na(dirty) || !nzchar(dirty)) {
    character()
  } else {
    strsplit(dirty, "\n", fixed = TRUE)[[1L]]
  }
  provenance <- list(
    bundle_name = name,
    built_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC",
      usetz = TRUE),
    version = version,
    git_commit = git_field(root, c("rev-parse", "HEAD")),
    git_branch = git_field(root, c("rev-parse", "--abbrev-ref", "HEAD")),
    git_state = if (is.na(dirty)) {
      "unknown (no git)"
    } else if (!length(dirty_paths)) {
      "clean"
    } else {
      paste0("DIRTY (", length(dirty_paths), " path(s))")
    },
    git_dirty_paths = dirty_paths,
    source_digest = current_digest,
    r_version = R.version.string,
    platform = R.version$platform,
    max_file_mb = arguments$max_file_mb
  )
  say("  commit  : ", substr(provenance$git_commit, 1L, 12L), " (",
    provenance$git_state, ")")
  say("  R/ digest: ", short_digest(current_digest))

  contents <- bundle_contents(root)
  staged <- stage_bundle(root, contents, staging,
    arguments$max_file_mb * 1024^2)

  # The Haxby receipts are individually optional because their names have
  # moved; the set is not. Fail if the arm brought nothing.
  receipts <- grep("^exemplars/haxby2001/results/", staged$included,
    value = TRUE)
  if (length(receipts) < 4L) {
    stop("only ", length(receipts), " Haxby receipt(s) resolved; expected the ",
      "committed smoke set", call. = FALSE)
  }
  if (length(staged$missing)) {
    say("  optional files absent from the tree (", length(staged$missing),
      "): ", paste(basename(staged$missing), collapse = ", "))
  }

  writeLines(running_tests_page(root, version),
    file.path(staging, "RUNNING-TESTS.md"))
  write_binding(staging, certification_binding(root, current_digest),
    current_digest)

  manifest <- write_manifest(staging, provenance, staged$oversized)

  # COPYFILE_DISABLE keeps macOS from writing `._` AppleDouble members that a
  # reviewer on Linux would see as junk files the manifest does not list.
  status <- system2("tar",
    c("-czf", shQuote(tarball), "-C", shQuote(output_dir), shQuote(name)),
    env = "COPYFILE_DISABLE=1")
  if (!identical(as.integer(status), 0L)) {
    stop("tar failed with status ", status, call. = FALSE)
  }

  tarball_bytes <- file.info(tarball)$size
  say("")
  say("  files    : ", length(manifest$files))
  say("  payload  : ", format(round(manifest$bytes / 1024^2, 2), nsmall = 2),
    " MB uncompressed")
  say("  tarball  : ", format(round(tarball_bytes / 1024^2, 2), nsmall = 2),
    " MB")
  if (!is.null(staged$oversized)) {
    say("  excluded for size: ", nrow(staged$oversized), " file(s) (see ",
      "MANIFEST.txt)")
  }
  if (tarball_bytes > 25 * 1024^2) {
    say("  WARNING: tarball exceeds the 25 MB budget.")
  }

  if (!arguments$keep_staging) {
    unlink(staging, recursive = TRUE)
  } else {
    say("  staging  : ", staging, " (kept)")
  }

  say("")
  say("  verify with:")
  say("    Rscript ", file.path("exemplars", "review-bundle", "build-bundle.R"),
    " --verify=", tarball)
  invisible(TRUE)
}

main()
