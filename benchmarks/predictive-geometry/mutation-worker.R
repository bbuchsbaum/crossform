args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: mutation-worker.R <repo> <mutation-id> <output-rds>")
repo <- normalizePath(args[1], mustWork = TRUE); id <- args[2]; output <- args[3]
pkgload::load_all(repo, quiet = TRUE)
source(file.path(repo, "benchmarks/predictive-geometry/mutations.R"))
source(file.path(repo, "benchmarks/provenance.R"))
catalog <- pg_mutation_catalog()
index <- which(vapply(catalog, `[[`, "", "id") == id)
if (length(index) != 1L) stop("Unknown mutation")
mutation <- catalog[[index]]
path <- file.path(repo, "tests/testthat", paste0("test-", mutation$file, ".R"))
source_before <- .crossform_source_tree_digest(repo)
testthat::set_max_fails(Inf)
failures <- function(results) {
  rows <- lapply(results, function(test) {
    bad <- Filter(function(x) inherits(x, "expectation_failure") || inherits(x, "expectation_error"), test$results)
    if (!length(bad)) return(NULL)
    data.frame(test = test$test, type = vapply(bad, function(x) class(x)[1], ""),
      message = vapply(bad, conditionMessage, ""), stringsAsFactors = FALSE)
  })
  value <- do.call(rbind, rows)
  if (is.null(value)) data.frame(test = character(), type = character(), message = character()) else value
}
run <- function(replacements = NULL) {
  instrument <- NULL
  if (!is.null(replacements)) {
    instrument <- pg_mutation_instrument(replacements)
    do.call(testthat::local_mocked_bindings, c(instrument$replacements,
      list(.package = "crossform", .env = environment())))
  }
  results <- testthat::test_file(path, reporter = "silent", stop_on_failure = FALSE)
  list(failures = failures(results), tests = vapply(results, `[[`, "", "test"),
    expectations = sum(vapply(results, function(x) length(x$results), 0L)),
    hits = if (is.null(instrument)) integer() else instrument$counter$hits)
}
baseline <- run()
if (nrow(baseline$failures)) stop("Unmodified primary test file fails; mutation result would be uninterpretable")
replacements <- mutation$transform()
originals <- stats::setNames(lapply(names(replacements), get, envir = asNamespace("crossform")), names(replacements))
control <- run(originals)
if (nrow(control$failures)) stop("No-op instrumentation changes the primary tests")
mutated <- run(replacements)
primary <- mutated$failures[grepl(mutation$primary, mutated$failures$test, fixed = TRUE), , drop = FALSE]
killed <- nrow(primary) > 0L && all(mutated$hits > 0L)
result <- list(id = id, description = mutation$description, symbols = names(replacements),
  primary = mutation$primary, test_file = sub(paste0(repo, "/"), "", path, fixed = TRUE),
  source_digest = source_before, test_digest = digest::digest(file = path, algo = "sha256"),
  original_body_digests = vapply(originals, function(x) digest::digest(body(x), algo = "sha256"), ""),
  mutated_bodies = lapply(replacements, function(x) deparse(body(x), width.cutoff = 120L)),
  baseline = baseline, instrumented_control = control, mutated = mutated,
  primary_failures = primary, killed = killed)
stopifnot(identical(source_before, .crossform_source_tree_digest(repo)))
saveRDS(result, output)
cat(id, if (killed) "KILLED" else "SURVIVED", "primary failures", nrow(primary),
  "actual production calls", sum(mutated$hits), "\n")
if (!killed) quit(status = 2L)
