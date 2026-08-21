# Relation source capabilities -----------------------------------------------
#
# Layer 2 (values). What a source can do is a fact about the source, so this
# file owns both halves of it: the `source_capabilities()` record and its
# validator, and the admission check that refuses a relation whose sources
# cannot be read in bounded blocks.
#
# Both halves arrived here from elsewhere, and for the same reason each time.
# The admission check lived in `R/compiler.R` only because the compiler was
# its first caller, which made every value that admits a relation appear to
# reach into the executor. The record lived in `R/receipt.R` only because the
# receipt is where a capability is finally written down -- but a receipt is a
# record of what an execution used, and `relation.R`, `relation-fit.R`, and
# `study-facts.R` all have to *state* capabilities long before any receipt
# exists. Defining the value here lets them say so without calling the record
# that quotes them back.

#' Declare source execution capabilities
#'
#' `source_capabilities()` states what an out-of-memory relation source can
#' actually do, so the compiler can choose a bounded execution route and
#' record the source revision in the receipt. Declare it when supplying a
#' custom source such as [file_matrix_source()].
#'
#' @param block_read Whether bounded feature-block reads are supported.
#' @param reopenable Whether a fresh read-only handle can be opened safely.
#' @param thread_safe Whether concurrent reads through one handle are supported.
#' @param stable_revision A strong immutable source revision or checksum.
#' @return An `effect_source_capabilities` value with the three logical flags
#'   `$block_read`, `$reopenable`, `$thread_safe`, and the
#'   `$stable_revision` identifier bound into the execution receipt.
#' @seealso [file_matrix_source()], which carries these capabilities, and
#'   [compute_policy()], which is checked against them before execution.
#' @family relation planning and fitting
#' @examples
#' # A block-readable, reopenable source identified by a content checksum.
#' revision <- paste0("sha256:", paste(rep("a", 64), collapse = ""))
#' capabilities <- source_capabilities(
#'   block_read = TRUE, reopenable = TRUE, stable_revision = revision
#' )
#' capabilities$block_read
#' capabilities$thread_safe
#'
#' # The revision must be a strong identifier: a mutable label such as a file
#' # modification time is refused, because receipts must stay verifiable.
#' weak <- try(
#'   source_capabilities(TRUE, stable_revision = "2026-08-15"), silent = TRUE
#' )
#' conditionMessage(attr(weak, "condition"))
#' @export
source_capabilities <- function(block_read, reopenable = FALSE,
                                thread_safe = FALSE, stable_revision) {
  flags <- list(
    block_read = block_read,
    reopenable = reopenable,
    thread_safe = thread_safe
  )
  if (!all(vapply(flags, function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }, logical(1)))) {
    .input_error("Source capability flags must each be TRUE or FALSE.")
  }
  if (!.strong_sha256(stable_revision)) {
    .input_error(paste0(
      "`stable_revision` must be a sha256 identifier with 64 hexadecimal ",
      "digits."
    ))
  }
  value <- structure(
    c(flags, list(stable_revision = stable_revision)),
    class = "effect_source_capabilities"
  )
  .validate_source_capabilities(value)
}

.validate_source_capabilities <- function(value) {
  expected <- c("block_read", "reopenable", "thread_safe", "stable_revision")
  if (!inherits(value, "effect_source_capabilities") ||
      !identical(names(value), expected)) {
    .input_error("Source capabilities are missing required canonical fields.")
  }
  flags <- value[c("block_read", "reopenable", "thread_safe")]
  if (!all(vapply(flags, function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }, logical(1)))) {
    .input_error("Source capability flags must each be TRUE or FALSE.")
  }
  if (!.strong_sha256(value$stable_revision)) {
    .input_error(
      "Source capability revision must be one strong sha256 identifier."
    )
  }
  structure(as.list(value), class = "effect_source_capabilities")
}

.relation_source_capabilities <- function(x) {
  if (is.null(x$capabilities)) {
    .input_error(paste0(
      "Opaque relation sources require explicit `source_capabilities()` ",
      "before execution."
    ))
  }
  capabilities <- lapply(x$capabilities, .validate_source_capabilities)
  if (!all(vapply(capabilities, function(value) isTRUE(value$block_read),
    logical(1)))) {
    .input_error("All relation sources must support bounded block reads.")
  }
  capabilities
}
