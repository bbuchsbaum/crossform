# Relation source capabilities -----------------------------------------------
#
# Layer 2 (values). Admitting a relation's sources is a statement about the
# sources, not an act of compilation: it validates each declared
# `source_capabilities()` record and refuses anything that cannot be read in
# bounded blocks. It lived in `R/compiler.R` only because the compiler was the
# first caller, which made every value that admits a relation appear to reach
# into the executor.

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
