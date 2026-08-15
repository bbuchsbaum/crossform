# External first-moment compiler conformance -------------------------------

.require_adapter_version <- function(package, supported) {
  if (!requireNamespace(package, quietly = TRUE)) {
    .capability_refusal(
      sprintf("Adapter package `%s` is not installed.", package),
      capability = "installed_compiler_adapter",
      namespace = "relation_compiler",
      reasons = paste0("Missing optional package: ", package, "."),
      remedies = paste0("Install the tested `", package, "` version ",
        supported, ".")
    )
  }
  observed <- as.character(utils::packageVersion(package))
  if (!identical(observed, supported)) {
    .capability_refusal(
      sprintf("Adapter `%s` has not been certified for version %s.",
        package, observed),
      capability = "supported_compiler_version",
      namespace = "relation_compiler",
      reasons = paste0("Certified version: ", supported,
        "; installed version: ", observed, "."),
      remedies = c(
        paste0("Use `", package, "` ", supported, "."),
        "Add and certify a version-specific adapter before granting conformance."
      )
    )
  }
  observed
}

#' Inspect first-moment compiler conformance
#'
#' The conformance court validates the portable receipt required by
#' `first-moment-relation-v1`. It reports construction evidence rather than
#' trusting an adapter name or package version.
#'
#' @param x A [plan_relation()] result.
#' @return A data frame with one row per partition and one boolean per required
#'   conformance field.
#' @export
compiler_conformance <- function(x) {
  x <- .validate_relation_plan(x)
  values <- lapply(x$partitions, function(partition) {
    receipt <- .validate_design_receipt(x$design_receipts[[partition]])
    source <- x$study$observations$capabilities[[partition]]
    data.frame(
      partition = partition,
      semantic_identity = grepl(
        "^(design-model|raw-design-model)-sha256:[[:xdigit:]]{64}$",
        receipt$design_model_id
      ),
      regressor_axis = identical(
        receipt$coefficient_axis, colnames(receipt$design)
      ),
      condition_lowering = identical(
        receipt$coefficient_axis, colnames(receipt$lowered_target)
      ),
      row_lineage = is.data.frame(receipt$row_lineage) &&
        nrow(receipt$row_lineage) == nrow(receipt$design),
      rank_and_aliases = is.integer(receipt$rank) &&
        is.character(receipt$aliases),
      censor_accounting = is.list(receipt$censoring) &&
        identical(receipt$censoring$retained_rows,
          as.integer(receipt$row_lineage$source_row)),
      solver_diagnostics = is.character(receipt$solver) &&
        length(receipt$solver) == 1L,
      whitening_provenance = is.list(receipt$observation_whitener) &&
        .strong_sha256(receipt$observation_whitener$signature),
      source_revision = .strong_sha256(source$stable_revision),
      portable_receipt = isTRUE(
        receipt$capabilities$portable_design_receipt
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, values)
}
