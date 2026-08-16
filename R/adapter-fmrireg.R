# fmrireg execution adapter ------------------------------------------------

#' Execute a relation plan with fmrireg
#'
#' This installed-consumer adapter independently executes the OLS point
#' relation using `fmrireg::fmri_ols_fit()`. It deliberately returns no
#' residual error channel; analytic second-moment uncertainty therefore
#' refuses rather than borrowing unsupported standard-error output.
#'
#' @param x A [plan_relation()] result using an OLS observation model.
#' @return An `effect_relation_fit` with point-relation capabilities only.
#' @export
fmrireg_relation <- function(x) {
  version <- .require_adapter_version("fmrireg", "0.1.2")
  plan <- .validate_relation_plan(x)
  if (!identical(plan$observation_model$kind, "ols")) {
    .capability_refusal(
      "The certified fmrireg adapter supports OLS relation plans only.",
      capability = "supported_observation_model",
      namespace = "relation_compiler",
      reasons = paste0(
        "The plan declares `", plan$observation_model$kind,
        "`; fmrireg 0.1.2 parity is certified only without whitening."
      ),
      remedies = c(
        "Execute the plan with `estimate()`.",
        "Add and certify a GLS-specific external adapter."
      )
    )
  }
  planned <- .planned_observation_sources(plan)
  extractors <- lapply(plan$partitions, function(partition) {
    design <- plan$design_receipts[[partition]]$design
    target <- plan$design_receipts[[partition]]$lowered_target
    coefficient_operator <- fmrireg::fmri_ols_fit(
      diag(nrow(design)), design
    )$beta
    effect_extractor(
      target %*% coefficient_operator,
      effects = plan$lowered_effects[[partition]]$effect_space,
      estimator = "fmrireg::fmri_ols_fit",
      diagnostics = list(
        adapter_version = version,
        residual_df = nrow(design) - qr(design)$rank
      )
    )
  })
  names(extractors) <- plan$partitions
  relation_value <- relation(
    planned$sources,
    extract = extractors,
    source_dims = planned$dimensions,
    partitions = plan$partitions,
    domain = plan$study$observations$domain,
    capabilities = planned$capabilities,
    provenance = list(
      contract_version = "first-moment-relation-v1",
      relation_plan_id = plan$relation_plan_id,
      adapter = "crossform::fmrireg_relation",
      adapter_package = "fmrireg",
      adapter_version = version,
      design_receipt_ids = vapply(
        plan$design_receipts, `[[`, character(1), "design_receipt_id"
      )
    )
  )
  relation_fit(
    relation_value,
    error_models = NULL,
    provenance = list(
      contract_version = "first-moment-relation-v1",
      relation_plan_id = plan$relation_plan_id,
      adapter = "crossform::fmrireg_relation",
      adapter_version = version,
      analytic_error_channel = "withheld"
    )
  )
}
