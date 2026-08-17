# fmrireg execution adapter ------------------------------------------------

#' Execute a relation plan with fmrireg
#'
#' This installed-consumer adapter independently executes the OLS point
#' relation using `fmrireg::fmri_ols_fit()`. It deliberately returns no
#' residual error channel; analytic second-moment uncertainty therefore
#' refuses rather than borrowing unsupported standard-error output.
#'
#' @param x A [plan_relation()] result using an OLS observation model.
#' @return An `effect_relation_fit` with point-relation capabilities only: its
#'   `$error_models` are all `NULL`, so `residual_blocks`,
#'   `effect_covariance`, and `residual_df` are `FALSE` for every partition.
#'   Its `$provenance` records the `relation_plan_id`, the adapter, and the
#'   pinned `adapter_version`.
#' @family relation planning and fitting
#' @seealso [estimate_relation()], the in-package execution verb that does
#'   return a residual channel, and [relation_fit_capabilities()] to see the
#'   difference.
#' @examples
#' # This adapter is certified against exactly one fmrireg version and refuses
#' # any other, so the example runs only under that version.
#' if (requireNamespace("fmrireg", quietly = TRUE) &&
#'     identical(as.character(utils::packageVersion("fmrireg")), "0.1.2")) {
#'   set.seed(1)
#'   domain <- abstract_domain(3L, id = "fmrireg-example")
#'   index <- observation_index(paste0("scan-", 1:4), "run-1")
#'   facts <- study(observations(
#'     list(`run-1` = matrix(rnorm(12), 4L, 3L)), list(`run-1` = index), domain
#'   ))
#'   design <- cbind(face = c(1, 0, 1, 0), body = c(0, 1, 0, 1))
#'   rownames(design) <- paste0("scan-", 1:4)
#'   target <- rbind(`face-body` = c(1, -1))
#'   colnames(target) <- colnames(design)
#'   plan <- plan_relation(
#'     facts, raw_design_model(list(`run-1` = design)), raw_effect_map(target),
#'     observation_model("ols", sampling_unit = "scan")
#'   )
#'
#'   # An independent execution of the same point relation: use it to check
#'   # parity against the in-package estimator.
#'   external <- fmrireg_relation(plan)
#'   print(all.equal(
#'     relation_block(external, "run-1", 1:3),
#'     relation_block(estimate_relation(plan), "run-1", 1:3)
#'   ))
#'
#'   # Parity is for point estimates only; no residual channel is claimed.
#'   print(relation_fit_capabilities(external)$residual_blocks)
#'
#'   # A plan outside the certified OLS slice is refused, not approximated.
#'   gls <- plan_relation(
#'     facts, raw_design_model(list(`run-1` = design)), raw_effect_map(target),
#'     observation_model(
#'       "fixed_gls", sampling_unit = "scan", whitener = diag(4)
#'     )
#'   )
#'   print(catch_refusal(fmrireg_relation(gls))$capability)
#' }
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
        "Execute the plan with `estimate_relation()`.",
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
