# fmrireg execution adapter ------------------------------------------------

# What an external executor has to build for itself, built here from the
# published surface only: the plan's documented `$study` and `$retained_rows`,
# the receipts `relation_plan_receipts()` returned, and
# `source_capabilities()`, one of the five sanctioned developer entry points.
#
# A censored plan does not read the rows it dropped, so the source an executor
# hands to `relation()` is the planned source -- the study's reader restricted
# to the rows the plan retained -- and not the study's own. That restriction is
# a *different* source from the one the study declared, which is why the
# capabilities below are re-declared rather than copied: the derived source is
# not reopenable, is not thread safe, and its revision has to bind the row
# selection and the design receipt that produced it, or two plans over one
# study would claim the same content.
#
# The recipe is the one `estimate_relation()` uses, deliberately, so that the
# two executions of a plan agree about what they read as well as what they
# compute.
.fmrireg_planned_sources <- function(plan, receipts) {
  observations <- plan$study$observations
  named <- function(values) stats::setNames(values, plan$partitions)
  sources <- named(lapply(plan$partitions, function(partition) {
    source <- observations$sources[[partition]]
    rows <- plan$retained_rows[[partition]]
    force(source)
    force(rows)
    function(features) source$read(features)[rows, , drop = FALSE]
  }))
  dimensions <- named(lapply(plan$partitions, function(partition) {
    c(length(plan$retained_rows[[partition]]), observations$n_features)
  }))
  capabilities <- named(lapply(plan$partitions, function(partition) {
    declared <- observations$capabilities[[partition]]
    source_capabilities(
      block_read = declared$block_read,
      reopenable = FALSE,
      thread_safe = FALSE,
      stable_revision = .sha256_signature(list(
        source_revision = declared$stable_revision,
        retained_rows = plan$retained_rows[[partition]],
        design_receipt_id = receipts[[partition]]$design_receipt_id
      ))
    )
  }))
  list(sources = sources, dimensions = dimensions, capabilities = capabilities)
}

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
  version <- adapter_version_certificate("fmrireg", "0.1.2")
  # `relation_plan_receipts()` is the public verb that makes a plan prove
  # itself: it runs the plan validator -- rebuilding the plan from its own
  # inputs and refusing anything that does not come back identical -- and
  # hands over the validated per-partition receipts this adapter then reads.
  # The plan that survives that call is `x`, unchanged, by definition.
  receipts <- relation_plan_receipts(x)
  plan <- x
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
  planned <- .fmrireg_planned_sources(plan, receipts)
  extractors <- lapply(plan$partitions, function(partition) {
    design <- receipts[[partition]]$design
    target <- receipts[[partition]]$lowered_target
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
        receipts, `[[`, character(1), "design_receipt_id"
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
