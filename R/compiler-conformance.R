# External first-moment compiler conformance -------------------------------

#' Certify the installed version of an adapter's upstream package
#'
#' `vignette("crossform-extending")` places two obligations on an adapter
#' author, and this is the first of them: certify against the *installed*
#' version of the package you adapt, and record what you certified against.
#' The obligation is part of the protocol, so the refusal that discharges it
#' is too --- an adapter outside this package raises the same two refusals, in
#' the same namespace, as the two adapters shipped here.
#'
#' Call it first, before touching the upstream package, and record the
#' returned string in the provenance of whatever you build. A version other
#' than the certified one is refused rather than attempted: an untested
#' upstream release that still runs is the failure mode this exists to
#' prevent.
#'
#' @param package Name of the upstream package the adapter is certified
#'   against.
#' @param supported The one version string the adapter has been tested with.
#' @return The installed version of `package`, as a character string equal to
#'   `supported`. Anything else refuses:
#'   `installed_compiler_adapter` when `package` is not installed, and
#'   `supported_compiler_version` when the installed version is not the
#'   certified one. Both refusals carry the `relation_compiler` namespace and
#'   name the certified version in their remedies.
#' @family relation planning and fitting
#' @seealso [compiler_conformance()], which checks the receipt an adapter
#'   produced rather than the version that produced it, and
#'   [catch_refusal()] to inspect either refusal.
#' @examples
#' # crossform certifies its own fmridesign adapter against exactly one
#' # version; this is the call that enforces it.
#' if (requireNamespace("fmridesign", quietly = TRUE)) {
#'   installed <- as.character(utils::packageVersion("fmridesign"))
#'   print(identical(
#'     adapter_version_certificate("fmridesign", installed), installed
#'   ))
#' }
#'
#' # An uncertified version is refused, not attempted. The refusal names the
#' # missing capability, so a caller can branch on the cause.
#' refusal <- catch_refusal(
#'   adapter_version_certificate("stats", "0.0.0-never-released")
#' )
#' refusal$capability
#' refusal$remedies
#'
#' # A package that is not installed at all refuses differently.
#' catch_refusal(
#'   adapter_version_certificate("crossformNotAPackage", "1.0.0")
#' )$capability
#' @export
adapter_version_certificate <- function(package, supported) {
  .check_string(package, "package")
  .check_string(supported, "supported")
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
#' @return A data frame with one row per partition: `partition` plus the
#'   logical columns `semantic_identity`, `regressor_axis`,
#'   `condition_lowering`, `row_lineage`, `rank_and_aliases`,
#'   `censor_accounting`, `solver_diagnostics`, `whitening_provenance`,
#'   `source_revision`, and `portable_receipt`.
#' @family relation planning and fitting
#' @seealso [relation_plan_receipts()] for the underlying receipts,
#'   [plan_relation()] for the plan, and [study_capabilities()] for the
#'   equivalent report on the study facts.
#' @examples
#' set.seed(1)
#' domain <- abstract_domain(3L, id = "conformance-example")
#' index <- observation_index(paste0("scan-", 1:4), "run-1")
#' facts <- study(observations(
#'   list(`run-1` = matrix(rnorm(12), 4L, 3L)), list(`run-1` = index), domain
#' ))
#' design <- cbind(face = c(1, 0, 1, 0), body = c(0, 1, 0, 1))
#' rownames(design) <- paste0("scan-", 1:4)
#' target <- rbind(`face-body` = c(1, -1))
#' colnames(target) <- colnames(design)
#' plan <- plan_relation(
#'   facts, raw_design_model(list(`run-1` = design)), raw_effect_map(target),
#'   observation_model("ols", sampling_unit = "scan")
#' )
#'
#' # Every field is evidence read back off the receipt, so this holds for an
#' # external compiler exactly as it does for a hand-built design.
#' conformance <- compiler_conformance(plan)
#' t(conformance)
#' all(as.matrix(conformance[-1L]))
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
