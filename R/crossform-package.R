#' crossform: cross-generalized effect geometry for task fMRI
#'
#' `crossform` answers one question about task-fMRI data: where does an
#' experimental effect, or a whole representational geometry, reproduce across
#' independent runs? You declare the condition effects you fitted, the spatial
#' measurements you want answers at, and which runs must generalize. The
#' package compiles that request once and reads contrasts, representational
#' dissimilarity matrices (RDMs), representational similarity analysis (RSA)
#' coefficients, and crossnobis distances from the same compiled object rather
#' than from separate pipelines.
#'
#' @section The mental model:
#'
#' Three declarations compile into one plan, and every scientific result is a
#' query against that plan:
#'
#' ```
#' relation + frame + pairing  --plan_geometry()->  geometry plan  ->  views
#' ```
#'
#' - **relation** --- run-by-condition effect estimates. Build one from beta
#'   matrices with [relation()], or fit one from scan responses with
#'   [plan_relation()] and [estimate_relation()].
#' - **frame** --- where results are reported. Choose [searchlights()],
#'   [regions()], [voxelwise()], or [whole_brain()], and compile it against a
#'   neural domain with [compile_frame()]. A *measurement* is one such spatial
#'   unit; every result has one row per measurement.
#' - **pairing** --- which partitions must generalize, declared with
#'   [cross_partitions()]. This is bound into the plan's identity, so
#'   generalizing over runs and over sessions are different estimands.
#' - **geometry plan** --- what [plan_geometry()] returns. It records the
#'   estimand (the quantity you mean to estimate) before any neural values are
#'   read, and keeps that identity separate from execution details such as
#'   block size or storage.
#' - **views** --- [contrast_energy()], [rdm()], [rsa()], and [crossnobis()],
#'   all reading the same plan. Asking a second question costs a query, not a
#'   refit.
#'
#' @section Start here:
#'
#' - `vignette("introduction")` --- the guided entry point, from a ready
#'   relation through contrasts, RDMs, RSA, and standard errors.
#' - [example_fmri_effects()] --- a generated fixture with known truth, so the
#'   whole workflow runs after installation and can check its own answer.
#' - [plan_geometry()] --- compile a relation, frame, and pairing into a plan.
#' - [contrast_energy()] --- where a named contrast reproduces across runs,
#'   split into the part carried by the local mean pattern (`coherent`) and the
#'   reproducible spatial pattern beyond it (`configuration`).
#' - [rdm()] --- crossvalidated squared distances between conditions, for all
#'   pairs or only the pairs you name.
#' - [rsa()] --- fixed linear RSA coefficients against named model RDMs.
#' - [crossnobis()] --- the crossvalidated Mahalanobis distance under a fixed
#'   neural metric.
#' - [rdm_sampling_covariance()] --- within-measurement analytic standard
#'   errors, available only when the fit retains the required residual channel.
#' - [catch_refusal()] --- capture a refused operation and read the unmet
#'   requirements and their remedies as data.
#'
#' @section Scope:
#'
#' - Experimental. Exported names may still change.
#' - Sequential execution; there is no parallel backend.
#' - No preprocessing, registration, masking, or hemodynamic response
#'   modeling. Bring fitted effects, or scan responses with a design.
#' - Within-participant results only. There is no group-inference path.
#' - Crossvalidated estimates are retained when negative; they are not clipped.
#' - Uncertainty is refused rather than approximated when the fit cannot
#'   support the admitted analytic law.
#'
#' @keywords internal
#' @name crossform-package
#' @aliases crossform
"_PACKAGE"
