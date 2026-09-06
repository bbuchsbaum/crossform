#' crossform: cross-generalized effect geometry for task fMRI
#'
#' Crossform gives univariate contrasts, multivariate distances and fixed linear
#' RSA, and population summaries one declared crossvalidated bilinear geometry,
#' then decomposes each scale's reproducible effect into coherent and
#' configurational estimand components.
#'
#' You declare condition effects, spatial measurements, independent partitions,
#' and --- for a population --- realized transports and a group model. The
#' package compiles that estimand before reading outcomes. Supported contrasts,
#' RDMs, fixed linear RSA, crossnobis distances, and population summaries then
#' share the same geometry. Coherent and configurational name additive estimand
#' components, not separate biological mechanisms.
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
#'   relation through contrasts, RDMs, RSA, and standard errors. Nine vignettes
#'   ship with the package; `vignette(package = "crossform")` lists them. They
#'   are present only if the package was installed with its vignettes built ---
#'   `remotes::install_local(".", build_vignettes = TRUE)`, or `R CMD build`
#'   followed by `R CMD INSTALL` --- and are otherwise readable at
#'   <https://bbuchsbaum.github.io/crossform/>.
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
#' - [crossform_conditions] --- the four classes every crossform failure
#'   carries, so a caller can branch on the cause instead of the prose.
#'
#' @section Scope:
#'
#' - Experimental. Exported names may still change.
#' - Sequential execution; there is no parallel backend.
#' - No preprocessing, registration, masking, or hemodynamic response
#'   modeling. Bring fitted effects, or scan responses with a design.
#' - Population point estimates, classical and HC3 pointwise intervals,
#'   null-imposed wild bootstrap, decomposition, prevalence, heterogeneity,
#'   influence, and support diagnostics are implemented under their declared
#'   conditional targets.
#' - No simultaneous, transport-marginal, arbitrary cross-node, or universally
#'   calibrated population inference is implied.
#' - Crossvalidated estimates are retained when negative; they are not clipped.
#' - Uncertainty is refused rather than approximated when the fit cannot
#'   support the admitted analytic law.
#' - Every failure is a classed condition. A wrong argument raises
#'   `effect_input_error`; two objects that do not belong together raise
#'   `effect_contract_error`; an impossible internal result raises
#'   `effect_invariant_error` and asks for a bug report; an interpretation
#'   that cannot be earned from the supplied objects raises
#'   `effect_capability_refusal`. See [crossform_conditions].
#'
#' @keywords internal
#' @useDynLib crossform, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @name crossform-package
#' @aliases crossform
"_PACKAGE"
