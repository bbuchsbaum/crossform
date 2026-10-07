# crossform: cross-generalized effect geometry for task fMRI

`crossform` answers one question about task-fMRI data: where does an
experimental effect, or a whole representational geometry, reproduce
across independent runs? You declare the condition effects you fitted,
the spatial measurements you want answers at, and which runs must
generalize. The package compiles that request once and reads contrasts,
representational dissimilarity matrices (RDMs), representational
similarity analysis (RSA) coefficients, and crossnobis distances from
the same compiled object rather than from separate pipelines.

## The mental model

Three declarations compile into one plan, and every scientific result is
a query against that plan:

    relation + frame + pairing  --plan_geometry()->  geometry plan  ->  views

- **relation** — run-by-condition effect estimates. Build one from beta
  matrices with
  [`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md),
  or fit one from scan responses with
  [`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
  and
  [`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md).

- **frame** — where results are reported. Choose
  [`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md),
  [`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md),
  [`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md),
  or
  [`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md),
  and compile it against a neural domain with
  [`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md).
  A *measurement* is one such spatial unit; every result has one row per
  measurement.

- **pairing** — which partitions must generalize, declared with
  [`cross_partitions()`](https://bbuchsbaum.github.io/crossform/reference/cross_partitions.md).
  This is bound into the plan's identity, so generalizing over runs and
  over sessions are different estimands.

- **geometry plan** — what
  [`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
  returns. It records the estimand (the quantity you mean to estimate)
  before any neural values are read, and keeps that identity separate
  from execution details such as block size or storage.

- **views** —
  [`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
  [`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
  [`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md),
  and
  [`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md),
  all reading the same plan. Asking a second question costs a query, not
  a refit.

## Start here

- [`vignette("introduction")`](https://bbuchsbaum.github.io/crossform/articles/introduction.md)
  — the guided entry point, from a ready relation through contrasts,
  RDMs, RSA, and standard errors. Nine vignettes ship with the package;
  `vignette(package = "crossform")` lists them. They are present only if
  the package was installed with its vignettes built —
  `remotes::install_local(".", build_vignettes = TRUE)`, or
  `R CMD build` followed by `R CMD INSTALL` — and are otherwise readable
  at <https://bbuchsbaum.github.io/crossform/>.

- [`example_fmri_effects()`](https://bbuchsbaum.github.io/crossform/reference/example_fmri_effects.md)
  — a generated fixture with known truth, so the whole workflow runs
  after installation and can check its own answer.

- [`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
  — compile a relation, frame, and pairing into a plan.

- [`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
  — where a named contrast reproduces across runs, split into the part
  carried by the local mean pattern (`coherent`) and the reproducible
  spatial pattern beyond it (`configuration`).

- [`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) —
  crossvalidated squared distances between conditions, for all pairs or
  only the pairs you name.

- [`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) —
  fixed linear RSA coefficients against named model RDMs.

- [`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md)
  — the crossvalidated Mahalanobis distance under a fixed neural metric.

- [`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
  — within-measurement analytic standard errors, available only when the
  fit retains the required residual channel.

- [`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)
  — capture a refused operation and read the unmet requirements and
  their remedies as data.

- [crossform_conditions](https://bbuchsbaum.github.io/crossform/reference/crossform_conditions.md)
  — the four classes every crossform failure carries, so a caller can
  branch on the cause instead of the prose.

## Scope

- Experimental. Exported names may still change.

- Sequential execution; there is no parallel backend.

- No preprocessing, registration, masking, or hemodynamic response
  modeling. Bring fitted effects, or scan responses with a design.

- Within-participant results only. There is no group-inference path.

- Crossvalidated estimates are retained when negative; they are not
  clipped.

- Uncertainty is refused rather than approximated when the fit cannot
  support the admitted analytic law.

- Every failure is a classed condition. A wrong argument raises
  `effect_input_error`; two objects that do not belong together raise
  `effect_contract_error`; an impossible internal result raises
  `effect_invariant_error` and asks for a bug report; an interpretation
  that cannot be earned from the supplied objects raises
  `effect_capability_refusal`. See
  [crossform_conditions](https://bbuchsbaum.github.io/crossform/reference/crossform_conditions.md).

## See also

Useful links:

- <https://bbuchsbaum.github.io/crossform/>

- <https://github.com/bbuchsbaum/crossform>

- Report bugs at <https://github.com/bbuchsbaum/crossform/issues>

## Author

**Maintainer**: Bradley Buchsbaum <bradley.buchsbaum@gmail.com>
