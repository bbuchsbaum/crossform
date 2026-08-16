# crossform elite pass — plan

Date: 2026-08-15. Branch: `elite-pass`. Baseline: commit 3480ec3 (0 ERR / 0 WARN / 1 NOTE; 2402 tests; covr 87.7%).

Source: eight-lens audit (check, API, numerics, architecture, tests, docs, DX, science). Goal: a
package that is compelling in the first hour and elite at every level — correctness, API, docs,
tests, architecture, presentation, scientific positioning.

## Waves (each ends with tests + R CMD check green, then a commit)

0. **Hygiene.** Untrack `.mote/`; move root design docs to `design/`; move novelty / failure
   gallery / correlation policy to `vignettes/articles/`; fix README benchmark numbers; reconcile
   certification-report contradictions; DESCRIPTION fields (Language, Roxygen, Config/Needs/website,
   `>=` pins); roxygen-managed NAMESPACE; `inst/WORDLIST`; R-CMD-check + pkgdown workflows.
1. **Confirmed bugs.** Plug-in sampling covariance signal term (isotropic surrogate → real residual
   covariance) + Monte Carlo tests at two normalizations with a non-spherical metric + disclose Δ
   bias; pin digest `serializeVersion` via one helper; clamp searchlight offset enumeration to the
   domain bounding box; `numerical_agreement("scheduling")` bitwise.
2. **First hour.** Rename masking exports (`contrast`→`contrast_energy`, `estimate`→
   `estimate_relation`, `events`→`observation_events`, `voxels`→`voxelwise`, `geometry`→
   `materialize_geometry`, `pair_contrast`→`coupling_contrast`); print/format for every
   user-facing class and
   the refusal condition; `?crossform` page; `@examples` + `@family` on every export; errors that
   report what was received; `rsa()` intercept message; runnable BYO-data path; README rewrite with
   a diagram; vignette order + neuroim2 article.
3. **Tests.** Bind certification artifacts to source SHA and ship gate files so they run under
   check; wire up or delete dead oracle helpers; make sampling-law tests exercise the package;
   condition-class assertions for refusal paths; dedicated tests for operations/query-structured/
   coupling-plan/task.
4. **Architecture.** Extract leaf primitives (break kernel↔views); move executor out of compiler;
   break plan↔compiler recursion; `.check_*` guard helpers; error taxonomy; delete scheduler/
   checkpoint dead machinery; collapse coupling reducers.
5. **Science docs.** cvMANOVA / Walther / Nili / Kriegeskorte / Framed-RSA citations and a
   "relation to cvMANOVA" section; unbiasedness statement; coherence-fraction caveats
   (sphere size, frame weighting, selection); `correlation_rdm()` if the contract permits.

Out of scope for this pass (follow-ups): group inference path, methods-paper simulation study,
six-subject Haxby decomposition figure, making the GitHub repo public / enabling Pages (user action).
