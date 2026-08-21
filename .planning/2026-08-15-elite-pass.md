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

## Outcome (2026-08-17)

Branch `elite-pass`, 14 commits ahead of `main` (11 from this pass, 3 from a concurrent session that
added an Rcpp fused pair-query / packed-form kernel). Every pass commit landed with
`R CMD check --as-cran` at 0 ERROR / 0 WARNING (the only WARNING now is the local Homebrew-clang
`-Wfixed-enum-extension` from R's own headers, introduced with the compiled code) and a green suite
(2,402 -> ~4,840 expectations under `testthat::test_local()`; the certified
figure is `PASS 4112` under `R CMD check --as-cran`, which skips 10 blocks —
see `design/certification-report.md`, "Recorded certification metrics").

Fresh-eyes first-hour acceptance: run 1 revelation 8 / friction 6 -> run 2 9 / 8 -> run 3 **9 / 8,
ACCEPT** ("every executable claim is true; every remaining defect was a documentation promise" —
all fixed in the final polish commit).

Real defects found and fixed along the way (all with Monte Carlo or oracle proof):
- plug-in sampling covariance signal term (isotropic surrogate; off by ~support size);
- null-target SE finite-sample bias (raw plug-in tr Sigma^2; ~8% at radius-1.5, ~100% at 500 voxels);
- searchlight radius cubic hang; digest serializeVersion; pattern-matrix frames; bitwise agreement;
- contrast_energy() 7.25 s -> 0.86 s at 20k searchlights (fused first-moment pass).

Left for the maintainer:
- Make the GitHub repo public and enable Pages (workflows are in place; all pkgdown URLs 404 until then).
- `shard-admission.rds` stays UNBOUND (installed `shard` 0.2.1 vs recorded 0.2.0) — decide whether to re-run.
- Certify `fmrireg` 0.2.0 for the adapter (currently refuses non-0.1.2).
- Decide whether the native kernel's summation reordering (last-ulp drift) is within the package's
  bitwise "scheduling" guarantee; benchmarks/README headline numbers need refreshing after the
  native kernel (recertify agent has the new values).
- Follow-ups outside this pass: group-inference path, parallel backend, methods-paper simulation
  study, six-subject Haxby decomposition figure, `correlation_rdm()` view.
