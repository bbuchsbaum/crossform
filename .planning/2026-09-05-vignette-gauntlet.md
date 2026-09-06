# Vignette pedagogy gauntlet

User request: improve the vignettes' balance of pedagogy and explication;
make them clear, lucid, clean and tidy. This follow-up does not include
committing, pushing or publishing the revised documentation.

Baseline: commit `4484c7fc6cbb63295e051ad7670113de4cb43b8e`.
Source archive and rendered baseline:
`/private/tmp/crossform-vignette-gauntlet-94nfvarn/`.
Mote: `bd-01M1T8A9RCF18NXNSQPS4RCB3N`.

Review criteria:

- A concrete reader question, prerequisites and an early interpreted result.
- A coherent progression from core workflow to optional deeper explanation.
- Scientific precision: signed estimation, fixed queries, latent PSD
  description and independent predictive gain remain distinct.
- Executable evidence supports the numerical claims near it.
- Readable rendered code, tables, equations, figures and cross-links.
- Enough explanation to understand the result, without repetitive history,
  implementation bookkeeping or unsupported claims.

All fifteen vignettes are in scope. The root builder handles introduction,
interpretation, predictive geometry, common-geometry derivation and novelty.
Two bounded builders handle the other ten guides. A separate critic received
the requirements and actual baseline artifacts with fresh context, and then
the actual revised artifacts without builder self-assessments. Full-page
comparison is independent but not blinded; no human comprehension study is
claimed. Rendering is serial, in fresh R processes.

The baseline review identified overstatements about null sign probabilities,
observed noise floors and independence; delayed first results; literal
unrendered mathematics; wide output and evidence tables; development-history
digressions; and a figure asset referring outside the built site. These are
being revised while retaining scientific checks and the source/test/artifact
certification boundary. Final review and verification results follow below.

## Completed revision

All fifteen guides were revised. The articles navigation in `_pkgdown.yml`
was also reordered around reader tasks; its reference configuration is unchanged.
No production R, tests, reference manuals or certification artifacts changed.
The pre-existing exemplar edits were preserved.

| Guide | Main improvement |
|---|---|
| Introduction | Earlier interpreted result, explicit independence assumptions, optional simulation setup and task-directed next steps. |
| Interpreting results | Correct null and radius interpretations, clear selection caveats and readable capability diagnostics. |
| From observations | First fit before optional derivations; reuse one dataset across equivalent routes. |
| neuroim2 inputs | Precise membership, overlap and backprojection interpretation. |
| Migrating from rmvpa | Explicit equivalences and differences; distinguish scan residuals from beta uncertainty. |
| Failure gallery | Executable refusals with supported remedies. |
| Matched interpretability | Portable figure asset, clearer caption and explicit scale and availability checks. |
| Conservative frames | Coherent walkthrough; quiet executable assertions with results and conservation gaps visible. |
| Predictive geometry | Model directions explained, rank and penalty distinguished, signed held-out gain interpreted. |
| Evidence pairing | Covariance before correlation, explicit orientation and rank qualifications. |
| Population forms | First group result before variants; optional helper details and bounded calibration claims. |
| Correlation and distance | Concrete comparison of magnitude and pattern correlation, with normalization assumptions. |
| Extension guide | Current protocol, reproducible source identity and readable diagnostics. |
| Common geometry | Correctly rendered mathematics, compact numerical evidence and a three-column estimator table. |
| Novelty | Shorter claims, preserved prior attribution, expandable evidence ledger and explicit signed/latent boundary. |

The independent critic reviewed baseline and revised source, rendered HTML,
mathematics, tables, links and image artifacts. Iterations addressed both
scientific overstatements and presentation defects. Final verdict: accepted;
all fifteen improve over baseline, with no remaining must-fix identified.
The comparison was not blinded and was not a human comprehension study.

## Verification

- All 15 revised articles rendered successfully in fresh R processes. The
  final three touch-ups were rendered again after the last edits.
- Focused documentation, adapter, predictive example, evidence ledger and
  terminology tests: **318 passed; 0 failures, errors, test warnings or skips**.
- HTML audit: **359 article/anchor links and 6 local image references**, with
  no missing targets, MathML errors, duplicate chunk labels, empty chunks or
  unintended control characters.
- All **11 digest-bound certification artifacts** still bind to current R
  source; `shard-admission.rds` remains intentionally unbound.
- Source digest:
  `sha256:10bfc4506bcf9c110f8ac10a0bfba425fb0f3eb4c5975e608bd1cd095690bf12`.
- `git diff --check` passes. This documentation-only change did not require
  rerunning the full numerical suite or regenerating certification evidence.
- Live viewport inspection was unavailable: the product browser could not
  start. No user Chrome profile was used. The final browser audit found no
  automated top-level browser processes.

Local render and verification artifacts are in
`/private/tmp/crossform-vignette-gauntlet-94nfvarn/`: `revised-site/articles/`,
`final-render.log`, `final-touchups-render.log`, `tests.log`,
`test-results.rds`, `binding.log` and `check-html.py`.
These are local review artifacts, not a published website. API reference pages
were not rebuilt; the link audit covers articles and their anchors.

The review concluded with changes uncommitted. On 2026-09-06, the user
authorized committing and pushing this verified documentation slice.
