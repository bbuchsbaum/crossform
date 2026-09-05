# Predictive geometry documentation decisions

## Research brief

The primary reader has condition-by-feature estimates from independent runs
and external representational models. A secondary reader builds analysis
software and needs storage, provenance and refusal contracts. Crossform
compiles declared cross-partition geometry into spatial measurements; this
addition learns a model-supported prediction and evaluates it independently.

The three ideas are: model spectra impose directional preferences; fitting
creates a biased PSD prediction from signed geometry; held-out gain measures
improvement over zero with the fitted norm charged once. The ordinary route
is `model_basis()` -> `fit_geometry()` -> `score_geometry()`. The model input
rank and neural response rank have separate meanings.

Actual source-origin lineage, effect meanings, neural metric, frame and
generalization axis constrain valid evaluation. Origins declare a sampling
assumption and cannot be inferred from numeric equality or partition labels.
Memory/block persistence is explicit. The current executor admits independent
same-condition undirected forms, identity or fixed SPD neural metrics, and
Frobenius loss. New-condition transfer and generic GLS remain unsupported.

The inspected tests confirm dense numerical parity, signed mode evidence,
single-offset component readout, origin overlap refusal, persistent records,
bounded workspace and training-only selection. The larger source-bound
performance and certification gates belong to G20/G22. No unrun gate is
described as current certification.

## Reader journey and page contracts

| Page | Reader question | Evidence | Next route |
|---|---|---|---|
| README advanced paragraph | Can I evaluate a learned model-supported geometry? | Names the actual three verbs and estimand | Predictive geometry vignette |
| `predictive-geometry` vignette | How do I fit on two runs and interpret independent evidence on two others? | Executable four-condition example and gain checks | Fit/score reference and standalone script |
| `model_basis`, `fit_geometry`, `score_geometry` reference | What exactly can I pass, retain and read? | Parameter and result contracts; runnable examples | Vignette and related verbs |
| Exemplar README/script | Can I run and adapt the complete workflow without downloaded data? | Kernel/features/squared-RDM equivalence, modes, rank-zero and refusal checks | Vignette and storage parameters |
| Migration design | How do earlier model-coordinate APIs change? | Explicit retained estimands and changed admission | Normative contract |
| Introduction and rMVPA migration vignette | When should I move from fixed RSA to adaptive geometry prediction? | Separate the coefficient and predictive-gain questions | Predictive geometry vignette |
| Interpreting results | What do amplitude, signed mode evidence and gain mean? | Explicit units, null expectation, selection boundary and worked-example link | Predictive geometry vignette |
| Conservative frames | Does a rank budget establish independent prediction? | Executable clipped-negative/truncated-positive mass accounting; frozen component evidence distinction | Predictive geometry vignette |
| Common-geometry equivalence | Where does nonlinear learning rejoin the fixed-query algebra? | Frozen query plus offset and conditional risk identity, outside the fixed RSA theorem | Predictive geometry vignette |
| Evidence pairing | Which observations must be independent? | Distinguish factors inside each cross-product from training/evaluation separation | Four-run predictive example |
| Novelty ledger | Which model-coordinate and predictive claims have evidence? | Keep descriptive/projector trace energy separate from penalized squared-error gain and its source-bound validation | Predictive certification and public guide |

The reference index groups the three predictive verbs together. The existing
descriptive `model_geometry()` reader stays in its own advanced section.
Internal fold scheduling, cache accounting and certification records stay in
design documents; they are not the newcomer's first route.

## API friction and decisions

| Workflow and evidence | Reader cost | Decision and compatibility |
|---|---|---|
| Independent score requires `provenance$observation_origins` | A sampling assumption needs several explicit fields | Show one complete manifest. Keep declaration strict; a future source-ingestion helper could construct it without changing the contract. |
| Fit and descriptive model reader both learn PSD forms | Similar names can suggest interchangeable outputs | Keep each estimand, separate navigation, and state trace energy versus squared-error gain. Ordinary RSA does not replace nonnegative additive fitting. |
| Model input rank and fitted rank coexist | One word denotes two scientific choices | Document rank at each verb; no implicit tuning on final test outcomes. |
| Group views formerly reordered or combined repeated requested rows | Downstream joins and tabulation become surprising | Preserve requested row occurrences and stable empty schemas; numerical gain is unchanged. |
| Nested selection and cache orchestration are internal | Users cannot assume a stable high-level tuning API | Document training/inner/test separation without advertising private functions as public. A public selection vocabulary needs a separate usability decision. |

The public fixed-split workflow is the release-facing example. Broad claims
about adaptive tuning interfaces, condition transfer or platform performance
would still exceed the public contract; they are intentionally absent.

## Vignette follow-up verification, 2026-09-05

The predictive guide now executes squared-RDM input parity and the frozen
coherent/configuration readout with one prediction cost. Its mode table shows
amplitudes 2/1, independent evidence 2/0 and gains 4/-1; the total gain is 3.
The conservative-frames guide executes a rank-one descriptive projection and
checks clipped negative mass, truncated positive mass and their sum separately.
The introduction, interpretation, common-geometry equivalence, evidence-pairing,
rMVPA migration and novelty guides now lead readers to the same public
workflow, with its scientific boundaries stated at the relevant decision.
The old novelty entry remains a projector-energy result and is explicitly
distinguished from predictive gain.

All eight edited articles rendered serially through pkgdown in fresh R
processes using the installed package from the final implementation check.
Every evaluated vignette assertion passed. The predictive article was rendered
again after its final mode-table formatting change. HTML inspection checked
the mode table, spectral formula, MathML error nodes, internal anchors, and
article/reference routes against the source inventory. This was local rendering;
the site was not published and remote URLs were not revalidated.

The focused source tests (`predictive-example`, `common-geometry-vignette`,
`conservative-frames-vignette`) passed 73 assertions with zero failures, errors,
test warnings or skips. The final render also executes the predictive checks
after the table-formatting change. `git diff --check` passed. No production
R code, test files or certification artifacts were changed in this follow-up;
all 11 digest-bound artifacts still match
`sha256:10bfc4506bcf9c110f8ac10a0bfba425fb0f3eb4c5975e608bd1cd095690bf12`.

Local rendered articles, render/test scripts and verification logs are under
`/private/tmp/crossform-vignette-followup-bbd2grmt/`. The Mote follow-up is
`bd-01M1RWJFD1G73VMGTVGAC6FP0D`. Changes remain uncommitted.
