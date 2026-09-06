# Evidence status and non-goals

Status: normative claim ledger

Ledger version: `evidence-status-v1`

Date: 2026-08-21

This ledger states what each headline Crossform claim has earned. The canonical
rows are stored in
`inst/extdata/certification/evidence-status-ledger.csv`; this document defines
their vocabulary and promotion rules. The companion
`inst/extdata/certification/evidence-claim-registry.csv` assigns exactly one
owner and current status to each claim, while
`design/evidence-promotion-history.md` preserves dated status changes. A claim
may have several evidence rows but exactly one registry row.
For example, the common-geometry identity has a proof, an internal oracle, and
one external parity case. Those rows do not combine into a broader claim than
their stated boundaries.

## 1. Non-goals

Crossform does not provide or claim to unify these tasks:

1. **Preprocessing.** It does not perform motion correction, distortion
   correction, normalization, registration, masking, smoothing, denoising, or
   quality control for raw images.
2. **Universal HRF or GLM modeling.** The ingestion layer accepts declared
   designs and a bounded set of adapters. It is not a universal HRF language,
   first-level modeling suite, BIDS workflow, or substitute for checking model
   misspecification.
3. **Classification and prediction.** A fixed contrast geometry is not a
   classifier, decoding accuracy, a trained decision rule, or a prediction
   benchmark.
4. **Nonlinear or adaptive distances.** Correlation and cosine normalization,
   rank transforms, learned kernels, locally selected features, and embeddings
   fitted on evaluated data are outside the fixed bilinear theorem.
5. **Every RSA method.** The common-geometry claim includes fixed linear RSA
   with a declared full-rank design. It does not include every similarity
   measure, model-fitting objective, noise ceiling, condition-generalization
   procedure, or inference method associated with RSA.
6. **Transport learning.** Crossform consumes a realized, identified transport.
   It does not learn anatomical registration, hyperalignment, functional
   correspondence, or other subject-to-group mappings.
7. **Universal downstream inference.** Implemented conditional HC3 and wild
   bootstrap procedures do not imply calibrated inference for every design,
   multiplicity control, cross-node covariance, random fields, arbitrary
   missingness, transport uncertainty, or causal interpretation.

These are scope boundaries, not judgments that the omitted methods are
unimportant. An upstream or extension package may supply them, but it must
preserve their identities and training/evaluation provenance.

## 2. Evidence classes

The CSV uses exactly these eight classes.

| Evidence class | What must exist | What the class permits |
|---|---|---|
| `algebraic_theorem` | A derivation with explicit assumptions and boundary conditions. | Claim the identity under those assumptions. It says nothing by itself about code or data. |
| `internal_oracle` | An implementation independent of the production lowering, mutation-sensitive tests, and declared numerical tolerances. | Claim computational agreement in the tested regimes. |
| `external_parity` | A named external implementation and version, an exactly mapped estimand, pinned inputs, tolerances, results, and source/artifact drift detection. | Claim parity only for the mapped case. |
| `matched_simulation` | A declared generator, known target, reproducible seed/replications, comparison methods, and pass criteria. | Claim recovery or operating behavior in the simulated regimes. |
| `existing_illustration` | A versioned real dataset and executable retrospective analysis with provenance and limitations. | Show that the workflow runs and expose descriptive diagnostics. It cannot confirm a hypothesis chosen after seeing the data. |
| `prospective_protocol` | A dated, hashed protocol frozen before eligible outcomes are inspected, including data eligibility, exclusions, estimands, analysis, thresholds, and failure handling. | Describe what will be tested. It cannot be called demonstrated, completed, or replicated. |
| `completed_real_data_result` | An eligible prospective protocol plus immutable inputs, execution receipts, results, and a deviation log. | Report the protocol's real-data result at the prespecified scope. |
| `independent_replication` | A completed result repeated through an eligible external dataset or analyst route that did not construct the original claim, with its own frozen protocol and public artifacts. | Claim replication of the prespecified target. Same-team reruns and retrospective examples do not qualify. |

The classes are not a single ladder. A theorem cannot be “promoted” into a
simulation, and external parity does not become scientific validation. Each
row answers a different question.

## 3. Promotion and refusal rules

- A source implementation earns `internal_oracle` only when the oracle does
  not call the production lowering and a deliberate mutation makes the court
  fail.
- A comparison earns `external_parity` only after versions, estimand mapping,
  normalization, partitions, tolerance, and a source/artifact manifest are
  fixed. A new external version creates a new evidence row.
- A simulation earns `matched_simulation` only after the target, generator,
  regimes, seeds or replication rule, comparison methods, and pass criteria
  are recorded. Exploring a simulation and then choosing a favorable endpoint
  remains exploratory.
- An `existing_illustration` never changes class merely because it uses public
  data or is reproducible. It can motivate a later protocol.
- A `prospective_protocol` becomes a `completed_real_data_result` only if its
  hash predates access to eligible outcomes, the run satisfies its eligibility
  rules, and every deviation is reported. If those facts cannot be proved, the
  result remains an illustration.
- A completed result becomes an `independent_replication` only through the
  external route defined above. Reanalysis of the same subjects, a new seed,
  another Crossform contributor, or another method on the same exploratory
  result does not qualify.
- A prospective row has `completed = FALSE` and cannot use “demonstrated,”
  “validated,” “confirmed,” or “replicated” as its status. The test court
  enforces this rule.
- Claims about intervals, prevalence, or power require matched-simulation
  calibration for their declared regime. Formula parity alone licenses
  “implemented,” not “calibrated.”

Claim `CF-H10` is a completed `matched_simulation`. Its 48-seed paired court
covers three planted organizations, three noise regimes, two sample sizes,
three SNRs, and four scales. All prespecified recovery, false-separation,
conservation, ambiguity, and negative-control gates pass in
`inst/extdata/certification/matched-interpretability-metrics.csv`. The
finite-sample coherence share is the signed crossvalidated coherent component
divided by the signed crossvalidated total; no positivity screen is applied.
This row licenses interpretation of the planted regimes only. It is not
real-data evidence and does not calibrate population intervals, prevalence,
power, or empirical neuroscience claims.

## 4. Current real-data boundary

The Haxby 2001 artifacts are `existing_illustration` and `external_parity`
evidence. They show that the matched point estimator, decomposition, adapter,
and refusal paths run on one public subject (and that the region-level
population slice runs on its recorded subjects). They are retrospective
regression fixtures. They do not establish a face/house neuroscience result,
population calibration, or independent replication.

The OpenNeuro ds003745 population-slice2 artifacts are also
`existing_illustration`. They provide valuable execution, transport,
commutation, and diagnostic evidence on public data. The analysis was not a
frozen prospective test. Its transport-efficiency ratio failed its own
interpretation audit, so the artifacts support no claim that functional
transport outperforms anatomical transport.

There are currently no `completed_real_data_result` or
`independent_replication` rows. Adding either class requires a new ledger row,
the prerequisite artifacts in Section 2, and review of every public document
listed for the claim.

## 5. Headline-document coverage

The `headline_documents` column maps each claim to the public pages that use
it. The executable ledger test requires coverage for these headline sources:

- `README.md`;
- `vignettes/introduction.Rmd`;
- `vignettes/from-observations.Rmd`;
- `vignettes/evidence-pairing.Rmd`;
- `vignettes/conservative-frames.Rmd`;
- `vignettes/population-form.Rmd`;
- `vignettes/common-geometry-equivalence.Rmd`;
- `vignettes/novelty.Rmd`;
- `design/unification-contract.md`; and
- `design/population-estimand-contract.md`.

This is coverage of headline scientific claims, not an assertion that every
sentence in those documents is a separate hypothesis. A new headline claim
must receive a claim ID, at least one evidence row, a boundary, and its public
document mapping before publication.

## 6. Ownership and review trigger

Every claim has exactly one accountable owner in the companion registry. The
owner is responsible for keeping its strongest artifact, limitation, current
status, headline-document mapping, and promotion history consistent. Any
change to a certification threshold or headline claim triggers review of the
ledger, registry, promotion history, public documents, and affected checksum
manifests. Protocol and external-parity rows keep their explicit boundaries
during that review; neither class is silently promoted.
