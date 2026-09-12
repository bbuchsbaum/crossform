# Agent-operated scientific state: crossform and rMVPA

Contract: `scientific-state-v1`  
Status: normative design; the inspection protocol and new adapters below are **not implemented by this documentation revision**.  
Reviewed against: crossform `b7cffa9` (`elite-pass`), rMVPA `7c93f8d` (`master`).  
Date: 2026-09-12

## 1. The system, not a collection of functions

An agent should operate a **typed, inspectable graph of scientific state**. It should not have to reconstruct the meaning of a result by reading several thousand lines of source, or guess whether a completed computation licenses a scientific claim.

The graph connects a question to identified observations, a declared estimand, a frozen learning procedure, numerical artifacts, independent evidence, and bounded claims. Existing R objects remain the computational objects. The graph is a small manifest and receipt layer over them, not a replacement workflow engine, database, autonomous scientist, or third modeling package.

The operating loop is: understand the current state; identify the cheapest eligible action that resolves the current uncertainty; preview its scientific and resource consequences; execute that action; inspect its receipt; retain the new evidence and limitations. A failed experiment can be a successful action when it cheaply rules out an approach. A faster computation is not a success when it changes the scientific question without permission.

There are two users of this interface: an analysis agent choosing scientific operations and an implementation agent changing the software. Both need identities, bounded scope, executable checks, and an accurate account of what remains unknown.

## 2. A tower of linked abstractions

| Level | Authoritative object | Question it answers | Must not be inferred from |
|---|---|---|---|
| 0 | Source/domain manifest | Which observations, features, units and anatomical coordinates exist? | Filenames, equal matrix shapes or positional train/test prefixes |
| 1 | Question/estimand declaration | What quantity, population, baseline and generalization axis are requested? | A convenient metric or default method name |
| 2 | Analysis/selection plan | Which observations may train, tune, calibrate and evaluate which operators? | Fold count or an assertion that a fit is frozen |
| 3 | Fitted relation and operators | What was learned, in which coordinates, with which numerical diagnostics? | Unmatched component column numbers |
| 4 | Prediction and measurement | What is predicted or measured, using which spatial access and metric? | An undifferentiated importance map |
| 5 | Evidence and uncertainty | What reproduced, under which null/error/multiplicity assumptions? | Sparsity, PSD projection, selected rank or algorithm completion |
| 6 | Claim and reusable knowledge | What assertion is supported within what boundary, and what next action is justified? | A high score, an old benchmark or another agent's confidence |

These are dependency levels, not automatic promotions. A fit may support descriptive prediction while failing an inference gate. A theorem, oracle, simulation and real-data study are different evidence types, not rungs on one confidence ladder. Preserve the vocabulary and claim ownership in `evidence-status-ledger.md`.

rMVPA owns observation-level fitting, CV, forward patterns, decoding, restricted predictors and its confirmation/loading-level procedures. Crossform owns effect relations, geometry/measurement queries, supported geometry scoring and population-form procedures. The interoperability contract translates objects and evidence obligations; neither package impersonates the other's private classes.

## 3. One meaning per quantity

Every readout carries a `quantity_id`, units, axes, baseline, conditioning set, aggregation and uncertainty scope. The initial controlled vocabulary distinguishes:

- Forward loading, calibrated-score Haufe pattern, decoding weight, signal standard deviation, and Gaussian model conditional information.
- Whole-brain predictive loss, local-restricted loss, local-adapted loss, independently fitted ROI loss, and geometry predictive gain.
- Fitted PSD geometry, signed cross-partition evidence, descriptive PSD projection, confirmed loading, selected rank, subspace stability and supported-rank evidence.

`rank_mean` is the current rMVPA metric averaging fold-selected ranks. It is neither an integer model rank nor significant rank. Keep existing public names; add explicit meaning rather than cosmetic renaming.

Observation weights, target/feature-set weights, neural metrics, spatial-frame weights, partition-pair weights and population weights are distinct typed roles. Transformations must record their order. In the inspected rMVPA implementation training/tuning can use observation weights while reported prediction metrics remain unweighted, and confirmation rejects nonuniform weights. An agent must see that difference before comparing outputs.

A numerical zero, a screened feature, absent coverage, rank deficiency, an unsupported request and a failed computation are distinct statuses. Do not represent all of them as zero or an unexplained NA.

## 4. Identity is the basis for trust and reuse

Use separate identities for scientific meaning, estimator specification, actual input revisions, implementation/environment and execution. The scientific ID includes effect/feature identities, units, target transformation, neural metric and composition, frame normalization, generalization axis, baseline and population target. An estimator ID adds model family, penalties, selection rules, tolerances and numerical conventions. Input ancestry includes upstream training dependencies, not only arrays consumed by the final call.

An artifact content digest establishes byte identity, not scientific validity. Matching scientific IDs establishes comparability of questions, not permission to reuse a fitted object. Cache reuse additionally requires identical eligible training inputs, transformations, estimator and relevant implementation. Hashes must be calculated from canonicalized contents; never trust a caller-supplied hash without binding it to those contents.

Every artifact declares parents, source revisions, feature and target order, origin ancestry, numerical status and available readouts. Orthogonal/oblique coordinate changes receive representation IDs while retaining an equivalence relation only after a tested transport proves the represented operator unchanged. Comparisons match explicit identifiers, never just row counts.

Protocol upgrades negotiate versions and capabilities. Unknown major versions refuse. A migration creates a new artifact with an explicit source link; it does not relabel old bytes as current. Scientific identity must not change merely because a tile size changes, but finite-tolerance execution and RNG details remain recorded in the execution identity.

## 5. Eligibility and outcome exposure

Each operator has a dependency set and role: anatomical fixed, externally fixed, training-learned, calibration-learned, or evaluation-derived. Composing operators unions their dependencies. Source aliases, renamed runs, cached projections, normalization statistics, nuisance models and target transforms retain ancestry.

Before confirmation, freeze the entire selection procedure, endpoints, exclusions, null, correction family, stopping rule and budget. Numerical hashes alone do not establish that a protocol preceded outcome access. Record the actual outcome-exposure history, including visual inspection, intermediate metrics and agent decisions.

A confirmation candidate becomes exposed as soon as an outcome-dependent result is read for choosing the analysis. That exposure propagates through derived artifacts. The same data can remain useful for exploration; it cannot regain untouched-confirmation status by changing a seed, renaming a file, rerunning a model or opening a new agent session. Exploring many ROIs or components on evaluation data is also selection.

A cross-partition product requires the stated error independence; a frozen predictive score additionally requires independence of evaluation from training/selection. Do not convert overlapping CV training estimates into independent runs. Ancestry checks can detect contradictions but cannot verify physical independence or correct temporal whitening. User assertions remain explicitly asserted, not certified.

## 6. The agent's compact state view

Supply one versioned machine-readable view of existing objects, mirrored by concise human output. The first response should fit a declared context budget and answer: what is this; what is known; what is missing; what may be claimed; what can be done next; what will it cost?

Required sections are identity, availability/status, scientific meaning, dependencies, evidence references, capabilities/refusals, diagnostics and bounded next actions. The accompanying `agent-protocol.schema.json` defines a minimum interchange shape, not a new public R API. Implement package-local adapters over existing print, validation, plan, result and receipt functions; do not add seven parallel workflow runners.

Progressive disclosure has three levels: a compact brief; a task-specific slice such as `locality` or `inference`; and exact source/artifact references for audit. Arrays, full histories and unbounded logs are not in the brief. Read-only inspection must not fit, materialize a geometry, scan all images or launch workers. Metadata uncertainty is printed as unknown rather than resolved through hidden expensive work.

Public examples must distinguish existing R calls from proposed protocol operations. Exported-symbol/signature checks, examples and method routing should be generated from one capability registry once implemented. Static guidance is a routing index, not an alternative source of runtime truth.

## 7. Actions are bounded scientific transitions

An action descriptor contains the existing executable entry point, arguments or an argument-schema reference, input artifact IDs, required capabilities, declared reads/writes, resource envelope, expected outputs and postconditions. It also states whether it fits, retunes, exposes outcomes, changes the estimand, mutates shared state or requires human authorization. A runnable descriptor cannot point to an unimplemented function.

Separate proposed, runnable, refused and completed actions. A refusal gives a stable reason code, the exact failing condition, relevant artifact IDs, and minimal remedies with costs and scientific consequences. Remedies are data, not automatically executed commands. Never repair an inference refusal by silently downgrading to a different null, global shuffle, descriptive map or different target.

A plan diff separates scientific changes from implementation-only changes: changing ROI access, neural metric, response weights, correction family or the training baseline changes meaning; changing chunk size need not. Rank and penalty searches are training actions even when cheap.

Side-effect-free preview is the default. Writes to source repositories, shared libraries, remote compute and persistent data require existing authorization. Never install into a user's R library, rewrite source data or modify a benchmark threshold merely to make an action pass.

## 8. Economical execution without hidden scientific shortcuts

Budget the full dependency graph: data reads, copies, target compression, covariance pilot, outer/inner fits, candidate path, retained folds, regional queries, resampling and serialization. Report estimated peak memory as a range with assumptions, separately from measured RSS; an R workspace estimate is not a hard process-memory guarantee.

Use a cost ladder: metadata validation; tiny independent oracle; representative training-only pilot; scoped fit; full evaluation. Scaling probes must not consume reserved confirmation outcomes. Rank-zero/baseline and unpenalized comparisons are explicit useful actions, not embarrassing special cases. Do not assume spatial sparsity improves prediction.

Cache only when profiling and repeated reuse justify it. Current rMVPA does not need a speculative global fit cache to remain usable. Cache anatomical geometry broadly; cache learned statistics under exact training/transform keys. Label-permuted fits cannot reuse supervised intermediates unless their transformation equivalence is proved. Reuse does not mean averaging independently rotated solutions.

Warm starts are within compatible penalty paths. Cross-rank warm starts are not enabled merely because a lower-rank artifact exists; the historical rMVPA implementation removed them after convergence problems. Record any approximate reuse and its tolerance separately from exact reuse.

Checkpoint atomic completed units with plan/input/implementation IDs, RNG stream, dimensions, solver phase and numerical diagnostics. Resume must validate the checkpoint and preserve selection bookkeeping. A changed plan forks a new execution, never edits a completed receipt. Concurrent workers/agents use leases or compare-and-swap and immutable outputs; incompatible results are not merged by last-writer-wins.

Parallelism has one owned level unless a measured nested policy is explicitly admitted. Stop on budget exhaustion with retained valid outputs, unavailable cells and a resumable state; do not report partial evaluation as a complete benchmark.

## 9. Evidence that accumulates rather than fossilizes

Retain small experience records: problem fingerprint, attempted plan, source/implementation IDs, observed behavior, uncertainty, failed assumption, alternative tested, test/receipt references, applicability boundary, invalidation triggers and review owner. These are not conversation transcripts or executable instructions from untrusted artifacts.

Before proposing a new action, retrieve the nearest applicable records and explain matches and mismatches. A benchmark is not portable merely because p and n match: signal organization, noise, blocks, target rank, missingness, BLAS and budget matter. Conflicting records coexist with their scopes. New evidence supersedes a recommendation, not the historical observations.

A useful discovery becomes one of: a regression test, a scoped benchmark row, a refusal, a recipe with preconditions, or a documented unresolved question. Promotion to a default requires reviewed evidence across declared regimes and a versioned decision. No self-modifying numerical policy or silent automatic threshold relaxation.

Crossform's existing evidence classes/claim registry remain authoritative. Agent usability evidence is separately identified; fewer tool calls cannot license stronger neuroscience claims. Historical documentation reports are `recorded_unverified` until their artifacts are checked against the relevant source, not automatically recertified by this revision.

## 10. Interpretation and loss accounting

Keep the low-rank relation and its coordinate transforms primary. A rotated loading display is a view, not a newly fitted model or a new confirmatory basis. Preserve coupling terms for oblique factors. A sum of component predictions has cross-terms in its squared norm; additive spectral-mode gain must not be copied onto arbitrary rotated components.

For a frozen predicted form F, geometry gain is 2<F,G_test> - ||F||^2. Its null expectation for nonzero F is negative, not zero. Coherent/configuration inner products may add before subtracting the prediction cost once. A geometry score is not classification accuracy, explained variance, or an automatic rank test.

A signal-SD map summarizes task covariance in original measurement units. Conditional-information maps describe the working Gaussian model, not causal necessity or an empirical categorical-information estimator. Under x independent of y given deterministic t=C'y, information about x carried through y equals that through t; a Gaussian surrogate for categorical labels does not establish that identity for empirical class information.

## 11. Locality and population are contracts

Distinguish a graph regularizer from a spatial measurement frame. Distinguish restriction of an existing precision metric from inversion of a restricted covariance. Local-restricted prediction uses (Psi_RR)^-1, not (Psi^-1)_RR; test preprocessing must not read outside-region measurements. Local-adapted covariance and independently fitted ROI models get different identities.

Group geometry and signed group loadings are different operations. Current rMVPA loading pooling requires exactly compatible target subspaces and declared one-to-one feature mappings; current crossform population forms have their own transport, coverage and metric gates. Neither pathway supplies missing cross-feature covariance for arbitrary interpolation. New adapters cannot bypass these restrictions by coercing classes.

## 12. Agent acceptance is a testable product requirement

An agent starting without chat history must identify the checked-out state, find the owning contract, select the correct existing entry point, diagnose a refusal, estimate cost, preserve untouched confirmation, and hand off a minimal reproducible state. Evaluate this with the fixed tasks and budgets in `agent-evaluation-plan.md`, against current documentation as a baseline.

Safety/correctness violations have zero tolerance on the defined fixtures. Resource and context improvements are measured, not promised. Record runtime, bytes, tool calls, context, unnecessary reads, recomputation, corrections and task success independently. A task suite does not establish universal agent competence.

## 13. Authority and rollout

Existing mathematical/numerical contracts continue to govern their estimands; this contract adds the control, identity and handoff layer. Conflicts are diagnosed and resolved through a versioned decision, not silently won by whichever file an agent read last. The active rMVPA plan owns its implementation tasks; `pattern-model-interop-plan.md` owns cross-package tasks. Historical plans and receipts are archived intact and clearly labeled.

Start with metadata-only inspection and accurate source/status indexing. Next add frozen-prediction interchange and refusal conformance. Add persistent execution caching/checkpoints only where measured workloads justify them. Keep local precision, TV support envelopes, new sequential rank inference, oblique inference and joint hierarchical fitting outside this revision's implementation commitment.

No model code, runtime adapter, new inferential capability or scientific certification is delivered by these documents. The design enables subsequent agents to implement those changes without having to rediscover the system.
