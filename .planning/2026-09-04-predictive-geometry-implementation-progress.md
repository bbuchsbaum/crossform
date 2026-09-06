# Predictive geometry implementation progress

Objective: implement all G01–G22 and T01–T62 in the work-package/test plans.
The source of live work-item IDs, dependency state and per-test evidence is
`2026-09-04-predictive-geometry-implementation.json` in this directory.

Started 2026-09-04, actor `codex-predictive-geometry-implementation`.
Mote epic `bd-01M1QMFAQJB9BVJZGKXR72FRYT` has all 22 bounded children and
their declared dependency edges. No prior implementation was discarded.

The original dirty source/tests/docs/evidence are preserved in
`/private/tmp/crossform-predictive-baseline-dz2i2luh/baseline.tar.gz`.
Initial R digest: `sha256:24afa8ce392d71efa930c2dcb0afd9ac0951d100a424965c5a89a949355f1a8d`.

G01 migration decisions are recorded in `design/predictive-geometry-migration.md`.
The full goal remains active until feature, validation, mutation and release
evidence have passed requirement-by-requirement scrutiny.

G01–G04 complete. The v1 contract freezes PSD admission, post-truncation
normalization, positive-support restriction, tie policy and independent gain.
The base-R oracle passes its 14 fixtures, 16 wrong alternatives, and dense
partition/minimum-norm/completed-square checks. Explicit kernel input now
passes 41 new admission assertions and 246 existing model-basis assertions.
G05 is active. Broader T01–T62 obligations remain pending until their full
production paths exist; these counts are not release certification.

G05 complete: 95 predictive-kernel and 246 existing model-basis assertions
pass. Model rank-input fixture now has a resolved spectral gap; T11 owns
the changed positive-tie refusal. Normalization uses retained trace and
records numerical roots separately. G06 now admits overlap and pools kernels.

G06 complete. Overlap is admitted in the basis and pooled form; the old
coefficient reader refuses it before execution. Duplicate/split-weight,
simplex-vertex and near-overlap checks pass. Five targeted files pass, with
one existing simulation skipped under the on-CRAN test policy. G09 proceeds
from its satisfied dependencies while provenance G07/G08 remain pending.

G09 complete. The pure estimator restricts to positive pooled support,
shifts by the inverse kernel spectrum and reuses one shifted decomposition
for the rank path. Analytic/noncommuting/singular/zero/tie, minimum-norm
factor optimization and controlled-conditioning tests pass. Existing latent
tests pass with two designed on-CRAN skips; architecture passes. Only the
private rank-PSD helper now admits zero; public latent admission is unchanged.
G07 is active: actual observation ancestry and fitting dependencies.

G07 complete. relation(provenance = list(observation_origins = ...))
canonicalizes a declared parent manifest and attaches ancestry to sources.
Subsets/aliases preserve it. The higher-layer collector records positive-
weight endpoints, source revisions, extractors and upstream fitting supports.
44 ancestry checks, relation/lowering regression tests and architecture pass.
The graph checker initially counted local names as calls; descriptive local
names remove those false edges without weakening the acyclic-layer rule.
G08 is active: metadata admission and reduced execution.

G08 complete. Reduced preparation excludes zero/unused endpoints, preserves
the fixed target and uses existing materialize_geometry() without the old
trace bank. 93 dense/read/admission assertions pass for regional/searchlight
frames, identity/dense metrics and non-square extractors. Explicit effect-space
meaning is now bound in the model declaration. The executor currently needs
undirected symmetric full forms and SPD metrics; directed, singular, learned
and whitened routes have explicit preflight refusals. These executor limits
are recorded in the v1 contract. A legacy hash-mutation test could replace
0 with 0; it now guarantees a changed byte. Regression and architecture
checks pass. Next: G10 fitted records and G11/G12 independent frozen scoring.

G10 completed: frozen factor/spectrum storage and public fit_geometry; predictive-fit, API-surface and architecture tests pass. Temporary load_source roxygen generation misclassified base print/format generics; restored their S3 registrations, confirmed export count 122. Additive-frame fixture uses the documented members route for stable IDs. G11 started: target and independence preflight.

G11 complete: target match and origin preflight, 22 assertions plus regression courts. G12 started: frozen gain and signed mode/group evaluation.

G12 complete: frozen independent scoring, mode/group evidence, coherent/configuration single-offset readout, 63 assertions plus regressions. Numeric payload and semantic-byte immutability are checked separately from R JIT changes to serialized closure bytecode. G13 started: memory admission and failure cleanup.

G13 complete: reserve predictive memory before compiler admission, bounded row readers, owned failure cleanup; 44 assertions plus regressions. The existing fixed-metric complete executor has no block route; explicit typed refusal added. G14 begins the executable public example.

G14 complete: public fixed-split example runs under a fresh isolated package install; 10 exact assertions. The fixture caught whole_brain default local normalization; example explicitly declares none to match the analytic feature-sum geometry and penalty units. G15 begins prespecified statistical calibration.

G15 complete: 87,000 independent production datasets, all 12 conditional-risk bias/precision gates and four generator moments pass; null leakage control detected. Pilot and production seeds/counts are frozen and separate. Paired model/isotropic/mismatched comparisons reported without universal superiority claims. R source remained unchanged throughout the run; later integration needs G22 rebinding. G16 is active.

G16 complete: 52 assertions, manual weighted edge folds, retained fitted predictions, separate norm costs and block round trips. Cross-fit records describe procedure risk at fold training sizes, not an all-data refit. G17 begins nested selection and its complete-pipeline calibration.

G17 deterministic work: 54 nested-selection assertions and architecture pass. Both scopes, exact ties/candidate order, manual weighted inner fits, outer-test perturbation, full stochastic-generator/public-route parity and serialization/mutation checks pass. Adaptive calibration pilot has started; G17 remains doing until pilot-fixed production calibration passes.

G17 complete: 54 deterministic assertions plus architecture; full adaptive calibration passes on 2000 null / 5000 signal datasets with independently piloted frozen counts. Null normalized identity error -0.0000943; signal -0.00287449, both five-MCSE and two-percent equivalence gates pass. Invalid outer-test selection has detectable optimism. G18 started: isolated training caches with explicit budget and invalidation.

G18 design in progress: keep separate training and inner-validation cache objects; no global adaptive cache. Keys must cover actual support/revisions/upstream transforms, ordered target (including measurement order), original compute policy, model declaration, pooled weights, penalty and signed-row content for spectral reuse. A context change clears neural forms/paths; data-independent pools may survive. Candidate family is declared up front so cache memory has a conservative bound. Reserve every retained reduced form/eigensystem and concurrent partner cache before compiler admission. Use a full-union rank buffer bound for cached runs so rank candidates share one geometry execution budget/receipt. Score-cache admission similarly bounds frozen input by full-union factor width. Reuse fit/score executors via optional private cache arguments; public verbs remain unchanged. The cache file must never call fit/score/selection files, avoiding file-level cycles. Nested selection creates one training and one validation cache per inner fold plus a final training cache; return counters separately from scientific selection identity. G18 tests must count source materialization and neural spectral calls, prove parity and invalidation for source/fold/effect/component/frame/order/metric/model/normalization changes, and reject undeclared cache candidates rather than silently growing beyond admission.
G19 follow-up found during review: as.data.frame(score, view="modes" or "groups", rows=integer()) currently needs an explicit empty-result schema or typed refusal; default summary already supports empty rows. Add regression while completing views/docs. Final G22 must rerun both new calibration runners against the final source; their pilot counts are frozen and should not be tuned to production outcomes.

G18 complete: 77 cache assertions plus architecture pass; all targeted source/fold/effect/extractor/upstream/frame/order/metric/model/normalization/compute changes invalidate correctly. Rank candidates reuse forms and spectral paths; cached and uncached predictions, scores and selection identities agree. Separate evaluation caches never fit modes. Conservative memory and undeclared-recipe admission precede neural reads. G19 starts final public workflow, views and migration documentation.

G19 complete: runnable predictive-geometry vignette, fitted/scored examples, same-condition scope, model spectral identity and migration/navigation align. Score empty/repeated/reordered views have stable schemas; 80 score, 47 fit, 13 example, 26 API, 12 architecture assertions pass. Rendered guide math/output/links/reference inventory pass at /private/tmp/crossform-g19-site-1604c78786d15. Site build exposed seven existing exported population topics missing from index; added them to the existing section. R cache redirected to /private/tmp for sandboxed site build; no published site changes. G20 begins matched tiny/medium/large runtime/memory/read/parity evidence.

G20 complete: 39 focused assertions plus all seven serial benchmark routes pass. Artifacts benchmark-results/predictive-geometry-performance bind source 10bfc4506bcf9c110f8ac10a0bfba425fb0f3eb4c5975e608bd1cd095690bf12 and all four harness files. Large blocked workload: 65.949s, 145792120 planned bytes, 1335558144 incremental RSS bytes; maximum dense error 3.24e-14. Report explicitly distinguishes logical workspace admission from OS RSS: 512 MiB is not a hard process cap. Tiny pilot full baseline initially read unused parent partitions; subset original relation to identical endpoints before comparison, after which reads match. ps_loadavg supplies observed host load despite sysctl shell refusal. Existing descriptive trace/unpenalized form parity passes. Medium six-candidate cache: 12->2 materializations, 12000->4000 spectral paths; 10.718s->9.479s warm instrumented time, unchanged gains. G21 begins real production-function mutations in disposable test scopes.

G21 complete: all 20 real production-closure mutants killed by named primary T-tests. Baseline and no-op instrumented controls pass; every replaced function was exercised. New seven-assertion witness file proves public signed-before-PSD compression (1 vs wrong 1.5) and invalidation when only edge weights change at fixed endpoints/shapes. Source R unchanged at 10bfc4506bcf9c110f8ac10a0bfba425fb0f3eb4c5975e608bd1cd095690bf12; performance artifacts remain current. Mutation bodies, hit counts, semantic failures, baseline/control records and source/test/harness bindings retained. G22 begins full integration, T01-T62 evidence audit, R CMD check and final source-bound recertification.

G22 complete. Full source suite: 11851 pass, 0 fail/error/warning, 7 documented optional skips. Final source-tooling additions: 3 oracle-binding and 15 frozen-count restoration assertions pass; focused certification/admission/parity courts pass. All 62 requirements audited and linked in design/predictive-geometry-test-evidence.md. Nine legacy runners, both fixed-count predictive calibrations (87000 fixed-parameter + 7000 selected datasets) and pinned R/Python parity pass on unchanged R digest 10bfc4506bcf9c110f8ac10a0bfba425fb0f3eb4c5975e608bd1cd095690bf12. All 20 production mutants killed with exercised positive controls. Compact 27476-byte receipt and 11 source-bound shipped artifacts pass binding; only the explicitly excluded shard artifact remains unbound. Promoted CSV summaries match artifact fields. Initial R CMD check 0/1/1 exposed a Homebrew compiler versus R-header warning; an isolated Apple-compiler override first omitted C++17, then was corrected. Final R CMD check 0 errors, 0 warnings, 1 explained incoming development/URL note; 9538 installed assertions pass, 119 designed source-only/optional skips. No global compiler setting or production R source changed during final verification. Receipt CSV LF correction was regenerated through the complete external pipeline. Final report: design/predictive-geometry-certification.md. No commit or publication; pre-existing measurement/coupling defect and missing-cache Haxby receipt remain separately labelled.
