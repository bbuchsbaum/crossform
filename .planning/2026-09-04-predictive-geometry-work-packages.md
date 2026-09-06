# Predictive geometry: granular implementation plan

Date: 2026-09-04. Status: proposed implementation program, not implemented
functionality. Planning task: `bd-01M1QFP2456NYHAT6WYRATZZ5Z`.

**Build one regularized representational-form estimator and one independent
score, using crossform's existing effect-space algebra.** A model specifies
preferred directions through its kernel, the fit learns their neural
amplitudes, and the frozen score measures predictive improvement.

This is an independent design plan. The other coder's M2–M9 work is a useful
partially integrated baseline, not an instruction to continue its API or
ticket sequence. Its tests certify its earlier contracts. Preserve that
evidence and the dirty worktree; reuse components where they satisfy this
objective. Do not commit, delete or rewrite that work as part of planning.

This document supersedes the broad P1–P7 sequence and provisional migration
choice in the [earlier review](2026-09-04-predictive-geometry-plan.md).
The [test specification](2026-09-04-predictive-geometry-test-specification.md)
defines T01–T62. Each package below has a bounded output, dependency list,
test obligations and a completion gate. IDs G01–G22 are plan identifiers;
they have not been created as 22 live Mote issues.

**Architecture to retain, and decisions to change**

- Reuse centered model factor construction and extractor-based lowering
  after independent tests establish their contracts. Reuse existing geometry
  execution, packed codecs and blocked storage.
- Preserve signed source geometry. Fit
  `A = U [U' S U - lambda D^-1]_(+,s) U'` on the positive support of
  `J_alpha = sum(alpha_i R_i R_i')`. Score
  `2 <A,S_test> - ||A||_F^2`. These are different operations and objects.
- Permit overlapping models for invariant form prediction. Do not require
  a unique coefficient for every redundant feature.
- Keep one data-independent model declaration and two verbs as the target
  surface: provisionally `model_basis()`, `fit_geometry()`, `score_geometry()`.
  The constructor name is reusable because it already carries model factors
  and spectra. The old four-structure reader is not part of the minimal
  predictive API; its final disposition is an explicit G01 decision.
- Treat positive-penalty full-span models as regularized predictors. Retain
  a labelled unregularized full-span baseline. Compare model preference
  against an isotropic kernel on the same span.
- First support a fixed neural metric, fixed conditions/measurement target,
  and independent partition subsets of one declared parent relation. Hold
  out condition transfer, learned metric schedules and GLS until their own
  estimands and validation contracts exist.

Suggested user workflow (new arguments/verbs are proposals):

```r
basis <- model_basis(
  kernels = list(visual = K_visual, semantic = K_semantic),
  conditions = relation$effect_space,
  rank = c(visual = 8L, semantic = 5L),
  normalize = "trace"
)
fit <- fit_geometry(train_plan, basis,
  weights = c(visual = 0.5, semantic = 0.5),
  rank = 4L, penalty = chosen_penalty)
evidence <- score_geometry(fit, test_plan)
```

The rank, weights and penalty in this example are declared or chosen within
training. Scoring cannot adjust them. The two plans state neural partition
generalization; the fit record separately states what trained the model.

**Dependency map**

| Work package | Depends on | Main output |
|---|---|---|
| G01 | — | Scope, reuse and migration decision |
| G02 | G01 | Mathematical, numerical and validation contract |
| G03 | G02 | Independent oracle and exact fixtures |
| G04 | G02 | Explicit kernel admission |
| G05 | G04 | Truncation and normalization |
| G06 | G05 | Overlapping factors and pooled support |
| G07 | G02 | Actual training-support records |
| G08 | G06, G07 | Reduced-plan preparation |
| G09 | G03, G06 | Pure spectral estimator and rank path |
| G10 | G08, G09 | Fitted prediction object |
| G11 | G07, G08 | Evaluation admission and target matching |
| G12 | G03, G10, G11 | Frozen predictive score |
| G13 | G10, G12 | Bounded storage and failure handling |
| G14 | G13 | First complete fixed-split example |
| G15 | G14 | Statistical calibration |
| G16 | G14 | Declared cross-fit schedule |
| G17 | G15, G16 | Inner selection and outer evaluation |
| G18 | G16, G17 | Training-local reuse and cache identity |
| G19 | G01, G14, G17 | Final API, explanatory views and migration |
| G20 | G13, G18 | Matched performance evidence |
| G21 | G03, G14, G17 | Mutation-based test adequacy |
| G22 | G15, G19, G20, G21 | Integrated verification and certification |

G03, G04 and G07 can progress independently after G02. This is a dependency
observation, not an instruction to start agents or edit overlapping paths.
The first complete internal vertical slice is G14. G22 is the completion
gate for the whole planned program. A smaller release omitting tuning must
explicitly narrow its contract and apply G15/G19/G20/G21/G22 to that scope;
G14 alone is not release certification.

**G01 — Decide scope and migration from the existing partial work**

Dependencies: none.
Tests: T59, T60.

Inventory the uses of `model_basis`, `model_geometry`, latent rank projection,
the old result classes and their tests. Classify each as retain, generalize,
internalize, or replace, with its scientific question and caller migrations.
Use `design/api-tiers.md`, not the current export count, to judge the surface.
Keep ordinary RSA's fixed regression contract. Default to excluding the
four-structure taxonomy from new predictive examples; retain separate
nonnegative/additive hypotheses only where they earn a clear user workflow.
Output a migration table and a final naming decision. Do not let an old
contract assertion silently determine a new estimand.

Exit: every affected public behavior has a disposition and every retained
mathematical guarantee has a destination. The old working tree remains
recoverable; any later removal has an explicit migration test.

**G02 — Freeze the objective, admissible inputs and score interpretation**

Dependencies: G01.
Tests: T04, T07, T11, T16, T22, T30, T39, T55.

Write a versioned predictive-geometry contract. Pin normalization after
truncation, named simplex weights, rank zero, rank clamping, pooled-support
tolerance, positive rank-cut ties, units and target effect coordinates.
Recommend numerical-only PSD repair by default; material model projection
must be an explicit transformation with recorded mass. Declare that penalty
shifting can rotate modes, source spectra stay signed, and gain has squared
geometry units. Specify conditional unbiasedness and actual-origin support,
plus exact refusal reasons for unsupported cases. Freeze numerical-error and
simulation equivalence margins before production results exist.

Exit: the contract determines the outcome of every fixture and every
unsupported-input example without consulting an implementation. This step
owns specifications; downstream packages own the corresponding code assertions.

**G03 — Establish an independent mathematical court**

Dependencies: G02.
Tests: T12, T13, T14, T15, T19, T20, T33, T34, T38, T62.

Promote reviewed versions of the 14 exact fixtures into a base-R oracle,
with no package numerical helper imports. Add dense raw partition-product
construction, completed-square identities and minimum-norm feature-factor
checks. Keep analytic 2x2 solutions independent of the production eigensolver.
The earlier exploratory script remains supporting evidence, not the sole
oracle. Include deliberately wrong local variants and record which exact
assertions reject them.

Exit: the oracle runs alone and the fixture arithmetic is checked; package
comparison entry points are ready for G09/G12. Numeric optimization is labelled
a supplementary comparison rather than proof of global optimality.

**G04 — Admit explicitly typed model kernels**

Dependencies: G02.
Tests: T01, T02, T03, T04.

Extend or replace the relevant constructor portion in `R/model-basis.R`.
Support declared kernels, features and squared-Euclidean RDMs without guessing
matrix meaning. Align both kernel/RDM axes by name. Check shared effect units,
scales, centering convention, symmetry, finiteness and PSD policy. Record the
source representation and any authorized transformation. Validation precedes
neural reads. Preserve the explicit centered-subspace construction.

Exit: all three input routes produce the same kernel on an exact fixture;
invalid or ambiguous representations produce typed, actionable errors.

**G05 — Make model truncation and scale reproducible**

Dependencies: G04.
Tests: T05, T06, T07, T10, T11, T23.

Apply per-model rank selection in the centered subspace, then normalize the
retained kernel under the declared convention. Record original and retained
trace, requested/effective rank, numerical dropped roots and material dropped
mass separately. Implement tied-cut policy and zero-model refusal. Do not
silently select the truncation from neural responses. Stress the existing
poor-conditioning regression before reusing its tolerance machinery.

Exit: input rescaling and baseline shifts obey the stated laws; F10 and
rank/tie boundaries behave exactly as declared.

**G06 — Admit overlap and pool inside fixed union coordinates**

Dependencies: G05.
Tests: T08, T09, T18, T21, T23.

Allow rectangular `R` in the model declaration and rederive all affected
validators. Construct `J_alpha` from the per-model coordinate kernels and
factor only its positive support. Zero weights remove model support without
letting an unpenalized null space enter the fit. Move requirements for unique
old coefficients into their surviving readers. Record union dimension,
pooled dimension and conditioning separately.

Exit: duplicates with split weights give identical predictions; F02 refuses
the pseudoinverse-only shortcut; zero-penalty/full-span and positive-penalty
interpretations are correctly distinguished.

**G07 — Represent actual fitting dependencies**

Dependencies: G02.
Tests: T26, T27, T28, T29, T30, T32.

Define a sealed training-support record for actually used observation origins,
partition products, source revisions and upstream transforms/selections.
Carry origin ancestry through lowering and any partition subsetting; do not
derive it from display names or matrix content alone. For the first version,
admit only known disjoint subsets of a common declared origin manifest.
The manifest records caller/ingestion assertions; no code can prove hidden
independence or detect dishonest provenance from arbitrary numeric arrays.
Unknown cross-dataset origin is a refusal, not an independence flag bypass.

Likely seam: a layer-2 record/validator plus a higher-layer collector that
reads plans, avoiding a relation↔support dependency cycle. Reuse existing
ingestion identities where they actually carry observation ancestry.

Exit: copied and renamed known observations remain overlapping, while equal
numeric values from distinct declared independent origins are not conflated.

**G08 — Prepare the reduced geometry without a second execution design**

Dependencies: G06, G07.
Tests: T24, T25, T26, T31, T58.

Factor shared admission/lowering preparation out of the current reader where
useful. Compose the model pullback with extractors and use ordinary
`plan_geometry`/`materialize_geometry`. Preserve frame, fixed metric,
pairing, compute policy and lineage. Retain nonlinear/unsupported schedule
refusals. A predictive fit does not request the current reader's additional
centered/orthogonal trace bank. The shared helper must not call back into
either high-level reader.

Exit: dense raw algebra and reduced execution agree; instrumentation proves
that only the declared partition effects and reduced forms are read.

**G09 — Implement the pure regularized spectral fit**

Dependencies: G03, G06.
Tests: T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23.

Implement the small objective on admitted symmetric `S` and pooled support.
Reuse the rank-PSD numerical primitive if it passes the independent court;
extend its zero-rank handling deliberately. Store the shifted spectrum and
fitted amplitudes with distinct meanings. One eigendecomposition supports
all rank budgets at fixed kernel/penalty. No full-space projection, iterative
model rotation, diagonal/block optimizer or generic GLS belongs in this step.

Exit: analytic optima, numerical constraints, metamorphic laws and stress
residuals pass. No source mutation and no unsupported full residual claim.

**G10 — Create a frozen fitted prediction value**

Dependencies: G08, G09.
Tests: T16, T29, T39, T41, T42, T58.

Add `fit_geometry()` and a sealed `effect_geometry_fit`. Store reduced modes
and amplitudes, model/support specification, measurement IDs, fit diagnostics,
actual training dependencies and parent execution receipts. Keep the signed
training form as a block-store reference when retained. Separate prediction
semantics from factor orientation and execution identity. Do not inherit
raw-geometry capabilities. Avoid claiming rotation-invariant byte hashes of
arbitrary floating-point factors.

Exit: a fit can be serialized/reopened and scored without rereading training
outcomes; mutation of a derived or identity-bearing field is detected.

**G11 — Admit evaluation only for a matching, independent target**

Dependencies: G07, G08.
Tests: T27, T28, T29, T30, T31, T32, T40, T55.

Build the score preflight. Check actual train/test support and both declared
independence requirements; compare target effect coordinates, units, metric,
component and frame separately from pairing/plan identity. Establish a unique
measurement-ID mapping before reading evaluation values. Reject unsupported
transports, unknown origins and contradictory records with typed reasons.
Introspection and metadata integrity checks must not be treated as model
training on every source stored in a parent relation.

Exit: valid held-out plans admit; each wrong-target or overlap counterexample
fails before expensive evaluation work.

**G12 — Implement frozen gain and mode/group evidence**

Dependencies: G03, G10, G11.
Tests: T32, T33, T34, T35, T36, T37, T38, T39, T40, T41, T42.

Add `score_geometry()` using reduced test forms and existing isometric packed
coordinates. Return inner product, squared prediction norm and gain, plus
optional mode/group evidence. Preserve negative evidence/gain. Apply the
same frozen operator when exposing coherent/configuration evidence and count
its offset once. Rank or penalty is never selected here. Score records carry
both origins and a separate validation record, not invented standard errors.

Exit: F05–F09 give exact results; dense/pulled-back/packed routes agree;
reordering measurement rows does not misapply predictions.

**G13 — Bound memory and handle partial execution honestly**

Dependencies: G10, G12.
Tests: T25, T41, T43, T44, T57.

Consume reduced stored forms in bounded row blocks. Retain O(r*s) prediction
factors per measurement, not O(d*d) coefficients for redundant model features.
Account for input, temporary spectral buffers, output stores and worker count
in memory admission. Use existing blocked storage and cleanup mechanisms.
Read failures, interrupts and incomplete files cannot produce a complete fit
or score receipt. Cover all admitted routes; refuse unsupported ones.

Exit: guarded stores prove bounded reads, supported execution paths agree,
and deliberate partial failures cannot masquerade as completed predictions.

**G14 — Deliver the fixed-split vertical slice**

Dependencies: G13.
Tests: T01, T24, T31, T33, T38, T40, T42, T55, T60.

Create a small executable example from synthetic effect partitions with known
origins and an admitted spatial frame. Use one training pair and one disjoint
evaluation pair, fixed model weights/rank/penalty, and all three model-input
routes. Present fitted amplitudes, signed replicated evidence and gain. Show
the zero predictor and one harmful overestimated mode. Exercise the complete
public path and its refusal examples without external data or API credentials.

Exit: another developer can run the example from a clean package installation
and explain every displayed value using the fixture algebra.

**G15 — Calibrate the fixed-fit predictive claim statistically**

Dependencies: G14.
Tests: T38, T45, T46, T47, T54.

Implement the scheduled independent-data generator and dense known-signal
target. Include pure noise, planted low-rank signal, model mismatch,
heterogeneous partition noise and admitted fixed metrics. Record risk-identity
error across independent datasets, not individual dependent pairing edges.
Run planned leakage controls. Compare same-span isotropic and model kernels
with matched settings. Use prespecified precision/margins from the test spec.

Exit: measured bias and its uncertainty satisfy the declared contract, or
the report identifies insufficient precision/model-assumption failure.
Do not change thresholds or seeds to obtain a pass.

**G16 — Reuse the fitter/scorer over a declared cross-fit schedule**

Dependencies: G14.
Tests: T26, T30, T48, T49.

Extract the useful pairing-driven schedule from the old cross-fitted energy
path, keeping its estimand distinct. Each fold calls the same fit and frozen
score. Preserve admitted edge orientations, nonuniform weights and actual
training origins. Record fold-specific ranks, predictions and offsets.
Declare what the aggregate estimates: the training procedure at those fold
training sizes, not the performance of a later all-data refit.

Exit: explicit manual per-edge computations match; unsupported schedules
refuse; offsets aggregate as sums of squared norms, not the norm of a mean.

**G17 — Add inner selection with untouched outer evaluation**

Dependencies: G15, G16.
Tests: T32, T36, T50, T51, T52, T55.

Implement internal orchestration or a narrowly scoped example first, before
inventing another public model-selection framework. Accept declared candidate
tuples and a feasible nested split. Select using weighted inner validation
gain, resolve exact ties deterministically, refit on outer training data and
score the untouched outer test. Record selection scope: global or per
measurement. Include all model/preprocessing choices in the training graph.

Exit: changing outer-test outcomes cannot change any learned choice, the
eight-partition manual fixture agrees, and full-pipeline calibration passes.
Four partitions are not advertised as sufficient for arbitrary nested tuning.

**G18 — Reuse training products across compatible candidates**

Dependencies: G16, G17.
Tests: T18, T32, T49, T53.

Cache reduced forms within one declared training fold, pooled-kernel support
per weight specification, and shifted eigensystems per penalty. Rank paths
reuse spectra. Keys include origins/revisions, model truncations and
normalization, effect target, component, fixed metric, frame and fold.
Start by treating changed model-rank specifications as separate preparations;
only add maximal-basis reuse after a measured need. No global mutable cache
that can mix training and evaluation states.

Exit: counters prove reuse where valid, every invalidation case passes, and
cached/uncached predictions are numerically equivalent.

**G19 — Finish the user workflow and migrate the partial API**

Dependencies: G01, G14, G17.
Tests: T35, T36, T37, T39, T41, T55, T59, T60.

Implement concise printing/data-frame views, executable documentation and the
G01 migration map. Default summaries show measurement, gain and effective
rank; detailed views separate training amplitudes, evidence and cost. State
same-condition partition prediction versus condition transfer. Remove or
internalize provisional APIs if the scoped migration supports it; do not
retain them to preserve 121 exports. Conversely, do not pretend ordinary RSA
replaces nonnegative additive-kernel fitting when retiring a reader.

Exit: the documented workflow is executable, the public surface has one
clear job per verb, and every changed old test names its revised contract.

**G20 — Establish matched performance and memory evidence**

Dependencies: G13, G18.
Tests: T25, T43, T44, T53, T56, T57.

Implement the three benchmark sizes in the test spec. Compare the proposed
route with full-space dense reference where feasible and the current view
where its outputs overlap. Match metric, folds, source reads, requested
outputs, rank and numerical precision. Report incremental RSS, wall time,
stored bytes, neural read counts, kernel/eigen calls and parity. Run the heavy
cases serially; the user-reported earlier memory-pressure failure reinforces
this requirement, but is not a measured regression in this implementation.

Exit: memory budgets have a byte model and measured evidence; any speed
claim states its configuration and matched work. Baselines are recorded
before setting later regression ratios; no invented universal speedup target.

**G21 — Demonstrate that the test suite rejects wrong algorithms**

Dependencies: G03, G14, G17.
Tests: T05, T09, T12, T14, T15, T25, T26, T27, T28, T32, T33, T34, T36, T37, T40, T43, T50, T53, T62.

Run the critical mutation matrix against disposable variants of the real
production paths. Include numerical errors, premature projection, wrong
packing, missing score offsets, source-origin shortcuts, leakage, stale
caches and wrong measurement alignment. Each mutant must fail a targeted
invariant; parser failures alone do not count. Reduce uncovered failures
into a small deterministic regression case, then retain it.

Exit: all critical listed mutants are detected. Document any untested
failure class explicitly; a raw coverage percentage is not the gate.

**G22 — Verify the integrated feature and bind its evidence**

Dependencies: G15, G19, G20, G21.
Tests: T23, T44, T45, T46, T47, T52, T54, T56, T57, T58, T59, T60, T61, T62.

Run focused tests during development. Once the selected implementation scope
is stable, run the full suite, R CMD check, all admitted execution routes and
scheduled numerical/statistical courts. Register files and exports through
the existing architecture/API governance. Re-run persisting certification
runners serially per `benchmarks/RECERTIFY.md`, using an isolated install where
required; retain the source digest before/after, promote only matching
artifacts, and check summary/artifact consistency. Bind new evidence to its
oracle, seed and configuration as well as production code.

Exit: no unreviewed architecture relaxation, no unexplained failures or
stale-evidence skips, and an exact source/evidence record for the completed
scope. Existing optional-data follow-ups stay separately labelled. A commit
or publication is a subsequent action, not performed by this planning task.

**Later programs, with explicit entry criteria**

Condition generalization requires a training-defined feature/centering map
and cross-kernel, an explicit target centering on held-out conditions, and
independent-evidence admission. Its first oracle is the minimum-norm
continuation `F_BB = L F_AA L'`, not a new fitted eigensystem on test conditions.
Existing M10 is relevant context, not a dependency that silently grants this
capability to the first release.

Observation-weighted/GLS loss requires an observation map, frozen precision
and a different fitting solver. The score-adjoint law remains reusable.
Existing M11 is relevant context; its covariance capability must be checked
against current admitted designs.

Learned neural metrics and genuinely learned feature transformations need
their own training dependency graph. The fixed-feature factor objective in
the supplied proposal is already covered by G09; do not implement a second
optimizer merely to satisfy old M12 wording.

**Evidence from this planning task**

The current source and existing tests were inspected; the earlier focused
test run and exploratory numerical checks remain recorded in the previous
review. This turn added and ran exact fixture witnesses: 14 fixtures, with
16 deliberately wrong alternatives separated by their assertions. It did
not implement G01–G22, rerun the historical full suite, or alter production
source/certification artifacts. Work-package/test references and the
dependency graph are checked for completeness and cycles before delivery.
Validation found 22 packages, 62 unique assigned test specifications and 38
dependency edges; the table agrees with the individual records, the graph
is acyclic, and local artifact links resolve. The rechecked production digest
remains `sha256:24afa8ce392d71efa930c2dcb0afd9ac0951d100a424965c5a89a949355f1a8d`.
