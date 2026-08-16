# First-moment relation compiler and ingestion contract

Date: 2026-08-15

Mote epic: `bd-01M03BM44QBBNN50PWB1PPZ6EF`

## Outcome

Extend effectagram one level upstream so that raw observations become an
identified relation family through the same discipline already used for
second-moment geometry:

\[
\text{event and observation facts}
\longrightarrow
\operatorname{plan\_relation}
\longrightarrow
\text{design receipt}
\longrightarrow
B_r
\longrightarrow
\operatorname{plan\_geometry}
\longrightarrow
\text{views}.
\]

The estimand-bearing first-moment object is the composition

```r
plan_relation(
  study,
  model = design_model(...),
  effects = effect_map(...),
  observation_model = observation_model(...)
)
```

It is not `design_model()` alone. For observations \(Y\), compiled design
\(X\), condition-space functional lowered to \(T\), and fixed whitener \(L\),

\[
E=T(LX)^+L,
\qquad
B=EY.
\]

The existing `effect_extractor` is the concrete \(E\). `estimate(plan)` emits
the existing `effect_relation_fit`, so `plan_geometry()` and every downstream
view keep their current contract.

## Boundary prerequisite

The current version 0.1 worktree must be reviewed, committed, and certified as
an immutable baseline before first-moment implementation begins. Contract work
may proceed now, but version 0.2 implementation must not be mixed into the
uncommitted closure-calculus changes.

The present-tense estimability-message repair remains independent work under
`bd-01KZY9HQ0QSV0WEAR4NV7Z4VRJ`. It should land immediately and must not wait
for this architecture.

## Five distinctions

Every contract clause must place information in exactly one category.

1. **Facts:** observation indexes and source revisions; event timing and
   schema; observation-level confounds; partition hierarchy; neural domain.
2. **Requests:** design terms; condition-space effect functionals; the later
   `generalizes_over` choice.
3. **Assumptions:** observation/noise model, sampling unit, independence and
   exchangeability, fixed versus learned status.
4. **Receipts:** realized design matrices, row lineage, censor realization,
   ranks, aliases, solver, lowering route, whitening diagnostics, runtime and
   memory.
5. **Earned results:** fitted relation, residual/error channel, effect
   covariance, residual degrees of freedom, and capability record.

Dataset hierarchy supplies a typed vocabulary for a generalization request; it
never chooses that request. It also does not prove independence or
exchangeability. `sampling_unit` belongs to `observation_model()` and is
cross-checked against declared axes rather than inferred from event rows.

## Typed core

### Observation facts

`observations()` binds a lazy observation-by-feature source to one exact neural
domain and one typed observation index per partition. The index records stable
observation identifiers, time coordinates and units when available, and source
revision. Timing is a capability rather than a requirement for already reduced
or trial-level inputs.

### Event and confound facts

`events()` records onsets, durations, units, experimental variables, item
identifiers, and partition metadata. Columns carry schema, not permanent model
roles. A variable may be a target in one design and nuisance in another.

Observation-level motion, physiology, censor flags, and acquisition covariates
belong to a separate table on the observation axis. They are not forced into an
event table.

### Study binding

`study()` binds observation clocks, event clocks, confound rows, partitions,
and neural domains. It does not equijoin events to scans. The design compiler
projects event intervals onto observation times and records exact row lineage.

Misalignment refuses by class and names the partition, axes, units, observed
coverage, and remedy.

### Design model

`design_model()` is a pure declaration owned by a conforming first-level
compiler. It identifies HRF bases, term construction, modulation, drift,
nuisance and censor policies. Fixed HRF convolution is linear; selecting a
different HRF basis changes model identity. Adaptive censor or noise procedures
are declared learned operations and record their realized state in the receipt.

### Condition-space effect map

`effect_map()` is declared against named experimental conditions or other
semantic model coordinates, never positionally against compiled matrix
columns. The compiler lowers it to a concrete \(T\) for the chosen coding.

Treatment and cell-means codings of the same declared functional therefore
share one effect identity without comparing floating-point operators. Coding
and lowering details belong to the design receipt. Changing the condition
functional, units, scale, or targeted HRF component changes effect identity.

Raw \((X,T)\) input remains a degenerate adapter. Its parameterization is bound
into identity and it cannot claim symbolic coding invariance.

### Observation model

`observation_model()` owns OLS or GLS semantics, whitening, autocorrelation
policy, sampling unit, and fixed versus learned provenance. Conditional on a
compiled design, effect map, and fixed observation model, the neural-data map
\(Y\mapsto B\) is linear.

A fixed whitener may grant the admitted separable-GLM error channel. A learned
AR or GLS whitener records its training source and weaker capabilities; it does
not inherit the fixed-whitener analytic law by string declaration.

### Relation plan, receipt, and fit

`plan_relation()` composes the study, design declaration, condition-space
effect map, and observation model. `estimate()` is the only verb required to
read neural values; metadata and design compilation may run before it.

The identity ladder is:

- `design_model_id`: semantic event-to-regressor declaration;
- `effect_map_id`: condition-space scientific functional;
- `relation_plan_id`: composed scientific request and declared assumptions;
- `design_receipt_id`: realized design, rows, route, diagnostics and compiler;
- `relation_fit_id`: plan plus source revisions and realized receipt.

Changing QR to SVD, block size, storage, or another numerically equivalent
execution route does not change `relation_plan_id`. Receipts distinguish those
routes. Different realized censor masks under one declared policy retain the
request identity but produce different receipt and fit identities.

No identity comparison uses a floating-point tolerance. Numerical invariance
is a tested law with a declared tolerance, not an identity algorithm.

## External compiler conformance

Effectagram owns the protocol, identities, capabilities, and refusals. It does
not own every event grammar, HRF basis, censoring implementation, or fitting
engine.

A conforming compiler must return a versioned, portable record containing:

- semantic model identity and exact compiler version;
- named regressor axes and condition-space lowering;
- exact observation-row lineage and censor accounting;
- per-partition design matrices or stable revisions;
- ranks, aliases, estimability diagnostics and residual degrees of freedom;
- solver and numerical tolerances;
- observation-model and whitening provenance;
- fixed versus learned declarations with training revisions;
- error-channel capabilities; and
- source revisions sufficient to reproduce the fit identity.

`fmridesign` should own the event/HRF grammar, `fmrireg` the estimation and
noise engines, and effectagram the conformance protocol. BIDS is an adapter
into observation, event, confound and partition facts rather than the core data
model.

## Route-stability and refusal laws

The first law suite must establish:

1. changing only a compilation or numerical route preserves plan identity;
2. receipts distinguish routes and fit identity binds the realized receipt;
3. two supported codings of one condition-space effect share its semantic and
   plan identities;
4. changing the HRF target, effect functional, units, observation model, or
   sampling assumption changes plan identity;
5. hierarchy constrains but never infers `generalizes_over`;
6. unbiased cross-partition consumers require an explicit independence
   declaration or a capability that proves it;
7. fixed and learned observation models expose different capabilities;
8. raw \((X,T)\) adapters cannot claim coding invariance;
9. non-estimability names the partition, alias structure, failed functional,
   and remedies; and
10. misalignment, missing timing, unearned analytic uncertainty, and
    nonportable compiler receipts fail as classed capability refusals.

All numerical-route tests use floating-point fixtures and declared tolerances.
Exact-integer fixtures are insufficient evidence for route invariance.

## Delivery phases

### Phase 1 — Normative contract

Mote: `bd-01M03BNA49FF208RSWWFCK6EN0`

Write `ingestion-contract.md` with the algebra, object schemas, identity
ladder, capability implications, refusal vocabulary, conformance schema,
positive laws, negative fixtures, and explicit non-goals. Freeze public names
only after the laws can be stated without implementation exceptions.

### Phase 2A — Condition-space effects

Mote: `bd-01M03BNBGMWGBK41QD4NPVW94P`

Implement symbolic, axis-bound effect functionals and compiler lowering. Prove
coding invariance without numerical operator equivalence. Preserve raw \((X,T)\)
as the narrower-capability adapter.

### Phase 2B — Observation and study facts

Mote: `bd-01M03BNC0FBBA6KPHQKVK746BZ`

Implement observation indexes, event schema, observation-level confounds,
partition hierarchy, clock binding, row lineage and alignment refusals. Make
generalization, sampling, independence and exchangeability declarations
explicit at their consuming layers.

Phases 2A and 2B may proceed independently after the contract lands.

### Phase 3 — Relation plan and estimation

Mote: `bd-01M03BNCBXY031N8GEJ17MFSSK`

Implement `plan_relation()`, `observation_model()`, design receipts, fit
identity and `estimate()`. Adapt `lm_relation_fit()` rather than breaking the
existing downstream relation contract. Establish route-stable identity and
fixed-versus-learned capabilities.

### Phase 4 — External conformance

Mote: `bd-01M03BND69JEWY0JAF8ETZB9RC`

Publish the compiler protocol and version-pinned `fmridesign`, `fmrireg`, and
BIDS adapters. Validate them as installed consumers, not only from source
checkouts.

### Phase 5 — End-to-end evidence

Mote: `bd-01M03BNDXGEBXN0TFP5F3MV1CN`

Run one multi-partition raw-observation problem through facts, relation plan,
fit, geometry, RDM, RSA and admitted uncertainty. Require direct
linear-algebra oracles, legacy `lm_relation_fit()` parity, coding-invariance
evidence, unequal-run and censor cases, fixed and learned observation models,
negative fixtures, and map-scale resource receipts.

### Phase 6 — Journey and exact-artifact certification

Mote: `bd-01M03BNEZXF36GN5CCKSZY04Q3`

Ship a question-first vignette and API reference. Certify the exact commit with
the full suite, installed cross-package integration, R CMD check, pkgdown,
hosted scale gates, and reproducible benchmark receipts. Close the epic only
when those artifacts and all child work are complete.

## Evidence gates

The epic is complete only when all of the following are demonstrated:

1. **Algebraic parity:** `estimate(plan_relation(...))` reproduces the existing
   fixed-design `lm_relation_fit()` relation, residuals, effect covariance and
   residual degrees of freedom.
2. **Coding invariance:** two supported codings of one semantic effect have
   one plan identity and agree numerically against an independent oracle.
3. **Route stability:** QR/SVD and block/whole routes preserve plan identity,
   differ in receipt identity, and agree within declared tolerance.
4. **Realized-sample identity:** changing a censor realization preserves the
   declared request identity but changes receipt and fit identity.
5. **Declaration discipline:** hierarchy never silently selects
   generalization, sampling units, independence or exchangeability.
6. **Capability discipline:** learned observation models cannot acquire fixed
   analytic uncertainty without an explicit uncertainty model.
7. **Interoperability:** version-pinned installed adapters reproduce their
   native compiler and fitter outputs.
8. **Scientific vertical slice:** the resulting `effect_relation_fit` supports
   point geometry and, when earned, analytic uncertainty without a special
   downstream path.
9. **Operational evidence:** every new public execution path has correctness,
   runtime and memory gates, and the exact checked artifact is recoverable from
   version control.

## Explicit non-goals

- Reimplement every HRF, first-level engine, preprocessing step, or BIDS tool
  inside effectagram.
- Infer independence, exchangeability, generalization, or sampling units from
  dataset shape.
- Treat event columns as permanently target or nuisance variables.
- Identify scientific equality by approximate matrix comparison.
- Grant fixed-design covariance laws to learned whitening or adaptive
  preprocessing by default.
- Delay the current estimability-message repair until this architecture ships.
