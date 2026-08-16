# First-moment relation compiler contract

Status: normative architecture contract

Contract version: `first-moment-relation-v1`

Date: 2026-08-15

This document freezes the scientific and software semantics for
`bd-01M03BM44QBBNN50PWB1PPZ6EF`. It governs the path from event and observation
facts to the relation families consumed by `plan_geometry()`. The independent
oracles in `tests/testthat/test-ingestion-contract.R` certify the initial
algebraic and identity laws without calling the production implementation.

The first-moment layer applies the same discipline already used by the
second-moment layer: an estimand-bearing plan is distinct from its execution
receipt, axes are named and bound rather than positional, capabilities are
earned by construction, and unsupported interpretations refuse explicitly.

## 1. Scope and category statement

The layer identifies and estimates a family of experimental--neural relations

\[
B_r\colon \mathcal E\longrightarrow\mathcal N,
\]

one for each declared partition \(r\). Its public chain is

```text
event and observation facts
  -> plan_relation()
  -> design receipt
  -> estimate_relation()
  -> effect_relation_fit
  -> plan_geometry()
  -> views
```

`plan_relation()`, not `design_model()`, is the first-moment analog of
`plan_geometry()`. The scientific request exists only after composing the
study, model, condition-space effects, and observation model:

```r
plan_relation(
  study,
  model = design_model(...),
  effects = effect_map(...),
  observation_model = observation_model(...)
)
```

The contract does not turn crossform into an event grammar, HRF library,
BIDS indexer, or general GLM engine. Those facilities may be supplied by
conforming packages. Crossform owns the identities, conformance boundary,
capability implications, refusals, and the relation object emitted downstream.

## 2. The linear relation and its conditional scope

For one partition, let

- \(Y\in\mathbb R^{n\times p}\) be observations by neural features;
- \(X\in\mathbb R^{n\times k}\) be the compiled design;
- \(L\in\mathbb R^{n\times n}\) be the declared observation whitener; and
- \(T\in\mathbb R^{q\times k}\) be the lowering of a named condition-space
  effect map into the compiled coefficient basis.

The estimator is the existing linear extractor

\[
\boxed{E=T(LX)^+L},
\qquad
\boxed{B=EY}.
\]

Conditional on the compiled design, lowered effects, and a fixed observation
model, the map \(Y\mapsto B\) is linear. This statement is deliberately
narrow:

- fixed HRF convolution is linear, but selecting the basis changes the model
  request;
- a fixed censor mask selects rows linearly, but adaptive censoring produces a
  realized receipt;
- a fixed whitener may support the admitted separable-GLM covariance law, but
  estimating an AR or GLS model from the observations is learned and has
  weaker capabilities; and
- fitting and selecting among nonlinear design or noise models is outside the
  conditional linear claim.

`effect_extractor` remains the concrete representation of \(E\).
`estimate_relation()` must emit the existing `effect_relation_fit` contract so no
second-moment API needs a parallel ingestion-specific relation type.

## 3. Five disjoint information classes

Every field belongs to exactly one class.

| Class | Contents | May affect |
|---|---|---|
| Facts | source revisions, observation indexes and clocks, event schema and timing, confound rows, partition hierarchy, neural domain | study and fit identity |
| Requests | model terms, HRF targets, named condition-space functionals, units and scale; later, a `generalizes_over` selection | semantic plan identity |
| Assumptions | observation/noise law, sampling unit, fixed or learned status, independence and exchangeability declarations | plan identity and earned capabilities |
| Receipts | compiled matrices, coding, lowering, row lineage, realized censor mask, ranks, aliases, solver, tolerances, diagnostics, runtime and memory | receipt and fit identity only |
| Earned results | relation estimates, residual/error channel, effect covariance, residual degrees of freedom, capability record | downstream admissibility |

Facts do not silently become assumptions. A partition hierarchy supplies the
typed vocabulary for a later generalization request; it does not choose that
request and does not prove that partitions are independent or exchangeable.
Likewise, event rows do not determine a sampling unit. `sampling_unit` is an
explicit declaration in `observation_model()` and is checked against available
study axes.

## 4. Typed fact objects

### 4.1 `observations()`

An observation record binds, per partition:

- a lazy observation-by-feature source and stable source revision;
- one exact neural-domain identity;
- ordered, unique observation identifiers;
- optional time coordinates and their physical units; and
- acquisition facts needed to interpret observation-level rows.

Rows are never joined by position alone. A source without stable revision may
be explored, but cannot earn portable fit identity. Timing is a capability:
already reduced trial-level observations may be valid without a scan clock,
whereas an event-to-scan compiler may require one and must refuse if absent.

### 4.2 `observation_events()`

An event record contains typed schema and factual values: onset, duration,
units, experimental variables, item identifiers, and partition membership.
Columns do not carry permanent target/nuisance roles. Role is model-relative:
the same variable can define an effect in one `design_model()` and a nuisance
term in another.

### 4.3 `observation_confounds()`

Motion, physiology, acquisition covariates, spike regressors, and censor flags
live in a separate table bound to the observation axis. An event table is not
used as a container for scan-level values. Confound rows preserve stable
observation identifiers and source revisions.

### 4.4 Partition hierarchy

The study declares named axes and nesting, for example trials within runs,
runs within sessions, and sessions within subjects. Axis identity includes its
ordered levels and parent relationship. Equal labels at different hierarchy
positions are not interchangeable axes.

The hierarchy validates requests such as `run_axis(study)` but never inserts a
default `generalizes_over`. Cross-run and cross-session requests therefore
remain distinct even when both happen to compile to the same number of edges.

### 4.5 `study()`

`study()` binds observations, events, confounds, partitions, clocks, and neural
domains. It is a correspondence declaration, not an equijoin of event onsets
to scan times. A conforming compiler projects event intervals and model bases
onto the observation grid and records exact input-to-design row lineage.

Binding validates, before neural values are read:

- partition membership and uniqueness;
- clock units and origins;
- observation/confound row agreement;
- event coverage relative to observation windows; and
- neural-domain and source-revision consistency.

Ambiguous clocks, missing required timing, uncovered rows, or mismatched
confounds refuse with the partition, both axes, units, observed coverage, and
specific remedies.

## 5. Request and assumption objects

### 5.1 `design_model()`

`design_model()` is a pure semantic declaration supplied by a conforming
compiler. It identifies at least:

- compiler protocol and model-language versions;
- event-to-regressor terms and their semantic coordinates;
- HRF or temporal basis and targeted components;
- interactions and parametric modulation policy;
- drift and nuisance construction;
- censor policy; and
- any fixed preprocessing that changes the mean model.

Changing any of these semantic choices changes `design_model_id`. A change in
parameter coding, column order, QR pivot, block size, or solver does not change
the model request when the compiler proves it is a route for the same semantic
mean structure; it changes the design receipt.

### 5.2 `effect_map()`

An effect map is a finite, named linear functional on a semantic model space,
normally a condition-by-basis space. It binds:

- the complete ordered semantic coordinates;
- named output effect coordinates;
- numerical functional weights expressed in that semantic vocabulary;
- units and scale; and
- the targeted temporal-basis component when the model has more than one.

It is never declared positionally against compiled design columns.

Let \(\mu\) be semantic condition parameters, let a parameterization satisfy
\(\mu=C\beta\), and let a declared functional be \(F\mu\). The compiler lowers
it to

\[
T=FC.
\]

If two full-rank codings are related by \(X_2=X_1A\), their lowered targets
must satisfy \(T_2=T_1A\). Consequently

\[
T_2(LX_2)^+L=T_1(LX_1)^+L
\]

for the estimable functional. Coding invariance is thus a property of the plan
language and lowering proof, not a tolerant comparison of hashed floating
point operators. `effect_map_id` hashes the symbolic functional, its bound
semantic axes, units, and scale; it never hashes a lowered \(T\).

Changing the functional, condition ordering, units, scale, or HRF target
changes `effect_map_id`. Relabeling an effect without changing its bound
identity is forbidden.

### 5.3 Raw `(X, T)` adapter

Already compiled design and target matrices remain supported as a degenerate
input route. Their values, column order, parameterization, and revisions are
part of request identity because no semantic compiler proof is available.
This route cannot claim:

- `symbolic_effects`;
- `coding_invariant`;
- event-to-design row lineage that was not supplied; or
- semantic equivalence to another raw parameterization.

Numerically equal fitted effects do not upgrade those capabilities.

### 5.4 `observation_model()`

The observation model declares:

- OLS, WLS, GLS, or another admitted error law;
- sampling unit and its bound study axis;
- fixed versus learned whitening/noise status;
- independence or covariance assumptions actually consumed by estimation;
- estimation/training provenance for learned components; and
- residual degrees-of-freedom semantics.

`"independent"` is not an implicit default. Any downstream unbiased
cross-partition interpretation must consume either an explicit independence
declaration or a capability earned from a stronger design contract.

A fixed whitener can earn `fixed_observation_model`. A learned whitener records
its training sources and realized parameters in the receipt. Calling a learned
matrix "fixed" after fitting does not retroactively earn the fixed-whitener
sampling law; it may instead earn a conditional or frozen-learned capability
whose limitations remain visible.

## 6. `plan_relation()` and the identity ladder

`plan_relation(study, model, effects, observation_model)` is an immutable
scientific request. Creating it may validate facts and compile metadata, but it
must not read neural values. `estimate_relation(plan)` is the data-touching verb.

The identity ladder is:

1. `study_id`: factual axes, schemas, domains, and declared source revisions;
2. `design_model_id`: semantic event-to-mean declaration;
3. `effect_map_id`: semantic condition-space functional;
4. `observation_model_id`: declared error law and sampling assumptions;
5. `relation_plan_id`: the composed request `(study, model, effects,
   observation_model)`;
6. `design_receipt_id`: the realized compiler route and diagnostics; and
7. `relation_fit_id`: plan identity plus exact observation source revisions and
   realized receipt identity.

`relation_plan_id` derives only from canonical semantic declarations. It is
stable under QR versus SVD, whole versus blocked fitting, design-column order,
supported coding changes, and equivalent fixed computation routes. It is not
stable under changes to event facts, semantic terms, requested effects, units,
sampling assumptions, or the fixed/learned observation-model declaration.

Receipts distinguish every execution route. Different realized censor masks
under the same declared adaptive policy retain the plan identity but produce
different receipt and fit identities. A fit cannot be substituted for another
fit merely because its plan identity matches.

No scientific identity is defined by a numerical tolerance. Floating-point
route equality is a law tested at a declared tolerance; it is not an identity
algorithm.

## 7. Design receipt conformance

A conforming compiler returns one portable, versioned receipt per plan and
partition. At minimum it contains:

- protocol name and version;
- compiler package, version, and implementation revision;
- `study_id`, `design_model_id`, `effect_map_id`, and
  `observation_model_id` consumed;
- exact ordered design-regressor axis and semantic provenance for every
  column;
- condition-space-to-coefficient lowering with named axes;
- observation-row lineage from source identifiers to retained design rows;
- censor policy and realized inclusion/exclusion with reasons;
- realized design revision or portable matrix representation;
- rank, singular values or equivalent diagnostics, alias groups, tolerance,
  and per-effect estimability;
- solver and factorization route;
- whitening/noise provenance, including training source revisions when
  learned;
- residual degrees of freedom and their derivation;
- source revisions needed for fit identity;
- capabilities granted and the evidence for each; and
- optional runtime, memory, blocking, and storage diagnostics.

The receipt must be finite, axis-bound, canonically serializable, and free of
live closures, external pointers, and session-only addresses. A stable digest
is computed over the complete portable record. Omitting a required field does
not fall back to an "unknown" string; conformance refuses and identifies the
missing field.

An adapter may store large compiled matrices outside the receipt only through
a stable, verified source revision plus dimensions, codec, and content digest.

## 8. Earned fit and capability implications

The fit contains one named relation matrix per partition and preserves the
existing `effect_relation_fit` interface. It may also carry a residual source,
effect covariance, residual degrees of freedom, sampling unit, and conformance
receipt.

Capabilities describe construction guarantees, not observed numerical shape.
The initial implications are:

```text
symbolic_effects + valid_lowering
  => coding_invariant

aligned_observations + complete_row_lineage + estimable_effects
  => identified_relation

identified_relation + stable_source_revision + portable_design_receipt
  => portable_relation_fit

residual_channel + valid_df + fixed_observation_model + separable_glm_law
  => analytic_effect_covariance
```

None of the reverse implications holds. In particular:

- a numerically unchanged relation does not prove coding invariance;
- residuals alone do not grant an analytic covariance law;
- an invertible realized whitener does not prove that it was fixed; and
- a relation family does not prove cross-partition independence.

Capability records name both granted and missing guarantees. Learned
observation models must expose their training provenance and the weaker law
under which downstream uncertainty is conditional, approximate, resampled, or
unavailable.

## 9. Generalization remains a downstream request

`plan_relation()` produces partition-indexed first-moment relations. It does
not select which partition axis a later analysis generalizes over.

The later request remains explicit, for example:

```r
cross_partitions(
  fit,
  generalizes_over = run_axis(study),
  independence = declared_independence(...)
)
```

The exact public spelling may evolve before release, but these semantics may
not: `generalizes_over`, independence, and exchangeability cannot be inferred
from fold count, labels, nesting, or a silent default. They participate in the
identity of the downstream estimand that consumes them.

## 10. Required refusals

Contract-level failures use `effect_capability_refusal` and carry
`capability`, `namespace`, all `reasons`, and concrete `remedies`. Shape and
type errors that do not claim a scientific interpretation may remain ordinary
validation errors.

The initial refusal vocabulary is:

| Capability | Refuses when | Required context |
|---|---|---|
| `aligned_observations` | clocks, partitions, confounds, or row identifiers do not bind | partition, axes, units, coverage, offending rows |
| `timing_resolved` | an event compiler requires timing that facts do not supply | partition, required clock, available alternatives |
| `estimable_effects` | a named functional lies outside a partition design's row space | partition, rank, regressor count, aliases, failed effects, remedies |
| `valid_effect_lowering` | semantic coordinates cannot be lowered to the compiled model | missing coordinates, model term, requested functional |
| `declared_sampling_unit` | a consumer needs a sampling unit that was not declared or does not bind | requested law, available axes, declaration remedy |
| `declared_independence` | unbiased cross-partition use lacks an admitted independence guarantee | paired partitions, requested target, remedies |
| `fixed_observation_model` | an analytic law is requested from learned whitening provenance | training source, available conditional/resampling alternatives |
| `portable_design_receipt` | a compiler omits required lineage, provenance, version, or digest fields | compiler, protocol version, missing fields |
| `stable_source_revision` | a portable fit is requested from mutable or unidentified sources | partition and source |
| `analytic_effect_covariance` | residual, df, separability, or fixed-model premises are absent | every unmet premise and valid alternatives |

Validation occurs as early as the needed facts allow. Alignment, schema,
effect-axis, and conformance failures occur before neural values are read.
Estimability occurs at design compilation, before fitting neural blocks.

## 11. Normative laws and independent fixtures

The implementation must satisfy all laws below.

### Identity laws

1. QR versus SVD, blocked versus whole fitting, and supported coefficient
   codings preserve `relation_plan_id`.
2. Those route changes produce distinct receipt identities where their
   realized routes differ.
3. `relation_fit_id` binds the exact receipt and source revisions.
4. Changing the effect functional, HRF target, units, observation law,
   sampling unit, or fixed/learned declaration changes plan identity.
5. Raw `(X, T)` parameterization changes raw-route plan identity even when a
   floating-point oracle finds the same fitted values.

### Numerical laws

6. Cell-means and treatment codings of one declared condition-space
   functional yield the same extractor and fitted relation to a declared
   tolerance on a non-integer fixture.
7. QR and SVD routes yield the same estimable relation to a declared tolerance
   on a floating-point fixture.
8. Unequal partition lengths compile and fit independently while preserving
   one shared effect-space identity.
9. `estimate_relation(plan)` agrees with direct \(T(LX)^+LY\) and with the legacy
   `lm_relation_fit()` route for their common fixed-design subset.
10. Residual and effect-covariance outputs agree with an independent separable
    GLM oracle only when all advertised premises hold.

### Discipline laws

11. Event schema never silently assigns target or nuisance roles.
12. Hierarchy constrains but never selects `generalizes_over`.
13. Sampling unit and independence remain explicit declarations.
14. Learned and fixed observation models expose different capabilities.
15. Every required failure in Section 10 is classed, names all unmet reasons,
    and offers at least one scientifically valid remedy.
16. Semantic and conformance failures that can be known from metadata occur
    before the first neural source read.

Independent fixtures certify the algebra and identity model, while the public
vertical-slice court exercises every production type through the downstream
geometry and covariance APIs. The independent oracle remains separate so the
production implementation never becomes its own judge.

## 12. External compiler adapters

The conformance protocol, not a particular package, is the core boundary.

- `fmridesign` may own formulas, event expansion, HRF bases, and design
  construction.
- `fmrireg` may own fitting and learned noise engines.
- BIDS adapters may map `events.tsv`, imaging metadata, and confounds into typed
  facts.
- Other compilers may participate if they emit the same portable receipt and
  pass the conformance court.

An adapter is a morphism into the typed core. The core is not BIDS-shaped and
does not make a dependency mandatory merely because one adapter supports it.
Integration certification uses installed package versions and records them in
receipts; source-checkout agreement alone is insufficient evidence.

## 13. Compatibility and migration

The current `effect_extractor()`, `lm_extractor()`, and `lm_relation_fit()` APIs
remain valid for explicit linear designs. They are the degenerate ingestion
path and keep their honest, narrower capabilities. The 0.2 layer may implement
them through shared internals but must not silently grant symbolic provenance
or change their numerical estimand.

`plan_geometry()` continues to consume the existing relation-fit output. No
`ingested_geometry`, `first_level_rsa`, or parallel downstream engine is
introduced.

The pairing default records independence as undeclared. A public path that
needs cross-partition independence must require an explicit declaration or an
earned guarantee, and the choice remains part of estimand identity.

## 14. Explicit non-goals for version 0.2

Version 0.2 does not promise:

- every BIDS derivative or modality;
- a new HRF formula language inside crossform;
- nonlinear or Bayesian first-level inference;
- automatic discovery of exchangeability or independence;
- semantic equivalence testing for arbitrary user-supplied `(X, T)` pairs;
- unconditional analytic uncertainty after learned noise estimation;
- replacement of `fmridesign` or `fmrireg`; or
- population-level inference over subjects or conditions.

The vertical slice is successful when one conforming raw-observation workflow
produces an identified, residual-bearing relation family; matches independent
linear-algebra and legacy oracles; preserves condition-functional identity
across supported codings; refuses invalid interpretations; and proceeds
unchanged through geometry, RDM, RSA, and admitted uncertainty.
