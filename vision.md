# Vision for `effectagram`

Status: working vision
Date: 2026-08-13 (calculus framing noted 2026-08-15)

> Framing update, 2026-08-15: the public claim has since been sharpened into
> the **evidence-pairing calculus** recorded in `novelty.md` — self- and
> cross-experimental axes, self- and cross-neural measurements, and spatial
> frames as boundary closures of one typed second-order pairing, with the
> executable estimand contract as the mechanism that keeps those closures
> intact through execution. This document's rMVPA-successor mission is
> unchanged; the calculus is its organizing mathematics.

The name is pronounced “effect-a-gram.” An effectagram is the
cross-generalized geometry of experimental effects across a declared spatial
frame. The package name is broader and more memorable than the mathematical
class name; durable objects therefore retain precise names such as
`effect_relation`, `effect_geometry`, and `effect_view`.

`effectagram` will become the foundation for task-fMRI effect analysis in R and the eventual successor to `rMVPA`. It will replace the scientific role of `rMVPA`, not reproduce its source code, public classes, or every historical feature. Researchers should be able to perform the important analyses for which they now turn to `rMVPA`, while gaining a clearer account of what was estimated, where it was measured, across which partitions it generalized, and how the reported number was derived.

The project will achieve this through unification rather than a smaller catalogue of methods. Voxelwise, searchlight, regional, surface, and whole-brain analyses will use common relations, spatial measurements, generalization declarations, and result geometry. Contrasts, representational dissimilarities, RSA, predictive summaries, and later population models will be interpretations or extensions of those shared objects. The system can grow broad without making every new analysis another engine, registry entry, result class, and combining path.

A beautiful abstraction should produce a beautiful execution plan. Conceptual
unity, statistical validity, and computational efficiency are coequal design
requirements: the package is not successful if its cleanest estimand is too
slow or too large to use on a brain.

Unification must improve execution as well as explanation. Common algebra exposes work that can be shared, moved, factorized, or avoided: sufficient statistics can be reused across locations, low-rank queries can be contracted before large forms are built, sparse supports can replace dense operators, and adjoint identities can select the cheaper of equivalent contraction orders. The project will pursue such exact reductions aggressively. It will not accept inherited method boundaries as computational necessities when a new factorization, sparse representation, matrix-free action, or streaming algorithm can remove the repeated work.

Speed never licenses a change in the estimand. Every optimized lowering must agree with an independent reference implementation under a declared numerical contract. An approximation, reduced-precision calculation, or altered estimator must be named and validated as a different scientific option; it may not silently replace the exact path. Performance claims must report their dimensional regime, runtime, peak memory, and materialization costs on realistic workloads. The aim is to make the principled analysis the practical analysis: fast, memory-efficient, and exactly faithful to the question the user declared.

The long-term goal is therefore larger than a compact version 0.1 package:

> Build a better system for the work `rMVPA` set out to support, with effect geometry as its stable scientific and computational foundation.

## The future system

A fixed bilinear analysis will follow this irreducible sequence:

\[
\text{partitioned relation}
\xrightarrow[\text{pairing}]{\text{spatial frame}}
\text{effect geometry}
\xrightarrow{\text{scientific view}}
\text{reported values}.
\]

The full system will grow around that sequence:

```text
response data, effects, designs, and external estimators
                         |
                  input adapters
                         |
                         v
                  effect_relation
                         |
             spatial frame + pairing
                         |
                         v
                  effect_geometry
                         |
        +----------------+----------------+
        |                |                |
   fixed views     learned methods   population and
   and maps        with held-out     calibration models
                   evaluation
        |                |                |
        +----------------+----------------+
                         |
                  reports and tools
```

These are dependency boundaries, not necessarily names of future packages. The core must not depend on the workflow, learning, inference, or presentation layers above it. Those layers may grow as separate packages or focused modules when their assumptions and dependencies warrant it.

In R, the foundational workflow should read approximately as follows:

```r
rel <- effect_relation(
  sources = list(run1 = B1, run2 = B2, run3 = B3),
  effects = effect_space(
    c("face", "house", "object"),
    basis_id = "condition-means-v1"
  ),
  domain = volume_domain(mask)
)

g <- geometry(
  rel,
  at = searchlights(
    domain(rel),
    radius = 6,
    units = "mm",
    normalization = "local"
  ),
  over = cross_partitions(rel)
)

result <- geom_contrast(g, c(face = 1, house = -1))
```

The exact function names may change while the contracts are tested. The sequence should not: define the experimental–neural relation, spatial measurement, and generalization requirement before asking for an interpretation.

## What replacement means

`rMVPA` proved that one R system can support regional and searchlight decoding, RSA, cross-domain analysis, specialized representational models, real-data demonstrations, parallel execution, diagnostics, and command-line workflows. It also accumulated separate model families, execution modes, engine-selection rules, fast paths, custom callbacks, result types, and adapters. A change in scientific method can therefore require changes across several software ontologies.

`effectagram` will replace that organizing model. It will aim for comparable or better coverage of important scientific jobs, but it will not seek one-for-one API or feature parity. Each inherited workflow will receive one of four explicit dispositions:

1. **Core geometry:** the analysis is a fixed relation, frame, pairing, and geometry query.
2. **Disciplined extension:** the analysis trains a model or transform, freezes it, and evaluates it on independent relations through shared contracts.
3. **Adjacent system:** the task concerns estimation, calibration, population modeling, orchestration, or presentation and belongs in a focused layer around the core.
4. **Deliberate retirement:** the behavior depends on ambiguous semantics, duplicated machinery, or an unrestricted callback that defeats the common contracts.

Replacement will be measured by user journeys, not exported-symbol counts. A researcher should be able to move from an experimental design or set of effect estimates to a regional, searchlight, or global result; compare representational models; test generalization across runs or experimental domains; inspect diagnostics; and render interpretable outputs without returning to `rMVPA`. The route need not use the old function names, and some routes should make different scientific choices explicit rather than preserve old defaults.

## Interpretive advantage

The unification is valuable because it changes what the software lets researchers see.

`effectagram` will not treat “univariate” and “multivariate” as two kinds of neural signal. For additive spatial measurements, it will retain both the coherent regional effect and reproducible spatial configuration. Researchers can inspect what a regional mean explains and what remains in the spatial pattern without destructively demeaning the data or comparing separately fitted pipelines.

Generalization will be part of the estimand rather than an execution detail. Explicit pairings will state whether an effect must reproduce across runs, sessions, tasks, stimulus sets, or experimental domains. Cross-validation will no longer be identified merely by a fold count.

Every estimated effect has a sampling law that is conceptually distinct from
its value. The system will state whether that law is available. Generalization
edges state where an effect must reproduce; they are not automatically
independent observations for calibration. When a fitted relation carries an
identified error channel, uncertainty can be transported over the same
evidence queries without changing their estimands. When only precomputed
effects are supplied, unavailable within-participant uncertainty will be stated
as a capability boundary rather than guessed from edge spread.

Geometry will remain queryable before scalarization. A contrast energy, RDM coefficient, eigenvalue, classifier summary, and map may discard different information from the same underlying relation. Preserving the parent geometry as a compiled operator makes those losses visible and allows several views without refitting, while dense packed storage remains an explicit materialization rather than a prerequisite.

Finite-sample cross-product estimates will remain signed and may be indefinite. Negative dissimilarities and eigenvalues will not be silently repaired. Information, effective-rank, and manifold summaries will require a model that supplies the positive-semidefinite object and units they assume.

The result will not be a universal theory that makes every method equivalent. Fixed bilinear geometry, learned prediction, nonlinear representation, and calibrated population inference have different assumptions. The system will connect them without erasing those distinctions.

## One core, broad capabilities

The numerical core should remain small enough for a maintainer to audit completely. Breadth will come from stable contracts around it.

### Inputs and estimation

Relations may be backed by raw response data, beta or contrast images, condition means, trialwise estimates, latent representations, or external estimators. Input adapters will establish orientation, experimental basis, units, feature identity, partitions, and provenance. Rich design, HRF, preprocessing, and noise-estimation systems may construct extractors, but they will not become hidden branches in the geometry kernel.

### Spatial measurement

One frame abstraction will cover points, smoothing kernels, searchlights, regions, surfaces, multiscale bases, and whole-brain measurements when their metrics can be stated explicitly. Dense or learned metrics may use later representations, but they must not create separate regional and searchlight engines.

### Learned analyses

Classification, local covariance estimation, feature selection, alignment, and other adaptive methods cannot be disguised as fixed geometry queries. A learning layer will declare its training relations, return a frozen transform or model, and evaluate it on independent relations. It may derive predictive summaries, but it must preserve training provenance and expose the geometry or representation that supports interpretation when one exists.

### Inference and population analysis

Permutation tests, uncertainty estimates, hierarchical models, prevalence, geometry transport, and multiplicity control will attach to declared geometry or views. They will not be baked into every result type. Population models should consume geometry before reducing each participant to accuracy or one RSA coefficient whenever the scientific question requires the fuller object.

### Presentation and workflow

Neuroimaging adapters, renderers, reports, command-line tools, and interactive applications can provide the convenience expected of a mature analysis system. They will compile to the same explicit objects and must not introduce a second hidden analysis API.

## Rules that prevent renewed sprawl

The following rules apply to the whole system, not only version 0.1:

- `effect_relation` and `effect_geometry` are the durable scientific intermediates. New modules exchange these objects or small values with equally explicit contracts.
- Spatial scope is a frame property. Searchlight, regional, surface, and global computations do not become top-level method families.
- Generalization is an explicit pairing or sampling declaration, not an implicit loop hidden inside a model.
- Scalarization occurs as late as practical. A scalar result must retain a route back to the relation, geometry, assumptions, and units that produced it.
- Learned operations separate training, freezing, and held-out evaluation. Training data identities and fitted transformations are part of provenance.
- Calibration remains distinct from effect estimation. A p-value or interval does not define the effect, and cross-validation does not by itself provide inference.
- Sampling covariance is itself a queryable operator. It should be transported directly to the requested contrast, RDM model, or other view and materialized only when the complete covariance is scientifically required.
- Extensions may not use an unrestricted per-location callback as their principal integration mechanism. They must state their inputs, outputs, invariants, and failure conditions.
- Optimization must remove work rather than alter the scientific question. Every optimized lowering requires an independent reference implementation, equivalence tests, and realistic runtime and memory gates. Approximate or randomized algorithms require an explicit error contract and a distinct estimator identity; performance never justifies silent semantic drift.
- Compatibility adapters may read common artifacts and help users migrate, but the project will not reproduce `rMVPA`'s class hierarchy or dispatch behavior inside the core.
- A feature that cannot reuse the common contracts belongs in an adjacent module, requires a deliberate extension of the algebra, or should not be added.

## How we set the stage now

The first release should remain narrow, but its contracts must anticipate the larger system. That means preserving experimental-space identity, metric units, partition semantics, feature-domain identity, provenance, estimator status, and resource requirements from the beginning. Retrofitting those concepts after users depend on loosely labelled matrices would recreate the very problem the project is meant to solve.

Development should proceed along six tracks:

1. **Prove and accelerate the core.** Build an independent dense oracle and executable algebraic laws, then derive exact sparse, streamed, query-fused, and matrix-free lowerings from them. Set realistic brain-scale runtime and memory gates before broad workflows. Prefer eliminating algebraic work to distributing avoidable work.
2. **Maintain a replacement map.** Track important `rMVPA` user journeys and classify each as core geometry, disciplined extension, adjacent system, or deliberate retirement. Record current support and the evidence needed to certify a replacement.
3. **Build comparative exemplars.** Reproduce selected regional, searchlight, RSA, and cross-domain analyses on synthetic and public datasets. Compare matched estimands, scaling, and sampling rules—not merely similar-looking outputs or runtime.
4. **Design migration without inheritance.** Provide adapters and guides for common matrices, `neuroim2` domains, effect images, masks, atlases, and design outputs. Explain conceptual translations from `rMVPA`; do not emulate its model specifications and result classes.
5. **Make extensions prove reuse.** Before adding a learning, population, surface, or workflow module, demonstrate that it consumes or produces the common intermediates and does not require a parallel engine hierarchy.
6. **Document tasks and interpretations.** Documentation should begin with research questions and real workflows, show the estimand and assumptions, and then reveal the common algebra. Users should experience a better system, not be required to study an abstract compiler before doing an analysis.

This work should be staged rather than collapsed into the first package release. The additive geometry kernel earns the architecture. Practical volume analyses and migration exemplars earn adoption. Learned, population, surface, and workflow layers earn the claim that the architecture can support a complete successor system.

## Version 0.1 responsibility

Version 0.1 will establish the constitutional core: validated dense and block-backed relations, explicit linear extractors, additive nonnegative sparse frames, normalized off-diagonal pairings, one packed geometry result, dense and streamed implementations, and a small set of contrast and representational views. It will support abstract and volume domains and demonstrate voxel, searchlight, regional, and global measurements through one kernel.

Version 0.1 will not attempt to reproduce the full `rMVPA` experience. It will omit classifier registries, arbitrary local callbacks, formula and HRF languages, local learned covariance, population inference, distributed orchestration, and method-specific result families. These omissions protect the foundation; they are not declarations that the future system will never support prediction, inference, surfaces, or convenient workflows.

## Definition of success

The core succeeds when its optimized results reproduce an independent brute-force calculation, its units and generalization semantics are explicit, voxel/searchlight/region/global computations differ only through declared measurement frames, and its fast paths meet stated runtime and memory gates at realistic scale without changing the estimand.

The broader system succeeds when:

- the important scientific journeys now served by `rMVPA` have clearer effectagram-based replacements or documented reasons for retirement;
- adding a new analysis usually means adding an extractor, frame, pairing, query, learner, or downstream model rather than another execution engine;
- researchers can inspect coherent effect, spatial configuration, generalization structure, geometry, and scalar summaries without rerunning unrelated pipelines;
- predictive and adaptive methods expose how they were trained and how held-out evidence was evaluated;
- public-data exemplars and matched-estimand tests substantiate scientific and numerical claims;
- realistic benchmarks show that shared algebra, sparse representations, query fusion, and streaming eliminate avoidable work without sacrificing exactness;
- practical adapters and documentation make the improved semantics usable rather than merely theoretically attractive; and
- new analyses no longer need `rMVPA` as their default implementation substrate.

The project should call itself the successor to `rMVPA` only after it has earned that status across real user journeys. The destination is explicit now so that every early contract supports it: a broader, more interpretable, and more maintainable analysis system built on one trustworthy geometry foundation.
