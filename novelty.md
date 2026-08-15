# What is novel in effectagram?

`effectagram` is an estimand compiler for second-order neuroimaging analysis.
Its primary architectural contribution is an **evidence-pairing calculus** in
which familiar analyses arise by closing different experimental and neural
boundaries. Its software contribution is the **executable estimand contract**
that prevents those scientific objects from changing silently during
execution.

That is the claim being developed. This page distinguishes the parts already
demonstrated from the stronger interpretations that still have to be earned.

## The category difference

The calculus separates two choices that are usually bundled together: whether
the experimental spaces are the same, and whether the neural measurements are
the same.

Each cell carries its evidence status from the ledger below, so the table
claims territory only at the strength actually earned.

| | Same neural measurement | Different neural measurements |
|---|---|---|
| **Same experimental space** | retained signed marginals; crossvalidated contrast energy; ordinary representational geometry *(demonstrated)* | cross-node effect coupling *(implemented; small-node contract)*; normalized connectivity only when its sampling contract is met |
| **Different experimental spaces** | rectangular ER-RSA or cross-task geometry *(implemented; real exemplar pending)* | rectangular cross-domain, cross-region forms *(prospective)* |

Signed activation is a retained first-moment marginal, not a bilinear statistic.
Its reproducible energy, squared-distance geometry, and the other entries in
the table are second-order closures. Keeping both channels in one plan lets a
single fitted relation report the signed effect together with the geometry that
asks whether it reproduces.

The irreducible second-order observable is

\[
\mathscr E_{\Gamma}(H,K)
=
\sum_{r,s}\Gamma_{rs}
\operatorname{tr}\!\left(
H^\top B_{L,r}K B_{R,s}^\top
\right).
\]

Here the relations \(B_{L,r}\) and \(B_{R,s}\) bind named experimental and
neural spaces, \(H\) asks an experimental-side question, \(K\) specifies a
neural-side metric or bridge, and \(\Gamma\) declares which independently
estimated relations must reproduce one another. Closing the neural boundaries
gives an effect form

\[
F_K=\sum_{r,s}\Gamma_{rs}B_{L,r}K B_{R,s}^\top,
\]

while closing the experimental boundaries gives the adjoint neural form

\[
Q_H=\sum_{r,s}\Gamma_{rs}B_{L,r}^\top H B_{R,s}.
\]

The scalar pairing can be read from either side:

\[
\mathscr E_{\Gamma}(H,K)
=\langle H,F_K\rangle_F
=\langle Q_H,K\rangle_F.
\]

This is more than saying that one kernel computes several distances. It gives
one typed construction for self- and cross-forms, square and rectangular
experimental axes, and local and cross-location neural measurements.

## Relationship to RSA and `rsatoolbox`

RSA, crossnobis distance, noise-precision estimation, searchlight RSA, model
fitting, noise ceilings, and inference over subjects and conditions are
established work. [`rsatoolbox`](https://rsatoolbox.readthedocs.io/) is a
sophisticated reference implementation of an RDM-centered workflow: it
supports multiple dissimilarities including crossnobis, residual-based noise
precision, searchlights, fixed and flexible models, and inferential procedures
for different generalization targets. Its methods paper explicitly treats
generalization over subjects and conditions
([Schütt et al., 2023](https://doi.org/10.7554/eLife.82566)). `effectagram`
does not claim any of those ingredients as inventions.

Nor is a shared interface across ROIs, searchlights, whole-brain maps, and
multiple multivariate measures new by itself. The Decoding Toolbox
([Hebart et al., 2015](https://doi.org/10.3389/fninf.2014.00088)),
CoSMoMVPA
([Oosterhof et al., 2016](https://doi.org/10.3389/fninf.2016.00027)), and
PyMVPA
([Hanke et al., 2009](https://doi.org/10.3389/neuro.11.003.2009)) established
powerful unifying software abstractions. The claim here is the exact algebraic
relationship among the outputs, not the fact that one package can dispatch
several analyses.

The common fixed-linear subset is simple. If a square self-form is
\(G=BKB^\top\), squared-distance extraction is a linear map \(d=\mathcal D(G)\).
A fixed linear RSA readout \(\beta=Cd\) can therefore be compiled into a query
\(H_\beta\) such that

\[
\beta=\langle H_\beta,G\rangle_F.
\]

This observation licenses a controlled parity comparison; it is not by itself
a scientific novelty result. A public `rsatoolbox` parity analysis is still an
open evidence gate. The current external comparison is with rMVPA.

The intended distinction is that an RDM is an optional view of an effect form,
not the mandatory intermediate representation. The principal `rsatoolbox`
`RDMs` abstraction stores dissimilarities over one shared pattern axis as a
vector or square matrix. `effectagram` instead keeps the identified relation,
left and right axes, generalization operator, measurement frame, and error
channel, when the relation was fitted with one, available until the requested
query has been compiled. This is an architectural comparison, not a claim that
every nonlinear dissimilarity or inferential method in `rsatoolbox` is
contained in the bilinear core.

In particular, Pearson correlation distance is outside that core because its
per-pattern norm is nonlinear in the fitted patterns. The package gives that
boundary a named policy rather than a quiet escape hatch; see
[the correlation-distance policy](correlation-distance-policy.md).

## What is distinctive

### 1. Boundary-closure unification

Crossvalidated contrast energy, squared-distance RDMs, fixed linear RSA,
ordered cross-domain similarity, and neural-side effect coupling are queries or
partial materializations of the same evidence pairing. This is implemented and
covered by independent small-matrix law tests. Normalized connectivity has
additional repeated-variation and covariance requirements; a numerical
off-diagonal block does not acquire that interpretation from its shape alone.

The neural-side construction has close scientific neighbors in
[representational connectivity](https://pmc.ncbi.nlm.nih.gov/articles/PMC2605405/),
[informational connectivity](https://pmc.ncbi.nlm.nih.gov/articles/PMC3566529/),
and multidimensional-connectivity methods
([Basti et al., 2020](https://doi.org/10.1016/j.neuroimage.2020.117179)). The
distinctive claim is not ownership of connectivity; it is that neural coupling
is the adjoint-side materialization of the same identified experimental-neural
pairing, with stronger normalized views admitted only by their own contracts.

### 2. Exact, non-destructive coherent/configuration accounting

For an admitted spatial metric, `effectagram` partitions the metric into a
rank-one coherent mode and its configuration complement. Pulling the partition
through the relation gives

\[
G^{\mathrm{total}}
=G^{\mathrm{coherent}}+G^{\mathrm{configuration}}.
\]

Every fixed linear downstream query inherits the same additive partition. The
defensible contribution is therefore not the elementary matrix identity, nor a
claim to have invented the distinction between mean activation and pattern
structure. It is an exact, crossvalidated, query-preserving accounting that
keeps both components instead of destructively demeaning the data and calling
the remainder "purely multivariate." The implementation and numerical laws
exist; the one-plan flagship public demonstration is still an evidence gate.

This claim sits beside a substantial literature on differences between
univariate and multivoxel analyses
([Davis et al., 2014](https://doi.org/10.1016/j.neuroimage.2014.04.041)),
[mean-pattern subtraction](https://pmc.ncbi.nlm.nih.gov/articles/PMC3786542/),
pattern-component modelling
([Diedrichsen et al., 2018](https://pubmed.ncbi.nlm.nih.gov/28843540/)), and
second-moment accounts of RSA and encoding models
([Diedrichsen and Kriegeskorte, 2017](https://doi.org/10.1371/journal.pcbi.1005508)).

### 3. Rectangular, identified cross-domain forms

The left and right experimental spaces may be different, unequal, and
directed. An encoding-by-retrieval form, for example, need not be forced into a
square condition axis. `plan_geometry(x, at, over, right = )` compiles the
rectangular plan publicly; axis-bound `pair_query()`s execute against it
query-first, and `geometry()` materializes a rectangular form that satisfies
the exact algebraic identity `total = coherent + configuration`. The public
test verifies that identity numerically to `1e-12`. The engine, independent
oracles, and public constructor exist. A real scientific analysis with missing
or unequal items and pair-space covariates has not yet been shipped, so the
stronger empirical claim remains gated.

### 4. Spatial frames with a qualified conservation law

Voxels, regions, searchlights, and whole-brain summaries enter as measurement
frames rather than separate analysis engines. Under a conservative,
feature-additive frame,

\[
\sum_x L_x^\top L_x=M_\Omega
\quad\Longrightarrow\quad
\sum_x B L_x^\top L_x B^\top=B M_\Omega B^\top.
\]

This law is implemented and tested. Today it should be read as structure the
framework provides, not yet as a headline scientific contribution. The public
demonstration still needs to show which overlap-accounting error the law
prevents, compare against a mass-preserving global measurement, and state its
feature-additive fixed-metric scope. Conservative feature weights do not by
themselves certify arbitrary non-diagonal or learned metrics.

### 5. Query-first execution

Selected contrasts, distance edges, and fixed linear RSA coefficients compile
without requiring the complete RDM as the public intermediate: every RDM edge
is the rank-one operator \((e_i-e_j)(e_i-e_j)^\top\), and the kernels evaluate
it as two row differences and a Hadamard product instead of materializing a
dense packed query. The recorded large-\(q\) gate measures the consequence at
100 conditions over 1,080 searchlights: one hundred selected pairs in 0.23 s,
the full fused 4,950-coordinate RDM in 4.6 s against 18.3 s for
materialize-then-project, with a direct-oracle error of `4.4e-16`. Incremental
peak RSS was 361 MB (345 MiB), below the 512 MiB gate, while structured
execution avoided a separate ~200 MB (191 MiB) dense query allocation.

### 6. Generalization bound to estimand identity

Existing RSA inference already recognizes that generalization over
measurements, subjects, and conditions changes the scientific claim. The
distinctive formalization here is narrower: \(\Gamma\) participates in the
identity of **every** effect-form estimand, not only the final inferential
procedure. The axis is typed, not inferred from labels:
`cross_partitions(relation, generalizes_over = "run")` and the same call with
`"session"` produce distinct plan identities even under identical generic
partition names and identical fold counts, and identity tests enforce it.
Runs, sessions, tasks, item sets, sites, and ordered cross-domains can all be
represented by named pairing relations.

## The contract is the proof mechanism

The calculus defines the scientific objects. The contract is how those objects
survive real execution:

- plan identity records the relation, queries, metric, frame, units, and
  generalization relation;
- receipt identity records storage, tiling, execution path, and numerical
  diagnostics without redefining the plan;
- capabilities such as symmetry, self-form status, positive semidefiniteness,
  fixed-metric status, and retained uncertainty are construction guarantees;
- optimized paths are required to carry independent numerical oracles; and
- refusals are first-class results when an interpretation has not been earned.

The retained error channel is an important consequence. For the admitted
fixed-metric equal-partition model, an RSA coefficient is a fixed linear
functional of the RDM, so its covariance transports exactly from the analytic
RDM covariance of
[Diedrichsen, Provost, and Zareamoghaddam (2016)](https://doi.org/10.48550/arXiv.1607.01371).
Precomputed effects without an identified error channel are not
reverse-engineered from the spread of their edges. Refusals are classed
conditions: `catch_refusal()` returns the missing capability, every unmet
reason, and remedies as data, and `sampling_capabilities()` answers the
admission question before it is provoked. The executable
[failure gallery](failure-gallery.md) shows five realistic errors that the
package guards against. Callable unsupported interpretations return classed
refusals; clipping is absent rather than offered as a biased option; and
changes in generalization produce distinct estimand identities.

## Evidence ledger

Every claim carries one of four statuses — **established algebraically**
(proved and independently law-tested as a mathematical identity),
**implemented** (reachable through the public package path and executed by
compiled plans), **demonstrated** (validated in a substantive comparison,
generative recovery, or realistic example), or **prospective** — and every
row links the machine-checkable artifact that certifies it, so the ledger can
be audited rather than believed.

| Claim | Status | Evidence artifact and boundary |
|---|---|---|
| Two-sided evidence-pairing laws | **Established algebraically** | Forward, adjoint, scalar, rectangular, reversal, decomposition, and tomography identities against independent bounded oracles: [`helper-evidence-pairing-laws.R`](tests/testthat/helper-evidence-pairing-laws.R), [`helper-effect-form-laws.R`](tests/testthat/helper-effect-form-laws.R), [`test-tomography.R`](tests/testthat/test-tomography.R). Algebraic and software evidence, not a scientific benchmark. |
| Crossnobis point parity | **Demonstrated** | The Haxby 2001 exemplar agrees with an independent loop to `1.33e-15` and rMVPA to `8.88e-16` over 577 VT searchlights: [`exemplars/haxby2001`](exemplars/haxby2001/). Matched crossvalidated squared-Euclidean/crossnobis estimand, not correlation distance. |
| Error-bearing refit and linear uncertainty transport | **Demonstrated under an admitted model** | Refit reproduces the point RDM to `4.44e-16`; a fixed linear RSA coefficient consumes factorized analytic covariance: [`exemplars/haxby2001`](exemplars/haxby2001/), [`test-evidence-sampling-kernel.R`](tests/testthat/test-evidence-sampling-kernel.R), [`test-evidence-sampling-generative.R`](tests/testthat/test-evidence-sampling-generative.R). Boundary: the declared equal-partition, fixed-metric, separable plug-in model. |
| Query-first execution at scale | **Demonstrated** | Recorded gate artifact at q = 100 over 1,080 searchlights: selected 100 pairs 0.23 s, full fused RDM 4.6 s vs 18.3 s materialize-then-project, oracle error `4.4e-16`, incremental peak RSS 361 MB (345 MiB): [`benchmarks/run-query-first-scale.R`](benchmarks/run-query-first-scale.R), [`benchmark-results/query-first-scale-gate.rds`](benchmark-results/), [`test-query-first-scale.R`](tests/testthat/test-query-first-scale.R). |
| Coherent/configuration family | **Demonstrated** | One plan yields the signed contrast, the three energies with exact recomposition, the RDM, the RSA coefficient, and the admitted analytic SE, with planted-effect recovery (signal carried by configuration; null regions centred on zero): the [introduction vignette](vignettes/introduction.Rmd) with executable checks, plus [`test-integrity-guards.R`](tests/testthat/test-integrity-guards.R) and [`test-measurement-decomposition.R`](tests/testthat/test-measurement-decomposition.R). A real-data decomposition narrative on Haxby remains desirable but is no longer the gate. |
| Rectangular cross-domain analysis | **Implemented** | Public constructor, query-first pair queries, oracle parity, and materialized recomposition: [`test-rectangular-plan.R`](tests/testthat/test-rectangular-plan.R). A real match/control, pair-covariate analysis is Gate 3. |
| Adjoint coupling from the plan vocabulary | **Implemented** | `coupling(plan, between, by)` compiles the adjoint closure against the plan's own frame and pairing, small-node contract enforced: [`test-coupling-views.R`](tests/testthat/test-coupling-views.R). |
| Generalization axis in estimand identity | **Implemented** | Typed `generalizes_over` bound into task and metric-pairing identity; cross-run vs cross-session distinct under identical generic labels: [`test-generalization-axis.R`](tests/testthat/test-generalization-axis.R). |
| Route-stable view identity | **Implemented** | Fused query-first and materialize-then-project executions of one estimand carry one scientific id with distinct execution receipts: [`test-route-identity.R`](tests/testthat/test-route-identity.R), [`test-effect-form-certification.R`](tests/testthat/test-effect-form-certification.R). |
| Spatial conservation | **Implemented** | `frame_conservation()` diagnostic plus public tests with the unnormalized global comparator and the total-only boundary stated: [`test-integrity-guards.R`](tests/testthat/test-integrity-guards.R). The overlap-accounting demonstration is Gate 4. |
| Refusal discipline | **Implemented** | Classed `effect_capability_refusal` conditions with all-reasons reporting, `catch_refusal()`, `sampling_capabilities()`, and the executable [failure gallery](failure-gallery.md): [`test-capability-refusals.R`](tests/testthat/test-capability-refusals.R). |
| Cross-package speed advantage | **Not demonstrated** | The recorded matched-estimand comparison is internal (fused vs materialized routes of one plan). No cross-package speedup is claimed; that would require matched estimands, independent parity, warm-up, repeated timings, hardware, and map-scale budgets against the other implementation. |

Reproducible Haxby scripts and qualifications are in
[`exemplars/haxby2001`](exemplars/haxby2001/). Versioned performance evidence
lives under [`benchmarks`](benchmarks/).

## Claim-promotion gates

The stronger novelty statement is a target, not a documentation shortcut. Its
sections are promoted only as these gates land:

1. **`rsatoolbox` specialization:** reproduce a fixed crossnobis plus linear-RSA
   result with a version-pinned environment and an independent oracle. *Open.*
2. **Inherited decomposition family:** one plan yields signed, coherent,
   configuration, total, RDM, and RSA outputs, with additive recomposition and
   planted-effect interpretation. **Landed:** the
   [introduction vignette](vignettes/introduction.Rmd) demonstrates the full
   family from one plan with executable recomposition and planted-truth
   checks.
3. **Real rectangular analysis:** distinct unequal axes, match/control coding,
   pair-space covariates, missing items, and explicit operation order. *Open;
   the public rectangular constructor and query path landed.*
4. **Operational spatial accounting:** overlapping searchlights, coverage
   correction, a mass-preserving global comparator, and an interpretive
   consequence. *Open.*
5. **Large-\(q\) query-first execution:** selected queries without a full RDM,
   with numerical, memory, and runtime receipts. **Landed:** see the
   [query-first scale gate](benchmarks/run-query-first-scale.R) and its
   recorded result artifact.
6. **Failure gallery:** executable safeguards for correlation-distance
   normalization, negative clipping, destructive demeaning claims,
   generalization identity changes, and learned-metric leakage. **Landed:**
   see the [failure gallery](failure-gallery.md). Callable unsupported
   interpretations use classed refusals; the gallery also documents absent
   transformations and estimand-identity guards.

The decisive comparison is not "which RSA toolbox has more methods?" It is:
reproduce an established RDM-first result as a specialization, then obtain a
scientifically meaningful result—such as inherited coherent/configuration
accounting or a real rectangular cross-domain analysis—that requires retaining
structure the square RDM view has already collapsed.

## What is not claimed

`effectagram` does not claim to have invented RSA, crossnobis, analytic RDM
covariance, searchlights, noise ceilings, condition/subject generalization, or
connectivity. It does not yet claim empirical superiority to `rsatoolbox`, a
matched-estimator speed advantage, or a complete population-inference layer.
Nor does it claim that every nonlinear RDM comparison belongs in the bilinear
core.

The current, defensible novelty is the conjunction of a broader typed
evidence-pairing architecture and an implementation that makes its scientific
identity, uncertainty preconditions, and refusal boundaries executable. The
scientific importance of its strict extensions will be upgraded only when the
public gates above supply the evidence.
