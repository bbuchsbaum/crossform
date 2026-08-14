# Evidence-pairing and two-sided tomography contract

Status: normative architecture contract

Contract version: `evidence-pairing-v1`

Date: 2026-08-13

This document freezes the semantic foundation for
`bd-01KZWD1P7Y7DA5HY779MXGQ731`. It is normative for the evidence-pairing,
measurement-form, coupling, and tomography epics. The independent numerical
oracles in `tests/testthat/test-evidence-pairing-contract.R` certify the laws
stated here without calling production effectagram functions.

The contract extends, and does not replace, `effect-form-v1`. No production API
is introduced by this document. In particular, connectivity materialization is
deferred until its capability gates exist.

`evidence-pairing-v1` defines an estimand, not its sampling distribution. When
relations are estimated from data, covariance among repeated evidence
estimates is governed separately by `evidence-sampling-v1`. The statistical
contract consumes the same identified queries and may reuse the same compiled
sufficient statistics, but it neither changes the irreducible observable nor
adds uncertainty capabilities to a pure relation.

## 1. Four identified spaces

Let

\[
B_L\in\mathbb R^{q_L\times p_L},
\qquad
B_R\in\mathbb R^{q_R\times p_R}
\]

be a left and right experimental--neural relation. They bind four independent
spaces:

```text
experimental_left  D_L, dimension q_L
experimental_right D_R, dimension q_R
neural_left        N_L, dimension p_L
neural_right       N_R, dimension p_R
```

Every space identity includes its ordered coordinates or basis, units, scale,
and provenance signature. Equal dimensions, coordinate labels, or numerical
values are not identity. A query or measurement leg is compatible only when
its complete bound identity agrees with the relation boundary it closes.

An experimental pair query is

\[
H\in\mathbb R^{q_L\times q_R},
\]

bound in this order to \((\mathcal D_L,\mathcal D_R)\). A neural pair query or
bridge is

\[
K\in\mathbb R^{p_L\times p_R},
\]

bound in this order to \((\mathcal N_L,\mathcal N_R)\). Shape checking is
necessary but never sufficient.

## 2. The irreducible observable

The scalar evidence pairing is

\[
\boxed{
\mathscr E_{LR}(H,K)
=
\operatorname{tr}\!\left(H^\top B_L K B_R^\top\right).
}
\]

It is defined when all four identities match, all dimensions conform, and all
participating values are finite. It is bilinear in \(H\) and \(K\). It need not
be symmetric, positive, or nonnegative.

This scalar is the smallest completed scientific observable in the framework.
Effect forms and measurement forms are useful partial materializations of the
same pairing; neither is a second mathematical engine.

## 3. Adjoint second-order form transports

The relation pair induces the forward transport

\[
\mathcal P_{LR}(K)=B_LKB_R^\top
\]

and its Frobenius adjoint

\[
\mathcal P_{LR}^{*}(H)=B_L^\top H B_R.
\]

For the Frobenius inner product
\(\langle X,Y\rangle_F=\operatorname{tr}(X^\top Y)\),

\[
\boxed{
\mathscr E_{LR}(H,K)
=\left\langle H,\mathcal P_{LR}(K)\right\rangle_F
=\left\langle\mathcal P_{LR}^{*}(H),K\right\rangle_F.
}
\]

The precise interpretation is therefore:

> A pair of experimental--neural relations induces an adjoint pair of
> second-order form transports. Effect geometry and neural coupling are the
> forward and adjoint pictures of that transport.

They are not described as two loosely related partial contractions. Their
equality is an adjoint law with named domains and codomains.

## 4. The conceptual second-order lift

Using column-major vectorization,

\[
\boxed{
\operatorname{vec}\!\left[\mathcal P_{LR}(K)\right]
=(B_R\otimes B_L)\operatorname{vec}(K).
}
\]

The conceptual lifted relation is

\[
\mathbb B_{LR}=B_R\otimes B_L,
\]

and its Euclidean adjoint is

\[
\mathbb B_{LR}^{*}=B_R^\top\otimes B_L^\top.
\]

This identity is explanatory and testable on small fixtures. The compiler must
not materialize the Kronecker matrix by default, and no plan, result, or receipt
may store it as the ordinary representation of a task. Execution contracts
matrix products and factorizations directly. Explicit lift materialization is
permitted only as a bounded diagnostic or test operation carrying a deliberate
opt-in and size guard.

## 5. Open boundaries and semantic materializations

The pairing has four boundaries: \(\mathcal D_L,\mathcal D_R,\mathcal N_L,\)
and \(\mathcal N_R\). Closing different pairs produces distinct semantic
objects.

| Boundary state | Materialization | Mathematical object |
|---|---|---|
| Neural pair closed; experimental pair open | `effect_form` | \(F_K=B_LKB_R^\top\) |
| Experimental pair closed; native neural pair open | neural evidence form | \(Q_H=B_L^\top H B_R\) |
| Experimental pair closed; measured neural pair open | `measurement_form` | \(C_H=\mathcal L_LQ_H\mathcal L_R^\top\) |
| Both pairs closed | scalar evidence | \(\mathscr E_{LR}(H,K)\) |
| All four boundaries open | conceptual transport | \(\mathcal P_{LR}\), never a dense default result |

These public result types remain distinct because their axes, valid views, and
scientific capabilities differ. Internally they must be partial
materializations of one typed contraction plan. The architectural rule is:

> Unification at the compiler; specialization at the semantic boundary.

An output's open boundaries are part of task identity. A result cannot be
relabelled from one materialization into another merely because its matrix has
the same dimensions.

## 6. Identified measurement legs

A left measurement leg at node \(x\) and a right measurement leg at node \(y\)
are

\[
L_{L,x}\in\mathbb R^{d_x\times p_L},
\qquad
L_{R,y}\in\mathbb R^{d_y\times p_R}.
\]

They bind their source neural-space identities and ordered output measurement
spaces. For fixed \(H\), the measured cross-form is

\[
\boxed{
C_{xy,H}
=L_{L,x}B_L^\top H B_RL_{R,y}^\top
=L_{L,x}Q_HL_{R,y}^\top.
}
\]

The metric \(M_x=L_x^\top L_x\) is sufficient for a self-contracted local
geometry, but it is not a replacement for \(L_x\) when oriented output
coordinates or cross-location blocks matter. Two legs can induce the same
metric while differing by an output rotation; their raw cross-block entries
then differ. Plans that may produce measurement forms must preserve identified
legs, not only metrics.

For a query

\[
D\in\mathbb R^{d_x\times d_y}
\]

on the measured output coordinates, the induced neural query is

\[
K_D=L_{L,x}^\top D L_{R,y}.
\]

The complete four-leg law is

\[
\boxed{
\left\langle H,
B_LL_{L,x}^\top D L_{R,y}B_R^\top
\right\rangle_F
=
\left\langle D,
L_{L,x}B_L^\top H B_RL_{R,y}^\top
\right\rangle_F.
}
\]

This law is the execution equivalence between pushing a measured neural query
forward and pulling an experimental query backward.

## 7. Reversal law

Side reversal exchanges every left and right identity and transposes every
pair query:

\[
(B_L,B_R,H,K,L_{L,x},L_{R,y})
\longmapsto
(B_R,B_L,H^\top,K^\top,L_{R,y},L_{L,x}).
\]

It obeys

\[
\mathscr E_{RL}(H^\top,K^\top)=\mathscr E_{LR}(H,K),
\]

\[
\mathcal P_{RL}(K^\top)=\mathcal P_{LR}(K)^\top,
\qquad
\mathcal P_{RL}^{*}(H^\top)=\mathcal P_{LR}^{*}(H)^\top,
\]

and

\[
C_{yx,H^\top}=C_{xy,H}^\top.
\]

Reversal is a structural operation, not a boolean metadata flag.

## 8. Direct sums, decompositions, and factorized contraction

If the neural spaces are direct sums and \(K\) is block diagonal, forward
transport is additive over those blocks. This recovers the existing additive
spatial-frame law as a specialization.

If a measurement leg is vertically decomposed,

\[
L_x=
\begin{pmatrix}
L_{x,1}\\ \vdots \\ L_{x,k}
\end{pmatrix},
\]

then every edge form inherits the crossed block decomposition

\[
C_{xy}^{a,b}=L_{x,a}Q_HL_{y,b}^\top.
\]

This is the decomposition-lifting law. It applies equally to coherent versus
configuration components, frequency bands, cortical depths, anatomical versus
functional bases, and shared versus subject-specific modes. A subspace can be
canonical without a particular basis in that subspace being canonical.
Consequently raw entries require identified oriented bases; singular values,
Frobenius norms, canonical correlations, subspace angles, and geometry
alignment are invariant under appropriate within-subspace rotations.

Low-rank queries must be contractible without forming either complete form. If

\[
H=UV^\top,
\qquad K=RS^\top,
\]

then

\[
\boxed{
\mathscr E_{LR}(H,K)
=\left\langle U^\top B_LR,\;V^\top B_RS\right\rangle_F.
}
\]

This is both a required algebraic law and a compiler optimization opportunity.

## 9. Exact effect-form specialization

Closing the neural boundary first gives

\[
F_K=B_LKB_R^\top,
\]

which is exactly the logical object governed by `effect-form-v1`. A factorized
measurement bridge with legs \(L_{left}\) and \(L_{right}\) supplies

\[
K=L_{left}^\top L_{right}.
\]

The existing identity-bridge, ordered-partition-edge, reversal, storage-codec,
and query-fusion laws remain unchanged. In a self construction with a common
leg, the familiar geometry

\[
G=BL^\top LB^\top
\]

is the self-adjoint positive specialization. Cross-generalized self forms can
remain symmetric but indefinite exactly as specified by `effect-form-v1`.

No evidence-pairing implementation may weaken an existing effect-form
capability or silently change its canonical operation order.

## 10. Construction guarantees and observed diagnostics

Capabilities are justified by construction provenance. Diagnostics are facts
observed from a finite numerical realization. They are stored separately.

Construction capabilities include, where applicable:

```text
query_role: effect | variation
positive_query: true | false
self_construction: true | false
symmetric: true | false
guaranteed_psd: true | false
joint_covariance: true | false
sampling_axis: time | trial | subject | session | ...
normalization_axis: neural_features | experimental_samples |
                    partition_pairs | form_entries
normalization_policy: explicit named policy
```

Diagnostics include, under a declared tolerance and algorithm:

```text
effective_rank
condition_number
observed_min_eigenvalue
observed_max_eigenvalue
zero_variance_coordinates
regularization_used
```

A positive sampled eigenvalue cannot create a `guaranteed_psd` capability. A
matrix that happens to be symmetric cannot create a symmetric construction
capability. Effective rank is a diagnostic; the role and sampling provenance
that make repeated variation scientifically meaningful are structural.

Every result receipt must preserve both classes and state which claims are
guarantees versus diagnostics.

### Evidence capability is not calibration capability

For an estimated relation, the construction may additionally declare whether
an identity-bound error channel is available. This declaration remains
separate from the relation value and from the evidence result:

```text
error_channel: relation_fit | external | absent
sampling_covariance: available | unavailable
metric_status: fixed | learned
metric_uncertainty: not_applicable | ignored | propagated
partition_model: equal | heterogeneous
```

A pure precomputed relation is fully capable of point evidence while lacking
within-participant sampling covariance. Numerical shape, edge count, or a
crossvalidated estimator cannot manufacture the missing capability. Detailed
sampling and refusal laws belong to `evidence-sampling-v1`.

## 11. Effect coupling is broader than connectivity

For every conformable \(H\),

\[
C_{xy,H}=L_xB_L^\top H B_RL_y^\top
\]

is an algebraically valid `effect_coupling`. The name makes no covariance or
sampling claim.

A `covariance_coupling` additionally requires a declared repeated-variation
axis and a valid positive covariance construction. `canonical_coupling`
additionally requires compatible positive self-blocks, a coherent joint
covariance construction, and an explicit regularization policy. Gaussian
mutual information additionally requires a coherent joint Gaussian covariance
model. A validating `connectivity()` convenience may accept only inputs whose
receipts establish the required capabilities.

### The rank-one semantic guard

For a single contrast in a self relation,

\[
H=cc^\top,
\qquad u=B^\top c,
\qquad Q_H=uu^\top.
\]

After two measurements, with \(a=L_xu\) and \(b=L_yu\),

\[
C_{xy}=ab^\top,
\quad C_{xx}=aa^\top,
\quad C_{yy}=bb^\top.
\]

Whenever both sides are nonzero, unregularized linear geometry alignment is
one and the only nonzero canonical correlation is one. Those perfect values
are algebraic degeneracy, not estimated connectivity.

The governing semantic rule is:

> Do not contract away the axis whose repeated variation is supposed to define
> normalized coupling.

The rank-one block remains a valid effect coupling. It must not acquire a
connectivity capability from its numerical value.

## 12. Axis-aware normalization and ordered stages

The word `correlation` is not a complete operation. Every normalization records
the axis on which moments are estimated:

- pattern correlation normalizes over `neural_features`;
- functional-connectivity correlation normalizes over
  `experimental_samples`;
- a partition-pair standardization normalizes over `partition_pairs`;
- geometry alignment or entry transforms act on `form_entries`.

These operations can have identical matrix shapes and different meanings.
Axis identity is mandatory and belongs in plan identity.

The certified effect-form specialization retains its exact order:

```text
ordered partition-pair product
  -> neural-feature spatial normalization
  -> edge transform
  -> partition reduction
  -> experimental pair query
```

Future evidence tasks carry an explicit ordered stage ledger. At minimum it
records raw sufficient-statistic construction, each axis-qualified
normalization, each nonlinear edge transform, partition aggregation, boundary
closure, and requested materialization. A policy may deliberately normalize
within each partition pair before aggregation or after aggregation from
separately accumulated positive self-moments. Those are different estimands
and therefore different task identities.

Normalization, Fisher or rank transformation, partition aggregation, and
boundary closure may be fused or reordered only under a certified algebraic
law. In particular,

\[
\operatorname{normalize}\!\left(\sum_e w_e C_e\right)
\ne
\sum_e w_e\operatorname{normalize}(C_e)
\]

in general. Crossvalidated or otherwise indefinite self-blocks must not be
silently used as correlation denominators.

The stage ledger describes how edge contributions form an estimand. It does
not imply that edge contributions are independent sampling replicates. In
particular, unordered cross-partition products share endpoints even when the
underlying partition estimates are independent. Any standard error, interval,
or test must use a declared sampling law rather than the number or empirical
spread of reducer edges. `evidence-sampling-v1` governs this distinction.

## 13. Two-sided block calculus and tomography

Experimental sameness and neural sameness are independent axes:

| | Same neural measurement | Different neural measurements |
|---|---|---|
| Same experimental space | ordinary effect geometry | neural coupling |
| Different experimental spaces | ER-RSA / cross-domain effect form | cross-domain, cross-region form |

Thus ER-RSA reads experimental off-diagonal blocks while coupling reads neural
off-diagonal blocks. Ordinary geometry lies on both diagonals.

For fixed \(H\), stack identified measurement legs as

\[
\mathcal L_L=\begin{pmatrix}L_{L,1}\\ \vdots \\ L_{L,m}\end{pmatrix},
\qquad
\mathcal L_R=\begin{pmatrix}L_{R,1}\\ \vdots \\ L_{R,n}\end{pmatrix}.
\]

The complete measured block form is

\[
\boxed{
\mathcal C_H=\mathcal L_LQ_H\mathcal L_R^\top.
}
\]

Its diagonal blocks are within-location evidence and its off-diagonal blocks
are between-location evidence. If both stacked frames have full column rank,
let

\[
S_L=\mathcal L_L^\top\mathcal L_L,
\quad S_R=\mathcal L_R^\top\mathcal L_R,
\]

and canonical duals

\[
\widetilde{\mathcal L}_L=\mathcal L_LS_L^{-1},
\qquad
\widetilde{\mathcal L}_R=\mathcal L_RS_R^{-1}.
\]

Then

\[
\boxed{
Q_H=\widetilde{\mathcal L}_L^\top
\mathcal C_H
\widetilde{\mathcal L}_R.
}
\]

For Parseval frames, \(S_L=S_R=I\). The corresponding reconstruction uses the
ordinary stacked legs. This is the spatial completion law: local node blocks
plus edge blocks can be a lossless coordinate representation of the global
neural evidence operator. Local-geometry conservation preserves total
evidence; complete node--edge tomography can preserve the operator itself.

The same dual-frame law applies on the experimental side to frame
representations of \(F_K\). The grand structural statement is therefore:

> The framework is a two-sided tomography of one second-order relation.
> Experimental forms and neural forms are dual frame representations of the
> same evidence pairing.

Tomography claims require validated frame rank, recorded tolerances, explicit
dual construction, and numerical reconstruction evidence. A collection of
local maps alone must never claim operator completeness.

## 14. Required executable law court

Before this contract can close, independent tests must cover:

1. the forward/adjoint/scalar equality;
2. the column-major vectorization identity on explicitly bounded small lifts;
3. refusal to materialize a Kronecker lift by default;
4. four-leg measured-query equivalence;
5. effect-form specialization;
6. reversal and transpose laws;
7. direct-sum and low-rank factorization laws;
8. same-shaped but differently identified boundary rejection;
9. separation of construction capabilities from numerical diagnostics;
10. rank-one normalized-coupling degeneracy and capability rejection;
11. axis-aware normalization identity and noncommutation with aggregation;
12. deterministic randomized and adversarial-scale fixtures.

All comparison tolerances must be explicit. Randomized fixtures must use fixed
seeds and report enough dimensions to reproduce a failure. Tests in this stage
are a contract court, not production conformance tests.

## 15. Scope boundary

For estimated relations, this contract's observable is the estimand and
`evidence-sampling-v1` is its statistical continuation. If

\[
\widehat e_a=\widehat{\mathscr E}(H_a,K_a),
\]

then the adjacent contract governs

\[
\mathcal V_{ab}=\operatorname{Cov}(\widehat e_a,\widehat e_b).
\]

That separation is constitutional: the point-evidence contraction remains
testable without a statistical apparatus; a `relation_fit` supplies a distinct
error channel; and calibration plans may reuse but never redefine the evidence
plan. The first admitted implementation is the fixed-metric, equal-partition
crossvalidated-distance covariance in Diedrichsen et al. (2016), with exact
query-first actions. Learned metrics, heterogeneous partitions, and spatial
cross-location covariance require later sampling laws.

This contract authorizes the later epics to build one private contraction-plan
IR and distinct semantic materializations. It does not itself authorize:

- a public `connectivity()` API;
- automatic all-pairs searchlight materialization;
- dense default Kronecker construction;
- inferred sampling semantics from matrix shape;
- CCA, correlation, or mutual-information claims without capability evidence;
- adaptive or data-learned measurement legs without explicit fitting
  provenance and cross-fitting where scientific validity requires it;
- replacement of identified measurement legs by metrics when orientation
  matters;
- inference from partition-edge spread or edge count without an admitted
  sampling law;
- calibration of precomputed relation values that lack an identified error
  channel.

Those boundaries are part of `evidence-pairing-v1`, not temporary omissions.
