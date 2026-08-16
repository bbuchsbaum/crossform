# Universal effect-form contract

Status: normative architecture contract

Contract version: `effect-form-v1`

Date: 2026-08-12

This document freezes the semantics used by the effect-form implementation
epic. It is normative for every child of
`bd-01KZW4WEZSJYMDK9XKB78N28FR`. The independent numerical oracles in
`tests/testthat/helper-effect-form-laws.R` never call production crossform
functions; `tests/testthat/test-effect-form-laws.R` runs the laws stated here
against the package, and `tests/testthat/helper-geometry-component-oracle.R`
supplies the first-principles total/coherent/configuration oracle used by the
view tests. This contract does not itself change production code
or the supported 0.1 API.

## 1. The logical object

An effect form is an ordered bilinear form

\[
F\colon \mathcal E_L\times\mathcal E_R\longrightarrow\mathbb R,
\qquad F\in\mathbb R^{q_L\times q_R},
\]

between a named `left_space` and a named `right_space`. Each axis is bound by
the complete `effect_space` identity: ordered coordinates, basis, units,
scale, provenance, and signature. Equal dimensions or labels are not identity.
Changing the order on one axis changes that space and requires the data and
query on that axis to be permuted with it.

A pair query is a finite matrix

\[
H\in\mathbb R^{q_L\times q_R},\qquad
\langle H,F\rangle_F=\operatorname{tr}(H^\top F).
\]

It is bound to both ordered space identities. The direct scalar primitive is
`t(b_left) %*% H %*% b_right`; the complete-form and direct-query paths must
produce the same scalar.

Direction is structural. Reversing a form exchanges the spaces, relation
families, partition endpoints, bridge legs, and query axes and produces
`t(F)`. A boolean called `directed` is not part of the new logical object.

## 2. Capabilities are guarantees

The following capabilities describe construction, not an observation made by
examining floating-point values:

- `self_form`: `left_space` and `right_space` are the same complete identity.
- `symmetric`: the construction guarantees `F == t(F)`. This implies
  `self_form`; a numerically symmetric rectangular result does not acquire the
  capability.
- `guaranteed_psd`: the construction guarantees
  `t(c) %*% F %*% c >= 0` for all `c`. This implies `symmetric` and
  `self_form`.

Cross-generalized self forms are symmetric when their ordered-edge/reducer
construction is self-adjoint, but they are not guaranteed PSD. Negative
eigenvalues and negative query values therefore remain valid estimates.
Self-products with a nonnegative spatial metric and identical measurement legs
may be guaranteed PSD; sampled numerical tolerance does not define that
capability.

Capability implication is thus

```text
guaranteed_psd  =>  symmetric  =>  self_form
```

and never the reverse.

## 3. Logical coordinates and physical codecs

The universal rectangular codec is column-major vectorization. For
`F[i, j]`, the one-based coordinate is

\[
k(i,j)=i+(j-1)q_L.
\]

The left coordinate therefore varies fastest, followed by the right
coordinate. A rectangular query uses the identical order, so
`sum(vec(H) * vec(F)) == sum(H * F)`.

A form carrying the `symmetric` capability may use the existing isometric
packed codec. It traverses the lower triangle column by column,

```text
(1,1), (2,1), ..., (q,1), (2,2), ..., (q,q),
```

and multiplies off-diagonal entries by `sqrt(2)`. For symmetric `H` and `F`,

\[
\operatorname{svec}(H)^\top\operatorname{svec}(F)
=\langle H,F\rangle_F.
\]

Rectangular and symmetric-packed are storage codecs, not distinct scientific
objects. Codec identity belongs in a storage manifest and execution receipt,
not in the scientific estimand. A packed codec is forbidden without the
`symmetric` guarantee.

`complete_form` and `query_only` are result capabilities, not codecs. A
complete result can answer any later compatible query and must be able to read
every logical coordinate. A query-only `effect_view` contains values for its
declared query and cannot claim completeness or be upgraded by relabeling,
even when a chosen set of queries happens to span the form.

## 4. Ordered partition edges and the self specialization

The kernel primitive for one ordered partition edge is never symmetrized. For
left relation `B_l` and right relation `B_r`, it forms the ordered outer
product after applying the declared spatial measurement or bridge:

\[
P_e=B_{l(e)}K_eB_{r(e)}^\top .
\]

An ordered edge table contains at least `left`, `right`, and a finite reducer
weight. Endpoint order has semantic meaning. Duplicate handling, if admitted,
is an explicit compiler operation and preserves canonical order and total
weight.

The existing undirected `cross_partitions()` estimator remains an exact
self-adjoint specialization. Each current unordered edge `(a, b, w)` expands
to

```text
(a, b, w / 2)
(b, a, w / 2)
```

and the linear reducer sums the two ordered products. Consequently

\[
\frac w2 B_aB_b^\top+\frac w2 B_bB_a^\top
=w\,\operatorname{sym}(B_aB_b^\top),
\]

which is the existing estimator. The old stored orientation remains
scientifically irrelevant. This expansion, not implicit kernel
symmetrization, is what grants the new result its `symmetric` capability.

## 5. Measurement bridges

Relations on distinct neural spaces cannot be multiplied by dimension alone.
A cross-space task requires a fixed `measurement_bridge` with

- a left leg `L_left` of shape `k x p_left`;
- a right leg `L_right` of shape `k x p_right`;
- one exact common-measurement-space identity of dimension `k`;
- the two source neural-space identities; and
- portable provenance describing how the fixed legs were obtained.

The bridge induces

\[
K=L_{left}^\top L_{right},
\]

but arbitrary dense `K` is not a public substitute for the factorized bridge.
For relation matrices with features in columns,

\[
B_LKB_R^\top
=(B_LL_{left}^\top)(B_RL_{right}^\top)^\top.
\]

The bridge identity covers both source identities, ordered legs and their
values, the common-space identity, and provenance. Reversing a bridge exchanges
the legs and source identities and transposes the induced form. Incompatible
leg dimensions or common-space identities fail before source reads.

When both sides have the exact same neural-space and feature identity, no
explicit bridge means the canonical identity bridge. It is numerically equal
to explicit identity legs. Merely same-shaped but differently identified
spaces have no implicit bridge and must fail before source reads.

Because both legs enter one common space, a joint self-product constructed from
the measured left and right relations is PSD. This fact does not make an
ordered cross block by itself symmetric or PSD.

## 6. Canonical operation order

Every plan has the following semantic order, with identity stages represented
explicitly in plan identity:

```text
ordered partition-pair product
  -> spatial normalization
  -> edge transform
  -> partition reduction
  -> pair query
```

Stages may be fused only when an algebraic law proves the same result. The
compiler must not silently move normalization, centering, ranking, Fisher
transformation, or a query across partition reduction.

The initial spatial normalizers are defined for positive finite spatial
weights `w_v`, mass `a = sum(w)`, left feature vectors `x_v`, and right feature
vectors `y_v`:

\[
T=\sum_v w_vx_vy_v^\top,\quad
u_L=\sum_v w_vx_v,\quad
u_R=\sum_v w_vy_v,
\]

\[
C=T-u_Lu_R^\top/a.
\]

- `inner_product()` returns `T`.
- `covariance()` returns the population-weight covariance `C / a`. Arbitrary
  frame weights do not imply a frequency-weight or unbiased degrees-of-freedom
  correction.
- `cosine()` divides `T[i,j]` by the corresponding uncentered weighted norms.
- `correlation()` divides `C[i,j]` by the corresponding centered weighted
  standard deviations.

For correlation, zero variance is never converted to `NA`. A declared
`zero_variance = "error"` policy fails and identifies the edge and axis
coordinate. A declared `zero_variance = "zero"` policy returns exactly zero
for every entry touching a zero-variance coordinate. The default is `"error"`.
The policy is semantic and enters plan identity.

`fisher_z()` is the elementwise `atanh` edge transform and is valid only for a
correlation-valued input. With the default `boundary = "error"`, any value with
absolute magnitude at least one fails. An explicit `boundary = "clip"` requires a finite
`delta` in `(0, 1)` and evaluates `atanh(pmin(1-delta,
pmax(-1+delta, r)))`. The boundary policy and `delta` enter plan identity.
There is no silent clipping and no infinite durable result.

`rank_edges()` ranks within each declared edge coordinate set before
partition reduction. The default `ties = "average"` assigns the mean occupied
rank. Other tie methods, if admitted, require distinct explicit values; a
stable-order method uses canonical coordinate order, never task completion
order. Missing values are not admitted by this contract. The tie policy enters
plan identity.

Correlation and Fisher transformation are nonlinear. In general,

```text
normalize(reduce(edges)) != reduce(normalize(edges))
fisher(reduce(correlations)) != reduce(fisher(correlations))
```

The unqualified estimand is edge-first: normalize, transform, then reduce. An
aggregate-first estimand, if later supported, requires a distinct public name,
plan identity, and receipt.

## 7. Reduction, centering, and queries

The partition reducer is a declared operation over ordered edge results. Its
weights are finite, and any normalization of those weights is explicit. A
linear weighted sum is the compatibility reducer. It does not imply that
partition edges are statistically independent.

The `independence` declaration on a pairing concerns the compatibility of its
endpoint relation estimates for cross-products. It does not declare the edge
products independent of one another. In an all-unordered-pairs estimator, two
edges sharing a partition also share estimation noise. Reducer rows are
contributions to the effect estimand, not sampling replicates for a
spread-across-edges standard error. Sampling covariance and calibration are
governed by `evidence-sampling-v1`.

The centered sufficient statistic is exact:

\[
C=T-u_Lu_R^\top/a
=\sum_vw_v(x_v-\bar x)(y_v-\bar y)^\top,
\]

with `bar(x) = u_L / a` and `bar(y) = u_R / a`. Rectangular coherent geometry
is `u_L u_R^T / a`; configuration is `T - coherent`. A singleton positive-mass
frame therefore has zero configuration. Implementations may retain
`T`, `u_L`, `u_R`, and `a` and need not materialize centered features.

Linear fusion is lawful. For fixed `H`,

\[
\left\langle H,\sum_ew_eP_e\right\rangle_F
=\sum_ew_e\,b_{L,e}^\top H b_{R,e}.
\]

This permits a final pair query to be pushed into an otherwise linear kernel.
It does not permit a query to cross a nonlinear normalization or transform.

A pair query cancels arbitrary additive left and right effect baselines exactly
only when both marginals of `H` vanish:

\[
H\mathbf 1_R=0,\qquad \mathbf 1_L^\top H=0.
\]

This is a property of the final query after all weighting and covariate
adjustment. A nominally double-centered matrix `Q` does not preserve both
properties under arbitrary weighting: `diag(a) Q` retains zero row sums but
generally loses zero column sums. Constructors must diagnose the final `H` and
must not overclaim baseline invariance.

## 8. Noise-unbiasedness of cross-partition products

Unbiasedness is a consequence of the pairing, not a property the analyst
confers. `independence = "independent"` records that the precondition below is
*asserted*; it never establishes it. This section states the precondition, what
follows from it, and where it breaks.

### Model

Write the relation estimated within partition `r` as

\[
\widehat B_r=B+\Xi_r,
\]

with `B` the common latent target, shared by every partition, and \(\Xi_r\) the
partition's estimation error. Assume

- **(A1)** \(\mathbb E[\Xi_r]=0\) for every partition contributing to the
  pairing;
- **(A2)** \(\mathbb E[\Xi_rM\Xi_s^\top]=0\) for every fixed feature-space `M`
  and every
  ordered pair with \(\Gamma_{rs}\neq0\) — the estimation errors at the two
  endpoints of each paired edge are uncorrelated. Independence of the two
  errors implies this; the theorem needs only the weaker statement; and
- **(A3)** the metric `K` and the frame weights are fixed: neither is a
  function of the data entering the products.

The pairing satisfies \(\Gamma_{rr}=0\) and \(\sum_{r\neq s}\Gamma_{rs}=1\).

### Theorem

Under (A1)–(A3),

\[
\mathbb E\!\left[\widehat G\right]
=\mathbb E\!\left[\sum_{r\neq s}\Gamma_{rs}\widehat B_rK\widehat B_s^\top\right]
=BKB^\top=G .
\]

*Proof.* Expand one edge:
\(\widehat B_rK\widehat B_s^\top
=BKB^\top+BK\Xi_s^\top+\Xi_rKB^\top+\Xi_rK\Xi_s^\top\).
The two first-order terms have zero expectation by (A1) and the fixity of `K`
and `B`; the second-order term has zero expectation by (A2), which applies
because \(r\neq s\) on every edge. Summing with weights that total one gives
`G`. No \(\Xi_rK\Xi_r^\top\) term — the one term whose expectation is a
noise variance rather than zero — is ever formed, which is why no
multiplicative bias correction is required. ∎

### Corollary: the components inherit it

For a frame row `w` with mass `a`, the coherent and configuration results are
the same cross-products read under two different **fixed** metrics,
\(ww^\top/a\) and \(D(w)-ww^\top/a\). Both depend only on the declared frame,
so (A3) holds for each, and the theorem applies to each separately:

\[
\mathbb E[\widehat G^{\mathrm{coherent}}]=G^{\mathrm{coherent}},\qquad
\mathbb E[\widehat G^{\mathrm{configuration}}]=G^{\mathrm{configuration}},
\]

and their sum is the unbiased total. The same argument covers every fixed pair
query, since \(\langle H,\cdot\rangle\) is linear:
\(\mathbb E[\langle H,\widehat G\rangle]=\langle H,G\rangle\). It does **not**
cover nonlinear readouts. `coherence_fraction` is a ratio of two unbiased
estimates and is not itself unbiased; that is one reason it is reported only
where the components form a nonnegative partition.

### Where it fails

- **Run-shared nuisance regressors.** A confound, baseline, or drift basis
  estimated jointly across partitions correlates \(\Xi_r\) with \(\Xi_s\), so (A2)
  fails and the residual bias is the \(\Gamma\)-weighted sum of
  \(\mathbb E[\Xi_rK\Xi_s^\top]\), which is not zero.
- **Temporal autocorrelation spanning partitions.** Noise correlated across a
  partition boundary breaks (A2) at exactly the edges the estimator uses.
  Partition boundaries must be chosen so that this correlation is negligible;
  the declaration does not make it so.
- **An extractor sharing `X` across partitions.** A design, filter, or
  whitening matrix fitted on pooled data and then applied per partition leaves
  a shared fitted object inside every \(\Xi_r\). Partition-wise estimation from
  disjoint data is the requirement.
- **A learned or estimated metric.** If `K` is estimated from the paired data,
  (A3) fails and \(\mathbb E[\widehat B_r\widehat K\widehat B_s^\top]\neq
  BKB^\top\) in general. Under a training-excluded policy, with \(\widehat K\)
  estimated from data disjoint from both endpoints, the theorem holds
  conditionally: \(\mathbb E[\widehat G\mid\widehat K]=B\widehat KB^\top\).
  That is unbiased for the estimand *defined by the realized metric*, not for a
  fixed-metric estimand, and the two differ whenever the metric is random.
- **Unequal partitions.** Unequal partition sizes do not break the theorem —
  it needs only weights summing to one — but the sampling-covariance laws in
  `evidence-sampling-v1` additionally assume the equal-partition model and are
  not licensed by this section.

### What this section does not claim

Unbiasedness of a point estimate is not a standard error, a calibration, or a
test. It also implies that individual estimates may be negative: negative
crossvalidated values are the visible cost of removing the noise term, and
clipping them at zero reintroduces exactly the bias the pairing removed. That
is why no clipping option exists.

## 9. Algebraic laws required of implementations

All production paths must preserve these laws within the declared numerical
contract:

1. **Self specialization:** the two half-weight ordered orientations equal the
   existing symmetrized undirected estimator.
2. **Transpose:** reversing sides, endpoints, bridge legs, and query axes
   transposes the complete form and preserves the scalar query value.
3. **Direct sum:** concatenating disjoint common-measurement coordinates adds
   their forms. This is the feature-blocking law.
4. **Packed/rectangular agreement:** `svec` and `vec` contractions agree for
   every symmetric self form and symmetric query.
5. **Linear fusion:** direct and late linear queries agree.
6. **Final-query baseline invariance:** it holds exactly when the final `H`
   has zero row and column marginals.
7. **Centered sufficient statistic:** `T - u_L u_R^T / a` equals explicit
   weighted feature centering.
8. **Bridge:** common-coordinate evaluation equals evaluation with
   `t(L_left) %*% L_right`; bridge reversal transposes.
9. **No bridge:** implicit and explicit identity legs agree for one exact
   neural identity, while distinct identities fail even at equal shape.

Negative fixtures are equally normative: correlation and Fisher
transformation do not commute with partition reduction, and covariate-weighted
`diag(a) Q` can lose column balance.

## 10. Identity and receipts

The scientific plan identity is a canonical hash of all semantic inputs:

- ordered left and right effect-space identities;
- ordered left and right relation-family identities;
- ordered partition endpoints, semantic edge weights, and edge selection;
- neural-space identities and the bridge identity or canonical same-space
  identity bridge;
- spatial frame and normalization, including zero-variance policy;
- edge transform and its boundary/tie parameters;
- partition reducer, weight convention, and stage order;
- result component and `complete_form` versus `query_only` capability; and
- the axis-bound query `H` for a query-only plan (or an explicit absent-query
  marker for complete materialization).

Changing any item above changes the scientific plan ID. Storage location,
tile sizes, worker count, progress reporting, and other execution choices do
not change it.

The execution receipt carries the scientific plan ID and separately identifies
the two neural/domain signatures, stable source revisions, bridge revision,
kernel version, task partition, reduction plan, storage codec, precision,
numerical contract, completion state, and observed resource facts. A codec or
execution change can therefore change receipt identity without changing the
scientific estimand. A receipt cannot repair or broaden the capabilities of its
result.

## 11. Public compatibility and deprecation

The lift is a conservative generalization:

- `materialize_geometry()` remains the public complete self-form constructor.
- `bilinear_query()` remains the symmetric same-space compatibility
  constructor and compiles to an axis-bound pair query.
- `cross_partitions()` retains its current unordered public meaning and exact
  numerical estimator through ordered-edge expansion.
- eligible self forms retain symmetric-packed memory and block storage.
- existing `effect_geometry` views remain available when their required
  capabilities are present.
- query-only execution continues to return `effect_view`, never a partial
  `effect_geometry`.

No existing call is silently reinterpreted as a rectangular, directed,
aggregate-first, normalized, bridged, or query-complete analysis. Any future
retirement of a compatibility constructor requires a documented replacement,
an ordinary deprecation cycle with warnings, and differential evidence first.
The old self kernel may be removed only after the new specialization passes the
full compatibility court.
