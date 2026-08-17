# Conservative geometry contract — detection versus attribution

Status: normative architecture contract

Contract version: `conservative-geometry-v1`

Date: 2026-08-17

This document freezes the semantics of *conservative* spatial frames: what
they estimate, what they conserve, what they do not conserve, and what a
population layer may assume about them. It is normative for WS-D and is the
interface WS-E consumes. It does not itself change production code or the
supported 0.1 API.

It sits below `design/effect-form-contract.md` (`effect-form-v1`), which
already fixes the bilinear-form algebra, the packed codec, the canonical
operation order, and the noise-unbiasedness theorem. Every statement here is
a statement about a *choice of frame normalization* made inside that algebra;
nothing here relaxes an `effect-form-v1` law.

Three oracle scripts under `design/oracles/` are independent of the production
contraction paths wherever a first-principles construction is possible; the
numbers quoted in this document are their measured output. Section 10 maps
every numbered claim to a test id or an oracle path.

Line numbers below are as of `elite-pass`, 2026-08-17; the accompanying
function name is the durable reference.

---

## 1. Two instruments: the detection map and the attribution map

A spatial frame is a nonnegative measurement-by-feature operator `W`, row `x`
holding the weight each domain feature contributes to node `x`. `effect-form-v1`
§6 defines the total component of a node as the weighted sum of ordered
per-feature products. Writing the two paired relation estimates as `B^(L)` and
`B^(R)` with feature columns `b_v`, the node total is

\[
G_x^{\mathrm{total}}=\sum_v w_{xv}\,b_v^{(L)}b_v^{(R)\top}.
\]

The *only* difference between the two instruments is which margin of `W` is
normalized.

### 1.1 `normalization = "local"` — the detection map

Rows sum to one: \(\sum_v w_{xv}=1\). Node `x`'s value is a weighted **average**
of per-feature products over its own support. Its estimand is

\[
\theta_x^{\mathrm{loc}}
=\Big\langle H,\ \textstyle\sum_v \tilde w_{xv}\,b_v^{(L)}b_v^{(R)\top}\Big\rangle,
\qquad \textstyle\sum_v\tilde w_{xv}=1,
\]

an intensive quantity: the mean evidence density inside the neighbourhood.
Because the normalizer is the node's own mass, values are directly comparable
across nodes of different size, and a node's value does not shrink merely
because it is large. This is the right instrument for the question *where is
there evidence?* — a detection map, read as a statistical image, thresholded
and reported as peaks. It is what `searchlights()`, `regions()` and
`whole_brain()` give by default.

A local frame's node values are **not** contributions to anything. Overlapping
neighbourhoods double-count shared features, so summing a detection map is
meaningless; `frame_conservation()` on a local frame returns
`conserved = FALSE` and reports the per-feature mass deviation
(`test-integrity-guards.R` "conservative frames conserve total evidence; local
frames do not").

### 1.2 `normalization = "conservative"` — the attribution map

Columns sum to one: \(\sum_x w_{xv}=1\), implemented as `W D(1/colSums W)`
(`.normalize_frame()`, `R/frame.R:425-440`). Every domain feature's unit of
evidence is *partitioned*
among the nodes that see it. Node `x`'s value is an extensive quantity: its
**share of one fixed global budget**. Its estimand is

\[
\theta_x^{\mathrm{cons}}=\big\langle H,\ G_x\big\rangle,
\qquad \sum_x\theta_x^{\mathrm{cons}}=\big\langle H,\ G_\Omega\big\rangle,
\]

where `G_Ω` is the geometry of the whole domain under the **unnormalized**
operator `whole_brain("none")`. This is the right instrument for the question
*how is the total effect distributed?* — an attribution map, read as a ledger,
reported as shares of a whole and aggregated over regions, networks or
subjects.

Nodes of a conservative frame are not comparable to one another as densities:
a large node holds more of the budget than a small one for the trivial reason
that it covers more features. Conversely, a detection map cannot be
aggregated. **Choosing between them is a scientific decision about the
question, not a preference about scaling, and the two must never be compared
node-for-node.**

`voxelwise()` defaults to `"conservative"` (a point frame is trivially both);
every other constructor defaults to `"local"`. Conservative must be asked for.

---

## 2. The conservation theorem

**Claim 2.** For a column-normalized frame, the total component conserves
exactly:

\[
\sum_x G_x^{\mathrm{total}}=G_\Omega^{\mathrm{total}},
\qquad\text{hence}\qquad
\sum_x\langle H,G_x\rangle=\langle H,G_\Omega\rangle\ \ \text{for every fixed }H.
\]

*Proof.* Exchange the sums:
\(\sum_x\sum_v w_{xv}b_v^{(L)}b_v^{(R)\top}
=\sum_v\big(\sum_x w_{xv}\big)b_v^{(L)}b_v^{(R)\top}
=\sum_v b_v^{(L)}b_v^{(R)\top}=G_\Omega\). The query statement follows from
linearity of \(\langle H,\cdot\rangle\) (`effect-form-v1` §7). ∎

Two preconditions are load-bearing and are the usual source of a spurious
"conservation failure":

1. **The comparator must be unnormalized.** `whole_brain()` defaults to
   `"local"`, which divides by the feature count; the comparator is
   `whole_brain("none")`.
2. **The law covers `total` only.** See §4.

Measured (`design/oracles/conservative-multiscale-ledger.R`): for point,
radius-1.01 and radius-2.01 conservative frames on a 12-feature line domain,
`max |Σ_x G_x − G_Ω|` is `0`, `1.78e-15` and `1.78e-15` respectively.

---

## 3. The conservative total is a smoothed univariate voxel ledger

**Claim 3.** Because the total is \(\sum_v w_{xv}b_v^{(L)}b_v^{(R)\top}\), the
per-voxel outer products form a **fixed ledger** that the frame merely
reallocates. The frame appears only as a nonnegative linear map applied *after*
the multivariate content has been fixed per voxel. Consequently the total
field at any scale is a spatially smoothed version of one univariate map (one
per query coordinate), and carries no cross-voxel information whatsoever.

This is visible in the implementation: `.effect_form_feature_task()`
(`R/task.R:123-190`) forms one atom row per feature — the ordered outer product
of that feature's two relation columns — and the frame enters only at the
spatial contraction `weight_slice %*% atom_slice` in the streamed kernel
(`R/kernel.R:354-377`). Nothing between the two steps mixes features.

Measured: reconstructing the ledger from first principles (no crossform
contraction) and contracting it with the frame reproduces
`geometry_component(., "total")` to `1.1e-16 – 2.2e-16` at every scale, and its
column sums equal `G_Ω` exactly (`0e+00`).

### 3.1 Corollary: per-scale totals under an α-weighted family are fixed by α

Let \(W_1,\dots,W_S\) each be column-normalized and let
\(\alpha_1,\dots,\alpha_S\ge 0\) with \(\sum_s\alpha_s=1\). Stacking
\(\alpha_s W_s\) gives a frame whose columns still sum to one, so it is
conservative — and, block by block,

\[
\sum_{x\in\text{scale }s}G_{s,x}=\alpha_s\,G_\Omega .
\]

Measured (three scales, `α = (0.2, 0.5, 0.3)`, 12-feature domain): worst
absolute deviation `8.88e-16`, worst relative deviation `3.43e-16`. For the
query `H = I`, with `⟨H, G_Ω⟩ = 4.077044941591`:

| scale | α_s | E_s = Σ_x ⟨H, G_{s,x}⟩ | α_s·⟨H, G_Ω⟩ |
|---|---|---|---|
| point | 0.2 | 0.815408988318 | 0.815408988318 |
| searchlight (r = 1.01) | 0.5 | 2.038522470795 | 2.038522470795 |
| wide (r = 2.01) | 0.3 | 1.223113482477 | 1.223113482477 |

**Normative consequence.** A multiscale panel plotting *total energy per scale*
is a plot of the analyst's own α vector. It is not a finding, and no such
panel may be presented as evidence about spatial scale. Any flagship figure
built on the conservative frame family must make this explicit.

### 3.2 The informative object is the coherence spectrum

What is not fixed by α is the **split** of each scale's fixed budget into
coherent and configuration parts. The coherent share

\[
\phi_s=\frac{\sum_{x\in s}\langle H,G^{\mathrm{coh}}_x\rangle}
             {\sum_{x\in s}\langle H,G_x\rangle}
\]

depends on the data through cross-voxel products and is the only scale-resolved
quantity in the table above that is not proportional to α. The same holds
location-wise: the coherent share as a function of node is a genuine map.

Measured on the same fixture:

| scale | E_s (fixed by α) | E_s^coh | coherent share φ_s |
|---|---|---|---|
| point | 0.815408988 | 0.815408988 | 1.000000 |
| searchlight (r = 1.01) | 2.038522471 | 1.186543279 | 0.582060 |
| wide (r = 2.01) | 1.223113482 | 0.301364694 | 0.246391 |

The left column is α times a constant; the right column is not. (The point
scale reads 1.000000 because a singleton positive-mass frame has exactly zero
configuration — `effect-form-v1` §7, `test-effect-form-laws.R` "a singleton
positive-mass frame has exactly zero configuration".)

**Normative:** the conservative frame family's scientific product is the
*coherence spectrum* — coherent share versus scale and versus location — not a
scale-resolved energy budget.

---

## 4. Coherent does not conserve, and coherent budgets are frame-relative

**Claim 4.** \(\sum_x G_x^{\mathrm{coh}}\neq G_\Omega^{\mathrm{coh}}\), for any
frame.

The coherent component of a node is the rank-one part of the node's form along
that node's **own** weighted common mode: the metric-general construction is
\(K_{\mathrm{coh}}=aa^\top/(a^\top K^{-1}a)\) with `a` the node's normalized
frame weights (`metric_components()`, `R/metric.R:915-996`; the identity-metric specialization
`(Bw)(Bw)^\top/\sum w` is what `src/packed-form.cpp:144-222` computes from
first moments). Different nodes have different `a`; there is no linear operator
that turns a collection of node-local rank-one projections into the global one.
Nothing in the construction could make them sum.

Measured: with `α`-weighted scales, `E_s^coh` above differs from
`α_s·⟨H, coh_Ω⟩ = (0.235591844, 0.588979610, 0.353387766)` in every row. Per
scale, `max |Σ_x G^coh_x − G^coh_Ω|` is `7.41`, `2.70` and `1.27` on a fixture
where `max |G_Ω|` is order 1.

**Normative consequences.**

1. `frame_conservation()` reports `component = "total"` and must never be read
   as covering the components. Its documentation already says so
   (`frame_conservation()`, `R/frame.R:355-423`).
2. \(\sum_x G^{\mathrm{coh}}_x\) is **not a global quantity**. Any fraction,
   share or percentage computed against it is a fraction of *that frame's*
   coherent mass and must be labelled with the frame identity. Two frames
   produce two incomparable denominators.
3. A coherence fraction may be reported only where the components form a
   nonnegative partition. `crossform` already masks rather than clamps
   (`R/result.R:859-863`): `coherence_fraction` is `NA` unless
   `total > 0 && coherent >= 0 && configuration >= 0`, and the mask ships
   alongside as `coherence_fraction_valid`.

---

## 5. The metric composition law, and the `composition =` choice

### 5.1 The native composition

`crossform` composes a node metric with the frame by symmetric congruence in
the square-root weights (`.compose_frame_metric()`, `R/metric.R:814-895`, provenance
`law = "D(sqrt(w)) K D(sqrt(w))"`):

\[
K_x^{\mathrm{native}}=D(\sqrt{w_x})\,Q\,D(\sqrt{w_x}).
\]

Summing over a conservative frame gives

\[
\sum_x K_x^{\mathrm{native}}=S\circ Q,
\qquad S_{uv}=\sum_x\sqrt{w_{xu}w_{xv}},
\]

an entrywise (Hadamard) reweighting of `Q`. Column normalization forces
\(S_{uu}=\sum_x w_{xu}=1\) but says nothing about \(S_{uv}\) for `u ≠ v`.

**Claim 5.**

- **Identity and diagonal metrics conserve.** If `Q = D(q)` then only the
  diagonal of `S` is used, and \(\sum_x K_x = D(q) = Q\). The implementation
  takes the exact diagonal specialization `diag(w * diag(K))`
  (`R/metric.R:831-836`), and `.metric_additive_frame()` folds it into the
  frame (`R/metric.R:590-637`).
- **Dense fixed metrics do not conserve.** The off-diagonal `S_uv ≠ 1`
  survives into the sum. `.metric_frame_conservation()`
  (`R/metric.R:998-1083`) already refuses to certify a non-diagonal schedule:
  `feature_additive = FALSE`, `global_metric_kind = "support_pair_operator"`,
  reason "Conservative frame weights alone do not conserve a non-diagonal
  metric schedule".
- **A learned local precision breaks it by construction**, and additionally
  breaks (A3) of `effect-form-v1` §8, so its estimand is metric-conditional in
  the first place.

Measured (`design/oracles/conservative-metric-composition.R`, 9-feature line
domain, radius-1.01 conservative searchlights, `Q = crossprod(A)/9 + I`):

| composition | max relative deviation of Σ_x G_x from G_Ω(Q) | signed trace error |
|---|---|---|
| implicit identity | `1.11e-16` | `−0.0000%` |
| diagonal `D(q)`, `q ~ U(0.5, 2.5)` | `2.63e-16` | `−0.0000%` |
| dense `Q`, native | `2.38e-01` | `+21.2009%` |
| dense `Q`, whitened (§5.2) | `1.53e-15` | `−0.0000%` |

The algebraic law is verified against the executed package result to
`8.88e-16`: `Σ_x G_x = B(S∘Q)B^⊤` exactly, with `max |S_uu − 1| = 2.2e-16` and
off-diagonal `S_uv` spanning `[0.0000, 0.8165]`.

**Correction to `.planning/2026-08-17-feedback-assessment.md` §Part 2(3).** That
document reports the dense-metric failure as `−6.6%`. That number is one draw,
not a property. Sweeping twelve draws of `crossprod(A)/p + I` on the same frame
gives signed trace errors of
`+15.70, +38.42, +155.36, −34.47, −23.58, −20.73, −0.85, −17.41, +60.83, +188.35, +3.06, −11.37`
percent — range `[−34.5%, +188.4%]`, both signs, and one draw within 1% of zero.
**The contract claim is the algebraic law, not a percentage**, and in particular
a small measured deviation on one dataset is not evidence that a dense metric
conserves.

### 5.2 The whitened composition, and why it must be an explicit choice

The alternative composition places the frame in whitened coordinates:

\[
K_x^{\mathrm{wh}}=Q^{1/2}D(w_x)Q^{1/2},
\qquad
\sum_x K_x^{\mathrm{wh}}=Q^{1/2}D\Big(\sum_x w_x\Big)Q^{1/2}=Q^{1/2}IQ^{1/2}=Q .
\]

This conserves exactly, for any SPD `Q`. Measured: `1.53e-15` relative.

It is **a different estimand**, not a bug fix. Under the native composition a
node weights *features* and then measures them in the `Q` geometry; under the
whitened composition it weights *whitened coordinates*, which are spatially
delocalized whenever `Q` is dense — a node's support is no longer its support.
Measured per-node disagreement on the fixture above: `6.01e-01` against
`max |native| = 4.46e+00`, i.e. ≈13.5% of the largest node value.

**Normative:** a conservative crossnobis-style analysis must expose
`composition = c("native", "whitened")`, defaulting to `"native"` (the current
behaviour), carrying the choice into the scientific plan identity, and emitting
a conservation certificate. The switch must never be applied silently to
"repair" a failed conservation check. `crossform` has no `"whitened"` code path
today; §10 records the oracle that defines what it must compute.

### 5.3 The diagonal-metric fold, and what it must keep

`.metric_additive_frame()` folds a native-diagonal metric into the frame
weights and rebuilds with `normalization = "none"` (`R/metric.R:590-637`).
That declaration is correct rather than lossy: after the fold the columns sum
to `q_v`, not to one, so declaring `"conservative"` would make the frame
validator refuse a numerically correct operator.

`.planning/2026-08-17-feedback-assessment.md` §Part 2(3) recorded this as a
latent provenance inconsistency ("discarding the declared normalization from
provenance even though the numbers survive"). **That has since been fixed on
`elite-pass`** and the fix is now part of this contract's baseline: the folded
frame carries `$metric_folded` with `declared_normalization`, `metric_kind`,
`metric_signature`, `composition = "diagonal_metric_fold"` and a
`reference_mass` equal to the folded metric diagonal, and
`frame_conservation()` certifies against `reference_mass` rather than against
one, additionally returning `declared_normalization` and `metric_folded`
(`R/frame.R:396-423`).

Measured: for a conservative radius-1.01 frame folded with
`q ~ U(0.5, 2.5)`, `normalization` reads `"none"`,
`$metric_folded$declared_normalization` reads `"conservative"`,
`|reference_mass − q| = 0`, and `|colSums(folded) − q| = 0`.

**What remains required.** The fold still returns a frame without `$index` or
`$specification`, so node labels do not survive it. That is the same gap §7.1
records for the family route, and the family constructor must close both.

---

## 6. The signed estimation layer and the latent PSD descriptive layer

**Claim 6.** Crossvalidated contributions are **signed**, and the arithmetic
that treats a node's value as part of a nonnegative whole is invalid on that
layer.

`effect-form-v1` §8 establishes that the cross-partition estimator is unbiased
precisely because no \(\Xi_r K\Xi_r^\top\) term is ever formed, and that
negative values are the visible cost of removing the noise term. Conservation
(§2) is a statement about a signed sum. It does not make the summands
nonnegative.

Measured (pure-noise fixture, 8 nodes): per-node totals
`+1.81 +1.76 +0.11 +0.59 −0.65 −0.05 −0.92 +0.37`, three of eight negative, with
conserved budget `+3.0246`. The implied "shares" `x / Σx` span `[−0.304,
0.599]` — outside `[0,1]`. Clipping to nonnegativity raises the total from
`+3.0246` to `+4.6413`, an inflation of `+53.45%`, which both reintroduces the
bias the pairing removes and destroys conservation.

**Two layers, named separately.**

| | signed estimation layer | latent PSD descriptive layer |
|---|---|---|
| what it holds | the unbiased crossvalidated estimates | a declared nonnegative projection of them |
| unbiased? | yes (`effect-form-v1` §8) | no — clipping is a bias |
| conserves? | yes (§2) | no |
| admits fractions, cumulative curves, `n_eff`, entropy, "top-k explains X%" | **no** | yes, *within the layer* |
| how it is obtained | the estimator | an explicit, recorded projection |

**Normative:**

1. Fractions, cumulative-contribution curves, effective-node counts and any
   other functional requiring a nonnegative partition live **only** on the
   latent layer, and are labelled with the projection that produced them.
2. PSD clipping, eigenvalue truncation, or any nonnegativity projection is
   **never silent**: it is a named operation, it enters plan identity, and the
   receipt records how much mass it moved.
3. Coherence fraction is already handled correctly by masking
   (`R/result.R:859-863`); that discipline is the template, not an exception.
4. Transport and OLS are linear operators on the estimates, so unbiasedness
   survives them by the same linearity argument `effect-form-v1` §8 uses for
   fixed queries and components. A nonlinear latent projection applied
   *before* aggregation does not survive, and is therefore forbidden upstream
   of a population fit.

---

## 7. Transport readiness: what a population layer needs

**Claim 7.** WS-E can consume a conservative result if and only if the
following seven items are supplied. Items marked **missing** do not exist
today and are WS-D deliverables.

### 7.1 Per-row frame metadata — **missing**

A transport operator maps *native nodes* to *group nodes*. Every row of a
conservative frame must therefore be self-describing:

| field | meaning | today |
|---|---|---|
| `family` | which frame in a multi-frame family the row came from | absent |
| `scale` | the row's own scale parameter (e.g. searchlight radius) | frame-wide on `$specification$radius`, not per row |
| `center` | the row's anchor feature or coordinate | recoverable only for `searchlights` via `$support_index` |
| `label` | the row's scientific identity | `$index$measurement`, present only on **compiled** frames |
| `alpha` | the row's family weight, when stacked | absent |

Measured: a compiled `searchlights(1.01, "conservative")` frame has
`$index` columns `measurement` only and `$specification` fields
`kind, normalization, radius`. Node labels *do* reach the result
(`as.data.frame(view)$measurement` reads `vox1, vox2, vox3`). But a two-scale
family built the only way that works today — `rbind()` the weight matrices and
call `additive_frame(..., "conservative")` — conserves correctly
(`frame_conservation()$conserved == TRUE`, 16 rows) while carrying **neither**
`$index` nor `$specification`, so its result's `measurement` column degrades to
positions `1, 2, 3, 4, …`. The family route is numerically sound and
provenance-blind. This is the single largest gap between today's frame and a
transportable one.

### 7.2 A Frobenius-consistent packed codec — **exists**

The subject-Gram trick requires that Euclidean geometry on stored rows *be*
Frobenius geometry on forms. `symmetric_packed` (`effect-form-v1` §3; √2 on
off-diagonals, `crossform:::.svec_symmetric`) satisfies

\[
\langle\operatorname{svec}(A),\operatorname{svec}(B)\rangle=\langle A,B\rangle_F ,
\]

measured to `1.78e-15` worst relative over 200 random symmetric pairs, at
storage width `q(q+1)/2` instead of `q²`.

### 7.3 The subject-Gram eigen trick — **exists as a consequence**

Geometry-space covariance over `N` subjects has rank at most `N−1`, so its
nonzero spectrum is recoverable from the `N × N` subject Gram. Measured with
`N = 6`, `P = 15`: `rank(Σ) = 5 = N−1`; the top 5 eigenvalues of the `P × P`
covariance and of the `N × N` Gram agree to `2.67e-15`; the remaining `P−N+1`
eigenvalues are `≤ 8.14e-16`. For realistic `q`, this is the difference between
a feasible and an infeasible population covariance.

### 7.4 Row-stochastic transport with a required sink node — **missing**

A transport `P` with `P ≥ 0` and unit row sums preserves each subject's budget
exactly: measured `|Σ(Pᵀc) − Σc| = 0`. The sink is **not optional**. Measured:
dropping three of ten native nodes from the atlas — a perfectly ordinary
partial-coverage situation — silently loses `30.0%` of the subject's budget,
because the rows are no longer stochastic. With an explicit sink column the
budget closes to `0e+00` and the unmapped mass is *visible* in the sink rather
than absent from the ledger. Conservation is what makes the loss
detectable at all: a detection map has no budget to check against, so the same
omission leaves no numerical trace.

### 7.5 Budget versus density semantics — **missing, and must be declared**

`P` preserves budgets, so a group node receives mass in proportion to how many
native nodes map to it. Two subjects with the same total evidence but different
native frame resolutions therefore contribute unequally per territory. Measured
with four group nodes and two subjects each carrying budget 1:

| subject | native nodes | budget into group nodes | density into group nodes |
|---|---|---|---|
| fine | 10 | `0.20 0.30 0.20 0.30` | `0.10 0.10 0.10 0.10` |
| coarse | 4 | `0.25 0.25 0.25 0.25` | `0.25 0.25 0.25 0.25` |

Under budget semantics the fine subject puts `0.30` into group node 2 and the
coarse one `0.25`, purely because its native frame is finer. Neither column is
correct in general; the choice is a declared field on the transport object and
enters plan identity.

### 7.6 Labelling of transported components — **normative**

Transported coherent and configuration values are a **ledger of native-node
coherence carried to a group location**. They are not the coherence of the
group node's own geometry: no group frame defines a common mode there, and §4
says node-local rank-one projections do not sum. Result labels and print
methods must say so.

### 7.7 What crossform must not be asked to supply

Registration and functional-transport *learning* stay outside the package. WS-E
accepts a typed, sparse, provenance-bearing transport as an input. A transport
efficiency `η_transport` must be cross-fitted on independent runs or tasks, may
be negative on held-out data, and is reported as-is — clipping it would repeat
the §6 error at the population level.

---

## 8. What exists today, and what this contract requires

| item | today | this contract requires |
|---|---|---|
| `normalization = "conservative"` | `W D(1/colSums W)` in `.normalize_frame()` (`R/frame.R:425-440`); checked, not applied, by `additive_frame()` (`R/scope.R:430-436`) | unchanged |
| conservation diagnostic | `frame_conservation()`: `component = "total"`, `normalization`, `declared_normalization`, `metric_folded`, `feature_mass`, `reference_mass` | unchanged; keep the `total`-only scope explicit |
| default normalization | `voxelwise()` conservative; `searchlights()`, `regions()`, `whole_brain()` local | unchanged; documented as detection-by-default |
| frame **families** (multi-scale, α-weighted) | no API. `rbind(...)` + `additive_frame(..., "conservative")` works end to end but is undocumented | a named constructor carrying per-row `family`, `scale`, `center`, `label`, `alpha` |
| per-row metadata | `$index$measurement` on compiled frames only; radius frame-wide on `$specification` | per-row metadata that survives stacking |
| `additive_frame()` | drops `$index` and `$specification` | a family constructor that does not |
| `contribution(by = region/network)` | does not exist | aggregation over a declared grouping of rows, budget semantics |
| coherence spectrum | derivable from `contrast_energy()` per frame; no named object | a named per-scale / per-location coherent-share object |
| metric composition | native `D(√w) Q D(√w)` only | explicit `composition = c("native", "whitened")` in plan identity, with a conservation certificate |
| metric conservation certificate | `.metric_frame_conservation()` / `.require_metric_conservation()`, internal | surfaced wherever a conservative claim is printed |
| diagonal-metric frame fold | folds to `normalization = "none"` and records `$metric_folded` provenance (`R/metric.R:590-637`); the assessment's "discards provenance" finding is **fixed** | also carry `$index` / `$specification` through the fold |
| latent PSD layer | does not exist; `coherence_fraction` masking is the only instance of the discipline | a named layer object; no silent clipping anywhere |
| packed codec | `symmetric_packed`, √2 off-diagonals, Frobenius consistent | unchanged; it is the population storage format |
| transport / population | does not exist (only `design/archive/searchlight-conversation-ledger.md:982-1053`) | WS-E, consuming the above |
| feature-space bridge | `measurement_bridge()` between named spaces | unchanged; it is *not* node transport and must not be confused with it |
| sampling covariance | `sampling_covariance(..., operation = "transport")` transports the covariance through an output-by-evidence linear map | unchanged; no named `contrast_energy` sampling route yet |

---

## 9. Numerical contract

All statements above hold to the following measured tolerances on the oracle
fixtures. They are the tolerances tests must assert, not aspirational figures.

| law | tolerance |
|---|---|
| `Σ_x G_x = G_Ω`, identity or diagonal metric | `1e-12` absolute (measured `≤ 1.8e-15`) |
| per-scale `Σ_{x∈s} G_{s,x} = α_s G_Ω` | `1e-12` absolute (measured `8.9e-16`) |
| total equals the frame contraction of the voxel ledger | `1e-12` (measured `2.2e-16`) |
| `Σ_x K_x^wh = Q` (whitened composition) | `1e-12` relative (measured `1.5e-15`) |
| `⟨svec A, svec B⟩ = ⟨A,B⟩_F` | `1e-12` relative (measured `1.8e-15`) |
| subject-Gram spectrum vs `P×P` covariance spectrum | `1e-12` absolute (measured `2.7e-15`) |
| dense-metric non-conservation | asserted as `> 1%` relative, **never** as an equality |

---

## 10. Test and oracle index

| claim | statement | evidence |
|---|---|---|
| 1 | local = detection (intensive, row-normalized); conservative = attribution (extensive, column-normalized); non-conservation of local frames | `tests/testthat/test-integrity-guards.R` — "conservative frames conserve total evidence; local frames do not" (`frame_conservation()` TRUE vs FALSE); `tests/testthat/test-frame.R` — "conservative frames partition global feature mass" |
| 2 | `Σ_x G_x = G_Ω` for a column-normalized frame | `tests/testthat/test-workflow.R` — "public frame laws preserve point decomposition and global total" (1e-13); `tests/testthat/test-integrity-guards.R` — "conservative frames conserve total evidence; local frames do not" (query level, 1e-10); `design/oracles/conservative-multiscale-ledger.R` §O1.a |
| 3a | total is the frame contraction of a fixed per-voxel ledger | `tests/testthat/test-conservative-geometry-contract.R` — "a conservative total is the frame contraction of a voxel ledger"; `design/oracles/conservative-multiscale-ledger.R` §O1.b |
| 3b | per-scale totals under an α-weighted family equal `α_s G_Ω` | `tests/testthat/test-conservative-geometry-contract.R` — "an alpha-weighted frame family conserves total scale by scale"; `design/oracles/conservative-multiscale-ledger.R` §O1.c |
| 3c | the coherence spectrum, not the energy spectrum, is the informative object | `design/oracles/conservative-multiscale-ledger.R` §O1.d (share column varies; `E_s` column does not) |
| 4a | coherent does not conserve | `tests/testthat/test-integrity-guards.R` — "conservative frames conserve total evidence; local frames do not" (`expect_gt(abs(Σ coherent − global coherent), 1e-8)`); `design/oracles/conservative-multiscale-ledger.R` §O1.d |
| 4b | coherence fraction is masked, never clamped | `tests/testthat/test-views.R` — "contrast returns one exact decomposition and signed marginal"; `R/result.R:859-863` |
| 4c | singleton frames have zero configuration (why the point scale reads share 1) | `tests/testthat/test-effect-form-laws.R` — "a singleton positive-mass frame has exactly zero configuration" |
| 5a | identity and diagonal metrics conserve; dense metrics do not | `tests/testthat/test-conservative-geometry-contract.R` — "conservation survives a diagonal metric and fails for a dense one"; `design/oracles/conservative-metric-composition.R` §O2.a–c |
| 5b | the failure obeys `Σ_x G_x = B(S∘Q)Bᵀ` with `S_uu = 1` | `tests/testthat/test-conservative-geometry-contract.R` — "conservation survives a diagonal metric and fails for a dense one" (its `overlap` assertions); `design/oracles/conservative-metric-composition.R` §O2.c |
| 5c | the assessment's `−6.6%` is fixture-specific; sign and size both vary | `design/oracles/conservative-metric-composition.R` §O2.c′ |
| 5d | the whitened composition `Q^{1/2} D(w) Q^{1/2}` conserves and is a different estimand | `tests/testthat/test-conservative-geometry-contract.R` — "the whitened composition conserves where the native one does not"; `design/oracles/conservative-metric-composition.R` §O2.d |
| 5e | the package refuses to certify a non-diagonal schedule | `tests/testthat/test-metric.R` — "conservative identity conservation is capability-gated"; `design/oracles/conservative-metric-composition.R` §O2.e |
| 5f | the diagonal-metric fold declares `"none"` but keeps `declared_normalization` in `$metric_folded`, and `frame_conservation()` certifies against `reference_mass` | `design/oracles/conservative-metric-composition.R` §O2.f; `R/frame.R:396-423`, `R/metric.R:590-637`. Node labels still do not survive the fold — see §7.1 |
| 6 | contributions are signed; shares and clipping are invalid on the estimation layer | `tests/testthat/test-conservative-geometry-contract.R` — "contribution shares are undefined on the signed estimation layer"; `design/oracles/conservative-transport-readiness.R` §O3.e; `design/effect-form-contract.md` §8 |
| 7a | `symmetric_packed` is Frobenius consistent | `tests/testthat/test-effect-form-laws.R` — "the packed symmetric codec equals the oracle svec" and "querying a form equals the oracle Frobenius pairing"; `design/oracles/conservative-transport-readiness.R` §O3.a |
| 7b | geometry covariance has rank ≤ N−1; the subject Gram recovers its spectrum | `design/oracles/conservative-transport-readiness.R` §O3.b |
| 7c | row-stochastic transport preserves budgets; a missing sink silently loses mass | `design/oracles/conservative-transport-readiness.R` §O3.c |
| 7d | budget and density semantics differ; neither is a default | `design/oracles/conservative-transport-readiness.R` §O3.d |
| 7e | frame families are provenance-blind today | `design/oracles/conservative-transport-readiness.R` §O3.f |

Oracle scripts are plain R and run standalone:

```sh
Rscript design/oracles/conservative-multiscale-ledger.R
Rscript design/oracles/conservative-metric-composition.R
Rscript design/oracles/conservative-transport-readiness.R
```
