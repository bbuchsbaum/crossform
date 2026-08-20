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

**Reviewer note (2026-08-20).** All source line numbers in this document were
re-checked against commit `7ec6e37` (the contract's own commit) and resolve
correctly there. They are already stale against the working tree: `R/metric.R`
has since grown by 72 lines ahead of line 590, so every `R/metric.R:NNN`
citation below is off by `+72` in an uncommitted checkout (`.metric_additive_frame`
590→662, `.compose_frame_metric` 814→886, `metric_components` 915→987,
`.metric_frame_conservation` 998→1070). The function names resolve in both. The
line numbers have deliberately **not** been rewritten to the working tree, which
is mid-edit by other work streams; readers should resolve by name.

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

**Claim 3** (identity or diagonal metric only — see the precondition below).
Because the total is \(\sum_v w_{xv}b_v^{(L)}b_v^{(R)\top}\), the
per-voxel outer products form a **fixed ledger** that the frame merely
reallocates. The frame appears only as a nonnegative linear map applied *after*
the multivariate content has been fixed per voxel. Consequently the total
field at any scale is a spatially smoothed version of one univariate map (one
per query coordinate), and carries no cross-voxel information whatsoever.

**Metric precondition (load-bearing).** The expression
\(G_x=\sum_v w_{xv}b_v^{(L)}b_v^{(R)\top}\) is the *feature-additive* branch of
the native composition of §5.1: it is \(BD(\sqrt{w_x})QD(\sqrt{w_x})B^\top\)
specialized to \(Q=I\) or \(Q=D(q)\), where \(D(\sqrt{w})QD(\sqrt{w})=D(wq)\)
stays diagonal. For a **dense** \(Q\) the node form is
\[
G_x=\sum_{u,v}\sqrt{w_{xu}w_{xv}}\;Q_{uv}\,b_u^{(L)}b_v^{(R)\top},
\]
which pairs *distinct* voxels inside the node. There is then no per-voxel
ledger, the total is genuinely multivariate across space, and neither the
smoothing picture nor §3.1 applies. Claim 3 is therefore scoped to schedules
for which `metric_capabilities()$feature_additive` is `TRUE`; §5 is where the
general metric enters, and it is exactly the case where the ledger dissolves
and conservation fails together. Both the oracle and the cited test exercise the
identity-metric branch, so the evidence matches the scoped claim.

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

Two preconditions are load-bearing here and must be enforced by the family
constructor (D2/D3), not assumed:

1. **Every scale must be column-normalized *on its own*.** The per-scale
   identity needs \(\sum_{x\in s}w^{(s)}_{xv}=1\) for every feature `v`
   separately, which is strictly stronger than the stacked frame conserving.
   A scale that leaves any feature uncovered cannot be column-normalized at
   all (`.normalize_frame()` refuses), so the constructor must reject a family
   member whose support does not cover the domain rather than renormalizing the
   stack as a whole — the stack can conserve while no individual block does.
2. **α must sum to one**, or the block identity reads \(\alpha_s G_\Omega\)
   against a family that does not conserve. Whether the constructor
   renormalizes α, or refuses, is a D2 decision; this contract requires only
   that it be one of the two and that the applied α be recorded per row
   (§7.1).

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

**Reviewer note (2026-08-20): the statement is stronger than "not proportional
to α" — φ_s is *exactly invariant* to α_s, and that is what makes the coherence
spectrum well-posed.** Both components are homogeneous of degree one under a
rescaling \(w_x\mapsto\alpha w_x\) of a single row: the total is linear in `w`;
and the coherent part is \(K_{\mathrm{coh}}=aa^\top/(a^\top K_x^{-1}a)\) with
\(a=w_x/\!\sum w_x\) *invariant* to α while \(K_x^{-1}\) scales as \(1/\alpha\),
so \(K_{\mathrm{coh}}\) scales as α (identity-metric form
\((Bw)(Bw)^\top/\!\sum w\): manifestly degree one). The ratio therefore cancels
α exactly. Consequence for D5: the coherence spectrum is a property of the
frame family's *geometry*, not of the analyst's α vector, so it may be reported
without disclosing α — the exact opposite of the energy panel §3.1 forbids.

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
whose global coherent component has `max |G^coh_Ω| = 1.21`. The failure is
therefore several times the size of the quantity it is meant to reproduce, not
a tolerance effect.

*(Corrected 2026-08-20. This paragraph previously compared the coherent
deviation to `max |G_Ω|` and called it "order 1"; measured, `max |G_Ω| = 8.62`
on that fixture, so the comparison was both numerically wrong and against the
wrong object. The comparator for a coherent deviation is the global coherent
component.)*

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
   (`.new_effect_contrast_view()`, `R/result.R:859-863`): `coherence_fraction`
   is `NA` unless
   `is.finite(total) && total > 0 && coherent >= 0 && configuration >= 0`, and
   the mask ships alongside as `coherence_fraction_valid`. (The `is.finite`
   conjunct was missing from this contract's first draft; it is part of the
   guard.)

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

#### 5.2.1 The root is part of the estimand — pin it

**Claim 5.2.1 (added 2026-08-20).** The conservation argument above never uses
symmetry of the root. For *any* `R` with \(RR^\top=Q\),

\[
\sum_x R\,D(w_x)\,R^\top=R\Big(\sum_x D(w_x)\Big)R^\top=RIR^\top=Q ,
\]

so conservation is root-invariant. **The per-node values are not.** Measured
(`conservative-metric-composition.R` §O2.d′, the same 9-feature fixture as §5):
the symmetric PSD root and the lower Cholesky factor both conserve —
`1.53e-15` and `2.86e-16` relative — while their per-node forms differ by
`6.33e-01` against `max |symmetric| = 4.03e+00`, i.e. **15.7% of the largest
node value**, comparable to the native-versus-whitened gap of 13.5% the section
opens with. An independent draw during review gave `28.6%`; as in §5, the
invariant claim is the algebraic one, and the percentage is fixture-specific.

Writing the law as \(Q^{1/2}D(w_x)Q^{1/2}\) with the same symbol on both sides
silently selects the symmetric PSD root; the general form is
\(R\,D(w_x)\,R^\top\), and "whitened" alone does not name an estimand.

**Normative.** `composition = "whitened"` means the **symmetric positive
semidefinite root** \(Q^{1/2}\) — the choice the oracle
(`conservative-metric-composition.R` §O2.d) and the test
(`test-conservative-geometry-contract.R` "the whitened composition conserves
where the native one does not") both already take. The root identity must be
recorded in the plan identity alongside the composition, so that a future
Cholesky or ZCA variant is a distinguishable estimand rather than a silent
substitution. A conservation certificate is **not** sufficient evidence that
two whitened analyses computed the same thing.

#### 5.2.2 Whitening is not a support-local operation

A second consequence, for implementation rather than for semantics:
\(Q^{1/2}D(w_x)Q^{1/2}\) is dense on the whole domain even though `D(w_x)` is
supported on the node. It therefore cannot be produced through
`.compose_frame_metric()`, whose contract requires
`identical(metric$support, node_value$support)` (`R/metric.R:814-895`, HEAD) and
which returns an operator on the node's own support. The feasible route is to
whiten the relation once — \(\tilde B=BQ^{1/2}\) — and then run the ordinary
identity-metric conservative pipeline on \(\tilde B\), which reproduces
\(\tilde B D(w_x)\tilde B^\top=BQ^{1/2}D(w_x)Q^{1/2}B^\top\) exactly. D6 should
budget for a relation-level transform, not a metric-schedule variant.

**Normative:** a conservative crossnobis-style analysis must expose
`composition = c("native", "whitened")`, defaulting to `"native"` (the current
behaviour), carrying the choice **and the root identity of §5.2.1** into the
scientific plan identity, and emitting a conservation certificate. The switch
must never be applied silently to "repair" a failed conservation check.
~~`crossform` has no `"whitened"` code path today; §10 records the oracle that
defines what it must compute.~~ *(True until D6; see immediately below.)*

**Delivered (D6, 2026-08-20).** `plan_geometry(..., metric = , composition =
c("native", "whitened"))` implements this section. `"native"` is the default
and no existing plan identity moved. `"whitened"` is admitted for a fixed
positive-definite domain-wide metric only — a learned recipe gets a
`whitened_metric_composition` capability refusal, because a per-support
operator has no single global root — and it takes the seam §5.2.2 prescribed: a
relation-level transform. `plan_geometry()` forms `B̃ = BQ^{1/2}` once at plan
time with the symmetric root of §5.2.1, compiles the task on `B̃`, and attaches
a `whitened_metric_before_frame` schedule that lowers to the ordinary
implicit-identity `additive_contraction`, so nothing in the compiler, executor,
or kernel branches on the composition. `composition` and the root convention
string `"symmetric_psd_root"` enter the schedule's semantic digest and so
`$scientific_plan_id`, alongside the metric signature; they are absent from the
other schedule kinds rather than carried as `"native"`, which is what keeps
every pre-D6 identity byte-identical. `.metric_frame_conservation(frame,
composition = "whitened", metric = )` certifies the law by its exact algebraic
residual `Σ_x K_x − Q = Q^{1/2}D(m−1)Q^{1/2}` (`m` the column mass), and the
native refusal now names the alternative without implying it is a repair.
Evidence: `tests/testthat/test-composition.R`; design note
`design/crossform-execution-design.md` §16.8.

**One consequence left open by D6.** Because the plan's relation *is* the
whitened one, the residual-channel sampling routes
(`rdm_sampling_covariance(plan, fit, …)`) refuse a whitened plan: the fit the
caller holds is not identity-bound to the plan's relation. That is a refusal,
not a wrong answer, and it is the correct conservative behaviour — the sampling
covariance of a whitened estimand needs the residual second moment transformed
too, `Σ̃ = Q^{1/2}ΣQ^{1/2}`, and which estimand that names is a decision this
contract has not made. The refusal message is generic rather than explaining
the whitened case. **Owed:** either a targeted capability refusal naming the
composition, or the whitened calibration law itself (a D8 question, since §7.6a
and gap G8 already own "conservation gives no uncertainty").

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

**What does survive transport** (E2 must state this, because a blanket warning
invites the wrong conclusion that transported components are unusable): `P` is
linear and the decomposition `total = coherent + configuration` is an identity
per native node, so it holds after transport as well —
\(P^\top c_{\mathrm{coh}}+P^\top c_{\mathrm{cfg}}=P^\top c_{\mathrm{tot}}\)
exactly. Transported components are a valid *additive decomposition of the
transported budget*. What fails is only the reinterpretation of the coherent
part as a group-node common mode, and — by §4 — any *fraction* taken against a
summed coherent denominator.

**Required labelling.** The transported coherent and configuration fields must
not reuse the bare names `coherent` / `configuration`, which would let them be
read as group-node geometry. E2 fixes explicit transported names (e.g.
`transported_coherent`) and a print line naming the native frame family the
ledger came from. This contract does not fix the spelling; it fixes that the
names must differ and that the native provenance must be printed.

### 7.6a Conservation is a point-estimate law — **normative**

The conservation theorem (§2) is a statement about estimates. It carries **no**
implication for their uncertainty: \(\sum_x\widehat\theta_x=\widehat\theta_\Omega\)
does not give \(\sum_x\operatorname{SE}_x=\operatorname{SE}_\Omega\), because
node estimates of an overlapping frame are strongly positively correlated and
variances do not add. The variance of the conserved budget requires the full
cross-node sampling covariance, not a per-node margin.

A conserved ledger therefore buys no free error bars, and a population layer
must not build one by summing per-node standard errors. This is exactly what
D8 (`contrast_energy` sampling route) has to supply, and until it does, a
conservative attribution map is a point ledger reported without inference.

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
| frame **families** (multi-scale, α-weighted) | **delivered (D2):** `frame_family(..., alpha, normalization = "conservative")` (`R/frame.R`). Refuses `abs(sum(α) − 1) > tol` rather than renormalizing (G1) and validates column normalization **per member** rather than on the stack (G2). The undocumented `rbind(...)` + `additive_frame(..., "conservative")` route still works and is still provenance-blind (`test-frame-family.R`, "the bare rbind route is still provenance blind") | met |
| per-row metadata | **delivered (D2)** for the family route: `$index` carries `measurement` (`"<family>::<node>"`), `family`, `node`, `scale`, `center`, `alpha`, one row per measurement, and reaches the result's `$index` (`test-frame-family.R`, "family metadata survives a geometry evaluation"). A single compiled frame's `$index` is unchanged: still `measurement` only, radius still frame-wide on `$specification` | met for families; a single frame's own `$index` is unchanged |
| `additive_frame()` | still drops `$index` and `$specification` — deliberately, since a declared frame has no generator. `frame_family()` carries both instead, so the gap is closed on the family route only; the **metric fold** route of §5.3 still drops them | met for families; the metric fold still owes it (§5.3) |
| `contribution(by = region/network)` | **delivered (D4):** `contribution(x, by, using = NULL)` (`R/views.R`) aggregates an `effect_contrast_view` or a query-only `effect_view` **by row** over a grouping — a column of a per-measurement metadata table such as a family's `$index`, or a label vector the analyst supplies. Refusals: a locally normalized frame (§1.1), a view that does not record its frame normalization, and the non-additive readouts (`effect_spectrum_view`, `effect_rdm_view`, `effect_rsa_view`, `effect_crossnobis_view`, packed geometry, an unevaluated plan). `total` is budget-exact; `coherent` / `configuration` are labelled `frame_relative` (§4); `signed` is **masked**, because it is the local weighted *mean* and a density does not add over a territory; the group coherence fraction is recomputed from the aggregated components under the §4 mask, never averaged. Provenance in `$metadata$aggregation` (`aggregated_by`, `frame_relative`, `budget_exact`, `masked`, `overlap_split = FALSE`); the receipt carries a derived `scientific_plan_id` (`test-contribution.R`) | met |
| coherence spectrum | derivable from `contrast_energy()` per frame; no named object | a named per-scale / per-location coherent-share object |
| metric composition | **delivered (D6):** `plan_geometry(..., metric = , composition = c("native", "whitened"))`. Default `"native"`, unchanged and identity-stable. `"whitened"` is a `whitened_metric_before_frame` schedule carrying `composition` and `root = "symmetric_psd_root"` in `$scientific_plan_id`; the transform is relation-level (`B̃ = BQ^{1/2}` at plan time), per §5.2.2, so the compiler, executor, and kernel are untouched. Fixed positive-definite domain-wide metrics only; learned recipes refused (`whitened_metric_composition`) | met |
| metric conservation certificate | `.metric_frame_conservation()` / `.require_metric_conservation()`, internal; **D6** added `composition`, `root`, and `max_deviation` to the certificate, a `"whitened"` branch certifying `Σ_x Q^{1/2}D(w_x)Q^{1/2} = Q` by its exact residual, and a separate `.require_metric_conservation()` message for a certificate that conserves without being feature additive | met for the composition law; still not surfaced wherever a conservative claim is *printed* |
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
| `Σ_x K_x^wh = Q` (whitened composition, symmetric root) | `1e-12` relative (measured `1.5e-15`) |
| `Σ_x R D(w_x) Rᵀ = Q` for any root `RRᵀ = Q` | `1e-12` relative (measured `1.5e-15` symmetric, `2.9e-16` Cholesky) |
| root-dependence of *node* values under whitening | asserted as `> 1%` relative, **never** as an equality (measured `15.7%`, and `28.6%` on a second draw) |
| `⟨svec A, svec B⟩ = ⟨A,B⟩_F` | `1e-12` relative (measured `1.8e-15`) |
| subject-Gram spectrum vs `P×P` covariance spectrum | `1e-12` absolute (measured `2.7e-15`) |
| dense-metric non-conservation | asserted as `> 1%` relative, **never** as an equality |

---

## 10. Test and oracle index

| claim | statement | evidence |
|---|---|---|
| 1 | local = detection (intensive, row-normalized); conservative = attribution (extensive, column-normalized); non-conservation of local frames | `tests/testthat/test-integrity-guards.R` — "conservative frames conserve total evidence; local frames do not" (`frame_conservation()` TRUE vs FALSE); `tests/testthat/test-frame.R` — "conservative frames partition global feature mass" (column margin) **and "local normalization produces unit measurement mass" (row margin — the intensive half of the claim, uncited in the first draft)** |
| 2 | `Σ_x G_x = G_Ω` for a column-normalized frame | `tests/testthat/test-workflow.R` — "public frame laws preserve point decomposition and global total" (1e-13); `tests/testthat/test-integrity-guards.R` — "conservative frames conserve total evidence; local frames do not" (query level, 1e-10); `design/oracles/conservative-multiscale-ledger.R` §O1.a |
| 3a | total is the frame contraction of a fixed per-voxel ledger, **for a feature-additive (identity or diagonal) metric only** | `tests/testthat/test-conservative-geometry-contract.R` — "a conservative total is the frame contraction of a voxel ledger" (identity metric); `design/oracles/conservative-multiscale-ledger.R` §O1.b. The dense-metric counterexample that bounds the claim is claim 5b's `B(S∘Q)Bᵀ` law — **no test asserts the ledger picture *fails* for a dense `Q`**; D6 should add one |
| 3b | per-scale totals under an α-weighted family equal `α_s G_Ω` | `tests/testthat/test-conservative-geometry-contract.R` — "an alpha-weighted frame family conserves total scale by scale"; `design/oracles/conservative-multiscale-ledger.R` §O1.c |
| 3c | the coherence spectrum, not the energy spectrum, is the informative object | `design/oracles/conservative-multiscale-ledger.R` §O1.d (share column varies; `E_s` column does not) |
| 4a | coherent does not conserve | `tests/testthat/test-integrity-guards.R` — "conservative frames conserve total evidence; local frames do not" (`expect_gt(abs(Σ coherent − global coherent), 1e-8)`); `design/oracles/conservative-multiscale-ledger.R` §O1.d |
| 4b | coherence fraction is masked, never clamped | **Primary:** `tests/testthat/test-conservative-geometry-contract.R` — "contribution shares are undefined on the signed estimation layer", which asserts `any(!coherence_fraction_valid)` *before* `all(is.na(...))` and so exercises the masked branch. **Secondary (formula only):** `tests/testthat/test-views.R` — "contrast returns one exact decomposition and signed marginal"; `R/result.R:859-863`. *Citation corrected 2026-08-20: in `view_geometry_fixture()` both nodes are valid (`total` 5, 3; `coherent` 0.25, 0.75; `configuration` 4.75, 2.25), so that file's `expect_true(all(is.na(fraction[!valid])))` reduces to `all(is.na(numeric(0)))` and passes vacuously. It verifies the fraction on valid nodes, not the mask.* |
| 4c | singleton frames have zero configuration (why the point scale reads share 1) | `tests/testthat/test-effect-form-laws.R` — "a singleton positive-mass frame has exactly zero configuration" |
| 5a | identity and diagonal metrics conserve; dense metrics do not | `tests/testthat/test-conservative-geometry-contract.R` — "conservation survives a diagonal metric and fails for a dense one"; `design/oracles/conservative-metric-composition.R` §O2.a–c |
| 5b | the failure obeys `Σ_x G_x = B(S∘Q)Bᵀ` with `S_uu = 1` | `tests/testthat/test-conservative-geometry-contract.R` — "conservation survives a diagonal metric and fails for a dense one" (its `overlap` assertions); `design/oracles/conservative-metric-composition.R` §O2.c |
| 5c | the assessment's `−6.6%` is fixture-specific; sign and size both vary | `design/oracles/conservative-metric-composition.R` §O2.c′ |
| 5d | the whitened composition `Q^{1/2} D(w) Q^{1/2}` conserves and is a different estimand | `tests/testthat/test-conservative-geometry-contract.R` — "the whitened composition conserves where the native one does not" (uses the symmetric eigen root, lines 181-184); **`tests/testthat/test-composition.R` (D6)** — "both compositions equal their hand-computed dense algebra" (both compositions against hand-computed dense algebra at `1e-12`, and the estimand gap asserted `> 1%`) and "a dense metric conserves under whitening and not natively" (whitened `< 1e-12`, native `> 1%`), now against the executed package path rather than first principles; `design/oracles/conservative-metric-composition.R` §O2.d |
| 5g | conservation holds for **any** root `RRᵀ = Q`, but node values are root-dependent, so `"whitened"` alone does not name an estimand (§5.2.1) | `design/oracles/conservative-metric-composition.R` §O2.d′ (added by the 2026-08-20 review: both roots conserve to `1.5e-15` / `2.9e-16`, node values differ by `15.7%`). **Test added (D6, gap G5 closed):** `tests/testthat/test-composition.R` — "whitening means the symmetric root, not any root that conserves" asserts both roots conserve at `1e-12`, that their node values differ by `> 1%`, and that the executed package values equal the *symmetric*-root values at `1e-12` while differing from the Cholesky ones by `> 1%`. A Cholesky implementation passes every conservation assertion in that test and fails this one |
| 5h | whitening is not support-local, so it cannot go through `.compose_frame_metric()` (§5.2.2) | source-level: `.compose_frame_metric()` requires `identical(metric$support, node_value$support)` (`R/metric.R`, HEAD). **Discharged by construction (D6, gap G6 closed):** the whitened route never reaches `.compose_frame_metric()` — it is a relation-level transform feeding the identity path — and `tests/testthat/test-composition.R` — "a whitened plan declares the transform it performed" pins that, asserting the plan lowers to `additive_contraction` and that `.metric_additive_frame()` refuses the whitened schedule outright |
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

**Reviewer note (2026-08-20): the oracles are not executed by anything.**
`grep -rl design/oracles` over the repository matches only the three scripts
themselves and this document — no test, no `Makefile`, and none of the four
workflows in `.github/workflows/` runs them. All three do run clean today
(exit 0, verified 2026-08-20) and every number quoted in this contract
reproduces exactly. But claims **3c, 5c, 5g, 7b, 7c, 7d and 7e** have *no other
evidence*, so their only protection against silent rot is a human remembering
to run a script. That is below the standard the rest of this contract sets, and
it will not survive the API churn WS-D is about to introduce.

Recommended (D2, since it is the first ticket to touch this material): add a
skipped-by-default `test-conservative-oracles.R` that shells out to the three
scripts under `testthat::skip_on_cran()` and asserts exit status plus a few
anchor numbers, or promote the six oracle-only claims into
`test-conservative-geometry-contract.R`. The second is preferable — the oracles
should stay readable derivations, not become the test suite.

---

## 11. Fresh-context review (2026-08-20)

Independent review of `conservative-geometry-v1` at commit `7ec6e37`, against
`HEAD` = `a27b246`. Scope: re-derivation of every numbered claim in §1–§7,
citation fidelity, evidence audit of the §10 index, execution of the three
oracles, and gap analysis for WS-D / WS-E.

**Verdict: accept with amendments.** The mathematics is sound. Every claim in
§1–§7 re-derives correctly, all three oracles run clean and reproduce every
quoted figure exactly, and every file:line citation resolves at the contract's
own commit. The amendments below are applied in place; none of them overturns a
claim, but two of them (§3's missing metric precondition, §5.2.1's root
ambiguity) change what an implementer must build, so this document should not
be treated as frozen until D6 confirms §5.2.1.

### 11.1 Claims verified

Re-derived independently and confirmed: **§2** (exchange of summation; both
preconditions are real and are the actual failure modes), **§3.1** (per-scale
`α_s G_Ω` follows from block-wise column normalization), **§3.2** (the coherent
share is the only data-dependent column), **§4** (node-local rank-one
projections have no linear operator summing them to the global one), **§5.1**
(`Σ_x K_x^native = S ∘ Q` with `S_uv = Σ_x √(w_xu w_xv)`; `S_uu = 1` by column
normalization and `S_uv ∈ [0,1]` by Cauchy–Schwarz, consistent with the measured
`[0, 0.8165]`), **§5.2** (whitened conservation, subject to §5.2.1), **§5.3**
(the `"none"` declaration is correct, not lossy), **§6** (conservation is a
signed-sum statement and does not make summands nonnegative; the arithmetic in
the measured block is internally consistent), **§7.2–7.5** (Frobenius isometry,
rank ≤ N−1 subject-Gram trick, row-stochastic budget preservation, budget vs
density).

**Oracles.** All three exit 0. Every number quoted in §2, §3.1, §3.2, §4, §5,
§5.2, §5.3, §6, §7.1 and §7.2–§7.5 matches the printed output to the digits
shown — including the twelve-draw sweep in §5, reproduced element for element.
Oracle fidelity is the strongest part of this contract.

**Evidence audit.** All fifteen cited test names exist. Fourteen of fifteen
citations test the claim they are cited for, non-vacuously; the exception is
claim 4b (§11.3). The three most load-bearing tests —
"an alpha-weighted frame family conserves total scale by scale",
"conservation survives a diagonal metric and fails for a dense one" (whose
`overlap` block asserts the `B(S∘Q)Bᵀ` law exactly as claim 5b describes), and
"the whitened composition conserves where the native one does not" — are precise
and independent of the production contraction path.

### 11.2 Corrections applied

1. **§3 was stated without its metric precondition** (over-statement, now
   fixed). "The total field at any scale … carries no cross-voxel information
   whatsoever" is true only on the feature-additive branch. For a dense `Q` the
   native composition gives `G_x = Σ_{u,v} √(w_xu w_xv) Q_uv b_u b_vᵀ`, which
   pairs distinct voxels inside a node; there is no per-voxel ledger and §3.1
   does not apply. Claim 3 is now scoped to
   `metric_capabilities()$feature_additive`. This mattered: §3 is the basis of
   the "multiscale energy panels are not findings" prohibition, and §5 is
   precisely about the case where its premise fails.
2. **§5.2 named an estimand that is not determined** (new §5.2.1). Conservation
   holds for *any* root `RRᵀ = Q`, not only the symmetric one, but the per-node
   values are root-dependent. `conservative-metric-composition.R` gained §O2.d′
   this review: symmetric PSD root and lower Cholesky factor both conserve
   (`1.5e-15`, `2.9e-16` relative) while their node forms differ by **15.7%** of
   the largest node value — comparable to the native-vs-whitened gap of 13.5%
   that §5.2 leads with, and `28.6%` on an independent draw. `composition =
   "whitened"` is now normatively defined as the symmetric PSD root, with the
   root identity required in plan identity. A conservation certificate cannot
   distinguish two roots, so it is not sufficient evidence of a reproduced
   analysis.
3. **§4's fixture comparison was numerically wrong** (now corrected). The
   coherent deviations `7.41 / 2.70 / 1.27` were compared to "`max |G_Ω|` … order
   1"; measured, `max |G_Ω| = 8.62`. The correct comparator is the global
   *coherent* component, `max |G^coh_Ω| = 1.21`, which supports the point better.
4. **§4 normative 3 misstated the mask guard** (now corrected): the code
   conjunct is `is.finite(total) & total > 0 & coherent >= 0 & configuration >= 0`;
   `is.finite` was omitted.
5. **§3.2 understated its own result** (strengthened, flagged as a reviewer
   note). φ_s is not merely "not proportional to α" — it is *exactly invariant*
   to α_s, because total and coherent are both homogeneous of degree one under a
   row rescaling (`a = w_x/Σw_x` is α-invariant while `K_x⁻¹` scales as `1/α`).
   This is what makes the coherence spectrum well-posed and reportable without
   disclosing α, and D5 should be built on the stronger statement.
6. **§3.1 gained two enforcement preconditions**: each family member must be
   column-normalized *on its own* (the stack can conserve while no block does),
   and the α-sum rule must be either enforced or refused by the constructor,
   never assumed.
7. **§7.6 gained what survives transport, and §7.6a is new.** The prohibition on
   reading transported coherent as group-node geometry invited the wrong
   inference that transported components are unusable; `total = coherent +
   configuration` survives `P` exactly by linearity, and E2 should say so.
   §7.6a records that conservation is a point-estimate law with no implication
   for standard errors.

### 11.3 Citation fixes

- **Claim 4b's primary citation was vacuous.** `test-views.R` "contrast returns
  one exact decomposition and signed marginal" asserts
  `all(is.na(fraction[!valid]))`, but in `view_geometry_fixture()` both nodes are
  valid (`total` 5, 3; `coherent` 0.25, 0.75; `configuration` 4.75, 2.25), so the
  subset is empty and the assertion passes on `all(is.na(numeric(0)))`. It
  verifies the fraction formula, not the mask. Repointed: the non-vacuous
  evidence is `test-conservative-geometry-contract.R` "contribution shares are
  undefined on the signed estimation layer", which asserts
  `any(!coherence_fraction_valid)` first. `test-views.R` is retained as the
  formula check.
- **Claim 1's row-margin half was uncited.** Added `test-frame.R` "local
  normalization produces unit measurement mass".
- **All other line numbers verified at `7ec6e37`**, including the exact ones:
  `R/result.R:859-863` is precisely the mask block; `R/scope.R:430-436` is
  precisely the conservative column check; `R/kernel.R:354-377` is the
  `weight_slice %*% atom_slice` contraction; `src/packed-form.cpp:144-222` is
  `coherent_effect_form_atoms_cpp`, which does compute `(Bw)(Bw)ᵀ/Σw` from first
  moments with `√2` off-diagonal packing; `R/frame.R:355-423`, `:396-423`,
  `:425-440` and all four `R/metric.R` spans resolve. `design/effect-form-contract.md`
  §8 (A3) exists and says what §5.1 attributes to it. No stale citation found —
  only the working-tree drift noted in the preamble.
- **Claim 5e is slightly under-cited**: the cited `test-metric.R` test asserts
  the refusal and `feature_additive = FALSE` but not the
  `global_metric_kind = "support_pair_operator"` string that §5 quotes; only the
  oracle covers that. Left as-is, noted here.
- **The oracles are wired to nothing** — see the note at the end of §10. Six
  claims are oracle-only and unprotected against rot.

### 11.4 Gaps for WS-D / WS-E

Each with a one-line proposed resolution and owning ticket.

| # | gap | proposed resolution | ticket |
|---|---|---|---|
| G1 | α-weight normalization across a family is assumed, never specified: renormalize, or refuse? May α be per-row rather than per-scale? | Constructor refuses `abs(sum(alpha) - 1) > tol` and records the applied α per row; per-row α permitted only as a within-scale reweighting that preserves the block column sum | D2 |
| G2 | Each family member must individually cover the domain for §3.1's block identity; the stack conserving does not imply it | Validate column normalization per block at construction, not on the stack | D2 |
| G3 | "Coherent share versus location" is underdetermined when a location sits in several scales — it is a function of (location, scale), not a number | Define the spectrum on `(center, scale)` and require any location-wise collapse to be a declared, named reduction (α-weighted mean, argmax-scale, …) | D5 |
| G4 | ~~`contribution(by = region)` granularity is unspecified: grouping rows *by center* preserves the budget exactly, splitting an overlapping node's mass across regions requires a second partition and can double-count~~ **CLOSED (D4).** Grouping is by row: every row belongs to exactly one group, so the group sums partition the budget exactly, and `$metadata$aggregation$overlap_split` records `FALSE` rather than leaving it implied. Overlap-splitting is not offered; if it is ever wanted it arrives as its own declared, separately certified reduction. Two things the gap did not anticipate: `signed` had to be **masked** rather than summed (it is the local weighted mean, already divided by the node's own frame mass, so it is a density), and the guard needs the frame's *declared* normalization, since a diagonal-metric fold leaves `normalization = "none"` on a conservative frame (§5.3) — `.execution_metadata()` now records `declared_normalization` alongside it | Default to grouping by row center (budget-exact); offer overlap-splitting only as a declared, separately certified reduction | D4 |
| G5 | The whitened root is not part of the estimand as written (§5.2.1) | Pin the symmetric PSD root and record root identity in plan identity; oracle §O2.d′ added this review, the matching test is still owed. **Closed (D6):** the schedule carries `root = "symmetric_psd_root"` in its semantic digest, and `test-composition.R` — "whitening means the symmetric root, not any root that conserves" is the owed test | D6 ✔ |
| G6 | Whitening is not support-local and cannot go through `.compose_frame_metric()` (§5.2.2) | Implement as a relation-level transform `B̃ = BQ^{1/2}` feeding the identity-metric path; budget accordingly. **Closed (D6):** implemented exactly so, at the `plan_geometry()` seam; the resident bytes are recorded in `$execution_hints` and enforced against a declared workspace budget before the first read | D6 ✔ |
| G7 | "Nonnegativity projection" is underdetermined — a per-node total clamp and an eigenvalue truncation of the form are different operators moving different mass | Latent layer takes a named projection from a closed set, and the receipt records moved mass per projection kind | D7 |
| G8 | Conservation gives no uncertainty; nothing currently supplies cross-node sampling covariance for `contrast_energy` (§7.6a) | Named `contrast_energy` sampling route returning the cross-node covariance, not per-node margins | D8 |
| G9 | Between-subject budget heterogeneity is unaddressed: `⟨H, G_Ω⟩` differs per subject, so a group ledger sums incommensurable budgets | Population contract declares a per-subject budget normalization (none / unit-budget / precision-weighted) as a plan-identity field | E1 |
| G10 | §7.6 fixes that transported components must be renamed but not to what | E2 fixes the names and the print line naming the native family | E2 |
| G11 | Node labels still do not survive the diagonal-metric fold (§5.3) or the family route (§7.1) | Family constructor carries `$index`/`$specification` through both | D2 |
| G12 | Oracle-only claims have no regression protection; nothing in the repo runs `design/oracles/` | Promote 3c, 5c, 5g, 7b–7e into `test-conservative-geometry-contract.R` | D2 |

### 11.5 Over-statements found

Two, both now corrected in place: §3's unscoped "no cross-voxel information
whatsoever" (§11.2 item 1) and §5.2's under-determined "the whitened
composition" (§11.2 item 2). §4's fixture comparison was a third, of a smaller
kind (§11.2 item 3). Everything else in §1–§7 is stated at, not above, the
strength its evidence supports — and §3.2 is stated *below* it.
