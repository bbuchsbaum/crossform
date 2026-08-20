# Population form contract — mass-preserving transport of conservative geometry

Status: normative architecture contract

Contract version: `population-form-v1`

Date: 2026-08-20

This document freezes the semantics of the WS-E **population layer**: how each
subject's conservative geometry is carried to a common group node set, what the
carrying operation preserves, what it destroys, and what the group form
estimates. It is normative for WS-E and consumes
`design/conservative-geometry-contract.md` (`conservative-geometry-v1`) as its
input interface — in particular §7 (transport readiness) and gaps **G9**
(per-subject budget normalization) and **G10** (labelling of transported
components), which this contract closes.

It sits below `design/effect-form-contract.md` (`effect-form-v1`), which fixes
the bilinear-form algebra, the packed codec, the canonical operation order and
the noise-unbiasedness theorem. Nothing here relaxes an `effect-form-v1` law,
and nothing here relaxes a `conservative-geometry-v1` law: the population layer
is a *linear* layer stacked on top of a conservative result, and every
non-linear step it adds (a ratio, a PSD projection, a fraction) is named and
confined to the latent descriptive layer of `conservative-geometry-v1` §6.

Three oracle scripts under `design/oracles/` support the numbered claims. They
are **pure matrix algebra** — the `crossform` package is not loaded and no
production contraction path is exercised — so they are independent evidence
rather than a re-execution of the code they govern. Every number quoted below
as *measured* is their output at the commit that introduced them. Three places
are **not** measurements and say so where they appear: §4.1's equal-precision
degeneracy and §6.1's exact identity are definitional, and §9.1 is a stipulation
of the required interface transcribed into executable predicates. §12 maps each
claim to its oracle section.

```sh
Rscript design/oracles/population-transport-contract.R      # claims 1-4
Rscript design/oracles/population-geometry-split.R          # claims 5-6
Rscript design/oracles/population-transport-diagnostics.R   # claims 7-9
```

**Nothing in this contract exists in `R/` today.** Verified by grep over `R/`
and `NAMESPACE` at `elite-pass`: no `location_transport`, no `plan_population`,
no `population_form`, no `transport_object`. §10 is the full gap table. This
document is therefore a specification, not a description; every "the
implementation must" is an unmet obligation on E2 and later.

---

## Notation

| symbol | meaning |
|---|---|
| `i = 1…N` | subjects |
| `x = 1…n_i` | subject `i`'s **native** conservative frame nodes |
| `j = 1…m` | **group** nodes; `⊥ = m+1` is the sink |
| `q` | query coordinates (conditions); forms are `q × q` symmetric |
| `p = q(q+1)/2` | packed width under `symmetric_packed` |
| `G_{i,x}` | subject `i`'s form at native node `x` (signed, crossvalidated) |
| `Y_i ∈ R^{n_i × p}` | svec'd native forms, one row per native node |
| `c_i ∈ R^{n_i}` | a *queried* native ledger, `c_{i,x} = ⟨H, G_{i,x}⟩` |
| `P_i ∈ R^{n_i × (m+1)}` | subject `i`'s transport |
| `T_i = Σ_x c_{i,x}` | subject `i`'s native budget (`= ⟨H, G_Ω^{(i)}⟩` by §2 of `conservative-geometry-v1`) |
| `X ∈ R^{N × k}` | the group design matrix |

---

## 1. The transport object

**Claim 1.** A transport is a typed, sparse, provenance-bearing value carrying
exactly six things. It is an *input* to `crossform`, never something
`crossform` estimates (§9).

### 1.1 Required structure

`P_i` is a sparse nonnegative matrix with `n_i` rows and `m + 1` columns, whose
**last column is the sink**, satisfying

\[
P_i \ge 0,\qquad \sum_{j=1}^{m+1} (P_i)_{xj} = 1 \ \ \text{for every native row } x .
\]

Row-stochasticity is asserted *including* the sink. The sink column is
**required and always materialized**, even when it is identically zero: its
presence is what makes partial coverage a visible number rather than a silent
loss (`conservative-geometry-v1` §7.4). A transport with `m` columns and no
sink is not a degenerate transport; it is ill-formed, and the constructor
refuses it.

Measured (`population-transport-contract.R` §P1.a) on a 12-row, 4-group-node
fixture with one row splitting `0.6 / 0.4` across two group nodes, one row
`0.7 / 0.3` between a group node and the sink, and one row entirely in the
sink: `max |rowSums − 1| = 0`, 14 nonzeros of 60 (76.7 % sparse). Four
well-formedness failures are each detected by the predicate set — sink column
dropped (`has_sink_column` and `row_stochastic` both `FALSE`, worst row deficit
`1`), rows summing to `0.8`, a negative entry, and a functional transport
arriving without cross-fit provenance.

### 1.2 Required fields

| field | content | why it is required |
|---|---|---|
| `matrix` | the sparse `n_i × (m+1)` operator | the operator itself |
| `native_index` | one row per native node: `family`, `scale`, `center`, `label`, `alpha` | §9; a transport maps *named* rows, and `conservative-geometry-v1` §7.1 records that today's frames cannot name them |
| `group_index` | `m` group-node labels and coordinates, plus the sink | displacement and coverage diagnostics (§7) need group coordinates |
| `semantics` | `"budget"` or `"density"` | §1.3; changes the estimand |
| `row_mass` | declared `μ_i ∈ R^{n_i}_{>0}` | §1.3; the denominator of density semantics |
| `provenance` | `kind ∈ {"anatomical","functional"}`, the construction method, and — for `"functional"` — `cross_fit` naming the runs/tasks used | §7; without it a circular transport is indistinguishable from an honest one |

`semantics`, `row_mass` and the transport's content hash all enter the group
**scientific plan identity**, not merely the physical signature: two group
analyses that differ in `P` estimate different things, and a receipt that hides
that is a false receipt.

### 1.3 Budget and density are two linear maps built from one `P`

**Budget semantics.** The group value is the transported sum,

\[
b_{ij} = (P_i^\top c_i)_j .
\]

A group node receives mass in proportion to how many native nodes map into it,
so a subject with a finer native frame weighs more per territory
(`conservative-geometry-v1` §7.5, measured there).

**Density semantics — defined exactly.** Density is a **declared row-mass
normalization**: the transport carries a positive vector `μ_i` giving each
native row's own territory measure (default `μ ≡ 1`, the native node count),
and the group value is transported budget per unit transported mass,

\[
d_{ij} = \frac{(P_i^\top c_i)_j}{(P_i^\top \mu_i)_j},
\qquad
d_{ij} := \mathrm{NA}\ \text{when}\ (P_i^\top \mu_i)_j = 0 .
\]

With `μ ≡ 1` this is exactly the arithmetic of
`conservative-transport-readiness.R` §O3.d. A group node reached by no native
mass returns `NA`, never `0`: zero is a measurement and absence is not.

Three normative consequences.

1. **Density is still linear in the data.** `d = D(1/(P^⊤μ)) P^⊤ c` is a fixed
   linear operator, because `P` and `μ` are declared, not estimated. Measured:
   the explicit linear map reproduces the ratio to `8.3e-17`
   (§P1.b). Everything §3 proves about commutation therefore holds for density
   as well. What density gives up is **conservation**, not commutation —
   measured on the same fixture, the budget columns sum to `−3.750403187968`
   (the native total minus the sink, §2) while the density columns sum to
   `−1.361403033504`, which satisfies no conservation law at all.
2. **The sink is always reported in budget units**, under either semantics. It
   is a mass-accounting column, not an estimate of anything at a location; a
   "density" of unmapped territory has no referent.
3. **Neither is a default.** `semantics` has no default value; the constructor
   requires it.

### 1.4 Provenance, and the cross-fit obligation

`kind = "anatomical"` means `P` was built from a registration or a parcellation
that never saw the response data; the provenance records the warp or atlas
identity. `kind = "functional"` means `P` was built from data, and then
`cross_fit` is **required**: it names the runs, sessions or tasks used to build
`P`, and the plan refuses to evaluate on any partition listed there. §7 measures
what happens when this is not enforced.

### 1.5 `P` is part of the estimand

The group form is not "the population geometry"; it is the population geometry
*as resolved by `P`*. Two transports give two different estimands, both valid,
and the contract's position is the same as `conservative-geometry-v1` §5.2 takes
on `composition`: the choice is recorded in plan identity and never applied
silently.

---

## 2. Budget preservation

**Claim 2.** Under budget semantics, each subject's transported total over the
group nodes equals its native total minus the sink mass, exactly.

*Proof.* `P_i` is row-stochastic, so `P_i 1 = 1` and

\[
\sum_{j=1}^{m+1}(P_i^\top c_i)_j = \mathbf 1^\top P_i^\top c_i
= (P_i\mathbf 1)^\top c_i = \mathbf 1^\top c_i = T_i .
\]

Splitting the left sum at the sink gives
`Σ_{j≤m} b_{ij} = T_i − b_{i⊥}`. ∎

The argument uses only `P_i 1 = 1`; it never uses `c_i ≥ 0`. **This is a
signed-sum law, not a mass law** — the word "mass" throughout this contract is
bookkeeping vocabulary for a signed budget, exactly as in
`conservative-geometry-v1` §6.

Measured (§P2): on the fixture of §1.1 with 8 of 12 native entries negative,
`|group + sink − native| = 8.9e-16` and `|group − (native − sink)| = 4.4e-16`;
over 500 random transports with fractional rows and random sink mass, the worst
deviation is `7.1e-15`. Deleting the sink column from the *same* `P` silently
loses `48.1 %` of the native total, because 1.3 rows of native row mass were
pointed at the sink.

**Normative.** A population plan asserts Claim 2 as a certificate on every
subject at fit time. The tolerance is **scaled by the ledger's `L¹` norm**,
`1e-12 · Σ_x |c_{i,x}|`, not by the total: the total is a signed sum and may sit
arbitrarily close to zero while the summands are large, so a relative-to-total
tolerance is unbounded exactly where the estimate is weakest. A subject failing
the certificate has an ill-formed transport, and the plan refuses rather than
renormalizing.

---

## 3. Commutation, and where it stops

**Claim 3.** Let the query `⟨H,·⟩`, the transport `P_i^⊤`, and the group fit act
on the **coordinate**, **node** and **subject** axes respectively. If the fit's
weight operator is constant across the node axis and the coordinate axis, all
three are linear maps on distinct tensor factors and every evaluation order
gives the same answer.

*Proof sketch.* Write the data as `A_i ∈ R^{n_i × p}`. Transport is
`A_i ↦ P_i^⊤ A_i` (left multiplication, node axis); query is `M ↦ M h` with
`h = svec(H)` (right multiplication, coordinate axis); the fit is, for each
`(node, coordinate)` cell, `Y_{·jk} ↦ (X^⊤ΩX)^{-1}X^⊤Ω Y_{·jk}` (subject axis).
Left multiplication, right multiplication, and a cellwise map along a third
index pairwise commute; the order-independence follows. ∎

**Scope of the tensor argument.** With subject-specific `P_i` on heterogeneous
`n_i` there is no common node axis, and *fit-then-transport is not even
definable* — the fit needs the subjects stacked at a shared node. The argument
as stated therefore licenses query↔fit and query↔transport for any `P_i`, and
licenses **transport↔fit only on a common native node set** (the slice-1
geometry, where every subject carries the same parcellation). §P3.a exhibits
the first two by putting transport first in every order; the `2.2e-16` row of
§P3.b is what exhibits transport↔fit. This is a real restriction, not a
technicality: with heterogeneous native frames the transport *must* precede the
fit, and the compiler has no ordering freedom there to exploit.

Measured (§P3.a), `N = 8` subjects with heterogeneous native node counts,
`q = 4`, `m = 5`, an intercept-plus-covariate design and signed indefinite
forms:

| comparison | max `|Δβ|` |
|---|---|
| query→transport→fit vs transport→fit→query | `8.9e-16` |
| query→transport→fit vs transport→query→fit | `6.7e-16` |
| GLS with a common `Ω`, order A vs order B | `1.3e-15` |

**The commuting class is larger than OLS.** GLS with a *subject-constant* `Ω`
(the same `N × N` operator at every node and coordinate) commutes too. What
breaks commutation is weight *variation along an axis a later operator mixes*.

### 3.1 Counterexample: weights varying along the node axis

Transport mixes nodes. If the fit's weights vary per node, fitting then
mapping ≠ mapping then fitting, and no choice of aggregated weight repairs it.
Measured (§P3.b) on a slice-1 geometry (six shared native regions aggregated
into three networks, so both orders are defined): with weights constant across
nodes the two orders agree to `2.2e-16`; with per-node weights they differ by
`3.55e-01`, **25.2 %** of the largest coefficient.

The exact condition is that the two orders agree **iff the per-node WLS
operators coincide on each group node's support** — `W^{(x)} = W^{(j)}` for
every native `x` with `P_{xj} > 0`. §P3.b tests one aggregation rule (the
mass-weighted mean weight); the claim is that algebraic condition, not that one
rule happened to fail.

### 3.2 Counterexample: weights varying along the coordinate axis

The query mixes coordinates. If the fit's weights vary per packed coordinate —
which is exactly what per-cell precision weighting does — then querying the
fitted coefficients ≠ fitting the queried values. Measured (§P3.c): coefficients
`(+3.866810, −2.471789)` versus `(+3.183510, −1.487006)`, a gap of `9.85e-01` or
**25.5 %** of the largest coefficient; the same code with coordinate-constant
weights agrees to `8.9e-16`.

### 3.3 Normative

1. **OLS is the default group fit.** GLS with a subject-constant `Ω` is an
   admitted commuting variant.
2. **Per-node and per-coordinate weighting are explicit modes that give up
   query and aggregation commutation.** They are not forbidden — a
   precision-weighted group fit is a legitimate estimand — but selecting one
   changes what the compiler is allowed to do.
3. **The sink column is fitted but never reported as a group estimate.** It is
   carried through the group fit (§P3.a fits all `m+1` columns) because its
   trajectory across subjects and covariates is a diagnostic — systematic
   covariate-linked sink mass means the transport is failing differentially.
   It is excluded from `Z` wherever a geometry is formed (§7.1) and is never
   returned as a value at a location.
4. **When commutation fails, the evaluation order becomes part of the
   estimand** and is recorded in the scientific plan identity. The compiler may
   choose the order freely *only* in the commuting class; outside it, the order
   is a declared field and the receipt states it.
5. This refines, and does not contradict,
   `.planning/2026-08-17-feedback-assessment.md` Part 1 correction 8. That
   correction says "per-node or per-coordinate weighting breaks the query↔
   regression commutation". Precisely: **per-coordinate** weights break
   query↔fit; **per-node** weights leave query↔fit intact but break
   aggregation↔fit (and therefore transport↔fit). Both are real, they are
   different failures, and a contract that conflates them would let a
   per-group-node weighted fit ship with a false commutation certificate.

---

## 4. Per-subject budget normalization (closes gap G9)

**Claim 4.** Subjects have incommensurable native budgets `T_i`, so a group
ledger that sums or averages them without a declared normalization sums
incommensurable quantities. The normalization is a closed set of three, it is
part of plan identity, and the three give materially different answers — up to
and including a different argmax.

### 4.1 The closed set

Each mode is a **group functional written out in full**. There is deliberately
no single `Σ_i s_i c_i / Σ_i s_i` template: `none` and `precision_weighted` are
weighted means of the ledgers, while `unit_budget` is the *unweighted* mean of
*unit-normalized* ledgers, and no common divisor covers all three. (An earlier
draft of this section did state a single template; it is wrong —
`Σ_i (c_i/T_i) / Σ_i (1/T_i)` sums to the harmonic mean of the `T_i`, not to 1.)

| mode | group functional `g` | estimand | conserves | requires |
|---|---|---|---|---|
| `none` (default) | `N^{-1} Σ_i c_i` | the **mean subject ledger** in native evidence units; subjects with more evidence dominate | yes: `Σ_j g_j = mean_i T_i` | nothing |
| `unit_budget` | `N^{-1} Σ_i c_i / T_i` | the **mean subject attribution *share*** — dimensionless, every subject contributes budget 1 | yes: `Σ_j g_j = 1` exactly | `T_i` bounded away from zero; see §4.3 |
| `precision_weighted` | `Σ_i π_i c_i / Σ_i π_i` | the **precision-weighted mean ledger**; `Σ_j g_j = Σ_i π_i T_i / Σ_i π_i` | yes, against that total | a declared external precision; see §4.5 |

**All three conserve** — each is a per-subject scalar rescaling followed by an
affine combination whose weights sum to one, and the node-sum identity of
`conservative-geometry-v1` §2 survives both. What differs is *which total* is conserved. Measured (§P4):
`0`, `0` and `1.1e-16` respectively.

### 4.2 They are not small differences

Measured on six subjects whose native totals span `4.0 … 0.4` and whose spatial
profiles differ (three high-budget subjects favouring `g1`, three low-budget
subjects favouring `g4`):

| mode | g1 | g2 | g3 | g4 | total |
|---|---|---|---|---|---|
| `none` | `+0.810833` | `+0.521667` | `+0.345833` | `+0.321667` | `+2.000000` |
| `unit_budget` | `+0.255000` | `+0.203333` | `+0.200000` | `+0.341667` | `+1.000000` |
| `precision_weighted` | `+0.050050` | `+0.071175` | `+0.122019` | `+0.298586` | `+0.541829` |

`argmax` is group node **1** under `none` and **4** under both alternatives.
Maximum relative difference `none` vs `unit_budget` is **68.55 %**, `none` vs
`precision_weighted` **93.83 %**. `precision_weighted` reduces to `none` exactly
when the precisions are equal — a definitional identity, pinned in the oracle
(`0`) rather than measured.

**Normative:** the mode is not a display option. It is a plan-identity field
with no silent default beyond `none`, and the printed result names it.

### 4.3 `unit_budget` has two failure modes, and they are declared

1. **It equalizes subjects regardless of evidence.** Measured: weakening one
   subject's total from `0.40` to `0.02` leaves the `unit_budget` ledger
   *bit-identical* (`max |Δ| = 0`). A subject with essentially no evidence
   carries the same weight as the strongest subject. That is the definition of
   the estimand, not a bug — but it must be stated where the number is printed.
2. **Its divisor is a signed estimate.** By `conservative-geometry-v1` §6 the
   ledger is signed, so `T_i` can be near zero or negative. Measured: a subject
   with ledger `(+0.30, +0.25, −0.50, −0.30)` has `T_i = −0.25` and
   `unit_budget` shares `(−1.20, −1.00, +2.00, +1.20)` — the two nodes holding
   **positive** evidence receive **negative** shares, and `|share|` reaches
   `2.0`.

**Normative.** `unit_budget` (a) must be computed against a `T_i` estimated
from partitions independent of `c_i`, or else declared descriptively as a
same-data ratio in the receipt; and (b) must refuse, per subject, when `T_i` is
not bounded away from zero by a declared criterion, marking that subject
`NA` rather than emitting a divergent share. **The refusal threshold is an
open maintainer decision (§14.2)** — it needs a per-subject standard error,
which `conservative-geometry-v1` gap G8 records does not exist yet.

### 4.4 The normalization commutes with transport and with the query

`s_i` is a scalar on the subject axis, so applying it to the native ledger and
applying it after transport give identical results. Measured (§P4) on six
subjects with native node counts `7/8/9/10/11/12` and a genuine sparse
transport per subject: `max |scale-then-transport − transport-then-scale| =
1.1e-16`. (An earlier draft cited a check that rescaled the already-transported
table and so exercised no transport at all — the vacuous-citation failure
`conservative-geometry-v1` §11.3 caught for its own claim 4b.) It
does **not** commute with the fit: rescaling responses per subject changes the
fitted coefficients, which is precisely why it is an estimand choice and not an
implementation detail.

### 4.5 `precision_weighted` is gated

`π_i = 1/\operatorname{Var}(\widehat T_i)` requires the variance of a *conserved
budget*, which by `conservative-geometry-v1` §7.6a needs the full cross-node
sampling covariance and not a sum of per-node margins — the object D8 / gap G8
owes. Until D8 lands, `precision_weighted` accepts `π` only as an **externally
supplied, provenance-bearing vector**, and the plan refuses to synthesize one.

---

## 5. The subject-Gram trick

**Claim 5.** Geometry-space covariance across `N` subjects has rank at most
`N − 1`, and the `N × N` Gram of svec'd transported forms recovers its nonzero
spectrum and its modes exactly. The `√2` off-diagonal weighting of
`symmetric_packed` is load-bearing.

*Proof sketch.* Let `Y ∈ R^{N×D}` stack the svec'd transported forms
(`D = m·p`) and `C = Y − 1\bar y^\top`. Then `Σ = C^⊤C/(N−1)` and
`Γ = CC^⊤/(N−1)` share every nonzero eigenvalue, and `Γu = λu ⇒ Σ(C^⊤u) =
λ(C^⊤u)`. Centering is across subjects, so `1_N^⊤C = 0`: the all-ones vector
lies in the left null space and `rank C ≤ N−1`. The identification of Euclidean
geometry on packed rows with Frobenius geometry on forms is `effect-form-v1`
§3, and holds only with `√2` off-diagonals. ∎

Measured (`population-geometry-split.R` §P5), `N = 5`, `q = 6` (`p = 21`),
`m = 4`, `D = 84`:

| quantity | measured |
|---|---|
| `rank(Σ)` | `4` (bound `N−1 = 4`) |
| `max |λ_k(Σ) − λ_k(Γ)|`, `k ≤ N−1` | `1.07e-14` |
| `max |λ_k(Σ)|`, `k > N−1` | `7.39e-15` |
| `\|Σ v₁ − λ₁ v₁\|` for `v₁ = C^⊤u₁/‖·‖` | `2.67e-15` |
| `|⟨v₁, first eigenvector of Σ⟩|` | `1.000000000000000` |
| `Γ − Σ_j Γ_j` (node-additivity) | `3.55e-15` |
| eigen cost ratio `(D/N)³` | `4.74e+03` |

**The `Σ` and `Γ` of this section are plug-in objects.** Formed from
partition-averaged forms they estimate `Σ_B + Σ_ε/R`, not `Σ_B` — §6.2 measures
that inflation at **+62.7 %**. Claim 5 is a statement about *linear algebra*
(rank, spectrum, mode recovery, streaming), and it holds for whatever `N × N`
matrix is fed to it. **Which matrix a population plan is allowed to feed it is
fixed by §6.4, and it is the cross-fitted Gram, not the plug-in one.**

Two consequences the implementation must use.

1. **The `D × D` covariance is never formed.** Modes come back as `C^⊤u`, and
   unpack to symmetric `q × q` forms (measured `max |M − M^⊤| = 0`).
2. **The Gram is additive over group nodes** — `Γ = Σ_j Γ_j` — so it accumulates
   node by node while streaming and the `N × D` stack is never materialized.
   This is the population-layer analogue of the streamed spatial contraction in
   `R/kernel.R`.

**The `√2` is not cosmetic.** Measured over 200 random symmetric pairs: the
`symmetric_packed` codec reproduces `⟨A,B⟩_F` to `3.11e-15` worst relative,
while a naive `upper.tri` packing is off by `2.96e+00` — an `O(1)` error, not a
tolerance. On the same fixture the naive packing's top "covariance" eigenvalue
is `18.846001934` against the correct `24.734361948`. A population spectrum
computed from a non-isometric packing is not a spectrum of geometries.
`crossform`'s codec is already correct (`crossform:::.svec_symmetric`,
`conservative-geometry-v1` §7.2, measured `1.78e-15` against the package);
this claim exists so that a future storage change cannot break the population
layer silently.

---

## 6. The consensus / heterogeneity split

**Claim 6.** The population second moment of transported forms decomposes
exactly as `V = Q^C + Q^H`; the plug-in estimator of that decomposition is
biased by within-subject sampling noise; the cross-fitted estimator is unbiased
and may be indefinite; and the subject loadings on the heterogeneity modes are
read off the Gram of §5.

### 6.1 The identity

With `Y_i ∈ R^D` the svec'd transported forms of subject `i` and
`\bar Y = N^{-1}Σ_i Y_i`,

\[
V=\frac1N\sum_i Y_iY_i^\top
=\underbrace{\bar Y\bar Y^\top}_{Q^{C}}
+\underbrace{\frac1N\sum_i (Y_i-\bar Y)(Y_i-\bar Y)^\top}_{Q^{H}} ,
\]

mean square equals squared mean plus variance, with the `1/N` divisor that
makes the identity exact.

**Relation to the archived sketch — a substitution, not a specialization.**
`design/archive/searchlight-conversation-ledger.md:982-1053` decomposes
`Q_A^W = E[B_i^⊤AB_i] = \bar B^⊤A\bar B + E[U_i^⊤AU_i]`, which is quadratic in
the subject *patterns* `B_i` and lives in `q × q`. The identity above is
quadratic in the subject *forms* and lives in `D × D`. It is the **same
algebraic identity applied one level up**, and the objects genuinely differ —
the contract's `Q^C = \bar Y\bar Y^⊤` is always rank one, while the sketch's
`\bar B^⊤A\bar B` is not. The substitution is *forced*: node transport gives a
common node space but never a common feature space, so `\bar B` does not exist
at a group node (this is the same fact that makes a group-node common mode
undefined, §8). What is inherited verbatim is the estimator recipe of §6.2.
`Q^C` is population consensus and `Q^H` reliable heterogeneity, as in the
sketch.
Measured (§P6.a): `max |V − (Q^C + Q^H)| = 4.44e-16`; `tr V = 70.101024857210 =
15.261496149356 + 54.839528707855`; consensus share `0.217707176`.

**Naming hazard.** `Q^H` here uses the `1/N` divisor, which the exact identity
requires; the *estimator* discussed in §6.2 uses `1/(N−1)`. They are two
different objects sharing one symbol. An implementation must name them apart
(`Q_H_identity` versus `Q_H_plugin`) rather than trusting the divisor to be
inferred from context.

### 6.2 The identity is exact; the plug-in estimator is not

Under `Y_{ir} = μ + U_i + ε_{ir}` with independent partitions `r`, the plug-in
`Q^H` estimates `Σ_B + Σ_ε/R` — within-subject sampling noise is booked as
between-subject heterogeneity. The cross-fitted split avoids it exactly as
`effect-form-v1` §8 avoids the `Ξ_r K Ξ_r^⊤` term:

\[
\widehat Q^{W}=\frac1{2N}\sum_i\big(Y_{i1}Y_{i2}^\top+Y_{i2}Y_{i1}^\top\big),
\qquad
\widehat Q^{C}=\frac1{2N(N-1)}\sum_{i\neq i'}\big(\bar Y_i\bar Y_{i'}^\top+\bar Y_{i'}\bar Y_i^\top\big),
\]
\[
\widehat Q^{H}=\widehat Q^{W}-\widehat Q^{C}.
\]

`E[\widehat Q^W] = μμ^⊤ + Σ_B` because the two partitions are independent;
`E[\widehat Q^C] = μμ^⊤` because distinct subjects are independent; hence
`E[\widehat Q^H] = Σ_B`. **Cross-partition products kill within-subject noise;
cross-participant products kill the consensus-squared term.**

Measured (§P6.b/c), 2 000 Monte Carlo replications, `N = 12`, `D = 21`,
`tr Σ_B = 10.823228`, `tr Σ_ε = 13.440000` (halved to `6.72` by averaging two
runs):

| estimator | mean `tr \widehat Q^H` | bias |
|---|---|---|
| plug-in (unbiased divisor) | `+17.610196` | `+6.786968` (predicted `+6.72`) |
| cross-fitted | `+10.883837` | `+0.060609` (MC se `0.0824`) |

The plug-in inflates the heterogeneity trace by **+62.7 %**. The cross-fitted
estimator's bias is within one Monte Carlo standard error of zero.

**`\widehat Q^H` is routinely indefinite**: measured in **100 %** of the 2 000
replications, worst minimum eigenvalue `−2.079761`. This is the population-level
face of `conservative-geometry-v1` §6 and it is reported as-is.

### 6.3 Subject loadings

The eigenvectors of the centered `N × N` Gram are the subject loadings on the
heterogeneity modes. Measured (§P6.d) with a planted rank-2 heterogeneity on
the **plug-in** Gram: spectrum `11.5805, 4.9475, 1.4569, 1.0067`, and `R²` of
the Gram modes on the planted loading plane `0.958967` (mode 1), `0.932292`
(mode 2), `0.011605` (mode 3, noise). The rank-2 planted structure is recovered
in the top two modes and nowhere else.

**Which Gram an implementation may use is not free.** The plug-in Gram recovers
the *directions* well — noise inflates the spectrum roughly isotropically, so
the leading modes survive — but its *eigenvalues* are inflated, so it may not
carry a variance-explained number. Measured on the cross-fitted Gram of §6.4
the same recovery is `0.980363` and `0.938227`. **Normative: loading directions
may be read from either Gram provided the source is named; any eigenvalue,
spectrum, variance-explained figure or `n_eff` must come from the cross-fitted
Gram.**

### 6.4 The cross-fitted subject Gram — one `N × N` matrix behind §5, §6 and §7

§5 gives a computational trick and §6.2 gives an unbiased estimator; the object
that is *both* is the **cross-fitted subject Gram**. For independent partitions
`A, B` of the same data,

\[
\widehat\Gamma_{ii'}=\tfrac12\big(\langle Z_i^{(A)},Z_{i'}^{(B)}\rangle
+\langle Z_i^{(B)},Z_{i'}^{(A)}\rangle\big).
\]

`E[\widehat\Gamma] = \Gamma_{\text{true}}`, the Gram of the *noiseless* subject
vectors — **on the diagonal too**, because the two partitions are independent.
It is still `N × N`, still node-additive, and it is the matrix a population plan
feeds to §5.

Measured (`population-geometry-split.R` §P6.f), 2 000 replications:

| quantity | cross-fitted | plug-in |
|---|---|---|
| mean diagonal bias vs the true Gram | `+0.038398` (MC se `0.0537`) | `+6.750411` (predicted `+6.7200`) |
| mean centered trace vs planted `tr Σ_B = 10.8232` | `+10.8071` (sd `3.5754`) | `+17.5327` (sd `3.5188`) |
| negative eigenvalues, one draw | `5` of 12 | `1` |
| `R²` of mode 1 on the planted loading plane | `0.980363` | `0.958967` |

Two further measured facts fix the unification.

1. **`mean(diag \widehat\Gamma)` *is* `V^W` of §7.1** — measured `0e+00`, an
   algebraic identity, not an approximation. `mean(offdiag \widehat\Gamma)` is
   the cross-participant consensus energy. **So the consensus share `R(P)` of
   §7 is the off-diagonal-versus-diagonal contrast of this one matrix**, and
   `η_transport`, the heterogeneity spectrum and the subject loadings are all
   read from a single `N × N` object per transport.
2. **It is node-additive** — measured `1.42e-14` — so it accumulates while
   streaming group nodes, exactly as the plug-in Gram does (§5).

`\widehat\Gamma` is **indefinite by construction** (measured 5 negative
eigenvalues of 12 on one draw), and a single draw of its centered trace is not
an estimate of `tr Σ_B`: measured sd `3.5754` against a mean of `10.8071`. Both
facts are reported, never repaired.

**Normative.** A population plan computes `\widehat\Gamma` cross-fitted, records
which partitions formed `A` and `B`, and refuses to report any spectrum,
variance-explained figure or `n_eff` from the plug-in Gram. The plug-in Gram
remains admissible for **directions only**, named as such (§6.3).

### 6.5 Latent-layer discipline (mirrors gap G7)

`Q^C` and `Q^H` as *estimated above* are signed objects on the estimation
layer. Any statement of the form "the consensus mode is a geometry", any
variance-explained percentage, any `n_eff` over modes, and any nonnegative
spectrum requires a projection, and therefore lives on the latent PSD
descriptive layer of `conservative-geometry-v1` §6.

Measured (§P6.e) for one cross-fitted `\widehat Q^H`: 10 of 21 eigenvalues
negative, minimum `−1.348081`; signed trace `+14.664075533`, PSD-clipped trace
`+19.801714850`, **moved mass `5.137639317` = +35.04 %** of the signed trace.

**Normative.** (a) The projection is drawn from a closed named set — at
minimum `eigenvalue_clip` and `per_node_clamp`, which are *different operators
moving different mass*; (b) its name enters plan identity; (c) the receipt
records the moved mass per projection kind; (d) no fraction, cumulative curve
or effective-count is computed on the signed layer. This is
`conservative-geometry-v1` §6 verbatim, applied one level up.

---

## 7. `η_transport` and the required diagnostics slice

**Claim 7.** `η_transport` is a cross-fitted, held-out difference of consensus
shares between two transports. It may be negative, it is reported as-is, and it
is uninterpretable without the coverage diagnostics reported beside it.

### 7.1 Definition

For a transport `P`, a held-out data set, and two independent partitions
`A, B` of that held-out data (runs, sessions or tasks — *not* the partitions
used to build `P`), with `Z_i^{(r)}` the concatenated svec'd transported forms
over the `m` group nodes (the sink excluded):

\[
V^{W}(P)=\frac1N\sum_i\big\langle Z_i^{(A)},Z_i^{(B)}\big\rangle,
\qquad
V^{C}(P)=\frac1{N(N-1)}\sum_{i\neq i'}\big\langle \bar Z_i,\bar Z_{i'}\big\rangle,
\qquad
R(P)=\frac{V^{C}(P)}{V^{W}(P)} .
\]

`V^W` is cross-partition (unbiased for `‖μ‖² + tr Σ_B`), `V^C` is
cross-participant (unbiased for `‖μ‖²`), and `R(P)` is the **consensus share of
reproducible energy** — dimensionless, and invariant to any global rescaling of
`P`. Then

\[
\eta_{\mathrm{transport}} = R(P^{F}) - R(P^{A}) .
\]

Because both terms are cross-fitted, either share may fall outside `[0,1]` and
`η` may be negative. `V^W` is itself a cross-partition product and can be zero
or negative when nothing reproduces; the share is then **`NA`, not a large
number**. The plan requires `V^W` to be bounded away from zero by a declared
criterion before `R(P)` is formed, and reports `V^W` regardless so the reader
sees why.

### 7.2 Cross-fit provenance is required, not advisory

Measured (`population-transport-diagnostics.R` §P7.a/b): 16 subjects whose
functional bump is displaced from anatomy by `δ_i ∈ {−2…2}` group nodes, three
independent runs, run 1 used to estimate the functional transport and runs 2–3
held out.

| transport | `V^C` | `V^W` | share `R` | `η` vs `P^A` |
|---|---|---|---|---|
| anatomical `P^A` (identity) | `+3.9750` | `+11.1251` | `+0.357300` | — |
| functional `P^F`, fitted on run 1 | `+4.4698` | `+10.9665` | `+0.407587` | **`+0.050286`** |
| fractional (partial-volume) `P^F` | `+4.1714` | `+9.5427` | `+0.437129` | `+0.079829` |
| oracle `P` (true shifts) | `+5.7926` | `+11.0164` | `+0.525813` | `+0.168512` (ceiling) |
| **circular** `P^F`, fitted on run 2 | `+5.6775` | `+11.0061` | `+0.515852` | **`+0.158552`** |

The circular transport — the *same estimator*, differing only in that it was
fitted on a partition that is also used for evaluation — reports **3.15×** the
honest gain and lands essentially at the oracle ceiling. In a result object the
two numbers are indistinguishable. Hence: **`cross_fit` provenance is a
required field (§1.2), and the plan refuses to evaluate `η` on any partition
named there.**

### 7.3 Negative `η` is a result

Measured (§P7.c): a null transport (random shifts) gives `η = −0.019578`; over
200 random null transports the mean `η` is `−0.048047` and **89.5 %** are
negative.

**The same 200 draws are the reference distribution `η` needs, and reporting
them is required.** Measured: null `sd = 0.036059`, `q95 = +0.012540`,
`max = +0.043867`; the honest `η = +0.050286` **exceeds all 200 null draws**.
That is what makes `+0.05` mean something. A bare `η` printed without its null
band is exactly the uninterpretable headline §7.4 warns about, so the null
band — obtained by permuting or randomizing the transport, with the method
recorded — joins §7.5's required slice. Clipping `η` at zero would convert a demonstrated alignment failure
into a null result, repeating the `conservative-geometry-v1` §6 error one level
up. `η` is reported signed, and both `V^C` and `V^W` are reported separately for
each transport so the ratio is auditable rather than being the only number on
the page.

### 7.4 `η` alone is not interpretable

A transport can raise the consensus share by **discarding** the nodes that
disagree. Measured (§P7.d): a transport that keeps only the two most central
group nodes and sinks everything else reports `η = +0.166938` — **3.32×** the
honest functional gain — while sending **83.3 %** of native territory to the
sink.

**Normative: `η_transport` may not be printed, plotted or returned without the
diagnostics of §7.5 in the same object.**

### 7.5 The required diagnostics slice

Slice 2 (and any `η` report) must carry all six, with these definitions.

| diagnostic | definition | data-free? |
|---|---|---|
| `V_C`, `V_W` per transport | §7.1, reported separately, not only as their ratio | no |
| `eta_null_band` | the distribution of `η` under randomized transports (§7.3), with the randomization method recorded; report `sd`, `q95` and the honest `η`'s rank within it | no |
| `eta_transport` | `R(P^F) − R(P^A)`, signed | no |
| `displacement` | per native row, renormalize the group-node part of the row to sum 1, take `‖center(x) − Σ_j \tilde p_{xj}\,center(j)‖`; report median, p90, max and the row-mass-weighted mean. Rows with zero group mass are **excluded and counted separately** | yes |
| `entropy` / `perplexity` | per native row, `H_x = −Σ_j \tilde p_{xj}\log \tilde p_{xj}` in nats over the renormalized group columns; `perplexity = exp(H_x)` is the effective number of group nodes a native node spreads into. All-sink rows are `NA`. **Both are reported as row means**, and `mean(exp H) ≠ exp(mean H)` — the two summaries below are not consistent with each other by design, and an implementation must label which it reports | yes |
| `sink` | **two numbers**: unmapped native *territory* `Σ_x μ_x P_{x⊥}/Σ_x μ_x` (a property of `P` alone, computable before any data), and sink *budget* `Σ_x P_{x⊥}c_{i,x}` per subject | territory yes, budget no |
| `group_node_subject_coverage` | number of subjects contributing nonzero mass to each group node; report the minimum and the count of nodes below a declared floor | yes |

Measured on the §7.2 fixture:

| transport | displacement (median / p90 / mass-wt mean) | entropy (nats) | perplexity | unmapped territory | all-sink rows | min subject coverage |
|---|---|---|---|---|---|---|
| `P^A` anatomical | `0.000 / 0.000 / 0.000` | `0.0000` | `1.0000` | `0.00 %` | `0` of 192 | `16` of 16 |
| `P^F` functional | `1.000 / 3.000 / 1.308` | `0.0000` | `1.0000` | `11.98 %` | `23` of 192 | `10` of 16 |
| `P^F` fractional | `1.000 / 3.000 / 1.309` | `0.8943` | `2.5515` | `13.18 %` | `0` of 192 | `16` of 16 |
| sink-heavy (adversarial) | `0.000 / 0.000 / 0.000` | `0.0000` | `1.0000` | `83.33 %` | `160` of 192 | `0` of 16, **10 nodes below 2 subjects** |

The adversarial transport is caught by three of the six diagnostics, and by
`η` **not at all** — `η` rises. Its `V^W` does fall from `+11.1251` to `+3.1049`
and its `V^C` from `+3.9750` to `+1.6277`, so the two components are themselves
visible traces; that is exactly why §7.3 requires them reported separately and
not only as their ratio. A hard-assignment transport has entropy `0` by
construction; a realistic partial-volume warp does not, which is why entropy is
reported rather than assumed.

**Normative.** A group node whose subject coverage is below a declared floor is
not a group estimate; the result marks it. **The floor itself is an open
maintainer decision (§14.3)** — the contract requires the number be reported and
the marking mechanism exist, and does not legislate the threshold.

---

## 8. Labelling of transported components (closes gap G10)

**Claim 8.** For a feature-additive (identity or diagonal) metric, the
transported *total* is the group node's own total, exactly. The transported
*coherent* is a ledger of native-node coherence carried to a group location and
is a different object from any group-node common mode — under every metric.

*Proof.* The total is linear in the frame weights, so with the group weight
vector `w^G_j = Σ_x P_{xj} w_x`,
`Σ_x P_{xj}\,B D(w_x)B^\top = B D(w^G_j)B^\top` — the group node's own total.
The coherent component is `a a^\top/(a^\top K^{-1}a)` with `a` the node's
*normalized* weights, identity-metric specialization `(Bw)(Bw)^\top/Σw`
(`conservative-geometry-v1` §4); it is a ratio, homogeneous of degree one but
not additive in `w`, so no such identity exists, and by §4 of that contract
node-local rank-one projections have no linear operator summing them. ∎

The proof of the total half uses the feature-additive branch of
`conservative-geometry-v1` §3 (identity or diagonal metric), which is also what
the oracle exercises. Under a **dense** metric the node total is
`B D(√w_x) Q D(√w_x) B^⊤`, whose entries carry `√(w_{xu}w_{xv})` and are
therefore *not* linear in `w`; so `Σ_x P_{xj}G_{i,x}` is not
`B D(√w^G_j) Q D(√w^G_j) B^⊤` in general. **A dense-metric transport does not
inherit the exact group-total identity**, and E2 must scope the certificate to
`metric_capabilities()$feature_additive`, exactly as `conservative-geometry-v1`
scopes its own Claim 3.

Measured (§P8) on a fixture where a group frame genuinely exists (20 template
voxels, radius-1 conservative searchlights as native nodes, four contiguous
5-voxel regions as group nodes, so `w^G_j` is well-defined and itself
conservative, `max |colSums − 1| = 0`):

| group node | `|transported total − group-node own total|` | ledger identity `|coh + cfg − total|` | transported vs group-node coherent |
|---|---|---|---|
| 1 | `4.44e-16` | `4.44e-16` | **201.50 %** (`tr` `2.4287` vs `0.7920`) |
| 2 | `4.44e-16` | `2.22e-16` | **28.84 %** (`5.2133` vs `3.4269`) |
| 3 | `1.78e-15` | `3.55e-15` | **104.75 %** (`6.3050` vs `2.7959`) |
| 4 | `4.44e-16` | `8.88e-16` | **79.96 %** (`5.5130` vs `3.5528`) |

Two things are proved at once: the additive decomposition **does** survive
transport exactly (`conservative-geometry-v1` §7.6, "what does survive"), and
the coherent part **is not** the group node's coherent part even where the
latter is computable. In `crossform` it is not computable at all — transport
maps *nodes*, there is no group frame and no group-node `a` — so the group
quantity these names would suggest does not exist.

### 8.1 The names (normative; E2 implements)

| field | meaning |
|---|---|
| `transported_total` | `Σ_x P_{xj} c^{tot}_{i,x}` — a genuine group-node total **for that subject** (proof above). The induced group weight `w^G_j = Σ_x P_{xj}w_x` is built from subject `i`'s own native weights in subject `i`'s own feature space, so there is no one group-node frame across subjects — the same reason there is no group-node `a` |
| `native_coherent_ledger` | `Σ_x P_{xj} c^{coh}_{i,x}` — native-node coherent evidence carried to group node `j` |
| `native_configuration_ledger` | `Σ_x P_{xj} c^{cfg}_{i,x}` |
| `sink_budget` | the `⊥` column of the transported ledger, always in budget units |
| `native_coherence_fraction` | permitted **only** on the latent PSD layer, labelled with the native frame family identity |

**Forbidden:** the bare names `coherent` and `configuration` on a transported
result, and any coherence fraction taken against a summed coherent denominator
on the estimation layer (`conservative-geometry-v1` §4 normative 2 — such a sum
is not a global quantity, so the fraction has no denominator).

**Required print line.** Every transported component view prints the native
frame family it is a ledger of (family, scale, normalization) and the transport
identity (id, `semantics`, `provenance$kind`, cross-fit status). A reader must
never have to infer which frame a "coherent" number belongs to.

**Required documentation sentence**, so that the prohibition does not invite
the opposite error: *`native_coherent_ledger + native_configuration_ledger =
transported_total` holds exactly; what fails is reading the coherent ledger as a
group-node common mode.*

---

## 9. The interface: what is required, what is refused

**Claim 9.** WS-E consumes a conservative result if and only if it supplies the
eleven items below, and supplies none of the four things `crossform` refuses to
compute.

### 9.1 Required of the conservative result (cross-references §7 of `conservative-geometry-v1`)

The "supplied by" column separates two things the first draft of this table
conflated: code that produces the item, versus a contract that merely legislates
it. Only the first can be consumed by a plan.

| item | source | supplied today by |
|---|---|---|
| `per_row_family`, `per_row_scale`, `per_row_center`, `per_row_label`, `per_row_alpha` | §7.1, gap G11 | **code — delivered by D2.** `frame_family(..., alpha, normalization = "conservative")` carries `measurement` (`"<family>::<node>"`), `family`, `node`, `scale`, `center`, `alpha` one row per measurement (`.frame_family_member_index()`, `R/frame.R`) |
| `packed_codec_frobenius` | §7.2 | code (`crossform:::.svec_symmetric`) |
| `normalization_declared` | §1.2, §5.3 | code (`frame_conservation()$declared_normalization`) |
| `conservation_certificate` | §2, §5 | code, internal (`.metric_frame_conservation()`); must be surfaced |
| `composition_and_root` | §5.2.1, gap G5 | **nothing — owed by D6** |
| `signed_layer_declared` | §6 | **contract only.** No `R/` object declares a layer; `conservative-geometry-v1` §6 legislates it and `coherence_fraction_valid` is the single instance of the discipline |
| `latent_projection_named` | §6, gap G7 | **nothing — owed by D7** |

**This section is a stipulation, not a measurement.** §P9 transcribes the table
into executable predicates so the required set is machine-checkable; it does not
interrogate the package. On the current tree the `frame_family()` route fails
**2 of 11** items (`composition_and_root`, `latent_projection_named`) while the
bare `rbind()` + `additive_frame(..., "conservative")` route — which still works
and is still provenance-blind — fails **7 of 11**. **A population plan refuses a
frame that fails any row.** Two nuances an implementer needs: `scale` is
`NA_real_` when a member has no `$specification$radius`, and `center` is
`NA_character_` for region and whole-brain rows.

### 9.2 Refused by `crossform`, required as typed input

1. **Image registration** of any kind.
2. **Functional-transport learning** — estimating `P^F` from data. `crossform`
   accepts a `P^F` bearing cross-fit provenance and evaluates it (§7); it does
   not fit it.
3. **Resampling or interpolation of subject images.** Transport acts on nodes,
   not on voxels.
4. **Any group frame over group features**, and therefore any group-node
   coherent component (§8). Should a future group frame be defined, it takes
   the name `coherent` and the ledger names remain distinct.

---

## 10. What exists today, and what this contract requires

Verified by grep over `R/` and `NAMESPACE` at `elite-pass`: **no WS-E symbol
exists.** `location_transport`, `plan_population`, `population_form` and
`transport_object` return no matches.

| item | today | this contract requires |
|---|---|---|
| transport object | does not exist | sparse row-stochastic `P` with a required sink, `semantics`, `row_mass`, `native_index`/`group_index`, provenance with cross-fit, entering scientific plan identity (§1) |
| budget preservation certificate | does not exist | asserted per subject at fit time, `1e-12` (§2) |
| density semantics | does not exist | declared row-mass ratio `D(1/(P^⊤μ))P^⊤`, `NA` on zero mass, sink always in budget units (§1.3) |
| `plan_population()` | does not exist | OLS default via one QR; commuting class declared; order pinned in identity outside it (§3) |
| budget normalization | does not exist; gap G9 open | closed set `none` / `unit_budget` / `precision_weighted`, plan-identity field, `none` default (§4) |
| subject-Gram population covariance | does not exist | **cross-fitted** `N × N` Gram accumulated node-by-node; modes as `C^⊤u`; plug-in Gram admissible for directions only (§5, §6.4) |
| `symmetric_packed` codec | exists, Frobenius consistent | unchanged; it is the population storage format, and the `√2` is load-bearing (§5) |
| consensus / heterogeneity split | only the archived sketch (`searchlight-conversation-ledger.md:982-1053`) | cross-fitted `\widehat Q^W`, `\widehat Q^C`, `\widehat Q^H`; plug-in refused for inference (§6) |
| latent PSD layer | does not exist (gap G7) | named projection from a closed set; moved mass in the receipt (§6.5) |
| `η_transport` | does not exist | cross-fitted held-out difference of consensus shares; signed; never printed without §7.5 (§7) |
| transport diagnostics | do not exist | six diagnostics with the definitions of §7.5 |
| transported component names | do not exist; gap G10 open | `transported_total`, `native_coherent_ledger`, `native_configuration_ledger`, `sink_budget`; bare `coherent`/`configuration` forbidden (§8) |
| per-row frame metadata | **delivered (D2)** on the `frame_family()` route: `family`, `node`, `scale`, `center`, `alpha` per row; a single compiled frame's `$index` is still `measurement` only, and the metric fold still drops it (`conservative-geometry-v1` §5.3) | required input; plan refuses without it (§9.1) |
| `measurement_bridge()` | exists but is **internal** (not in `NAMESPACE`; its own examples call it as `crossform:::measurement_bridge()`) — a **feature-space** bridge | unchanged; it is *not* node transport and must not be confused with it |
| per-subject `SE` for `T_i` | does not exist (gap G8) | needed by `precision_weighted` and by the `unit_budget` refusal rule; both gated (§4.3, §4.5) |

---

## 11. Numerical contract

Measured on this contract's own oracle fixtures. These are the tolerances tests
must assert.

| law | tolerance |
|---|---|
| `Σ_j (P^⊤c)_j = Σ_x c_x` (budget preservation, incl. sink) | `1e-12 · Σ_x|c_x|` (measured `≤ 7.1e-15` over 500 draws) |
| `Σ_{j≤m}(P^⊤c)_j = Σ_x c_x − sink` | `1e-12 · Σ_x|c_x|` (measured `4.4e-16`) |
| density equals the fixed linear map `D(1/(P^⊤μ))P^⊤` | `1e-12` (measured `8.3e-17`) |
| commutation of query / transport / OLS, common weights | `1e-12` absolute (measured `8.9e-16`) |
| commutation under GLS with a subject-constant `Ω` | `1e-12` absolute (measured `1.3e-15`) |
| non-commutation under per-node weights | asserted as `> 1 %` relative, **never** as an equality (measured `25.2 %`) |
| non-commutation under per-coordinate weights | asserted as `> 1 %` relative, **never** as an equality (measured `25.5 %`) |
| each budget normalization conserves against its own total | `1e-12` absolute (measured `≤ 1.1e-16`) |
| `precision_weighted` with equal precisions equals `none` | `1e-12` (measured `0`) |
| per-subject scalar commutes with transport | `1e-12` (measured `1.1e-16`, over a real sparse transport) |
| `⟨svec A, svec B⟩ = ⟨A,B⟩_F` (`√2` codec) | `1e-12` relative (measured `3.1e-15`) |
| failure of the naive packing | asserted as `> 1 %` relative (measured `2.96e+00`) |
| subject-Gram spectrum vs `D × D` covariance spectrum | `1e-12` absolute (measured `1.1e-14`) |
| `rank(Σ) ≤ N − 1`; residual eigenvalues | `1e-12` (measured `7.4e-15`) |
| `Γ = Σ_j Γ_j` (node-additive Gram) | `1e-12` (measured `3.6e-15`) |
| `V = Q^C + Q^H` (exact identity, `1/N` divisor) | `1e-12` (measured `4.4e-16`) |
| `mean(diag \widehat\Gamma) = V^W` | `1e-12` (measured `0`) |
| cross-fitted Gram is node-additive | `1e-12` (measured `1.4e-14`) |
| cross-fitted Gram diagonal unbiased for the true Gram | asserted as `< 1` MC se (measured `+0.038` vs se `0.054`, 2 000 reps) |
| cross-fitted centered Gram trace unbiased for `tr Σ_B` | asserted as `< 1` MC se (measured `+10.8071` vs planted `10.8232`, sd `3.5754`) |
| plug-in centered Gram trace biased | asserted as `> 10 %` (measured `+17.5327` vs `10.8232`, `+62 %`) |
| cross-fitted `tr \widehat Q^H` unbiased for `tr Σ_B` | asserted as `< 1` Monte Carlo se (measured `+0.061` vs se `0.082`, 2 000 reps) |
| plug-in `tr \widehat Q^H` biased | asserted as `> 10 %` inflation (measured `+62.7 %`, predicted bias `+6.72` vs measured `+6.79`) |
| `\widehat Q^H` indefiniteness | asserted as "occurs", **never** clipped (measured 100 % of 2 000 reps) |
| transported total equals the group node's own total, **feature-additive metric only** (§8) | `1e-12` (measured `1.8e-15`) |
| `native_coherent_ledger + native_configuration_ledger = transported_total` | `1e-12` (measured `3.6e-15`) |
| transported coherent vs group-node coherent | asserted as `> 1 %` relative, **never** as an equality (measured `28.8 %` – `201.5 %`) |
| circular vs cross-fitted `η` | asserted as "circular is materially larger" (measured `3.15×`) |
| `η` for a null transport | asserted as "may be negative", **never** clipped (measured `−0.048` mean, 89.5 % negative over 200 draws) |
| honest `η` against its null band | asserted as a rank, not a p-value (measured `+0.050286` exceeding 200 of 200 nulls; null `sd 0.036059`, `q95 +0.012540`) |

---

## 12. Test and oracle index

No `crossform` test covers any claim in this document, because no WS-E code
exists. Every row's `test` column is an obligation on E2 (implementation) or E3
(slice 1), not a citation.

| claim | statement | oracle | test owed |
|---|---|---|---|
| 1a | a well-formed transport is nonnegative, row-stochastic including the sink, with declared semantics, row mass and provenance; four failure modes detected | `population-transport-contract.R` §P1.a | E2 |
| 1b | budget and density are two linear maps from one `P`; density is a declared row-mass ratio; the sink is always in budget units | §P1.b | E2 |
| 2 | transported group total = native total − sink, exactly, for signed ledgers; deleting the sink loses mass silently | §P2 | E2 |
| 3a | query, transport and the fit commute under a subject-constant weight operator (OLS and common-`Ω` GLS) | §P3.a | E2 |
| 3b | per-node weights break aggregation↔fit commutation | §P3.b | E2 |
| 3c | per-coordinate weights break query↔fit commutation | §P3.c | E2 |
| 4 | `none` / `unit_budget` / `precision_weighted` all conserve, give different estimands, and can differ in argmax; `unit_budget`'s two failure modes | §P4 | E2 |
| 5a | the `√2` codec is Frobenius consistent; the naive packing is not and changes the spectrum | `population-geometry-split.R` §P5.a | E2 |
| 5b | `rank(Σ) ≤ N−1`; the `N × N` Gram recovers the spectrum exactly | §P5.b | E2 |
| 5c | modes recovered as `C^⊤u`; the `D × D` covariance is never formed | §P5.c | E2 |
| 5d | the Gram is additive over group nodes (streaming) | §P5.d | E2 |
| 6a | `V = Q^C + Q^H` is an exact identity | `population-geometry-split.R` §P6.a | E2 |
| 6b | the plug-in split is biased by within-subject noise | §P6.b | E2 |
| 6c | the cross-fitted split is unbiased for `Σ_B` and is routinely indefinite | §P6.c | E2 |
| 6d | subject loadings are read off the Gram modes | §P6.d | E3 |
| 6e | a PSD projection moves mass; the amount must be recorded | §P6.e | E2 |
| 6f | the cross-fitted Gram is unbiased for the noiseless subject Gram, node-additive, indefinite by construction, and its diagonal/off-diagonal contrast is the `V^W`/`V^C` split | §P6.f | E2 |
| 7a | the cross-fitted consensus share and `η_transport` | `population-transport-diagnostics.R` §P7.a | E3 |
| 7b | a circular transport inflates `η` by `3.15×`; cross-fit provenance is required | §P7.b | E2 (refusal), E3 (number) |
| 7c | `η` may be negative and is reported as-is; the null band is what calibrates a positive `η` | §P7.c | E2 |
| 7d | a sink-heavy transport buys consensus by discarding nodes; `η` alone is uninterpretable | §P7.d | E3 |
| 7e | the six diagnostics and their definitions | §P7.e | E2 |
| 8 | transported total = group-node own total exactly; transported coherent ≠ group-node coherent; the ledger identity survives | §P8 | E2 |
| 9 | the eleven readiness items; `frame_family()` fails 2, the bare `rbind` route fails 7 | §P9 — **a stipulation, not a measurement** (§9.1) | E2 |

**These oracles are wired to nothing.** `conservative-geometry-v1` §10 records
the same defect for its three scripts and gap G12 assigns the fix to D2. The
same fix applies here: **E2 must promote claims 2, 3a–3c, 4, 5a–5d, 6a–6f and 8
into a `test-population-form-contract.R`**, rather than leaving them protected
only by a human remembering to run a script. The oracles should remain readable
derivations.

---

## 13. Self-review against `.planning/2026-08-17-feedback-assessment.md`

Part 2's ten corrections, each checked against this document.

| # | correction | honoured where | notes |
|---|---|---|---|
| 1 | the conservative total is a smoothed univariate ledger; the informative object is the coherence spectrum, not per-scale energy | §8 (the transported total is exactly the group node's own total, so the smoothing picture carries through transport — **on the feature-additive branch only**, scoped in §8); §7.1 (the population headline is a *share*, `R(P)`, not an energy) | No panel in this contract reports transported energy per scale. `V^C`/`V^W` are reported as inputs to a share, never as findings on their own — the same discipline `conservative-geometry-v1` §3.1 imposes on α. |
| 2 | coherent budgets are frame-relative; fractions only in the latent layer, labelled with frame identity | §8.1 (`native_coherence_fraction` is latent-layer only and carries the native family identity; fractions against a summed coherent denominator are forbidden) | **Fixed during review**: an earlier draft of §8.1 listed a plain `coherence_fraction` field. |
| 3 | conservation fails for dense metrics; `composition` must be explicit, with the root pinned | §9.1 (`composition_and_root` is a required input item, and a frame lacking it is refused) | This contract inherits rather than re-derives; the population layer never composes a metric itself. |
| 4 | budget vs density is a declared choice; the sink is required, not optional | §1.1, §1.3 | Extended: density is given an exact definition via a declared row mass, and the sink is required to be *materialized* even when zero. |
| 5 | transported coherent is a ledger of native-node coherence | §8 in full | Extended with a measured demonstration on a fixture where both objects exist (`28.8 %` – `201.5 %` apart), plus the positive statement that the additive identity survives. |
| 6 | signed vs latent layers; transport and OLS are linear so unbiasedness survives; PSD clipping never silent | §2 (signed-sum law), §3 (linearity), §6.5 (named projections, moved mass measured at `+35.04 %`) | **Caught during review**: §6's `Q^C`/`Q^H` are *quadratic*, so unbiasedness does **not** follow from the linearity argument. §6.2 supplies the cross-fitted estimator instead and measures the plug-in bias. This is the one place where correction 6's licence does not extend, and the contract says so. |
| 7 | `η_transport` cross-fitted, may be negative, reported as-is; transport learning stays outside | §7.1–7.3, §9.2 | Extended: §7.4 shows `η` is gameable by sinking nodes, so §7.5 makes the diagnostics mandatory beside it. |
| 8 | **Part 2 correction 8**: the Gram trick needs a Frobenius-consistent packed codec, and `symmetric_packed` (√2 off-diagonals) already is | §5 ("the √2 is not cosmetic"); §10 codec row | Extended: the naive packing is measured to be off by `2.96e+00` — an `O(1)` error — and to change the population spectrum (`18.846` vs `24.734`). The claim exists so a future storage change cannot break the population layer silently. **A first draft of this table put Part *1* correction 8 in this row and left Part 2's correction 8 unaudited.** |
| 8′ | **Part 1 correction 8** (audited here because §3 rests on it): per-node weighting breaks commutation; OLS default, weighted variants explicit | §3.1–3.3 | **Refined**: per-*coordinate* weights break query↔fit; per-*node* weights break aggregation↔fit. Conflating them would let a per-group-node weighted fit ship with a false commutation certificate (§3.3 item 5). |
| 9 | frame families, per-row metadata, `contribution(by=)`, population code do not exist | §9.1, §10 (grep-verified), §12 | **Partly overtaken by D2**, which landed `frame_family()` with per-row `family`/`node`/`scale`/`center`/`alpha` while this contract was being written; `conservative-geometry-v1` §8 now records those rows as delivered. §9.1 and §10 are written against the current tree: `frame_family()` fails 2 of 11 readiness items, the bare `rbind` route still fails 7. The population half of correction 9 stands — no WS-E symbol exists (grep-verified §10). |
| 10 | slice 1 = region-level nodes on six Haxby subjects; slice 2 = normalized dataset with searchlight transport | §7 is written so slice 1 needs only §§1–6 (region-level transport is `P` with trivial rows and an empty sink), and only slice 2 needs `η` and the diagnostics | This contract does not choose the slice-2 dataset; that is E12's ticket. |

Four corrections applied to this document in the course of its fresh-context
review are marked above (rows 2, 6, 8 and 9), and three more are marked in the
body: §4.1's group functional (a single weighted-mean template does not cover
`unit_budget`), §4.4's vacuous commutation citation, and §6.4's addition of the
cross-fitted Gram — without which §5 would have specified the very plug-in
object §6.2 refutes. Row 6 remains the substantive one: the population second moment
is quadratic in the estimates, so the population layer is **not** purely linear,
and the contract now names the exact point where `conservative-geometry-v1` §6's
"transport and OLS are linear, so unbiasedness survives" stops applying.

---

## 14. Open maintainer decisions

Flagged, not decided, because each needs either a missing estimator or a
scientific judgement outside this contract's scope.

1. **§4.5 — `precision_weighted` is gated on D8.** Until a cross-node
   `contrast_energy` sampling route exists (gap G8), `π_i` can only be supplied
   externally. Alternative: ship `precision_weighted` with a
   `precision_source = "external"` requirement and no internal path.
2. **§4.3 — the `unit_budget` refusal threshold.** "`T_i` bounded away from
   zero" needs a per-subject standard error, which does not exist. Interim
   options: refuse `unit_budget` entirely until D8; accept it with a declared
   absolute floor; or accept it descriptively with a receipt warning. The
   contract requires the refusal mechanism, not the constant.
3. **§7.5 — the group-node subject-coverage floor.** A group node fed by two of
   twenty subjects is not a group estimate, but the threshold is a study-design
   judgement. The contract requires the number to be reported and the marking
   mechanism to exist.
4. **§7.1 — whether `η` or the two shares are the headline.** This contract
   requires both `V^C` and `V^W` to be reported per transport so the ratio is
   auditable. Whether `η` appears in the printed summary line at all is a
   presentation decision for E3.
5. **§7.3/§7.5 — how the `η` null band is generated.** The contract requires a
   null band and requires the randomization method to be recorded. Which
   randomization is right is a scientific choice: randomizing the transport
   (what §P7.c does), permuting subject labels, or permuting the held-out
   partition assignment test different nulls. E3 picks one per slice and
   records it; the contract does not legislate it.
6. **§1.3 — the default `row_mass`.** `μ ≡ 1` (native node count) is the stated
   default. For frames whose nodes differ greatly in size, the node's own frame
   mass `Σ_v w_{xv}` is arguably the better default. This contract fixes the
   *mechanism* (a declared positive vector) and the current default, and leaves
   the default to D2's per-row metadata work, which is what would make the
   alternative computable.
