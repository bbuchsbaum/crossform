# Model-coordinate geometry contract — reduced-rank model fitting as a lowering

Status: normative architecture contract

Contract version: `model-coordinate-v1`

Date: 2026-09-04

Last amended: 2026-09-04 (fresh-context review, §15; implementation amendments, §15.1a; §11 status)

Epic: `bd-01M1P2X7BB4TP1D2JXA2HAK3EA` (mote, tag `model-coordinate`)

This document freezes the semantics of **model-coordinate geometry**: how a
family of model representational dissimilarity matrices (RDMs), or model
feature matrices, becomes a declared basis of the effect axis; how the
cross-generalized geometry is read in that basis; which readings stay on the
signed estimation layer and which are latent projections; and what a
reduced-rank fit of the geometry to the model family does and does not test.
It is normative for every child of the epic above. The independent executable
court is `design/oracles/model-coordinate-geometry.R`; it uses only base
matrix algebra and never loads crossform. A child ticket compares the public
package route against it in `tests/testthat/test-model-coordinate-contract.R`.
This document does not itself change production code or the supported API.

It consumes `effect-form-v1` (the form, its codecs, the baseline-invariance
diagnosis of §7 and the fixedness premise **(A3)** of §8),
`conservative-geometry-v1` §6 and `population-form-v1` §6.5 (the signed
estimation layer versus the latent PSD descriptive layer), and the
partition-disjoint training discipline of `metric_training_policy()`.

## Why this contract exists

The proposal it answers is a two-sided reduced-rank RSA: a latent signal
geometry \(G = T C T^\top\) with a rank-limited factor \(T_i\) per model, a
positive semidefinite (PSD) \(C\) of rank at most \(s\), and full-rank
measurement noise kept apart from the signal. Its mathematics is sound. Its
integration sketch was a bank of \(r(r+1)/2\) fixed bilinear queries followed
by a per-measurement eigensolver. That sketch is valid, but it is strictly
dominated by a construction the package already has: the model basis is an
**extractor**, the lowered relation is an ordinary relation, and the
compressed form \(Q^\top G_x Q\) is a complete geometry on which every
existing reader applies. This contract fixes that construction and the
discipline around it, so that the four hypothesis classes of the proposal
become four readings of one form rather than four algorithms.

## Notation

- \(n\) conditions; the effect space of the relation, ordered and signed by
  `effect_space()`. \(H = I_n - \mathbf 1\mathbf 1^\top / n\).
- \(B_r \in \mathbb R^{n \times p}\): the relation estimated within partition
  \(r\); \(K_x\): the fixed measurement operator at node \(x\); \(\Gamma\): the
  declared ordered pairing and reducer. The total geometry is
  \(G_x = \sum_{r,s} \Gamma_{rs} B_r K_x B_s^\top\) (`unification-v1` §2.2).
- \(D_i\): model \(i\)'s squared-Euclidean RDM over the same \(n\)
  conditions; \(K_i = -\tfrac12 H D_i H\) its centered Gram.
- \(T_i \in \mathbb R^{n \times r_i}\): the rank-revealing factor of \(K_i\);
  \(T = [T_1, \ldots, T_k]\); \(T = QR\) with \(Q^\top Q = I_q\),
  \(q = \operatorname{rank}(T)\); \(R\) is square and invertible exactly when
  \(T\) has full column rank (§1.5).
- \(S_x = Q^\top G_x Q \in \mathbb R^{q \times q}\): the **compressed form**,
  the geometry in model coordinates. \(P = QQ^\top\).
- \([S]_{+,s}\): the matrix keeping the \(s\) largest positive eigenvalues of
  a symmetric \(S\) and setting every other eigenvalue to zero.

---

## 1. The model basis is a declared effect map

**Claim 1.** A model family enters crossform as one data-independent value,
the **model basis**: an orthonormal, centered basis \(Q\) of the joint model
span together with the triangular factor \(R\) that recovers each model's own
coordinates, the per-model spectra, the requested and effective ranks, the
negative model mass dropped, the span fraction, and a signature. It is
constructed before any neural data is read and it enters plan identity.

### 1.1 Construction

For each model \(i\):

\[
K_i = -\tfrac12 H D_i H = U_i \Lambda_i U_i^\top,
\qquad
T_i = U_{i, r_i^\ast}\,\Lambda_{i, r_i^\ast}^{1/2},
\]

where \(r_i^\ast = \min(r_i, \operatorname{rank}_\tau K_i)\) and
\(\operatorname{rank}_\tau\) counts eigenvalues above a relative tolerance
\(\tau\) of the largest absolute eigenvalue. Then \(T = [T_1, \ldots, T_k]\) is
factored by a **rank-revealing** decomposition (thin SVD, or a pivoted QR with
the same tolerance) into \(Q\) and \(R\).

**Normative.** (a) Every truncation is rank-revealing; a thin, unpivoted QR of
\(T\) is forbidden. (b) The requested rank \(r_i\) and the effective rank
\(r_i^\ast\) are both recorded; clamping is reported, never silent. (c) When
model features \(F_i\) are supplied instead of an RDM, \(K_i = H F_i F_i^\top H\)
and the same rules apply; a learned readout \(T_i = F_i L_i\) is **not**
admitted by this version (§14). (d) Every decomposition runs in an explicit
orthonormal basis \(V\) of \(\mathbf 1^\perp\): \(K_i\) is formed as
\(-\tfrac12 V^\top D_i V\) (or \((V^\top F_i)(V^\top F_i)^\top\)), its
eigenvectors are lifted as \(V \tilde U_i\), and \(T\) is factored through
\(V^\top T\) so that \(Q = V \tilde U\). Eigenvalues within \(\tau\) of zero
count on neither side: they are not kept and are not dropped mass.

Measured (fresh-context review of M2, 2026-09-04): a direct
eigendecomposition of \(-\tfrac12 H D H\) leaves the eigenvector of a small
kept root \(\lambda_k\) uncentered by \(\epsilon\,\lambda_{\max}/\lambda_k\),
which is \(\epsilon\) times the *square* of the condition number of \(R\); on a
two-feature model with a \(10^{-4}\) scale ratio that reached
\(|Q^\top\mathbf 1| = 1.1\times10^{-6}\), and an ordinary model tripped the
centering assertion. Under rule (d) the same grid stays at
\(6\times10^{-16}\). Rule (d) exists because of this measurement.

Measured (oracle §O2, the trap): the literal factor
\(U_{i,r}\Lambda_{i,r}^{1/2}\) with \(r = 2\) on a rank-1 category RDM has a
zero column; a thin QR of it returns a second basis vector with
\(|Q^\top \mathbf 1| = 1.10\), and under a condition-independent baseline shift
the compressed form moves by \(8.12\). The same trap was reproduced through
the package route on the README fixture: a thin QR of the two-column literal
factor of the rank-1 category RDM gave \(Q^\top\mathbf 1 = 1.15\) on its
second column, and the first draft of this contract quoted a trace split
measured on that basis (§15, correction 1). Rule (a) exists because of this
measurement.

### 1.2 Admission

**Normative.** (a) A model RDM is admitted to the Gram construction only under
an explicit `distance = "squared_euclidean"` declaration; `rsa()`'s admission
rule (any symmetric zero-diagonal matrix, `R/views.R:385-492`) is a different
rule for a different reading and is not inherited. (b) Negative eigenvalues of
\(K_i\) are dropped and their total mass is recorded as
`dropped_negative_mass`; when that mass exceeds a declared share of the
positive mass the constructor signals an `effect_capability_refusal` with
capability `"euclidean_model_geometry"`, namespace `"model_coordinate"`,
reasons `"model_gram_indefinite"`, and the remedy of supplying a Euclidean
embedding or the features that generated the RDM. (c) Condition labels are
aligned to the relation's effect coordinates by name, exactly as `rsa()`
aligns them; unlabeled models are admitted only when their dimension equals
the effect dimension.

### 1.3 Centering is a property of the basis

`effect-form-v1` §7 requires a constructor to diagnose whether the final
query cancels additive left and right baselines, which holds only when both
marginals of the operator vanish. For the model basis every induced operator
is \(\operatorname{sym}(q_j q_k^\top)\), so both marginals vanish if and only
if \(Q^\top \mathbf 1 = 0\).

**Normative.** The constructor asserts \(Q^\top \mathbf 1 = 0\) to a fixed
tolerance (\(10^{-10}\sqrt n\)) and records `baseline_invariant = TRUE`.
Failure is the package's invariant error (`effect_invariant_error`), not a
refusal and not a contract error between two objects: by §1.1(d) a correctly
built basis cannot fail it, so a failure is a construction bug.

Measured (oracle §O2): a random baseline shift moves the raw form by
\(17.2\) and the compressed form by \(1.14\times10^{-13}\).

### 1.4 Identity

**Normative.** The basis signature is a SHA-256 over the ordered condition
coordinates, the model names, kinds and their declared distances, \(\tau\),
the requested and the effective ranks (so that two families with one
stacked factor but different per-model attributions are different objects),
and the unrounded \(Q\) and \(R\). It is carried into the
effect space of the lowered relation (§2.1) and therefore into every plan,
receipt, and population signature downstream. Two bases with equal
dimensions and labels but different signatures are different objects.

### 1.5 Identification of model coordinates

\(C\) is read through \(A = R C R^\top\). When \(T\) has full column rank,
\(R\) is square and invertible and \(C = R^{-1} A R^{-\top}\) is identified.
When model spans overlap, \(q < \sum_i r_i^\ast\), \(R\) is \(q \times
\sum_i r_i^\ast\), and \(C\) is identified only modulo the null space of
\(R\): the shared and block metrics and the diagonal weights \(c\) are then
not unique, although \(A\) and the fitted geometry are.

**Normative.** (a) This version admits only bases with
\(q = \sum_i r_i^\ast\). A basis whose stacked factor is rank-deficient
signals an `effect_capability_refusal` with capability
`"identified_model_coordinates"`, namespace `"model_coordinate"`, reasons
`"model_span_overlap"`, and the remedies: reduce the ranks; drop the
redundant model; or merge the overlapping models into one. (b) The basis
records `span_overlap_rank = \sum_i r_i^\ast - q` (zero when admitted) and
`basis_condition`, the condition number of \(R\), so that near-overlap is
visible before a fit is read. (c) A later version may admit overlap by
reporting \(A\) only, under its own name.

---

## 2. Lowering, not querying

**Claim 2.** The geometry in model coordinates is obtained by lowering the
relation, not by querying the form. The lowered relation
\(\tilde B_r = Q^\top B_r\) is an ordinary relation on a new effect space; the
ordinary plan and executor produce the complete \(q \times q\) form
\(S_x = Q^\top G_x Q\).

### 2.1 The lowered relation

crossform defines a relation as \(B = E Y\) with a declared extractor \(E\)
(`R/extractor.R:3-8`), and builds the identity extractor when betas are
supplied directly (`R/relation.R:299`). The model basis is therefore

- on the betas path, an `effect_extractor()` with map \(Q^\top\) and effect
  space "model coordinates of \(\{i\}\)" whose identity is the basis
  signature; and
- on the ingestion path, an `effect_map()` with weights \(Q^\top\) over the
  study's `condition_space()`, lowered through the design parameterization
  exactly as any other declared effect functional (`R/relation-plan.R:121`).

**Normative.** (a) The model-coordinate effect space carries the basis
signature in its identity, so that two subjects lowered through the same
basis share one effect-space signature and two different bases never do.
(b) The lowered form is a `complete_form` self form with the `symmetric`
capability, stored in the packed codec of width \(q(q+1)/2\), and it is
**not** `guaranteed_psd`. (c) Nothing on the lowered plan is a new estimand
kind: under a fixed metric, `plan_geometry()`, `materialize_geometry()`,
`geometry_spectrum()`, `latent_geometry()`, the sampling covariance, and
`plan_population()` consume it unchanged. Learned local metrics are the
exception stated in §2.3.

### 2.2 Commutation, and where it stops

For every stage that is linear on the effect axis, lowering commutes with the
estimator:

\[
Q^\top\Big(\sum_{r,s}\Gamma_{rs} B_r K_x B_s^\top\Big) Q
= \sum_{r,s}\Gamma_{rs}\,(Q^\top B_r)\,K_x\,(Q^\top B_s)^\top ,
\]

and likewise for the coherent component
\(\Gamma_{rs}(B_r w_x)(B_s w_x)^\top / a_x\) and hence for configuration.

Measured, package route (2026-09-04, working tree, README fixture: 4
conditions, 4 runs, 280 searchlights; extractor route against projecting the
materialized full form):

| component | max abs difference |
|---|---|
| total | \(8.9\times10^{-16}\) |
| coherent | \(8.9\times10^{-16}\) |
| configuration | \(7.2\times10^{-16}\) |

Measured, oracle §O1 (6 conditions, 4 partitions, weighted node): total
\(5.7\times10^{-14}\), coherent \(1.8\times10^{-15}\), configuration
\(2.8\times10^{-14}\).

**Where it stops.** `effect-form-v1` §6 forbids moving a query across a
nonlinear normalization or transform. Under `correlation()`, `cosine()`,
`fisher_z()`, or `rank_edges()`, the lowered and projected estimands differ.

**Normative.** Model-coordinate lowering is admitted only for plans whose
edge normalizer and transform are linear on the effect axis
(`inner_product()` and `covariance()` with the identity transform). Any other
plan signals an `effect_capability_refusal` with capability
`"model_coordinate_lowering"`, namespace `"model_coordinate"`, reasons
`"nonlinear_edge_transform"`, and the remedy of reading the full form and
projecting after the transform. The refusal is not a claim that
lowered-then-correlated is meaningless: a correlation of model-coordinate
patterns across partitions is a well-defined estimand, but it is a different
one, it has no name or identity in this version, and refusing it keeps it
from being produced silently under the model-coordinate name. It may be
admitted later under its own name.

**Where the gate lives.** Edge normalizers and transforms are declared on the
evidence stage plan (`.evidence_stage_plan()`, `R/evidence-task.R`, layer 3),
not on `plan_geometry()`. The stage plan does not see the relations, so the
gate is `.model_coordinate_lowering_gate()`, called from
`.new_evidence_task()`, the one constructor every evidence task (compiled,
reversed, measurement) passes through with both relations and the stage plan
in hand. Today the
normalizers `inner_product()`, `covariance()`, `correlation()`, `cosine()`,
`fisher_z()` and `rank_edges()` are unexported (`R/operations.R:42-111`) and
every caller uses the default, so the gate is dormant: it guards a path no
public entry can reach yet, and M3 tests it against the internal stage plan
until a normalizer is public.

### 2.3 Why not a query bank

The engine accepts a packed \(q(q+1)/2 \times V\) bank
(`R/compiler.R:62-67`; production use at `R/population-driver.R:99-106` and
`:891`), and one operator of the proposal's bank reproduced the projected
entry to \(4.2\times10^{-16}\). The bank is nevertheless the wrong primary
route, for reasons that are properties of the package rather than of taste:

1. A query-only result "cannot claim completeness or be upgraded by
   relabeling, even when a chosen set of queries happens to span the form"
   (`effect-form-v1` §3). The lowered form is complete.
2. The population layer requires identical effect-space signatures across
   subjects (`R/population-plan.R:232-251`). A basis compiled once from a
   shared condition space satisfies that; a bank produces views that the
   population form cannot pool.
3. Under a fixed metric the lowered relation is materialized by the same
   packed kernel as any relation. Under a learned local metric **neither**
   route reaches the complete form: the compiler refuses full materialization
   (`scheduled_metric_full_materialization_not_admitted`,
   `R/compiler.R:219-224`) and contracts one signed contrast per call
   (`scheduled_metric_requires_signed_contrast`, `:274-290`). Model-coordinate
   readings under a learned local metric are therefore limited to contrast
   energies \(c^\top S_x c\) and are outside this version (§14). The Gram
   route serves structured pair-difference queries only
   (`R/kernel.R:531-532`) and bears on neither route.
4. Baseline invariance is one assertion on \(Q\) (§1.3) rather than a
   diagnosis of \(r(r+1)/2\) operators.

**Normative.** The bank remains admissible as a secondary route for a
query-first plan that must not be re-lowered; it produces an `effect_view`
and is subject to every restriction of a query-only result.

---

## 3. The two layers

**Claim 3.** The compressed form is a fixed linear query of the
cross-generalized form and inherits its noise-unbiasedness; the rank-\(s\)
fit is a latent projection of it and inherits none.

### 3.1 The compressed form is signed and unbiased

\(Q\) is fixed before any neural data is read, so premise **(A3)** of
`effect-form-v1` §8 holds and, by its Corollary ("the same argument covers
every fixed pair query"), \(\mathbb E[S_x] = Q^\top G_x^\star Q\). Negative
eigenvalues of \(S_x\) are valid estimates and are never clipped on this
layer.

Measured (oracle §O7, 600 replications): under pure noise the mean
addressable energy \(\operatorname{tr} S_x\) is \(0.0049\) (se \(0.0281\)); with
a planted signal the bias of \(\operatorname{tr} S_x\) against
\(\operatorname{tr}(Q^\top G^\star Q) = 361.79\) is \(-0.877\) (se \(0.445\)).

### 3.2 The rank-\(s\) fit is a named latent projection

**Theorem 4 of the proposal**, restated in the package's terms. For the
Frobenius objective \(\|G_x - T C T^\top\|_F^2\) over PSD \(C\) of rank at most
\(s\),

\[
\|G_x - QAQ^\top\|_F^2
= \underbrace{\|G_x - Q S_x Q^\top\|_F^2}_{\text{constant in } A}
+ \|S_x - A\|_F^2,
\qquad
A^\star = [S_x]_{+,s},
\qquad
C^\star = R^{-1} A^\star R^{-\top}.
\]

Measured (oracle §O3): the split holds to \(10^{-10}\); the closed form
attains \(15{,}219.4\) against a best of \(57{,}504.6\) over 2 000 random
feasible \(C\) and \(15{,}220.1\) over 2 000 local perturbations of
\(A^\star\).

This is `latent_geometry()`'s eigenvalue truncation (`R/latent.R:655`) with
a rank budget added. It therefore inherits the three numbered decisions of
`R/latent.R:22-44` verbatim, the projection is named from a closed set,
clipping is never silent, and the signed source is not touched, and it
prints `"latent descriptive layer; not for inference"` (`R/latent.R:63`).

**Normative.** (a) The projection kind is `psd_projection` with a declared
`rank`; the rank enters the latent scientific identity. (b) The receipt
records moved mass in **two** parts that are reported separately:
`clipped_negative_mass`, the negative roots set to zero, and
`truncated_positive_mass`, the positive roots beyond rank \(s\) set to
zero. They are different reasons for moving mass and are never summed
into one number without both being visible. (c) No fraction, cumulative
curve, effective count, or "variance explained" is computed on \(S_x\)
itself.

Measured (oracle §O7): under pure noise the rank-1 latent energy
\(\operatorname{tr}[S_x]_{+,1}\) has mean \(0.6446\) (se \(0.0167\)), while the
signed addressable energy is centered on zero. That gap is the bias this
clause names.

### 3.3 No noise floor is subtracted

The proposal's §7 fits \(V = N + TWW^\top T^\top\) and subtracts an estimated
\(\sigma^2\) from the retained eigenvalues. That model describes a
within-partition second moment. crossform's form is a cross-partition
product: no within-partition term is ever formed, the noise floor never
enters \(G_x\), and negative roots are its sampling error rather than a
floor to remove.

**Normative.** No estimator in this contract subtracts a noise variance from
\(S_x\) or its spectrum. Noise-aware fitting is admitted only through the
sampling covariance of the compressed coordinates (§8).

---

## 4. Hypothesis classes are readings of one form

**Claim 4.** With \(T = QR\), every structure the proposal defines is a
function of \(S_x\) and of model-only matrices. One neural pass serves all
four.

| structure | \(C\) | what is needed from the neural side | what it tests |
|---|---|---|---|
| isotropic | \(\operatorname{blockdiag}(\alpha_i I)\), \(\alpha_i \ge 0\) | \(S_x\) | the truncated model geometries with nonnegative weights |
| diagonal | \(\operatorname{diag}(c)\), \(c \ge 0\) | \(S_x\) | the model axes with unequal, possibly sparse weights |
| block | \(\operatorname{blockdiag}(C_i)\), \(C_i \succeq 0\), \(\operatorname{rank} C_i \le s_i\) | \(S_x\) | an independent learned metric per model |
| shared | \(C \succeq 0\), \(\operatorname{rank} C \le s\) | \(S_x\) | a learned rank-\(s\) readout of the joint model span |

### 4.1 Diagonal

\[
\tfrac12\|G_x - \sum_j c_j t_j t_j^\top\|_F^2
= \text{const} - g^\top c + \tfrac12 c^\top M c,
\qquad
g_j = t_j^\top G_x t_j = R_{\cdot j}^\top S_x R_{\cdot j},
\qquad
M_{jk} = (t_j^\top t_k)^2 .
\]

\(M\) is model-only and computed once; \(g\) is read from the compressed form.
The fit is a nonnegative quadratic program in \(\sum_i r_i^\ast\) variables.

Measured (oracle §O5): the objective identity holds to \(1.5\times10^{-11}\)
and \(g\) from \(S_x\) agrees with \(g\) from \(G_x\) to \(4.6\times10^{-13}\).

### 4.2 Block

With \(Q_i\) an orthonormal basis of \(\operatorname{col}(R_{\cdot,\text{block } i})\)
in \(Q\)-coordinates, the exact block update is
\(A_i \leftarrow [Q_i^\top (S_x - \sum_{j\ne i} Q_j A_j Q_j^\top) Q_i]_{+,s_i}\);
each exact block update is the minimizer over its block with the others
held fixed, so the objective never increases. Because the rank constraint is
nonconvex, a limit point is a blockwise minimizer only when every block
minimizer is unique, which eigenvalue ties break; the descent is therefore
reported with its sweep count and convergence flag rather than as a global
optimum. The trace-penalized relaxation soft-thresholds the eigenvalues and
is admitted as a named alternative.

Measured (oracle §O8, 2026-09-04, block ranks \((1, 2)\)): the objective fell
from \(120{,}449\) to \(1{,}261.9\) in 14 sweeps with a largest increase of
\(4.5\times10^{-13}\); one further exact update of every block moved the
iterate by \(3.6\times10^{-10}\), so the limit point is a blockwise
minimizer; and the nesting of the feasible sets shows in the objectives in
compressed coordinates, shared \(\le\) block \(\le\) diagonal
(\(1.1\times10^{-26} \le 1{,}261.9 \le 24{,}257\)).

### 4.3 Shared

The shared fit **is** §3.2 and is not implemented twice.

**Normative.** (a) The returned object for block and shared structures is
\(C\) (or its spectrum and eigenvectors), never a factor \(W\): \(W \mapsto WO\)
is unidentified, and under overlapping model spans \(C\) itself is identified
only modulo the null space of \(R\), which is why §1.5 refuses overlap. (b) Every fit records its structure and the sentence in the
last column above; the shared structure additionally records that it tests
the joint span and not the supplied geometries. (c) The isotropic structure
is the only one that tests the model geometries as supplied; the receipt
says so.

---

## 5. The trace decomposition

**Claim 5.** The proposal's model-addressable versus model-orthogonal split
is stated here in trace form, because a Frobenius norm of \(G_x\) is not a
fixed linear query and a trace is.

\[
\operatorname{tr}(H G_x H)
= \underbrace{\operatorname{tr}(P G_x)}_{\text{model-addressable}}
+ \underbrace{\operatorname{tr}\big((H - P) G_x\big)}_{\text{model-orthogonal}},
\qquad
\operatorname{tr}(P G_x) = \operatorname{tr}(S_x).
\]

Both terms are fixed symmetric bilinear queries with vanishing marginals;
both are signed and unbiased; the identity is exact at every measurement.
It mirrors `coherent + configuration = total` and is orthogonal to it: a
measurement is characterized along the spatial axis by that split and along
the representational axis by this one.

Measured (oracle §O6): \(371.79 = 358.10 + 13.69\), and
\(\operatorname{tr}(PG_x) - \operatorname{tr}(S_x) = 0\) to machine precision.
Package route at README searchlight 144, with the rank-revealing basis
(\(q = 1\) for the category model): \(4.10545 = 4.10297 + 0.00248\), and
\(\operatorname{tr}(S_x) = 4.10297\).

**Normative.** (a) The model-orthogonal term is read from the **original**
plan with the operator \(H - P\); it is not available from the lowered form
alone and the view says so. (b) The rank-\(s\) share of the addressable
term, \(\operatorname{tr}[S_x]_{+,s}\), lives on the latent layer of §3.2 and
is never reported as a fraction of a signed denominator. (c) The
terminology rows of §9 name the two components.

---

## 6. Saturation is a diagnostic and, at the limit, a refusal

**Theorem 3 of the proposal.** If \(\operatorname{col}(T) = \mathbf 1^\perp\),
every centered PSD geometry of rank at most \(s\) is representable as
\(TWW^\top T^\top\), and the shared fit reduces to the best rank-\(s\) PSD
approximation of \(H G_x H\) irrespective of which models produced \(T\).

Measured (oracle §O4): with a basis of all of \(\mathbf 1^\perp\) the fit
reproduces a random rank-2 target to \(1.2\times10^{-14}\); with the
three-dimensional model span the residual is \(6.16\). The shared fit's
degrees of freedom are \(qs - s(s-1)/2\) with \(s\) clamped to \(q\): \(5\) against \(9\) for the
unconstrained rank-2 form in that fixture.

**Normative.** (a) The basis records `span_fraction = q / (n - 1)` and every
fit records `df_fit = qs - s(s-1)/2` beside `df_free = (n-1)s - s(s-1)/2`.
(b) A shared fit with `span_fraction == 1` signals an
`effect_capability_refusal` with capability `"model_geometry_test"`,
namespace `"model_coordinate"`, reasons `"model_span_saturated"`, and the
remedies: reduce the model ranks; use the diagonal or block structure; or
evaluate the learned readout under the cross-fit of §7. (c) Above a declared
warning fraction the fit is admitted and labelled `near_saturated`. (d)
Isotropic and diagonal structures are never refused on this ground, because
they do not learn a rotation; their receipts still carry the fraction.

---

## 7. Cross-fitting the learned readout

**Claim 7.** The learned model metric \(C_x\) is the dual of the learned
neural metric \(K_x\), and the discipline that keeps a learned neural metric
honest applies to it unchanged.

### 7.1 The duality

\[
B_a\,K\,B_b^\top \quad\longleftrightarrow\quad T\,C\,T^\top .
\]

\(K\) declares which distinctions among neural features count; \(C\) declares
which distinctions among model coordinates count. When \(K\) is learned from
residuals, `metric_training_policy("exclude_evaluation")` trains each edge's
metric on partitions other than that edge's two endpoints
(`R/metric-learning.R:367-382`), errors at plan time, from `plan_geometry()` (`R/geometry-plan.R:681`),
when an edge is left with no training partitions (`:323-348`), and records
the policy in plan identity. When \(C\) is learned from the geometry, the same policy
applies with the same unit, the partition edge.

### 7.2 The edge-disjoint cross-fit

For an evaluation edge \((a, b)\), learn the rank-\(s\) eigenbasis
\(V_{(a,b)}\) of the compressed form built from edges disjoint from
\(\{a, b\}\); evaluate \(\operatorname{tr}\big(V_{(a,b)}^\top S_x^{(a,b)}
V_{(a,b)}\big)\); reduce over evaluation edges with the declared weights.
Because \(V_{(a,b)}\) is a function only of data independent of the products
it weights, the learned-metric clause of `effect-form-v1` §8 applies edge by
edge: the cross-fitted rank-\(s\) energy is unbiased for the estimand defined
by the realized readout, conditionally on that readout, and is zero in
expectation under pure noise. It is not unbiased for a fixed-readout
estimand, and §O7 shows only the null.

Measured (oracle §O7): under pure noise the cross-fitted rank-1 energy has
mean \(-0.0415\) (se \(0.0249\)), against the positive plug-in latent energy
of §3.2.

**Normative.** (a) The evaluation edges and the training edges are the
**declared pairing's** edges: for evaluation edge \((a, b)\) the readout is
learned on the pairing's edges disjoint from \(\{a, b\}\), with their
weights renormalized, never on a product the pairing did not declare. The
cross-fit therefore needs a pairing over at least four partitions in which
every edge has a disjoint edge, the same floor `heterogeneity()` enforces
(`R/population-heterogeneity.R:86-101`); a pairing that fails it (fewer
partitions, only self products, a star) makes the reader signal an
`effect_capability_refusal` with capability `"cross_fitted_model_energy"`,
namespace `"model_coordinate"`, reasons `"insufficient_disjoint_edges"`.
(b) The reading records, in `$cross_fit$edges`, which edges evaluated and
which partitions and how many pairing edges trained each readout, together
with the pairing's declared independence; the execution receipt's schema is
sealed and carries none of it. (c) The cross-fitted energy is signed and is
reported on the estimation layer; the plug-in rank-\(s\) energy is reported
on the latent layer; the two are never presented as the same number. (d)
This version defines the cross-fit for the shared structure only; the
isotropic, diagonal and block plug-in energies are fitted quantities too and
are not unbiased under the null, and their cross-fit is deferred (§14).

### 7.3 Condition hold-out is a different claim

No abstraction for generalizing over conditions or effects exists in the
package; `generalizes_over` names a partition axis only (`R/pairing.R:59-63`).
Holding out conditions tests whether a learned readout transfers to new
stimuli, which is a different scientific claim from reproducibility of the
readout in independent data for the same stimuli. The Nyström extension
\(T_B = K_{BA} U_r \Lambda_r^{-1/2}\) is the correct tool for it, with the
centering origin fixed at the training-condition mean.

**Normative.** Condition hold-out is deferred to its own object with its own
plan identity (§14, M10). This version neither implements nor emulates it.

---

## 8. GLS in compressed coordinates (admitted, deferred)

`sampling_covariance(x, queries = )` reduces a sampling covariance to a bank
of zero-sum contrast vectors, that is, rank-one operators
(`.sampling_query_bank()`, `R/evidence-sampling-product.R:679-717`; entry
`:1154-1180`). Because \(Q^\top\mathbf 1 = 0\), the bank
\(\{q_j\} \cup \{q_j + q_k\}_{j<k}\) is admissible, and every off-diagonal
entry follows by polarization,
\(q_j^\top G q_k = \tfrac12\big[(q_j+q_k)^\top G (q_j+q_k) - q_j^\top G q_j -
q_k^\top G q_k\big]\). The sampling covariance of
\(\operatorname{svec}(S_x)\) is therefore a fixed linear image of that bank's
covariance, at dimension \(q(q+1)/2\), and the proposal's §8 generalized
least squares fit can be posed there rather than over \(n(n-1)/2\) distances.

**Normative.** This route is admitted under `evidence-sampling-v1` and is
deferred (§14, M11). Any such fit is a latent reading unless cross-fitted
under §7.

---

## 9. Vocabulary

Rows to be added to `design/terminology.md` (M8). Nouns and adjectives
follow the existing convention.

| Preferred prose | API token | Meaning | Deprecated shorthand | Forbidden interpretation |
|---|---|---|---|---|
| model basis | model basis | Declared, centered, rank-revealing basis of a model family's joint span; data-independent | model kernel | a neural subspace |
| model-coordinate geometry | lowered form | The complete cross-generalized form read in model coordinates, \(Q^\top G_x Q\) | projected RDM | a denoised geometry |
| model-addressable component | addressable | Signed energy of the geometry within the model span, \(\operatorname{tr}(PG_x)\) | explained energy | variance explained by the model |
| model-orthogonal component | orthogonal | Signed energy of the centered geometry outside the model span | residual energy | noise |
| latent rank-\(s\) model energy | latent_rank | PSD rank-\(s\) projection of the addressable component; latent layer | low-rank fit | an unbiased estimate |
| cross-fitted model energy | cross_fit | Rank-\(s\) energy with the readout learned on disjoint edges; estimation layer | held-out fit | condition generalization |
| span fraction | span_fraction | \(q/(n-1)\); one at saturation | model rank | goodness of fit |

---

## 10. Layering and naming

Per `design/architecture.md`, every new file is registered in
`tests/testthat/test-architecture.R` and may call downward or sideways only.

| piece | layer | sits beside | proposed token |
|---|---|---|---|
| model basis value: construction, admission, centering, saturation, signature; produces an extractor and an effect map | 2 | `effect-map.R`, `metric.R` | `model_basis()` |
| lowering entry and the gate of §2.2 | 2 / 3 | `relation.R`, `evidence-task.R` | lowering by `relation(..., extract = basis)` or `effect_map`; gate in `.new_evidence_task()`, dormant until a normalizer is public |
| rank budget on the latent projection | 5 | `latent.R` | `latent_geometry(x, rank = s)` |
| model-geometry view: trace split, isotropic / diagonal / block fits, shared delegating to the latent layer, print / format / as.data.frame | 5 | `latent.R`, `views.R` | `model_geometry()` |
| cross-fitted rank-\(s\) energy | 5 | `population-heterogeneity.R` (precedent, `heterogeneity()`) | `model_geometry(..., training = )` |

**Normative.** (a) `rsa()` is not modified. (b) The word *geometry* keeps
its meaning, the neural cross-generalized form; the model-side value is a
*basis*, not a geometry. (c) The view's structure fits are readers of a
materialized lowered form; nothing in layer 5 executes a plan except by
re-planning through the layer-3 entries, as `population-heterogeneity.R`
does with `plan_geometry()` and `plan_population()`.

---

## 11. What exists today, and what this contract requires

| requirement | status | evidence |
|---|---|---|
| a declared linear map on the effect axis applied at read time | **exists** | `effect_extractor()` (`R/extractor.R:34`); identity extractor at `R/relation.R:299`; `effect_map()` (`R/effect-map.R:158`) |
| complete form, spectrum, latent projection on a lowered relation, fixed metric | **exists** | `design/oracles/model-coordinate-package-route.R` §R3: `geometry_spectrum()` and `latent_geometry()` ran unchanged on the lowered form |
| complete lowered form under a learned local metric | **not reachable** | `R/compiler.R:219-224`, `:274-290`; deferred (§14) |
| population pooling of lowered forms | **exists by construction** | signature equality check `R/population-plan.R:232-251`; requires §1.4 |
| sampling covariance of compressed coordinates | **exists by polarization** | zero-sum contrast bank, `R/evidence-sampling-product.R:679-717`; §8 |
| partition-disjoint training discipline | **exists** | `metric_training_policy()` and its records, `R/metric-learning.R:259-382` |
| rank-revealing model basis with admission, centering, saturation, overlap refusal | **exists** (2026-09-04) | `model_basis()`, `R/model-basis.R`; `tests/testthat/test-model-basis.R` |
| compiler gate for nonlinear edge transforms | **exists**, dormant (2026-09-04) | `.model_coordinate_lowering_gate()` in `.new_evidence_task()`; `tests/testthat/test-model-coordinate-lowering.R` |
| rank budget and two-part moved mass on the latent layer | **exists** (2026-09-04) | `latent_geometry(rank =)`, `.latent_rank_psd_form()`; `tests/testthat/test-latent-geometry.R` |
| trace split and structure fits | **exists** (2026-09-04) | `model_geometry()`, `R/model-geometry.R`; `tests/testthat/test-model-geometry.R` |
| edge-disjoint cross-fit of the readout | **exists** (2026-09-04) | `model_geometry(training =)`; `tests/testthat/test-model-geometry.R` |
| contract test against the oracle | **exists** (2026-09-04) | `tests/testthat/test-model-coordinate-contract.R` |
| vocabulary, novelty ledger, prior-work entry | **exists** (2026-09-04) | `design/terminology.md`, `vignettes/novelty.Rmd`, `design/relation-to-prior-work.md`, README |

---

## 12. Numerical contract

- Rank tolerances are relative to the largest absolute eigenvalue or
  singular value and are recorded.
- The lowered form agrees with the projected full form to \(10^{-12}\)
  relative on every component for linear edge transforms (§2.2).
- The Frobenius split of §3.2 holds to \(10^{-10}\) relative.
- The trace identity of §5 holds to \(10^{-10}\) relative.
- Under pure noise, the mean of \(\operatorname{tr} S_x\) and of the
  cross-fitted energy lie within four Monte Carlo standard errors of zero;
  the plug-in latent energy lies above four standard errors of zero.

---

## 13. Test and oracle index

| claim | statement | evidence |
|---|---|---|
| 1 | rank-revealing factor; the thin-QR trap leaks an uncentered direction | oracle §O2 (trap: \(|Q^\top\mathbf 1| = 1.10\), leakage \(8.12\)); package route on the README fixture; **test M7** |
| 1.3 | baseline invariance iff \(Q^\top\mathbf 1 = 0\) | oracle §O2 (\(17.2\) raw versus \(1.1\times10^{-13}\) compressed); **test M7** |
| 1.5 | overlapping spans leave \(C\) unidentified; refusal | **test M2** (refusal `model_span_overlap`); unmeasured by the oracle, whose fixture has \(q = \sum_i r_i^\ast = 3\) |
| 2.2 | lowering commutes for total, coherent, configuration | package route §R1 \(8.9\times10^{-16}\) / \(8.9\times10^{-16}\) / \(7.2\times10^{-16}\); oracle §O1; **test M7** |
| 2.2 | lowering does not commute across nonlinear edge transforms | **test M3** (correlation-normalized plan refused with `nonlinear_edge_transform`) |
| 2.3 | the bank reproduces one compressed entry | package route §R2 \(8.9\times10^{-16}\) |
| 3.1 | \(S_x\) is unbiased | oracle §O7 (\(0.0049\), se \(0.0281\); signal bias \(-0.877\), se \(0.445\)) |
| 3.2 | Theorem 4 and the Pythagorean split | oracle §O3 (split error \(7.3\times10^{-12}\); \(15{,}219.4\) versus \(57{,}504.6\) random and \(15{,}220.1\) local); **test M4** |
| 3.2 | the plug-in latent energy is biased | oracle §O7 (\(0.6446\), se \(0.0167\)) |
| 4.1 | diagonal objective identity; \(g\) from the compressed form | oracle §O5 (\(1.5\times10^{-11}\); \(4.6\times10^{-13}\)); **test M5** |
| 4.2 | block descent is non-increasing and reaches a blockwise minimizer; shared \(\le\) block \(\le\) diagonal | oracle §O8 (max increase \(4.5\times10^{-13}\); fixed-point move \(3.6\times10^{-10}\); \(1.1\times10^{-26} \le 1{,}261.9 \le 24{,}257\)); **test M5** |
| 5 | addressable + orthogonal = centered total; addressable = \(\operatorname{tr} S_x\) | oracle §O6; package route §R4 at searchlight 144 (\(4.10545 = 4.10297 + 0.00248\)); **test M5** |
| 6 | saturation; degrees of freedom | oracle §O4 (\(1.2\times10^{-14}\) versus \(6.16\); df \(5\) versus \(9\)); **test M2** (refusal `model_span_saturated`) |
| 7.2 | the edge-disjoint cross-fit is unbiased under noise | oracle §O7 (\(-0.0415\), se \(0.0249\)); **test M6** |
| 7.2 | four-partition floor | **test M6** (refusal `insufficient_disjoint_edges`) |

Two scripts back this index. The oracle is plain R and never loads the
package; the package-route witness loads it and must run from the
repository root:

```sh
Rscript design/oracles/model-coordinate-geometry.R
Rscript design/oracles/model-coordinate-package-route.R
```

Both ran clean on 2026-09-04. Every oracle number quoted here reproduces
from the first under seed `20260904`; every package-route number (§1.1's
\(1.15\), the §2.2 table, §2.3, §5, and the §R rows above) reproduces from the
second on the README fixture. As with the other oracles, nothing executes
either automatically; M7 promotes their laws into
`tests/testthat/test-model-coordinate-contract.R` rather than shelling out
to the scripts.

---

## 14. Scope boundary

Admitted by this version: fixed model RDMs or fixed model features; the four
structures of §4; the trace split of §5; the plug-in latent fit of §3.2 and
the edge-disjoint cross-fit of §7.2 for the shared structure; per-subject
and population-level lowered forms.

Deferred, each to its own ticket and plan identity:

- **M10** condition hold-out with the Nyström extension (§7.3);
- **M11** GLS in compressed coordinates under `evidence-sampling-v1` (§8);
- **M12** feature-learned model factors \(T_i = F_i L_i\), which test a
  learned readout of the features rather than a fixed geometry and need the
  cross-fit of §7 before they are admitted at all;
- complete lowered forms under a learned local metric (§2.3), which the
  compiler cannot produce today;
- lowered-then-correlated readings as a named estimand (§2.2);
- bases with overlapping model spans, reported through \(A\) only (§1.5);
- the edge-disjoint cross-fit of the isotropic, diagonal and block fits
  (§7.2(d)).

Not claimed: that this reading unifies every representational model; that
the shared structure tests the supplied RDMs; that any latent quantity in
this contract is an estimate; or that a response envelope in the sense of
envelope regression is identifiable from the form (the proposal's Theorem 5
shows it is not).

---

## 15. Fresh-context review (2026-09-04)

A reviewer with no prior context ran both scripts, checked every quoted
number, every file:line citation, the mathematics, and the consumed
contracts. Verdict: accept with amendments. Every oracle number matched; the
package-route numbers, the file:line citations, and the mathematics were
correct except as corrected below.

### 15.1 Corrections applied

1. **§5, §13.** The first draft's package-route trace split
   (\(\operatorname{tr} S_x = 4.09647\), orthogonal \(0.00898\)) was measured
   on the rank-deficient thin-QR basis that §1.1 forbids; on it \(H - P\) is
   not a projector. Re-measured with the rank-revealing basis
   (`design/oracles/model-coordinate-package-route.R`, \(q = 1\)):
   \(4.10545 = 4.10297 + 0.00248\). The linear identities of §2.2 and §2.3
   held on either basis and were re-measured on the compliant one.
2. **§2.1(c), §2.3, §11, §14.** Under a learned local metric the compiler
   refuses full materialization and contracts one signed contrast per call,
   so the lowered form is not an escape from the bank's one-column limit;
   both routes are limited alike, and complete lowered forms under a learned
   metric are deferred.
3. **§1.5 (new), Notation, §4.** \(R\) is invertible only when \(T\) has full
   column rank; overlapping model spans leave \(C\) identified only modulo the
   null space of \(R\). This version refuses overlap with
   `model_span_overlap` and records `span_overlap_rank` and
   `basis_condition`.
4. **§3.2.** `R/latent.R` states three numbered decisions at `:22-44` and the
   reading line at `:63`; the citation was corrected.
5. **§7.2, §10.** The exported reader is `heterogeneity()`; there is no
   `population_heterogeneity()`.
6. **§8, §11.** The sampling-covariance bank admits zero-sum contrast vectors
   only; \(\operatorname{svec}(S_x)\) follows by polarization because
   \(Q^\top\mathbf 1 = 0\).
7. **§3.1, §7.2.** Unbiasedness of the compressed form rests on the Corollary
   of `effect-form-v1` §8; the cross-fitted energy is unbiased conditionally
   on the realized readout, per that contract's learned-metric clause, and the
   oracle shows only the null.
8. **§2.2, §10.** Edge normalizers live on the evidence stage plan (layer 3),
   are unexported, and are never set by a public entry; the gate is dormant,
   and the refusal now says why lowered-then-correlated is refused rather
   than named.
9. **§6.** With a full-span basis the shared fit approximates \(H G_x H\),
   not \(G_x\); the degrees-of-freedom count clamps \(s\) to \(q\).
10. **§4.2, §13.** The block-descent convergence sentence is qualified and
    marked unmeasured.
11. **§13.** The package-route script is now committed beside the oracle;
    the closing paragraph says which script reproduces which numbers.
12. Refusal fields are `reasons` and `remedies` (plural); `.preflight_metric_training()`
    raises an input error at plan time rather than a refusal; `R/compiler.R:469-476`
    is a contract error and the user-facing gates are the `scheduled_metric_*`
    refusals; `population-heterogeneity.R` re-plans through layer-3 entries;
    `.validate_rdm_models()` begins at `R/views.R:385`.

### 15.1a Implementation review amendments (2026-09-04, M2 and M4)

Fresh-context reviews of the M2 and M4 implementations amended the
contract in three places, each recorded above in its section: §1.1(d), the
centered-frame construction and the two-sided tolerance on the model
spectrum; §1.3, the class of the centering failure; §1.4, the effective
ranks in the signature. Two package facts the reviews established are
recorded here so that they are not rediscovered. The `$extractor` product of
a basis lowers a labelled beta block by row name and refuses a mismatch,
but an unlabelled block by position; and every latent identity changed when
the rank budget entered the latent schema (`latent-sha256` schema 2), which
nothing in the repository pinned except one frozen exemplar artifact
(`exemplars/haxby2001/results/conservative-geometry.rds`, regenerated by
M9).

### 15.2 Gaps carried into tickets

| gap | statement | ticket |
|---|---|---|
| G1 | overlap identification: admission rule, `span_overlap_rank`, `basis_condition`, refusal test | M2 |
| G2 | the §2.2 gate is dormant; test against the internal stage plan until a normalizer is public | M3 |
| G3 | block-descent convergence is unmeasured | M5 — closed 2026-09-04 by oracle §O8 (§4.2, §13) |
| G4 | polarization bank for the compressed sampling covariance | M11 |
| G5 | complete lowered forms under a learned local metric | deferred, §14 |
