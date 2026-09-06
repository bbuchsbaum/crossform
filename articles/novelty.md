# What is novel in crossform?

`crossform` keeps condition effects, spatial weights and partition
pairings in one declared geometry. Contrast energies, distances and
fixed RSA coefficients read different aspects of that geometry.
Retaining its structure also makes it possible to separate local-mean
from spatial-pattern energy and to describe exactly what must generalize
across runs or sessions.

This guide explains what that organization adds to established methods
and which claims have supporting evidence. For the analysis itself,
start with the
[introduction](https://bbuchsbaum.github.io/crossform/articles/introduction.md);
for the underlying equalities, read [the common-geometry
derivation](https://bbuchsbaum.github.io/crossform/articles/common-geometry-equivalence.md).

## The category difference

The calculus separates two choices that are usually bundled together:
whether the experimental spaces are the same, and whether the neural
measurements are the same.

Each cell carries its evidence status from the ledger below, so the
table claims territory only at the strength actually earned.

|  | Same neural measurement | Different neural measurements |
|----|----|----|
| **Same experimental space** | retained signed marginals; crossvalidated contrast energy — a cvMANOVA-class statistic under a fixed metric ([Allefeld and Haynes, 2014](https://doi.org/10.1016/j.neuroimage.2013.11.043)); ordinary representational geometry *(demonstrated)* | cross-node effect coupling *(implemented; small-node contract)*; normalized connectivity only when its sampling contract is met |
| **Different experimental spaces** | rectangular ER-RSA or cross-task geometry *(demonstrated on a designed simulation)* | rectangular cross-domain, cross-region forms *(prospective)* |

Signed activation is a retained first-moment marginal, not a bilinear
statistic. Its reproducible energy, squared-distance geometry, and the
other entries in the table are second-order closures. Keeping both
channels in one plan lets a single fitted relation report the signed
effect together with the geometry that asks whether it reproduces.

The irreducible second-order observable is

``` math
\mathscr E_{\Gamma}(H,K)
=
\sum_{r,s}\Gamma_{rs}
\operatorname{tr}\!\left(
H^\top B_{L,r}K B_{R,s}^\top
\right).
```

Here the relations $`B_{L,r}`$ and $`B_{R,s}`$ bind named experimental
and neural spaces, $`H`$ asks an experimental-side question, $`K`$
specifies a neural-side metric or bridge, and $`\Gamma`$ declares which
independently estimated relations must reproduce one another. Closing
the neural boundaries gives an effect form

``` math
F_K=\sum_{r,s}\Gamma_{rs}B_{L,r}K B_{R,s}^\top,
```

while closing the experimental boundaries gives the adjoint neural form

``` math
Q_H=\sum_{r,s}\Gamma_{rs}B_{L,r}^\top H B_{R,s}.
```

The scalar pairing can be read from either side:

``` math
\mathscr E_{\Gamma}(H,K)
=\langle H,F_K\rangle_F
=\langle Q_H,K\rangle_F.
```

This is more than saying that one kernel computes several distances. It
gives one typed construction for self- and cross-forms, square and
rectangular experimental axes, and local and cross-location neural
measurements.

## Relation to cvMANOVA

The closest existing statistic to crossvalidated contrast energy is the
cross-validated MANOVA of [Allefeld and Haynes
(2014)](https://doi.org/10.1016/j.neuroimage.2013.11.043). Their pattern
distinctness is a leave-one-run-out trace that combines a contrast
estimated from the held-out run with one estimated from the training
runs, under an error-covariance metric: the same *shape* of object as
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
and it should be read as the direct ancestor rather than as a neighbor.
The differences are specific. cvMANOVA normalizes the contrast by its
estimated error term in a Bartlett–Lawley–Hotelling trace, corrects the
resulting bias with an explicit multiplicative factor, and recommends
dividing the map by $`\sqrt p`$ so that the null variance does not grow
with searchlight size. `crossform` does none of the three: the metric is
declared rather than estimated inside the statistic, unbiasedness comes
from the pairing itself because every product multiplies estimates from
two different partitions, and frame normalization is a statement about
what a measurement is rather than a null-variance correction. cvMANOVA
also supplies distributional theory and permutation inference, which
this package does not have at all.

Allefeld and Haynes’ Figure 6 already runs the voxel-axis control:
cvMANOVA on searchlight-mean-only data beside cvMANOVA on mean-removed
data. Those are two re-analyses of two modified datasets, and their
results do not sum to the unmodified analysis. That gap is where this
package’s own voxel-axis claim lives; it is stated under *What is
distinctive* below. The full attribution, including crossnobis,
`rsatoolbox`, pattern-component modeling, and Framed RSA, is in
[`design/relation-to-prior-work.md`](https://github.com/bbuchsbaum/crossform/blob/main/design/relation-to-prior-work.md);
the noise-unbiasedness theorem and its failure modes are §8 of
[`design/effect-form-contract.md`](https://github.com/bbuchsbaum/crossform/blob/main/design/effect-form-contract.md).

## Relationship to RSA and `rsatoolbox`

RSA, crossnobis distance, noise-precision estimation, searchlight RSA,
model fitting, noise ceilings, and inference over subjects and
conditions are established work. Searchlight information mapping
([Kriegeskorte et al., 2006](https://doi.org/10.1073/pnas.0600244103)),
the RSA toolbox ([Nili et al.,
2014](https://doi.org/10.1371/journal.pcbi.1003553)), the reliability
advantage of crossvalidated Mahalanobis distances over correlation
distance and classification accuracy ([Walther et al.,
2016](https://doi.org/10.1016/j.neuroimage.2015.12.012)), and whitened
unbiased RDM similarity ([Diedrichsen et al.,
2021](https://doi.org/10.51628/001c.27664)) are prior art that this
package neither reimplements nor claims.
[`rsatoolbox`](https://rsatoolbox.readthedocs.io/) is a sophisticated
reference implementation of an RDM-centered workflow: it supports
multiple dissimilarities including crossnobis, residual-based noise
precision, searchlights, fixed and flexible models, and inferential
procedures for different generalization targets. The inference it
implements is developed in “Statistical inference on representational
geometries” ([Schütt et al.,
2023](https://doi.org/10.7554/eLife.82566)), which treats generalization
over subjects and conditions explicitly. `crossform` does not claim any
of those ingredients as inventions.

Nor is a shared interface across ROIs, searchlights, whole-brain maps,
and multiple multivariate measures new by itself. The Decoding Toolbox
([Hebart et al., 2015](https://doi.org/10.3389/fninf.2014.00088)),
CoSMoMVPA ([Oosterhof et al.,
2016](https://doi.org/10.3389/fninf.2016.00027)), and PyMVPA ([Hanke et
al., 2009](https://doi.org/10.3389/neuro.11.003.2009)) established
powerful unifying software abstractions. The claim here is the exact
algebraic relationship among the outputs, not the fact that one package
can dispatch several analyses.

The common fixed-linear subset is simple. If a square self-form is
$`G=BKB^\top`$, squared-distance extraction is a linear map
$`d=\mathcal D(G)`$. A fixed linear RSA readout $`\beta=Cd`$ can
therefore be compiled into a query $`H_\beta`$ such that

``` math
\beta=\langle H_\beta,G\rangle_F.
```

This observation licenses a controlled parity comparison; it is not by
itself a scientific novelty result. That comparison has now been run
against a version-pinned Python `rsatoolbox`, in addition to the rMVPA
comparison in the Haxby exemplar: see
[`exemplars/rsatoolbox-parity`](https://github.com/bbuchsbaum/crossform/tree/main/exemplars/rsatoolbox-parity).

The intended distinction is that an RDM is an optional view of an effect
form, not the mandatory intermediate representation. The principal
`rsatoolbox` `RDMs` abstraction stores dissimilarities over one shared
pattern axis as a vector or square matrix. `crossform` instead keeps the
identified relation, left and right axes, generalization operator,
measurement frame, and error channel, when the relation was fitted with
one, available until the requested query has been compiled. This is an
architectural comparison, not a claim that every nonlinear dissimilarity
or inferential method in `rsatoolbox` is contained in the bilinear core.

In particular, Pearson correlation distance is outside that core because
its per-pattern norm is nonlinear in the fitted patterns. That is a
statement about the bilinear core, not about what can be reproduced: the
parity exemplar does recover `rsatoolbox`’s within-sample correlation
distance to `4.9e-15`, reading `G_ij = (G_ii + G_jj - d_ij)/2` from
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) and
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
downstream of a guaranteed-PSD self form with strictly positive
diagonals – the case the correlation-distance policy already licenses.
What stays refused is `rdm(normalize=)` over signed cross-generalized
diagonals, which may be zero or negative. The package gives that
boundary a named policy rather than a quiet escape hatch; see [the
correlation-distance
policy](https://bbuchsbaum.github.io/crossform/articles/correlation-distance-policy.md)
([`vignette("correlation-distance-policy", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/correlation-distance-policy.md)
offline).

## What is distinctive

### 1. Boundary-closure unification

Crossvalidated contrast energy, squared-distance RDMs, fixed linear RSA,
ordered cross-domain similarity, and neural-side effect coupling are
queries or partial materializations of the same evidence pairing. This
is implemented and covered by independent small-matrix law tests.
Normalized connectivity has additional repeated-variation and covariance
requirements; a numerical off-diagonal block does not acquire that
interpretation from its shape alone.

The neural-side construction has close scientific neighbors in
[representational
connectivity](https://pmc.ncbi.nlm.nih.gov/articles/PMC2605405/),
[informational
connectivity](https://pmc.ncbi.nlm.nih.gov/articles/PMC3566529/), and
multidimensional-connectivity methods ([Basti et al.,
2020](https://doi.org/10.1016/j.neuroimage.2020.117179)). The
distinctive claim is not ownership of connectivity; it is that neural
coupling is the adjoint-side materialization of the same identified
experimental-neural pairing, with stronger normalized views admitted
only by their own contracts.

### 2. An exact voxel-axis partition that every fixed linear query inherits

For a frame row $`w`$ with positive mass $`a=\mathbf 1^\top w`$, the
admitted spatial metric splits as

``` math
D(w)=\frac{ww^\top}{a}+\left[D(w)-\frac{ww^\top}{a}\right],
```

with both terms positive semidefinite. Pulling the split through the
relation gives

``` math
G^{\mathrm{total}}
=G^{\mathrm{coherent}}+G^{\mathrm{configuration}},
```

and because each part is itself a *fixed* metric, every fixed linear
query $`\langle H,\cdot\rangle`$ of the form inherits the same additive
partition from one execution: the contrast energy, each RDM edge, and a
fixed linear RSA coefficient all decompose the same way, and each
component is separately noise-unbiased under the pairing (§8 of the
[effect-form
contract](https://github.com/bbuchsbaum/crossform/blob/main/design/effect-form-contract.md)).

That is the claim, and it is small enough to state at theorem strength:
*the two pieces are the images of the frame metric under complementary
$`D(w)`$-orthogonal projectors — $`P=\mathbf 1w^\top/a`$ satisfies
$`P^2=P`$ and $`D(w)P=ww^\top/a`$, which is symmetric — so the split is
exact, both parts are PSD, and every fixed linear query of the form
inherits it additively.* What is not claimed is the distinction between
regional-mean and pattern structure, which has a substantial literature:
[Davis et al., 2014](https://doi.org/10.1016/j.neuroimage.2014.04.037),
pattern-component modeling ([Diedrichsen et al.,
2018](https://doi.org/10.1016/j.neuroimage.2017.08.051)), and
second-moment accounts of RSA and encoding models ([Diedrichsen and
Kriegeskorte, 2017](https://doi.org/10.1371/journal.pcbi.1005508)) — nor
the elementary matrix identity itself.

One nearby literature has to be separated rather than absorbed, because
it concerns a *different axis*. [Garrido et
al. (2013)](https://doi.org/10.3389/fnins.2013.00174) argue against
subtracting the cocktail mean — each voxel’s mean across conditions — on
the grounds that it is contaminated by condition-specific signal and
introduces dependencies between conditions. That mean runs along the
condition axis; the partition above runs along the voxel axis, and they
are different operations. `crossform` subtracts neither. A zero-sum
contrast cancels any pattern shared additively by all conditions
exactly, without estimating it, and the voxel-axis common mode is
retained as `coherent` rather than removed.

The nearest published construction is **Framed RSA** ([Taylor and
Kriegeskorte, 2025](https://doi.org/10.1101/2025.07.10.664257)), which
restores the population-mean information ordinary RSA discards by
augmenting the pattern set with two reference patterns, the origin and a
uniform constant pattern. It is a genuine treatment of the same problem,
and it differs from the partition here in kind: augmentation changes the
analyzed set and the model comparison, whereas the split above changes
nothing, sums exactly, is computed once, and travels into every
downstream fixed linear query. [Allefeld and Haynes’
(2014)](https://doi.org/10.1016/j.neuroimage.2013.11.043) Figure 6
reaches the same scientific question by re-analyzing mean-only and
mean-removed data, two results that do not sum to the original. A search
of this literature did not find the exact partition, with PSD components
and inheritance by every fixed linear query, published in this form.
That is recorded as *apparently unpublished* — the strength a literature
search can actually support — rather than as a priority claim.

The one-plan family — signed contrast, the three energies with exact
recomposition, the RDM, the RSA coefficient, and the admitted analytic
standard error — is demonstrated on a generated planted-truth fixture in
the [introduction
vignette](https://bbuchsbaum.github.io/crossform/articles/introduction.md)
([`vignette("introduction", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/introduction.md)
offline). That is Gate 2 in the ledger below.

It now also runs on real data. On Haxby 2001 subject 1 — 12 runs, 8
categories, 577 ventral-temporal searchlights at 11.25 mm — the face
minus house contrast recomposes to `5.6e-17`, its total is positive at
576 of 577 searchlights, and at the peak searchlight the three numbers
are coherent 1.139, configuration 0.090, total 1.229. Over the 536
searchlights where the components form a nonnegative partition, the
coherent share has median 0.53 and interquartile range \[0.28, 0.77\]:
in this subject’s VT, roughly half of the reproducible face/house energy
is carried by the searchlight’s own common spatial mode. Read at one
whole-VT region from the same plan, coherent is 0.212 and configuration
0.224. The retained signed marginal is what makes those energies
interpretable — it is negative at 568 of 577 searchlights, so the common
mode runs *house above face*, a direction no energy can report. This is
one subject, condition means rather than GLM betas, and no inference; it
is a decomposition narrative, not a result about faces. See
[`exemplars/haxby2001`](https://github.com/bbuchsbaum/crossform/tree/main/exemplars/haxby2001).

### 3. Rectangular, identified cross-domain forms

The left and right experimental spaces may be different, unequal, and
directed. An encoding-by-retrieval form, for example, need not be forced
into a square condition axis. `plan_geometry(x, at, over, right = )`
compiles the rectangular plan publicly; axis-bound
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)s
execute against it query-first, and
[`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md)
materializes a rectangular form that satisfies the exact algebraic
identity `total = coherent + configuration`. The public test verifies
that identity numerically to `1e-12`. The engine, independent oracles,
and public constructor exist, and one analysis with unequal axes,
missing matches, match/control coupling, and a pair-space covariate has
now run end to end in
[`exemplars/er-rsa`](https://github.com/bbuchsbaum/crossform/tree/main/exemplars/er-rsa),
recovering three planted regional structures against a closed-form
ground truth. That analysis is a designed simulation, so what is
demonstrated is that the machinery estimates what it claims to estimate
on a design whose answer is known; the empirical claim awaits a public
encoding-retrieval dataset with per-trial betas.

### 4. Spatial frames with a qualified conservation law

Voxels, regions, searchlights, and whole-brain summaries enter as
measurement frames rather than separate analysis engines. Under a
conservative, feature-additive frame,

``` math
\sum_x L_x^\top L_x=M_\Omega
\quad\Longrightarrow\quad
\sum_x B L_x^\top L_x B^\top=B M_\Omega B^\top.
```

The [conservative-frames
guide](https://bbuchsbaum.github.io/crossform/articles/conservative-frames.md)
demonstrates this accounting with overlapping searchlights and an
unnormalized whole-domain comparator. The interpretation depends on a
conservative, feature-additive frame; arbitrary dense or learned metrics
do not inherit that guarantee.

### 5. Query-first execution

Selected contrasts, distance edges, and fixed linear RSA coefficients
compile without requiring the complete RDM as the public intermediate:
every RDM edge is the rank-one operator $`(e_i-e_j)(e_i-e_j)^\top`$, and
the kernels evaluate it as two row differences and a Hadamard product
instead of materializing a dense packed query. The recorded
large-condition benchmark compares selected pairs and complete RDMs
under matched estimands. It supports the practical benefit of computing
only requested quantities, with numerical and memory checks. Its timings
compare execution routes within crossform, not speed against another
package; R-heap measurements also differ from process resident memory.
The exact workloads and receipts are in the [certification
report](https://github.com/bbuchsbaum/crossform/blob/main/design/certification-report.md).

### 6. Generalization bound to estimand identity

Existing RSA inference already recognizes that generalization over
measurements, subjects, and conditions changes the scientific claim. The
distinctive formalization here is narrower: $`\Gamma`$ participates in
the identity of **every** effect-form estimand, not only the final
inferential procedure. The axis is typed, not inferred from labels:
`cross_partitions(relation, independence = "independent", generalizes_over = "run")`
and the same call with `"session"` produce distinct plan identities even
under identical generic partition names and identical fold counts, and
identity tests enforce it. Runs, sessions, tasks, item sets, sites, and
ordered cross-domains can all be represented by named pairing relations.

### 7. Learned forms with independent fixed-query evidence

[`model_basis()`](https://bbuchsbaum.github.io/crossform/reference/model_basis.md)
retains model directions and their strengths,
[`fit_geometry()`](https://bbuchsbaum.github.io/crossform/reference/fit_geometry.md)
learns a regularized low-rank PSD form on training observations, and
[`score_geometry()`](https://bbuchsbaum.github.io/crossform/reference/score_geometry.md)
evaluates the frozen result on independent observations. The fitting
problem is nonlinear. Its evaluation is an ordinary fixed bilinear
readout plus a known prediction cost:

``` math
\Delta=2\langle F,G_{\mathrm{test}}\rangle_F-\|F\|_F^2.
```

The model’s positive eigenvalues enter an inverse-kernel trace penalty,
so an anisotropic full-span model still imposes directional preferences
when the penalty is positive. At zero penalty only the admitted span
matters; a full-span fit is labelled a generic baseline. Compression
transforms each partition’s effects before geometry accumulation, and a
small spectral solve learns the form without a complete neural RDM.

Low-rank factorization, PSD projection and kernel regularization are not
claimed as new ingredients. The package contribution is their
integration with typed effect coordinates, observation-origin admission
and a frozen predictive readout. Conditional on training, the expected
gain equals improvement in squared geometry error when the test form is
unbiased for the target signal. It is not a trace energy: under zero
signal its expectation is $`-\|F\|_F^2`$, even though the signed test
inner product is centered on zero.

The [predictive geometry
guide](https://bbuchsbaum.github.io/crossform/articles/predictive-geometry.md)
([`vignette("predictive-geometry")`](https://bbuchsbaum.github.io/crossform/articles/predictive-geometry.md)
offline) demonstrates the public fixed-split workflow. Its evidence
concerns independent runs on the same conditions, with identity or fixed
SPD neural metrics and Frobenius geometry loss. New-condition transfer,
generic GLS, learned neural metrics and generic inference for selected
eigenvalues are outside this public contract. The local statistical and
numerical validation is recorded separately in the [predictive
certification
report](https://github.com/bbuchsbaum/crossform/blob/main/design/predictive-geometry-certification.md);
it is not an empirical demonstration or a claim of universal speed
advantage.

## The contract is the proof mechanism

The calculus defines the scientific objects. The contract is how those
objects survive real execution:

- plan identity records the relation, queries, metric, frame, units, and
  generalization relation;
- receipt identity records storage, tiling, execution path, and
  numerical diagnostics without redefining the plan;
- capabilities such as symmetry, self-form status, positive
  semidefiniteness, fixed-metric status, and retained uncertainty are
  construction guarantees;
- optimized paths are required to carry independent numerical oracles;
  and
- refusals are first-class results when an interpretation has not been
  earned.

The retained error channel is an important consequence. For the admitted
fixed-metric equal-partition model, an RSA coefficient is a fixed linear
functional of the RDM, so its covariance transports exactly from the
analytic RDM covariance of [Diedrichsen, Provost, and Zareamoghaddam
(2016)](https://doi.org/10.48550/arXiv.1607.01371). Precomputed effects
without an identified error channel are not reverse-engineered from the
spread of their edges. Refusals are classed conditions:
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)
returns the missing capability, every unmet reason, and remedies as
data, and
[`sampling_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/sampling_capabilities.md)
answers the admission question before it is provoked. The executable
[failure
gallery](https://bbuchsbaum.github.io/crossform/articles/failure-gallery.md)
([`vignette("failure-gallery", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/failure-gallery.md)
offline) shows six realistic errors that the package guards against.
Callable unsupported interpretations return classed refusals. Raw
crossvalidated estimates remain signed; explicitly named PSD projections
produce separate descriptive objects. Changes in generalization produce
distinct estimand identities.

## Evidence ledger

The canonical ledger uses eight evidence classes because proofs,
internal oracles, external parity, simulations, retrospective
illustrations, prospective protocols, completed real-data results, and
independent replications answer different questions. The full ledger
below is generated from the package certification artifact; expand it
when you need the status of an individual governed claim. The
definitions and promotion rules are in
[`design/evidence-status-ledger.md`](https://github.com/bbuchsbaum/crossform/blob/main/design/evidence-status-ledger.md).
The companion claim registry assigns exactly one owner and current
status to each claim, and the promotion history records why that status
changed. The strongest new population evidence is matched simulation:
the interval court covers the declared Gaussian, heteroskedastic,
heavy-tailed, and informative coverage regimes, while the hierarchical
interpretability court recovers its planted component ordering and scale
profiles. Neither is empirical evidence. The discovery and replication
protocols are frozen and rehearsed, but current readiness is `BLOCKED`;
no completed real-data result or independent replication exists.

Full governed claim ledger

| Claim | Headline claim | Evidence class | Boundary |
|:---|:---|:---|:---|
| CF-H01 | Fixed linear first- and second-moment queries share one typed cross-generalized geometry. | algebraic_theorem | Fixed H, K, and Gamma only; no nonlinear or adaptive statistic. |
| CF-H01 | Fixed linear first- and second-moment queries share one typed cross-generalized geometry. | internal_oracle | Independent base-matrix oracle and package property courts; this is computational evidence. |
| CF-H01 | Fixed linear first- and second-moment queries share one typed cross-generalized geometry. | external_parity | One balanced fixed-crossnobis and fixed-OLS-RSA fixture, not every RSA method. |
| CF-H02 | Total geometry decomposes exactly into coherent and configuration components and every fixed linear query inherits the split. | algebraic_theorem | Exact arithmetic does not make either finite-sample crossvalidated component nonnegative. |
| CF-H02 | Total geometry decomposes exactly into coherent and configuration components and every fixed linear query inherits the split. | internal_oracle | Independent finite-matrix construction plus production recomposition tests. |
| CF-H02 | Total geometry decomposes exactly into coherent and configuration components and every fixed linear query inherits the split. | matched_simulation | Generated planted truth validates interpretation only in the declared fixture. |
| CF-H02 | Total geometry decomposes exactly into coherent and configuration components and every fixed linear query inherits the split. | existing_illustration | Retrospective one-subject Haxby illustration; no frozen hypothesis, population inference, or neuroscience replication. |
| CF-H03 | Conservative frames preserve the declared total budget and frame-family alpha fixes each scale total. | algebraic_theorem | Total geometry under the admitted additive metric; coherent energy is not conserved across overlapping nodes. |
| CF-H03 | Conservative frames preserve the declared total budget and frame-family alpha fixes each scale total. | internal_oracle | Independent conservation and overlap-accounting courts. |
| CF-H04 | Typed observations can be compiled to an identified relation and then to geometry without changing the declared estimand. | matched_simulation | One versioned linear-model fixture; not universal BIDS, HRF, censoring, or GLM coverage. |
| CF-H05 | Query-first execution preserves the estimand while avoiding full geometry materialization in the certified regime. | internal_oracle | Recorded q=100 searchlight regime and R-heap receipt; not a cross-package speed or OS-RSS claim. |
| CF-H06 | Crossform reproduces mapped crossnobis outputs from independent implementations. | external_parity | Haxby subject 1 regression evidence and rMVPA parity on one matched squared-Euclidean estimand. |
| CF-H06 | Crossform reproduces mapped crossnobis outputs from independent implementations. | external_parity | Pinned synthetic rsatoolbox case with source binding; no correlation distance or inference parity. |
| CF-H07 | Rectangular cross-domain queries reuse the same bilinear geometry. | matched_simulation | Designed encoding-retrieval simulation with planted truth; no empirical encoding-retrieval result. |
| CF-H08 | Fixed population queries commute with declared linear transport and group modeling. | algebraic_theorem | Conditional on the realized transport; no transport-estimation uncertainty is propagated. |
| CF-H08 | Fixed population queries commute with declared linear transport and group modeling. | internal_oracle | Independent matrix oracle plus package route; no cross-node covariance estimator. |
| CF-H08 | Fixed population queries commute with declared linear transport and group modeling. | existing_illustration | Retrospective ds003745 execution/diagnostic evidence; no frozen hypothesis and no scientific transport-superiority result. |
| CF-H09 | Population point estimates, HC3 intervals, wild bootstrap, prevalence, and heterogeneity are available only under their declared conditional targets. | internal_oracle | Formula and product-oracle evidence; calibration and broad operating-characteristic claims require PE-C/PE-D simulations. |
| CF-H09 | Population point estimates, HC3 intervals, wild bootstrap, prevalence, and heterogeneity are available only under their declared conditional targets. | matched_simulation | 500 paired datasets per regime quantify classical, HC3, and wild-bootstrap behavior; transport remains realized and informative coverage supports no marginal population claim. |
| CF-H10 | Coherent/configuration mixtures can be interpreted as distinct planted regimes at fixed total energy. | matched_simulation | 48 paired seeds across three planted organizations, three noise regimes, two sample sizes, and three SNRs; synthetic line-domain evidence only, with no empirical interpretation, population coverage, or power claim. |
| CF-H10 | Coherent/configuration mixtures can be interpreted as distinct planted regimes at fixed total energy. | matched_simulation | 200 paired 24-subject hierarchical replications recover planted population ordering and profiles conditional on realized transport; informative coverage remains an unsupported failure regime and no empirical claim is licensed. |
| CF-H11 | Functionally informed transport improves population recovery relative to anatomical transport. | existing_illustration | The retrospective ds003745 eta ratio failed its interpretation audit; no superiority claim is supported. |
| CF-H12 | A future prospective real-data analysis could establish a bounded neuroscience result. | prospective_protocol | Discovery and replication specifications are frozen and synthetically rehearsed, but readiness-current is BLOCKED; no eligible real-data result or independent replication exists. |
| CF-H13 | Unsupported interpretations fail as typed refusals rather than changing estimands silently. | internal_oracle | Covers the executable refusal gallery and named gates, not every possible misuse. |
| CF-H14 | Cross-node covariance can support precision weighting or joint spatial inference. | prospective_protocol | Future capability only; current population precision and joint cross-node inference refuse. |

The following overview groups those claims by the reader’s question.
Linked guides and receipts carry the numerical details and admission
conditions.

| Topic | Evidence | What it establishes and where to read it |
|----|----|----|
| Observation workflow | Implemented and checked | An executable [observation-to-relation workflow](https://bbuchsbaum.github.io/crossform/articles/from-observations.md) agrees with direct linear-model oracles. It covers an admitted linear design, not general preprocessing or GLM coverage. |
| Shared bilinear algebra | Established algebraically | The [common-geometry guide](https://bbuchsbaum.github.io/crossform/articles/common-geometry-equivalence.md) derives and checks fixed contrast, RDM and RSA mappings. Algebra alone provides no sampling or population guarantee. |
| External point parity | Demonstrated on matched estimands | [Haxby and rsatoolbox comparisons](https://github.com/bbuchsbaum/crossform/blob/main/design/relation-to-prior-work.md) reproduce specified distances and fixed coefficients. Matching one estimand does not equate complete toolboxes. |
| Analytic uncertainty | Validated under admitted models | Residual-bearing fits support fixed-metric, within-measurement covariance and fixed linear transport. See [uncertainty interpretation](https://bbuchsbaum.github.io/crossform/articles/interpreting-results.html#uncertainty-two-targets-and-what-refused-means). |
| Selected-query execution | Measured within crossform | [Source-bound benchmarks](https://github.com/bbuchsbaum/crossform/blob/main/design/certification-report.md) compare complete and selected readouts with numerical, timing and memory receipts. No cross-package speed advantage is claimed. |
| Coherent/configuration decomposition | Demonstrated | The [introduction](https://bbuchsbaum.github.io/crossform/articles/introduction.md) checks exact recomposition and recovers two planted spatial organizations. Haxby illustrates the same reading in one participant, without inference. |
| Rectangular cross-domain forms | Demonstrated in simulation | The [encoding–retrieval exemplar](https://github.com/bbuchsbaum/crossform/blob/main/exemplars/er-rsa/README.md) handles unequal item sets, controls and covariates. Its recovered ground truth is synthetic; fixed neural precision is not admitted on that rectangular path. |
| Adjoint coupling | Implemented, small-node scope | The [evidence-pairing guide](https://bbuchsbaum.github.io/crossform/articles/evidence-pairing.md) demonstrates node/edge forms and their adjoint. It does not establish brain-scale tomography or generic connectivity inference. |
| Generalization and route identity | Implemented and tested | [Identity tests](https://github.com/bbuchsbaum/crossform/blob/main/tests/testthat/test-generalization-axis.R) distinguish sampling axes; execution-route tests preserve an estimand across fused and materialized readings. |
| Spatial conservation | Demonstrated for admitted frames | The [conservative-frame example](https://bbuchsbaum.github.io/crossform/articles/conservative-frames.md) checks that a spatial ledger adds to its global budget. This does not apply to arbitrary measurement densities or metrics. |
| Refusal discipline | Implemented and executable | The [failure gallery](https://bbuchsbaum.github.io/crossform/articles/failure-gallery.md) shows missing capabilities and unsupported requests. A refusal can check a declaration, but cannot discover hidden dependence in the observations. |
| Descriptive model coordinates | Implemented and checked | [Model-coordinate oracles](https://github.com/bbuchsbaum/crossform/blob/main/design/model-coordinate-geometry-contract.md) check lowering, the signed trace split and latent fits. Cross-fitted projector energy is distinct from predictive gain; an unregularized full-span fit does not test model shape. |
| Regularized geometry prediction | Implemented and locally validated | The [predictive guide](https://bbuchsbaum.github.io/crossform/articles/predictive-geometry.md) learns a form on training runs and scores it independently. Calibration concerns the declared same-condition, Frobenius-loss designs; new-condition transfer and generic eigenvalue inference remain outside the interface. |

Reproducible Haxby scripts and qualifications are in
[`exemplars/haxby2001`](https://github.com/bbuchsbaum/crossform/tree/main/exemplars/haxby2001).
Versioned performance evidence lives under
[`benchmarks`](https://github.com/bbuchsbaum/crossform/tree/main/benchmarks).

## What would strengthen the scientific case?

Current algebra, software tests and matched simulations support the
specific claims above. Stronger empirical claims require different
evidence:

- A prospective real-data analysis conducted under a protocol frozen
  before its outcomes are inspected, followed by an independent
  replication.
- A demonstrated interpretive benefit on that data, beyond reproducing a
  conventional distance or coefficient. Examples include meaningful
  spatial accounting or a justified comparison across unequal
  experimental axes.
- For any cross-package speed claim, a matched-estimand comparison with
  independent numerical parity and recorded hardware, runtime and
  memory.

The [evidence-status
rules](https://github.com/bbuchsbaum/crossform/blob/main/design/evidence-status-ledger.md)
keep these claims separate. A passing simulation or executable example
cannot substitute for prospective or independently replicated evidence.

## What is not claimed

`crossform` does not claim to have invented RSA, crossnobis,
cross-validated MANOVA, analytic RDM covariance, searchlights, noise
ceilings, condition/subject generalization, the regional-mean versus
pattern distinction, or connectivity. It does not claim empirical
superiority to `rsatoolbox`, a matched-estimator speed advantage,
marginal inference under informative coverage, transport-superiority
from the retrospective ds003745 illustration, cross-node covariance, or
completed prospective or independently replicated evidence. Nor does it
claim that every nonlinear RDM comparison belongs in the bilinear core.

The current contribution is a common representation whose scientific
identity, uncertainty preconditions and allowed readouts remain explicit
through execution. Its value as an analysis tool can be explored in the
[worked
introduction](https://bbuchsbaum.github.io/crossform/articles/introduction.md),
[matched spatial
simulation](https://bbuchsbaum.github.io/crossform/articles/matched-interpretability.md)
and [predictive geometry
guide](https://bbuchsbaum.github.io/crossform/articles/predictive-geometry.md).
