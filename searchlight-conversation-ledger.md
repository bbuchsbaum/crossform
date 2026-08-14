# Searchlight theory conversation ledger

Purpose: preserve the evolution, mathematical core, publication strategy, and unresolved questions of the conversation about searchlight analysis. This is a cumulative research record, not merely a polished summary.

## Recording protocol

- Preserve every incoming part as a chronological conceptual snapshot.
- Separate what the conversation **asserts** from what has been independently verified.
- Maintain a living synthesis without rewriting away earlier formulations.
- Record later refinements, contradictions, changes in scope, and abandoned ideas explicitly.
- Keep exact equations when they carry the argument; compress repetition, not substance.
- Distinguish theorem-level identities, interpretive claims, computational consequences, empirical demonstrations, and novelty claims.

## Evolution log

### 2026-08-11 — Part 1: the publishable idea

#### The conceptual turn

The initiating claim is that a searchlight should not be understood primarily as a classifier repeatedly moved through space. It is more generally:

1. a spatially localized constructor of an observation-space Gram matrix; followed by
2. a statistical readout of that matrix.

This reframing exposes a dividing line that is obscured by the usual procedural description. Some searchlights genuinely couple features and must be evaluated locally. Others merely accumulate independently computable voxelwise evidence. For the latter class, repeated neighborhood extraction and repeated fitting are algebraically unnecessary.

The memorable formulation is:

> A broad class of searchlights is exactly a voxelwise quadratic statistic followed by spatial smoothing.

#### Foundational construction: a field of local Gram matrices

Let

\[
Z=[z_1,\ldots,z_V]\in\mathbb R^{n\times V},
\]

where each column \(z_v\) is the response vector for voxel \(v\) across \(n\) observations. At centre \(s\), spatial weights \(w_{sv}\ge 0\) define

\[
D_s=\operatorname{diag}(w_{s1},\ldots,w_{sV}).
\]

The weighted local kernel is

\[
K_s=ZD_sZ^\top=\sum_v w_{sv}z_vz_v^\top.
\]

Each voxel contributes a rank-one positive-semidefinite "Gram atom"

\[
K_v=z_vz_v^\top,
\]

and the searchlight pools these atoms. On a regular grid with \(w_{sv}=h(s-v)\), every element of the Gram field is a spatial convolution:

\[
[K_s]_{ij}=\big[h*(Z_i\odot Z_j)\big](s).
\]

Thus a searchlight generates a PSD-matrix-valued spatial field, or equivalently one filtered spatial channel for every observation pair.

#### Central theorem: collapse or commutation

If the readout is linear in the local Gram matrix,

\[
T_s=\langle A,K_s\rangle_F=\operatorname{tr}(AK_s),
\]

then trace cyclicity gives

\[
T_s
=\operatorname{tr}(AZD_sZ^\top)
=\operatorname{tr}(D_sZ^\top AZ)
=\sum_v w_{sv}z_v^\top A z_v.
\]

Define voxelwise quadratic evidence

\[
q_v=z_v^\top A z_v.
\]

Then

\[
T_s=\sum_vw_{sv}q_v,
\]

and, for translation-invariant weights,

\[
T=h*q.
\]

The substantive conclusion is that spatial pooling commutes with a linear statistical readout. An analysis in this class needs one voxelwise evidence calculation and one spatial filtering operation, not a refit at every centre.

#### Cross-validation is compatible with collapse

Fixed train/test contrasts can be embedded in observation space. With folds \(\ell=1,\ldots,L\), define

\[
A_{\mathrm{CV}}
=\frac{1}{2L}\sum_\ell
\left(c_\ell^{\mathrm{tr}}c_\ell^{\mathrm{te}\top}
+c_\ell^{\mathrm{te}}c_\ell^{\mathrm{tr}\top}\right).
\]

A crossvalidated bilinear contrast energy is then

\[
T_s=\operatorname{tr}(A_{\mathrm{CV}}K_s).
\]

The fold structure changes the observation-space operator \(A\); it does not itself prevent collapse. The actual obstructions are feature coupling, adaptive estimation, or a nonlinear readout.

#### Immediate method correspondences

**Squared-distance RSA.** Squared Euclidean distances are linear in \(K_s\), so the local RDM obeys

\[
d_s=\sum_vw_{sv}d_v.
\]

With a fixed RSA regression design \(R\), coefficient estimation is also linear:

\[
\hat\beta_s=(R^\top R)^{-1}R^\top d_s
=\sum_vw_{sv}\hat\beta_v.
\]

The claim extends to multiple model and nuisance RDMs and to fixed-fold crossvalidated squared distances. It excludes nonlinear normalizations, ranks, locally estimated denominators, and dense local whitening.

**Fixed-metric MANOVA / pattern distinctness.** A crossvalidated effect numerator can be represented by \(A\). The searchlight collapses if the feature precision is globally fixed or represented by fixed diagonal weights. Standard local cvMANOVA generally does not fully collapse because each searchlight re-estimates a dense residual precision \(E_s^{-1}\). On this account, the distinctly multivariate ingredient is the locally varying covariance geometry, not cross-validation.

**HSIC, MMD, and alignment.** Linear-kernel HSIC uses \(A=HLH\), where \(L\) is a fixed target kernel and \(H\) centres observations:

\[
T_s=\operatorname{tr}(K_sHLH)
=\sum_vw_{sv}z_v^\top HLHz_v.
\]

Linear-kernel two-sample MMD has the same form with a group-defined operator. These methods, RSA, and bilinear contrast tests differ primarily in the scientific operator \(A\) applied to the same Gram field.

**Classifiers and other nonlinear readouts.** Ridge, SVM, logistic regression, and thresholded accuracy generally have

\[
T_s=F(K_s)
\]

with nonlinear \(F\), so they do not collapse to a scalar convolution. The Gram representation still offers computational structure: neighboring kernels differ through rank-one additions and removals, enabling low-rank inverse updates, or the Gram field can be constructed by matrix-valued convolution before applying \(F\).

#### Proposed taxonomy: feature-additive versus feature-interactive

The conversation rejects "univariate versus multivariate" as the decisive distinction. An analysis may pool many weak effects, encode complex contrasts, and use cross-run products while remaining feature-additive.

A feature-additive readout obeys, with suitable treatment of the affine constant,

\[
F\!\left(\sum_vw_vK_v\right)=\sum_vw_vF(K_v).
\]

A feature-interactive readout does not. The proposed pairwise interaction defect is

\[
\Delta_F(K_u,K_v)
=F(K_u+K_v)-F(K_u)-F(K_v)+F(0).
\]

Affine readouts have zero defect; nonzero defects measure algorithmic curvature under feature combination. Such curvature can arise from covariance inversion, regularization, normalized similarity, kernel nonlinearities, or thresholding. It must not automatically be interpreted as biological interaction or a population code.

#### Deeper duality: the evidence graph

Replace the diagonal window with a symmetric local feature-space operator \(M_s\):

\[
K_s=ZM_sZ^\top.
\]

Define the whole-brain evidence graph

\[
Q_A=Z^\top AZ.
\]

Then

\[
T_s
=\operatorname{tr}(AZM_sZ^\top)
=\operatorname{tr}(M_sZ^\top AZ)
=\langle M_s,Q_A\rangle_F.
\]

This separates the scientific question \(A\), the data \(Z\), and the anatomical localization/metric \(M_s\). Expanding the trace gives an exact node-edge decomposition:

\[
T_s
=\sum_vM_{s,vv}Q_{vv}
+2\sum_{u<v}M_{s,uv}Q_{uv}.
\]

Here \(Q_{vv}\) is voxelwise evidence and \(Q_{uv}\) is pairwise co-evidence. Diagonal spatial metrics inspect nodes only; dense metrics inspect nodes and edges. For local metrics the full \(V\times V\) graph need not be materialized—only spatially co-occurring pairs are required.

#### Conservation through coverage correction

Overlapping searchlights ordinarily count the same voxel multiple times. If voxel coverage is

\[
c_v=\sum_s w_{sv},
\]

define partition-of-unity weights

\[
\pi_{sv}=\frac{w_{sv}}{c_v},\qquad \sum_s\pi_{sv}=1.
\]

Then coverage-corrected local kernels satisfy

\[
\sum_sK_s^\pi=ZZ^\top=K_{\mathrm{global}},
\]

and every linear score satisfies

\[
\sum_sT_s^\pi=T_{\mathrm{global}}.
\]

This is proposed as a conservative searchlight: local geometries become additive attributions of global representational geometry rather than overlapping summaries. It addresses redundant accounting, not the separate localization blur or displacement induced by extended windows.

#### Radius reinterpreted as spatial scale

For additive evidence, a kernel family yields

\[
T_\tau=h_\tau*q.
\]

With a cortical graph Laplacian \(L\), heat-kernel localization gives

\[
T_\tau=e^{-\tau L}q,
\qquad
\frac{\partial T_\tau}{\partial\tau}=-LT_\tau.
\]

The semigroup property turns radius into diffusion time and makes multiscale searchlights one progressive scale-space analysis rather than unrelated analyses at several arbitrary radii.

#### Computational program

For the collapsed class:

\[
Z\longrightarrow q_v=z_v^\top A z_v\longrightarrow Wq.
\]

If \(A=U\Lambda U^\top\) has rank \(r\),

\[
q_v=\sum_{k=1}^r\lambda_k(u_k^\top z_v)^2,
\]

giving cost \(O(Vnr)\) before spatial filtering. Candidate filters include separable Gaussian convolution, integral images, sparse neighborhood matrices, and graph heat kernels.

For nonlinear readouts, the full Gram field requires \(n(n+1)/2\) product channels. If \(Z\approx UC\) has low observation rank, local kernels can instead be represented through \(CD_sC^\top\), reducing the filtered field to \(r(r+1)/2\) latent products.

#### Publication form proposed in Part 1

Candidate titles:

- *When Does a Searchlight Need a Searchlight? A Gram–Convolution Identity for Local Multivariate Mapping*
- *Design–Space Duality in fMRI Searchlight Analysis* (more conservative)

Five proposed formal results:

1. Gram-convolution identity.
2. Linear-readout commutation/collapse theorem.
3. Feature-additivity criterion.
4. Evidence-graph duality.
5. Coverage-normalized conservation corollary.

One proposed four-panel figure:

1. a commuting diagram for pooling and linear readout;
2. numerical identity between a conventional RSA or fixed-metric cvMANOVA searchlight and the convolution implementation;
3. a simulation separating additive evidence from covariance-sensitive/nonlinear interaction;
4. runtime scaling of repeated fitting versus quadratic-evidence-plus-filtering.

The intended paper is a short methodological note whose force comes from an exact, simple identity with unexpectedly wide scope—not from introducing another arbitrary searchlight variant.

#### Novelty position as stated, not yet independently established

The conversation acknowledges that its ingredients already exist separately: searchlights, localization critiques, second-moment unifications of representational methods, feature additivity of squared Euclidean distance, trace forms for cvMANOVA, and kernel-distance correspondences.

The asserted synthesis is novel in combining these ingredients into:

1. a PSD-matrix-valued spatial convolution;
2. a general commutation/collapse result;
3. a feature-additive versus feature-interactive taxonomy;
4. observation-space/feature-space trace duality framed as an evidence graph;
5. a partition-of-unity conservation law for local representational geometry.

Formal priority remains an open bibliographic question.

### 2026-08-11 — Part 2: the univariate confound is a projection problem

#### How Part 2 changes the project

Part 1 identified when a searchlight is a spatial filter of voxelwise quadratic evidence. Part 2 uses that operator view to replace the traditional univariate/multivariate opposition with explicit decompositions. The computational identity becomes an interpretive theory: apparent disagreements between mass-univariate and multivoxel analyses can often be located in projections, order of aggregation, feature metrics, and population estimands rather than in two ontologically different kinds of neural signal.

Five questions previously conflated under “univariate versus multivariate” are separated:

1. Is there a condition effect in the regional mean?
2. Is there a condition effect at any individual voxel?
3. Is there a reproducible spatial configuration of voxel effects?
4. Is information genuinely joint and absent from all voxelwise marginals?
5. Is the effect direction shared across subjects or only reproducible within each subject?

Mean subtraction answers only the first question by excluding one feature-space direction.

#### Exact level/configuration decomposition

For a searchlight with normalized weights \(\mathbf w\), \(\mathbf 1^\top\mathbf w=1\), let \(\Omega=\operatorname{diag}(\mathbf w)\). Given independent contrast estimates \(\hat{\boldsymbol\delta}^{(a)}\) and \(\hat{\boldsymbol\delta}^{(b)}\), total crossvalidated effect energy is

\[
E_{\mathrm{total}}
=\hat{\boldsymbol\delta}^{(a)\top}\Omega
\hat{\boldsymbol\delta}^{(b)}.
\]

The feature metric decomposes as

\[
\Omega=\mathbf w\mathbf w^\top+
(\Omega-\mathbf w\mathbf w^\top),
\]

and hence

\[
E_{\mathrm{total}}=E_{\mathrm{level}}+E_{\mathrm{configuration}},
\]

where

\[
E_{\mathrm{level}}
=(\mathbf w^\top\hat{\boldsymbol\delta}^{(a)})
(\mathbf w^\top\hat{\boldsymbol\delta}^{(b)})
\]

is reproducible evidence in the weighted regional mean, and

\[
E_{\mathrm{configuration}}
=\hat{\boldsymbol\delta}^{(a)\top}
(\Omega-\mathbf w\mathbf w^\top)
\hat{\boldsymbol\delta}^{(b)}
\]

is reproducible evidence in departures from that mean.

For equal weights,

\[
\Omega-\mathbf w\mathbf w^\top
=\frac1p\left(I-\frac1p\mathbf1\mathbf1^\top\right)
=\frac1pP_\perp.
\]

Demeaning therefore removes exactly the constant spatial mode \(\operatorname{span}\{\mathbf1\}\). It does not remove voxelwise marginal effects. A zero-mean pattern such as \((1,-1,1,-1)\) survives untouched although every voxel has an ordinary condition mean difference; a constant pattern such as \((1,1,1,1)\) is removed completely despite potentially being reliable and meaningful.

The scientifically accurate description is:

> The analysis is restricted to spatial contrasts orthogonal to the regional-average response.

It is not purified of “univariate information.”

#### Gram-kernel version

For local pattern matrix \(X\), define

\[
K_{\mathrm{total}}=X\Omega X^\top,
\qquad
\mathbf m=X\mathbf w.
\]

Then

\[
K_{\mathrm{total}}
=\underbrace{\mathbf m\mathbf m^\top}_{K_{\mathrm{level}}}
+\underbrace{X(\Omega-\mathbf w\mathbf w^\top)X^\top}_{K_{\mathrm{configuration}}}.
\]

Both summands are PSD. Consequently every linear Gram readout decomposes exactly:

\[
T_{\mathrm{total}}=T_{\mathrm{level}}+T_{\mathrm{configuration}}.
\]

An especially concise identity is

\[
T_{\mathrm{configuration}}
=\sum_vw_vx_v^\top Ax_v-(X\mathbf w)^\top A(X\mathbf w).
\]

Configuration evidence is the average voxelwise evidence minus evidence carried by the average voxel. For squared-Euclidean or suitable fixed-metric crossvalidated RSA, the same units-preserving decomposition applies to the RDM itself:

\[
\mathrm{RDM}_{\mathrm{total}}
=\mathrm{RDM}_{\mathrm{level}}
+\mathrm{RDM}_{\mathrm{configuration}}.
\]

#### Mean removal creates feature coupling

The equal-weight centered metric

\[
\frac1pP_\perp=\frac1pI-\frac1{p^2}\mathbf1\mathbf1^\top
\]

has negative off-diagonal entries. It is a dense all-to-all contrast: each voxel is defined relative to all other voxels. Thus demeaning is not a neutral deletion that otherwise leaves voxel evidence untouched. It imposes a relational spatial code through preprocessing.

This yields a crucial interpretive warning:

\[
\text{relative spatial code imposed by the analyst}
\ne
\text{condition-dependent neural interaction}.
\]

#### The univariate–multivariate gap as spatial covariance

Let \(\mathcal S_s f=\sum_vw_{sv}f_v\). Compare the additive searchlight numerator

\[
E_{\mathrm{SL}}(s)
=\mathcal S_s(\hat{\boldsymbol\delta}^{(a)}\odot
\hat{\boldsymbol\delta}^{(b)})
\]

with evidence in the locally averaged contrast,

\[
E_{\mathrm{mean}}(s)
=\mathcal S_s\hat{\boldsymbol\delta}^{(a)}\,
\mathcal S_s\hat{\boldsymbol\delta}^{(b)}.
\]

Their exact difference is

\[
E_{\mathrm{SL}}(s)-E_{\mathrm{mean}}(s)
=\operatorname{Cov}_{w_s}
(\hat{\boldsymbol\delta}^{(a)},
\hat{\boldsymbol\delta}^{(b)}).
\]

With identical rather than independent maps, this becomes the weighted variance/Jensen gap:

\[
\mathcal S_s(\delta^2)-(\mathcal S_s\delta)^2
=\operatorname{Var}_{w_s}(\delta).
\]

This locates a major source of apparent disagreement. Local averaging rewards a spatially coherent signed effect; additive effect energy retains stable positive and negative voxel effects. Demeaned additive mapping extracts their reproducible spatial heterogeneity. The contrast can be summarized as “pool then square” versus “square then pool,” whose difference is spatial variance—not necessarily a distinct neural ontology.

#### Smooth-field implication: demeaned searchlights can act as edge detectors

For a smooth true contrast field and a symmetric window with spatial covariance \(C_w\), a local Taylor expansion gives, to leading order,

\[
E_{\mathrm{configuration}}(s)
\approx
\nabla\delta^{(a)}(s)^\top C_w\nabla\delta^{(b)}(s).
\]

With an isotropic window \(C_w=\tau^2I\), expected configuration evidence is proportional to \(\|\nabla\delta(s)\|^2\). Demeaning can therefore transform an ordinary smooth activation blob into low evidence at its flat centre and a high-valued ring or shell at its boundary. It removes the zero-order/DC component but preserves gradients, curvature, gain variation, and higher spatial frequencies.

#### Confound control belongs in the correct space

Feature-space mean removal is not a general control for reaction time, difficulty, arousal, motion, or other experimental nuisance variables. Such nuisances will generally have heterogeneous voxel gains, so demeaning removes only their constant component.

The general statistic

\[
T(A,M)=\operatorname{tr}(AXMX^\top)
\]

separates two intervention sites:

- \(A\), in observation/design space, defines which relationships among trials or conditions count as target evidence.
- \(M\), in feature/voxel space, defines which spatial directions and voxel relationships count as evidence.

For nuisance design \(N\), observation-space residualization can use

\[
R_N=I-N(N^\top N)^{-1}N^\top,
\qquad
A_{\mathrm{target}\mid N}
=R_NA_{\mathrm{target}}R_N.
\]

Feature invariances such as removing the regional mean or prescribed spatial basis functions instead modify \(M\). Concisely: voxel demeaning is a right-side projection; reaction-time regression is a left-side projection. They are not interchangeable.

#### Metric-aware “above regional mean” decomposition

Ordinary demeaning is orthogonal only in Euclidean geometry. If the analysis uses feature metric \(M\), for example a precision matrix, a nuisance feature subspace \(U\) should be decomposed in that same geometry:

\[
M_U=MU(U^\top MU)^{-1}U^\top M,
\qquad
M_{\perp U}=M-M_U.
\]

Then

\[
\boldsymbol\delta_a^\top M\boldsymbol\delta_b
=\boldsymbol\delta_a^\top M_U\boldsymbol\delta_b
+\boldsymbol\delta_a^\top M_{\perp U}\boldsymbol\delta_b.
\]

For \(U=\mathbf1\), this is the proper level/configuration split in the classifier’s own noise geometry. Euclidean mean subtraction followed by a differently whitened analysis can mix the components again.

#### Why raw-minus-demeaned accuracy is not attribution

Classification accuracy is nonlinear and the classifier is retrained after projection. Therefore a change such as 70% to 66% cannot be interpreted as four percentage points of “univariate contribution.” Signals may be redundant; removing a noisy mean may improve prediction; regularization changes; and accuracy saturates.

For nonlinear readout \(F\), the relevant algorithm-specific interaction is

\[
\Delta_F
=F(K_{\mathrm{level}}+K_{\mathrm{configuration}})
-F(K_{\mathrm{level}})
-F(K_{\mathrm{configuration}})
+F(0).
\]

This can characterize synergy or redundancy under that readout but not, without further argument, neural interaction. The preferred primary decomposition uses an additive crossvalidated evidence scale, with accuracy as a secondary predictive summary.

#### Across-subject decomposition: consensus versus individuality

Let subject-specific true contrast patterns follow

\[
\boldsymbol\delta_i=\boldsymbol\mu+\boldsymbol\eta_i,
\qquad
\operatorname{Cov}(\boldsymbol\delta_i)=\Sigma_B.
\]

A mass-univariate group contrast primarily tests coordinates of the shared signed mean \(\boldsymbol\mu\). A within-subject multivariate energy targets \(\boldsymbol\delta_i^\top M\boldsymbol\delta_i\). Across subjects,

\[
E_i[\boldsymbol\delta_i^\top M\boldsymbol\delta_i]
=\boldsymbol\mu^\top M\boldsymbol\mu
+\operatorname{tr}(M\Sigma_B).
\]

The first term is consensus/shared evidence; the second is subject-specific heterogeneity energy. Thus reliable within-subject decoding can coexist with zero population-mean topography when subject effects differ in sign or location. This discrepancy does not by itself imply a higher-dimensional neural code.

#### The two-by-two scientific decomposition

Combining \(M=M_{\mathrm{level}}+M_{\mathrm{configuration}}\) with the consensus/heterogeneity split yields four distinct targets:

| Spatial component | Shared across subjects | Subject-specific |
|---|---|---|
| Regional level | \(\boldsymbol\mu^\top M_{\mathrm{level}}\boldsymbol\mu\) | \(\operatorname{tr}(M_{\mathrm{level}}\Sigma_B)\) |
| Spatial configuration | \(\boldsymbol\mu^\top M_{\mathrm{configuration}}\boldsymbol\mu\) | \(\operatorname{tr}(M_{\mathrm{configuration}}\Sigma_B)\) |

These mean, respectively: shared regional activation; shared multivoxel topography; reliable but sign/magnitude-variable regional effects; and individually reliable but non-aligned spatial topographies.

#### Estimating the four components

With independent within-subject contrast estimates,

\[
\widehat W_M
=\frac1N\sum_i
\hat{\boldsymbol\delta}_i^{(a)\top}M
\hat{\boldsymbol\delta}_i^{(b)}
\]

estimates total within-subject energy. In a common aligned feature space, the cross-subject U-statistic

\[
\widehat C_M
=\frac{1}{N(N-1)}\sum_{i\ne j}
\hat{\boldsymbol\delta}_i^\top M
\hat{\boldsymbol\delta}_j
\]

estimates consensus energy \(\boldsymbol\mu^\top M\boldsymbol\mu\). Their difference

\[
\widehat H_M=\widehat W_M-\widehat C_M
\]

estimates heterogeneity energy \(\operatorname{tr}(M\Sigma_B)\). Applying this separately to level and configuration metrics produces four maps.

This suggests a new criterion for functional alignment: successful alignment should transfer energy from idiosyncratic configuration to shared configuration while approximately preserving total within-subject information.

#### Configuration is not the same as genuinely joint information

A mean-zero spatial contrast can remain a set of marginal voxelwise mean effects. Strictly joint information requires matched voxelwise marginals across conditions but a difference in the joint distribution, such as condition-dependent cross-voxel covariance or higher-order structure.

Part 2 proposes the following hierarchy:

1. regional-level information in the constant spatial mode;
2. multiple marginal voxel effects;
3. configuration information from nonconstant marginal effects;
4. metric-coupled information, where a dense geometry changes their combination;
5. marginally silent joint information in covariance, copula, or higher-order structure.

Mean removal moves attention from level 1 toward level 3. It does not establish level 5.

#### From a binary label to a spatial information spectrum

Let \(L_s\) be a local spatial graph Laplacian with eigenvectors \(u_k\), ordered by eigenvalue, with \(u_0\propto\mathbf1\). Effect energy decomposes spectrally:

\[
\hat{\boldsymbol\delta}^{(a)\top}
\hat{\boldsymbol\delta}^{(b)}
=\sum_k
(u_k^\top\hat{\boldsymbol\delta}^{(a)})
(u_k^\top\hat{\boldsymbol\delta}^{(b)}).
\]

The modes distinguish DC/regional level, smooth gradients, mesoscale topography, and fine spatial texture. Mean removal discards only \(k=0\). The proposed scientific question becomes:

> At what spatial modes, under what feature metric, and with what degree of population consensus is the experimental information expressed?

#### Unified workflow proposed in Part 2

1. Estimate independent runwise or split-half contrast patterns.
2. Compute fixed-metric crossvalidated effect energy rather than beginning with accuracy.
3. Decompose total energy into level and configuration in the same metric.
4. Optionally split configuration into smooth, mesoscale, and fine spectral bands.
5. Decompose each spatial component into consensus and subject-specific energy.
6. Add nonlinear classification as a secondary readout and investigate any increment as metric weighting, regularization, nonlinear decision structure, variance sensitivity, or genuine joint information.
7. Control behavioral and design confounds in observation space and spatial invariances in feature space.

The integrated target is

\[
\text{total evidence}
=\text{shared level}
+\text{shared configuration}
+\text{idiosyncratic level}
+\text{idiosyncratic configuration},
\]

with configuration optionally resolved by spatial frequency.

#### Publication centre proposed in Part 2

Part 2 proposes a stronger, more interpretive title:

*The Univariate–Multivariate Gap Is a Spatial Covariance: An Exact Decomposition of fMRI Searchlight Evidence into Regional Level, Spatial Configuration, and Population Heterogeneity.*

Its first central identity is the weighted product/covariance decomposition:

\[
\mathcal S_s(\hat{\boldsymbol\delta}^{(a)}\odot
\hat{\boldsymbol\delta}^{(b)})
=\mathcal S_s\hat{\boldsymbol\delta}^{(a)}\,
\mathcal S_s\hat{\boldsymbol\delta}^{(b)}
+\operatorname{Cov}_{w_s}
(\hat{\boldsymbol\delta}^{(a)},
\hat{\boldsymbol\delta}^{(b)}).
\]

Its second is the population second-moment decomposition:

\[
E_i[\boldsymbol\delta_i^\top M\boldsymbol\delta_i]
=\boldsymbol\mu^\top M\boldsymbol\mu
+\operatorname{tr}(M\Sigma_B).
\]

Together they explain disagreement along two axes: level versus configuration within brains, and consensus versus individuality across brains.

#### Decisive conclusion of Part 2

“Univariate” and “multivariate” usually name different projections, aggregations, metrics, nonlinearities, and population estimands applied to an underlying contrast field—not two cleanly separated kinds of neural signal. The proposed replacement is

\[
\text{level}
\oplus
\text{configuration}
\oplus
\text{population heterogeneity}
\oplus
\text{joint-distribution interaction}.
\]

### 2026-08-11 — Part 3: from statistical maps to Unified Effect Geometry

#### Threshold crossed

Part 3 explicitly enlarges the ambition. The project is no longer framed chiefly as a reformulation of searchlights. It proposes a common conceptual architecture for task-fMRI analysis, provisionally named **Unified Effect Geometry (UEG)**.

The historical diagnosis is that the GLM unified experimental designs while retaining the voxelwise statistical map as the basic spatial observable. Searchlights relaxed locality but positioned information mapping against activation mapping. cvMANOVA brought factorial design logic into multivoxel analysis, and second-moment methods unified encoding, PCM, and RSA. UEG proposes the missing synthesis: experimental design, spatial measurement, and population/generalization structure become coequal operators acting on a common effect object.

Its foundational claim is:

> An fMRI analysis is not inherently univariate, multivariate, local, global, fixed-effect, or random-effect. It measures a common effect geometry using specified experimental, spatial, and generalization operators.

The familiar categories become presets rather than distinct analysis species.

#### Primitive object: the effect matrix

For participant \(i\) and independent partition \(r\), retain

\[
B_{ir}\in\mathbb R^{q\times V},
\]

where rows span design coordinates and columns span spatial features. Conventional pipelines immediately collapse this matrix in different directions: SPM selects a row contrast and inspects columns, RSA pools columns and studies row geometry, searchlights repeatedly subset columns, classifiers learn spatial projections, and group analysis collapses participants afterward. UEG retains \(B_{ir}\) as a queryable object.

Let \(a=(i,r)\) index participant-partition units and introduce:

\[
A\in\mathbb R^{q\times q},\qquad
M\in\mathbb R^{V\times V},\qquad
\Gamma\in\mathbb R^{m\times m}.
\]

Their meanings are:

- \(A\): which experimental/design relationships count as evidence;
- \(M\): which spatial patterns, scales, supports, or feature metrics count;
- \(\Gamma\): between which runs, subjects, tasks, sites, or other units the effect must generalize.

The atomic cross-unit comparison is

\[
\mathcal E_{ab}(A,M)
=\operatorname{tr}(AB_aMB_b^\top).
\]

The proposed master statistic is

\[
T(\Gamma,A,M)
=\sum_{a,b}\Gamma_{ab}\operatorname{tr}(AB_aMB_b^\top).
\]

If the \(B_a\) are stacked vertically into \(\mathbf B\),

\[
\boxed{
T(\Gamma,A,M)
=\operatorname{tr}\left[(\Gamma\otimes A)
\mathbf B M\mathbf B^\top\right].
}
\]

UEG calls \(A\), \(M\), and \(\Gamma\) the design, brain, and generalization geometries. This is the new master equation around which all earlier identities are reorganized.

#### Evidence-operator duality generalized

Trace cyclicity defines the spatial effect operator

\[
Q_{\Gamma,A}
=\mathbf B^\top(\Gamma\otimes A)\mathbf B,
\]

so that

\[
T(\Gamma,A,M)=\langle M,Q_{\Gamma,A}\rangle_F.
\]

For fixed \(\Gamma\) and \(A\), voxelwise activation, smoothing, ROI summaries, searchlights, spatial modes, and global linear analyses become different measurements of this spatial operator. The dual participant-partition-design Gram object is

\[
K_M=\mathbf B M\mathbf B^\top,
\]

and

\[
T(\Gamma,A,M)
=\langle M,Q_{\Gamma,A}\rangle_F
=\langle\Gamma\otimes A,K_M\rangle_F.
\]

The conceptual claim is that spatial activation geometry and representational geometry are right- and left-Gram views of the same effect matrix.

#### First and second moments in one homogeneous effect state

The quadratic master equation does not by itself contain a signed linear activation contrast. Part 3 therefore embeds first- and second-order population information into a homogeneous moment matrix. With \(b_i=\operatorname{vec}(B_i)\), define

\[
\mu=E[b_i],\qquad G=E[b_ib_i^\top],
\]

and

\[
\mathbb H
=E\left[
\begin{pmatrix}1\\b_i\end{pmatrix}
\begin{pmatrix}1&b_i^\top\end{pmatrix}
\right]
=
\begin{pmatrix}
1&\mu^\top\\
\mu&G
\end{pmatrix}.
\]

This PSD effect state contains the signed mean \(\mu\), total second moment \(G\), and population covariance \(\Sigma_B=G-\mu\mu^\top\). Any degree-two polynomial query is a linear functional of \(\mathbb H\). A group contrast reads the mean block, effect-energy/RSA/MVPA queries read structured parts of the second moment, and standard errors or heterogeneity derive from their difference.

For linear query \(\ell\),

\[
t(\ell)
=\frac{\sqrt N\,\ell^\top\mu}
{\sqrt{\ell^\top\Sigma_B\ell}}.
\]

A voxelwise test uses \(\ell=c\otimes e_v\); a multivariate energy test uses a higher-rank quadratic operator such as \(M\otimes A\). The slogan is that SPM and MVPA read different blocks of the same parent effect state.

#### “Univariate” is coordinate-dependent

Under an invertible spatial change of basis \(B'=BR\), the same scientific geometry is preserved by the contravariant metric transformation

\[
M'=R^{-1}MR^{-\top},
\qquad
B'M'B'^\top=BMB^\top.
\]

An effect occupying many voxel coordinates can occupy one rotated coordinate, and vice versa. Therefore “univariate” is not a basis-invariant property of the underlying effect. Anatomical coordinates remain scientifically privileged when locality matters, but fundamental descriptions should concern invariant or explicitly geometry-relative properties such as rank, effective rank, generalized eigenvalues, norms, angles, subspace overlap, cortical smoothness, spatial support, and population consensus.

#### Activation–representation singular-value duality

For \(A,M\succeq0\), define

\[
Z=A^{1/2}BM^{1/2}=U\Sigma V^\top.
\]

Then

\[
ZZ^\top
=A^{1/2}BMB^\top A^{1/2}
=U\Sigma^2U^\top
\]

is the design-side/representational geometry, while

\[
Z^\top Z
=M^{1/2}B^\top ABM^{1/2}
=V\Sigma^2V^\top
\]

is the spatial activation geometry. They have identical nonzero eigenvalues, with paired experimental and spatial modes \((u_k,v_k)\) of strength \(\sigma_k\).

Traditional contrast mapping preselects an experimental direction; RSA preselects a spatial metric and studies design geometry. UEG proposes reporting their paired modes.

For a binary linear mean contrast \(\delta=c^\top B\), the spatial signal operator \(Q_c=\delta^\top\delta\) has rank one. Thus distributing a two-condition mean effect across many voxels does not make its signal subspace multidimensional. Higher-dimensional signal requires multiple independent design contrasts, heterogeneity spanning multiple modes, higher-order moment effects, or nonlinear feature lifting.

#### Coherent versus incoherent local integration

Part 3 renames and strengthens Part 2’s covariance decomposition. For normalized nonnegative spatial kernel \(h\), with \(D_h=\operatorname{diag}(h)\),

\[
D_h=hh^\top+(D_h-hh^\top),
\qquad D_h-hh^\top\succeq0.
\]

A crossvalidated smoothed local mean effect is

\[
E_{\mathrm{SPM}}
=(h^\top\delta^{(a)})(h^\top\delta^{(b)})
=\delta^{(a)\top}hh^\top\delta^{(b)},
\]

whereas additive searchlight energy is

\[
E_{\mathrm{SL}}
=\delta^{(a)\top}D_h\delta^{(b)}.
\]

Therefore

\[
E_{\mathrm{SL}}
=E_{\mathrm{coherent}}+E_{\mathrm{configuration}}.
\]

The proposed physical analogy is:

- SPM-like smoothing performs coherent spatial integration: signed effects align before accumulating.
- An additive searchlight performs incoherent energy integration: reproducible positive and negative effects contribute without cancellation.

The distinction is coherent local effect versus total local effect energy, not univariate signal versus multivariate information. Demeaning selects the configuration component rather than removing every coordinatewise effect.

For a latent or denoised effect vector, Part 3 introduces local coherence fraction

\[
\kappa_h
=\frac{(h^\top\delta)^2}{\delta^\top D_h\delta},
\qquad 0\le\kappa_h\le1.
\]

Its complement is the fraction of energy in spatial contrasts. Estimation from noisy crossvalidated quantities will require care because ratios of signed/near-zero estimates need not inherit the population bounds.

#### Activation blob and boundary as spectral components of one field

The Taylor expansion from Part 2 is retained and reinterpreted spectrally. Demeaned configuration energy of a smooth field has leading term

\[
\nabla\delta^{(a)}(s)^\top C_h\nabla\delta^{(b)}(s),
\]

so the boundary highlighted by a demeaned searchlight can be the spatial derivative energy of the coherent activation field. SPM emphasizes its low-frequency common mode; demeaning highlights gradients and curvature. These are bands of one spatial effect, not competing neural ontologies.

#### Universal activation spectrum

Given signal operator \(Q\) and noise/variability operator \(N\succ0\), define

\[
R(w)=\frac{w^\top Qw}{w^\top Nw}.
\]

The generalized eigenproblem

\[
Qw_k=\lambda_kNw_k
\]

or spectrum of \(N^{-1/2}QN^{-1/2}\) defines the proposed **activation spectrum**.

Traditional analyses become constrained queries of this spectrum:

| Analysis | Spatial constraint or summary |
|---|---|
| Unsmoothed mass-univariate | Coordinate direction \(w=e_v\) |
| Smoothed mass-univariate | Predetermined local mode \(w=h_s\) |
| Local MVPA | Optimize \(w\) with support in neighborhood \(\mathcal N_s\) |
| Whole-brain MVPA | Optimize over a global structured class |
| MANOVA/distinctness | Sum or transform several generalized eigenvalues |

Possible spectral summaries are trace-like \(\sum\lambda_k\), log-determinant \(\sum\log(1+\lambda_k)\), bounded explained-energy \(\sum\lambda_k/(1+\lambda_k)\), or largest-root statistics. The claimed unification is that voxelwise SPM, smoothed SPM, searchlight MVPA, global MVPA, and MANOVA measure the same signal–noise geometry under different constraints and summaries.

#### Decomposing “multivariate gain”

Part 3 proposes replacing a vague performance advantage with nested relaxations:

1. **Pooling gain:** pointwise to predetermined coherent local mean.
2. **Cancellation-avoidance gain:** \(hh^\top\) to diagonal energy \(D_h\).
3. **Relational/covariance gain:** diagonal to dense feature metric.
4. **Adaptive-direction gain:** fixed metric to a learned rank-one discriminant operator.
5. **Scale gain:** local to distributed/global support.
6. **Nonlinear-order gain:** first-order features to covariance or higher-order feature lifts.

Each relaxation must be trained on training data and evaluated independently. The intended result is a precise statement of which freedom produced predictive improvement.

#### Decoder and activation modes

In \(Qw=\lambda Nw\), \(w\) is the discriminative mode shaped by noise suppression. For rank-one signal \(Q=aa^\top\), the forward effect mode satisfies \(a\propto Nw\). UEG therefore makes the Haufe-style distinction structural: a decoder weight should be accompanied by its corresponding forward activation/effect mode. They are paired objects within the same signal–noise geometry, not interchangeable maps.

#### Population geometry becomes a first-class operator

Let

\[
B_i=\bar B+U_i,
\qquad E[U_i]=0.
\]

For \(A\succeq0\), total reproducible within-person geometry decomposes as

\[
Q_A^W
=E[B_i^\top AB_i]
=\underbrace{\bar B^\top A\bar B}_{Q_A^C}
+\underbrace{E[U_i^\top AU_i]}_{Q_A^H}.
\]

Thus

\[
Q_A^W=Q_A^C+Q_A^H,
\]

where \(C\) denotes population consensus and \(H\) reliable heterogeneity. This is the operator-valued extension of mean square equals squared mean plus variance.

With noisy independent partitions \(B_{ir}=B_i+\varepsilon_{ir}\), proposed cross-fitted estimators use symmetrized cross-partition products for \(\widehat Q_A^W\), cross-participant products for \(\widehat Q_A^C\), and their difference for \(\widehat Q_A^H\). Unbiased crossvalidated estimates can be indefinite in finite samples; eigenvalue truncation should not be used silently for inference, though separate shrinkage/PSD estimates may support descriptive mode visualization.

#### Group inference as a constrained population eigenproblem

For one contrast with population mean \(\mu\) and between-subject covariance \(\Sigma_B\),

\[
t_v^2
=N\frac{e_v^\top\mu\mu^\top e_v}
{e_v^\top\Sigma_Be_v}.
\]

This is a coordinate-restricted Rayleigh quotient. Its regularized global counterpart is \(N\mu^\top\Sigma_B^{-1}\mu\); a group searchlight restricts the direction to a neighborhood. The asserted unification is that voxelwise group SPM, local population MVPA, and global population MVPA ask the same consensus-versus-heterogeneity question along coordinate, fixed regional, learned local, or structured global directions.

#### Generalization as an explicit geometry

UEG replaces procedural fold language with \(\Gamma\), whose vertices may be runs, sessions, participants, tasks, stimulus sets, scanners, sites, age groups, or time points, and whose weighted edges state where agreement is required.

Examples include within-person cross-run stability, cross-session reliability, cross-participant consensus, cross-task abstraction, cross-item generalization, and cross-site robustness. Crossed forms can be represented through tensor products such as

\[
\Gamma_{\mathrm{subject}}
\otimes\Gamma_{\mathrm{item}}
\otimes\Gamma_{\mathrm{session}}.
\]

This changes the reported estimand from “fivefold accuracy” to effect geometry generalizing across named scientific axes. It also offers a route to treating participants and stimuli as crossed sampling dimensions.

#### Covariates, longitudinal structure, and “second level”

With subject model

\[
B_i=B_0+x_iG_x+U_i,
\]

the signed covariate effect is \(G_x\), and its geometry is \(G_x^\top AG_x\). Alternatively, a nuisance-residualized subject kernel \(\Gamma_x\propto\tilde x\tilde x^\top\) can encode the covariate query. Richer participant–occasion operators cover longitudinal effects. In UEG, the second level is not a separate species of pipeline but an operator on a population/occasion axis of the effect tensor.

#### Complete local-population-scale description

Combining population geometry with \(M_{\mathrm{total}}=M_{\mathrm{level}}+M_{\mathrm{configuration}}\) yields the four Part 2 components as operator queries: shared level, shared configuration, idiosyncratic level, and idiosyncratic configuration. Further decomposing configuration into graph-Laplacian, wavelet, or cortical-frequency bands yields a result indexed by

\[
\text{experimental effect}
\times
\text{population component}
\times
\text{spatial mode/scale}.
\]

The hoped-for reporting language becomes quantitative—for example, the fraction of reproducible energy in a shared low-frequency configuration mode—rather than separate “activation” and “decoding” declarations.

#### Conserved local/global measurement and activation tomography

For PSD spatial operators satisfying a resolution of identity,

\[
\sum_jM_j=I,
\]

linearity gives

\[
\sum_jT(\Gamma,A,M_j)=T(\Gamma,A,I).
\]

Orthogonal bases, tight frames, spectral bands, coverage-normalized parcels, and corrected searchlights can therefore provide local or multiscale measurements that sum to global evidence.

If \(Q\succeq0\) is normalized by positive trace, then \(p_j=\operatorname{tr}(M_jQ)\) forms a conserved distribution over the measurement frame. This connects the framework to resolution-of-identity and signal-tomography mathematics. Part 3 names a possible method **activation tomography**: estimate a common global effect operator and obtain controlled local projections, potentially using a rich operator frame capable of recovering node, edge, and subspace structure rather than merely painting centre scores.

This probability-like interpretation is conditional on PSD signal operators; crossvalidated or contrast-weighted effect operators can be indefinite.

#### Spatial scale and global regularization

Graph heat kernels \(H_\tau=e^{-\tau L}\) replace arbitrary radius with diffusion scale, while graph-spectral projectors \(P_b\) partition conserved evidence across common, broad-gradient, mesoscale, and fine modes. Local and global MVPA become members of a continuous spatial-operator family: a searchlight imposes hard support, whereas a whole-brain model uses soft smoothness, sparsity, anatomical, or low-rank regularization. A global decoder selects a rank-one operator \(ww^\top\); regularized learning optimizes a quotient such as

\[
\frac{w^\top Qw}{w^\top Nw+\lambda w^\top Lw}.
\]

#### Native-space population comparison and alignment

If participant \(i\) has native feature dimension \(V_i\), a subject-specific map \(P_i\in\mathbb R^{V_i\times d}\) can connect a shared anatomical, functional, or latent space to native features. Cross-person geometry can then be evaluated through terms of the form

\[
\operatorname{tr}
\left(AB_iP_iM_\star P_j^\top B_j^\top\right)
\]

without forcing identical native voxel coordinates. Alignment is defined as transferring reproducible energy from heterogeneity geometry \(Q_H\) into consensus geometry \(Q_C\) while preserving within-person geometry \(Q_W\). This becomes a testable objective rather than an informal expectation.

#### Joint information is indexed by moment order

A spatially nonconstant first-order contrast remains a vector of marginal mean effects. To target joint information with matched means, lift observations to second order:

\[
\phi_2(X)=\operatorname{vec}(XX^\top),
\qquad
\Delta C
=E[XX^\top\mid Y=1]-E[XX^\top\mid Y=0].
\]

A covariance-only effect has \(\operatorname{diag}(\Delta C)=0\) but \(\Delta C\ne0\). Higher-order effects use \(\phi_k(X)=X^{\otimes k}\), with kernels providing implicit lifts. The proposed replacement descriptors are

\[
\text{moment order}
\times
\text{spatial support}
\times
\text{spatial rank}
\times
\text{adaptivity}
\times
\text{population generalization}.
\]

#### Activation, representation, and connectivity

For time-by-space data \(Y\), spatial operators of the form \(Y^\top AY\) vary with observation/time operator \(A\): a task-effect rank-one \(A\) produces an outer product of a task contrast, while temporal centering \(A=H_T\) produces an unnormalized temporal covariance. Condition-weighted, localized, lagged, or frequency-selective operators can represent task-modulated or dynamic connectivity. Cross-participant \(\Gamma\) can express intersubject coupling.

The proposed unification is that activation and covariance/connectivity arise from different observation-space operators on common data. Normalized correlations and other nonlinear summaries remain additional readouts rather than pure linear operator queries.

#### The first-level GLM composes into the design operator

If first-level estimation is linear, \(B=LY\), then

\[
B^\top AB=Y^\top L^\top ALY.
\]

Design, nuisance regression, filtering, prewhitening, HRF bases, and contrasts can therefore compile into an observation-space operator \(A_Y=L^\top AL\). Beta matrices remain convenient sufficient summaries under appropriate models, but first-level estimation and later geometry form one composed operator chain. This makes preprocessing choices explicit parts of the estimand: temporal nuisance operations act on the left, while smoothing, demeaning, and spatial whitening act on the right.

#### Generative counterpart

Part 3 proposes a separable variance-component model

\[
\operatorname{Cov}[\operatorname{vec}(B)]
=\sum_k\theta_k
(\Gamma_k\otimes A_k\otimes M_k),
\]

estimated through REML, hierarchical Bayes, score tests, or crossvalidated components. This supplies a possible generative and inferential counterpart to the descriptive quadratic calculus, including uncertainty, shrinkage, and covariance-component tests.

#### Computational feasibility

The formal \(V\times V\) spatial operator can remain matrix-free:

\[
Q_{Gamma,A}x
=\mathbf B^\top(\Gamma\otimes A)(\mathbf Bx).
\]

Its rank is bounded by the rank of \(\Gamma\otimes A\), hence by \(\operatorname{rank}(\Gamma)\operatorname{rank}(A)\), subject also to the rank and dimensions of \(\mathbf B\). Low-rank design queries permit dual eigensolutions in participant-design space. Diagonal operators retain the convolution shortcut; low-rank spatial operators use projections; local dense metrics can exploit updates; latent and wavelet representations reduce dimension or sparsify multiscale computation.

#### Proposed UEG output

For each experimental effect, the primary result becomes an effect state rather than a thresholded map:

\[
\mathcal A_A=
(\bar B,Q_A^W,Q_A^C,Q_A^H,
N_{\mathrm{measurement}},N_{\mathrm{population}}).
\]

Derived summaries include:

- signed consensus direction \(\bar B\);
- reproducible within-person, consensus, and heterogeneity geometries;
- local coherence \(\kappa\);
- consensus fraction \(\pi_C(M)=\operatorname{tr}(MQ_A^C)/\operatorname{tr}(MQ_A^W)\);
- effective dimensionality \((\sum_k\lambda_k)^2/\sum_k\lambda_k^2\) for nonnegative spectra;
- spatial scale energies \(\operatorname{tr}(P_bQ_A)\);
- paired experimental/spatial singular modes;
- explicit cross-run, cross-session, cross-person, cross-task, cross-stimulus, and cross-site generalization profiles.

#### What UEG would and would not replace

UEG does not discard motion correction, temporal modeling, the GLM, permutation inference, or existing software backends. It seeks to replace the coordinate-first ontology that treats voxelwise, searchlight, ROI, representational, whole-brain, and second-level workflows as separate top-level analysis types.

Practices targeted for demotion or correction include generic mean subtraction as “univariate removal,” classifier accuracy as primitive effect size, decoder weights as activation maps, one arbitrary radius, separate local/global paradigms, ordinary group tests on chance-bounded accuracies, searchlight-only claims of multidimensionality, and mandatory exact voxel correspondence.

The framework explicitly does not solve design validity, causality, vascular/motion/behavioral confounding, HRF misspecification, alignment validity, unstable learned modes, or sampling bias. It names where those assumptions enter: \(A\), \(M\), \(\Gamma\), the noise models, and the study design.

#### Five-proposition mature framework

Part 3 condenses UEG to:

1. **Universal operator form:** a broad fixed-metric linear/quadratic class is \(T(\Gamma,A,M)\).
2. **Activation–representation duality:** left and right Gram geometries share nonzero spectra and paired singular modes.
3. **Coherence decomposition:** additive local energy equals coherent regional effect plus configuration energy.
4. **Population decomposition:** within-person geometry equals consensus plus heterogeneity geometry.
5. **Conservation:** spatial resolutions of identity make local/multiscale evidence sum to global evidence.

#### Publication program after Part 3

The conversation now distinguishes two papers:

- A focused note, *The Univariate–Multivariate Gap Is a Spatial Covariance*, covering coherent/incoherent integration, exact level/configuration decomposition, gradient-energy interpretation, and the convolution shortcut.
- A larger framework paper, *Activation Is an Operator: Unified Effect Geometry for fMRI*, centred on the master equation, homogeneous effect state, dual Gram geometry, activation spectrum, population/generalization operators, and conserved spatial measurement.

The decisive empirical demonstration is proposed as five regimes:

1. a coherent activation blob;
2. a zero-mean gradient composed of ordinary marginal voxel effects;
3. a shared multiscale configuration;
4. individually reproducible but mutually rotated participant patterns;
5. a covariance-only effect with matched voxelwise marginals.

Conventional SPM, demeaned searchlights, local MANOVA, and global decoding should appear discordant, while UEG should recover them as predictable allocations of effect geometry under explicit operator choices.

### 2026-08-11 — Part 4: no forward, no reverse, no common grid

#### New primitive: the brain–experiment relation

Part 4 argues that even the effect matrix is a coordinate realization of a deeper object: a relation between experimental space and native neural space. Encoding and decoding are opposing conditional factorizations of this relation, while activation and representation are its two Gram geometries. This move also changes population inference: invariant experimental geometry can be compared across participants without first equating native voxel bases.

The proposed conceptual chain is:

\[
\text{conditional density relation}
\rightarrow
\text{effect operator}
\rightarrow
\text{relation spectrum}
\rightarrow
\text{gauge-invariant geometry}
\rightarrow
\text{native-space population field}.
\]

#### Exact direction-free evidence identity

For experimental variable \(D\), neural response \(Y\), and nuisance/context \(Z\), define

\[
r(d,y\mid z)=\frac{p(d,y\mid z)}{p(d\mid z)p(y\mid z)}.
\]

Bayes’ rule gives the pointwise identity

\[
\frac{p(y\mid d,z)}{p(y\mid z)}
=\frac{p(d\mid y,z)}{p(d\mid z)}
=r(d,y\mid z),
\]

and therefore

\[
\log\frac{p(y\mid d,z)}{p(y\mid z)}
=\log\frac{p(d\mid y,z)}{p(d\mid z)}.
\]

Its expectation is conditional mutual information \(I(D;Y\mid Z)\). Encoding and decoding information gains are thus the same evidence ratio under one coherent joint model. This is a statement about statistical association/factorization, not causal direction: randomized design can still cause neural response.

For a fixed design, \(p(d\mid z)\) is interpreted as the explicit design measure over trials, conditions, items, or features rather than a metaphysical claim of random treatment labels.

#### Canonical score: held-out log information gain

Part 4 proposes scoring independently evaluated conditional models by log predictive density improvement over a nuisance-only baseline:

\[
\mathcal I^{Y\leftarrow D}_{a\to b}
=\frac1{n_b}\sum_{t\in b}
[\log q_a(y_t\mid d_t,z_t)-\log q_{0,a}(y_t\mid z_t)],
\]

with an analogous decoder score. If the two conditionals and marginals come from the same coherent joint model, their pointwise gains agree. Units are nats per observation, or bits after division by \(\log2\).

Accuracy, AUC, correlation, \(R^2\), distance, likelihood-ratio, and test statistics become derived transforms or diagnostics rather than separate analysis traditions. Accuracy is singled out as especially lossy because thresholding discards confidence, priors, and magnitude of near misses.

Approximate encoder and decoder models trained independently need not yield equal empirical scores; coherence of the fitted joint model is a substantive condition.

#### Effect, generalization, and calibration

Part 4 separates three axes that classical statistics and MVPA combine differently:

\[
\text{effect}\oplus\text{generalization}\oplus\text{calibration}.
\]

- **Effect:** what relation was estimated, at what moment order, scale, rank, and support?
- **Generalization:** across which runs, sessions, people, items, tasks, scanners, sites, or groups must it reproduce? This is encoded by \(\Gamma\).
- **Calibration:** what uncertainty or incompatibility with a null is attached through intervals, hierarchical posteriors, randomization/permutation tests, prevalence, or multiplicity control?

A p-value is calibration attached to an effect estimate; cross-validation specifies generalization. Neither is an alternative to the other.

#### Linear-Gaussian realization: relation operator and adjoint

Let common experimental space be \(\mathcal D\), participant-native neural space be \(\mathcal N_i\), and

\[
F_i:\mathcal D\to\mathcal N_i,
\qquad y_i=F_id+\epsilon_i.
\]

With \(d\sim\mathcal N(0,S_D)\) and \(\epsilon_i\sim\mathcal N(0,N_i)\), define the normalized relation

\[
R_i=N_i^{-1/2}F_iS_D^{1/2}.
\]

Its SVD pairs experimental and native neural modes:

\[
R_i=L_i\operatorname{diag}(s_{ik})V_i^\top.
\]

In whitened coordinates, the posterior-mean linear decoder is

\[
E[d_w\mid y_w]
=(I+R_i^\top R_i)^{-1}R_i^\top y_w,
\]

which applies modewise shrinkage \(s/(1+s^2)\). More generally, for design and neural inner products \(G_D\) and \(M_i\), the metric adjoint is

\[
F_i^*=G_D^{-1}F_i^\top M_i.
\]

The decoder is a regularized inverse built from the adjoint, not a second scientific object. This places the forward-pattern/backward-weight distinction inside the relation geometry itself.

#### Pullback representation and pushforward activation

The normalized relation produces two canonical Gram operators:

\[
G_i=R_i^\top R_i
=S_D^{1/2}F_i^\top N_i^{-1}F_iS_D^{1/2}
\]

in experimental space, and

\[
Q_i=R_iR_i^\top
=N_i^{-1/2}F_iS_DF_i^\top N_i^{-1/2}
\]

in native neural space. They share nonzero eigenvalues \(s_{ik}^2\). Part 4’s compact statement is:

\[
\text{representation}=R^\top R,
\qquad
\text{activation}=RR^\top.
\]

Representation describes experimental distinctions as expressed neurally; activation describes the participating neural directions. They are pullback and pushforward geometries of one morphism.

#### Relation spectrum unifies statistical and predictive summaries

For the Gaussian model,

\[
I(D;Y)
=\frac12\log\det(I+R^\top R)
=\frac12\sum_k\log(1+s_k^2).
\]

With canonical correlations \(\rho_k^2=s_k^2/(1+s_k^2)\),

\[
I(D;Y)=-\frac12\sum_k\log(1-\rho_k^2).
\]

Canonical \(R^2\), signal-to-noise roots, Wilks’ lambda, Gaussian information, rank-one \(t^2/F\), likelihood ratios, and linear prediction all summarize or calibrate this same spectrum in different ways. The proposed primary report is the modal spectrum

\[
\{s_k^2,\rho_k^2,\iota_k\},
\qquad \iota_k=\tfrac12\log(1+s_k^2),
\]

together with generalization domain and uncertainty. Exact finite-sample equivalence to a named \(t\) or \(F\) statistic depends on design, degrees of freedom, nuisance estimation, and covariance assumptions.

Cross-partition products

\[
\widehat G_i^{\mathrm{cv}}
=\frac1{R(R-1)}\sum_{r\ne s}
\operatorname{sym}(\widehat R_{ir}^\top\widehat R_{is})
\]

remove the positive self-product noise bias under independence and zero-mean estimation error. Such estimates may be indefinite; unbiased inference and PSD descriptive reconstruction must remain distinct.

#### Native voxel systems as gauge choices

Under an invertible native coordinate change

\[
y_i'=T_iy_i,
\quad F_i'=T_iF_i,
\quad N_i'=T_iN_iT_i^\top,
\]

the experimental pullback geometry is invariant:

\[
S_D^{1/2}F_i'^\top N_i'^{-1}F_i'S_D^{1/2}=G_i.
\]

Native activation coordinates are gauge-dependent; the induced experimental geometry is gauge-invariant, provided the coordinate change is invertible and the neural metric/noise transforms coherently. This does not cover information destroyed or added by noninvertible measurement transformations.

If \(R_i=U_iR_\star\) with subject-specific isometries, all participants share \(R_\star^\top R_\star\) and the same information spectrum even when the coordinatewise average \(N^{-1}\sum_iR_i\) cancels. The relation is shared although its embedding is not fixed-coordinate shared.

#### Population claims as a quotient hierarchy

The identifiable group object may be the equivalence class of native relation factors under metric-preserving transformations, represented by \(G_i=R_i^\top R_i\). Group inference should therefore model \(G_i\) before averaging voxel patterns or scalar accuracies.

Part 4 separates increasingly strong claims:

1. **Shared information magnitude:** similar \(\frac12\log\det(I+G_i)\).
2. **Shared representational geometry:** similar full \(G_i\).
3. **Shared transportable embedding:** biologically plausible transformations align relation factors.
4. **Shared fixed-coordinate topography:** patterns agree under anatomy-only correspondence.

Traditional individual decoding often addresses only the first; fixed-grid group analysis approaches the fourth. UEG is intended to fill the middle.

#### Bures geometry and alignment duality

For PSD geometries, Part 4 proposes Bures–Wasserstein distance

\[
d_B^2(G_i,G_j)
=\operatorname{tr}(G_i)+\operatorname{tr}(G_j)
-2\operatorname{tr}[(G_i^{1/2}G_jG_i^{1/2})^{1/2}].
\]

For suitable factors it equals the residual Frobenius discrepancy after optimal isometric/partial-isometric alignment. Thus direct comparison of invariant Gram geometries and Procrustes-style alignment of relation factors are dual views. The Bures barycenter

\[
G_\star=\arg\min_G\sum_i d_B^2(G,G_i)
\]

is proposed as consensus geometry modulo native neural rotations/embeddings.

This separates whether geometries are equivalent from whether their neural embeddings can be aligned by a plausible transport. Alignment becomes a secondary analysis of implementation after invariant consensus has been established.

Because raw crossvalidated geometry estimates may be indefinite, Bures modeling requires a latent PSD estimate or an inferential model that respects the distinction.

#### Do not scalarize before population modeling

RSA correlations and classification accuracies are many-to-one projections of \(G_i\). Equal scalars can conceal different residual geometries, dominant modes, spatial coherence, moment orders, nuisance exploitation, or simply rotated embeddings. The durable principle becomes:

> Model or aggregate the relation geometry before aggregating the scalar view generated from it.

#### Separate cortical base from functional fiber

For native cortical manifold \(\mathcal M_i\), local neural fiber \(\mathcal N_{i,x}\), and local relation

\[
F_i(x):\mathcal D\to\mathcal N_{i,x},
\]

define the invariant local geometry field

\[
G_i(x)
=S_D^{1/2}F_i(x)^\top N_i(x)^{-1}F_i(x)S_D^{1/2}.
\]

Two distinct correspondences are required:

- **Base transport \(\tau_i\):** where on the native cortex corresponds to a group cortical territory?
- **Fiber transport \(U_i\):** how do local functional axes map into a common neural embedding?

Geometry inference requires only base correspondence because every \(G_i(x)\) lives in common experimental space. With probabilistic anatomical correspondence \(\pi_i(dx\mid u)\), one can transport the low-dimensional field:

\[
\widetilde G_i(u)=\int G_i(x)\pi_i(dx\mid u).
\]

Raw BOLD and voxel patterns remain native. MNI is retained as reporting chart, visualization coordinate, or anatomical prior—not assumed to be the ontology of the neural effect.

#### Native-space local–global conservation

For participant-specific metric \(M_i\) and positive local operators satisfying

\[
\sum_\alpha M_{i\alpha}=M_i,
\]

define

\[
G_{i\alpha}
=S_D^{1/2}F_i^\top M_{i\alpha}F_iS_D^{1/2}.
\]

Then

\[
\sum_\alpha G_{i\alpha}=G_i^{\mathrm{global}}.
\]

The conservation law survives differing native dimensions, searchlight shapes, surface domains, and latent representations. A scalar searchlight is merely a query of the matrix-valued field \(\alpha\mapsto G_{i\alpha}\); nonlinear information views such as log determinants do not themselves add across the frame even when the matrices do.

#### Matrix-valued native-space searchlight and group model

The proposed replacement output is

\[
x\mapsto G_i(x),
\]

not \(x\mapsto\mathrm{accuracy}_i(x)\). Each local matrix retains modal strengths and angles, dimensionality, level/configuration decomposition, RSA projections, information spectrum, and residual geometry. Scalar hypothesis maps \(\operatorname{tr}[A G_i(x)]\) remain available on demand.

At group level, latent PSD geometry fields can be modeled using low-rank factors, Bures hierarchies, tangent-space regression, or participant relation factors. Outputs include group centre, geometry heterogeneity, information, covariate effects, prevalence, transportability, and fixed-coordinate coherence.

Part 4 distinguishes four group agreements: information magnitude, full geometry, transportability of embeddings, and fixed-coordinate topography. The combination of these identifies shared fixed topography, shared but rotated/transportable representation, abstract geometry without simple transport, true geometric difference, or scalar equivalence concealing divergent geometries.

#### Geometry versus implementation

Equal \(G_i\) establishes noise-normalized representational equivalence, not identical circuits or anatomy. Part 4 makes a ladder:

1. representational equivalence (gauge-invariant);
2. functional-embedding equivalence under a low-complexity transport;
3. anatomical implementation equivalence.

Claims must not jump between levels.

#### Nonlinear and higher-order extension

Mean relations live at first order. Covariance, nonlinear manifolds, and higher temporal moments require feature lifts such as

\[
\psi_2(y)=\operatorname{vec}(yy^\top),
\qquad \psi_m(y)=y^{\otimes m},
\]

and relation operators from experimental features into the lifted neural space. The direction-free density-ratio identity remains general, while the linear operator realization changes with the feature representation. The descriptors remain moment order, spatial scale, relation rank, population geometry, and generalization.

#### Experimental design remains central

Direction-free evidence does not weaken design. \(\mathcal D\) carries contrasts, factorial and continuous effects, stimulus features, nuisance subspaces, partition structure, and randomization/design weights. Conditional target information satisfies the same forward/reverse equality after conditioning on nuisance structure. At the linear level, the target is a conditioned or residualized relation subspace, with leakage avoided through appropriate cross-fitting.

#### Proposed canonical outputs

Participant-level output retains relation factor, invariant geometry, relation spectrum, information, generalization declaration, and uncertainty. Group output retains population geometry/barycenter, geometric heterogeneity, prevalence, transportability, and fixed-coordinate coherence. The key reporting improvement is to state how much information generalizes, which modes carry it, whether geometry is shared, and how its native embeddings relate.

#### Five grand theorems and paper

Part 4 proposes:

1. direction-free conditional evidence identity;
2. encoder/decoder adjoint or regularized-inverse duality;
3. activation–representation Gram duality;
4. gauge invariance of pullback relation geometry;
5. local–global conservation under a native measurement frame.

The proposed paper is *No Forward, No Reverse, No Common Grid: Gauge-Invariant Effect Geometry for Native-Space fMRI*. Its simulation crosses shared versus different geometries with fixed, rotated, or shifted native embeddings, and asks competing methods to distinguish shared information, shared geometry, transportable embedding, and fixed-coordinate topography.

### 2026-08-11 — Part 5: a geometry compiler for fMRI

#### Software thesis

Part 5 materializes the theory as a deliberately narrow **geometry compiler** positioned between effect estimation and scientific interpretation:

\[
\text{partitioned effect estimates}
\longrightarrow
\text{cross-generalized effect geometries}
\longrightarrow
\text{scientific views}.
\]

It should not own motion correction, HRF estimation, registration, generic classifiers, or visualization. It should define precise interfaces by which external systems provide or consume relation objects. The implementation maxim is:

> Relations enter. Geometry is the durable intermediate representation. Scalar maps are produced last as views.

#### One computational kernel: cross-generalized Gram

Let partitioned effects be \(B_{ir}:\mathcal D\to\mathcal N_i\). Factor a neural measurement metric as

\[
M_{i\alpha}=L_{i\alpha}^*L_{i\alpha},
\qquad L_{i\alpha}:\mathcal N_i\to\mathcal H_{i\alpha}.
\]

Then local/global measured geometry is

\[
G_{i\alpha}
=\sum_{r,s}\Gamma_{rs}
(L_{i\alpha}B_{ir})^*(L_{i\alpha}B_{is}).
\]

With a design-side projection \(C_\beta:\mathcal D\to\mathcal K_\beta\), define

\[
Z_{ir\alpha\beta}=L_{i\alpha}B_{ir}C_\beta^*,
\]

and the universal primitive

\[
\boxed{
G_{i\alpha\beta}
=\sum_{r,s}\Gamma_{rs}
Z_{ir\alpha\beta}^*Z_{is\alpha\beta}.
}
\]

Implementation stages are therefore neural measurement, design projection, cross-partition Gram formation, and graph-weighted accumulation. Full experimental geometry uses \(C=I\); one scalar contrast uses a one-dimensional \(\mathcal K\); large designs can use reduced or operator-valued geometry.

#### Minimal public vocabulary

The proposed user vocabulary has six concepts:

1. **Effects:** partitioned experimental-to-neural relations.
2. **Measure:** where/how native neural relations are observed.
3. **Generalize:** which independently estimated relations must agree.
4. **Geometry:** the resulting cross-generalized experimental-space relation.
5. **View:** a contrast, RDM, spectrum, information quantity, mode, or map derived from geometry.
6. **Population:** a model of geometry distribution across participants.

This vocabulary should serve notebooks, command-line tools, graphical workbenches, and backend execution alike.

#### Language-neutral core objects

**`Space`.** A named vector space with semantic identity, units, coordinate schema, basis identity, and optional domain—not merely an integer dimension. `DesignSpace` names experimental coordinates and transformations. `NeuralSpace` distinguishes domain (where measurements live) from basis (how signals are represented), preserving the Part 4 separation of cortical base and functional fiber.

**`EffectRelation<D,N>`.** A participant/partition-specific abstract linear map with source and target spaces, apply/adjoint operations, block access, estimator/independence metadata, and provenance. Storage may be dense, sparse, chunked, factored, latent, or lazy; public semantics must not depend on array orientation conventions.

**`NeuralMeasure<N,H>`.** A linear measurement \(L\) with apply/adjoint, support, location, scale, metric, and provenance. Point, weighted mean, diagonal searchlight, whitened local space, ROI basis, graph band, latent projection, and frozen learned decoder are all measures.

**`MeasurementFrame<N>`.** An indexed family \(\{L_\alpha\}\) with coverage, normalization, scale, and conservation declaration. A conservative frame satisfies

\[
\sum_\alpha L_\alpha^*L_\alpha=M_{\mathrm{global}},
\]

making local geometries sum to global geometry. Overlap and nonconservation must be explicit metadata rather than hidden behavior.

**`GeneralizationGraph`.** Sparse weighted edges over independently estimated relation nodes, with symmetry, normalization, semantic labels, and independence constraints. It replaces `folds = 5` with a declaration of what must reproduce. Tensor composition can express crossed run/session/item/participant generalization.

**`GeometryEstimate<D>`.** A symmetric design-space value/operator tagged by measurement, generalization, estimator type, pair count, uncertainty, and provenance. Types must distinguish an unbiased finite-sample `CrossGeometryEstimate`, which may be indefinite, from a modeled PSD `LatentGeometry` used for spectra and descriptive information.

**`GeometryField<Index,D>`.** The primary searchlight/multiscale result \(\alpha\mapsto G_\alpha\), logically shaped subject × measurement × design × design. Symmetry permits triangular storage; large design spaces can use projections, factors, matrix-free operators, or direct scalar queries without changing the interface.

**`DesignView`.** A pure, cheap interpretation of existing geometry: contrast energy, RDM, RSA model projection, eigen-spectrum, effective dimension, Gaussian information, coherence/configuration, paired modes, or local/global fraction. Switching views must not refit the analysis.

#### Immutable intermediate representation: `GeometryPlan`

Every frontend compiles to one content-addressed, serializable, inspectable, backend-independent plan containing effect source, canonical design space, measurement frame, generalization graph, design materialization, requested views, population model, calibration, and execution policy.

Familiar analyses become transparent presets over primitives:

- voxelwise geometry → point measurements;
- smoothed SPM-like query → weighted-mean measures;
- searchlight → local multidimensional frame;
- ROI → ROI frame;
- whole-brain MVPA → global or learned measure;
- RSA → geometry/model view;
- learned decoder → external learner returning a frozen measure.

No core method should be named `run_univariate_analysis`, `run_mvpa`, or `run_rsa`.

#### Materialization strategies

- **Full geometry:** retain all \(q\times q\) design relationships for many later views.
- **Projected geometry:** materialize a lower-dimensional contrast, factor, model, or PCA basis.
- **Direct scalar:** compute a one-dimensional hypothesis without building the full matrix.
- **Operator geometry:** expose apply, trace, leading-eigenpair, and approximate-logdet operations for very large design spaces.

Geometry is the abstraction; dense materialization is only one representation.

#### Compiler, not method dispatcher

The planner lowers every request to

\[
Z_{ir}=L_iB_{ir}C^*,
\qquad
G_i=\sum_{r,s}\Gamma_{rs}Z_{ir}^*Z_{is},
\]

then selects an algorithm from operator structure:

- diagonal measures → exact voxelwise quadratic filtering/collapse;
- rank-one measures → projected local effects;
- sparse measures → sparse multiplication/indexed gathering;
- low-rank measures → compute only \(LB\);
- global measures → small design Gram directly;
- translation-invariant windows → convolution/separable/FFT/integral algorithms;
- surface frames → sparse graph, geodesic, or heat-kernel filters;
- repeated views → cache geometry once.

The algebra should be closed under composition, adjoint, sum, scaling, direct sum, tensor product, restriction, factorization, and projection. New methods normally add an operator, graph, view, or population model—not a new analyzer subclass.

#### Native-space population interface

Participants contribute geometry fields in common design space despite different neural bases. A `GeometryCohort` combines design identity, subject fields, covariates, sampling uncertainty, and optional location correspondence. Population models may operate in the linear space of symmetric matrices, low-rank factors, tangent spaces, Bures geometry, or hierarchical shrinkage.

Two transport types must be unambiguously distinct:

- `LocationCorrespondence` transports geometry fields between native and reporting domains and may be probabilistic.
- `EmbeddingTransport` maps native neural bases and is required for signed patterns, decoder modes, or fixed-coordinate embedding comparisons.

Grouping geometry requires only the former; grouping basis-dependent neural patterns requires the latter. Implicit MNI resampling is forbidden.

The population rule is “model the relation/geometry object first; scalarize last.” Signed effects remain available through a separate `RelationField`; `effects.pattern(contrast)` is basis-dependent, whereas `geometry.energy(contrast)` is invariant under coherent basis changes.

#### Storage and provenance

An `EffectBundle` records manifests, spaces/domains, partitioned relations, measurement and generalization specifications, geometry fields, views, uncertainty, and provenance. Large arrays are chunked and lazy. Compact operator specifications replace large explicit matrices.

Stable identities derive from data, design, partitions, measures, generalization graph, software, precision, and estimator configuration. This supports exact caching, incremental recomputation, backend comparison, and auditing. Every number must answer:

> Which relation, measured how, generalized across what, and summarized in which way produced it?

#### Calibration stays orthogonal

Primary geometry types do not bake in p-values. A separate `Calibrator` declares target view/geometry, null transformation, exchangeability, resampling, and multiplicity. It may attach intervals, permutation/randomization distributions, posteriors, prevalence, or error control without changing the underlying effect/generalization object.

#### Adaptive analyses are train/freeze/evaluate

External learners receive training relations and return immutable `LearnedMeasure` objects recording fitted transform, training identities, regularization, and diagnostics. The geometry engine evaluates the frozen measure only on independent relations. Local covariance estimation, feature selection, classifiers, and whole-brain discriminants therefore obey the same leakage controls as fixed measures.

#### Exact decompositions, not reruns

Declared measure decompositions should expose total/coherent/configuration and DC/low/mesoscale/fine geometries with exact additive invariants. Users should query one geometry decomposition rather than run multiple pipelines and compare thresholded maps. This is both computationally cheaper and scientifically clearer.

#### Workbench and modules

A GUI becomes an editor/viewer over the same plan, with persistent Effects, Measure, Generalize, and Geometry panels and a continuously visible estimand. The proposed core modules are `spaces`, `relations`, `measures`, `generalization`, `geometry`, `population`, and `plan`; NIfTI/CIFTI/BIDS/SPM/zarr adapters remain outside the mathematical core.

#### Correctness laws

The highest-level software specification consists of algebraic/property tests:

1. adjoint law;
2. trace duality between design- and neural-side implementations;
3. gauge invariance under coherent coordinate/metric transformation;
4. measurement-frame conservation;
5. exactness of declared decompositions;
6. dense/streaming/sparse/convolution equivalence;
7. direct-query/full-materialization equivalence;
8. cross-null centering under independent zero-effect partitions;
9. rejection of implicit transport across incompatible neural bases.

These laws define backend substitutability more deeply than ordinary example tests.

#### Execution strategy

Because neural dimension usually dominates design dimension, the engine streams neural blocks, applies measurements while reading relations, and accumulates only small symmetric design matrices. It exploits low rank and sparsity, avoids \(V\times V\) matrices, parallelizes across participants/frame elements, uses deterministic reductions, and caches transformed blocks only when reused. Diagonal frames trigger the Part 1 filtering identity; graph frames use suitable approximations without changing public semantics.

#### Disciplined first release

The proposed first core includes partitioned linear relations, named design/native neural spaces, volume/surface domains, fixed point/mean/diagonal-searchlight/ROI/basis/global measures, run/session graphs, full/projected geometry, standard views, conservative frames, geometry-field storage, basic group modeling, explicit location correspondence, and matrix-free/out-of-core execution.

Learned measures, nonlinear lifts, covariance-only orders, Bures models, embedding transport, probabilistic atlases, distributed execution, and interactive visualization get extension points but need not be comprehensive initially. Preprocessing, GLM/HRF estimation, registration, generic classifiers/permutation engines, and publication graphics remain outside the core.

#### Smallest complete architecture

The implementation reduces to:

\[
\texttt{EffectRelation}
\to\texttt{NeuralMeasure}
\to\texttt{DesignView}
\to\texttt{GeneralizationGraph}
\to\texttt{CrossGram}
\to\texttt{GeometryView}
\to\texttt{PopulationModel}.
\]

The engine implements the algebra from which analyses are composed, not a catalogue of analyses.

### 2026-08-11 — Part 6: `effectgeom`, a clean-room successor concept

#### Design brief

The conversation now targets a new R package—not an rMVPA rewrite—that covers much of rMVPA’s scientific territory with a quarter or less of the code by replacing method- and engine-centred architecture with geometry algebra. Working name: **`effectgeom`**.

The package mental model is

\[
\texttt{effects}\xrightarrow[\Gamma]{W}\texttt{geometry}\xrightarrow{A}\texttt{view},
\]

followed optionally by geometry-level population modeling. \(W\) determines where/how brain features are measured, \(\Gamma\) declares which independent partitions must agree, and \(A\) or its compiled packed representation determines the scientific query.

The intended absences are architectural: no separate searchlight, ROI, or global engines; no RSA/model/result class families; no classifier registry; no classification-versus-regression branch. One numerical kernel should support fixed bilinear geometry queries.

#### Live rMVPA grounding at the time of Part 6/7

A read-only shallow checkout of `bbuchsbaum/rMVPA` at commit `3b12a8855b06f549e5772ccb9af070639a34752b` (2026-08-11) confirmed the comparison baseline:

- `DESCRIPTION` collated 82 R files and declared 23 imported packages.
- The README described `mvpa_dataset → mvpa_design → model_spec → engine → result`, with regional/searchlight engines and model families including RSA and cross-decoding.
- `workflow_api.R` exposed a `mode` choice among searchlight, regional, and global and branched to `run_searchlight`, `run_regional`, or `run_global`.
- `searchlight_engine.R` registered legacy, SWIFT, dual-LDA, and naive cross-decoding paths with eligibility/fallback logic.
- `data_roi_result.R` converted the newer `roi_result` back to a one-row tibble expected by legacy combiners.
- `mvpa_result.R`/`performance.R` maintained separate binary classification, multiclass, and regression result/performance behaviors.
- `dataset.R` included a Feature-RSA-specific multiple-feature validation in a general dataset constructor.
- `rsa_model.R` used RSA-specific design, model, fitting, merging, and output-schema machinery.

The clean-room conclusion is not that these features are individually ill-conceived; it is that algorithm and workflow names have become software ontology, forcing later fast paths and adapters to reconcile divergent execution and result systems.

#### Part 6 numerical normal form

For partition effects \(B_r\in\mathbb R^{p\times q}\), Part 6 first proposed one featurewise cross-generalized atom

\[
H_v=\sum_{r,s}\Gamma_{rs}b_{rv}b_{sv}^\top
\]

and an additive spatial frame \(W\in\mathbb R^{m\times p}\), yielding

\[
G_j^{\mathrm{total}}=\sum_vW_{jv}H_v.
\]

With symmetric vectorization \(\operatorname{svec}\) and \(h=q(q+1)/2\), stack atoms into \(Z_\Gamma(B)\in\mathbb R^{p\times h}\). The geometry field becomes

\[
\mathcal G=WZ_\Gamma(B).
\]

For \(k\) compiled linear queries \(C\in\mathbb R^{h\times k}\), output is

\[
Y=WZ_\Gamma(B)C.
\]

Spatial scope changes \(W\), generalization changes \(\Gamma\), and experimental interpretation changes \(C\); the engine remains fixed.

#### Built-in coherent/configuration compression

For frame row mass \(a_j=\sum_vW_{jv}\) and weighted partition relation \(\bar b_{jr}=a_j^{-1}\sum_vW_{jv}b_{rv}\), Part 6 defines

\[
G_j^{\mathrm{coherent}}
=a_j\sum_{r,s}\Gamma_{rs}\bar b_{jr}\bar b_{js}^\top,
\]

and

\[
G_j^{\mathrm{configuration}}
=G_j^{\mathrm{total}}-G_j^{\mathrm{coherent}}.
\]

The result simultaneously retains signed weighted effect, coherent regional energy, configuration energy, and total energy. Demeaning becomes a nondestructive component view, not preprocessing.

Local row normalization makes differently sized measurements comparable as averages; conservative column normalization makes local total geometries sum to global geometry. Both are explicit frame metadata.

#### RSA and classification treatment

Squared-distance RDMs are fixed linear transforms of packed geometry. Multiple-regression RSA compiles its distance transform and regression hat matrix into \(C_{\mathrm{RSA}}\), eliminating per-searchlight RDM regression and RSA-specific execution classes. Rank/correlation-normalized RSA remains a labeled nonlinear view.

The core does not reproduce a classifier zoo. A contrast query \(c^\top Gc\) is fixed-metric cross-generalized discriminability. Gaussian information or expected accuracy may be derived only under explicit noise-normalized, equal-covariance assumptions. Adaptive methods must return a frozen metric/measurement learned on independent training data; nonlinear methods must lift features upstream.

#### First object model proposed in Part 6

Part 6 initially proposed five durable concepts: `effect_set`, sparse `spatial_frame`, edge-table `pairing`, packed `geometry_field`, and compiled `geom_query`. Searchlight, ROI, voxel, and whole-brain calls differed only in frame constructor. One block-streamed engine accumulated feature atoms and weighted local relations, then obtained configuration by subtraction.

The first code budget was approximately 5,000 production R lines and a public API around effects, frames, pairings, geometry, views, simple transport, and population regression. Mathematical law tests—decomposition, conservation, query duality, RSA equivalence, streaming equivalence, invariance, null centering, and group/query commutation—were placed above example-based tests.

### 2026-08-11 — Part 7: the irreducible package, relation → frame → geometry

#### Revision of the primitive

Part 7 makes one decisive correction to Part 6: `effect_set` and mandatory beta estimates remain too concrete. The primitive becomes a lazy experimental–neural **relation** backed by raw time series, trialwise images, condition averages, precomputed effects, latent features, or an external estimator. The geometry kernel must not know which backing it received.

For partition \(r\), let response field be \(Y_r\in\mathbb R^{n_r\times p}\) and effect extractor be \(E_r\in\mathbb R^{q\times n_r}\). Then

\[
\boxed{B_r=E_rY_r}
\]

is the experimental–neural relation, but need never be materialized. The kernel requests feature blocks and applies \(E_r\) lazily.

For supplied design \(X_r\), target coefficient map \(C_r\), and observation precision \(P_r\), a GLS extractor is

\[
E_r=C_r(X_r^\top P_rX_r)^{-1}X_r^\top P_r,
\]

preferably constructed stably via whitening and QR rather than an explicit inverse. Precomputed effects are the special case \(Y_r=B_r,E_r=I\).

This yields the first foundational package law:

\[
\operatorname{geometry}(Y,E)
=\operatorname{geometry}(EY,I).
\]

OLS, GLS, ridge, FIR, LSA, LSS, condition averaging, temporal contrasts, and nuisance-adjusted effects differ only in the extractor once their design choices are fixed. The package accepts explicit extractors and offers at most a thin linear-model convenience; design/HRF packages can construct richer extractors externally.

#### Complete relation kernel

For neural measurement \(M_j\succeq0\) and pairing weights \(\Gamma\), the sole first-level estimand is

\[
\boxed{
G_j=\operatorname{sym}\sum_{r,s}\Gamma_{rs}B_rM_jB_s^\top
}
\]

or, from raw sources,

\[
G_j=\operatorname{sym}\sum_{r,s}\Gamma_{rs}
E_rY_rM_jY_s^\top E_s^\top.
\]

Raw versus beta changes \((Y,E)\); spatial scope changes \(M_j\); cross-run/session/task semantics change \(\Gamma\); contrast/RSA/MANOVA/information changes only the view.

The governing admission rule is **bilinear closure**: a core method must be expressible by this contraction followed by a declared query. Otherwise it must produce an extractor \(E\), measurement \(M\), transformed response \(Y\), or stay outside the core.

#### Additive-frame lowering and compiler normal form

The primary fast representation is \(M_j=D(w_j)\) with nonnegative sparse frame matrix \(W\). Feature atoms are

\[
H_v=\operatorname{sym}\sum_{r,s}\Gamma_{rs}b_{rv}b_{sv}^\top,
\]

packed into \(Z\), giving

\[
\mathcal G=WZ.
\]

Linear queries compile to \(C\), so the final normal form is

\[
\boxed{V=W\,\Phi_\Gamma(EY)\,C.}
\]

- \(E\) extracts experimental effects;
- \(\Gamma\) demands cross-partition agreement;
- \(W\) determines spatial measurement;
- \(C\) determines scientific interpretation.

Association can be chosen for efficiency: with few queries compute \(W(ZC)\) directly; full geometry is the same compiler with \(C=I_h\). These are optimization choices, not different scientific code paths.

At most one generic alternative lowering is proposed: factor frames \(M_j=L_j^\top L_j\) for fixed dense metrics/transforms. Additive diagonal and factor representations are storage/algebra choices, never named analysis engines.

#### Exact frame decomposition corrected for arbitrary mass

For nonnegative frame row \(w\) with mass \(a=\mathbf1^\top w\),

\[
D(w)=\frac{ww^\top}{a}
+\left[D(w)-\frac{ww^\top}{a}\right].
\]

Both terms are PSD. Pulling them through the relation gives total = coherent + configuration, with

\[
G^{\mathrm{coherent}}
=\frac1a\operatorname{sym}\sum_{r,s}\Gamma_{rs}u_ru_s^\top,
\qquad u_r=B_rw.
\]

A point frame has zero configuration. Searchlights/ROIs expose regional-mean and within-region configuration geometry simultaneously. Conservation remains a property of column-normalized frames.

#### Only two durable intermediate representations

Part 7 compresses the core to:

**`relation`.** Lazy response sources plus partition extractors, common experimental names, native feature/domain metadata, partition identities, and provenance. Its only computational protocol is

```r
relation_block(x, partition, features)
```

returning experimental dimension × requested features.

**`geometry`.** Signed local relation mean plus packed total and coherent geometries, frame index, experimental coordinates, pairing/metric/normalization metadata, diagnostics, and provenance. Configuration is derived lazily and exactly as total minus coherent.

Frames, pairings, and queries are small declarative values: sparse matrix plus index, edge table, and query matrix. Lightweight class tags may validate/print them, but there is no computational inheritance tree or `UseMethod()` cascade.

#### Refined API and boundaries

Canonical raw interface:

```r
rel <- relation(
  bold_runs,
  extract = lm_extractor(X_runs, effects = C, whiten = L_runs),
  domain = volume_domain(mask)
)
```

Precomputed effects omit the extractor and use identity. The package does not own events, HRFs, drift construction, motion correction, AR estimation, LSS design generation, registration, or file organization. It represents estimation without colonizing it and avoids a formula DSL.

Spatial scope is solely `at = voxels()`, `searchlights()`, `regions()`, or `whole_brain()`. Generalization is an edge table constructed by `cross_partitions()`, `cross_domains()`, or `pairing()`. Diagonal self-pairs are allowed only with explicit noise-bias labeling.

A contrast view returns signed, coherent, configuration, total, and coherence fraction together. RDM, compiled RSA, spectrum, effective rank, and conditional information views derive from the same geometry. Information is gated by metric/unit validity metadata.

Prediction no longer organizes the package. Learned transforms are frozen upstream and supplied as frames/metrics; nonlinear representations are supplied as transformed relations. No arbitrary “run any callback in every ROI” escape hatch is permitted because it defeats semantics and compiler optimization.

#### Plain R representation and error model

The proposed implementation uses validated lists, numeric matrices, `Matrix::dgCMatrix` frames, and data-frame edge tables. It avoids nested list columns, one-row tibble results, per-metric image objects, and neural-data-carrying fold objects.

Most errors are compile-time: extractor/source dimension mismatch, inconsistent experimental coordinates, frame/relation feature mismatch, zero-mass rows, invalid pairings/independence, and query dimension errors. Per-location status is a compact diagnostic table (`valid`, `support_size`, `mass`, `rank`, `reason`), not thousands of nested error objects.

Domains provide feature identity/coordinates, frame compilation, and rendering. They never encode analysis-specific constraints. A one-feature relation is valid; a frame or query requiring more features rejects it locally.

#### Group readiness without group-platform sprawl

Packed native geometry field \(H_i\) can be transported by sparse location correspondence \(P_i\):

\[
H_i^\star=P_iH_i.
\]

Raw patterns stay native. Day-one code need only preserve experimental-basis identity, packed coordinates, locations, uncertainty, and provenance; a later or separate group package can fit matrix linear/hierarchical models. Correct output makes the core group-ready without embedding another workflow platform.

#### Final API/code budget after Part 7

Proposed public core: `relation`, `effect_extractor`, `lm_extractor`; frame constructors; pairing constructors; `geometry`; and view functions (`contrast`, `rdm`, `rsa`, `spectrum`, `information`, `as_neuro`). Low-level extensions return algebraic values rather than registering plugins.

The source target is fewer than 3,000 executable R lines, with eight small modules (`relation`, `extractor`, `frame`, `pairing`, `kernel`, `geometry`, `views`, optional `neuroim2`) and only `Matrix` mandatory. Tests should exceed production code.

#### New foundational laws and prohibitions

Part 7 retains Part 6’s algebraic tests and adds raw/effect equivalence as the first law. Its constitutional prohibitions reject modes, data-bearing model specs, arbitrary callbacks, registries, result classes per method, destructive mean removal, fold objects carrying data, implicit resampling, mandatory betas, premature maps, workflow-platform features, and any new feature lacking algebraic reuse.

The irreducible architecture is:

\[
\boxed{
\text{Relation IR } B_r=E_rY_r
\quad\longrightarrow\quad
\text{Geometry IR }G_j
\quad\longrightarrow\quad
\text{late views}.
}
\]

### 2026-08-12 — Part 8: execution without inherited futures

#### Question

The mathematical and conceptual compression is not enough. rMVPA’s execution model is based on `future`/`furrr`, serial ROI preparation, batching heuristics, shared-memory exceptions, global plan state, and result/error reconciliation. The new package should improve execution just as radically as it improves the estimand.

#### Live execution audit

At current rMVPA commit `3b12a8855b06f549e5772ccb9af070639a34752b`:

- `future`, `future.apply`, and `furrr` are mandatory dependencies.
- `mvpa_iterate()` extracts ROI data into main-process batches before `future_pmap()`, with different memory/chunk heuristics for searchlight and regional modes.
- Worker model specs are manually stripped of the dataset to reduce serialization.
- Runtime behavior depends on `future::nbrOfWorkers()` and therefore the active global future plan.
- CLI and custom workflows temporarily mutate and restore that global plan.
- An experimental `shard` backend duplicates dataset-specific extraction and worker logic to avoid serial ROI extraction and serialization.
- Progress uses worker closures/`progressr`; remote shared-memory execution intentionally loses the local progress closure.
- Runtime errors become per-ROI/per-batch tibble rows, and partial failure handling is intertwined with result combination.
- The relevant iterator, shard, engine, fast-path, global, custom, CLI, and utility files together exceed 8,000 physical R lines.

This is not simply a poor choice of parallel library. It is the consequence of treating every ROI/searchlight fit as a job.

#### Central execution correction

In `effectgeom`, searchlights are rows of sparse frame \(W\), never tasks. The vectorized compiler

\[
V=W\Phi_\Gamma(EY)C
\]

changes the work graph:

1. stream a canonical neural-feature block;
2. lazily extract relation blocks \(B_{r,I}=E_rY_r[:,I]\);
3. form featurewise cross-generalized atoms or direct query values once;
4. contract the block through every spatial measurement using \(W[:,I]\);
5. reduce blocks in canonical order.

The package therefore removes repeated searchlight extraction/fitting, model serialization, prediction-table construction, per-searchlight RSA, and most scheduler overhead before adding any parallelism.

#### Explicit compute policy

Scientific semantics and execution resources stay separate. A small immutable `compute_policy` records workers, memory budget, feature-block size, process mode, optional checkpoint directory, progress callback, and threads per worker. It does not mutate session state or alter the scientific plan.

The default is deterministic one-process streaming. `workers > 1` is opt-in and valid only when the source access pattern can support it without silently duplicating the full dataset.

#### Compiled execution stages

**Compile.** Validate bases/dimensions, canonicalize pair edges, compile frame/query matrices, select direct-query versus full-geometry width, estimate memory/I/O/work, choose a feature block size, and inspect source capabilities.

**Feature task.** A pure task owns one disjoint feature interval and returns relation blocks plus packed atoms/direct query values and a small receipt. It receives no searchlight, ROI, model spec, output map, or arbitrary function.

**Coordinator reduction.** Completed blocks are consumed in canonical feature order. Sparse frame multiplication accumulates total geometry/direct views and local relation marginals. Workers never write shared scientific output and never return full-map partials.

**Certification.** Construct `effect_geometry` only after exact feature coverage, dimensions, provenance, and applicable algebraic laws pass. Partial fields do not masquerade as complete results.

#### Parallelism hierarchy

1. Parallelize participants first; they are genuinely independent and native-space separated.
2. If necessary, parallelize disjoint feature blocks within one participant.
3. Let compiled sparse/dense numerical kernels do local arithmetic.
4. Never combine multiple R processes with unbounded multithreaded BLAS; default to one numerical thread per worker.
5. Never schedule searchlight centres.

Parallelism may prove unnecessary for many analyses because the algebra eliminates the old redundant workload.

#### Source-aware locality

Each relation source declares whether it is seekable, reopenable, concurrent-read-safe, descriptor-serializable, or in memory.

- Coordinator-only/in-memory sources use bounded local streaming rather than copying the whole dataset to every process.
- Reopenable file-backed sources permit workers to receive descriptors plus feature indices and read locally.
- Future shared/memory-mapped sources implement the same protocol rather than creating dataset-specific execution engines.

There is no automatic shared-memory attempt followed by silent fallback.

#### Backend boundary

The core imports neither `future` nor `furrr`. Version 0.1 ships the sequential executor. If benchmarks justify it, an owned local-process executor may use base R `parallel` for worker-readable sources. The call creates, initializes, and closes only its own workers.

An internal pure task entry point permits later `future`, `mirai`, Slurm, or remote adapters. Such adapters only schedule fixed task manifests and return fixed task results; they cannot change scientific semantics or reduction.

#### Memory, determinism, and RNG

The compiler derives block size from a declared memory budget, accounting for fixed frame/output/local-relation storage plus active response/relation/atom blocks and buffers. It fails early with remedies if the budget is impossible. There is no analysis-type-dependent ROI batch knob and no reliance on `gc()` as scheduling.

Core geometry uses no RNG. Random frames, pairings, permutations, or bootstraps are realized upstream and become immutable, hashed inputs. Task completion order never determines reduction order. Worker count, block size, and backend must agree within declared numerical tolerance; a fixed ordered reduction supports reproducibility.

#### Failure and progress semantics

Compile-time structural problems fail before workers start. A runtime task failure is identified by task, stage, feature interval, source, condition, worker, and timing. Default response: stop new dispatch, close only owned resources, preserve explicitly requested checkpoints, and return no ordinary geometry.

No silent backend fallback, arbitrary retries, partial geometry views, or per-location error tibbles are allowed.

Workers return structured receipts rather than progress closures. The coordinator emits stable events based on processed features/frame nonzeros or planned work. Completed results retain one execution receipt with plan hash, workers/threads, source access, memory/block choices, bytes, task counts, timing, software/BLAS identity, and completion state.

#### Checkpointing

Later checkpointing writes atomic content-addressed feature-task shards. Resume validates plan/kernel/source identities and reruns the same canonical reduction. Checkpoints are operational artifacts, not a third scientific IR, and are deleted on success unless explicitly retained.

#### Revised implementation order

1. **E0:** deterministic sequential, query-aware feature streaming, memory planning, and execution receipt.
2. **E1:** pure task protocol, source capabilities, arbitrary-completion-order tests.
3. **E2:** owned local processes only if realistic benchmarks show benefit.
4. **E3:** atomic resume and optional external executor adapters.

The execution specification was initially saved as
`effectgeom-execution-design.md`; after the package naming decision it lives at
`effectagram-execution-design.md`.

### 2026-08-12 — Part 9: retain shard, remove the shard analysis path

#### Correction

Part 8 was too dismissive of shared memory. rMVPA's `shard` integration has a
real and important systems benefit: several R workers can read one large
immutable dataset without receiving one serialized copy each. That mechanism
was developed specifically around R's process-memory limitations and should not
be thrown away merely because rMVPA wrapped it in a parallel model-specific
execution branch.

The revised rule is:

> Do not reject `shard`; reject allowing shared memory to define a second
> scientific engine.

#### Live shard audit

The current `shard` 0.2.0 source at commit
`233c71186ebac0e2e98eba70f10671450fc5e1be` provides:

- immutable shared inputs through POSIX shared memory or memory mapping;
- lightweight descriptors that workers can attach to;
- explicit shared output buffers;
- persistent supervised PSOCK workers;
- recycling under RSS drift and worker failure;
- deterministic chunk-ordered reductions;
- copy/materialization and aggregate memory diagnostics;
- bounded table/result sinks and cleanup APIs.

This is closely matched to `effectgeom`'s remaining operational problem after
the algebra removes repeated searchlight fitting. It should be reused rather
than reimplemented with a thin, less capable base-`parallel` layer.

#### Architectural placement

`effectgeom` retains sole ownership of:

```text
relation + frame + pairing + query
feature-block task definition
cross-Gram kernel
canonical scientific reduction
effect_geometry result
```

`shard` may own:

```text
shared segment publication and attachment
worker pool and supervision
coarse task dispatch
worker recycling
bounded shared buffers
memory and copy telemetry
```

There is no `shard_model_spec`, `run_geometry.shard`, dataset-specific shared
extractor, or shard-specific result. The same pure feature-block task used by
the sequential compiler is the only worker kernel.

#### Two useful sharing points

For relation \(B_r=E_rY_r\), the adapter may explicitly stage either side of
the extraction boundary.

**Shared response.** Publish \(Y_r\) once. Each worker reads its requested
feature block and computes \(E_rY_r[:,I]\). This preserves lazy raw-data
extraction and avoids copying the largest source.

**Shared relation.** Compute \(B_r\) once and publish the \(q\times p\) relation.
When \(q\ll n_r\), this is much smaller than the raw response and can be reused
across frames, queries, and geometry plans. It is a transient compiler cache,
not a mandatory beta-image pipeline stage.

The first adapter should make staging explicit (`none`, `response`, or
`relation`). An automatic choice may be added only after benchmarks establish
a reliable cost model using source size, extraction cost, \(n_r/q\), and reuse.

#### Coarse tasks and bounded results

`shard_map()` should receive a small number of coarse feature blocks, never one
job per searchlight. Workers return canonically indexed relation/atom/query
blocks. `effectgeom` reduces them in feature order through \(W\).

If gathered blocks exceed the memory plan, workers may write disjoint atom or
query rows to a `shard::buffer()` and the coordinator may scan that buffer in
order. Workers do not write the final overlapping spatial geometry. Scientific
reduction remains in `effectgeom`; `shard_reduce()` is used only where its
associative combine and numerical order exactly match the required law.

#### Lifecycle is part of correctness

Shared memory introduces a stronger resource contract:

1. the coordinator alone owns and unlinks a segment;
2. workers attach read-only and only detach;
3. the owner outlives every worker attachment;
4. every dispatch generation has a distinct identity;
5. recycled or persistent workers reopen current descriptors and never reuse a
   stale cached handle;
6. success, error, interrupt, timeout, and retry exhaustion use one idempotent
   cleanup path;
7. a run with uncertain cleanup is failed, not successful with a warning.

One-worker execution uses the direct local loop by default. Shared staging is
only paid for in that case when the user explicitly wants a reusable staged
relation.

#### Evidence gate

The adapter must prove more than faster wall time. Benchmarks report:

- aggregate private RSS across coordinator and workers;
- shared mapped bytes separately;
- serialization volume;
- hidden materialization/copy events;
- staging cost and reuse break-even point;
- worker recycling and end-of-run memory return;
- cleanup after success, failure, and interrupt;
- numerical agreement with dense, streamed, and ordinary-source execution.

This preserves the main Part 8 insight—parallelism operates over feature blocks,
not analysis objects—while retaining the strongest part of the infrastructure
already developed for rMVPA.

### 2026-08-12 — Part 10: execution review closes the specification gaps

#### Verdict received

An independent architectural review approved the execution thesis but correctly
rejected the description “fully specified.” It confirmed the rMVPA diagnosis
and identified nine gaps: contraction temporaries, overbroad searchlight scope,
unbounded out-of-order buffering, inconsistent determinism language, ambiguous
direct-query results, contradictory version-0.1 parallel promises, mutable
callbacks inside `compute_policy`, incomplete checkpoint integrity, and missing
integration into the public package contract. It also identified an arbitrary-
orientation defect in undirected left/right pairing marginals.

#### Scope correction

The authoritative slogan is now:

> **Additive-frame searchlights with fixed bilinear queries are rows of a sparse
> matrix, not jobs submitted to workers.**

Locally trained classifiers, locally estimated covariance, nonlinear
normalizations, and generic factor frames do not generally collapse to
(W\Phi_\Gamma(EY)C). They may later receive separate compiled lowerings, but
cannot inherit the additive theorem by name.

#### Contraction-memory correction

Feature blocking alone does not bound the dense intermediate from
(W[:,I]Z_I). For large measurement count (m) and packed width
(h=q(q+1)/2), one multiplication can allocate gigabytes, and replacing the
accumulator can allocate another full object.

The canonical contraction is now tiled over feature block (I_t), measurement
rows (J_a), and packed coordinates (K_b):

\[
G[J_a,K_b]\mathrel{+}=W[J_a,I_t]Z_{I_t}[:,K_b].
\]

Full geometry may be in-memory or block-backed behind the same complete
`effect_geometry` contract. The memory planner includes result storage,
contraction temporaries/R copies, serialization overlap, bounded reordering,
checkpoint buffers, shared mappings, and per-task state. Its claim is
conservative planning validated against measured peak memory—not exact portable
prediction of R allocation.

#### Bounded scheduling and numerical contract

Dynamic scheduling now has hard `max_inflight` and `max_reorder_bytes` limits
derived from the memory budget. Dispatch pauses under backpressure; completed
tasks are released immediately after canonical reduction; checkpoints may spill
validated blocks but cannot authorize unbounded dispatch.

Version 0.1 promises scheduling invariance and tolerance-qualified agreement
across block/tile choices and platforms. It does not promise bitwise equality.
A later stronger mode may separate I/O blocks from fixed reduction microblocks
and use a canonical pairwise/compensated tree.

#### Honest public results

Two calls now have distinct contracts:

```r
geometry(..., materialize = "full", storage = "memory" | "block")
  -> effect_geometry

evaluate_geometry(..., query = geom_contrast(...))
  -> effect_view
```

Storage changes physical representation, not completeness. A direct query is
cheaper but cannot later answer an uncomputed query and never masquerades as a
complete geometry.

Version 0.1 accepts `workers = 1` only; any other value fails before source
access. `compute_policy` contains only immutable numeric execution choices.
Reporter callbacks and checkpoint destinations are separate, non-semantic
attachments with separate receipt fields.

#### Pairing correction

For an undirected edge stored once, left/right labels are arbitrary. The signed
summary is therefore the orientation-invariant endpoint marginal

\[
\bar B_j^{\mathrm{endpoint}}
=\frac12\sum_e\omega_e
\frac{(B_{l(e)}+B_{r(e)})w_j}{a_j}.
\]

Only genuinely directed pairings retain role-specific left/right marginals.

#### Failure, checkpoint, and verification corrections

Raised execution conditions carry partial receipts and cleanup state. Reporter
failure is isolated from scientific execution. Checkpoints require disk-space
preflight, restrictive permissions, schema/checksum metadata, same-filesystem
atomic rename, strong source revisions, and exact task/tile/precision/kernel/
reduction-plan identity for resume.

Verification is organized as contract, dense-oracle, metamorphic, adversarial,
regression, and performance test families. An in-process scheduling simulator
precedes all worker processes.

#### Mote implementation epics

The review was converted into three dependency-ordered epics with fifteen
children in the repository's new `.mote` store:

1. `bd-01KZTVQ3MATVABJA9ME77JS9P4` — specify effectgeom 0.1 semantic contracts;
2. `bd-01KZTVQDKFTD0NTAEKBMQ477XP` — build the memory-bounded sequential compiler;
3. `bd-01KZTVQDZ36XCT46WMXWHPW3N9` — harden scheduling, recovery, and optional executors.

The first ready implementation items are theorem scope
(`bd-01KZTVR1WZCXFKR3KESDXF91B8`), undirected marginal semantics
(`bd-01KZTVR2SNZFHAF853NZR0DG3N`), and the sequential-only version contract
(`bd-01KZTVR3H7K4F3911WQ2VZA1PH`). Epics are blocked by their children, and
later runtime work is blocked by the semantic and sequential foundations.

### 2026-08-12 — Part 11: `effectagram` becomes the canonical name

The package and project name is now **`effectagram`**, pronounced
“effect-a-gram.” An effectagram is the cross-generalized geometry of
experimental effects across a declared spatial frame. The name retains the
mathematical echo of a Gram matrix while also naming an interpretable record of
an effect, which better fits the intended system than the provisional
`effectgeom` name.

The rename is package-facing, not a replacement of the mathematical
vocabulary. Durable scientific objects remain `effect_relation`,
`effect_geometry`, and `effect_view`; “effect geometry” remains the underlying
concept. Earlier `effectgeom` references in this ledger are preserved as
historical snapshots. Current package metadata, diagnostics, tests, benchmarks,
vision, mission, and design artifacts use `effectagram`.

Preliminary checks on 2026-08-12 found `effectagram` available on CRAN,
Bioconductor, and GitHub and found no obvious exact-name scientific-software
collision. This establishes a practical working identity, not trademark
clearance.

### 2026-08-13 — Part 12: from evidence pairing to sampling-covariance transport

#### The statistical completion

The evidence-pairing theorem remains the irreducible point-evidence law:

\[
\mathscr E_{LR}(H,K)
=
\operatorname{tr}(H^\top B_LKB_R^\top).
\]

The new step is to distinguish an estimated pairing from its sampling law. For
a family of evidence queries \(a=(H_a,K_a)\), define

\[
\boxed{
\mathcal V_{ab}
=
\operatorname{Cov}
\left[
\widehat{\mathscr E}(H_a,K_a),
\widehat{\mathscr E}(H_b,K_b)
\right].
}
\]

This does not replace or enlarge the irreducible observable. It is the second
moment of the same observable under relation-estimation noise. The architecture
therefore becomes

\[
\boxed{
\text{relation fit}
\to
\text{evidence plan}
\to
\text{estimated evidence}
\to
\text{sampling-covariance plan}
\to
\text{calibrated view}.
}
\]

Effect, generalization, and calibration remain distinct. An interval does not
define an effect; a partition pairing does not define an interval.

#### Published law recovered

Diedrichsen, Provost, and Zareamoghaddam (2016), *On the distribution of
cross-validated Mahalanobis distances* (arXiv:1607.01371), derive the complete
covariance of the vector of crossvalidated distances. Their Eq. 13 is

\[
V
=
\left[
4\frac{\Delta\circ\Xi}{M}
+
2\frac{\Xi\circ\Xi}{M(M-1)}
\right]
\frac{\operatorname{tr}(\Sigma_R^2)}{P^2}.
\]

The diagonal gives the sampling variance of each distance; the off-diagonal
gives covariance between distances. Section 3.4 explicitly treats averaging
across cross-validation folds whose inner products reuse partitions. The first
term is signal-dependent and the second is signal-independent.

Their Eq. 10 estimator initially appeared different from effectagram's
all-unordered-pairs average. Algebra and numerical checks showed exact
equivalence:

\[
\frac1M\sum_m
\widehat\delta_m^\top
\left(\frac1{M-1}\sum_{n\ne m}\widehat\delta_n\right)
=
\frac1{\binom M2}\sum_{m<n}
\widehat\delta_m^\top\widehat\delta_n.
\]

The maximum numerical discrepancy in the verification was
\(2.3\times10^{-17}\). An independently derived diagonal formula agreed with
Eq. 13 to \(1.1\times10^{-16}\) across 24 checked cells and with empirical
simulation standard deviations to within 1.2 percent. This is a verification
and rediscovery of published theory, not a novelty claim.

#### The partition-edge correction

The pairing contract already stated that a linear edge reducer does not imply
independent edges. Part 12 supplies the sampling consequence:

> Independent partition estimates can form unbiased cross-products, while the
> cross-products themselves remain dependent because several edges reuse each
> partition.

Thus the rows returned by `cross_partitions()` are contributions to the
generalization estimator, not replicates for
`sd(edge_values) / sqrt(number_of_edges)`.

One scaling statement required correction. Comparing the signal term to a
known marginal edge variance divided by the nominal pair count yields the
familiar asymptotic factor \(\sqrt{M-1}\). That is not the exact factor for the
usual sample-SD-across-edges recipe. For pair fluctuations
\(y_{ij}=a_i+a_j\),

\[
\operatorname{Var}(\bar y)=\frac{4\sigma_a^2}{M},
\qquad
\mathbb E\left[
\frac{s_y^2}{\binom M2}
\right]
=
\frac{4\sigma_a^2}{M(M+1)}.
\]

The true variance divided by the expected naive squared standard error is
therefore \(M+1\), giving a standard-error factor \(\sqrt{M+1}\). Later prose
must name the naive estimator before naming its factor.

#### What the package already contains

On the `lm_relation_fit()` path, public accessors expose the ingredients of the
fixed-metric, equal-partition specialization:

```text
effect_covariance()        -> effect-coordinate covariance
residual_block()           -> whitened residual information
residual_df()              -> residual degrees of freedom
residual_pair_statistics() -> support-local residual cross-products
```

A public-accessor calculation recovered mean analytic standard errors of
0.00279, 0.02059, and 0.04094 for simulations whose empirical standard
deviations were 0.00269, 0.01909, and 0.04153. The corresponding nominal 95
percent coverages were 0.987, 0.967, and 0.940. These checks establish the
special-case ingredients and recipe; they do not establish a general exported
calibration API.

A pure `relation()` made from precomputed betas correctly refuses all four
error-model accessors. That path supports point evidence but cannot recreate
effect covariance, residual covariance, or residual degrees of freedom that
were discarded upstream. The prior error text made this look like a malformed
object rather than a capability boundary. The design correction is to preserve
`relation` as the pure algebraic object and make `relation_fit` or an explicit
external error channel the source of sampling capabilities.

Subject-level resampling may estimate population uncertainty when subjects are
available. It does not recover missing within-participant uncertainty for a
single precomputed-beta relation.

#### Boundaries retained

The verified four-line recipe is not universal. Its simple scalar form assumes
a common fixed metric, equal partition covariance, the declared separable GLM
construction, and its stated distance normalization. The complete theory
retains the metric-transformed residual trace term and provides more general
expressions for unequal partitions.

Learned local metrics \(K_{x,e}\) are random and may vary by location and
evaluation edge. Point evidence remains an evidence pairing with a
provenance-frozen, on-demand operator recipe, but Eq. 13 cannot be transferred
unchanged while claiming calibrated intervals. Metric-estimation uncertainty
and dependence must be propagated or the limitation must be explicit. The
existing certification report is therefore correct to export no interval or
LD-t for learned metrics.

Likewise, one local covariance \(V_x\) concerns distances or queries at a
location. It does not provide covariance between overlapping searchlights and
does not by itself solve spatial multiplicity.

#### Execution consequence

For \(D\) distances, a dense covariance has \(D^2\) entries. The same
query-first principle now applies to uncertainty: preserve an exact sampling
operator and compute `diagonal`, selected entries, `apply(a)`, or
\(a^\top Va\) directly. Full materialization is explicit and size-preflighted.
Structured Hadamard products, low-rank query transport, symmetry, sparsity, and
matrix-free action are preferred exact lowerings. Approximation remains a
separately named estimator with an error contract.

#### Normative disposition

Part 12 produced `evidence-sampling-v1`, a statistical contract adjacent to and
dependent on `evidence-pairing-v1`. The evidence-pairing contract receives only
the abstract boundary: it defines point evidence, while the statistical
contract defines covariance among estimated evidence values. The pure
effect-form laws remain free of Gaussian or residual-model assumptions.

The first implementation milestone is fixed-metric, equal-partition analytic
covariance with exact linear transport, explanatory refusal on a bare
precomputed relation, independent product-oracle tests, generative coverage,
and query-first scale gates. Heterogeneous partitions, learned-metric
uncertainty, cross-location covariance, LD-t, and population inference remain
separately gated extensions.

#### Mote implementation program

The statistical completion was converted into one dependency-ordered program
epic and three phases:

1. `bd-01KZYJH5MFBB2ZSM4TFCSF1SN8` — sampling-covariance transport and
   calibrated evidence;
2. `bd-01KZYJHYZKMPW86N8XTGW2KX97` — freeze sampling semantics and capability
   gates;
3. `bd-01KZYJHZP2JNDJXZ7JZ20G3SA2` — implement fixed-metric analytic
   covariance transport; and
4. `bd-01KZYJJ04RRCN63A7HE33DKKXV` — validate and teach calibrated evidence.

The program depends on the completed statistical and brain-scale evidence
pairing epic. The existing partition-dependence and precomputed-beta P1s are
children of the validation phase rather than duplicated findings. Unequal
partitions and learned-metric uncertainty are tracked as a non-blocking
research branch so the fixed-metric slice can close on its own evidence.

#### Part 12 implementation checkpoint: exact, calibrated, and scale-qualified

The first statistical vertical slice is now implemented through one compiled
sampling plan rather than a parallel RDM inference engine. The product path is

```text
lm_relation_fit
  -> fixed-metric geometry plan
  -> rdm point estimate
  -> factorized RDM sampling covariance
  -> exact covariance query
```

The public surface is deliberately two functions:

- `rdm_sampling_covariance()` constructs the local sampling law and requires an
  explicit `target = "plugin"` or `target = "null"`;
- `sampling_covariance()` reads its diagonal, selected entries, action,
  quadratic form, linear transport, or explicitly materialized dense form.

The covariance object stores row factors for the signal and design terms, not
the dense \(D\times D\) covariance and not one fourth-order tensor. Exact
queries use the identity

\[
v^\top[(LL^\top)\circ(RR^\top)]v
=
\|L^\top\operatorname{diag}(v)R\|_F^2.
\]

The corresponding action and transport use the same compact weighted
cross-product. This replaces a literal row-tensor rank-product loop with BLAS
contractions while preserving the declared covariance exactly. Dense
materialization is explicit and charges both output and working-row storage
before allocation.

The support-local residual path also changed during scale review. A one-node
query now reads residual blocks only for that support. Frame-wide sparse pair
statistics are an explicit batch strategy, not a hidden default. The same
searchlight support graph remains available to share residual cross-products
when many overlapping nodes are requested.

The statistical gate used 10,000 independent repetitions under the declared
fixed-metric, equal-partition matrix-normal model. Relative Frobenius error for
the complete covariance was 1.16 percent under nonzero signal and 1.96 percent
under the null. Coordinatewise analytic coverage ranged from 94.75 to 95.29
percent; a declared linear transport had variance ratio 0.9952 and coverage
95.06 percent. The partition-mean plug-in policy was mildly conservative at
97.70 to 97.99 percent and is reported as such. Treating dependent partition
edges as replicates reduced signal coverage to 72.84 to 73.97 percent.

The scale gate used a 60,000-feature volumetric domain and 60,000 searchlights
with mean support 112.3, together with \(q=120\) effects and \(D=7{,}140\)
distance coordinates. The spatial frame compiled in 6.94 seconds and occupied
193.2 MB. The exact factorized covariance occupied 16.6 MB; diagonal,
100-entry, quadratic-form, and eight-output transport queries each completed
in 0.18--0.30 seconds on the recorded machine. The corresponding dense
covariance payload is 407.8 MB and was correctly refused under a 64 MiB
budget. These measurements qualify the named fixtures and numerical contract;
they are not universal runtime promises.

This checkpoint realizes the mission's performance clause in its intended
form: mathematical elimination, factorization, support locality, query fusion,
and bounded action first; approximation only as a separately named future
estimator. It also preserves the vision's boundary. Precomputed beta relations
retain point evidence but refuse within-participant covariance because their
error channel is absent. Learned metrics remain point-estimate capable but
uncalibrated until metric uncertainty is propagated. Local covariance remains
local and is not relabelled as a spatial random field.

## Living synthesis after Part 1 (historical snapshot)

### Essence in one paragraph

A searchlight forms a local observation-space kernel by spatially pooling voxelwise rank-one Gram atoms. If the scientific score is linear in that kernel, the score can be pushed through the pooling operation: compute a quadratic evidence value at every voxel and spatially filter the resulting map. This is true even for fixed-fold crossvalidated bilinear statistics. A searchlight only needs genuinely local multifeature computation when the feature metric couples voxels, the procedure adapts to local data, or the kernel readout is nonlinear. General local feature operators are dual to scans of the evidence graph \(Z^\top AZ\), while coverage-normalized diagonal operators partition the global kernel exactly. The framework therefore supplies an algebraic taxonomy, computational shortcuts, a conservation principle, and a path to multiscale graph-based mapping.

### Argument dependency chain

1. Voxel columns generate rank-one observation-space Gram atoms.
2. A diagonal searchlight window linearly pools those atoms.
3. A trace readout is linear in the pooled kernel.
4. Therefore readout and spatial pooling commute.
5. Fixed cross-validation can be encoded in the trace operator \(A\).
6. Failures of collapse locate the true feature interaction: dense metrics, local adaptation, or nonlinear readouts.
7. Replacing diagonal windows by general \(M_s\) exposes an evidence graph with node and edge terms.
8. Normalizing overlapping windows to a partition of unity conserves the global kernel and all linear scores.
9. A family of smoothing kernels turns radius into scale-space or graph diffusion time.

### Boundary map

| Procedure | Collapse status asserted in Part 1 | Decisive condition |
|---|---|---|
| Squared-Euclidean RSA / fixed RSA regression | Exact | Fixed design and linear, unnormalized readout |
| Fixed-fold crossvalidated bilinear contrast | Exact | Fold contrasts compile into fixed \(A\) |
| Fixed-metric cvMANOVA numerator/effect | Exact or partial | Precision is global or fixed diagonal |
| Standard locally whitened cvMANOVA | Generally no | Dense precision is re-estimated locally |
| Linear-kernel HSIC or MMD | Exact | Fixed target/group kernel |
| Pearson/cosine/rank RSA | Generally no | Local normalization or ranking is nonlinear |
| Ridge/SVM/logistic fitting and accuracy | Generally no | Inversion, optimization, or thresholding |
| General dense \(M_s\) with trace readout | Linear graph scan, not node-only smoothing | Off-diagonal metric reads evidence edges |

### Important conceptual cautions already present

- "Feature-additive" does not mean scientifically trivial or conventionally univariate.
- A nonzero algorithmic interaction defect does not by itself demonstrate biological synergy.
- Conservation corrects double counting but not spatial blur or centre-assignment bias.
- The numerator of a statistic may collapse even when a locally estimated denominator does not.
- The novelty claim is provisional until a broader literature search is completed.

### Questions to carry into later parts

1. What exact assumptions on centering, nuisance projection, missing observations, and fold balance are required for each identity?
2. Should \(A\) be assumed symmetric without loss of generality, since only its symmetric part contributes to \(z_v^\top A z_v\)?
3. How should affine rather than strictly linear readouts and non-unit window sums be stated in the additivity criterion?
4. Does crossvalidated RSA regression remain a single fixed linear map when its noise covariance or variance normalization is estimated from local data?
5. What is the exact relationship between a fixed dense feature metric and the evidence-graph formulation? A dense \(M_s\) produces explicit pairwise terms but can still have a linear trace readout.
6. Which claims concern computational equivalence, which concern inferential equivalence, and which require identical permutation/null procedures?
7. How should boundary coverage, masked voxels with \(c_v=0\), and surface topology be handled in the conservation result?
8. Does conservation of signed crossvalidated evidence remain scientifically interpretable when local contributions can be negative?
9. What empirical example will best distinguish harmless additivity from scientifically useful feature interaction?
10. What prior work already uses convolution, sufficient-statistic maps, searchlight linearity, or partition-of-unity localization under different terminology?
11. Is "evidence graph" the best term, given that \(Q_A\) may be indefinite and its off-diagonal entries are signed?
12. Which result is the true centre of the note: collapse, design/feature duality, or conservation? A 2–3 page paper may need one primary theorem and the rest as corollaries/extensions.

## Living synthesis after Part 2 (historical snapshot)

### Essence in one paragraph

The framework now has two layers. Algebraically, a searchlight constructs a local Gram kernel; any fixed linear readout of that kernel collapses to spatial filtering of voxelwise quadratic evidence. Interpretively, the resulting energy can be decomposed in feature space into regional level and spatial configuration, and in population space into consensus and subject-specific heterogeneity. Mean removal is just projection away from the constant spatial mode: it leaves nonconstant marginal effects intact, introduces dense voxel coupling, and can turn the gradient of a smooth activation blob into an apparently “multivariate-only” boundary map. Actual behavioral confounds belong in observation-space operator \(A\); desired spatial invariances belong in feature metric \(M\). Configuration is not synonymous with genuinely joint information, which requires differences in the joint distribution despite matched marginals. The mature question is therefore not whether information is univariate or multivariate, but where it lies across spatial modes, metric coupling, population consensus, and marginally silent joint structure.

### Unified operator map

The central object is

\[
T(A,M)=\operatorname{tr}(AXMX^\top)=\langle M,X^\top AX\rangle_F.
\]

It organizes the whole argument:

| Object or operation | Mathematical location | Scientific role |
|---|---|---|
| Target contrasts, labels, folds, behavioral nuisance control | Observation-space \(A\) | Defines which relationships among observations count |
| Spatial window, mean projection, whitening, spatial basis | Feature-space \(M\) | Defines which voxel modes and relations count |
| Data | \(X\) or \(Z\) | Connects design and feature spaces |
| Evidence graph | \(X^\top AX\) | Stores nodewise and pairwise target co-evidence |
| Searchlight score | \(\langle M_s,Q_A\rangle_F\) | Scans the evidence graph with an anatomical metric |

For diagonal \(M_s\) and linear readout, the scan reduces to spatial filtering of node evidence. Dense \(M_s\) reads graph edges as well. Nonlinear \(F(K_s)\) adds algorithmic interactions beyond this bilinear operator.

### Integrated decomposition

For fixed-metric crossvalidated effect energy, the developing framework seeks an additive account along multiple axes:

\[
\begin{aligned}
\text{total evidence}
={}&\text{shared level}
+\text{shared configuration}\\
&+\text{idiosyncratic level}
+\text{idiosyncratic configuration}.
\end{aligned}
\]

Configuration can then be resolved into smooth, mesoscale, and fine graph-spectral modes. Separately, an analysis can ask whether information remains after all voxelwise marginals are matched, which is the stronger criterion for genuinely joint-distribution information. This final joint component is not automatically additive with the contrast-energy decomposition and may require different estimators.

### Strongest emerging claims

1. A large and useful class of searchlights does not require repeated local fitting.
2. Cross-validation itself does not block collapse; feature coupling, adaptive estimation, and nonlinear readouts do.
3. Demeaning removes one spatial direction, not “the univariate signal.”
4. The difference between local additive effect energy and energy of the local mean is exactly reproducible spatial covariance.
5. For smooth contrast fields, demeaned configuration energy has a leading gradient-energy interpretation.
6. Design confounds and feature invariances act on different sides of the same bilinear statistic.
7. Within-subject decoding energy combines shared population effects and between-subject heterogeneity.
8. Configuration, algorithmic feature interaction, and marginally silent joint information are distinct concepts.
9. Coverage correction can make local kernels conservative, while spatial spectral decomposition can replace arbitrary univariate/multivariate labels with a graded description.

### Publication architecture now visible

There are at least two plausible notes:

1. **Algebra/computation note:** the Gram-convolution identity, collapse criterion, evidence-graph duality, and conservation law.
2. **Interpretation/decomposition note:** the spatial-covariance identity, metric-aware level/configuration split, gradient/edge consequence, and consensus/heterogeneity split.

A single short note may be strongest if the Gram identity is introduced only as the machinery needed to derive the spatial-covariance and population decompositions. A longer methods paper could retain the full evidence-graph, conservation, scale-space, spectral, and nonlinear-interaction program. This scope decision remains open.

### New questions introduced by Part 2

1. Under what sampling assumptions is \(\widehat W_M-\widehat C_M\) unbiased for heterogeneity when each subject estimate itself contains noise and the consensus statistic does or does not use independent splits?
2. Must subject pairs be crossvalidated across runs as well as subjects to prevent measurement noise from entering \(\widehat C_M\)?
3. How do anatomical misregistration and functional alignment affect the interpretation of \(\Sigma_B\) versus true idiosyncrasy?
4. Is the weighted projector \(\Omega-\mathbf w\mathbf w^\top\) the unique natural configuration metric, and under what inner product is it a projector for unequal weights?
5. What normalization makes level and configuration comparable across searchlights of different size and weight concentration?
6. How accurate is the gradient-energy approximation near mask boundaries, anisotropic windows, curved surfaces, and fields with appreciable Hessian terms?
7. Can the edge-ring prediction be demonstrated with a simple analytic Gaussian blob and a realistic fMRI simulation?
8. How should signed crossvalidated covariance maps be inferred at the group level without converting them prematurely to unsigned accuracy?
9. Which forms of “joint information” can be tested with acceptable sample complexity once voxelwise marginals are controlled?
10. Does the term “SPM-like local activation effect” risk conflating smoothing-before-modeling, model-before-smoothing, and regional aggregation?
11. Should functional alignment be judged by conservation of within-subject energy, or can regularization and resampling legitimately change that energy?
12. How should observation-space nuisance projection be constructed inside cross-validation to prevent leakage?

## Living synthesis after Part 3 (historical snapshot)

### Essence

UEG treats the participant-by-partition effect matrix as the persistent scientific object and ordinary analyses as queries specified by three geometries: experimental operator \(A\), spatial operator \(M\), and generalization operator \(\Gamma\). Linear and quadratic queries can be written as contractions of these operators with a common effect tensor. Right- and left-Gram constructions expose activation and representation as dual geometries; a homogeneous moment matrix additionally joins signed first moments to quadratic effect energy and population variance. The earlier searchlight results become special cases: diagonal local operators collapse to spatial filters, coherent smoothing plus configuration covariance equals total local energy, and partitions of identity conserve global evidence. Population effects similarly split into consensus and heterogeneity operators. Rather than classifying analyses as univariate/multivariate or local/global, UEG describes moment order, support, spatial rank, adaptivity, scale, noise metric, and the axes along which effects generalize.

### The conceptual hierarchy now in force

1. **Data/effect representation:** retain \(B_{ir}\), rather than immediately reducing it to one map or score.
2. **Experimental query:** choose \(A\), including target contrasts and correctly cross-fitted nuisance control.
3. **Generalization query:** choose \(\Gamma\), naming the run, subject, item, session, task, or site relationships across which effects must reproduce.
4. **Spatial measurement:** choose \(M\), including locality, smoothing, metric, basis, scale, and invariance.
5. **Noise geometry:** choose measurement and population variability operators against which signal modes are normalized.
6. **Summary:** report signed directions, quadratic energy, spectral modes, consensus/heterogeneity, coherence, and scale; add nonlinear prediction only when scientifically required.
7. **Inference:** preserve the cross-fitting and sampling structure implied by the chosen operators.

### How the three parts now nest

| Part | Original centre | Role inside UEG |
|---|---|---|
| 1 | Searchlight Gram-convolution and collapse | Efficient evaluation of diagonal/fixed linear spatial queries; evidence-graph duality and conservation |
| 2 | Level/configuration and consensus/heterogeneity | Exact interpretable decompositions of spatial and population effect energy |
| 3 | Unified Effect Geometry | Parent architecture organizing design, space, generalization, first/second moments, spectra, and inference |

### Most important distinction added by Part 3

Part 2’s two-sided formula \(T(A,M)\) distinguished observation and feature space. Part 3 adds a third axis: generalization is not merely a procedural wrapper around a statistic but part of the estimand. This is arguably UEG’s most original architectural move. A score is incomplete unless it says both what effect was measured and across which independent units it had to reproduce.

### Scope boundary for the master equation

The equation exactly covers fixed bilinear/quadratic analyses once preprocessing and feature maps are fixed. It does not automatically absorb:

- parameters learned on the same test data;
- nonlinear normalization, ranking, thresholding, or accuracy;
- locally estimated covariance inverses unless training/test separation is explicit;
- arbitrary causal or generative claims;
- all higher moments without an explicit feature lift;
- inference validity without a sampling/noise model.

Such analyses may still fit UEG as cross-fitted operator selection followed by held-out evaluation, as nonlinear readouts of an effect operator, or as queries in a lifted feature space. The distinction must remain visible.

### Candidate empirical spine

The five-regime simulation is now the clearest common spine for both papers. It can test whether UEG correctly separates:

| Regime | Expected UEG allocation |
|---|---|
| Coherent blob | High level/coherence, shared low-frequency mode |
| Zero-mean gradient | Configuration energy from marginal mean effects, low DC coherence |
| Shared multiscale pattern | Consensus configuration distributed across spectral bands |
| Rotated individual patterns | High within-person and heterogeneity geometry, weak consensus |
| Covariance-only condition effect | Second-order off-diagonal effect with matched first-order marginals |

The demonstration should compare not only maps but recovered operator components, numerical conservation, bias under cross-fitting, and failure modes of accuracy-based summaries.

### Critical questions after Part 3

1. **Notation and symmetry:** What assumptions on symmetry and definiteness of \(A\), \(M\), and \(\Gamma\) are required? For cross-products, should every operator and estimate be explicitly symmetrized?
2. **Generalization operator design:** How are signs and normalizations of \(\Gamma\) chosen so that a crossvalidated statistic is unbiased and has a clear scale? Which \(\Gamma\) are PSD, indefinite, or graph-like?
3. **Changing queries:** There is one common effect matrix \(\mathbf B\), but not literally one fixed \(Q\) for every analysis; \(Q_{\Gamma,A}\) changes with the design and generalization query. The rhetoric should preserve this distinction.
4. **Linear activation extension:** The homogeneous moment state is an extension beyond the quadratic master equation. Should it be foundational from the start, with the master equation presented as a structured block query?
5. **Conventional SPM equivalence:** The coherent-energy expression models reproducibility of a smoothed contrast numerator, not automatically a complete SPM \(t\)- or \(F\)-statistic with estimated variance and preprocessing. Claims of exact equivalence must name the numerator/estimand.
6. **Learned operators:** When \(M\), \(w\), alignment \(P_i\), or spectral regularization is data-adaptive, what nested cross-fitting is required to retain unbiased generalization claims?
7. **Population estimators:** What are the variance, degrees of freedom, and limiting distributions of \(\widehat Q_W\), \(\widehat Q_C\), and \(\widehat Q_H\)? How are unbalanced partitions and missing cells handled?
8. **Indefiniteness:** Crossvalidated estimators and contrast-weighted \(\Gamma\) can yield indefinite \(Q\). Which UEG summaries require population PSD geometry, and which remain valid for signed evidence operators?
9. **Coherence ratio:** How should \(\kappa_h\) and consensus fractions be estimated without unstable ratios, positive bias, or values outside \([0,1]\) in finite samples?
10. **Rank-one binary claim:** The binary mean-difference signal operator is rank one by construction, but observed discriminability can include nuisance contrasts, covariance differences, nonlinearities, and multiple latent condition effects. The interpretation must be tied narrowly to the first-order two-condition estimand.
11. **Basis invariance versus anatomy:** A coordinate label is not invariant, but anatomical localization, measurement physics, topology, and smoothness make some bases scientifically privileged. UEG should reject ontological “univariate” claims without pretending all bases are equally meaningful.
12. **Tomographic identifiability:** Which spatial operator frames permit stable reconstruction of \(Q\), how many measurements are required, and what regularization is scientifically defensible under low-rank/noisy sampling?
13. **Connectivity extension:** When does \(Y^\top AY\) represent activation, covariance, connectivity, or merely an unnormalized cross-product? Directionality, temporal autocorrelation, normalization, and causal language require separate care.
14. **Separable generative model:** How restrictive is the Kronecker variance-component assumption, and can deviations be diagnosed without losing computational tractability?
15. **Sufficiency:** Under which first-level noise and HRF models is \(B=LY\) an adequate sufficient state for later UEG queries?
16. **Inference:** What null hypotheses and permutation schemes correspond to each \((\Gamma,A,M)\) query, especially with learned operators and crossed participant/item generalization?
17. **Novelty:** Which parts are repackagings of tensor regression, kernel machines, variance-component models, dual PCA/CCA, crossnobis/cvMANOVA, representational geometry, graph signal processing, or operator-valued frames, and which synthesis is genuinely new?
18. **Paper scale:** Can the larger framework be made falsifiable and useful rather than appearing to rename all quadratic methods? The five propositions need executable estimators and at least one result unavailable or obscure in existing formalisms.

### Recommended claim discipline preserved in the ledger

- Use “a broad class of fixed-metric linear and quadratic analyses,” not “all fMRI analyses.”
- Use “SPM-like coherent effect numerator” when the denominator and full inferential machinery are not included.
- Call \(Q\) an evidence/effect operator without assuming PSD unless its construction guarantees PSD.
- Treat nonlinear decoders as learned or lifted extensions, not automatic instances of the fixed master query.
- Separate population identities from finite-sample crossvalidated estimators.
- Present activation tomography as a proposed method requiring identifiability analysis, not an established consequence of conservation alone.
- Preserve the anatomical meaning of the voxel basis while rejecting basis-dependent labels as fundamental signal properties.

## Living synthesis after Part 7 (historical snapshot)

### Essence

The conversation has converged on a relation-first geometry compiler. A partition exposes a lazy experimental–neural relation \(B_r=E_rY_r\), regardless of whether it is backed by raw time series, betas, condition averages, latent features, or another estimator. A spatial frame measures that relation; a pairing graph states where independently estimated relations must agree; and scientific queries read the resulting experimental-space geometry. In the fast additive case, the entire package compiles to \(V=W\Phi_\Gamma(EY)C\). Coherent regional effect and within-region configuration are exact simultaneous components, local/global scope is a frame choice, RSA and contrasts are compiled views, and scalar maps appear only at the boundary. Native geometry is group-ready without voxel-basis alignment. The clean-room R package therefore needs only two durable IRs—`relation` and `geometry`—plus small frame, pairing, and query values and one block-streamed contraction.

### Architecture frozen provisionally

```text
response sources Y_r
        +
effect extractors E_r
        |
        v
lazy relation B_r = E_r Y_r
        |
        +-- spatial frame W or factors L_j
        +-- partition pairing Gamma
        |
        v
cross-generalized geometry G_j
        |
        +-- signed/coherent/configuration/total
        +-- contrast/RDM/RSA/spectrum/information views
        +-- location transport / later population model
```

### Clean-room relationship to rMVPA

The new package seeks similar scientific coverage, not feature parity. rMVPA remains the evidence base for which user tasks matter and which architectural pressures recur. No rMVPA code, class structure, or compatibility layer should be copied. The replacement compresses workflow branching into algebraic composition:

| rMVPA pressure | `effectgeom` response |
|---|---|
| dataset/design/model spec coupled into analysis | lazy relation with explicit source and extractor |
| regional/searchlight/global modes | one spatial frame abstraction |
| procedural cross-validation objects | sparse partition-pair graph |
| classifier/model registry | fixed geometry core; learners return frozen transforms |
| model-specific RSA and MANOVA paths | compiled linear geometry queries |
| per-method result and combiner classes | one packed geometry plus late views |
| fast-path eligibility/fallback | compiler lowering based on operator representation |
| analysis-specific dataset restrictions | operation-local validation |
| map/prediction scalar as durable result | geometry as durable IR |

### MVP decision

The strongest first release is narrower than Parts 5–6 initially suggested:

- matrix/list/file-backed response sources;
- identity and explicit linear extractors, plus stable `lm_extractor()`;
- abstract and volume domains first, surface as suggested adapter if it does not delay the kernel;
- additive sparse frames for points, searchlights, regions, and whole brain;
- explicit off-diagonal partition pairings;
- packed total/coherent geometry and signed weighted relation means;
- contrast bundle and squared-RDM/linear-RSA views;
- dense and streamed reference implementations;
- conservative/local frame normalization;
- `as_neuro()` at the edge.

Information spectra, factor frames, adaptive metrics, native geometry transport, surface adapters, and population models should follow only after the core laws and units are proven. This sequencing protects the sub-3,000-line goal.

### Highest-risk mathematical/API decisions

1. **Orientation:** Parts 3–6 used both \(q\times p\) and \(p\times q\) conventions. The public API must be orientation-free; the internal canonical convention should be fixed once and tested at adapters.
2. **Pair normalization:** `cross_partitions()` needs an exact rule for ordered versus unordered edges, symmetrization, and total weight so estimates have predictable units.
3. **Signed mean:** a cross-generalized geometry does not by itself determine a signed local effect. The `geometry$mean` field needs a declared partition aggregation rule and is not necessarily crossvalidated evidence.
4. **Indefiniteness:** cross-geometries can be indefinite. Spectrum/information functions must refuse or model them rather than silently clipping.
5. **Metric semantics:** additive weights are a measurement metric, but noise whitening and local dense metrics require factor frames and training provenance. They cannot be smuggled into arbitrary weights.
6. **Coherence fraction:** ratios can be unstable or outside population bounds in crossvalidated finite samples. It should initially be descriptive with uncertainty/validity warnings, not a default inferential statistic.
7. **Raw extractor validity:** `lm_extractor()` must use numerically stable QR/SVD, expose estimability/rank, and be clear about whether temporal precision is supplied, estimated elsewhere, or fixed.
8. **Cross-partition independence:** metadata can record and validate declared independence but cannot prove it from arbitrary sources. The API must say “declared/structurally checked,” not certified.
9. **RSA intercept/centering:** compiled regression views must reproduce an explicitly documented RDM vectorization, intercept, weighting, nuisance, and scaling convention.
10. **Information units:** \(\frac12\log\det(I+G)\) is valid only for appropriate PSD, noise-normalized relation geometry. It is not a generic transformation of beta-energy fields.
11. **Conservation:** matrix geometry and linear views conserve under a partition of identity; nonlinear information, accuracy, ranks, and ratios generally do not.
12. **Group transport:** averaging transported low-dimensional geometry is valid as a linear operation, but uncertainty in correspondence and nonlinear manifold summaries require later models.

### Package constitutional test

Before any feature enters core, ask:

1. Can it be expressed by a source/extractor relation, a spatial frame/factor, a pairing graph, or a geometry query?
2. Can it reuse the single cross-Gram kernel?
3. Does it preserve one `geometry` result type?
4. Can its semantics and units be stated independently of algorithm name?
5. If learned, can training be frozen and evaluation separated?

If not, it belongs upstream, downstream, or outside the package.

### Next concrete artifact

The right next step is not package scaffolding. It is a small executable mathematical reference containing:

- dense `relation_block` and identity extractor;
- `svec`/`unsvec` with trace-isometry tests;
- explicit edge-normalized cross atoms;
- dense frame contraction;
- total/coherent/configuration laws;
- direct-query/full-geometry equivalence;
- raw/effect equivalence;
- a brute-force searchlight comparison.

Only after this kernel passes simulations should public class names or neuroimaging adapters be stabilized.

## Living synthesis after Part 8 (historical snapshot)

### Essence

`effectgeom` now has a unified scientific and operational architecture. Scientifically, a lazy relation \(B_r=E_rY_r\) is measured by a spatial frame, generalized by an explicit partition pairing, and read through compiled geometry queries. Operationally, the compiler lowers this to canonical feature-block tasks and one deterministic sparse contraction. Searchlights, ROIs, voxels, and the whole brain are never execution modes or job types. The default executor is sequential and memory-bounded; optional concurrency operates only over disjoint feature blocks or participants, under explicit source-locality and resource contracts. No global scheduler state, implicit dataset serialization, backend fallbacks, per-location error objects, or premature result rows are required. A completed geometry carries one execution receipt proving how the plan ran.

### Full architecture

```text
Scientific layer
  response Y + extractor E -> lazy relation B
  relation + frame W + pairing Gamma + query C

Compiler layer
  validate -> choose full/direct query -> estimate memory
  -> canonical feature blocks -> sparse contraction

Execution layer
  sequential by default
  optional owned workers for worker-readable sources
  canonical reduction + structured events + receipt

Result layer
  complete effect_geometry
  -> signed/coherent/configuration/total/RDM/RSA/spectrum views
```

### The key improvement over rMVPA execution

rMVPA has had to optimize a large number of local fits. `effectgeom` removes the local fits for its admitted bilinear class. Consequently, execution engineering focuses on data movement and deterministic contraction rather than scheduling thousands of analysis objects. This is why the package can improve speed, memory, reproducibility, failure behavior, and code size simultaneously.

### New execution laws

1. Sequential, streamed, and parallel executions agree with the dense oracle.
2. Worker count and task completion order do not alter ordered reduction semantics.
3. A one-worker request uses a direct loop and creates no worker infrastructure.
4. The full source is never replicated implicitly per worker.
5. No global `future` plan, RNG state, or unrelated process state is mutated.
6. Every owned worker, handle, and temporary artifact is released on success, failure, cancellation, and interrupt.
7. Partial execution never returns an ordinary complete geometry.
8. Execution receipt and scientific plan identities are preserved separately.
9. Parallel execution ships only after end-to-end benchmarks against the vectorized sequential compiler.

### Scope impact

The execution design adds a small compute-policy value and internal operational manifest but does not add another scientific object or workflow platform. It may modestly increase the initial code budget; that cost must be offset by keeping parallelism sequential-only in E0 and deferring process execution/checkpoint adapters. The under-3,000-line goal remains plausible for the sequential compiler, not necessarily for every later executor.

### Remaining execution questions

1. Are relation extraction and atom construction expensive enough to benefit from R processes after vectorization?
2. Which first file-backed source format gives efficient feature-major concurrent reads without adding a heavy dependency?
3. Should owned local processes use base `parallel`, or should every parallel backend remain an optional adapter?
4. What deterministic summation method gives sufficient cross-platform reproducibility without unacceptable cost?
5. How accurately can peak R memory be predicted given copy-on-write, sparse multiplication temporaries, and BLAS behavior?
6. Should participants be parallelized only by callers, or should a later cohort API own that outer execution axis?
7. What is the minimal event/receipt schema that remains stable without becoming a logging subsystem?

## Living synthesis after Part 9 (historical snapshot)

### Essence

The scientific compression and the systems design now reinforce one another.
`effectgeom` has one relation/geometry compiler and one pure feature-block task.
Sequential streaming remains the executable oracle. When R process isolation
would otherwise duplicate large matrices, `shard` becomes the preferred
optional transport and execution substrate: publish a response or compiled
relation once, dispatch coarse feature blocks to supervised workers, and return
bounded canonically indexed intermediates. The geometry kernel, reduction,
result, and interpretation do not branch. We retain shared memory, supervision,
and telemetry while deleting rMVPA's model-specific shard backend.

### Execution stack

```text
Scientific IR
  relation B = EY + frame W + pairing Gamma + query C

Compiler
  validate -> select staging -> choose blocks -> compile one task

Execution
  sequential direct loop (oracle/default)
  or optional shard adapter
    shared immutable Y or B
    supervised coarse feature tasks
    bounded block results/buffers

Canonical effectgeom reduction
  W Phi_Gamma(EY) C in fixed feature order

Result
  one complete effect_geometry + execution receipt
```

### Authoritative decisions

1. Searchlights and ROIs are frame rows, never jobs.
2. `future` does not define package state or semantics.
3. `shard` is valuable and is the preferred first local parallel adapter after
   the sequential compiler is proven.
4. `shard` is optional and operational; it creates no new scientific class,
   method family, kernel, or result.
5. Shared staging may occur at \(Y\) or \(B=EY\); the choice is explicit and
   provenance-recorded.
6. Workers read immutable inputs and produce bounded disjoint intermediates;
   only `effectgeom` performs final overlapping spatial reduction.
7. Shared-resource ownership, generation identity, and cleanup are correctness
   laws.
8. A one-worker request bypasses process/shared-memory machinery by default.
9. Parallelism is admitted on aggregate memory, serialization, cleanup,
   numerical parity, and wall-time evidence—not speed alone.

### Remaining implementation questions

1. Does the first adapter gather a small number of block results, use a shared
   atom/query buffer, or require a small streaming-result addition to `shard`?
2. Which staging threshold predicts when sharing \(B\) beats sharing \(Y\)?
3. Should reusable caller-owned pools be supported initially, or deferred until
   ownership semantics are exhaustively tested?
4. Which `shard` telemetry fields become stable fields in the `effectgeom`
   execution receipt?
5. What platform/backing combinations form the supported matrix for the first
   adapter?

## Living synthesis after Part 10

### Essence

`effectgeom` is now specified as a deliberately narrow version-0.1 compiler,
not a universal promise about every multivariate workflow. Its exact core is an
additive diagonal spatial frame plus fixed bilinear generalization/query
geometry. Complete geometry and direct query have distinct result contracts.
The sequential compiler tiles feature, measurement, and packed-geometry axes,
and may use block-backed complete output; its memory plan includes contraction
copies, serialization, reorder, checkpoint, and storage costs and must be
validated empirically. Scheduling is bounded, numerical guarantees are
tolerance-qualified, undirected signed summaries are orientation invariant, and
version 0.1 rejects parallel execution. `shard` remains the preferred later
adapter after the sequential and scheduling gates are earned.

### Contract stack

```text
Exact scientific scope
  additive frame + fixed bilinear pairing/query

Public results
  complete effect_geometry (memory or block storage)
  direct effect_view (one declared query)

Sequential compiler
  feature blocks x measurement tiles x geometry-coordinate tiles
  conservative memory plan checked against measured peak

Scheduling contract
  canonical order + bounded inflight/reorder window + backpressure
  tolerance-qualified numerical agreement

Later execution
  in-process simulator -> benchmark gate -> optional shard adapter
```

### Immediate work order

1. Close theorem scope, endpoint marginal, result-type, public-API, and v0.1
   contracts.
2. Implement the dense oracle and 20 algebraic/result laws.
3. Implement tiled sequential contraction and block-backed complete output.
4. Measure actual allocation and peak memory across realistic regimes.
5. Prove bounded scheduling with an in-process adversarial simulator.
6. Add checkpoint integrity only after task artifacts stabilize.
7. Benchmark `shard` against the complete sequential baseline.

## Revision and contradiction log

### Part 2 relative to Part 1

- **Strengthens:** The feature-additive/feature-interactive taxonomy is embedded in a broader hierarchy. Feature additivity is no longer the endpoint; the framework now distinguishes regional level, nonconstant marginal configuration, metric coupling, nonlinear algorithmic interaction, and marginally silent joint structure.
- **Strengthens:** The abstract operator \(T(A,M)\) gains a practical interpretation: left-side operations control observation/design relationships, while right-side operations specify spatial invariance and feature geometry.
- **Refines:** Part 1’s statement that demeaning/normalization can prevent scalar collapse needs precision. Mean-centering with locally normalized weights is itself a linear but searchlight-dependent dense feature operator. It may not collapse to simple smoothing of a single voxelwise map, yet it remains exactly decomposable as additive total minus rank-one regional-level evidence.
- **Refines:** Dense feature metrics are feature-interactive in the node-only sense because they read off-diagonal evidence-graph terms, but the overall statistic can remain linear in the Gram/evidence graph. “Feature interaction” must distinguish metric coupling from nonlinear readout curvature.
- **Adds:** Part 2 introduces an across-subject axis absent from Part 1. Within-subject energy and population-mean evidence are different estimands, with their difference attributable to heterogeneity under stated assumptions.
- **Adds:** The smooth-field Taylor expansion predicts edge/ring maps after local mean removal, providing a decisive simulation and interpretive example.
- **Changes publication emphasis:** Part 1 centred the computational collapse theorem; Part 2 proposes the spatial-covariance identity and population decomposition as the more compelling scientific story. Whether these form one note or two remains unresolved.
- **No direct contradiction identified:** Part 2 extends and qualifies Part 1 rather than negating its core trace identities.

### Part 3 relative to Parts 1–2

- **Promotes:** The common effect matrix \(B_{ir}\), not the local Gram field, becomes the primitive object. Searchlight kernels and maps become derived spatial queries.
- **Generalizes:** \(T(A,M)\) becomes \(T(\Gamma,A,M)\), adding explicit generalization geometry across partitions, participants, stimuli, tasks, sessions, and sites.
- **Unifies first and second moments:** The homogeneous effect state \(\mathbb H\) is added so signed activation and quadratic geometry can be described as blocks of one moment object.
- **Strengthens duality:** Evidence-graph trace cyclicity is extended to singular-value pairing of experimental and spatial modes.
- **Renames:** Regional level versus configuration is reframed as coherent versus incoherent integration, with coherence fraction \(\kappa_h\) proposed as a quantitative summary.
- **Generalizes population decomposition:** Scalar consensus/heterogeneity energies become operator-valued geometries \(Q_C\) and \(Q_H\).
- **Extends conservation:** Coverage correction becomes a general resolution-of-identity framework and motivates the speculative activation-tomography program.
- **Extends scale:** Searchlight radius and local spectra become part of a conserved spatial-frequency measurement system.
- **Extends moment order:** Genuinely joint information is represented through explicit second- or higher-order feature lifts rather than inferred from voxel count or demeaning.
- **Extends scope:** Activation, representation, decoding, population inference, alignment, and some covariance/connectivity analyses are proposed as members of one calculus.
- **Qualifies Part 1 language:** Not every named method is simply a measurement of one already-estimated \(Q\); changing \(A\) or \(\Gamma\) changes the effect operator, while learned metrics and nonlinear summaries require cross-fitting or extensions.
- **Qualifies Part 2’s SPM analogy:** The exact decomposition concerns coherent crossvalidated effect energy. A conventional SPM statistic adds variance estimation and other choices that do not automatically share the same algebra.
- **Changes publication scale:** The short covariance/searchlight note remains viable, but a separate full framework paper is now proposed.
- **No algebraic contradiction identified:** The earlier trace, covariance, population, and conservation identities are retained as special cases of the expanded architecture.

### Part 4 relative to Part 3

- **Promotes:** The brain–experiment relation, not the effect matrix or spatial effect operator, becomes the foundational object.
- **Adds:** Encoding and decoding are equal conditional density ratios under a coherent joint model; effect, generalization, and calibration are separated.
- **Sharpens:** Activation and representation become pushforward/pullback Gram geometries of one noise-normalized relation operator.
- **Adds:** Native neural coordinate systems are treated as gauges; experimental pullback geometry supports group comparison without fiber alignment.
- **Adds:** A population hierarchy distinguishes shared information, shared geometry, transportability, and fixed-coordinate topography.
- **Adds:** Bures geometry and base-versus-fiber transport motivate native-space geometry fields.
- **Qualifies:** Direction-free statistical evidence does not erase causal experimental direction.
- **Qualifies:** Gauge invariance requires invertible coordinate changes and coherent metric transformation; it does not cover lossy measurement.
- **Qualifies:** Bures and information summaries require latent PSD geometry, whereas unbiased crossvalidated estimates may be indefinite.

### Part 5 relative to Part 4

- **Materializes:** The theory becomes a geometry compiler between external effect estimation and late scientific views.
- **Introduces:** `EffectRelation`, `NeuralMeasure`, `MeasurementFrame`, `GeneralizationGraph`, `GeometryEstimate/Field`, `DesignView`, and immutable `GeometryPlan` as language-neutral concepts.
- **Freezes learned analyses:** adaptive methods train externally and return immutable measures evaluated on held-out relations.
- **Separates:** location correspondence from embedding transport in both math and API.
- **Makes laws executable:** adjoint, trace, gauge, frame, decomposition, streaming, query, null-centering, and no-implicit-transport properties define backend correctness.

### Part 6 relative to Part 5

- **Narrows:** The broad compiler becomes a clean-room R package aimed at a small, fixed bilinear core rather than a general multi-backend platform.
- **Compresses:** Additive frames, packed geometry atoms, and compiled queries produce the normal form \(WZC\).
- **Demotes:** Classifier and RSA identities become views or external learners instead of primary classes.
- **Targets:** Similar scientific coverage to rMVPA in roughly one quarter of the code, with no compatibility requirement.

### Part 7 relative to Part 6

- **Replaces:** `effect_set` with a lazy relation \(B_r=E_rY_r\); beta images are no longer mandatory inputs.
- **Reduces:** Five proposed durable objects become two IRs (`relation`, `geometry`) plus small declarative frame/pairing/query values.
- **Adds:** Raw/effect equivalence as the first implementation law.
- **Corrects:** The coherence metric for arbitrary row mass uses \(ww^\top/a\).
- **Tightens:** Bilinear closure becomes the explicit rule guarding core scope.
- **Narrows further:** The production target drops from about 5,000 to under 3,000 R lines; group modeling becomes a later/separate concern.
- **Rejects:** A formula DSL, arbitrary per-ROI callbacks, workflow platform, mandatory beta stage, and engine-specific fast paths.
- **Supersedes:** Where Part 6 says `effect_set`, read Part 7’s `relation`; where Part 6 proposes five central objects, Part 7’s two-IR design is authoritative.

### Part 8 relative to Parts 6–7

- **Rejects a false substitution:** Replacing `future` with another parallel-map library would retain the wrong per-searchlight unit of work.
- **Promotes feature blocks:** Concurrency is expressed over canonical relation/atom blocks, while spatial measurements remain rows of one sparse contraction.
- **Adds an explicit operational layer:** `compute_policy`, task manifests, source capabilities, structured events, and execution receipts are operational values, not new scientific IRs.
- **Makes sequential primary:** The deterministic vectorized compiler ships before any parallel backend; parallelism must beat it empirically.
- **Separates source locality:** coordinator-read, worker-read, and explicit shared sources implement one protocol instead of separate dataset/engine hierarchies.
- **Strengthens reproducibility:** canonical reduction order replaces schedule-dependent reduction, and the core contains no RNG.
- **Strengthens failure semantics:** structural issues fail before work; runtime task failures stop the plan; incomplete output is not a normal geometry.
- **Strengthens ownership:** the package never inherits or mutates global `future` state and closes only resources it creates.
- **Defers:** checkpoint/resume, local-process execution, and external executor adapters until after the sequential compiler and task protocol are verified.

### Part 9 relative to Part 8

- **Corrects an overreach:** The problem is not shared memory or `shard`
  itself. rMVPA's shard path demonstrates real protection against per-worker
  duplication of large R objects.
- **Retains the mechanism:** Immutable shared inputs, supervised workers,
  recycling, bounded buffers, deterministic reduction support, and memory/copy
  telemetry are valuable infrastructure.
- **Rejects only the branching:** There remains no shard-specific scientific
  engine, dataset extractor hierarchy, model spec, or result type.
- **Promotes `shard`:** After the sequential oracle, it becomes the preferred
  first local executor/source adapter rather than a rejected or generic future
  possibility.
- **Adds two staging points:** Either share raw response matrices \(Y_r\), or
  materialize and share compact relations \(B_r=E_rY_r\) for reuse.
- **Strengthens lifecycle laws:** owner/attachment roles, generation identity,
  stale-handle prevention, idempotent cleanup, and aggregate-memory evidence
  become explicit correctness requirements.
- **Preserves Part 8's central claim:** Feature blocks are tasks; searchlights
  remain rows of \(W\), and final scientific reduction stays canonical.

### Part 10 relative to Parts 8–9

- **Scopes:** “Searchlights are not tasks” is now explicitly limited to
  additive diagonal frames with fixed bilinear queries.
- **Corrects memory:** Feature blocks do not bound dense sparse-product output;
  row/coordinate tiling and block-backed complete storage are now required
  options.
- **Replaces exactness language:** Portable exact memory accounting and bitwise
  policy invariance are replaced by conservative measured-memory planning and
  tolerance-qualified numerical contracts.
- **Bounds scheduling:** `max_inflight`, `max_reorder_bytes`, dispatch
  backpressure, prompt release, and spill rules make canonical reordering an
  actual bounded mechanism.
- **Separates results:** Full materialization returns `effect_geometry`; direct
  query returns `effect_view`.
- **Freezes release scope:** Version 0.1 rejects `workers != 1`; `shard` remains
  a benchmark-gated later adapter.
- **Purifies policy:** reporters and checkpoint paths are removed from immutable
  numeric `compute_policy` identity.
- **Hardens recovery:** failure conditions carry receipts; checkpoints add disk,
  checksum, permission, source-integrity, and exact resume-plan contracts.
- **Corrects marginals:** undirected edges use one endpoint marginal; left/right
  is reserved for directed pairings.
- **Operationalizes:** three live mote epics and fifteen dependency-ordered
  children now encode semantic, sequential-memory, and later-runtime work.

## Verification ledger

### Established within the conversation by direct algebra

- Gram atoms pool linearly under a diagonal spatial window.
- Trace readouts commute with that pooling.
- A fixed bilinear fold construction can be symmetrized into an observation-space operator.
- Coverage-normalized diagonal windows sum to the global Gram matrix when every included voxel has positive coverage.
- Trace cyclicity gives the feature-space/observation-space dual expression.
- The normalized diagonal feature metric splits into regional-level and centered-configuration PSD components.
- The additive searchlight product equals product of local means plus weighted cross-map covariance.
- Population mean energy plus metric-weighted between-subject covariance equals expected subject-specific quadratic energy.
- Encoding and decoding density ratios are pointwise identical when derived from one coherent conditional joint distribution.
- Pullback relation geometry is invariant under invertible native-coordinate changes when the neural metric/noise transforms coherently.
- Raw response plus linear extractor and materialized effect relations are algebraically equivalent inputs to the same geometry kernel.
- Additive spatial frames and packed symmetric queries yield the compiler normal form \(W\Phi_\Gamma(EY)C\).
- The arbitrary-mass coherent/configuration split follows from \(D(w)=ww^\top/a+[D(w)-ww^\top/a]\).

### Confirmed against live rMVPA source on 2026-08-11

- Audited public commit `3b12a8855b06f549e5772ccb9af070639a34752b` in a read-only shallow checkout.
- `DESCRIPTION` collated 82 R files and listed 23 imports; `R/` contained about 49,441 physical lines and `NAMESPACE` about 190 exports.
- README and workflow source confirmed dataset/design/model/engine/result framing and explicit searchlight/regional/global dispatch.
- Searchlight source confirmed legacy, SWIFT, dual-LDA, and naive cross-decoding engine policies with eligibility/fallback logic.
- ROI/result source confirmed conversion of structured `roi_result` values to the tibble format expected by older combiners.
- Dataset source confirmed Feature-RSA-specific multiple-feature validation in the general image dataset constructor.
- Execution source confirmed serial main-process ROI preparation followed by `future_pmap`, analysis-type-dependent batch/chunk rules, global-plan worker discovery, and worker-spec dataset stripping.
- The experimental shard path confirmed a second dataset-specific shared-memory preparation/extraction layer and a duplicated `run_future` implementation.
- CLI/custom paths confirmed temporary mutation/restoration of the global future plan; remote shard workers confirmed progress-closure limitations.
- The audited iterator/shard/engine/fast-path/global/custom/CLI/utilities execution surface exceeded 8,000 physical R lines.

### Confirmed against live shard source on 2026-08-12

- Audited public commit `233c71186ebac0e2e98eba70f10671450fc5e1be`
  (`shard` 0.2.0) in a read-only shallow checkout.
- Public APIs include immutable shared inputs, shared/memory-mapped segments,
  reopenable descriptors, explicit buffers, worker pools, `shard_map()`, and
  deterministic `shard_reduce()`.
- Dispatch code supervises worker health, requeues work after worker loss, and
  replays published worker state after restart/recycling.
- Reduction documentation and tests require chunk-ordered master reduction and
  worker-count reproducibility at fixed chunking.
- Diagnostics expose worker memory, copy/materialization, recycling, timing,
  and end-of-run behavior—the right evidence surface for an optional
  `effectgeom` adapter.
- This source inspection establishes capability, not yet fitness or stability
  for `effectgeom`; integration and lifecycle stress tests remain required.

### Still requiring independent mathematical or implementation checking

- Exact equivalence for every named RSA and cvMANOVA variant.
- Handling of scaling constants, degrees of freedom, intercepts, centering, and variance estimators.
- Equivalence of complete inferential workflows, including permutations and uncertainty estimates.
- Runtime and memory advantages in realistic regimes.
- Stability and usefulness of interaction-defect measures.
- Scientific utility of conservative kernel decomposition and scale-space persistence.
- Exact unbiasedness and variance of the proposed within-subject, cross-subject, and heterogeneity estimators under finite noisy samples.
- Metric-aware projection formulas under singular, regularized, estimated, or subject-varying precision matrices.
- Range of validity and boundary corrections for the smooth-field gradient-energy approximation.
- Whether the proposed four-way maps remain inferentially well calibrated after alignment and spatial localization.
- Exact finite-sample equivalence between relation-spectrum roots and named \(t/F\), prediction, and likelihood summaries under realistic nuisance/noise estimation.
- Statistical properties of Bures population models built from noisy or crossvalidated geometry estimates.
- Identifiability and uncertainty of base correspondence versus functional fiber transport.
- Canonical pairing normalization, directed marginal semantics, and variance estimates for the proposed package.
- Numerical behavior of packed streaming accumulation, rank-deficient extractors, and conservative sparse frames in a prototype.
- Whether the under-3,000-line implementation target survives documentation-quality validation and adapters.

### Still requiring literature verification

- Priority for matrix-valued convolution framing.
- Priority for the general collapse/commutation theorem in searchlight analysis.
- Prior feature-additive/feature-interactive taxonomies under other names.
- Prior evidence-graph or design/feature trace dualities applied spatially.
- Prior partition-of-unity conservation arguments for representational searchlights.
- Prior exact spatial-covariance accounts of the mass-univariate/additive-searchlight discrepancy.
- Prior level/configuration by consensus/heterogeneity decompositions and corresponding U-statistic estimators.
- Prior graph-spectral decompositions of local representational evidence.

## Implementation checkpoint: execution architecture (2026-08-12)

The execution design has now been reduced to one scientific kernel with a
separate runtime boundary:

- `crossgram_feature_task()` is the canonical pure feature-block task. It
  receives relation blocks and pairing data, and returns bounded relation and
  packed-Gram contributions without sources, frames, executors, or output
  objects.
- One ordered coordinator reducer is shared by sequential execution and any
  future executor. Completion order may vary, but scientific reduction order
  is canonical and numerical equivalence is tolerance-qualified.
- Reopenable source descriptors provide immutable, serializable metadata and
  strong content revisions. Worker handles own only their attachments; closing
  them never deletes the backing source.
- Complete geometries and direct query results remain distinct public
  contracts. Output tiling, bounded reordering, checkpoint integrity, and
  failure receipts are explicit parts of the execution design.

The optional-process decision was settled by measurement rather than taste. A
clean-room shard benchmark used the same canonical task and reducer for
sequential, shared-response, and shared-relation execution, including cold and
warm reuse. All modes were numerically equivalent and cleaned up normally, and
the injected-failure cleanup test passed. None met the predeclared 1.1x cold
speedup gate: observed shard speedups ranged from 0.258x to 0.794x, although
process-tree peak RSS was lower. Therefore shard is not admitted into the core
or version 0.1.

`future` is also rejected as a core-owned executor because its ambient global
plan conflicts with explicit executor ownership. `mirai` has a technically
compatible named-daemon lifecycle, but is deferred because no process backend
has yet earned admission for the measured workload. Batch schedulers remain a
coarse participant/plan orchestration concern, not an inner feature-block
engine.

The resulting release decision is intentionally conservative:

> Version 0.1 is sequential, memory-bounded, checkpointable, and executor-ready,
> but it does not ship an unearned parallel adapter.

This is not a retreat from the original goal. The package first removes the
redundant searchlight work algebraically; it parallelizes only the irreducible
remainder when a benchmark demonstrates an end-to-end benefit under the same
scientific estimand.
