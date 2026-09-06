# Relation to prior work

Status: positioning note, normative for public claim language

Date: 2026-08-16

This note states precisely which existing statistic each part of `crossform`
is a relative of, and what is left over. It exists so that public claims can
be checked against the literature rather than against enthusiasm. The
evidence ledger in `vignettes/novelty.Rmd` governs *status*; this note governs
*attribution*.

## The statistic

For relations `B_r` estimated within partitions `r`, a fixed neural metric
`K`, a fixed experimental query `H`, and a pairing weight `Γ`, the package's
irreducible observable is

    E(H, K) = Σ_{r,s} Γ_rs · tr(Hᵀ B_r K B_sᵀ),    Γ_rr = 0.

Everything public — crossvalidated contrast energy, the crossnobis RDM, a
fixed linear RSA coefficient, adjoint coupling — is a query of this form.

## Cross-validated MANOVA (Allefeld & Haynes, 2014)

This is the closest existing statistic to crossvalidated contrast energy and
had not previously been cited here. Their **pattern distinctness** `D̂` is a
leave-one-run-out trace `tr(H_l E_l⁻¹)` combining a contrast estimated from
the held-out run with one estimated from the training runs, under an error
covariance metric. It is the same *shape* of object: a quadratic form in a
contrast, crossvalidated across runs, evaluated under a noise metric. Four
things differ, and none of them is cosmetic.

1. **Normalization.** `D̂` is a Bartlett–Lawley–Hotelling-type trace: the
   contrast sum-of-squares is normalized by the estimated error term `E⁻¹`.
   `crossform`'s total is the unnormalized pairing `⟨H, B_r K B_sᵀ⟩` under a
   metric the analyst declares. With `K` set to an inverse noise covariance
   the two coincide in spirit, but the package does not build the error term
   into the statistic's definition.
2. **Bias correction.** cvMANOVA removes the bias of `E⁻¹` as an estimate of
   `(nΣ)⁻¹` with an explicit multiplicative factor,
   `((m−1)f_E − p − 1) / ((m−1)n)`. `crossform` has no multiplicative
   correction. Its unbiasedness route is structural: every product multiplies
   estimates from two *different* partitions, so no squared noise term is ever
   formed (see `design/effect-form-contract.md` §8).
3. **Voxel-count standardization.** cvMANOVA recommends mapping
   `D̂_s = D̂/√p`, because the null variance of `D̂` grows with the number of
   voxels. `crossform` standardizes by the measurement frame's own weight mass
   under `normalization = "local"`. That is a declaration about what a
   measurement is, not a null-variance correction, and it does not equalize
   null variance across searchlights of different sizes.
4. **The voxel-axis split.** Allefeld and Haynes' Figure 6 already runs the
   mean-versus-pattern control: cvMANOVA on searchlight-mean-only data, and on
   mean-removed data. Those are two re-analyses of two modified datasets, and
   their results do not sum to the unmodified analysis. `crossform` instead
   partitions the *metric* — `D(w) = wwᵀ/a + [D(w) − wwᵀ/a]`, both parts PSD —
   so `total = coherent + configuration` holds exactly from one computation and
   is inherited by every fixed linear query of the same form.

cvMANOVA also supplies distributional theory and permutation inference.
`crossform` supplies no inference layer, and says so.

## Crossnobis / linear discriminant contrast

The crossvalidated Mahalanobis (LDC) distance and its reliability advantage
over correlation distance and classification accuracy are established
(Walther et al., 2016), as is its sampling distribution under an equal-partition
model (Diedrichsen, Provost & Zareamoghaddam, 2016) and the second-moment
framework relating encoding models, PCM, and RSA (Diedrichsen & Kriegeskorte,
2017). `crossform`'s `rdm()` **is** that estimator — the Haxby exemplar matches
an independent implementation to `8.88e-16` — and none of it is claimed as new.
The only difference is architectural: the distance is one fixed linear query of
a retained form rather than the object the analysis is built around.

## RSA tooling and inference

Searchlight information mapping (Kriegeskorte, Goebel & Bandettini, 2006), the
RSA toolbox and its RDM comparison practice (Nili et al., 2014), whitened
unbiased RDM similarity (Diedrichsen et al., 2021), and inference on
representational geometries over subjects and conditions (Schütt et al., 2023)
are prior art that `crossform` neither reimplements nor claims. Multivariate
connectivity has its own review (Basti et al., 2020); the package's adjoint
coupling is a materialization of the same pairing, not a new connectivity
method.

## Framed RSA (Taylor & Kriegeskorte, 2025)

This is the nearest published construction to the coherent/configuration
split, and it postdates most of the mean-versus-pattern literature (Davis et
al., 2014). Framed RSA restores the regional-mean information that ordinary
RSA discards by **augmenting the pattern set** with two reference patterns:
the origin and a uniform constant pattern. `crossform` does not augment
anything. It splits the metric with complementary `D(w)`-orthogonal
projectors — `P = 1wᵀ/a` is idempotent and `D(w)P = wwᵀ/a` is symmetric — and
that differs in three ways. The decomposition is exact and additive rather
than a change of the analyzed set; it is inherited by every fixed linear query
of the form without recomputation; and both components are retained, so
nothing is chosen between.
The coherent term is the common mode under the *effective* metric weights, so
under a precision-weighted metric it is a precision-weighted common mode
rather than an arithmetic regional mean.

## Model-constrained low-rank geometry (Diedrichsen & Kriegeskorte, 2017; PCM; reduced-rank regression)

Diedrichsen and Kriegeskorte (2017) put encoding models, pattern-component
modelling (PCM) and RSA on one footing: every representational model is a
statement about the second moment of the condition patterns, `G = T C T'`
for a model factor `T` and a positive semidefinite `C`, and the three
methods differ in what they fix, marginalize, or estimate about `C`. PCM
fits `C` (or a parameterized family of it) by maximizing a Gaussian
likelihood of the observed patterns, which requires a full-rank noise model
and estimates the signal geometry as a within-sample second moment with the
noise variance subtracted. Feature-reweighted RSA (Kaniuth & Hebart, 2022;
Jozwik et al., 2016) learns nonnegative weights on the model's feature
dimensions before comparing RDMs, and reduced-rank regression and envelope
regression (Cook, 2018) constrain the rank of a coefficient matrix or of the
subspace on which the response varies.

crossform's descriptive model-coordinate reader (`model_basis()`, `model_geometry()`;
`design/model-coordinate-geometry-contract.md`) takes the second-moment
parameterization and changes where it is estimated. The model family enters
as a declared, centered, rank-revealing basis `Q` of its joint span with
`T = QR`, and the relation is *lowered* through it, `B~ = Q'B`; the ordinary
cross-partition estimator then yields the complete compressed form
`S_x = Q' G_x Q` at every measurement, which is a fixed linear query of the
crossvalidated form and inherits its noise unbiasedness. No noise floor is
subtracted, because no within-partition term is ever formed. The four model
structures of the second-moment framework -- one weight per model
(isotropic), one per model coordinate (diagonal, the feature-reweighting
case), one PSD block per model, one PSD rank-limited matrix over the joint
span (shared) -- are four readings of that one form: nonnegative quadratic
programs in the model-only Gram for the first two, exact block coordinate
descent for the third, and `latent_geometry()`'s rank-budgeted PSD
truncation for the fourth, which is the closed-form Frobenius minimizer
(Theorem 4 of the proposal the contract answers). Every fit returns `C`,
never a factor `W`, because `W` is identified only up to rotation; and every
fit is a latent projection reported apart from the signed
model-addressable energy `tr(S_x)` and the signed model-orthogonal energy
`tr((H - P) G_x)`, whose sum is the centered total exactly. The learned
readout of the shared structure is cross-fitted over partition edges under
the same `metric_training_policy("exclude_evaluation")` discipline that keeps
a learned neural metric honest, restoring a signed, zero-mean-under-noise
energy.

Two things are not claimed. A basis that spans every centered direction
reproduces the best rank-limited PSD approximation of the centered geometry
whatever models produced it (Theorem 3 of the proposal), so the shared fit
is refused at saturation and labelled near it. And a response envelope in
the sense of envelope regression is not identifiable from the
cross-generalized form: the form fixes the span of the signal second moment
but not a complementary subspace on which the response is invariant
(Theorem 5 of the proposal), so crossform reports spans and energies and
never an envelope.

The predictive workflow (`model_basis()` -> `fit_geometry()` ->
`score_geometry()`; `design/predictive-geometry-contract.md`) uses that same
lowering and retains each model's spectrum. Named simplex weights pool the
retained kernels. On their positive support `K = Q U D U' Q'`, the prediction
is `F = Q U [U' S U - penalty * D^-1]_{+,rank} U' Q'`. This is the exact
Frobenius geometry minimizer with inverse-kernel trace penalty and a PSD rank
budget, as established by completing the square in the contract. It is not
the solution to a generic RDM-weighted objective. Full-span positive-penalty
models remain meaningful predictive hypotheses; a zero-penalty full-span
fit is instead labelled a generic baseline.

The invariant prediction is the form `F`; its retained reduced factor is a
storage representation with arbitrary latent rotation. Independent scoring
freezes it and reads `2<F,G_test> - ||F||^2`. Conditional on training and an
independent unbiased test form, this estimates squared-error improvement
over zero. Its null expectation is `-||F||^2`, unlike the descriptive
reader's zero-mean cross-fitted projector energy. The code therefore keeps
these results as different classes and estimands. The contribution claimed
here is their composition through a reduced effect form and frozen readout,
with explicit fitting/evaluation dependencies; no new principle of spectral
truncation, kernel regularization or condition transfer is claimed.

## What is therefore claimed

Only this, at theorem strength: the voxel-axis partition is exact, its parts
are PSD, and every fixed linear query inherits it. Not claimed: the trace
statistic, cross-validation across runs, noise normalization, the
mean-versus-pattern distinction, or any inferential procedure.

## References

- Allefeld, C., & Haynes, J.-D. (2014). Searchlight-based multi-voxel pattern
  analysis of fMRI by cross-validated MANOVA. *NeuroImage*, 89, 345–357.
  <https://doi.org/10.1016/j.neuroimage.2013.11.043>
- Basti, A., Nili, H., Hauk, O., Marzetti, L., & Henson, R. N. (2020).
  Multi-dimensional connectivity: a conceptual and mathematical review.
  *NeuroImage*, 221, 117179.
  <https://doi.org/10.1016/j.neuroimage.2020.117179>
- Davis, T., LaRocque, K. F., Mumford, J. A., Norman, K. A., Wagner, A. D., &
  Poldrack, R. A. (2014). What do differences between multi-voxel and
  univariate analysis mean? How subject-, voxel-, and trial-level variance
  impact fMRI analysis. *NeuroImage*, 97, 271–283.
  <https://doi.org/10.1016/j.neuroimage.2014.04.037>
- Diedrichsen, J., Berlot, E., Mur, M., Schütt, H. H., Shahbazi, M., &
  Kriegeskorte, N. (2021). Comparing representational geometries using
  whitened unbiased-distance-matrix similarity. *Neurons, Behavior, Data
  Analysis, and Theory*, 5(3), 1–31. <https://doi.org/10.51628/001c.27664>
  (preprint: <https://doi.org/10.48550/arXiv.2007.02789>)
- Diedrichsen, J., & Kriegeskorte, N. (2017). Representational models: a
  common framework for understanding encoding, pattern-component, and
  representational-similarity analysis. *PLOS Computational Biology*, 13(4),
  e1005508. <https://doi.org/10.1371/journal.pcbi.1005508>
- Cook, R. D. (2018). *An introduction to envelopes: dimension reduction for
  efficient estimation in multivariate statistics*. Wiley.
  <https://doi.org/10.1002/9781119422976>
- Jozwik, K. M., Kriegeskorte, N., & Mur, M. (2016). Visual features as
  stepping stones toward semantics: explaining object similarity in IT and
  perception with non-negative least squares. *Neuropsychologia*, 83,
  201–226. <https://doi.org/10.1016/j.neuropsychologia.2015.10.023>
- Kaniuth, P., & Hebart, M. N. (2022). Feature-reweighted representational
  similarity analysis: a method for improving the fit between computational
  models, brains, and behavior. *NeuroImage*, 257, 119294.
  <https://doi.org/10.1016/j.neuroimage.2022.119294>
- Diedrichsen, J., Provost, S., & Zareamoghaddam, H. (2016). On the
  distribution of cross-validated Mahalanobis distances. *arXiv*:1607.01371.
  <https://doi.org/10.48550/arXiv.1607.01371>
- Kriegeskorte, N., Goebel, R., & Bandettini, P. (2006). Information-based
  functional brain mapping. *Proceedings of the National Academy of Sciences
  USA*, 103(10), 3863–3868. <https://doi.org/10.1073/pnas.0600244103>
- Nili, H., Wingfield, C., Walther, A., Su, L., Marslen-Wilson, W., &
  Kriegeskorte, N. (2014). A toolbox for representational similarity analysis.
  *PLoS Computational Biology*, 10(4), e1003553.
  <https://doi.org/10.1371/journal.pcbi.1003553>
- Schütt, H. H., Kipnis, A. D., Diedrichsen, J., & Kriegeskorte, N. (2023).
  Statistical inference on representational geometries. *eLife*, 12, e82566.
  <https://doi.org/10.7554/eLife.82566>
- Taylor, J. E., & Kriegeskorte, N. (2025). Framed RSA: representational
  comparisons that honor both geometry and population-mean response
  preferences. *bioRxiv*.
  <https://doi.org/10.1101/2025.07.10.664257>
- Walther, A., Nili, H., Ejaz, N., Alink, A., Kriegeskorte, N., &
  Diedrichsen, J. (2016). Reliability of dissimilarity measures for
  multi-voxel pattern analysis. *NeuroImage*, 137, 188–200.
  <https://doi.org/10.1016/j.neuroimage.2015.12.012>
