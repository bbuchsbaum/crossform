# Matched multiscale interpretability simulation contract

Status: normative simulation contract, PE-C

Contract version: `matched-interpretability-v1`

Date: 2026-08-21

This contract defines the generated truth used to ask whether equal total
effects with different spatial organization can be distinguished. The
generator is
`benchmarks/matched-interpretability/00-mixture-generator.R`. It uses base R
only and returns truth separately from any Crossform estimator output.

## 1. Fixed-total mixture

For `p` features, define the normalized regional-mean direction

\[
u_{\mathrm{mean}}=\mathbf 1/\sqrt p
\]

and a normalized configuration direction satisfying

\[
\mathbf 1^\top u_{\mathrm{config}}=0,
\quad
u_{\mathrm{mean}}^\top u_{\mathrm{config}}=0,
\quad
\lVert u_{\mathrm{mean}}\rVert_2
=\lVert u_{\mathrm{config}}\rVert_2=1.
\]

For fixed total magnitude `T > 0` and `theta` in `[0, pi/2]`, the planted
contrast pattern is

\[
b(\theta)=\sqrt T\{
  \cos(\theta)u_{\mathrm{mean}}
  +\sin(\theta)u_{\mathrm{config}}
\}.
\]

Under the unnormalized whole-territory identity metric, the coherent projector
is `P = 11' / p` and its exact complement is `I - P`. Therefore

\[
T_{\mathrm{total}}=b^\top b=T,
\qquad
T_{\mathrm{coherent}}=b^\top P b=T\cos^2(\theta),
\qquad
T_{\mathrm{configuration}}
=b^\top(I-P)b=T\sin^2(\theta).
\]

The endpoints are pure regional mean at `theta = 0` and pure configuration at
`theta = pi/2`. Interior values are prespecified mixtures, not post hoc labels.
The relation-ready effect matrix stores `condition_a = b/2` and
`condition_b = -b/2`, so the named contrast `(1, -1)` recovers `b` exactly.

## 2. Generator record

Every generated object records:

- the schema and basis identifiers;
- the mean and configuration vectors and their Gram matrix;
- the direction, frame, and coherent-projector normalization;
- `T`, `theta` in radians and degrees, and the seed;
- expected total, coherent, and configuration values and shares; and
- the rank-one neural spectrum.

The strict default rejects nonunit, nonorthogonal, or noncanonical directions.
Three explicit relaxation flags admit such inputs only for negative fixtures.
When a relaxation is active, the generator computes truth from the actual
projectors and records `strict_basis = FALSE`; it does not report the
`cos^2/sin^2` law for a basis that does not satisfy its assumptions.

## 3. Evidence boundary

## 3. Conservative multiscale extension

PE-C2 places the same fixed-total patterns on a 17-feature line and evaluates
four overlapping conservative searchlight frames at radii `0.01`, `1.01`,
`2.01`, and `4.01`. Equal alpha weights stack them into one frame family. The
fixture deliberately includes smaller boundary supports; each compiled weight
matrix remains sparse, and every nonpoint scale overlaps features across
nodes.

For conservative frame weights `W_s`, pattern `b`, and row mass `a_x`, the
independent expected scale summaries are

\[
T_s=\alpha_s\sum_{x,v}W_{s,xv}b_v^2=\alpha_s T,
\qquad
C_s=\alpha_s\sum_x\frac{(\sum_vW_{s,xv}b_v)^2}{a_x},
\qquad
Q_s=T_s-C_s.
\]

The coherent share `C_s/T_s` is independent of alpha. Three prespecified
scenarios use the C1 path: pure broad coherent (`theta = 0`), a 50/50
broad-plus-alternating mixture (`theta = pi/4`), and pure fine alternating
configuration (`theta = pi/2`). The acceptance thresholds are fixed before
noisy simulations:

- broad coherent share is one at every scale;
- fine configuration share exceeds `0.85` at the first nonpoint radius and
  peaks at the largest declared radius;
- mixed configuration share first exceeds `0.45` at radius `2.01`; and
- scale ordering, frame-family ordering, and family relabeling do not change
  values after joining on numeric scale.

These thresholds belong to this exact line-domain fixture. They are recovery
targets for later paired-noise simulations, not universal statements about
spatial frequency or searchlight radius.

## 4. Paired end-to-end observations

PE-C3 generates balanced two-condition observations in four independent run
partitions. The prespecified sample sizes are 6 and 24 trials per condition.
The SNR grid is `0`, `0.2`, and `0.8`, labeled null, low power, and recoverable.
SNR means contrast-pattern RMS divided by average feature-noise standard
deviation. Every noise covariance is normalized to unit average marginal
variance, so an SNR label has the same scale across regimes.

One seeded standard-normal array supplies all cells. Three fixed transforms
produce:

- iid Gaussian noise;
- diagonal heteroskedastic noise with feature SDs spanning 0.55 to 1.45 before
  RMS normalization; and
- spatially correlated noise with AR-like correlation `0.6^distance`.

Within a replicate, broad, mixed, and fine scenarios reuse the same design,
partitions, and transformed noise. SNR levels reuse the same noise. The smaller
sample is a condition-stratified prefix of the larger sample. This pairing
prevents an organization from receiving an easier draw and supports paired
comparisons in later children.

The generator stores four layers separately: immutable base noise and designs,
latent truth, observation cells, and estimator output returned only by
`matched_observation_fit()`. Its compact manifest identifies every cell by
scenario, noise regime, sample size, SNR, partition count, feature count, seed,
design, base-noise source, and truth row.

## 5. Conventional-summary ambiguity

PE-C4 compares the spectrum with two conventional quantities computed directly
from fitted partition contrast patterns. Neither comparator is weakened or
redefined for the demonstration.

1. **Regional activation** is the mean feature value of the signed contrast
   pattern after averaging partitions. It keeps first-moment sign and units.
2. **Aggregate multivariate magnitude** is the uniform mean over all unordered
   independent partition pairs of their Euclidean contrast-pattern inner
   product. In the noiseless fixture it is the global crossvalidated squared
   contrast magnitude. For two conditions it is also the sole nonredundant
   squared-Euclidean RDM entry.

For scenarios `a` and `b`, a comparator is operationally ambiguous when its
absolute truth difference is at most `1e-12` while the maximum absolute
difference between their configuration-share spectra is at least `0.2`. The
pure broad and 50/50 mixed scenarios meet this criterion for aggregate
multivariate magnitude by construction: both have the same total, but their
scale organization differs. Pure fine and mixed do as well.

The signed regional activation is a negative control. It correctly separates
the broad and mixed scenarios because their mean-mode coefficients differ. The
claim is therefore not that conventional summaries are always insufficient.
It is that an aggregate magnitude cannot identify how a fixed total is divided
between coherent and configuration structure.

The independent comparator uses base matrix arithmetic only. Tests compare it
with the package's whole-territory `signed` and `total` readings on fitted
Gaussian, heteroskedastic, and spatially correlated cells.

## 6. Canonical figure

The canonical reader-facing figure is generated by
`benchmarks/matched-interpretability/04-canonical-figure.R`. It uses 24 paired
Gaussian replications at 24 trials per condition and SNR 0.8. Six panels share
one scenario color mapping:

1. planted ground-truth patterns;
2. conventional signed regional activation;
3. aggregate crossvalidated multivariate magnitude;
4. coherent magnitude at the widest radius;
5. configuration magnitude at the widest radius; and
6. coherent share over all four scales.

Panels 4 and 5 use one y-axis. Analytic truth is marked separately from Monte
Carlo means. Isolated estimates use empirical 2.5%--97.5% whiskers; the
spectrum uses the same interval as a band. The underlying generated table is
versioned beside the PNG, so equal totals and scenario ordering are numerical
assertions rather than visual impressions.

Rendering tests check panel inventory, scenario and scale coverage, common
component axes, PNG dimensions, and file substance. They do not use a pixel
hash, font-specific coordinate comparison, or a subjective image score.

## 7. Certification

`benchmarks/matched-interpretability/05-certify.R` is the full certification
job. It runs 48 paired seeds over all 54 scenario/noise/sample-size/SNR cells.
Balanced-design OLS and frame contraction are evaluated by direct formulas;
the ordinary test suite separately exercises public-package smoke cells. The
full job records:

- the versioned parameter and threshold manifest;
- analytic global and scale truth;
- grouped estimates and empirical intervals;
- recovery, false-separation, conservation, ambiguity, and negative-control
  metrics; and
- source, contract, evidence-ledger, numeric-artifact, and article checksums.

The PNG is intentionally not byte-hashed. Its scientific input CSV is hashed,
while the plot itself is guarded by its six-panel contract, dimensions, and
semantic data checks. This avoids treating platform font rasterization as a
scientific result.

The predeclared gates are: conservation and recomposition at `1e-10`; paired
null false separation at `1e-12`; recoverable aggregate relative bias and
coherence-share MAE at `0.12`; broad-greater-than-mixed-greater-than-fine
ordering in at least `80%` of recoverable replications; paired aggregate
equivalence within `15%` of the true total; broad-minus-mixed wide-scale share
separation of at least `0.35`; and a regional-activation negative-control gap
of at least `0.15`.

The estimated coherence share is the signed crossvalidated coherent component
divided by the signed crossvalidated total. Crossvalidated component estimates
are unbiased but are not constrained to be nonnegative in finite samples.
Certification therefore does not truncate components or discard a replication
merely because either component crosses zero; doing so would turn the recovery
court into a positivity-selected analysis. A zero numerical denominator
remains undefined.

The default test court runs a two-seed full-grid smoke reconstruction. Setting
`CROSSFORM_FULL_MATCHED_SIMULATION=true` additionally rebuilds the 48-seed
numeric artifacts and compares them with the committed certification record.
The dedicated workflow runs that full gate on demand and weekly.

The checksum manifest binds this contract and evidence-ledger claim `CF-H10`
to the scenario definitions and thresholds. Changing any of them without an
explicit ledger update and certification rebuild fails the artifact court.

## 8. Evidence boundary

The generator and multiscale extension prove no broad empirical interpretation.
PE-C1 establishes analytic fixed-total truth; PE-C2 establishes its noiseless
scale expectations and public-route alignment; PE-C3 establishes paired,
reproducible observation generation and execution. The completed 48-seed court
adds comparator ambiguity, visualization, and prespecified recovery evidence,
so evidence-ledger claim `CF-H10` is a `matched_simulation`. Its scope remains
the declared synthetic line-domain regimes. It is not empirical evidence and
does not establish population coverage, power, or a neuroscience result.
