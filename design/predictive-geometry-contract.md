# Predictive geometry contract

Version: `predictive-geometry-v1`. Frozen for implementation on 2026-09-04.
Scope: G01–G22 and T01–T62 of the predictive-geometry work packages.
Migration: `design/predictive-geometry-migration.md`.

## Estimator and evaluation

A model declaration supplies fixed centered factors `T_i`. In a common
orthonormal union basis `Q`, write `T_i=Q R_i`. The neural estimation object
is the signed cross-partition form `G_x`, and its restriction is
`S_x=Q' G_x Q`. The estimator minimizes

```
0.5 * ||G_x - F||_F^2 + lambda * tr(K_alpha^+ F)
```

over PSD `F` of rank at most `s`, with range inside the positive support of
`K_alpha = sum(alpha_i T_i T_i')`. Weights are nonnegative and sum to one.
The positive support of `J_alpha = sum(alpha_i R_i R_i')` is `U D U'`.
The fitted reduced form is exactly

```
A_x = U [U' S_x U - lambda * D^-1]_(+,s) U'
F_x = Q A_x Q'
```

The range restriction is essential when a weight is zero. Subtracting the
pseudoinverse penalty in the entire union space is not this estimator.
No neural PSD projection precedes compression or the penalty shift.
The small eigensolver is exact for this Frobenius loss, not RDM regression
or arbitrary observation-weighted GLS.

For a frozen fitted form, independent evaluation returns

```
inner_product = <A_x, S_test,x>
prediction_norm_sq = ||A_x||_F^2
gain = 2 * inner_product - prediction_norm_sq
```

Conditional on all training and selection, unbiased `G_test` for `G*` gives
`E[gain] = ||G*||^2 - ||G*-F||^2`. Under zero signal this expectation is
`-||F||^2`, not zero. Scoring never selects a rank, clips test evidence or
discards a negative-gain mode. No generic standard errors are supplied.

## Input and coordinate contract

`model_basis()` accepts explicitly named kernel, feature and squared-Euclidean
RDM inputs. Kernel and RDM axes align independently by declared condition
names. Unlabelled matrices require a declared same-size condition axis.
Duplicate, missing and one-sided labels refuse. Effect coordinates must carry
shared units and scales. The declared centering is the arithmetic condition
mean, implemented in an explicit orthonormal frame of `1`-perpendicular.
Nonorthogonal effect rescaling changes the loss; no false invariance is claimed.

Kernels must be finite, square, symmetric and PSD under their declared repair
policy before centering. RDMs also require zero diagonals, nonnegative entries
and `distance="squared_euclidean"`. Features enter through their centered Gram.
All matrix products must remain finite; overflow is an input-scale error.

The default `negative_share=0` allows only numerical negative roots below
the relative spectral tolerance. Setting a larger limit is an explicit
positive-part model replacement, and input/centered negative mass and positive
truncation are recorded separately. A model with no positive centered
variation refuses; it is not the rank-zero response predictor.

Model ranks are positive integers or `NULL` (all numerical positive modes),
clamped to each model's available rank and recorded. `normalize="trace"`
scales each **retained** model kernel to unit trace; `normalize="none"`
preserves the declared model scale and records it. The default remains
`"none"` for the existing constructor, while predictive examples explicitly
declare trace normalization. Source representation, raw/retained scale,
normalization, PSD policy and centering enter the model signature.

The union factor `R` may be rectangular. Overlap is reported and admitted.
Readers that require unique redundant coefficients must refuse there.
Positive simplex weights are not estimated from neural data by the constructor.
All names must match the declared models exactly. Finite nonnegative weights
whose sum differs from one by more than `1e-12` refuse; weights are never
silently renormalized. A valid simplex vertex is admitted.

Response rank is a finite integer at least zero; oversized budgets clamp
to pooled support dimension and record the requested/effective budget.
`rank=0` is exactly the zero predictor with zero score. Penalty is a required
finite nonnegative scalar. With unit-trace model kernels it has the units of
the neural form. Gain has squared-geometry units. With effect unit `u` and
a dimensionless fixed neural metric these are `u^2` and `u^4`, respectively.

## Numerical and nonuniqueness policy

Default model/support relative eigenvalue tolerance is `1e-10`; callers may
declare a different finite tolerance in `(0,1)`. The factorization's relative
singular-value tolerance is separately recorded as that same declared numeric
value; eigenvalue and singular-value tests are not conflated. Root classes
and all losses follow the declared effective kernel rather than pretending
discarded numerical support was fitted.

Symmetry uses a scale-relative check, followed only by roundoff symmetrization.
Centering is asserted at `1e-10 * sqrt(n)`. Orthogonality, support residuals
and PSD checks use combined absolute/relative numerical tolerances. Stress
evidence reports conditioning and eigengaps; individual eigenvector accuracy
is never required at a vanishing eigengap.

A positive tie straddling an explicitly requested model or response rank cut
refuses with `ambiguous_rank_cut`; its gap is compared to the declared tie
tolerance times spectral scale. A rank-zero fit has no cut ambiguity.
All-positive modes retained, or clipping roots at zero, do not require choosing
a partial positive eigenspace. Ties wholly retained preserve the form, but
mode-wise evidence is reported for the tied group. Scientific prediction
identity is the form; exact coordinate/execution records may still differ
across floating-point eigendecompositions. No universal rotation-invariant
hash of arbitrary serialized factors is promised.

## Training support and matching targets

A training record contains actually used observation origins, pairing
products, source revisions, upstream transformations, selections and their
scope. The initial supported route uses disjoint subsets of one declared
origin manifest. Origin declarations come from ingestion or explicit caller
provenance; they are assumptions, not independence proved from arrays.
Known copies and renamings preserve origin ancestry. Distinct origins with
identical numeric values are not automatically dependent. Unknown external
lineage refuses independent scoring rather than accepting a boolean bypass.

The original plan may store unused partitions. Metadata validation and hashes
are not training on their outcomes. Numeric reads and actual training
dependencies are restricted to the training pairing. Parent plan identity
is distinct from the training footprint.

Evaluation requires both the independence needed for the neural cross-product
and the independence of those observations from every fitting/selection
dependency. Pairing endpoint counts alone are insufficient. Matching
targets require the same effect meaning, units, normalization, fixed metric,
component and spatial measurements. Training and evaluation pairings and
plan IDs may differ. Match measurement rows by stable IDs; missing/duplicate
IDs refuse. Origin scope and recorded assumptions remain visible in the score.

Execution admission refined by G08's package-route witnesses: the existing
complete symmetric-form executor requires an **undirected self-form pairing**.
Its complete coherent/configuration decomposition requires a **positive
definite fixed neural metric** (or implicit identity). Predictive preflight
refuses directed full-form plans, singular fixed metrics, whitened schedules,
and metrics marked `learned_frozen`; it does not silently symmetrize a directed
estimand or replace a learned metric's training policy. These are executor
admission limits, not limits on the spectral theorem. Reversing the spelling
of an undirected edge preserves its operator and weight. An explicitly supplied
effect-space signature is retained by the model declaration and checked before
lowering; matching labels and units alone cannot replace that declared meaning.

Fitting without sufficient provenance may produce a latent fit if that
limitation is explicitly recorded; it cannot acquire independent scoring
capability. No implementation claims to detect undeclared hidden dependence.
Fixed-effect geometry fitting does not admit learned/whitened metric schedules,
condition transfer, arbitrary spatial transport or GLS under this version.

## Records, decomposition, storage and selection

The public verbs are `fit_geometry()` and `score_geometry()` consuming the
model declaration. Fits and scores are separately sealed records, not
subclasses of raw signed geometry or latent descriptive projections. They
carry parent execution receipts and their own statistical records. Reduced
prediction factors use O(r*s) storage per measurement. Optional signed
training forms retain existing memory/block stores and bounded row reads.
Incomplete execution cannot emit a completed fit, score or reopenable artifact.

Predictive memory admission reserves model declarations, retained factors and
spectra, diagnostics, spectral scratch, bounded row buffers and serialization
overlap before delegating the remaining budget to the existing geometry
compiler. The result keeps both plans in `$workspace`; the numerical receipt
continues to describe the geometry execution alone. Admission can reduce the
requested row block, and refuses if one row cannot fit. This is a conservative
owned-workspace byte model, not a claim about total process RSS.

Only the existing sequential worker route is admitted. Complete fixed-metric
execution currently supports memory output only; predictive block storage
requires the implicit identity metric and otherwise refuses with
`fixed_metric_block_storage_not_admitted`. Failed or interrupted writes remove
only the newly created output directory. Completion RDS records are published
by a final rename after all rows and validators succeed.

One frozen prediction may read total, coherent and configuration evidence.
The linear evidence terms add. The squared prediction norm appears once in
the total score. Separately fitted nonlinear forms do not generally conserve.
The subspace residual is identified; an absolute whole-form residual or an
unbiased signal-energy denominator is not recovered from compression.

Cross-fitting uses only declared edges and their weights; each fold retains
its own fitted form and squared-norm offset. The aggregate is evidence for
the training procedure at those training sizes, not a final all-data refit.
Nested selection uses only inner validation gain with declared candidate
tuples and weights. Exact selection ties prefer lower rank, then larger
penalty, then a stable canonical candidate ID. Per-measurement versus global
selection is explicit. Outer evaluation cannot modify choices. Training-local
caches include source revisions, origins, folds, effect/metric/frame target,
component, model truncation/normalization, weights and penalty where relevant.

## Required evidence and statistical precision

T01–T62 are acceptance requirements in
`.planning/2026-09-04-predictive-geometry-test-specification.md`.
The base-R court cannot import production numerical helpers. Exact fixtures,
independent dense products, metamorphic laws, stress residuals, execution
parity and critical production-path mutants are all required.

For stochastic identity calibration use independent datasets as the MCSE unit
after dependent within-dataset folds are aggregated. Freeze seeds and arms.
Let `b=||G*||^2+E||G_test-G*||^2` from independently verified generator moments;
use the exact zero test if `b=0`. Otherwise practical margin is `delta=.02*b`.
Require `abs(mean(error)) <= 5*MCSE + numerical_tolerance` and
`abs(mean(error))+5*MCSE <= delta`. Freeze a Bonferroni familywise error budget
of `0.001`; the simultaneous normal multiplier is
`max(5, qnorm(1 - 0.001/(2 * number_of_arms)))` in place of five when larger.
An independently seeded pilot selects the fixed final replication count to
target half the equivalence margin; insufficient precision is not a pass.
Negative controls have separately prespecified detectable effects.

Heavy benchmarks and recertification run serially. Record actual source reads,
packed widths, temporary/retained bytes, incremental RSS, timing, configuration,
hardware and parity. Bind final evidence to source and oracle/harness versions.
No unexplained stale-evidence skip counts as passing certification.
