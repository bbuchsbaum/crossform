# effectagram

`effectagram` solves a task-fMRI problem: estimate where experimental effects
are reproducible, how their multivariate geometry changes across the brain,
and which local patterns generalize across independent runs without building a
different engine for voxels, ROIs, or searchlights. It turns partitioned
experimental effects into effect-space or measurement-space forms. When the
same named effects appear on both axes, an
effect form specializes to symmetric geometry. A voxel, searchlight, ROI, and
whole brain differ only in the spatial measurement applied to the same
relation. Contrast energy, squared-distance RDMs, regression RSA, node effects,
and cross-node coupling are partial views of one evidence pairing; they do not
trigger separate model-fitting pipelines.

The package is an experimental 0.1 implementation of the additive, fixed-
bilinear core. It runs sequentially, retains negative crossvalidated estimates,
and does not provide a classifier registry, raw fMRI preprocessing, or image
registration.

## Installation

The package is not yet published. From a local checkout, install it with:

```r
install.packages(".", repos = NULL, type = "source")
```

## One complete analysis

Suppose two independent runs already contain three condition effects over six
native features. The coordinates below stand in for voxel or surface
coordinates.

```r
library(effectagram)

effects <- effect_space(
  c("face", "house", "object"),
  basis_id = "condition-means:v1",
  units = "percent-signal"
)
domain <- abstract_domain(
  6,
  coordinates = cbind(x = 0:5, y = 0),
  feature_ids = paste0("feature", 1:6),
  id = "native:demo"
)

run1 <- rbind(
  face   = c(2, 1, 0.5, -0.5, 0, 1),
  house  = c(0, 1, 2, 1, 0.5, 0),
  object = c(0.5, 0.5, 1, 1.5, 1, 0.5)
)
run2 <- rbind(
  face   = c(2.2, 0.8, 0.4, -0.4, 0.1, 0.9),
  house  = c(0.1, 1.2, 1.8, 0.9, 0.6, 0.1),
  object = c(0.4, 0.6, 1.1, 1.4, 0.9, 0.6)
)

rel <- relation(list(run1 = run1, run2 = run2), effects = effects,
  domain = domain)
at <- compile_frame(searchlights(radius = 1.01, normalization = "local"), domain)
plan <- plan_geometry(rel, at = at, over = cross_partitions(rel))
plan
g <- geometry(plan)
g
```

`plan` is the reusable query-first object. This example explicitly materializes
it because the next sections ask several questions of the complete result. `g`
contains total and coherent geometry; configuration geometry is their exact
difference. No searchlight-specific model was fit.

```r
face_house <- contrast(g, c(face = 1, house = -1, object = 0))

round(cbind(
  signed = face_house$signed,
  coherent = face_house$coherent,
  configuration = face_house$configuration,
  total = face_house$total
), 3)
#>      signed coherent configuration total
#> [1,]  0.925    0.850         1.250 2.100
#> [2,]  0.133    0.017         2.083 2.100
#> [3,] -1.017    1.033         0.317 1.350
#> [4,] -1.117    1.244         0.189 1.433
#> [5,] -0.333    0.111         0.889 1.000
#> [6,]  0.200    0.037         0.488 0.525
```

The columns answer different questions:

| View | Meaning |
|---|---|
| `signed` | Weighted regional contrast, averaged across pairing endpoints |
| `coherent` | Reproducible energy in the weighted regional mean |
| `configuration` | Reproducible energy in spatial departures from that mean |
| `total` | Exact sum of coherent and configuration energy |

Cross-generalized components can be negative because they are unbiased around
zero. `coherence_fraction` is therefore reported only at locations where the
observed components form a nonnegative partition.

The same `g` supports other questions without recomputation:

```r
distances <- rdm(g)
roots <- geometry_spectrum(g)

category_rdm <- matrix(
  c(0, 0, 1,
    0, 0, 1,
    1, 1, 0),
  3, 3, byrow = TRUE,
  dimnames = list(rownames(run1), rownames(run1))
)
model_fit <- rsa(g, models = list(category = category_rdm))
```

RDMs are squared geometry distances and may also be negative under
cross-generalization. RSA compiles the distance transform and regression into
one linear query. `geometry_spectrum()` preserves negative roots instead of
silently converting the estimate into a positive-semidefinite description.

## Beyond self geometry

A general effect form has independent ordered left and right effect spaces.
`pair_query()` binds a fixed matrix to those exact axes. Higher-level
constructors including `match_control()` and `pair_lm_query()` compile
matched-control and pair-space regression coefficients to the same matrix
representation, including unequal encoding and retrieval axes. They report the
row and column balance of the final operator; they do not promise that later
external weighting will preserve baseline invariance.

The executable form compiler currently admits fixed inner-product evidence.
It records partition reduction order explicitly; where a public workflow
accepts `reduce_partitions()` or `aggregate_first()`, those values name
different estimands. Generic correlation, ranking, and Fisher edge pipelines
remain outside the public API until they have an executable and scale-qualified
path.

Two relations with different neural-domain identities require an explicit
fixed `measurement_bridge()`. Its left and right legs enter one named
`measurement_space()`; equal dimensions alone never imply alignment. The
package does not learn bridges or perform hyperalignment, pair-edge independence
corrections, or group inference.

## Read node effects and coupling from one form

`measurement_form()` closes the experimental side with `by` and leaves only
the neural edges named by `between` open. Its `left` and `right` arguments are
reserved for experimental-neural relation sides. Neural endpoints use
`edge_frame(from, to, ...)`, and every requested edge is explicit - there is no
automatic searchlight-by-searchlight graph.

```r
nodes <- measurement_frame(
  list(
    seed = matrix(c(1, 0, 0, 0, 0, 0), 1),
    target = matrix(c(0, 1, 0, 0, 0, 0), 1)
  ),
  domain = domain,
  id = "prespecified-nodes:v1"
)
between <- edge_frame(
  from = c("seed", "seed", "target", "target"),
  to = c("seed", "target", "seed", "target"),
  frame = nodes
)
by <- variation_query(
  (diag(3) - matrix(1 / 3, 3, 3)) / 2,
  effects,
  sampling_axis = "condition",
  construction = "joint_covariance"
)
within_partition_products <- pairing(
  rel$partitions,
  rel$partitions,
  directed = TRUE,
  self_pairs = "allow_biased",
  independence = "not_independent"
)
measured <- measurement_form(
  left = rel,
  between = between,
  by = by,
  over = within_partition_products
)
connectivity(measured, view = "correlation")$values
```

The raw form is always algebraically available through `effect_coupling()`.
Other views validate stronger interpretation contracts:

| View | Additional requirement |
|---|---|
| `covariance_coupling()` | repeated variation, effective rank above one, and a coherent joint covariance |
| `connectivity(..., "correlation")` | valid scalar self-variances and oriented rank-one measurements |
| `canonical_coupling()` | valid multivariate self-covariances and explicit ridge regularization |
| `geometry_alignment()` | nonzero self-geometries; static CKA/RV-like form-entry normalization |
| Gaussian information | an explicit joint Gaussian model and reported information units |

A rank-one contrast remains valid effect coupling but is rejected as normalized
connectivity: one effect direction does not provide repeated variation.
Crossvalidated self-blocks are not silently used as covariance denominators.
`geometry_alignment()` names static Gram-matrix alignment; it does not redefine
dynamic informational-connectivity methods based on learned discriminability
series.

For a frame-complete block form, `reconstruct_evidence()` uses Parseval or
canonical-dual reconstruction when both frames have full column rank. A
rank-deficient frame returns only an explicitly projected operator. Missing
edges, incompatible bases, unacceptable conditioning, and insufficient dense
workspace are rejected before a lossless claim.

This measurement/connectivity surface is explicitly a small-node dense path in
version 0.1. `measurement_frame()` rejects estimated dense frame and leg
payloads above 256 MiB before converting an additive frame, and tomography has
a 512 MiB default reconstruction budget. Brain-scale local work belongs on the
support-local geometry plan until matrix-free measurement frames pass their
scale gates.

See `vignette("evidence-pairing", package = "effectagram")` for signed scalar
correlation, multivariate coupling, coherent/configuration edge channels,
cross-domain forms, and Parseval reconstruction.

## Change spatial scope without changing the engine

All ordinary additive scopes compile to one sparse measurement frame:

```r
point_frame  <- compile_frame(voxels(), domain)
local_frame  <- compile_frame(searchlights(radius = 1.01), domain)
region_frame <- compile_frame(regions(c("A", "A", "A", "B", "B", "B")), domain)
global_frame <- compile_frame(whole_brain(), domain)
```

Local normalization makes every row sum to one, so measurements are local
averages. For the feature-additive metric path, including the current
Euclidean default `K = I`, conservative normalization makes every feature
distribute unit mass across overlapping measurements and local total
geometries sum to global geometry. This conservation is a compiled capability,
not a property inherited by a non-diagonal local metric. The normalization and
metric identity are part of the result metadata.

Frame weights and a support-local metric compose by the exact congruence

```text
K_effective = D(sqrt(w)) K D(sqrt(w)).
```

It recovers the existing additive metric exactly when `K = I` and preserves
positive semidefiniteness. Native-diagonal metrics retain the shared additive
contraction; non-diagonal metrics require a support-streamed pair contraction.
The latter is scientifically valid but does not inherit the claim that all
overlapping searchlight work has been collapsed into one sparse multiply.

## Plan once, then query what you need

The primary path is query-first. `plan_geometry()` validates the relation,
frame, pairing, source capabilities, and compute policy without reading neural
blocks or allocating packed geometry:

```r
plan <- plan_geometry(rel, at, cross_partitions(rel))
query <- bilinear_query(tcrossprod(c(1, -1, 0)), effects = effects)
only_total <- evaluate_geometry(plan, query = query)
```

The query is fused before spatial contraction, so a one-coordinate question
does not create the full packed effect form. The result is an `effect_view`
marked `query_only`; it never pretends that full geometry was materialized.

When several later views genuinely need every coordinate, materialize it
explicitly:

```r
g <- geometry(plan)
```

`geometry(rel, at, over)` remains a compatibility form, but it compiles and
executes this same plan. Complete materialization is size-preflighted and can
use memory or block-backed storage.

## What the package owns

`effectagram` accepts a lazy relation, not necessarily beta files. A supplied
linear extractor can map raw response rows into effects block by block:

```r
extractor <- lm_extractor(
  design = X,
  effects = C,
  observation_whitener = L
)
rel <- relation(response_runs, extract = extractor, domain = domain)
```

When later work needs a noise metric or within-participant calibration, retain
the statistical channel explicitly:

```r
condition <- factor(rep(c("face", "house", "object"), each = 12))
X <- model.matrix(~ 0 + condition)
colnames(X) <- levels(condition)
C <- diag(3)
rownames(C) <- colnames(X)

fit <- lm_relation_fit(
  response_runs,
  design = X,
  effects = C,
  effect_names = colnames(X),
  sampling_unit = "trial",
  observation_whitener = L,
  domain = domain
)
```

`fit$relation` is still the pure experimental-neural relation. The separate
error channel supplies lazy whitened residual blocks, residual degrees of
freedom, and the unscaled effect-coordinate covariance under a named
separable-GLM assumption. A precomputed `relation()` remains valid, but it does
not claim noise-estimation capabilities that its inputs cannot recover.

For a noise-normalized contrast, compile the error channel into a frozen
on-demand metric schedule and query it directly:

```r
crossnobis_plan <- plan_crossnobis(
  fit,
  at = compile_frame(searchlights(radius = 6), domain),
  over = cross_partitions(fit$relation),
  metric = shrinkage_precision(0.1),
  training = metric_training_policy("exclude_evaluation")
)
crossnobis_map <- crossnobis(crossnobis_plan, condition_contrast)
crossnobis_map
head(as.data.frame(crossnobis_map))
```

`plan_crossnobis()` stores the residual sufficient-statistic graph, estimator
identity, training assignment, regularization, and source revisions. It never
stores one covariance or precision matrix per searchlight and fold. Each local
solve is derived only when its support is evaluated; negative crossvalidated
estimates remain negative. The alternative
`all_partitions_residual_orthogonality` policy requires a written justification
and produces a different plan identity. Neither policy currently exports
confidence intervals or LD-t calibration, because uncertainty in the estimated
metric has not yet been propagated.

The training-only policy is the conservative default. In a 500-replication
Gaussian-GLM validation, training-only and explicitly justified all-run
residual reuse were equivalent within a predeclared +/-0.005 margin. When the
same analysis deliberately omitted AR(1) temporal structure, their signal
estimates differed by -0.00936 (95% Monte Carlo interval -0.01140 to -0.00732),
almost entirely because they learned different metric targets. The
metric-training policy is therefore scientific identity, not a speed setting.

### Ask for analytic RDM covariance when the metric is fixed

An all-partition-pairs RDM is one estimator. Its partition-pair contributions
are not independent replicates: two pairs that share a run also share that
run's estimation error. Do not use the spread of pair values divided by the
square root of the pair count as a standard error.

When the geometry plan uses one common fixed metric and `fit` has equal
partition error structures, effectagram implements the complete sampling
covariance from Diedrichsen, Provost, and Zareamoghaddam (2016), Eq. 13. The
calibration target is explicit because the signal-dependent term differs under
the fixed zero null and a plug-in estimate:

```r
known_precision <- noise_precision(
  diag(6), domain,
  covariance = diag(6),
  provenance = list(source = "fixed-before-evaluation")
)
fixed_plan <- plan_geometry(
  fit$relation,
  at = compile_frame(whole_brain(), domain),
  over = cross_partitions(fit$relation),
  metric = known_precision
)

distance_estimate <- rdm(fixed_plan)
distance_covariance <- rdm_sampling_covariance(
  fixed_plan, fit, target = "plugin"
)
distance_variance <- sampling_covariance(distance_covariance)
distance_se <- sqrt(distance_variance)
```

Use `target = "null"` when the scientific calculation is explicitly under the
zero-distance null. The returned object is the covariance law, not an automatic
confidence interval or test. It is stored as exact factors: diagonal
variances, selected entries, covariance actions, quadratic forms, and linear
transports are evaluated without constructing the dense distance-by-distance
matrix. `operation = "materialize"` is explicit and size-preflighted.

This path requires `lm_relation_fit()`. A relation made from precomputed beta
matrices has valid point estimates but no recoverable residual covariance or
effect-coordinate covariance, so effectagram refuses analytic within-person
uncertainty and explains the missing entry-path capability. Subject-level
resampling estimates population uncertainty; it does not reconstruct those
discarded first-level quantities. The current analytic law also refuses
learned metrics, unequal partition error structures, and cross-location field
inference rather than silently extending Eq. 13 beyond its assumptions.

In the committed 10,000-repetition validation under the declared fixed-metric,
equal-partition matrix-normal model, full-covariance relative error was 1.16%
under signal and 1.96% under the null; analytic coordinatewise coverage ranged
from 94.75% to 95.29%. The naive pair-spread interval covered only 72.84% to
73.97% under signal. These results validate that named model and target; they
are not a general inference claim. See the
[sampling-covariance contract](evidence-sampling-contract.md) and the
[validation artifact](benchmark-results/sampling-covariance-validation-summary.csv).

The named brain-scale gate uses 52,416 active voxels, 6.1 mm searchlights, and
one dense learned local shrinkage metric per searchlight. It completed the
fitted-relation, residual-statistics, crossnobis, and NeuroVol-mapping path in
123.5 seconds with 0.95 GB incremental peak RSS on the recorded development
machine. No pair-frame, pair-atom field, or node-by-edge factor table was
materialized. These figures qualify this exact fixture; they are not a promise
for arbitrary support sizes, fold counts, hardware, or metric learners.
This is a scale qualification, not a statistical-recovery result: the fixture
has 30 training residual degrees of freedom and supports up to 33 features, so
the declared shrinkage estimator supplies the invertible local metric. The
separate 500-replication study above carries the estimator-policy evidence.

The package does not construct `X`, choose an HRF, estimate temporal whitening,
correct motion, or register images. Function-backed sources must declare a
strong immutable revision with `source_capabilities()`; matrix-backed sources
receive one automatically. Execution validates the compute policy, source
capabilities, frame, pairing, query, and conservative memory plan before reading
neural blocks. `compute_policy(workspace_bytes = ...)` optionally imposes a
hard budget on effectagram-owned live workspace; process RSS is recorded
separately rather than folded into that limit.

For neuroim2 volumes, `neuroim2_volume_domain()` and
`neuroim2_searchlights()` provide a narrow conditional bridge to stable
full-volume indices. They do not import neuroim2 ROI iteration, data
extraction, result objects, or parallel state.

Use `as_neurovol(values, mask, domain)` to place one compact result value at
each exact full-volume index carried by `neuroim2_volume_domain()`. It preserves
the source NeuroSpace and performs no interpolation or smoothing.

Version 0.1 supports one owned sequential worker. For compiler-proven
feature-additive metrics it removes repeated searchlight work through the
additive contraction. Dense local metrics are support-streamed and carry their
larger pair and factorization cost explicitly. Process-level parallelism is
deferred until it can pass the package's bounded-memory and deterministic-
reduction gates.

See `vignette("introduction", package = "effectagram")` for the continuous
workflow and the generated function reference for exact argument contracts.
