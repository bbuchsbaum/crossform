# effectagram

`effectagram` turns partitioned experimental effects into local or global
cross-generalized geometry. A voxel, searchlight, ROI, and whole brain differ
only in the spatial frame used to measure the same relation. Contrast energy,
squared-distance RDMs, and regression RSA are views of the resulting geometry;
they do not trigger separate model-fitting pipelines.

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
g <- geometry(rel, at = at, over = cross_partitions(rel))
```

`g` is the durable result. It contains complete total and coherent geometry;
configuration geometry is their exact difference. No searchlight-specific
model was fit.

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

## Change spatial scope without changing the engine

All ordinary additive scopes compile to one sparse measurement frame:

```r
point_frame  <- compile_frame(voxels(), domain)
local_frame  <- compile_frame(searchlights(radius = 1.01), domain)
region_frame <- compile_frame(regions(c("A", "A", "A", "B", "B", "B")), domain)
global_frame <- compile_frame(whole_brain(), domain)
```

Local normalization makes every row sum to one, so measurements are local
averages. Conservative normalization makes every feature distribute unit mass
across overlapping measurements, so local total geometries sum to global
geometry. The normalization is part of the result metadata.

## Materialize geometry or evaluate one query

Use `geometry()` when several later views are likely. It always returns a
complete `effect_geometry`, whether stored in memory or in a block-backed
directory.

Use `evaluate_geometry()` when only a fixed bilinear question is needed:

```r
query <- bilinear_query(tcrossprod(c(1, -1, 0)), effects = effects)
only_total <- evaluate_geometry(rel, at, cross_partitions(rel), query)
```

This returns an `effect_view` marked `query_only`; it never pretends that the
full geometry was materialized.

## What the package owns

`effectagram` accepts a lazy relation, not necessarily beta files. A supplied
linear extractor can map raw response rows into effects block by block:

```r
extractor <- lm_extractor(design = X, effects = C, whiten = L)
rel <- relation(response_runs, extract = extractor, domain = domain)
```

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

Version 0.1 supports one owned sequential worker. It first removes repeated
searchlight work through the additive contraction; process-level parallelism is
deferred until it can pass the package's bounded-memory and deterministic-
reduction gates.

See `vignette("introduction", package = "effectagram")` for the continuous
workflow and the generated function reference for exact argument contracts.
