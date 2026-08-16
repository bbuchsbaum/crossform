# crossform

`crossform` answers one question about task-fMRI data: **where does an
experimental effect — or a whole representational geometry — reproduce across
independent runs?**

![Effects, frame, and pairing compile into one geometry plan, which four view functions then read.](man/figures/crossform-overview.svg)

You declare the condition effects you fitted, the spatial measurements you want
answers at, and which runs must generalize. `crossform` compiles those once into
a **geometry plan**, which records the quantity you are estimating before any
brain data is read. Contrasts, RDMs, RSA, and crossnobis distances are queries
against that one plan rather than four separate pipelines.

A **measurement**, used consistently below, is a spatial unit: a searchlight, a
region, a voxel, or the whole brain. Every result has one row per measurement.

## Install

The repository is not yet public and the package is on neither CRAN nor
R-universe, so today the route is a local checkout:

```r
install.packages(".", repos = NULL, type = "source")
```

Once published, `remotes::install_github("bbuchsbaum/crossform")` and the
R-universe repository `https://bbuchsbaum.r-universe.dev` become the supported
routes.

## Why it is different

- **One fitted geometry, many questions.** `contrast_energy()`, `rdm()`,
  `rsa()` and `crossnobis()` read the same compiled plan, so a contrast map and
  an RDM at the same searchlights come from the same estimates by
  construction, and a second question costs a query rather than a refit. With
  100 conditions, ask for three of the 4,950 pairs and the rest are never
  computed.
- **Generalization is part of what you estimated, not a fold count.**
  `cross_partitions(rel, generalizes_over = "run")` is bound into the plan's
  identity. Crossvalidating over runs and over sessions are different
  scientific quantities, and the plan says which one you ran; changing block
  size or storage changes only the execution record.
- **Negative estimates are kept; uncertainty is refused when it is not
  earned.** Crossvalidated energies and distances are centered on zero under
  the null, so about half go negative where there is no effect, and clipping
  them would bias the estimator upward. Conversely, if your relation is beta
  matrices alone, `rdm_sampling_covariance()` refuses to return standard errors
  and names the missing input instead of substituting the spread across run
  pairs.

## A complete first analysis

```r
library(crossform)

example <- example_fmri_effects()   # 4 conditions, 4 runs, 245 searchlights

plan <- plan_geometry(
  example$fit$relation,
  at   = example$frame,
  over = cross_partitions(
    example$fit$relation,
    independence     = "independent",
    generalizes_over = "run"
  )
)
plan
#> <effect_geometry_plan>
#>   effects:      4 x 4
#>   measurements: 245
#>   features:     245
#>   metric:       implicit identity
#>   generalizes:  6 partition pairs over run, endpoints independent
#>   lowering:     additive_contraction
#>   dense payload: 68.9 Kb
#>   state:        query-first; call materialize_geometry(plan) to materialize
```

Save `plan` with the analysis record: it names the conditions, domain, frame,
run pairs, metric, and units. Now ask where the animate-versus-inanimate
contrast reproduces:

```r
effect <- contrast_energy(plan, example$contrast)
peak   <- which.max(effect$total)

round(c(
  searchlight   = peak,
  signed        = effect$signed[peak],
  total         = effect$total[peak],
  coherent      = effect$coherent[peak],
  configuration = effect$configuration[peak]
), 3)
#>   searchlight        signed         total      coherent configuration
#>       166.000         1.381         4.223         1.906         2.318

peak %in% example$truth$signal_measurements
#> [1] TRUE
```

`effect$signed` is the ordinary GLM contrast at that searchlight. The three
energy columns are what is new. `total` is the contrast energy that reproduces
across runs; `coherent` is the part carried by the searchlight's average
pattern — the familiar univariate story; `configuration` is the reproducible
spatial pattern beyond that mean. Because `coherent + configuration = total`,
you can say how much of a searchlight's effect survives once its regional mean
is accounted for. Here configuration carries most of it, as the planted
multivariate signal should.

Change the question without recompiling the geometry:

```r
distances <- rdm(plan)
category  <- rsa(plan, models = list(category = example$model_rdm))

dim(distances$values)       #> [1] 245   6   — searchlights by condition pairs
dim(category$coefficients)  #> [1] 245   2   — intercept and category model

# Or ask for exactly the pairs you mean; the rest is never materialized.
rdm(plan, pairs = rbind(c("face", "house")))
```

For conditions `i` and `j`, `rdm()` reports
$d_{ij} = G_{ii} + G_{jj} - 2G_{ij}$: under cross-run pairing a signed
crossvalidated squared Euclidean distance, and squared Mahalanobis
(crossnobis) once you supply a fixed neural metric. It is deliberately not
`1 - Pearson correlation`; see the
[correlation-distance policy](https://bbuchsbaum.github.io/crossform/articles/correlation-distance-policy.html).

## Uncertainty, only when it is earned

```r
round(sqrt(sampling_covariance(
  rdm_sampling_covariance(plan, example$fit, target = "null", at = peak)
)), 4)
#>  face - body face - house  face - tool body - house  body - tool house - tool
#>       0.0145       0.0145       0.0145       0.0145       0.0145       0.0145
```

These are within-searchlight standard errors under a declared equal-partition,
fixed-metric, separable error model — not a random-field model, a confidence
interval, or group inference. They need the residual channel that
`lm_relation_fit()` and `estimate_relation()` retain; beta matrices alone cannot
supply it. `sampling_capabilities(plan, example$fit)` answers the admission
question before you provoke it, and `catch_refusal(expr)` returns the unmet
requirements and their remedies as data. The
[failure gallery](https://bbuchsbaum.github.io/crossform/articles/failure-gallery.html)
shows six realistic errors the package guards against.

## Bring your own data

If you already have one condition-by-voxel beta matrix per run, this is the
whole path — it runs as written:

```r
set.seed(1)
conditions <- c("face", "house", "cat", "chair")

mask <- array(TRUE, c(8, 8, 6))                        # substitute your mask
dom  <- volume_domain(mask, spacing = c(3, 3, 3))

# Substitute your betas: one condition-by-voxel matrix per run, rows named.
betas <- lapply(1:4, function(run) matrix(
  rnorm(4 * sum(mask)), 4, sum(mask), dimnames = list(conditions, NULL)
))
names(betas) <- paste0("run", 1:4)

rel      <- relation(betas, effects = effect_space(conditions), domain = dom)
frame    <- compile_frame(searchlights(radius = 6), dom)
own_plan <- plan_geometry(rel, at = frame, over = cross_partitions(
  rel, independence = "independent", generalizes_over = "run"
))

own_effect <- contrast_energy(own_plan, c(face = 1, house = -1, cat = 0, chair = 0))
```

Swap `searchlights(radius = 6)` for `regions(labels)`, `voxelwise()`, or
`whole_brain()`; swap `volume_domain()` for `abstract_domain(n_features)` when
the features are not a volume. For `neuroim2` users,
`neuroim2_volume_domain()` and `neuroim2_searchlights()` preserve stable
full-volume indices and `as_neurovol()` maps compact results back without
interpolation or smoothing — see `vignette("neuroim2-data")`.

If you have scan-by-feature responses rather than betas, start one step earlier
with `study()`, `plan_relation()`, and `estimate_relation()`.

## What is actually novel?

Not RSA, crossnobis, or the analytic covariance formula — those are
established. The contribution is architectural: contrast energy,
squared-distance RDMs, fixed linear RSA, and ordered cross-axis hypotheses are
queries against one typed cross-partition form, and the plan keeps that
scientific identity separate from how it was executed.
[What is novel in crossform?](https://bbuchsbaum.github.io/crossform/articles/novelty.html)
is the ledger separating what is demonstrated from what is still gated.

## Where next

- [Introduction](https://bbuchsbaum.github.io/crossform/articles/introduction.html) — the entry point: a ready relation to contrasts, RDMs, RSA, and standard errors.
- [Coming from rMVPA](https://bbuchsbaum.github.io/crossform/articles/from-rmvpa.html) — how the familiar objects map onto this vocabulary.
- [Fit condition effects from observations](https://bbuchsbaum.github.io/crossform/articles/from-observations.html) — scan responses, events, confounds, censoring.
- [neuroim2 data](https://bbuchsbaum.github.io/crossform/articles/neuroim2-data.html) — `NeuroVol` and `NeuroVec` in, brain maps out.
- [Evidence pairing](https://bbuchsbaum.github.io/crossform/articles/evidence-pairing.html) — results that relate two regions, and the contracts each connectivity view requires.
- [What is novel in crossform?](https://bbuchsbaum.github.io/crossform/articles/novelty.html) · [Failure gallery](https://bbuchsbaum.github.io/crossform/articles/failure-gallery.html) · [Correlation-distance policy](https://bbuchsbaum.github.io/crossform/articles/correlation-distance-policy.html)
- [`design/`](design/) (contracts) · [`exemplars/haxby2001/`](exemplars/haxby2001/) (public-data comparison) · [`benchmarks/`](benchmarks/) (runtime and memory records)

## Status and scope

`crossform` is experimental; exported names may still change. It runs
sequentially, and it does not preprocess fMRI, register images, build
hemodynamic response models, or perform group inference.

**Demonstrated.** On the [Haxby 2001 exemplar](exemplars/haxby2001/) (12
conditions, 577 ventral-temporal searchlights), `crossform` agrees with an
independent reference loop to `1.33e-15` and with `rMVPA` to `8.88e-16` on the
matched crossvalidated squared-Euclidean/crossnobis estimand; refitting the raw
responses to retain the error channel reproduces the point RDM to `4.44e-16`.
At 100 conditions over 1,080 searchlights, one hundred selected pairs run in
0.27 s and the fused full RDM in 5.16 s against 13.40 s for
materialize-then-project — a ratio of 0.39, with a `4.4e-16` oracle.

**Not demonstrated.** Those results show numerical parity and an integrated
uncertainty path. They do **not** show a matched-estimator speed advantage,
and say nothing about correlation-distance RSA. `rsatoolbox` parity, a real
rectangular cross-axis exemplar, and an operational conservation example remain
to be earned; map-scale runtime and storage claims are qualified under
[`benchmarks/`](benchmarks/).
