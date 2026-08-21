# Coming from rMVPA

> **This is a mapping guide, not a benchmark.** Nothing here compares
> the two packages on speed, accuracy, or scientific quality. It answers
> one question: if you already think in
> [rMVPA](https://github.com/bbuchsbaum/rMVPA)’s vocabulary, what is the
> corresponding `crossform` object, and what is deliberately missing?
> The one place the two packages have been measured against each other
> is a numerical-parity check, reported at the end with its caveats.
>
> The rMVPA calls below are shown but not executed — rMVPA is not on
> CRAN, so this vignette cannot depend on it. The `crossform` calls are
> executed.

## The one structural difference

rMVPA is organized around a **model** that is *run over* a set of
locations. You build a dataset, a design, and a model specification,
then hand the whole thing to `run_searchlight()` or `run_regional()`,
which returns a performance map.

`crossform` is organized around a **plan** that is *queried*. You
declare a fitted relation, where to measure it, and what must
generalize; that produces a plan. Contrasts, RDMs, and RSA are then
different reads of the same plan, and they do not refit anything:

``` r

# rMVPA: one run per question
run_searchlight(model_a, radius = 8)
run_searchlight(model_b, radius = 8)

# crossform: one plan, many questions
plan <- plan_geometry(rel, at = frame, over = pairing)
contrast_energy(plan, weights)
rdm(plan)
rsa(plan, models = list(category = model_rdm))
```

The practical consequence is where cross-validation lives. In rMVPA the
blocking variable is consumed by a cross-validation scheme that folds a
*model*. In `crossform` the generalization axis is bound into the
*estimand’s identity* by
[`cross_partitions()`](https://bbuchsbaum.github.io/crossform/reference/cross_partitions.md),
before any data is read, and it travels with the plan.

## Concept map

| rMVPA | crossform | Notes |
|----|----|----|
| `mvpa_dataset(train_data, mask =)` | `neuroim2_volume_domain(mask)` + `relation(list_of_run_matrices, ...)` | rMVPA holds the image; `crossform` splits *where* (domain) from *what was estimated* (relation). |
| `mvpa_design(df, y_train = ~y, block_var = ~block)` | `effect_space(condition_names)` for the conditions; run structure is the names of the list passed to [`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md) | The same table drives both routes in “Getting your own data in” below: the condition column becomes the effect space, the run column becomes the partitions. |
| `blocked_cross_validation(block_var)` | `cross_partitions(rel, independence = "independent", generalizes_over = "run")` | rMVPA folds train/test; `crossform` forms products between *different* partitions. Not the same object: see below. |
| `run_searchlight(spec, radius = 8)` | `compile_frame(searchlights(radius = 8), domain)` or `neuroim2_searchlights(mask, radius = 8)`, then [`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md) | Radius is millimeters in both. |
| `run_regional(spec, region_mask)` | `compile_frame(regions(labels), domain)` | [`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md) takes one label per **feature** (a vector), not a `NeuroVol`. |
| `rsa_design(~ model, data = list(model = D), block_var =)` | the `models =` list of [`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md), plus [`cross_partitions()`](https://bbuchsbaum.github.io/crossform/reference/cross_partitions.md) for the blocking | `crossform` aligns model RDMs by condition **name**. |
| `rsa_model(dset, des, distmethod =, regtype =)` | `rsa(plan, models = list(...))` | Two separate arguments: `distmethod` picks the distance, `regtype` the second-order statistic. See “distance conventions” below. |
| `contrast_rsa_model(dset, msreve_des, output_metric = "beta_delta")` | `contrast_energy(plan, weights)` | The nearest analogue, not an equivalent. rMVPA regresses the crossvalidated second-moment matrix on several competing contrast RDMs at once and multiplies each signed coefficient by the *center voxel’s* projection onto that contrast, giving one signed map per contrast. [`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md) scores one contrast on its own, returns that contrast’s crossvalidated energy for the whole measurement, and splits it into coherent and configuration parts. |
| `vector_rsa_model(dset, des, distfun = cordist())` | no equivalent | Trial-level; `crossform` aggregates to condition level. |
| `compute_crossvalidated_means_sl(..., estimation_method = "crossnobis")` + `compute_crossnobis_distances_sl()` | `rdm(plan)` under [`cross_partitions()`](https://bbuchsbaum.github.io/crossform/reference/cross_partitions.md), or [`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md) / [`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md) | This is the pairing where the two packages have been checked against each other. |
| `mvpa_model(model = "sda_notune", ...)`, classification accuracy | no equivalent | `crossform` has no classifier. |
| train/test folds | cross-partition pairing | A fold count is an execution detail in `crossform`; the *generalization axis* is part of plan identity. |
| performance metrics (`compute_performance`) | the views: [`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md), [`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md), [`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) | Energies and signed squared distances, not accuracies. |
| `feature_selector("FTest", "top_k", 100)` | no equivalent | See absences below. |
| `run_permutation_searchlight(...)`, `permutation_control()` | no equivalent | See absences below. |
| `future::plan(multisession)` (rMVPA uses furrr/future.apply) | [`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md), sequential only | `workers > 1` is refused, not silently ignored. |
| `save_results(x, dir, ...)` | [`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md) + your own writer | `crossform` writes no files. |

### The cross-validation mapping deserves a second look

`blocked_cross_validation(block_var)` and
`cross_partitions(rel, generalizes_over = "run")` both encode “runs are
the independent unit”, but they are not interchangeable objects.

A leave-one-run-out fold scheme produces training means that **overlap
across folds**. Cross-partition pairing forms products between two
*disjoint* run estimates, so the noise in the two factors is independent
and cancels in expectation. That is exactly why `crossform`’s estimates
are unbiased and may be negative. The Haxby exemplar records a concrete
instance of the difference mattering: an older rMVPA build derived
crossnobis folds from leave-one-run-out training means, whose overlap
made the distances positively biased, and the exemplar script asserts
fold semantics rather than trusting a version number.

## Side-by-side minimal example

The same question — where does an animate-versus-inanimate pattern
reproduce across runs? — written both ways.

### rMVPA

``` r

library(rMVPA)

dset <- mvpa_dataset(train_data = beta_vec, mask = vt_mask)

des <- mvpa_design(
  data.frame(
    y = factor(rep(conditions, times = n_run), levels = conditions),
    block = run_of
  ),
  y_train = ~ y,
  block_var = ~ block
)

rsa_des <- rsa_design(~ animacy, data = list(animacy = model_rdm),
                      block_var = run_of)
mspec <- rsa_model(dset, rsa_des, distmethod = "pearson", regtype = "pearson")
res <- run_searchlight(mspec, radius = 6, method = "standard")

# res$results$animacy is a NeuroVol of second-order correlations.
```

### crossform

``` r

example <- example_fmri_effects()

plan <- plan_geometry(
  example$fit$relation,
  at = example$frame,
  over = cross_partitions(
    example$fit$relation,
    independence = "independent",
    generalizes_over = "run"
  )
)

effect <- contrast_energy(plan, example$contrast)
category <- rsa(plan, models = list(category = example$model_rdm))

c(
  measurements = length(effect$total),
  peak = which.max(effect$total),
  rsa_columns = ncol(category$coefficients)
)
#> measurements         peak  rsa_columns 
#>          280          144            2
```

Two things differ beyond syntax. `plan` is reusable: `rdm(plan)` and any
number of further contrasts cost no refit. And
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md)
returns OLS coefficients with an intercept column by default —
[`ncol()`](https://rdrr.io/r/base/nrow.html) above is 2 for one model
because the intercept is retained. rMVPA reaches the same family with
`rsa_model(..., regtype = "lm")`; its default `regtype = "pearson"`
scores the match as a correlation between the two vectorized RDMs
instead.

``` r

distances <- rdm(plan)
dim(distances$values)
#> [1] 280   6
colnames(category$coefficients)
#> [1] "(Intercept)" "category"
```

## Getting your own data in: from a trial-level beta series

The example above starts from a fitted relation. If what you have is
what rMVPA wants — a trial-level beta series plus a design table with a
condition column and a run column — this section is the bridge.
Everything here runs.

``` r

set.seed(2026)
conditions <- c("face", "body", "house", "tool")
runs <- paste0("run", 1:4)

# One row per trial, five trials per condition per run. This is the
# `mvpa_design()` table, and `beta_series` is the matrix behind the
# `mvpa_dataset()` NeuroVec once the mask has been applied.
trial_design <- expand.grid(
  trial = 1:5, condition = conditions, run = runs, stringsAsFactors = FALSE
)
trial_design$condition <- factor(trial_design$condition, levels = conditions)
n_trial <- nrow(trial_design)
n_voxel <- 40L

beta_series <- matrix(rnorm(n_trial * n_voxel), n_trial, n_voxel)
planted <- 1:8
animate <- trial_design$condition %in% c("face", "body")
beta_series[animate, planted] <- beta_series[animate, planted] +
  rep(rep(c(1.5, -1.5), length.out = length(planted)), each = sum(animate))

dim(beta_series)
#> [1] 80 40
head(trial_design, 3)
#>   trial condition  run
#> 1     1      face run1
#> 2     2      face run1
#> 3     3      face run1
```

[`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md)
wants one **condition-by-voxel matrix per run**, not a trial list,
because the second-order geometry is defined on condition estimates and
the runs are the partitions that must generalize. There are two ways to
get there, and they are not equivalent in what they earn you.

### Route 1: average the betas within condition, per run

The direct translation. [`split()`](https://rdrr.io/r/base/split.html)
on the run column, then a column mean per condition.

``` r

domain <- abstract_domain(n_voxel, id = "beta-series-example")
by_run <- split(seq_len(n_trial), trial_design$run)

run_matrices <- lapply(by_run, function(rows) {
  t(vapply(conditions, function(condition) {
    colMeans(beta_series[rows[trial_design$condition[rows] == condition], ,
      drop = FALSE])
  }, numeric(n_voxel)))
})

averaged <- relation(run_matrices, effects = effect_space(conditions),
                     domain = domain)
plan_averaged <- plan_geometry(
  averaged,
  at = compile_frame(whole_brain(), domain),
  over = cross_partitions(averaged, independence = "independent",
                          generalizes_over = "run")
)
animacy <- c(face = 0.5, body = 0.5, house = -0.5, tool = -0.5)
energy_averaged <- contrast_energy(plan_averaged, animacy)
energy_averaged$total
#> [1] 0.4209631
```

Row names matter:
[`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md)
aligns partitions by effect **name**, so the matrices must be named by
condition. `t(vapply(...))` above produces that.

### Route 2: hand the trial betas over as observations

Better, when you can.
[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md)
takes the trials themselves plus a per-run design matrix and fits the
condition estimates inside `crossform`, which keeps the residual error
channel that the analytic sampling law needs. The design is the
condition indicator; `effects` says each condition *is* its own
coefficient.

``` r

designs <- lapply(by_run, function(rows) {
  X <- stats::model.matrix(~ 0 + condition,
                           data = trial_design[rows, , drop = FALSE])
  colnames(X) <- conditions
  X
})
sources <- lapply(by_run, function(rows) beta_series[rows, , drop = FALSE])
targets <- diag(length(conditions))
dimnames(targets) <- list(conditions, conditions)

fit <- lm_relation_fit(
  sources, design = designs, effects = targets,
  effect_names = effect_space(conditions),
  sampling_unit = "trial", domain = domain
)
plan_fitted <- plan_geometry(
  fit$relation,
  at = compile_frame(whole_brain(), domain),
  over = cross_partitions(fit$relation, independence = "independent",
                          generalizes_over = "run")
)
energy_fitted <- contrast_energy(plan_fitted, animacy)
```

With this balanced indicator design the two routes estimate the same
condition means, so they report the same energy:

``` r

c(route_1 = energy_averaged$total,
  route_2 = energy_fitted$total,
  max_difference = max(abs(energy_averaged$total - energy_fitted$total)))
#>        route_1        route_2 max_difference 
#>   4.209631e-01   4.209631e-01   1.110223e-16
```

What route 2 buys is the error channel. Only the fitted relation can be
asked for a within-measurement standard error:

``` r

sampling_capabilities(plan_fitted, fit)
#> <effect_sampling_capabilities>
#>   analytic sampling law: available 
#>   metric: fixed | partitions: equal | error channel: relation_fit

round(sqrt(sampling_covariance(rdm_sampling_covariance(
  plan_fitted, fit, target = "null", at = 1L
))), 4)
#>  face - body face - house  face - tool body - house  body - tool house - tool 
#>       0.0247       0.0247       0.0247       0.0247       0.0247       0.0247
```

Ask the averaged relation the same question and you get a refusal naming
exactly what is missing, rather than a number:

``` r

catch_refusal(
  rdm_sampling_covariance(plan_averaged, averaged, target = "null", at = 1L)
)
#> <effect_capability_refusal>
#>   capability:  sampling_covariance
#>   namespace:   evidence_sampling
#>   reasons:
#>     - missing_error_channel
#>   remedies:
#>     - Refit raw observations with `lm_relation_fit()`.
#>   state:       refused; no partial result was produced
```

Two practical notes. Runs need not have equal trial counts — `design`
accepts a named list with one matrix per partition, which is why it is
built inside [`lapply()`](https://rdrr.io/r/base/lapply.html) above. And
if your betas came from a GLM whose residuals are temporally
autocorrelated, pass `observation_whitener =` as well; see
[`vignette("from-observations", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/from-observations.md)
for what that argument is asserting.

## Distance and second-order conventions are not the same

This is the most likely source of surprise, and it is not a bug on
either side.

The two arguments of `rsa_model()` do separate jobs, and only one of
them is where the packages part company.

- **The distance.** `crossform`’s
  [`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md)
  reports signed **crossvalidated squared Euclidean** distances,
  `d_ij = G_ii + G_jj - 2 G_ij` (squared Mahalanobis under a fixed
  metric). They can be negative, and are not truncated. rMVPA’s
  `rsa_model()` takes `distmethod` in `{"pearson", "spearman"}` and
  nothing else, so its neural RDM is always a **correlation distance**.
  This is the difference that cannot be bridged. (rMVPA’s crossnobis
  distances live elsewhere, on
  `contrast_rsa_model(estimation_method = "crossnobis")` and the
  `compute_crossnobis_distances_sl()` helper — which is why the parity
  check below uses those and not `rsa_model()`.)
- **The second-order statistic.** This is a *separate* argument,
  `regtype`, taking `"pearson"`, `"spearman"`, `"lm"`, or `"rfit"`.
  `crossform`’s
  [`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) is
  OLS with an intercept, which is the same family as `regtype = "lm"`.
  rMVPA’s default is `regtype = "pearson"`, a correlation between the
  two vectorized RDMs; set `regtype = "lm"` and this half of the
  comparison lines up.

`crossform` cannot express correlation distance at all: dividing each
pattern by its own norm is not a fixed linear operation on the patterns,
so it is outside the bilinear core. Per-sphere mean-centering *is*
expressible as a metric; the unit-norm denominator is not. The boundary
and the conditions under which it could be crossed are recorded in the
[correlation-distance
policy](https://bbuchsbaum.github.io/crossform/articles/correlation-distance-policy.md)
([`vignette("correlation-distance-policy", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/correlation-distance-policy.md)
offline).

If you want the two packages to compute the same thing, the estimand
both can express is the condition-level crossvalidated squared Euclidean
(crossnobis) distance — which is the pairing used for the parity check
below.

## What crossform does not do (yet)

Each absence below was checked against the package’s actual exports
rather than recalled.

- **Classification and decoding.** There is no classifier, no accuracy,
  no AUC, no confusion matrix, no `MVPAModels` registry. `crossform`
  computes second-order geometry only.
- **Group inference.** No second-level model, no random-effects map, no
  across-subject test.
  [`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
  is a *within*-measurement, within-participant covariance law and says
  so.
- **Feature selection.** No `FTest`/`catscore` ranking, no `top_k`
  cutoff. The frame decides which features each measurement reads, and
  it is declared in advance.
- **Permutation and bootstrap.** No label shuffling, no null
  distributions, no FDR correction.
- **Parallel backends.** Execution is sequential. This is refused rather
  than silently downgraded:

``` r

compute_policy(workers = 4)
#> Error:
#> ! crossform 0.1 owns no process pool, so `workers` must be 1; received `4`. Sequential execution is a capability boundary, not a performance default: an executor that spawns workers has to pass memory and determinism gates first.
```

- **Preprocessing and I/O.** No registration, no smoothing, no file
  reading or writing.
  [`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md)
  hands you a `NeuroVol`; saving it is yours.

Two things rMVPA users may expect that `crossform` *does* have, in a
different shape: a crossnobis path
([`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md),
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md),
with
[`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)
and
[`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md)
for the whitening), and explicit refusals when an uncertainty claim is
not supported by the fit
([`sampling_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/sampling_capabilities.md),
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)).

## The one measured comparison

The [Haxby 2001
exemplar](https://github.com/bbuchsbaum/crossform/tree/main/exemplars/haxby2001)
runs both packages on eight stimulus categories estimated in twelve
runs, over 577 ventral-temporal searchlights at an 11.25 mm radius. On
the matched estimand — `crossform`’s
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) over
a cross-run geometry against rMVPA’s
`compute_crossvalidated_means_sl(estimation_method = "crossnobis")`
followed by `compute_crossnobis_distances_sl()` — the two agree to a
**maximum absolute difference of 8.88e-16** at every one of the 577
centers. Sphere membership was verified identical rather than assumed.

Read that number with the exemplar’s own caveats, which are recorded
there in full:

1.  **It is not `rsa_model`.** The exemplar originally named `rsa_model`
    for the rMVPA arm and could not use it: its distances are
    correlation distances, which `crossform` does not compute.
    `vector_rsa_model` is trial-level. `contrast_rsa_model` does build
    the cross-validated second-moment matrix internally, but what it
    returns is a signed per-contrast contribution evaluated at each
    searchlight’s center voxel, never the condition-by-condition
    distances themselves. The crossnobis helper is the only route, and
    `compute_crossnobis_distances_sl()` is not exported, so the script
    calls it with `:::`.
2.  **It is version-sensitive.** The comparison requires an rMVPA build
    whose crossnobis fold estimates are per-run condition means. The
    script asserts that semantics and stops with an explanatory error
    otherwise, rather than checking a version string.
3.  **The originally specified estimand was not achievable.** A
    correlation-distance RDM scored by Spearman cannot be expressed by
    `crossform` at all: `rsa_model()`’s `distmethod` offers only
    `"pearson"` and `"spearman"`, and `crossform` has no correlation
    distance to meet either one. The second-order statistic is not the
    obstacle — that is a separate argument, `regtype`, whose `"lm"`
    setting is the same OLS family as `crossform`’s
    [`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md).
    It is the distance that cannot be shared. The matched comparison is
    therefore a substitution, and the exemplar says so.

What the number demonstrates is numerical parity on a shared estimand.
It demonstrates nothing about speed, and nothing about
correlation-distance RSA.

## Where to go next

- [`vignette("introduction", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/introduction.md)
  — what the energies, RDMs, and RSA coefficients mean.
- [`vignette("neuroim2-data", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/neuroim2-data.md)
  — the input and output adapters, if your data is already in neuroim2
  objects.
- [Reading contrast energies, RDMs, and
  uncertainty](https://bbuchsbaum.github.io/crossform/articles/interpreting-results.md)
  ([`vignette("interpreting-results", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/interpreting-results.md))
  — how to read the numbers, and the three interpretive traps.
- [What is novel in
  crossform?](https://bbuchsbaum.github.io/crossform/articles/novelty.md)
  ([`vignette("novelty", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/novelty.md))
  — what is architecture and what is established statistics.
