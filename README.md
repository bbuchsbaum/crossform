# crossform

Ask a task-fMRI dataset four questions from one fit, and find out what *kind*
of effect you found.

![Effects, frame, and pairing compile into one geometry plan, which four view functions then read.](man/figures/crossform-overview.svg)

You declare the condition effects, the measurements you want answers at, and
which runs must generalize. That compiles once into a **geometry plan**, which
records the quantity you are estimating before any brain data is read.
Contrasts, RDMs, RSA, and crossnobis distances are queries against it.

A **measurement** is a spatial unit: a searchlight, a region, a voxel, or the
whole brain. Every result has one row per measurement.

## Five lines to a first result

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

effect <- contrast_energy(plan, example$contrast)
plot(effect,
  highlight       = example$truth$signal_measurements,
  highlight_label = "planted signal"
)
```

![Left: every searchlight in the coherent/configuration plane. Right: total contrast energy along the searchlight index, planted searchlights filled.](man/figures/readme-contrast-energy.png)

The fixture plants an animate-versus-inanimate pattern in 57 of the 245
searchlights. The analysis is not told which; `highlight` only colors the
answer afterwards.

The right panel is the whole map. The 55 largest energies in it are all planted
searchlights, and the largest energy anywhere else is `0.03`, against `4.22` at
the peak. The rest lies *on* zero rather than above it, 55% of it below,
because crossvalidated energy is centered at zero under the null and negatives
are kept rather than clipped. That band is the noise floor itself. The left
panel splits the same energy in two, and the planted searchlights climb the
vertical axis.

Every object prints as itself.

```r
effect
#> <effect_contrast_view>
#>   measurements: 245
#>   contrast:     face 0.5, body 0.5, house -0.5, tool -0.5
#>  measurement    signed   coherent configuration     total coherence_fraction
#>            1 -0.001094 -0.0010841    -0.0059598 -0.007044                 NA
#>            2  0.015774 -0.0014452    -0.0002452 -0.001690                 NA
#>            3  0.029483 -0.0015089    -0.0036478 -0.005157                 NA
#>            4  0.002694 -0.0007510    -0.0028487 -0.003600                 NA
#>            5 -0.028748 -0.0004546     0.0028346  0.002380                 NA
#>            6  0.004236 -0.0017620    -0.0032531 -0.005015                 NA
#>   ... 239 more measurements

peak <- which.max(effect$total)

round(c(
  searchlight   = peak,
  signed        = effect$signed[peak],
  total         = effect$total[peak],
  coherent      = effect$coherent[peak],
  configuration = effect$configuration[peak]
), 3)
#>   searchlight        signed         total      coherent configuration
#>       166.000         1.381         4.223         1.906         2.318
```

`coherence_fraction` is `NA` in those rows because the two components do not
form a nonnegative partition there, and a share of a negative quantity is not a
share of anything.

One line in that block deserves a second look. `cross_partitions()` pairs every
partition with every other one, once. `generalizes_over` names the axis those
partitions *are* — runs, sessions, tasks — and binds it into the plan's
identity, because cross-run and cross-session reproduction are different
quantities at the same fold count; if your partitions are sessions, pass
`"session"`. `independence = "independent"` declares that the partitions'
estimation errors are independent, which separate runs or sessions with
separate noise satisfy; omit it and you keep a point estimate but earn neither
cross-generalized nor analytic-uncertainty capabilities. See
`?cross_partitions`, and `?pairing` for nested or directed designs.

## Now change the question. No refit.

`plan` holds an identity, not a matrix.

```r
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

Save it with the analysis record; it names the conditions, domain, frame, run
pairs, metric, and units. Nothing is materialized until a view asks.

```r
distances <- rdm(plan)
plot(distances, measurement = peak)
```

![Signed squared distances at searchlight 166: the two animate conditions and the two inanimate conditions sit near zero apart, every animate-inanimate pair above four.](man/figures/readme-rdm.png)

At searchlight 166, face and body are `0.016` apart and house and tool
`-0.003`, while every animate-versus-inanimate pair is above `4`. Nothing asked
for that block structure.

```r
category <- rsa(plan, models = list(category = example$model_rdm))
dim(category$coefficients)
#> [1] 245   2
```

Linear RSA against a model RDM, from the same plan: an intercept and the model.

```r
one_pair <- rdm(plan, pairs = rbind(c("face", "house")))
dim(one_pair$values)
#> [1] 245   1
```

One of the six pairs; the other five are never materialized. With 100
conditions, three of the 4,950 pairs cost three pairs.

`rdm()` reports $d_{ij} = G_{ii} + G_{jj} - 2G_{ij}$: under cross-run pairing a
signed crossvalidated squared Euclidean distance, and squared Mahalanobis
(crossnobis) once you supply a fixed neural metric. It is deliberately not
`1 - Pearson correlation`; see the
[correlation-distance policy](https://bbuchsbaum.github.io/crossform/articles/correlation-distance-policy.html).

## What kind of effect is it?

That left panel is the part that is hard to get elsewhere. A searchlight's
reproducible energy splits in two:

- **`coherent`** is the energy carried by the searchlight's frame-weighted
  average pattern, which is the familiar univariate story;
- **`configuration`** is the reproducible spatial pattern beyond that mean;
- `coherent + configuration = total`, exactly, at every measurement, which is
  why the dashed guides are lines of constant total.

So you can say how much of an effect survives once its regional mean is
accounted for, without a second analysis and without discarding either part.
Here 96% of the planted searchlights carry more configuration than coherent
energy. The effect is a pattern, and the picture says so.

`signed`, back in that first table, is a first moment: the ordinary contrast of
the region-average pattern, in your effect's units and sign. The energies are
second moments, squared and crossvalidated. So `signed` sees only the coherent
half.

```r
round(as.data.frame(effect)[123, c("signed", "coherent", "configuration", "total")], 3)
#>     signed coherent configuration total
#> 123 -0.037        0         1.993 1.993
```

Searchlight 123 is planted. Its signed contrast is smaller than the largest
signed value anywhere outside the signal (`0.189`), so a signed map cannot tell
it from noise, yet it reproduces `1.993` of pattern energy. Seventeen of the 57
planted searchlights are like this.

Most toolboxes ask you to choose first: remove the regional mean and the
coherent half is gone, keep it and the two are summed into one number.
[Reading the results](https://bbuchsbaum.github.io/crossform/articles/interpreting-results.html)
covers the traps in the coherent share.

## Bring your own data

If you already have one condition-by-voxel beta matrix per run, this is the
whole path. It runs as written.

```r
set.seed(1)
conditions <- c("face", "house", "cat", "chair")

mask  <- array(TRUE, c(8, 8, 6))                       # substitute your mask
dom   <- volume_domain(mask, spacing = c(3, 3, 3))
frame <- compile_frame(searchlights(radius = 4), dom)

# Substitute your betas: one condition-by-voxel matrix per run, rows named.
# Twelve voxels here carry the same face-minus-house pattern in every run,
# alternating in sign so the regional mean stays near zero.
planted <- array(FALSE, dim(mask))
planted[3:4, 3:5, 3:4] <- TRUE
pattern <- rep(c(3, -3), length.out = sum(planted))

betas <- lapply(1:4, function(run) {
  b <- matrix(rnorm(4 * sum(mask)), 4, sum(mask), dimnames = list(conditions, NULL))
  b["face",  planted[mask]] <- b["face",  planted[mask]]  + pattern
  b["house", planted[mask]] <- b["house", planted[mask]] - pattern
  b
})
names(betas) <- paste0("run", 1:4)

rel      <- relation(betas, effects = effect_space(conditions), domain = dom)
own_plan <- plan_geometry(rel, at = frame, over = cross_partitions(
  rel, independence = "independent", generalizes_over = "run"
))

own_effect <- contrast_energy(own_plan, c(face = 1, house = -1, cat = 0, chair = 0))
plot(own_effect)
```

The same story on "your" data. Did the searchlights that actually touch the
planted voxels win?

```r
touching <- which(Matrix::rowSums(frame$weights[, planted[mask], drop = FALSE]) > 0)
ranked   <- order(own_effect$total, decreasing = TRUE)

c(
  searchlights          = nrow(frame$index),
  touching_signal       = length(touching),
  top_energies_touching = sum(ranked[seq_along(touching)] %in% touching),
  peak_energy           = round(max(own_effect$total), 2),
  largest_elsewhere     = round(max(own_effect$total[-touching]), 2)
)
#>          searchlights       touching_signal top_energies_touching
#>                384.00                 44.00                 44.00
#>           peak_energy     largest_elsewhere
#>                 27.48                  1.29
```

The 44 largest energies in the map are exactly the 44 searchlights that touch
the planted voxels. Every one of them is configuration-dominated, and 48% of
the rest fall below zero.

Swap `searchlights(radius = 4)` for `regions(labels)`, `voxelwise()`, or
`whole_brain()`, and `volume_domain()` for `abstract_domain(n_features)` when
the features are not a volume. `neuroim2` users have
`neuroim2_volume_domain()`, `neuroim2_searchlights()`, and `as_neurovol()`,
which keep stable full-volume indices and map results back without
interpolation; see `vignette("neuroim2-data")`. Starting from scan-by-feature
responses instead? Use `study()`, `plan_relation()`, and `estimate_relation()`.

### Reading a map without ground truth

`highlight` is optional — the call above omits it — and the picture still
works. In the right panel the band around zero *is* the noise floor, because
crossvalidated energy is centered there when nothing reproduces; candidates are
the measurements standing clear of it. In the left panel each candidate then
declares its kind: high on the vertical axis is a reproducible pattern, far
along the horizontal axis is a regional mean.

What the picture does not give you is a threshold or a p-value. `crossform`
does no spatial and no group inference (see [Status and scope](#status-and-scope)),
and a within-measurement standard error needs the residual channel, so fit with
`lm_relation_fit()` or the
[from-observations](https://bbuchsbaum.github.io/crossform/articles/from-observations.html)
route rather than importing betas.
[Reading a map without ground truth](https://bbuchsbaum.github.io/crossform/articles/interpreting-results.html#reading-a-map-without-ground-truth)
is the longer version.

## Uncertainty, only when it is earned

`example_fmri_effects()` was built by `lm_relation_fit()`, so it kept its
residual channel and can be asked for standard errors.

```r
round(sqrt(sampling_covariance(
  rdm_sampling_covariance(plan, example$fit, target = "null", at = peak)
)), 4)
#>  face - body face - house  face - tool body - house  body - tool house - tool
#>       0.0145       0.0145       0.0145       0.0145       0.0145       0.0145
```

These are within-searchlight standard errors under a declared equal-partition,
fixed-metric, separable error model, not a random-field model, a confidence
interval, or group inference.

Ask the same of the betas-only relation above and the answer is an object that
names what is missing and how to earn it.

```r
catch_refusal(
  rdm_sampling_covariance(own_plan, rel, target = "null", at = 1L)
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

Beta matrices hold no residuals and no residual degrees of freedom, and the six
run pairs share run estimates, so their spread is not a standard error either.
`sampling_capabilities(plan, example$fit)` answers the admission question
before you provoke it, and the
[failure gallery](https://bbuchsbaum.github.io/crossform/articles/failure-gallery.html)
shows six realistic errors the package guards against.

## Why it is different

- **One fitted geometry, many questions.** A contrast map and an RDM at the
  same searchlights come from the same estimates by construction, and a second
  question costs a query rather than a refit.
- **Generalization is part of what you estimated, not a fold count.** The
  declared axis is bound into the plan's identity, so the plan says which
  quantity you ran; block size and storage change only the execution record.
- **The novelty is architectural.** Not RSA, crossnobis, or the analytic
  covariance formula, which are established, but that contrast energy,
  squared-distance RDMs, fixed linear RSA, and ordered cross-axis hypotheses
  are all queries against one typed cross-partition form.
  [What is novel in crossform?](https://bbuchsbaum.github.io/crossform/articles/novelty.html)
  is the ledger separating what is demonstrated from what is still gated.

## Install

The repository is not yet public and the package is on neither CRAN nor
R-universe, so today the route is a local checkout:

```r
install.packages(".", repos = NULL, type = "source")
```

Once published, `remotes::install_github("bbuchsbaum/crossform")` and the
R-universe repository `https://bbuchsbaum.r-universe.dev` become the supported
routes.

## Where next

- [Introduction](https://bbuchsbaum.github.io/crossform/articles/introduction.html) — the entry point, from a ready relation to standard errors.
- [Reading the results](https://bbuchsbaum.github.io/crossform/articles/interpreting-results.html) — what each column means, and three interpretive traps.
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
