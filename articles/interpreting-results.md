# Reading contrast energies, RDMs, and uncertainty

Use this guide after the
[introduction](https://bbuchsbaum.github.io/crossform/articles/introduction.md)
to interpret the numbers returned by a geometry plan. Start with the
contrast columns, then follow the section relevant to your result:
coherent shares, an exploratory map, RDMs, fixed RSA coefficients,
predictive gain or uncertainty. Everything below runs on the generated
fixture, which plants two blocks of equal amplitude at opposite ends of
the volume — one whose sign alternates between neighboring voxels, one
that shifts every voxel the same way — so each half of the energy
decomposition has a known home.

``` r

example <- example_fmri_effects()
pairing <- cross_partitions(
  example$fit$relation,
  independence = "independent",
  generalizes_over = "run"
)
plan <- plan_geometry(example$fit$relation, at = example$frame, over = pairing)
effect <- contrast_energy(plan, example$contrast)
peak <- which.max(effect$total)

pattern <- example$truth$pattern_measurements
mean_block <- example$truth$mean_measurements
c(
  measurements = length(effect$total),
  pattern_searchlights = length(pattern),
  mean_searchlights = length(mean_block),
  overlap = length(intersect(pattern, mean_block))
)
#>         measurements pattern_searchlights    mean_searchlights 
#>                  280                   52                   52 
#>              overlap 
#>                    0
```

The two blocks have equal planted amplitude and size, and their
searchlights do not overlap. This controlled comparison lets us examine
spatial sign structure while keeping those factors fixed.

## The four columns of a contrast view

``` r

head(as.data.frame(effect), 4)
#>   measurement      signed     coherent configuration         total
#> 1           1 -0.10647918  0.007018165 -6.660201e-03  0.0003579634
#> 2           2 -0.01466327 -0.002234820  1.515714e-02  0.0129223176
#> 3           3  0.05821104  0.001297404 -3.916758e-03 -0.0026193537
#> 4           4  0.01998561 -0.000585909 -6.014912e-05 -0.0006460581
#>   coherence_fraction
#> 1                 NA
#> 2                 NA
#> 3                 NA
#> 4                 NA
```

- **`signed`** is the ordinary signed contrast of the frame-weighted
  mean pattern. It is a *first moment*: it has the units and the sign of
  your effect, and it is the only column that answers “which way?”. It
  is not crossvalidated in the reproducibility sense — it is a retained
  marginal.
- **`coherent`** is the crossvalidated energy carried by that
  frame-weighted mean pattern. It is a *second moment* — a squared
  quantity — so it is not the signed effect and not its square either.
- **`configuration`** is the crossvalidated energy in reproducible
  spatial *departures* from the weighted mean, orthogonal to `coherent`
  under the frame’s weights.
- **`total`** is their exact sum.

``` r

max(abs(effect$total - (effect$coherent + effect$configuration)))
#> [1] 2.220446e-16
```

The fixture separates the two energies by construction. Here is the
strongest searchlight in each block:

``` r

strongest <- function(measurements) {
  measurements[which.max(effect$total[measurements])]
}
at <- c(pattern = strongest(pattern), mean = strongest(mean_block))
block_rows <- as.data.frame(effect)[at, -1]
rownames(block_rows) <- names(at)
round(block_rows, 4)
#>          signed coherent configuration  total coherence_fraction
#> pattern -0.0525  -0.0007        4.1037 4.1030                 NA
#> mean     2.0029   4.0109        0.0039 4.0148              0.999
```

Equal amplitude, equal block size, near-equal `total` — and opposite
decompositions. The pattern block’s energy is almost all
`configuration`: a shape that reproduces across runs while the sphere’s
average barely moves, which is why its `signed` value is close to zero
and its signed regional-mean contrast is small. The mean block’s energy
is almost all `coherent`, and its `signed` value is large, because the
whole neighborhood shifts together.

Note the `NA` in the first row. That searchlight is the largest `total`
in the map and it reports no coherence fraction at all, because its
`coherent` estimate landed a hair below zero. That is the subject of the
section after next.

### Negatives are estimates, not errors

Each energy combines products of estimates from different runs. With
independent, unbiased partition estimates, its expectation is zero under
the null. Negative estimates are therefore legitimate. Zero expectation
alone does not imply that positive and negative values are equally
likely.

``` r

c(
  measurements = length(effect$total),
  total_below_zero = sum(effect$total < 0),
  coherent_below_zero = sum(effect$coherent < 0),
  configuration_below_zero = sum(effect$configuration < 0),
  min_total = min(effect$total)
)
#>             measurements         total_below_zero      coherent_below_zero 
#>             280.00000000              79.00000000              95.00000000 
#> configuration_below_zero                min_total 
#>              82.00000000              -0.01089941
```

Keep the negative values in summaries. Clipping them at zero introduces
positive bias, including where the true effect is absent. An average
over prespecified measurements or participants should use the signed
estimates; selecting only positive values changes the quantity being
summarized.

## The coherence fraction and its validity flag

`coherence_fraction` is `coherent / total`: the share of reproducible
energy carried by the local mean rather than by pattern shape. It is
reported only where the raw components actually form a nonnegative
partition — that is, where `total > 0`, `coherent >= 0`, and
`configuration >= 0`. Elsewhere it is `NA`, and a companion logical says
why you got `NA`.

``` r

c(
  measurements = length(effect$total),
  valid = sum(effect$coherence_fraction_valid),
  na = sum(is.na(effect$coherence_fraction))
)
#> measurements        valid           na 
#>          280          143          137
```

``` r

c(
  total_not_positive = sum(effect$total <= 0),
  coherent_negative = sum(effect$coherent < 0),
  configuration_negative = sum(effect$configuration < 0)
)
#>     total_not_positive      coherent_negative configuration_negative 
#>                     79                     95                     82
```

Outside this gate the ratio can leave \[0, 1\], and a denominator near
zero makes it unstable. The package therefore withholds it. `NA` here
means “this measurement has no interpretable coherent share”, not
“missing data”.

The gate is not a test for signal. The strongest pattern-block
measurement above has an `NA` share because its coherent estimate is
slightly negative, even though its total energy is large.

The three sections below are the caveats that matter once you start
reporting those fractions.

## Trap (a): the coherent share shrinks as the sphere grows

Changing the radius changes the spatial quantity being measured, even
when the neural data stay fixed. In this example, a larger sphere
dilutes a fixed block of mean-shift signal with surrounding voxels.

A simple noise-free case explains the effect. Suppose a fraction $`f`$
of a uniformly weighted, locally normalized support carries contrast
$`a`$, and the rest carries zero. Then

``` math
\mathrm{total}=a^2 f,\qquad
\mathrm{coherent}=a^2 f^2,\qquad
\mathrm{configuration}=a^2 f(1-f).
```

The coherent share is $`f`$. Enlarging the support lowers that share as
the signal occupies less of it. Configuration need not stay constant; it
is the remainder after the squared local mean is removed.

Extra noise dimensions can change variability and which noisy ratios
pass the validity gate. They do not, by themselves, add positive
expected crossvalidated energy. Read the radius sweep as a change in the
declared measurement, not as a count of neural dimensions.

``` r

radii <- c(2, 3, 4.5, 6)
sweep <- t(vapply(radii, function(r) {
  frame_r <- compile_frame(searchlights(radius = r), example$domain)
  effect_r <- contrast_energy(
    plan_geometry(example$fit$relation, at = frame_r, over = pairing),
    example$contrast
  )
  c(
    median_support = median(Matrix::rowSums(frame_r$weights != 0)),
    valid = sum(effect_r$coherence_fraction_valid),
    median_fraction = median(effect_r$coherence_fraction, na.rm = TRUE),
    pattern_block = median(effect_r$coherence_fraction[pattern], na.rm = TRUE),
    mean_block = median(effect_r$coherence_fraction[mean_block], na.rm = TRUE)
  )
}, numeric(5)))
data.frame(radius = radii, round(sweep, 3))
#>   radius median_support valid median_fraction pattern_block mean_block
#> 1    2.0              1   138           1.000         1.000      1.000
#> 2    3.0              6   143           0.181         0.090      0.347
#> 3    4.5             14   175           0.135         0.047      0.305
#> 4    6.0             22   224           0.084         0.047      0.289
```

The mechanism is visible directly. Hold the radius at 4.5 mm and
compare, over the searchlights whose support touches the mean-shift
block, the two energies against the squared weighted share of planted
voxels in each support:

``` r

weighted_share <- function(frame, features) {
  as.numeric(Matrix::rowSums(frame$weights[, features, drop = FALSE]) /
               Matrix::rowSums(frame$weights))
}

frame_45 <- compile_frame(searchlights(radius = 4.5), example$domain)
effect_45 <- contrast_energy(
  plan_geometry(example$fit$relation, at = frame_45, over = pairing),
  example$contrast
)
share <- weighted_share(frame_45, example$truth$mean_features)
touching <- which(share > 0)

c(
  searchlights = length(touching),
  coherent = round(cor(effect_45$coherent[touching], share[touching]^2), 3),
  configuration = round(
    cor(effect_45$configuration[touching], share[touching]^2), 3
  )
)
#>  searchlights      coherent configuration 
#>        80.000         0.998         0.408
```

The coherent energy is almost a deterministic function of the squared
planted share. The configuration energy is not. This is the dilution
mechanism predicted by the simple mean-shift example.

At radius 2 every support is a single voxel, so the weighted mean *is*
the pattern and the fraction is exactly 1 by construction — in both
blocks. It then falls. The decomposition keeps the two blocks in the
right order at every larger radius, the mean block several times above
the pattern block, so the ordering is real. But the mean block’s own
fraction still slides from 0.347 to 0.289 while nothing about the
planted signal changes. A sentence such as “roughly 30% of the
reproducible signal was carried by the mean response” is therefore a
statement about your radius as much as about the brain. Report the
radius and the support sizes alongside the fraction, and never compare
fractions across analyses with different frames.

## Trap (b): “coherent” means coherent under *your* weights

The decomposition is orthogonal in the frame-weighted inner product, so
the weights are part of the definition. Reweighting the *same* voxels
moves energy between the two components — and changes the total.

``` r

frame_uniform <- compile_frame(searchlights(radius = 4.5), example$domain)

# Same supports, Gaussian instead of uniform weights.
coords <- example$domain$coordinates
uniform_weights <- as.matrix(frame_uniform$weights)
centers <- match(frame_uniform$index$measurement, example$domain$feature_ids)
gaussian_weights <- uniform_weights
for (i in seq_len(nrow(uniform_weights))) {
  members <- which(uniform_weights[i, ] != 0)
  offsets <- coords[members, , drop = FALSE] -
    matrix(coords[centers[i], ], length(members), 3, byrow = TRUE)
  w <- exp(-rowSums(offsets^2) / (2 * 3^2))
  gaussian_weights[i, ] <- 0
  gaussian_weights[i, members] <- w / sum(w)
}
frame_gaussian <- additive_frame(
  gaussian_weights, normalization = "local", domain = example$domain
)

identical(which(uniform_weights != 0), which(gaussian_weights != 0))
#> [1] TRUE
```

Read the two frames at the searchlight sitting on each block’s center.
That measurement is chosen from geometry alone, not from the results, so
nothing below is a selected extreme.

``` r

effect_uniform <- contrast_energy(
  plan_geometry(example$fit$relation, at = frame_uniform, over = pairing),
  example$contrast
)
effect_gaussian <- contrast_energy(
  plan_geometry(example$fit$relation, at = frame_gaussian, over = pairing),
  example$contrast
)

centered_on <- function(features) {
  centroid <- colMeans(coords[features, , drop = FALSE])
  which.min(rowSums(
    (coords[centers, , drop = FALSE] -
      matrix(centroid, length(centers), 3, byrow = TRUE))^2
  ))
}
middle <- c(
  pattern = centered_on(example$truth$planted_features),
  mean = centered_on(example$truth$mean_features)
)

reweighted <- cbind(
  uniform_coherent = effect_uniform$coherent[middle],
  gaussian_coherent = effect_gaussian$coherent[middle],
  uniform_total = effect_uniform$total[middle],
  gaussian_total = effect_gaussian$total[middle]
)
rownames(reweighted) <- names(middle)
round(reweighted, 4)
#>         uniform_coherent gaussian_coherent uniform_total gaussian_total
#> pattern           0.2678            0.1320        4.0762         4.0969
#> mean              4.0671            4.0252        4.0604         4.0199
```

Identical voxels, identical data, different numbers — but not equally
different. The mean block’s coherent energy barely moves: a shift that
is the same in every voxel is the same shift under any nonnegative
weights. The pattern block’s coherent energy roughly halves because its
alternating signal cancels differently under the two weightings. Both
are valid weighted projections of the same spatial pattern. Both totals
move.

This is not instability: each frame defines a different, well-specified
estimand, and the plan records which one you asked for. It does mean
that “the coherent component” is not a property of the brain region. It
is a property of the region *and* the weights you chose to average it
with.

## Trap (c): reported fractions are a selected sample

Because the fraction exists only where the components form a nonnegative
partition, the set of measurements that report one is not a random
subset of the map. The gate depends on the same quantities that appear
in the ratio, so it is enriched for measurements with real signal — and,
among the rest, for whichever noise realizations happened to make both
components positive.

``` r

frame_small <- compile_frame(searchlights(radius = 3), example$domain)
effect_small <- contrast_energy(
  plan_geometry(example$fit$relation, at = frame_small, over = pairing),
  example$contrast
)
signal <- example$truth$signal_measurements
valid <- effect_small$coherence_fraction_valid

c(
  signal_measurements = length(signal),
  signal_reporting_a_fraction = sum(valid[signal]),
  other_measurements = length(valid) - length(signal),
  other_reporting_a_fraction = sum(valid[-signal])
)
#>         signal_measurements signal_reporting_a_fraction 
#>                         104                         100 
#>          other_measurements  other_reporting_a_fraction 
#>                         176                          43
```

``` r

round(c(
  pattern_block = mean(valid[pattern]),
  mean_block = mean(valid[mean_block]),
  elsewhere = mean(valid[-signal])
), 3)
#> pattern_block    mean_block     elsewhere 
#>         0.942         0.981         0.244
```

``` r

c(
  share_of_map_reporting = mean(valid),
  share_of_reporters_that_are_signal = mean(seq_along(valid)[valid] %in% signal)
)
#>             share_of_map_reporting share_of_reporters_that_are_signal 
#>                          0.5107143                          0.6993007
```

In this fixture, both blocks clear the gate at similar rates. That need
not hold for every signal shape or noise regime. Under a quarter of the
remaining measurements report a fraction, on noise alone. The result is
that the average fraction “over the map” is really an average over a set
that is 70% planted signal. That is fine if you say so, and misleading
if you do not. If you want a map-wide summary, define the denominator
explicitly — for instance the share of *all* measurements that report a
coherent share at all, the first number above — rather than averaging
the reported values and calling the result a property of the map.

The same warning applies to any post-hoc gate: thresholding on
`total > 0` and then summarizing `coherent` conditions the summary on
the numerator’s sibling.

## Reading a map without ground truth

Every figure so far could mark the right answer, because the fixture
carries `truth$signal_measurements`. Your data will not. This section is
about what the same picture tells you when nothing is marked.

``` r

plot(effect)
```

![The same contrast view with no highlight supplied. Nothing here
identifies the two planted blocks; the reading below uses only what is
on the page.](interpreting-results_files/figure-html/blind-map-1.png)

The same contrast view with no highlight supplied. Nothing here
identifies the two planted blocks; the reading below uses only what is
on the page.

### Zero is a reference, not a threshold

Independent, unbiased partition estimates give null energy an
expectation of zero. That provides a useful reference line. It does not
tell you which locations are null, how variable their estimates are, or
what threshold controls false positives across the map.

Here the generated truth lets us examine measurements without planted
signal:

``` r

quiet <- setdiff(seq_along(effect$total), example$truth$signal_measurements)
round(c(
  measurements = length(quiet),
  mean = mean(effect$total[quiet]),
  sd = stats::sd(effect$total[quiet]),
  fraction_below_zero = mean(effect$total[quiet] < 0),
  largest = max(effect$total[quiet])
), 4)
#>        measurements                mean                  sd fraction_below_zero 
#>            176.0000              0.0027              0.0085              0.4489 
#>             largest 
#>              0.0297
```

Without the truth, the median and MAD summarize the observed map. Their
meaning depends on its mixture of signal and noise; they are not
automatically estimates of a null distribution:

``` r

round(c(
  median_over_all = stats::median(effect$total),
  mad_over_all = stats::mad(effect$total),
  fraction_below_zero_over_all = mean(effect$total < 0)
), 4)
#>              median_over_all                 mad_over_all 
#>                       0.0088                       0.0198 
#> fraction_below_zero_over_all 
#>                       0.2821
```

In this fixture the blind median lies near the known null center, while
the MAD is more than twice the standard deviation of the known null
measurements. The two signal blocks occupy 37% of the map and affect
that summary. A different signal prevalence or spatial noise pattern
could change both comparisons.

### Candidates, not findings

Sorted, this map has a visible shoulder: rank 104 is `0.478` and rank
105 is `0.030`, while rank 1 is `4.10`. The shoulder falls exactly at
the edge of the planted set — all 104 planted searchlights outrank every
other measurement. That is a fact about this fixture, not a method.

``` r

round(sort(effect$total, decreasing = TRUE)[c(1, 10, 40, 104, 105, 150, 280)], 4)
#> [1]  4.1030  3.4728  1.7863  0.4779  0.0297  0.0067 -0.0109
```

A visible shoulder is useful for exploration, but this view supplies no
spatial threshold or selected-peak p-value. A population analysis has
its own declared target and assumptions; it does not retroactively
calibrate this map. Treat selected locations as candidates for further
study:

- **Describe its estimated composition.** High configuration attributes
  more estimated energy to departures from the local mean; high coherent
  energy attributes more to that mean. Selection can exaggerate either
  component, so the description does not establish replication.
- **Attach a within-measurement standard error at a measurement you
  chose in advance**, using
  [`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
  with `target = "null"` — and only if the relation kept its residual
  channel, which means fitting with
  [`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md)
  or
  [`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md)
  rather than importing betas.
- **Take the candidates somewhere else**: a held-out session, an
  independent dataset, or a group-level analysis performed outside this
  package.

### A pointwise error bar does not account for peak selection

If you find `peak <- which.max(effect$total)` and then report a standard
error at `peak`, the pointwise standard error does not account for
having chosen an extreme estimate. An interval formed from it need not
retain its nominal coverage after selection. The peak summaries here
illustrate the fixture; they do not perform inference for a selected
location.

Specify the measurement in advance, or treat the peak as a candidate and
estimate it somewhere you did not select it.

## Reading the RDM

`peak` is the pattern block’s strongest searchlight, so the two
within-category pairs should sit near zero and the four across-category
pairs near the planted separation.

``` r

distances <- rdm(plan)
distances$pairs
#>    left right
#> 1  face  body
#> 2  face house
#> 3  face  tool
#> 4  body house
#> 5  body  tool
#> 6 house  tool
round(distances$values[peak, ], 3)
#>  face - body face - house  face - tool body - house  body - tool house - tool 
#>       -0.004        4.078        4.334        3.881        4.124        0.009
```

For conditions *i* and *j*,
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md)
reports `d_ij = G_ii + G_jj - 2 G_ij`. Under cross-partition pairing
this is a **signed crossvalidated squared Euclidean** distance; with a
fixed neural precision metric it is a squared Mahalanobis distance.
Three consequences:

1.  **They are squared.** Do not take square roots to “get distances
    back” — the square root of a negative estimate is undefined, and the
    square root of a noisy positive one is biased.
2.  **They can be negative**, for the same reason the energies can, and
    for the same reason they are kept.
3.  **They are not `1 - Pearson correlation`.** Correlation distance
    divides each pattern by its own norm, which is not a fixed linear
    operation, so it cannot be a query against this geometry. The
    package refuses to substitute it silently; the reasoning and the
    conditions under which the boundary could be crossed are in the
    correlation-distance policy
    ([online](https://bbuchsbaum.github.io/crossform/articles/correlation-distance-policy.md),
    or
    [`vignette("correlation-distance-policy")`](https://bbuchsbaum.github.io/crossform/articles/correlation-distance-policy.md)
    offline).

Because the RDM is a geometry query and not a decomposition, it does not
tell the two blocks apart. The strongest searchlight in each block
returns much the same six distances, even though their energies split in
opposite directions:

``` r

round(rbind(
  pattern = distances$values[at[["pattern"]], ],
  mean = distances$values[at[["mean"]], ]
), 3)
#>         face - body face - house face - tool body - house body - tool
#> pattern      -0.004        4.078       4.334        3.881       4.124
#> mean          0.005        4.147       3.960        4.055       3.887
#>         house - tool
#> pattern        0.009
#> mean          -0.015
```

The total RDM alone cannot distinguish those spatial organizations. Use
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
for a contrast decomposition;
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) also
accepts a `component` when you want distances from a particular spatial
component.

## Reading RSA coefficients

[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md)
regresses the measured distances on your model RDMs. It is ordinary
least squares, not a rank correlation.

``` r

category <- rsa(plan, models = list(category = example$model_rdm))
colnames(category$coefficients)
#> [1] "(Intercept)" "category"
round(category$coefficients[peak, ], 4)
#> (Intercept)    category 
#>      0.0025      4.1017
```

There are two columns because an **intercept is fitted by default**. It
absorbs the overall level of signed distances, letting the model slope
describe their relationship after accounting for that level. Without an
intercept, the slope must explain both level and shape: uniformly large
distances can increase the coefficient of a positive-valued model even
when its shape matches poorly.

``` r

no_intercept <- rsa(
  plan, models = list(category = example$model_rdm), intercept = FALSE
)
c(
  with_intercept = category$coefficients[peak, "category"],
  without_intercept = no_intercept$coefficients[peak, "category"]
)
#>    with_intercept.category without_intercept.category 
#>                   4.101730                   4.104211
```

Drop the intercept only when you have a reason to believe the level is
meaningful and shared. Coefficients are in units of distance per unit of
model RDM, so they are comparable across measurements within one
analysis but not across analyses with differently scaled model RDMs.

## Reading a learned geometry prediction

[`fit_geometry()`](https://bbuchsbaum.github.io/crossform/reference/fit_geometry.md)
learns a model-supported PSD form on training data.
[`score_geometry()`](https://bbuchsbaum.github.io/crossform/reference/score_geometry.md)
evaluates that frozen form on independent runs using
$`2\langle F,G_{\mathrm{test}}\rangle_F-\|F\|_F^2`$. The result is a
signed gain in squared geometry units, relative to predicting zero.
Positive gain supports improved prediction in expectation under the
declared independence assumptions; negative gain means the prediction’s
amplitude costs more than its independent evidence earns.

For each fitted mode, keep three quantities apart:

| Quantity    | Meaning                                                     |
|-------------|-------------------------------------------------------------|
| `amplitude` | Strength learned on training data                           |
| `evidence`  | Signed geometry in that frozen direction on evaluation data |
| `gain`      | Twice amplitude times evidence, minus amplitude squared     |

A positive training amplitude does not establish replication. Even
positive test evidence can produce negative gain if the training
amplitude is too large. Do not drop such a mode after inspecting final
evaluation results: that would select a different prediction using the
test data. Rank, penalty, model weights and preprocessing belong to
training or inner validation.

Gain is neither a fraction explained nor a trace energy. Under zero
signal, its conditional expectation is minus the prediction’s squared
norm; the test inner product alone is centered on zero. Likewise, the
descriptive projection from `latent_geometry(rank = ...)` and the trace
split from
[`model_geometry()`](https://bbuchsbaum.github.io/crossform/reference/model_geometry.md)
do not measure independent predictive improvement.

The [predictive geometry
guide](https://bbuchsbaum.github.io/crossform/articles/predictive-geometry.md)
([`vignette("predictive-geometry", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/predictive-geometry.md)
offline) executes a case with mode gains 4 and -1, total gain 3, and an
exact rank-zero baseline. It also shows the observation-origin manifest
that keeps training and evaluation separate. These scores concern new
runs on the same conditions; they do not establish transfer to new
conditions or supply generic standard errors for learned ranks and
eigenvalues.

## Uncertainty: two targets, and what “refused” means

[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
requires an explicit `target`, because the two choices answer different
questions:

``` r

null_se <- sqrt(sampling_covariance(
  rdm_sampling_covariance(plan, example$fit, target = "null", at = peak)
))
plugin_se <- sqrt(sampling_covariance(
  rdm_sampling_covariance(plan, example$fit, target = "plugin", at = peak)
))
round(rbind(distance = distances$values[peak, ],
            null = null_se, plugin = plugin_se), 4)
#>          face - body face - house face - tool body - house body - tool
#> distance     -0.0043       4.0781      4.3337       3.8806      4.1245
#> null          0.0142       0.0142      0.0142       0.0142      0.0142
#> plugin        0.0209       0.2377      0.2487       0.2321      0.2426
#>          house - tool
#> distance       0.0092
#> null           0.0142
#> plugin         0.0266
```

- **`target = "null"`** evaluates the signal-dependent variance term at
  a fixed zero effect. It is **exact on the variance scale** under the
  declared model, and it is what you want when calibrating a test of no
  effect. Notice it is the same for every pair: with no signal assumed,
  the variance depends only on the noise.
- **`target = "plugin"`** substitutes the estimated partition-mean
  pattern for the unknown signal. It tracks the actual distances — large
  where the distance is large, and back down toward the null value for
  the two within-category pairs — which is what you want when reporting
  an error bar around a nonzero estimate. It is **biased upward**, by a
  term that shrinks like 1/M² in the number of partitions and is largest
  when noise dominates the true distances, so read it as mildly
  conservative. See
  [`?rdm_sampling_covariance`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
  for the exact expression.

Both targets share one estimated ingredient: the residual covariance
$`\Sigma_w`$. “Exact” above means exact *given* that covariance — and
because the noise term is quadratic in it, a raw plug-in would inflate
every standard error by $`\sqrt{1+(1+P_{\text{eff}})/\nu}`$. crossform
applies the Wishart-unbiased correction and reports the two numbers that
govern it: the residual degrees of freedom $`\nu`$, and
$`P_{\text{eff}}`$, the number of residual directions this measurement’s
support actually spends variance on.

``` r

covariance <- rdm_sampling_covariance(
  plan, example$fit, target = "null", at = peak
)
c(residual_df = covariance$source$residual_df,
  effective_dimension = round(
    covariance$source$residual_effective_dimension, 2
  ))
#>         residual_df effective_dimension 
#>              112.00                5.88
```

A large support bought with few runs runs out of residual information,
and when $`\nu`$ falls below $`P_{\text{eff}}`$ the call refuses rather
than reporting a confidently small number.

There is no default. Choosing one silently would mean the same function
reported two different quantities depending on context.

This is a **within-measurement, within-participant** covariance under an
equal-partition, fixed-metric, separable error model. It is not a
confidence interval, not a spatial random-field correction, and not
group inference.

### Refusal is an answer

If the relation cannot support the analytic law, the package says so
with a classed condition rather than returning a plausible-looking
number. Ask before provoking:

``` r

capabilities <- sampling_capabilities(beta_plan, beta_only)
capabilities$available
#> [1] FALSE
```

| reason | remedy |
|:---|:---|
| missing_error_channel | Refit raw observations with [`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md). |

A relation built only from precomputed beta matrices has no residuals
and no residual degrees of freedom, so the analytic law is unavailable.
Provoking it anyway yields a refusal you can inspect as data:

``` r

refusal <- catch_refusal(
  rdm_sampling_covariance(beta_plan, beta_only, target = "null", at = 1)
)
refusal$capability
#> [1] "sampling_covariance"
refusal$reasons
#> [1] "missing_error_channel"
```

Here the refusal identifies missing statistical information. Other
refusals can identify unsupported operations; read `$reasons` and
`$remedies` to tell which requirement failed. The point estimates —
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) —
remain perfectly valid on such a relation; only the error bar is
withheld.

## Summary

| Reading | Do | Do not |
|----|----|----|
| Negative energies or distances | keep and average them | truncate at zero |
| `signed` | read direction and units | square it and call it energy |
| `coherence_fraction` | report with radius and support sizes | compare across frames |
| `NA` fraction | read as “no interpretable share” | impute or drop silently |
| RDM values | treat as signed squared distances | take square roots, or read as `1 - r` |
| RSA coefficient | keep the intercept | compare across differently scaled models |
| Predictive gain | retain its sign and the frozen prediction cost | call it a fraction explained or a null-centered energy |
| Mode evidence | distinguish training amplitude, test evidence and gain | remove harmful modes using final test outcomes |
| `target = "null"` | calibrate a test of no effect | use as an error bar on a nonzero distance |
| `target = "plugin"` | error bar on an estimate, read as conservative | treat as exact |
| a refusal | read `$reasons` and `$remedies` | retry until something returns a number |
| a measurement above the floor | call it a candidate and say what kind it is | threshold it and report it as a finding |
| an error bar at [`which.max()`](https://rdrr.io/r/base/which.min.html) | name the measurement in advance instead | read it as uncertainty about the peak |

## See also

- [`vignette("introduction", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/introduction.md)
  — the workflow these results come from.
- [Correlation-distance
  policy](https://bbuchsbaum.github.io/crossform/articles/correlation-distance-policy.md),
  offline
  [`vignette("correlation-distance-policy")`](https://bbuchsbaum.github.io/crossform/articles/correlation-distance-policy.md)
  — why
  [`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md)
  will not give you `1 - r`.
- [Failure
  gallery](https://bbuchsbaum.github.io/crossform/articles/failure-gallery.md),
  offline
  [`vignette("failure-gallery")`](https://bbuchsbaum.github.io/crossform/articles/failure-gallery.md)
  — six realistic errors and the refusals that stop them.
- [`?example_fmri_effects`](https://bbuchsbaum.github.io/crossform/reference/example_fmri_effects.md)
  for the two planted blocks, and
  [`?contrast_energy`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
  [`?rdm`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
  [`?rsa`](https://bbuchsbaum.github.io/crossform/reference/rsa.md),
  [`?rdm_sampling_covariance`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
  for the exact contracts.
