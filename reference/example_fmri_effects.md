# Generate a small task-fMRI relation with known spatial truth

This is the executable newcomer fixture. Four conditions are observed in
independent runs over a small volume, and raw trial responses are
retained through
[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md)
so both point geometry and the admitted fixed-metric uncertainty path
are available.

## Usage

``` r
example_fmri_effects(
  seed = 20260814L,
  dimensions = c(8L, 7L, 5L),
  partitions = 4L,
  trials_per_condition = 8L,
  noise_sd = 0.6,
  spacing = c(3, 3, 3),
  searchlight_radius = 4,
  plant = c("pattern", "mean")
)
```

## Arguments

- seed:

  Nonnegative integer random seed. The caller's random-number state is
  restored on exit.

- dimensions:

  Three volume dimensions, each at least three. The default is the
  smallest volume in which the two planted blocks, and the searchlights
  that touch them, stay disjoint; smaller volumes still work but the
  blocks crowd each other.

- partitions:

  Number of independent runs, at least two.

- trials_per_condition:

  Trials per condition and run, at least two.

- noise_sd:

  Positive residual standard deviation.

- spacing:

  Three positive voxel spacings in millimeters.

- searchlight_radius:

  Positive searchlight radius in millimeters.

- plant:

  Which blocks to plant: `"pattern"`, `"mean"`, or both (the default).
  Dropping one leaves its feature and measurement sets empty; the noise
  draw is unchanged either way, so the two settings differ only in the
  planted signal.

## Value

An `effect_example_effects` list.

## Details

Two animate-versus-inanimate blocks are planted, at opposite ends of the
longest axis and never overlapping. Both carry the same per-voxel
amplitude; only the sign structure differs, so the two halves of the
energy decomposition are each demonstrated by one block:

- the **pattern** block alternates sign between neighboring voxels
  within each axial slice, so its frame-weighted average nearly cancels
  and its reproducible energy is almost entirely `configuration`;

- the **mean** block shifts every voxel the same way, so its energy is
  almost entirely `coherent`.

The returned object is generated, not empirical data. Use the Haxby
exemplar under `exemplars/haxby2001` for the public-data parity
workflow.

## Structure

The returned list has six public elements.

- `$fit`: an
  [`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md)
  carrying a residual channel, so
  [`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
  is admitted. `$fit$relation` is what
  [`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
  takes.

- `$domain`: the full
  [`volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/volume_domain.md)
  the fit was made over.

- `$frame`: a searchlight
  [`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md)
  over that domain, one measurement per voxel.

- `$contrast`: named animate-versus-inanimate weights over the four
  conditions.

- `$model_rdm`: a condition-by-condition category model for
  [`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md).

- `$truth`: what was planted, listed below.

`$truth` holds feature indices into the domain, measurement indices into
the frame, and the generating settings.

- `$planted_features`, `$mean_features`: the two disjoint voxel sets, as
  positions in the domain. `$planted_feature_ids` gives the same voxels
  as domain feature identifiers.

- `$pattern_measurements`, `$mean_measurements`: the searchlights that
  overlap each block. `$signal_measurements` is their sorted union,
  which is what a map-reading example should highlight.

- `$contrast_pattern`: the planted per-voxel contrast profile over the
  whole domain, zero outside the two blocks.

- `$condition_patterns`: the noiseless condition-by-voxel means.

- `$noise_sd`, `$seed`: the generating settings.

## See also

[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
for the next step, then
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md), or
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md); and
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md),
which the retained residual channel makes available.

Other geometry plans and views:
[`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md),
[`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md),
[`coherence_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/coherence_spectrum.md),
[`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md),
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md),
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md),
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md),
[`geometry_component()`](https://bbuchsbaum.github.io/crossform/reference/geometry_component.md),
[`geometry_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/geometry_spectrum.md),
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md),
[`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md),
[`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md),
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md),
[`plot_views`](https://bbuchsbaum.github.io/crossform/reference/plot_views.md),
[`query_geometry()`](https://bbuchsbaum.github.io/crossform/reference/query_geometry.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md),
[`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md)

## Examples

``` r
example <- example_fmri_effects()

# The planted contrast and the two blocks it was planted in are carried
# alongside the data, so any result can be checked against ground truth.
example$contrast
#>  face  body house  tool 
#>   0.5   0.5  -0.5  -0.5 
c(pattern = length(example$truth$planted_features),
  mean = length(example$truth$mean_features))
#> pattern    mean 
#>      19      19 

# Everything a second-moment question needs is already built: a fit with
# residuals, a compiled searchlight frame, and the cross-run pairing.
plan <- plan_geometry(
  example$fit$relation, example$frame,
  cross_partitions(example$fit$relation, independence = "independent")
)
distances <- rdm(plan)
dim(distances$values)
#> [1] 280   6

# Each block reproduces the half of the decomposition it was built from.
energy <- contrast_energy(plan, example$contrast)
round(c(
  pattern_configuration = max(energy$configuration[
    example$truth$pattern_measurements]),
  mean_coherent = max(energy$coherent[example$truth$mean_measurements]),
  elsewhere = max(energy$total[-example$truth$signal_measurements])
), 3)
#> pattern_configuration         mean_coherent             elsewhere 
#>                 4.104                 4.011                 0.030 
```
