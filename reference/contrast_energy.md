# Read a contrast as signed, coherent, configuration, and total evidence

This view never removes a mean pattern or reruns an analysis. It reads
the exact coherent/configuration decomposition already contained in `x`.

## Usage

``` r
contrast_energy(x, ...)

# Default S3 method
contrast_energy(x, weights, remove_univariate = FALSE, ...)
```

## Arguments

- x:

  An `effect_geometry_plan` or complete `effect_geometry`.

- ...:

  Passed to the method. The generic dispatches on `x`: a geometry plan
  or a complete geometry is read by the default method documented here,
  and an `effect_population_result` from
  [`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
  by the group-level method in
  [population_views](https://bbuchsbaum.github.io/crossform/reference/population_views.md).

- weights:

  One finite contrast weight per experimental effect. Named weights are
  reordered to the relation's declared effect order and must name every
  effect exactly once, so naming them is the safe form. Unnamed weights
  are accepted positionally, in the order given by
  `x$task$left_relation$effects`; the returned `$weights` always carries
  the effect names, so print it to confirm the alignment you intended.

- remove_univariate:

  Must be omitted or `FALSE`. Destructive demeaning is refused: the
  coherent/configuration/total decomposition already reports the common
  spatial mode and its orthogonal remainder as an exact, non-destructive
  partition.

## Value

An `effect_contrast_view` with one value per measurement in `$signed`
(the signed contrast of the local weighted mean), `$coherent`,
`$configuration`, and `$total = coherent + configuration`, plus
`$coherence_fraction` (reported only where the raw cross-generalized
components form a nonnegative partition, flagged by
`$coherence_fraction_valid`), the aligned `$weights`, `$index`, and
`$receipt`.

## Structure

Each value element holds one number per spatial measurement, in `$index`
order.

- `$signed`: the signed contrast of the local weighted mean. It keeps
  its sign, so it says which way the effect goes; the energies below
  cannot.

- `$coherent`: the part of the energy carried by the measurement's own
  weighted common spatial mode.

- `$configuration`: the orthogonal remainder, the pattern part.

- `$total`: `$coherent + $configuration`, the crossvalidated energy of
  the contrast. Cross-generalized values may be negative.

- `$coherence_fraction`: `$coherent / $total`, and `NA` wherever the raw
  components do not form a nonnegative partition.

- `$coherence_fraction_valid`: `TRUE` exactly where that fraction was
  reported.

- `$weights`: the contrast, reordered to the relation's effect order and
  named. Print it to confirm the alignment.

- `$index`: the measurement identifiers, one per value, carried from the
  frame's `$index$measurement`.

- `$receipt`: the execution receipt for the run that produced the
  values.

Any element not listed here is internal and may change.

## See also

[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
to build `x`,
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) and
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) for
the other named views, and
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md)
for the same total under a declared noise-precision metric.

Other geometry plans and views:
[`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md),
[`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md),
[`coherence_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/coherence_spectrum.md),
[`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md),
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md),
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md),
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md),
[`example_fmri_effects()`](https://bbuchsbaum.github.io/crossform/reference/example_fmri_effects.md),
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
# Where does an animate-versus-inanimate pattern reproduce across runs?
example <- example_fmri_effects()
plan <- plan_geometry(
  example$fit$relation, example$frame,
  cross_partitions(
    example$fit$relation,
    independence = "independent", generalizes_over = "run"
  )
)
effect <- contrast_energy(plan, example$contrast)

# The strongest searchlight falls inside the planted signal, and its
# energy splits exactly into coherent and configuration parts.
peak <- which.max(effect$total)
c(
  signed = effect$signed[peak],
  coherent = effect$coherent[peak],
  configuration = effect$configuration[peak],
  total = effect$total[peak],
  planted = peak %in% example$truth$signal_measurements
)
#>        signed      coherent configuration         total       planted 
#> -0.0524585968 -0.0007331921  4.1037034421  4.1029702500  1.0000000000 

# Demeaning the univariate signal away is refused: the decomposition
# already separates the common spatial mode from its remainder.
refusal <- catch_refusal(
  contrast_energy(plan, example$contrast, remove_univariate = TRUE)
)
refusal$capability
#> [1] "nondestructive_decomposition"
refusal$remedies
#> [1] "Read `$coherent`, `$configuration`, and `$total` from the returned contrast view."
```
