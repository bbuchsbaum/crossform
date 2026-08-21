# Aggregate edge sufficient statistics before normalization

This names a distinct estimand. Raw sufficient statistics are combined
by partition weights before normalization and transformation. It is the
default for
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md)
and
[`coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling.md),
where averaging raw partition products and normalizing once is the
intended estimator.

## Usage

``` r
aggregate_first()
```

## Value

An `effect_partition_reducer` recording `$kind`, `$weight_convention`,
and the stage `$order` (`"aggregate_first"`).

## See also

[`reduce_partitions()`](https://bbuchsbaum.github.io/crossform/reference/reduce_partitions.md)
for the edge-first order, and
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md),
which records this choice in its plan identity.

Other geometry plans and views:
[`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md),
[`coherence_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/coherence_spectrum.md),
[`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md),
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
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
# Aggregate-first: raw partition products are pooled, then normalized once.
reducer <- aggregate_first()
reducer$order
#> [1] "aggregate_first"
reducer$weight_convention
#> [1] "normalized_unit_mass"

# This is the default measurement-pipeline order, so it is what
# `measurement_form()` and `coupling()` record when `reducer` is omitted.
identical(reducer, eval(formals(measurement_form)$reducer))
#> [1] TRUE
```
