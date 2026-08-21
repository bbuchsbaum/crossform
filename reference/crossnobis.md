# Evaluate a signed local crossnobis contrast

The result is the query-first evidence `tr(c c^T B_a K B_b^T)`
aggregated over the plan's independent partition edges. This function is
a validating view over the ordinary geometry compiler; it does not
introduce a second numerical engine. Negative finite estimates are
retained.

## Usage

``` r
crossnobis(x, weights)
```

## Arguments

- x:

  An `effect_geometry_plan` carrying either an explicit fixed
  [`noise_precision()`](https://bbuchsbaum.github.io/crossform/reference/noise_precision.md)
  metric, or a provenance-frozen learned metric schedule compiled from a
  recipe by
  [`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md)
  or `plan_geometry(metric = )`.

- weights:

  One finite contrast weight per experimental effect.

## Value

An `effect_crossnobis_view` whose `$values` holds one signed
crossvalidated squared Mahalanobis value per spatial measurement, with
the aligned `$contrast`, the named `$estimand`, the `$metric` and
`$pairing` identities, `$index`, and the executed `$receipt`.

## Structure

One signed value per spatial measurement, alongside the declarations
that make it a Mahalanobis reading.

- `$values`: the signed crossvalidated squared Mahalanobis contrast, one
  per measurement, in `$index` order. Negative estimates are retained.

- `$contrast`: the weights, reordered to the relation's effect order and
  named.

- `$estimand`: the named estimand,
  `"crossvalidated_squared_mahalanobis_contrast"`.

- `$metric`: the signature of the noise-precision metric the values were
  read under.

- `$pairing`: the identity of the independent partition edges they were
  generalized over.

- `$index`: the measurement identifiers, one per value, carried from the
  frame's `$index$measurement`.

- `$receipt`: the execution receipt for the run that produced the
  values.

The `$metadata` block and any other element not listed here are internal
and may change.

## Refusal

A plan carrying the implicit identity metric, or a fixed metric that was
not built by
[`noise_precision()`](https://bbuchsbaum.github.io/crossform/reference/noise_precision.md),
signals an `effect_capability_refusal` with capability
`"declared_noise_metric"` in namespace `"geometry_views"`. Inspect it
with
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md).

## One estimand, two views

On a fixed-metric geometry plan, `crossnobis(x, weights)` is the named
Mahalanobis reading of exactly `contrast_energy(x, weights)$total`: the
same compiled estimand, exposed as a single signed value. When the
analysis also needs the signed endpoint marginals or the exact
coherent/configuration decomposition of the same quantity, call
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
on the same plan; every component it returns inherits the plan's fixed
metric.

## See also

[`noise_precision()`](https://bbuchsbaum.github.io/crossform/reference/noise_precision.md)
for the fixed metric a geometry plan needs,
[`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md)
for the learned-metric route, and
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
for the decomposed reading of the same estimand.

Other geometry plans and views:
[`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md),
[`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md),
[`coherence_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/coherence_spectrum.md),
[`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md),
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md),
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
domain <- abstract_domain(3, id = "crossnobis-example")
run1 <- rbind(a = c(1, 0, 0), b = c(0, 1, 0))
run2 <- rbind(a = c(1.1, 0, 0), b = c(0, 0.9, 0))
relation <- relation(list(run1 = run1, run2 = run2), domain = domain)
plan <- plan_geometry(
  relation,
  compile_frame(whole_brain(), domain),
  cross_partitions(relation, independence = "independent"),
  metric = noise_precision(diag(3), domain)
)
result <- crossnobis(plan, c(a = 1, b = -1))
result
#> <effect_crossnobis_view>
#>   measurements: 1
#>   contrast:     a 1, b -1
#>  measurement crossnobis
#>  whole_brain     0.6667
as.data.frame(result)
#>   measurement crossnobis
#> 1 whole_brain  0.6666667
```
