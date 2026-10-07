# Read one component of a complete geometry

`geometry_component()` returns the packed geometry rows themselves, for
the cases where a view is not enough. Reach for it only after
[`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md);
a query-only result from
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md)
has no stored component to read.

## Usage

``` r
geometry_component(x, component = "total", rows = NULL)
```

## Arguments

- x:

  A complete `effect_form` (including an `effect_geometry`).

- component:

  One of `total`, `coherent`, or `configuration`. `configuration` is
  computed exactly as `total - coherent`.

- rows:

  Optional measurement rows to read. Block-backed stores read only the
  requested rows.

## Value

A numeric matrix with one row per measurement and one column per packed
geometry coordinate (`svec` order for symmetric self forms: the lower
triangle by column, off-diagonal entries scaled by `sqrt(2)`).

## See also

[`query_geometry()`](https://bbuchsbaum.github.io/crossform/reference/query_geometry.md)
to apply a linear query instead of reading packed coordinates, and
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) or
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
for the named scientific views.

Other geometry plans and views:
[`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md),
[`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md),
[`coherence_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/coherence_spectrum.md),
[`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md),
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md),
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md),
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md),
[`example_fmri_effects()`](https://bbuchsbaum.github.io/crossform/reference/example_fmri_effects.md),
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
domain <- abstract_domain(4, id = "component-example")
relation <- relation(
  list(run1 = rbind(a = c(1, 0, 2, 1), b = c(0, 1, 1, 0)),
       run2 = rbind(a = c(1.1, 0.1, 1.9, 0.8), b = c(0.1, 0.9, 1.2, 0.2))),
  domain = domain
)
geometry <- materialize_geometry(plan_geometry(
  relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
  cross_partitions(relation, independence = "independent")
))

# One row per region; three packed coordinates for a 2-effect self form,
# namely G[1,1], G[2,1], and G[2,2].
geometry_component(geometry, "total")
#>      [,1]       [,2] [,3]
#> [1,] 0.55 0.07071068 0.45
#> [2,] 2.30 1.59099026 0.60

# The coherent/configuration split is an exact partition, not a fit.
all.equal(
  geometry_component(geometry, "configuration"),
  geometry_component(geometry, "total") -
    geometry_component(geometry, "coherent")
)
#> [1] TRUE

# Read a single measurement without touching the rest.
geometry_component(geometry, "total", rows = 2)
#>      [,1]    [,2] [,3]
#> [1,]  2.3 1.59099  0.6
```
