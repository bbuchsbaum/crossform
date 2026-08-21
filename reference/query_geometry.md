# Apply a linear query to a complete geometry

`query_geometry()` projects an already materialized geometry through a
fixed linear query. It answers the same question as
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md)
and carries the same view identity; use this form when the complete
geometry already exists and several queries will be read from it.

## Usage

``` r
query_geometry(x, query, component = "total", row_block = 1024L)
```

## Arguments

- x:

  A complete `effect_form` (including an `effect_geometry`).

- query:

  An axis-bound
  [`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md),
  a compatible
  [`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md),
  or a finite physical-coordinate-by-view matrix.

- component:

  Geometry component to query.

- row_block:

  Positive number of measurement rows read at once. This bounds
  packed-geometry memory for block-backed stores.

## Value

An `effect_view`: `$values` has one row per measurement and one column
per query column, alongside `$query`, `$component`, `$index`, and a
`$receipt` recording that this view was projected from the parent
estimand. It is a view, not another geometry.

## See also

[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md)
for the query-first route that never materializes geometry, and
[`geometry_component()`](https://bbuchsbaum.github.io/crossform/reference/geometry_component.md)
for the packed rows.

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
[`geometry_component()`](https://bbuchsbaum.github.io/crossform/reference/geometry_component.md),
[`geometry_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/geometry_spectrum.md),
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md),
[`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md),
[`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md),
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md),
[`plot_views`](https://bbuchsbaum.github.io/crossform/reference/plot_views.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md),
[`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md)

## Examples

``` r
domain <- abstract_domain(4, id = "query-example")
relation <- relation(
  list(run1 = rbind(a = c(1, 0, 2, 1), b = c(0, 1, 1, 0)),
       run2 = rbind(a = c(1.1, 0.1, 1.9, 0.8), b = c(0.1, 0.9, 1.2, 0.2))),
  domain = domain
)
geometry <- materialize_geometry(plan_geometry(
  relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
  cross_partitions(relation, independence = "independent")
))

# The squared cross-generalized distance between the two effects, one row
# per region.
contrast <- bilinear_query(tcrossprod(c(1, -1)))
distance <- query_geometry(geometry, contrast)
distance
#> <effect_view>
#>   measurements: 2
#>  measurement view1
#>           v1  0.90
#>           it  0.65
as.data.frame(distance)
#>   measurement view1
#> 1          v1  0.90
#> 2          it  0.65

# The same query reads the two orthogonal modes, which sum back exactly.
coherent <- query_geometry(geometry, contrast, component = "coherent")
configuration <- query_geometry(
  geometry, contrast, component = "configuration"
)
all.equal(distance$values, coherent$values + configuration$values)
#> [1] TRUE
```
