# Evaluate a fixed query without materializing complete geometry

`evaluate_geometry()` is the query-first route: the requested view is
contracted during execution, so complete packed geometry is never
allocated. Use it whenever a fixed query answers the question; use
[`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md)
only when several readouts must share one stored geometry.

## Usage

``` r
evaluate_geometry(
  x,
  at = NULL,
  over = NULL,
  query = NULL,
  component = c("total", "coherent", "configuration"),
  compute = NULL,
  reporter = NULL
)
```

## Arguments

- x, at, over, compute, reporter:

  As in
  [`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md).

- query:

  A fixed
  [`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md),
  an axis-bound
  [`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)
  for rectangular cross-axis plans, or a packed-coordinate query matrix.

- component:

  One of `total`, `coherent`, or `configuration`.

## Value

A query-only `effect_view`: `$values` has one row per measurement and
one column per query column, with `$component`, `$index`, and the
execution `$receipt`. It carries no packed geometry, so
[`geometry_component()`](https://bbuchsbaum.github.io/crossform/reference/geometry_component.md)
does not apply to it.

## See also

[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
to build the plan,
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md), and
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) for
named views of the same plan.

Other geometry plans and views:
[`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md),
[`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md),
[`coherence_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/coherence_spectrum.md),
[`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md),
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md),
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md),
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
domain <- abstract_domain(4, id = "evaluate-example")
relation <- relation(
  list(run1 = rbind(a = c(1, 0, 2, 1), b = c(0, 1, 1, 0)),
       run2 = rbind(a = c(1.1, 0.1, 1.9, 0.8), b = c(0.1, 0.9, 1.2, 0.2))),
  domain = domain
)
plan <- plan_geometry(
  relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
  cross_partitions(relation, independence = "independent")
)

# A fixed contrast read straight off the plan, with no packed geometry.
view <- evaluate_geometry(plan, query = bilinear_query(tcrossprod(c(1, -1))))
view
#> <effect_view>
#>   measurements: 2
#>  measurement view1
#>           v1  0.90
#>           it  0.65
as.data.frame(view)
#>   measurement view1
#> 1          v1  0.90
#> 2          it  0.65

# A query is mandatory: query-first execution will not guess one.
missing_query <- try(evaluate_geometry(plan), silent = TRUE)
conditionMessage(attr(missing_query, "condition"))
#> [1] "Query-first execution requires an explicit fixed `query`."
```
