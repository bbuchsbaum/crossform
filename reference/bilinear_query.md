# Describe a bilinear geometry query

This compatibility constructor describes a symmetric query on one effect
space. Use
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)
for distinct or unequal axes.

## Usage

``` r
bilinear_query(operator, fixed = TRUE, effects = NULL)
```

## Arguments

- operator:

  A finite square symmetric numeric matrix.

- fixed:

  Whether the query is fixed before local data are inspected.

- effects:

  Optional
  [`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md)
  binding for the operator axes.

## Value

An `effect_query` with `kind = "bilinear"`, carrying the symmetric
`$operator`, the `$fixed` flag, and an optional `$effect_space` binding.

## See also

[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md),
which reads a plan with this query, and
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)
for distinct or unequal axes.

Other geometry plans and views:
[`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md),
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
# A contrast read as a rank-one bilinear query: `t(c) G c` for c = a - b.
contrast <- c(1, -1)
query <- bilinear_query(tcrossprod(contrast))
query$operator
#>      [,1] [,2]
#> [1,]    1   -1
#> [2,]   -1    1

# Evaluate it against a plan without materializing complete geometry.
domain <- abstract_domain(3, id = "bilinear-example")
relation <- relation(
  list(run1 = rbind(a = c(1, 0, 2), b = c(0, 1, 1)),
       run2 = rbind(a = c(1.1, 0.1, 1.9), b = c(0.1, 0.9, 1.2))),
  domain = domain
)
plan <- plan_geometry(
  relation, compile_frame(whole_brain(), domain),
  cross_partitions(relation, independence = "independent")
)
evaluate_geometry(plan, query = query)$values
#>          view1
#> [1,] 0.8333333

# A bilinear query is symmetric by construction; asymmetry is refused.
asymmetric <- try(bilinear_query(matrix(c(1, 2, 3, 4), 2)), silent = TRUE)
conditionMessage(attr(asymmetric, "condition"))
#> [1] "`operator` must be symmetric."
```
