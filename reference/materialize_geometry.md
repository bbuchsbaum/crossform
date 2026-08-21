# Materialize a complete cross-generalized geometry

`materialize_geometry()` explicitly materializes complete total and
coherent packed geometry. Prefer
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
followed by
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md)
when a fixed query is sufficient. The relation compatibility form
compiles that same plan first; it does not use a second execution path.

## Usage

``` r
materialize_geometry(
  x,
  at = NULL,
  over = NULL,
  storage = c("memory", "block"),
  storage_path = NULL,
  compute = NULL,
  reporter = NULL
)
```

## Arguments

- x:

  An `effect_geometry_plan` or, for compatibility, an `effect_relation`.

- at:

  A compiled additive `effect_frame`, required only with a relation.

- over:

  An `effect_pairing`, required only with a relation.

- storage:

  Either `"memory"` or `"block"`.

- storage_path:

  Durable directory for block-backed geometry.

- compute:

  A sequential
  [`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md)
  when `x` is a relation. A compiled plan already owns this choice.

- reporter:

  Optional nonsemantic coordinator-side event reporter.

## Value

An `effect_geometry` holding the packed `$total` and `$coherent` stores,
the signed endpoint `$marginals`, the `$effects` names, the measurement
`$index`, and the execution `$receipt`. Read it with
[`geometry_component()`](https://bbuchsbaum.github.io/crossform/reference/geometry_component.md),
[`query_geometry()`](https://bbuchsbaum.github.io/crossform/reference/query_geometry.md),
or the named views.

`$result_capability` is `"complete_form"` here and `"query_only"` on a
view; it is the one field whose vocabulary is the same across every
result kind, so branch on it rather than on the displayed
`$completeness`, which reads `"full"` on an `effect_geometry` and
`"complete"` on an `effect_measurement_form`. Those two words are
frozen: `completeness` is hashed into the measurement contract signature
and into evidence-task identity, so reconciling them would invalidate
recorded identities.

## See also

[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
then
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md)
for the query-first route;
[`geometry_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/geometry_spectrum.md),
which requires a complete geometry.

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
[`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md),
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md),
[`plot_views`](https://bbuchsbaum.github.io/crossform/reference/plot_views.md),
[`query_geometry()`](https://bbuchsbaum.github.io/crossform/reference/query_geometry.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md),
[`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md)

## Examples

``` r
domain <- abstract_domain(4, id = "materialize-example")
relation <- relation(
  list(run1 = rbind(a = c(1, 0, 2, 1), b = c(0, 1, 1, 0)),
       run2 = rbind(a = c(1.1, 0.1, 1.9, 0.8), b = c(0.1, 0.9, 1.2, 0.2))),
  domain = domain
)
plan <- plan_geometry(
  relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
  cross_partitions(relation, independence = "independent")
)

# Materialize only when several different readouts will share the same
# packed geometry; a single fixed query does not need this step.
geometry <- materialize_geometry(plan)
geometry$effects
#> [1] "a" "b"
geometry_component(geometry, "total")
#>      [,1]       [,2] [,3]
#> [1,] 0.55 0.07071068 0.45
#> [2,] 2.30 1.59099026 0.60

# The complete object unlocks readouts that a query-only view cannot give,
# such as the signed eigenvalue spectrum.
geometry_spectrum(geometry)$values
#>          root1      root2
#> [1,] 0.5707107 0.42928932
#> [2,] 2.8600089 0.03999113
```
