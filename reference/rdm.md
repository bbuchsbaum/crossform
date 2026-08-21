# Read squared experimental distances from geometry

For effects `i` and `j`, this view applies the fixed contrast
`c = e_i - e_j` to the selected geometry component `G`: \$\$d\_{ij} =
c^T G c = G\_{ii} + G\_{jj} - 2G\_{ij}.\$\$ With cross-partition
geometry this is the signed crossvalidated squared Euclidean distance,
or squared Mahalanobis distance when the plan carries a fixed neural
metric. It is not `1 - Pearson correlation`: correlation distance
requires data-dependent diagonal normalization and is outside this
linear view.

## Usage

``` r
rdm(x, ...)

# Default S3 method
rdm(
  x,
  component = c("total", "coherent", "configuration"),
  pairs = NULL,
  normalize = NULL,
  ...
)
```

## Arguments

- x:

  An `effect_geometry_plan` or a complete effect form carrying the
  symmetric self-form capability.

- ...:

  Passed to the method. The generic dispatches on `x`: a geometry plan
  or a complete effect form is read by the default method documented
  here, and an `effect_population_result` from
  [`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
  by the group-level method in
  [population_views](https://bbuchsbaum.github.io/crossform/reference/population_views.md).

- component:

  One of `total`, `coherent`, or `configuration`.

- pairs:

  Optional two-column matrix selecting the effect pairs to report, by
  effect name or index. The default reports every unordered pair.
  Selected pairs execute without materializing the remaining geometry:
  the RDM is a view, not a mandatory intermediate object.

- normalize:

  Must be omitted. Correlation-style diagonal normalization of a signed
  cross-generalized form is refused: crossvalidated diagonals can be
  zero or negative, so `1 - r` here is not conventional Pearson
  distance. The boundary is documented in the correlation-distance
  policy.

## Value

An `effect_rdm_view`. `$values` has one row per spatial measurement and
one column per requested experimental pair, `$pairs` is the
`left`/`right` table naming those columns, and `$component`, `$index`,
and `$receipt` record what was read. Cross-generalized distances may be
negative.

## Details

The returned values are point estimates. For equal-weight all-partition-
pairs geometry with a common fixed metric,
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
can construct a separate analytic within-measurement covariance law when
the plan was built from
[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md)
and its residual error channel is still available. Pair rows that share
a partition are not independent replicates for a spread-across-pairs
standard error.

## Structure

The distances are one measurement-by-pair matrix; the other elements
name its axes and record what was read.

- `$values`: one row per spatial measurement, one column per requested
  pair, in `$pairs` row order.

- `$pairs`: a data frame whose `left` and `right` columns name the two
  effects behind each column of `$values`.

- `$component`: the geometry component the distances were taken from.

- `$index`: the measurement identifiers, one per row of `$values`,
  carried from the frame's `$index$measurement`.

- `$receipt`: the execution receipt for the run that produced the
  values.

Any element not listed here is internal and may change.

## Refusals

`normalize` signals an `effect_capability_refusal` with capability
`"guaranteed_psd"`, and a rectangular cross-axis plan or a non-symmetric
form signals capability `"symmetric_self_form"`, both in namespace
`"geometry_views"`. Branch on them with
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md).

## See also

[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) to
regress model RDMs on these distances,
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
for a single contrast, and
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
for the admitted analytic uncertainty law.

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
[`query_geometry()`](https://bbuchsbaum.github.io/crossform/reference/query_geometry.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md),
[`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md)

## Examples

``` r
# Three conditions over two regions, generalizing across two runs.
domain <- abstract_domain(4, id = "rdm-example")
run1 <- rbind(
  face = c(1, 0.2, 0, 0), house = c(0, 1, 0.1, 0), tool = c(0, 0, 1, 0.3)
)
run2 <- rbind(
  face = c(0.9, 0.3, 0, 0), house = c(0.1, 0.9, 0, 0), tool = c(0, 0.1, 1.1, 0.2)
)
relation <- relation(list(run1 = run1, run2 = run2), domain = domain)
plan <- plan_geometry(
  relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
  cross_partitions(
    relation, independence = "independent", generalizes_over = "run"
  )
)

# All three unordered pairs, one column each.
distances <- rdm(plan)
distances
#> <effect_rdm_view>
#>   measurements: 2
#>  measurement face - house face - tool house - tool
#>           v1         0.64        0.47        0.400
#>           it         0.00        0.58        0.525
distances$pairs
#>    left right
#> 1  face house
#> 2  face  tool
#> 3 house  tool

# Selecting pairs is a narrower view, not a post-hoc subset: the remaining
# geometry is never computed.
rdm(plan, pairs = cbind("face", "house"))$values
#>      face - house
#> [1,]         0.64
#> [2,]         0.00

# Correlation-style normalization is refused, because crossvalidated
# diagonals can be zero or negative.
refusal <- catch_refusal(rdm(plan, normalize = "correlation"))
refusal$capability
#> [1] "guaranteed_psd"
```
