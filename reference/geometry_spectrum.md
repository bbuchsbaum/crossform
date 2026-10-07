# Read the signed eigenvalue spectrum of cross-generalized geometry

Unlike the linear views, eigenvalues are not additive across components:
the spectra of the `coherent` and `configuration` components do not sum
to the spectrum of `total`, even though the underlying matrices do.
Compare spectra across components only as separate decompositions.

## Usage

``` r
geometry_spectrum(
  x,
  component = c("total", "coherent", "configuration"),
  row_block = 1024L
)
```

## Arguments

- x:

  A complete effect form carrying the symmetric self-form capability, as
  returned by
  [`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md).
  A query-only view has no geometry to decompose.

- component:

  Geometry component to decompose.

- row_block:

  Positive number of measurement rows read per block.

## Value

An `effect_spectrum_view`. `$values` has one row per measurement and one
column per eigenvalue (`root1` largest), with `$component`, `$index`,
`$receipt`, and `$indefinite_estimates_preserved = TRUE` recording that
negative eigenvalues are never truncated at zero.

## Structure

One signed spectrum per spatial measurement, with the labels that say
what was decomposed.

- `$values`: one row per measurement, one column per eigenvalue, ordered
  `root1` (largest) through `rootq`. Negative roots are retained.

- `$component`: the geometry component that was decomposed.

- `$index`: the measurement identifiers, one per row of `$values`,
  carried from the decomposed form.

- `$indefinite_estimates_preserved`: always `TRUE`, recording that no
  eigenvalue was truncated at zero.

- `$receipt`: the execution receipt of the form that was decomposed.

Any element not listed here is internal and may change.

## Refusal

Passing a query-first `effect_geometry_plan` signals an
`effect_capability_refusal` with capability `"complete_geometry"` in
namespace `"geometry_views"` and remedy `materialize_geometry(x)`; a
non-symmetric form signals capability `"symmetric_self_form"`. See
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md).

## See also

[`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md),
which produces the complete geometry this view requires, and
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) for
the linear distance view.

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
domain <- abstract_domain(4, id = "spectrum-example")
relation <- relation(
  list(run1 = rbind(a = c(1, 0, 2, 1), b = c(0, 1, 1, 0)),
       run2 = rbind(a = c(1.1, 0.1, 1.9, 0.8), b = c(0.1, 0.9, 1.2, 0.2))),
  domain = domain
)
geometry <- materialize_geometry(plan_geometry(
  relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
  cross_partitions(relation, independence = "independent")
))

# Signed eigenvalues per region, largest first. Small negative roots are
# retained: a cross-generalized form is not guaranteed positive.
spectrum <- geometry_spectrum(geometry)
spectrum
#> <effect_spectrum_view>
#>   measurements: 2
#>  measurement  root1   root2
#>           v1 0.5707 0.42929
#>           it 2.8600 0.03999
as.data.frame(spectrum)
#>   measurement     root1      root2
#> 1          v1 0.5707107 0.42928932
#> 2          it 2.8600089 0.03999113

# Eigenvalues are not additive across components, so read each spectrum as
# its own decomposition rather than summing them.
geometry_spectrum(geometry, component = "coherent")$values
#>         root1        root2
#> [1,] 0.551134 -0.001134025
#> [2,] 2.389712 -0.014711504
```
