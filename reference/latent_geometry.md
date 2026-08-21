# Construct the latent PSD descriptive layer of a signed geometry

Crossvalidated contributions are **signed**, and the arithmetic that
treats a node's value as part of a nonnegative whole is invalid on them
(`design/conservative-geometry-contract.md` section 6): the implied
shares can leave `[0, 1]`, and clipping to nonnegativity both
reintroduces the bias the cross-partition estimator removes and destroys
conservation. Effective counts, cumulative-contribution curves and
fractions are legal only on a declared nonnegative projection of the
estimates, and `latent_geometry()` is the named operation that builds
one.

## Usage

``` r
latent_geometry(
  x,
  method = c("psd_projection"),
  component = c("total", "coherent", "configuration"),
  row_block = 1024L
)
```

## Arguments

- x:

  A complete symmetric effect form from
  [`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md),
  or the signed `effect_spectrum_view`
  [`geometry_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/geometry_spectrum.md)
  reads from one. A query readout has already been contracted against
  fixed weights and carries no spectrum, so it is refused.

- method:

  The projection to apply, from the latent layer's closed set.
  `"psd_projection"` truncates the symmetric form's eigenvalues at zero.
  `"nearest_psd"` is a declared member of the set and is not
  implemented; it refuses rather than falling back, because it moves
  different mass.

- component:

  Geometry component to project. Fixed by the view when `x` is an
  `effect_spectrum_view`.

- row_block:

  Positive number of measurement rows read per block when `x` is a form.

## Value

An `effect_latent_geometry`.

## Details

The projection is biased and does not conserve. That is not a defect to
be worked around; it is what makes the layer's arithmetic well defined,
and it is why the moved mass is reported per measurement rather than
absorbed. Read the layer descriptively – which modes carry a form's
nonnegative mass and how concentrated they are – and read the signed
geometry for anything that has to be unbiased or conserved.

## Structure

One nonnegative spectrum per spatial measurement, with the functionals
that only a nonnegative partition admits and the record of what the
projection cost.

- `$spectrum`: the projected spectrum, one row per measurement, one
  column per root, ordered `root1` (largest) through `rootq`. Every
  entry is nonnegative.

- `$cumulative`: the cumulative contribution curve \\C(k) = \sum\_{i \le
  k} \lambda_i / \sum_i \lambda_i\\, columns `C1` through `Cq`, so `Cq`
  is `1` wherever the curve is defined.

- `$n_eff`: the participation ratio \\(\sum_i \lambda_i)^2 / \sum_i
  \lambda_i^2\\ of the projected spectrum, one per measurement.

- `$moved_mass`: the absolute mass the projection removed from each
  measurement – the sum of the magnitudes of its negative roots. Zero
  for a measurement whose source spectrum was already nonnegative.

- `$moved_share`: that mass as a fraction of the source spectrum's total
  absolute mass, so a measurement that moved a lot in a large form is
  not confused with one that moved a little in a small one.

- `$component`, `$method`: what was projected, and by which named
  operator.

- `$index`: the measurement identifiers, carried from the source.

- `$projection`: the projection receipt – the operator, the
  per-measurement moved mass and share, the counts of clipped and masked
  measurements, and the source identity it was derived from.

- `$receipt`: the execution receipt. Its `$scientific_plan_id` is
  derived from the source's and the projection's name, so a latent layer
  never shares an identity with the signed estimates behind it, and its
  `$task_partition_id` ends in `+psd_projection`.

Any element not listed here is internal and may change.

## Masking

A measurement whose projected spectrum sums to zero – every root at or
below zero – has no nonnegative partition, so its `$n_eff` and its whole
`$cumulative` row are `NA` rather than `0/0`. `$projection$masked`
counts them. This is the discipline the coherence fraction already uses:
a number that was not earned is withheld, not clamped.

The guard is exact (`> 0`), not a tolerance, and it is the same guard
every per-node fraction in the package uses. A form that is zero only to
within rounding therefore returns `$n_eff` and `$cumulative` computed
from that rounding: a row of roots at `1e-17` gives the same `$n_eff` as
a row a hundred million times larger, because the ratio is scale free.
`$moved_share` will read near `0.5` for such a row and `$spectrum` will
show the scale, so the object says what happened – but the effective
count of a form that is numerically zero is noise, and no guard here
will tell you so. Whether these guards should take a relative tolerance
is an open contract decision (`design/conservative-geometry-contract.md`
§11.4, gap G3, where a singleton scale's aggregated configuration was
measured at `-2.8e-17`); this layer inherits it rather than making a
private one.

`$moved_share` is the one quantity that is *not* masked when its
denominator vanishes: a form with no absolute mass at all moved none of
it, so the share is `0`. `NA` there would read as "we could not tell you
how much was clipped" about a measurement where the answer is known.

## Refusal

A query-first `effect_geometry_plan` signals an
`effect_capability_refusal` with capability `"complete_geometry"`; a
rectangular or non-symmetric form signals `"symmetric_self_form"`; a
query readout
([`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md),
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md),
or a
[`query_geometry()`](https://bbuchsbaum.github.io/crossform/reference/query_geometry.md)
view) signals `"latent_projection_source"`, because clamping a
contracted value is a per-node total clamp and a different projection.
`method = "nearest_psd"` signals `"nearest_psd_projection"`. See
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md).

## See also

[`geometry_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/geometry_spectrum.md)
for the signed spectrum this layer projects, which never truncates a
negative root, and
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md)
for the additive ledger reading that stays on the signed layer.

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
domain <- abstract_domain(4, id = "latent-example")
relation <- relation(
  list(
    run1 = rbind(face = c(1, 0, 2, 1), house = c(0, 1, 1, 0),
      tool = c(0.5, -1, 0, 1)),
    run2 = rbind(face = c(0.9, 0.2, 1.7, 0.6), house = c(0.2, 0.8, 1.3, 0.1),
      tool = c(-0.4, 1, 0.2, -0.8))
  ),
  domain = domain
)
geometry <- materialize_geometry(plan_geometry(
  relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
  cross_partitions(relation, independence = "independent")
))

# The signed spectrum keeps its negative root; both regions have one.
geometry_spectrum(geometry)$values
#>          root1      root2      root3
#> [1,] 0.5283291 0.32918218 -0.6075113
#> [2,] 2.6177260 0.03822887 -0.4059548

# The latent layer truncates it and says what that cost, per region.
latent <- latent_geometry(geometry)
latent
#> <effect_latent_geometry>
#>   measurements: 2
#>   component:    total
#>   projection:   psd_projection (eigenvalue truncation at zero)
#>   moved_mass:   1.013 moved, 2 of 2 measurements clipped, max share 0.415
#>   n_eff:        median 1.46 (range 1.03 to 1.9)
#>   reading:      latent descriptive layer; not for inference
#>  measurement n_eff moved_mass moved_share  root1   root2 root3
#>           v1 1.898     0.6075      0.4147 0.5283 0.32918     0
#>           it 1.029     0.4060      0.1326 2.6177 0.03823     0
#>   next:         x$cumulative, x$projection
latent$cumulative
#>             C1 C2 C3
#> [1,] 0.6161191  1  1
#> [2,] 0.9856064  1  1
latent$projection
#> <effect_latent_projection_receipt>
#>   method:     psd_projection
#>   operator:   eigenvalue truncation at zero
#>   component:  total
#>   clipped:    2 of 2 measurements
#>   moved_mass: 1.013 total, max share 0.415
#>   masked:     none
#>   source:     geometry-sha256:21c27c763dd0...

# A contracted readout has no spectrum to project, and clamping it would be
# a per-node total clamp: a different operator moving different mass.
effect <- contrast_energy(geometry, c(face = 1, house = -1, tool = 0))
catch_refusal(latent_geometry(effect))$capability
#> [1] "latent_projection_source"
```
