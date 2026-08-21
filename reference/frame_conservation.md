# Diagnose local-to-global conservation of a compiled frame

Under a conservative frame every domain feature carries total weight
mass one, so local `total` geometries sum exactly to the global
geometry: \$\$\sum_x G_x^{\mathrm{total}} =
G\_\Omega^{\mathrm{total}}.\$\$ Overlapping neighborhoods under
`normalization = "local"` double-count shared features, so their local
values are not contributions to a whole-brain quantity. Two
preconditions matter when checking the law: the global comparator must
be the unnormalized operator (`whole_brain("none")`), because default
local normalization divides by the feature count; and the law covers the
`total` component only — `coherent` is defined by each measurement's own
weighted common mode, and local coherent values do not sum to the global
coherent component.

## Usage

``` r
frame_conservation(x, tolerance = 1e-10)
```

## Arguments

- x:

  A compiled `effect_frame`.

- tolerance:

  Nonnegative absolute per-feature mass tolerance.

## Value

An `effect_frame_conservation` list reporting `conserved`, the covered
`component` (`"total"`), the frame `normalization`, its
`declared_normalization` together with `metric_folded` (whether a
diagonal metric has been folded into the weights) and the `composition`
law in force (`"none"` for a declared frame, `"diagonal_metric_fold"`
for a folded one), the maximum per-feature deviation from the conserving
`reference_mass`, and the per-feature mass vector. `reference_mass` is
one for a declared frame and the folded metric diagonal for a
metric-folded one, because that frame's global comparator is read under
the same metric. A
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md)
additionally reports `members`, one row per family member giving its
`alpha`, its measurement count, and the deviation of its own per-feature
mass from that `alpha` – the block-by-block half of the law, which the
whole family conserving does not imply.

## See also

[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md)
and the normalization argument of
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md),
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md),
[`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md),
and
[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md).

Other neural domains and frames:
[`abstract_domain()`](https://bbuchsbaum.github.io/crossform/reference/abstract_domain.md),
[`additive_frame()`](https://bbuchsbaum.github.io/crossform/reference/additive_frame.md),
[`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md),
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md),
[`edge_frame()`](https://bbuchsbaum.github.io/crossform/reference/edge_frame.md),
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md),
[`measurement_frame()`](https://bbuchsbaum.github.io/crossform/reference/measurement_frame.md),
[`neuroim2_searchlights()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_searchlights.md),
[`neuroim2_volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_volume_domain.md),
[`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md),
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md),
[`volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/volume_domain.md),
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md),
[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)

## Examples

``` r
domain <- abstract_domain(4, id = "conservation-example")

# A conservative frame gives every feature total weight mass one, so local
# `total` geometries sum exactly to the global one.
conservative <- compile_frame(voxelwise(), domain)
frame_conservation(conservative)$conserved
#> [1] TRUE

# Under "local" each measurement is rescaled instead, and the per-feature
# mass tells you by how much the accounting is off.
report <- frame_conservation(compile_frame(whole_brain(), domain))
report$conserved
#> [1] FALSE
report$max_deviation
#> [1] 0.75
```
