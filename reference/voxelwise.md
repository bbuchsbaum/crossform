# Specify point measurements

Declares one measurement per neural feature. This is the finest scope
and, under the default conservative normalization, the only spatial
scope whose local values sum exactly to the whole-brain value.

## Usage

``` r
voxelwise(normalization = "conservative")
```

## Arguments

- normalization:

  Explicit frame normalization.

## Value

An `effect_frame_spec` with `$kind` `"voxels"` and the requested
`$normalization`. Pass it to
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md)
to obtain measurements.

## See also

[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md),
and
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md),
[`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md),
or
[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)
for coarser scopes.

Other neural domains and frames:
[`abstract_domain()`](https://bbuchsbaum.github.io/crossform/reference/abstract_domain.md),
[`additive_frame()`](https://bbuchsbaum.github.io/crossform/reference/additive_frame.md),
[`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md),
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md),
[`edge_frame()`](https://bbuchsbaum.github.io/crossform/reference/edge_frame.md),
[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md),
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md),
[`measurement_frame()`](https://bbuchsbaum.github.io/crossform/reference/measurement_frame.md),
[`neuroim2_searchlights()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_searchlights.md),
[`neuroim2_volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_volume_domain.md),
[`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md),
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md),
[`volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/volume_domain.md),
[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)

## Examples

``` r
domain <- abstract_domain(4L, id = "voxelwise-example")

# One measurement per feature: the operator is the identity.
frame <- compile_frame(voxelwise(), domain)
dim(frame$weights)
#> [1] 4 4

# Conservative normalization is the default here, so local totals add up to
# the whole-brain total exactly.
frame_conservation(frame)$conserved
#> [1] TRUE
```
