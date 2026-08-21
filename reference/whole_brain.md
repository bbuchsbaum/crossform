# Specify a whole-brain additive measurement

Declares a single measurement covering every neural feature in the
domain. Use `normalization = "none"` when this is the global comparator
against which local measurements are checked, because the default
averages instead of summing.

## Usage

``` r
whole_brain(normalization = "local")
```

## Arguments

- normalization:

  Explicit frame normalization.

## Value

An `effect_frame_spec` with `$kind` `"whole_brain"` and the requested
`$normalization`. Pass it to
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md).

## See also

[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md),
[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md),
and
[`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md)
for several measurements instead of one.

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
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md)

## Examples

``` r
domain <- abstract_domain(6L, id = "whole-brain-example")

# A single measurement row spanning the whole domain.
frame <- compile_frame(whole_brain(), domain)
dim(frame$weights)
#> [1] 1 6

# The default "local" normalization averages over features (mass 1); the
# unnormalized operator is the one that sums, and is the correct global
# comparator when checking local-to-global conservation.
as.numeric(Matrix::rowSums(frame$weights))
#> [1] 1
as.numeric(Matrix::rowSums(
  compile_frame(whole_brain("none"), domain)$weights
))
#> [1] 6
```
