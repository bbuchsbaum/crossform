# Specify region measurements

Declares one measurement per distinct label, so an atlas or region-label
vector becomes a set of non-overlapping regional measurements.

## Usage

``` r
regions(labels, normalization = "local")
```

## Arguments

- labels:

  One region label per neural feature. Missing labels are excluded
  unless conservative normalization is requested.

- normalization:

  Explicit frame normalization.

## Value

An `effect_frame_spec` with `$kind` `"regions"`, the supplied `$labels`,
and `$normalization`. Pass it to
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md).

## See also

[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md),
and
[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)
for the single-region case.

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
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md),
[`volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/volume_domain.md),
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md),
[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)

## Examples

``` r
domain <- abstract_domain(12L, id = "region-example")
labels <- rep(paste0("roi-", 1:3), each = 4L)

# Measurement order follows first appearance of each label.
frame <- compile_frame(regions(labels), domain)
frame$index$measurement
#> [1] "roi-1" "roi-2" "roi-3"

# Unlabeled features are dropped: with "none" the row sums are the member
# counts, so the effect of an NA label is visible.
as.numeric(Matrix::rowSums(
  compile_frame(regions(labels, normalization = "none"), domain)$weights
))
#> [1] 4 4 4
labels[[1L]] <- NA
as.numeric(Matrix::rowSums(
  compile_frame(regions(labels, normalization = "none"), domain)$weights
))
#> [1] 3 4 4
```
