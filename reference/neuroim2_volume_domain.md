# Construct a crossform domain from a neuroim2 volume mask

This conditional adapter records the mask's stable full-volume indices,
physical coordinates, spacing, and full neuroim2 spatial metadata. It
does not extract response data or import neuroim2 result types.

## Usage

``` r
neuroim2_volume_domain(mask, id = "neuroim2-volume")
```

## Arguments

- mask:

  A three-dimensional neuroim2 `NeuroVol` mask.

- id:

  Stable domain identity.

## Value

An `effect_domain` with `$kind` `"volume"`, `$feature_ids` giving the
mask's full-volume indices, `$coordinates` in millimeters, and
`$metadata` carrying `dim`, `spacing`, `voxel` indices, the logical
`mask`, and a `neuroim2_space_sha256` hash of the full neuroim2 space.

## See also

[`volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/volume_domain.md)
for the same object from a plain array,
[`neuroim2_searchlights()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_searchlights.md)
for matching neighborhoods, and
[`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md)
to write compact results back out.

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
[`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md),
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md),
[`volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/volume_domain.md),
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md),
[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)

## Examples

``` r
if (requireNamespace("neuroim2", quietly = TRUE) &&
    utils::packageVersion("neuroim2") >= "0.19.0") {
  values <- array(FALSE, c(5L, 5L, 4L))
  values[2:4, 2:4, 2:3] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(
    values, neuroim2::NeuroSpace(c(5L, 5L, 4L), spacing = c(3, 3, 3))
  )
  domain <- neuroim2_volume_domain(mask, id = "subject-mask")

  # Compact features, addressed by their stable full-volume indices, so a
  # result vector can always be placed back in the original array.
  print(domain$n_features)
  print(identical(domain$feature_ids, which(values)))

  # Physical spacing and a hash of the full neuroim2 space are recorded;
  # any frame built later must agree with that geometry.
  print(domain$metadata$spacing)
  print(substr(domain$metadata$neuroim2_space_sha256, 1, 24))
}
#> [1] 18
#> [1] TRUE
#> [1] 3 3 3
#> [1] "sha256:ea44fad00b51f3d39"
```
