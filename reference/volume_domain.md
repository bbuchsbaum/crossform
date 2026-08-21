# Construct a native volumetric neural feature domain

Use this when features are the in-mask voxels of a 3D volume. The
compact feature order follows the full-volume index order of the mask,
so results can always be written back to the original array positions.

## Usage

``` r
volume_domain(
  mask,
  spacing = c(1, 1, 1),
  id = "native-volume",
  coordinate_units = "mm",
  metadata = list()
)
```

## Arguments

- mask:

  A three-dimensional logical or numeric mask. Finite nonzero entries
  are included.

- spacing:

  Three positive finite voxel spacings.

- id:

  Stable domain identity.

- coordinate_units:

  One physical coordinate unit or one per axis.

- metadata:

  Optional uniquely named list of extra facts about where the geometry
  came from, recorded alongside the grid facts the constructor derives.
  See *Provider metadata*.

## Value

An `effect_domain` with `$kind` `"volume"`, `$feature_ids` giving the
stable full-volume indices of the included voxels, `$coordinates`
holding their physical positions, `$metadata` carrying `dim`, `spacing`,
`voxel` indices and the logical `mask` followed by anything `metadata`
added, plus the usual `$geometry_signature` and `$reference`.

## Provider metadata

A domain built by an adapter usually knows something about the geometry
that the array itself does not carry — the native header it came from,
the file, the transform it was resampled under. `metadata` is where that
goes, and it is not decoration: the domain's `$geometry_signature`
covers it, so two domains that agree on every voxel but disagree about
their provenance are correctly *different* domains, and anything holding
a compact result vector against one of them refuses the other.

That is the point.
[`neuroim2_volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_volume_domain.md)
records the hash of the full `neuroim2` space this way, which is what
makes writing a result back to voxels safe. The grid facts the
constructor derives itself — `dim`, `spacing`, `voxel`, `mask` — cannot
be overridden.

## See also

[`abstract_domain()`](https://bbuchsbaum.github.io/crossform/reference/abstract_domain.md)
for non-volumetric features,
[`neuroim2_volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_volume_domain.md)
to build the same object from a `NeuroVol`, and
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md)
to place searchlights on it.

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
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md),
[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)

## Examples

``` r
# A 5 x 5 x 3 volume with one voxel excluded from the mask.
mask <- array(TRUE, c(5L, 5L, 3L))
mask[1, 1, 1] <- FALSE
domain <- volume_domain(mask, spacing = c(3, 3, 3), id = "example-volume")
domain$n_features
#> [1] 74

# Coordinates are millimeters, so a searchlight radius is physical too.
utils::head(domain$coordinates, 3)
#>      [,1] [,2] [,3]
#> [1,]    3    0    0
#> [2,]    6    0    0
#> [3,]    9    0    0
nrow(compile_frame(searchlights(4), domain)$weights)
#> [1] 74

# feature_ids are full-volume indices, which is how compact results are
# written back into the original array.
utils::head(domain$feature_ids, 3)
#> [1] 2 3 4

# A provider records where the geometry came from, and the record is part
# of what the domain is: the same voxels under a different provenance are a
# different domain, not the same one.
stamped <- volume_domain(mask, spacing = c(3, 3, 3), id = "example-volume",
  metadata = list(source_file = "sub-01_mask.nii.gz"))
stamped$metadata$source_file
#> [1] "sub-01_mask.nii.gz"
identical(stamped$reference, domain$reference)
#> [1] FALSE
```
