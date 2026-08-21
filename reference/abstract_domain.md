# Construct an abstract neural feature domain

Use this when the neural features are a plain ordered set — surface
nodes, parcels, electrodes, or columns of a beta matrix — and any
geometry you have is supplied directly as a coordinate matrix rather
than read from a volume.

## Usage

``` r
abstract_domain(
  n_features,
  coordinates = NULL,
  feature_ids = NULL,
  id = "abstract",
  coordinate_units = "arbitrary"
)
```

## Arguments

- n_features:

  Positive neural feature count.

- coordinates:

  Optional finite feature-by-coordinate matrix used by spatial frame
  builders.

- feature_ids:

  Optional unique feature identifiers.

- id:

  Stable nonempty domain identity.

- coordinate_units:

  One coordinate unit or one per coordinate axis.

## Value

An `effect_domain`: a list with `$id`, `$kind` (`"abstract"`),
`$n_features`, `$feature_ids`, the optional `$coordinates` matrix,
`$coordinate_units`, a `$geometry_signature`, a compact `$reference`
that other objects store instead of the full domain, and `$metadata`.

## See also

[`volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/volume_domain.md)
for a mask-derived domain and
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md)
to turn a domain into measurements.

Other neural domains and frames:
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
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md),
[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)

## Examples

``` r
# Twelve features with no geometry: enough for region or whole-brain scopes.
domain <- abstract_domain(12L, id = "from-observations:v1")
domain$n_features
#> [1] 12
domain$kind
#> [1] "abstract"

# Supply coordinates when a spatial scope such as searchlights() will be
# compiled; the radius is then read in these coordinate units.
grid <- as.matrix(expand.grid(x = 1:3, y = 1:4))
placed <- abstract_domain(nrow(grid), coordinates = grid, id = "planar-grid")
dim(compile_frame(searchlights(1.5), placed)$weights)
#> [1] 12 12
```
