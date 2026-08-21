# Compile a spatial specification against a neural domain

Voxel, searchlight, region, and whole-brain scopes all compile to the
same sparse additive-frame representation.

## Usage

``` r
compile_frame(specification, domain)
```

## Arguments

- specification:

  An additive frame specification.

- domain:

  An `effect_domain`.

## Value

An `effect_frame`: a list whose `$weights` are the sparse
measurement-by-feature operator, with `$normalization`, an `$index` data
frame naming each measurement, `$domain` (the domain reference),
`$domain_kind`, the originating `$specification`, and, for neighborhood
scopes, a `$support_index`.

## Structure

A compiled frame is the spatial operator together with the record of how
it was built.

- `$weights`: the sparse measurement-by-feature operator. Row `m` holds
  the weight each domain feature contributes to measurement `m`, in
  domain feature order, after normalization was applied.

- `$index`: one row per measurement, in `$weights` row order. Its
  `measurement` column names each measurement: domain feature
  identifiers for
  [`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md)
  and
  [`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md),
  the distinct labels in first appearance order for
  [`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md),
  and `"whole_brain"` for
  [`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md).
  Views carry these identifiers through as their `$index`.

- `$normalization`: the normalization that was applied, which is what
  [`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md)
  reports against.

- `$specification`: the `effect_frame_spec` the frame was compiled from,
  so the scope and its arguments travel with the operator.

- `$domain`: the neural domain the columns are bound to, carrying its
  `id`, `n_features`, and `feature_ids`.

Any element not listed here, including `$support_index` and the
representation flags, is internal and may change.

One specification is not one operator: a multiscale
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md)
request compiles to a
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md),
whose `$index` carries the extra per-row `family`, `node`, `scale`,
`center`, and `alpha` columns and whose `$specification` records every
member rather than one scope.

## See also

[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md),
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md),
[`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md),
[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)
for the specifications,
[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md)
to check normalization, and
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md),
which consumes the compiled frame.

Other neural domains and frames:
[`abstract_domain()`](https://bbuchsbaum.github.io/crossform/reference/abstract_domain.md),
[`additive_frame()`](https://bbuchsbaum.github.io/crossform/reference/additive_frame.md),
[`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md),
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
domain <- abstract_domain(9L, id = "compile-frame-example")

# Every scope compiles to one sparse measurement-by-feature operator; only
# the number of measurement rows differs.
scopes <- list(
  voxelwise(), regions(rep(c("a", "b", "c"), each = 3L)), whole_brain()
)
vapply(scopes, function(scope) nrow(compile_frame(scope, domain)$weights),
  integer(1))
#> [1] 9 3 1

# The compiled frame remembers how it was built, which is what downstream
# receipts record.
frame <- compile_frame(voxelwise(), domain)
frame$normalization
#> [1] "conservative"
frame$index$measurement
#> [1] 1 2 3 4 5 6 7 8 9
```
