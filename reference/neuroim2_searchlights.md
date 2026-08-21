# Compile neuroim2 searchlight indices into a crossform frame

The function calls only
[`neuroim2::searchlight_indices()`](https://bbuchsbaum.github.io/neuroim2/reference/searchlight_indices.html)
and maps its stable full-volume indices to the ordered compact feature
columns of `domain`.

## Usage

``` r
neuroim2_searchlights(
  mask,
  radius,
  domain = NULL,
  normalization = "local",
  nonzero = TRUE,
  weights = NULL
)
```

## Arguments

- mask:

  A three-dimensional neuroim2 `NeuroVol` mask.

- radius:

  Positive spherical radius in millimeters. Several radii request a
  multiscale family, one member frame per radius.

- domain:

  An exact domain from
  [`neuroim2_volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_volume_domain.md).
  When omitted, it is constructed from `mask`.

- normalization:

  One of `none`, `local`, or `conservative`. Several radii admit only
  `conservative`.

- nonzero:

  Passed to
  [`neuroim2::searchlight_indices()`](https://bbuchsbaum.github.io/neuroim2/reference/searchlight_indices.html);
  version 0.1 requires `TRUE` so every member belongs to the compact
  domain.

- weights:

  Family weights for a multiscale request: one positive weight per
  radius, summing to one, matched to the radii in order or by the
  `"radius-<r>"` names. `NULL` (the default) weights the radii equally.
  A single radius is one frame with no budget to divide, so `weights` is
  refused there rather than ignored.

## Value

An `effect_frame` whose `$weights` are the sparse measurement-by-feature
operator, with `$index$measurement` holding the full-volume center
indices, `$normalization`, a `$specification` recording the radius and
the pinned `upstream_commit`, and a `$support_index`. Several radii
return a
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md)
instead: its `$index` carries one row per measurement with that row's
`family`, `node`, `scale`, `center`, and `alpha`.

## Details

Several radii request one conservative frame per radius, stacked into a
multiscale
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md);
see *Multiscale families*.

## Multiscale families

`neuroim2_searchlights(mask, c(4, 8), normalization = "conservative")`
builds one conservative frame per radius from the same neuroim2
neighborhoods, names them `"radius-4"` and `"radius-8"`, and stacks them
with
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md)
under family weights `weights` (equal by default). Only conservative
normalization is admitted, for the reason
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md)
gives: locally normalized values are not contributions to any total, so
a family of them has no budget for `weights` to divide.

Per-scale totals of such a family are \\\alpha_s G\_\Omega\\ by
construction, so a total-energy-by-scale panel reports the analyst's own
`weights`, not the data. Only the coherent share of each block's fixed
budget varies informatively with scale
(`design/conservative-geometry-contract.md` sections 3.1 and 3.2).

## See also

[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md)
plus
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md)
for the built-in neighborhood builder,
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md)
for the family several radii compile to,
[`neuroim2_volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_volume_domain.md)
for the domain, and
[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md)
to check the normalization you chose.

Other neural domains and frames:
[`abstract_domain()`](https://bbuchsbaum.github.io/crossform/reference/abstract_domain.md),
[`additive_frame()`](https://bbuchsbaum.github.io/crossform/reference/additive_frame.md),
[`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md),
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md),
[`edge_frame()`](https://bbuchsbaum.github.io/crossform/reference/edge_frame.md),
[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md),
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md),
[`measurement_frame()`](https://bbuchsbaum.github.io/crossform/reference/measurement_frame.md),
[`neuroim2_volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_volume_domain.md),
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
  domain <- neuroim2_volume_domain(mask)

  # Neighborhoods come from neuroim2; crossform only maps their full-volume
  # member indices onto the ordered compact feature columns.
  frame <- neuroim2_searchlights(mask, radius = 4, domain = domain)
  print(dim(frame$weights))
  print(identical(frame$index$measurement, domain$feature_ids))

  # The pinned upstream geometry is part of the frame's specification.
  print(frame$specification$upstream_commit)

  # Several radii stack into a conservative family, one member per radius.
  family <- neuroim2_searchlights(mask, c(4, 6), domain = domain,
    normalization = "conservative", weights = c(0.4, 0.6))
  print(unique(family$index[, c("family", "scale", "alpha")]))

  # Each block carries exactly its weight, so per-scale energy is the
  # `weights` vector and only the coherent share is a finding.
  print(frame_conservation(family)$members)

  # A mask whose geometry differs from the declared domain is refused.
  moved <- neuroim2::LogicalNeuroVol(
    values, neuroim2::NeuroSpace(c(5L, 5L, 4L), spacing = c(2, 2, 4))
  )
  print(try(neuroim2_searchlights(moved, radius = 4, domain = domain)))
}
#> [1] 18 18
#> [1] TRUE
#> [1] "77b1ddb"
#>      family scale alpha
#> 1  radius-4     4   0.4
#> 19 radius-6     6   0.6
#>     family alpha measurements max_deviation conserved
#> 1 radius-4   0.4           18  5.551115e-17      TRUE
#> 2 radius-6   0.6           18  1.110223e-16      TRUE
#> Error : `mask` and `domain` have different volume geometry, so a compact index in one does not name the same voxel in the other. The mask is 5 x 5 x 4 with spacing 2 x 2 x 4 and 18 in-mask voxels; domain `neuroim2-volume` is 5 x 5 x 4 with spacing 3 x 3 x 3 and 18 features. Pass the mask the domain was built from.
#> [1] "Error : `mask` and `domain` have different volume geometry, so a compact index in one does not name the same voxel in the other. The mask is 5 x 5 x 4 with spacing 2 x 2 x 4 and 18 in-mask voxels; domain `neuroim2-volume` is 5 x 5 x 4 with spacing 3 x 3 x 3 and 18 features. Pass the mask the domain was built from.\n"
#> attr(,"class")
#> [1] "try-error"
#> attr(,"condition")
#> <effect_contract_error>
#>   `mask` and `domain` have different volume geometry, so a compact index in
#>   one does not name the same voxel in the other. The mask is 5 x 5 x 4 with
#>   spacing 2 x 2 x 4 and 18 in-mask voxels; domain `neuroim2-volume` is 5 x
#>   5 x 4 with spacing 3 x 3 x 3 and 18 features. Pass the mask the domain
#>   was built from.
```
