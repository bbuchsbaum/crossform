# Specify Euclidean searchlights

The default `"local"` normalization matches the conventional
center-assigned searchlight map, in which overlapping neighborhoods
repeatedly count shared features. Local values under this default cannot
be summed into a whole-brain quantity. Request
`normalization = "conservative"` for exact local-to-global accounting of
the `total` component, and check any frame with
[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md).

## Usage

``` r
searchlights(radius, normalization = "local", weights = NULL)
```

## Arguments

- radius:

  Positive radius in domain coordinate units. Several radii request a
  multiscale family, one member frame per radius.

- normalization:

  Explicit frame normalization. Several radii admit only
  `"conservative"`.

- weights:

  Family weights for a multiscale request: one positive weight per
  radius, summing to one, matched to the radii in order or by the
  `"radius-<r>"` names. `NULL` (the default) weights the radii equally.
  A single radius is one frame with no budget to divide, so `weights` is
  refused there rather than ignored.

## Value

An `effect_frame_spec`. With one radius its `$kind` is `"searchlights"`
and it carries the requested `$radius` and `$normalization`, exactly as
before. With several its `$kind` is `"searchlight_family"` and it
carries the `$radius` vector and the `$weights` that will be applied.
Pass either to
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md).

## Details

Several radii request one conservative frame per radius, stacked into a
multiscale
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md);
see *Multiscale families*.

## Multiscale families

`compile_frame(searchlights(c(4, 8, 12), "conservative"), domain)`
returns a
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md):
one conservative member frame per radius, named `"radius-4"`,
`"radius-8"`, `"radius-12"`, stacked with family weights `weights`
(equal by default). Every row of the compiled frame's `$index` carries
its own `family`, `scale` (the radius it came from), `center`, and
`alpha`, so a result can be grouped by scale after the fact.

Only `normalization = "conservative"` admits several radii. Locally
normalized values are not contributions to any total, so a family of
them has no budget for `weights` to divide and the per-scale law below
is undefined (`design/conservative-geometry-contract.md` section 3.1).

**What a multiscale family can and cannot show.** Because every member
is column normalized and the weights sum to one, the family conserves
block by block: the `total` component summed over the rows of scale `s`
is exactly \\\alpha_s G\_\Omega\\, the scale's weight times the
whole-brain total, whatever the data say. A panel of *total energy by
scale* is therefore a plot of the analyst's own `weights` vector and is
not a finding. What does vary informatively with scale is the split of
each block's fixed budget into coherent and configuration parts: the
coherent share is invariant to `alpha` and is the scale-resolved
quantity a multiscale family exists to report
(`design/conservative-geometry-contract.md` sections 3.1 and 3.2).

## See also

[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md),
[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md)
to check the normalization you chose,
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md)
for the family a multiscale request compiles to, and
[`neuroim2_searchlights()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_searchlights.md)
for neighborhoods built by neuroim2 instead.

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
[`volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/volume_domain.md),
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md),
[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)

## Examples

``` r
grid <- as.matrix(expand.grid(x = 1:4, y = 1:4))
domain <- abstract_domain(
  nrow(grid), coordinates = grid, id = "searchlight-example"
)

# One center-assigned neighborhood per feature; the radius is read in the
# domain's own coordinate units.
local <- compile_frame(searchlights(1.5), domain)
dim(local$weights)
#> [1] 16 16

# The default "local" normalization double-counts features shared between
# overlapping neighborhoods, so its local values are not contributions to a
# whole-brain total. Ask for "conservative" when they must be.
frame_conservation(local)$conserved
#> [1] FALSE
frame_conservation(
  compile_frame(searchlights(1.5, normalization = "conservative"), domain)
)$conserved
#> [1] TRUE

# Several radii request a multiscale family instead: one conservative
# member frame per radius, each row labelled with the scale it came from.
family <- compile_frame(
  searchlights(c(1.5, 2.5), "conservative", weights = c(0.25, 0.75)), domain
)
unique(family$index[, c("family", "scale", "alpha")])
#>        family scale alpha
#> 1  radius-1.5   1.5  0.25
#> 17 radius-2.5   2.5  0.75

# Each block carries exactly its weight, whatever the data: per-scale
# energy is the `weights` vector, so only the coherent share is a finding.
frame_conservation(family)$members
#>       family alpha measurements max_deviation conserved
#> 1 radius-1.5  0.25           16  5.551115e-17      TRUE
#> 2 radius-2.5  0.75           16  2.220446e-16      TRUE
```
