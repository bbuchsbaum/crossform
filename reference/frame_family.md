# Combine conservative frames into one alpha-weighted family

Stacks several conservative frames over one domain into a single
compiled frame whose rows are the members' rows scaled by family weights
`alpha`. Because each member is column-normalized on its own and the
weights sum to one, the stacked columns still sum to one, so the family
is conservative and the `total` component conserves both overall and
block by block: \$\$\sum\_{x \in s} G\_{s,x} =
\alpha_s\\G\_\Omega,\qquad \sum_x G_x = G\_\Omega.\$\$

## Usage

``` r
frame_family(
  ...,
  alpha = NULL,
  normalization = "conservative",
  tolerance = 1e-12
)
```

## Arguments

- ...:

  Two or more compiled `effect_frame`s over one domain, ideally named. A
  member's name becomes its `family` identity; unnamed members are named
  `frame1`, `frame2`, and so on by position. Pass a list of frames with
  `do.call(frame_family, c(frames, list(alpha = alpha)))`.

- alpha:

  One positive family weight per member, summing to one. Named weights
  are matched to member names; `NULL` (the default) weights the members
  equally. Weights are never renormalized for you, because the per-block
  law reads against the weight actually applied.

- normalization:

  Frame normalization of the family. Only `"conservative"` is defined:
  the per-block law needs column-normalized members, and a family of
  locally normalized frames has no conserved budget to weight.

- tolerance:

  Nonnegative absolute tolerance for the two construction checks: that
  `alpha` sums to one, and that each member's columns sum to one on its
  own.

## Value

An `effect_frame` usable anywhere a compiled frame is, carrying
`$weights` (the alpha-scaled row-bind), `$index` (one row per
measurement, see *Per-row metadata*), and a `$specification` recording
every member specification together with its applied weight.

## Details

The consequence is that per-scale *energy* is fixed by `alpha` alone and
is never a finding. What varies with the data is the split of each
block's fixed budget into coherent and configuration parts. See
`design/conservative-geometry-contract.md` section 3.

## Per-row metadata

A family row must be self-describing, because its own scale and
provenance can no longer be read off a frame-wide specification.
`$index` therefore carries one row per measurement, in `$weights` row
order:

- `measurement`: the row's identity, `"<family>::<node>"`. It is unique
  across the family, which the same node label appearing at several
  scales is not, and it is what reaches a result's `$index`.

- `family`: the member the row came from.

- `node`: the row's label inside its own member, exactly as that
  member's `$index$measurement` had it.

- `scale`: the member's scale parameter – the radius for a searchlight
  member, `NA` for a member that has no scale.

- `center`: the anchor feature identifier, for members whose rows are
  anchored at a feature (points and searchlights); `NA` otherwise.

- `alpha`: the family weight applied to the row.

Join a result back to this table by `measurement` to group values by
scale, by center, or by member.

## See also

[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md)
for the members,
[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md),
which certifies a family both overall and block by block, and
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md),
whose multiscale form is the shorthand for a family of conservative
searchlight frames at several radii.

Other neural domains and frames:
[`abstract_domain()`](https://bbuchsbaum.github.io/crossform/reference/abstract_domain.md),
[`additive_frame()`](https://bbuchsbaum.github.io/crossform/reference/additive_frame.md),
[`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md),
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md),
[`edge_frame()`](https://bbuchsbaum.github.io/crossform/reference/edge_frame.md),
[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md),
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
domain <- abstract_domain(
  9L, coordinates = cbind(seq_len(9L) - 1, 0), id = "frame-family-example"
)

# Two conservative scales over one domain, weighted one quarter and three
# quarters.
family <- frame_family(
  point = compile_frame(voxelwise("conservative"), domain),
  narrow = compile_frame(searchlights(1.01, "conservative"), domain),
  alpha = c(point = 0.25, narrow = 0.75)
)
head(family$index, 3L)
#>   measurement family node scale center alpha
#> 1    point::1  point    1    NA      1  0.25
#> 2    point::2  point    2    NA      2  0.25
#> 3    point::3  point    3    NA      3  0.25

# The stacked family conserves, and each block carries exactly its alpha.
report <- frame_conservation(family)
report$conserved
#> [1] TRUE
report$members
#>   family alpha measurements max_deviation conserved
#> 1  point  0.25            9             0      TRUE
#> 2 narrow  0.75            9             0      TRUE

# Weights that do not sum to one are refused rather than renormalized: the
# per-block law reads against the weight actually applied.
refused <- try(
  frame_family(
    point = compile_frame(voxelwise("conservative"), domain),
    narrow = compile_frame(searchlights(1.01, "conservative"), domain),
    alpha = c(1, 1)
  ),
  silent = TRUE
)
conditionMessage(attr(refused, "condition"))
#> [1] "Family weights must sum to one; `alpha` sums to 2, off by 1. The per-block law reads the weight that was actually applied, so `alpha` is never renormalized for you: pass weights summing to one, for example `alpha / sum(alpha)`."
```
