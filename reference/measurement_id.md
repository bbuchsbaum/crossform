# Select measurements by identifier or by position

Wrap the measurements you pass to `at =` in
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md),
or to `measurement =` and `highlight =` in
[plot_views](https://bbuchsbaum.github.io/crossform/reference/plot_views.md),
to say whether the numbers are *identifiers* or *row positions*.

## Usage

``` r
measurement_id(x)

measurement_position(x)
```

## Arguments

- x:

  For `measurement_id()`, one or more measurement identifiers, numeric
  or character, exactly as a view's `measurement` column prints them
  (equivalently, values of the frame's `$index$measurement`). For
  `measurement_position()`, one or more whole-number row positions.

## Value

A classed selector, `crossform_measurement_id` or
`crossform_measurement_position`, holding `x`. It is resolved against
the frame or view it is passed with; an identifier the frame does not
carry, or a position outside `1..n`, is an input error.

## Details

Every compiled frame names its measurements twice. The **identifier** is
the frame's `$index$measurement`, carried into each view as `$index` and
printed in the view's `measurement` column: domain feature ids for
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md)
and
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md)
(the full-volume voxel indices of a masked
[`volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/volume_domain.md)),
region labels for
[`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md),
and `"family::node"` for a
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md).
The **position** is the row of the frame's `$weights` and of the view's
values, `1..n`: what `which.max(view$total)` or `order(view$total)`
returns.

On an
[`abstract_domain()`](https://bbuchsbaum.github.io/crossform/reference/abstract_domain.md)
with default ids, or a volume without holes, the two coincide and a bare
number works either way. When the identifiers are whole numbers that
differ from the positions – any masked volume – a bare number is
ambiguous: `117` copied from a printed view is voxel 117, but as a
position it is some other voxel. crossform does not guess. It refuses a
bare number there, and asks for one of:

- `measurement_id(117)`, or the character identifier `"117"`: the
  measurement printed as `117`.

- `measurement_position(which.max(view$total))`: the row position.

Identifiers that are not whole numbers (region names, `"family::node"`)
cannot be mistaken for positions, so bare numbers stay positions there.

## See also

Other sampling uncertainty:
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md),
[`reconstruct_evidence()`](https://bbuchsbaum.github.io/crossform/reference/reconstruct_evidence.md),
[`sampling_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/sampling_capabilities.md),
[`sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/sampling_covariance.md)

## Examples

``` r
mask <- array(TRUE, c(4L, 4L, 2L))
mask[1:2, , 1] <- FALSE
domain <- volume_domain(mask, spacing = c(3, 3, 3), id = "masked")
frame <- compile_frame(searchlights(3.1), domain)

# Identifiers are full-volume voxel indices; positions are rows 1..n.
head(frame$index$measurement)
#> [1]  3  4  7  8 11 12

# The same measurement, named both ways.
measurement_id(11)
#> <measurement_id> 11
measurement_position(match(11, frame$index$measurement))
#> <measurement_position> 5
```
