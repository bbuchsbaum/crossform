# Map a compact result vector back to a neuroim2 volume

The compact values are inserted at the exact full-volume indices carried
by a crossform volume domain. Features outside the domain receive
`fill`. This is an output adapter only; it performs no interpolation,
smoothing, or coordinate reinterpretation.

## Usage

``` r
as_neurovol(values, ...)

# S3 method for class 'numeric'
as_neurovol(
  values,
  mask,
  domain = NULL,
  fill = NA_real_,
  label = "crossform result",
  ...
)

# Default S3 method
as_neurovol(
  values,
  mask,
  domain = NULL,
  fill = NA_real_,
  label = "crossform result",
  ...
)
```

## Arguments

- values:

  One finite numeric value per compact domain *feature* (voxel), in
  `domain$feature_ids` order. crossform result views carry one value per
  *measurement* instead, which coincides with the features only for a
  voxelwise or searchlight frame. For a coarser frame, expand first with
  the frame's membership pattern —
  `as.numeric(Matrix::crossprod(frame$weights != 0, values))` — as shown
  in the "Measurements are not features" section of
  [`vignette("neuroim2-data")`](https://bbuchsbaum.github.io/crossform/articles/neuroim2-data.md).

- ...:

  Arguments passed on to methods. The methods crossform ships take no
  further arguments and refuse any.

- mask:

  The three-dimensional neuroim2 `NeuroVol` whose geometry defined
  `domain`.

- domain:

  The exact domain from
  [`neuroim2_volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_volume_domain.md).

- fill:

  Finite value written outside the compact domain.

- label:

  Optional result-volume label.

## Value

A neuroim2 `NeuroVol` on the same space as `mask`, carrying `values` at
`domain$feature_ids` and `fill` everywhere else.

## Details

`as_neurovol()` is an S3 generic dispatching on `values`, so a package
that owns its own result type can write that type out without crossform
having to know about it. Register a method the ordinary way —
`S3method(as_neurovol, my_result)` in your NAMESPACE — and it receives
`mask`, `domain`, `fill`, and `label` unchanged; it is expected to
return a `NeuroVol` on the space of `mask`. The generic validates
nothing itself, so a method is free to require different arguments, or
none beyond the object.

crossform ships the numeric-vector method described here. The default
method behaves identically, so any numeric vector still writes out
whether or not it carries a class, and anything that is not a numeric
vector is refused.

## See also

[`neuroim2_volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_volume_domain.md)
for the domain whose `feature_ids` fix the output positions, and
[`geometry_component()`](https://bbuchsbaum.github.io/crossform/reference/geometry_component.md)
for one source of the compact vector.

Other neural domains and frames:
[`abstract_domain()`](https://bbuchsbaum.github.io/crossform/reference/abstract_domain.md),
[`additive_frame()`](https://bbuchsbaum.github.io/crossform/reference/additive_frame.md),
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
if (requireNamespace("neuroim2", quietly = TRUE) &&
    utils::packageVersion("neuroim2") >= "0.19.0") {
  values <- array(FALSE, c(5L, 5L, 4L))
  values[2:4, 2:4, 2:3] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(
    values, neuroim2::NeuroSpace(c(5L, 5L, 4L), spacing = c(3, 3, 3))
  )
  domain <- neuroim2_volume_domain(mask)

  # One number per compact feature, in domain feature order.
  statistic <- seq_len(domain$n_features) / domain$n_features
  volume <- as_neurovol(statistic, mask, domain, label = "example statistic")

  # Values land at exactly the mask indices; everything else stays `fill`.
  print(dim(volume))
  print(identical(as.numeric(volume[domain$feature_ids]), statistic))
  print(all(is.na(as.array(volume)[!values])))

  # The generic is the extension point: a package with its own result type
  # registers a method for it and delegates the writing back here.
  as_neurovol.crossform_example_map <- function(values, mask, ...) {
    as_neurovol(values$statistic, mask, ...)
  }
  boxed <- structure(list(statistic = statistic),
    class = "crossform_example_map")
  print(identical(
    as.numeric(as_neurovol(boxed, mask, domain)[domain$feature_ids]),
    statistic
  ))
}
#> [1] 5 5 4
#> [1] TRUE
#> [1] TRUE
#> [1] TRUE
```
