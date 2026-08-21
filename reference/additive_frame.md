# Describe an additive diagonal spatial frame

An additive frame is the only frame representation whose locations
collapse to rows of one spatial contraction in the version 0.1 compiler.

## Usage

``` r
additive_frame(
  weights,
  normalization = "none",
  domain_id = "abstract",
  domain = NULL,
  members = NULL,
  measurements = NULL,
  construction = list(),
  specification = NULL
)
```

## Arguments

- weights:

  A finite, nonnegative measurement-by-feature base or sparse `Matrix`
  matrix. Supply this or `members`, never both.

- normalization:

  One of `none`, `local` (row sums equal one), or `conservative` (column
  sums equal one).

- domain_id:

  Stable identity of the neural feature domain.

- domain:

  Optional exact `effect_domain` or internal domain reference. The
  `members` route requires an exact domain.

- members:

  A nonempty list of integer vectors, one per measurement, giving the
  domain feature *positions* that measurement is supported on. This is
  the neighborhood route described under *Two routes*.

- measurements:

  Identifiers for the rows of the `members` route, one per measurement,
  recorded as `$index$measurement`. Defaults to the measurement
  positions.

- construction:

  Named list recording how a `members` neighborhood was produced — the
  rule, the provider, its parameters — kept with the support pattern and
  folded into its identity.

- specification:

  Optional named list recording what generated a `members` frame, kept
  as `$specification` the way a
  [`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md)
  result keeps the scope it was compiled from.

## Value

An `effect_frame` with `representation = "additive_diagonal"`, carrying
the `$weights` matrix, its `$normalization`, and the `$domain` reference
the weights are bound to.

## Structure

A declared frame carries the operator and the domain it claims, and
nothing about how it was chosen.

- `$weights`: the measurement-by-feature operator exactly as supplied.
  Row `m` holds the weight each domain feature contributes to
  measurement `m`, in domain feature order.

- `$normalization`: the normalization asserted about those rows. It is
  checked, not applied.

- `$domain`: the neural domain reference the columns are bound to, and
  `$domain_id` its identity.

Unlike a
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md)
result, a frame declared from `weights` has no `$index` and no
`$specification`: nothing generated it, so views index its measurements
by position. Any other element is internal and may change.

## Two routes

A frame can be declared from an operator or built from neighborhoods,
and the difference is who applies the normalization.

- **`weights`** is the declaration: the operator is taken exactly as
  supplied and `normalization` is *checked* against it, so unnormalized
  rows are refused rather than silently rescaled.

- **`members`** is the neighborhood route, for a spatial provider that
  computes its own supports — searchlights from an external package,
  parcels from an atlas, any rule crossform does not implement. Each
  element of `members` lists the domain feature positions one
  measurement covers; `normalization` is then *applied* to that
  membership pattern, exactly as
  [`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md)
  applies it to the neighborhoods it computes itself. What comes back is
  a frame of the same shape a compiled one has, carrying `$index`, the
  `$specification` you recorded, and the support pattern that makes a
  locality-aware plan possible.

The neighborhood route is the reason an external searchlight provider
does not have to reproduce crossform's normalization law or its support
bookkeeping to hand back a frame the rest of the package accepts;
[`neuroim2_searchlights()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_searchlights.md)
is written this way and is the worked example.

## See also

[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md)
with
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md),
[`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md),
or
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md),
which build additive frames from a neural domain;
[`neuroim2_searchlights()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_searchlights.md)
for the `members` route in use; and
[`measurement_frame()`](https://bbuchsbaum.github.io/crossform/reference/measurement_frame.md),
which adapts one into oriented measurements.

Other neural domains and frames:
[`abstract_domain()`](https://bbuchsbaum.github.io/crossform/reference/abstract_domain.md),
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
# Two overlapping local averages over four features, declared directly.
# `"local"` asserts that each measurement's weights already sum to one, so
# a frame row is a weighted mean rather than a weighted sum.
domain <- abstract_domain(4, id = "demo:native:v1")
frame <- additive_frame(
  matrix(c(
    1 / 2, 1 / 2, 0, 0,
    0, 1 / 3, 1 / 3, 1 / 3
  ), 2, 4, byrow = TRUE),
  normalization = "local", domain = domain
)
dim(frame$weights)
#> [1] 2 4
rowSums(as.matrix(frame$weights))
#> [1] 1 1

# The assertion is checked, not applied: unnormalized rows are refused
# rather than silently rescaled.
unnormalized <- try(
  additive_frame(matrix(1, 1, 4), normalization = "local", domain = domain),
  silent = TRUE
)
conditionMessage(attr(unnormalized, "condition"))
#> [1] "Locally normalized frame rows must sum to one."

# The declared width must match the domain it claims.
wrong <- try(additive_frame(matrix(1, 1, 3), domain = domain), silent = TRUE)
conditionMessage(attr(wrong, "condition"))
#> [1] "The frame width must match its exact neural domain."

# The neighborhood route: hand over the supports a provider computed and
# the declared normalization is applied to them, not asserted about them.
neighborhoods <- additive_frame(
  members = list(1:2, 2:4), measurements = c("left", "right"),
  normalization = "local", domain = domain,
  construction = list(kind = "declared_neighborhoods", provider = "example")
)
rowSums(as.matrix(neighborhoods$weights))
#> [1] 1 1

# It comes back shaped like a compiled frame: the measurements are named,
# not merely positional.
neighborhoods$index$measurement
#> [1] "left"  "right"
```
