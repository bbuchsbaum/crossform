# Define an experimental coordinate space

An effect space gives experimental coordinates a semantic identity
beyond their display labels. Its stable signature covers coordinate
order, basis, units, scaling, and provenance.

## Usage

``` r
effect_space(
  coordinates,
  basis_id = "unspecified",
  units = "arbitrary",
  scale = 1,
  provenance = list()
)
```

## Arguments

- coordinates:

  Unique ordered coordinate identifiers.

- basis_id:

  One nonempty semantic basis identifier.

- units:

  One unit shared by all coordinates, or one per coordinate.

- scale:

  One positive finite scale shared by all coordinates, or one per
  coordinate.

- provenance:

  Optional design or contrast provenance as a list.

## Value

An `effect_space`: a list with `$coordinates` (the ordered names),
`$basis_id`, `$units` and `$scale` (both named by coordinate),
`$provenance`, and a `$signature` that changes whenever any of those
does. Treat the value as immutable.

## See also

[`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md)
for the mean-model vocabulary an
[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md)
is declared against, and
[`effect_extractor()`](https://bbuchsbaum.github.io/crossform/reference/effect_extractor.md)
for the map that produces coordinates in this space.

Other studies and effect maps:
[`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md),
[`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md),
[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md),
[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md),
[`lower_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/lower_effect_map.md),
[`observation_model()`](https://bbuchsbaum.github.io/crossform/reference/observation_model.md),
[`raw_design_model()`](https://bbuchsbaum.github.io/crossform/reference/raw_design_model.md),
[`raw_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/raw_effect_map.md),
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md),
[`study_axis()`](https://bbuchsbaum.github.io/crossform/reference/study_axis.md),
[`study_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/study_capabilities.md)

## Examples

``` r
# Name the coordinates a fit will produce, and the basis they live in.
space <- effect_space(
  c("face", "house"),
  basis_id = "condition-means:v1", units = "arbitrary-BOLD"
)
space$coordinates
#> [1] "face"  "house"
space$units[["face"]]
#> [1] "arbitrary-BOLD"

# The signature covers meaning, not display labels: the same two names in
# the default unspecified basis are a different space, so objects built
# against one cannot be silently combined with the other.
identical(space$signature, effect_space(c("face", "house"))$signature)
#> [1] FALSE
```
