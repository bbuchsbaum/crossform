# Define a semantic condition space

A condition space identifies the named mean-model coordinates against
which an
[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md)
is declared. It is independent of design-matrix coding and column order.

## Usage

``` r
condition_space(
  coordinates,
  basis_id = "condition-means",
  units = "arbitrary",
  scale = 1,
  provenance = list()
)
```

## Arguments

- coordinates:

  Unique ordered semantic coordinates, such as condition names or
  condition-by-basis names.

- basis_id:

  One identifier for the targeted mean-model basis.

- units:

  One unit shared by all coordinates, or one per coordinate.

- scale:

  One positive scale shared by all coordinates, or one per coordinate.

- provenance:

  Portable semantic provenance.

## Value

An `effect_condition_space`: a list with `$coordinates`, `$basis_id`,
`$units` and `$scale` (named by coordinate), `$provenance`, and a
`$signature`. Two objects bind the same condition space only when their
signatures agree.

## See also

[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md)
to declare effects in this vocabulary and
[`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md)
to bind it to a compiled design.

Other studies and effect maps:
[`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md),
[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md),
[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md),
[`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md),
[`lower_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/lower_effect_map.md),
[`observation_model()`](https://bbuchsbaum.github.io/crossform/reference/observation_model.md),
[`raw_design_model()`](https://bbuchsbaum.github.io/crossform/reference/raw_design_model.md),
[`raw_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/raw_effect_map.md),
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md),
[`study_axis()`](https://bbuchsbaum.github.io/crossform/reference/study_axis.md),
[`study_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/study_capabilities.md)

## Examples

``` r
# The scientific vocabulary: three condition means in BOLD units.
conditions <- condition_space(
  c("face", "body", "tool"),
  basis_id = "scan-level-condition-mean:v1", units = "arbitrary-BOLD"
)
conditions
#> condition_space<3; scan-level-condition-mean:v1>
conditions$coordinates
#> [1] "face" "body" "tool"

# These coordinates are semantic, not design columns: the same space can be
# reached by cell-means or treatment coding, and the shared signature is
# what lets an effect map and a design model be checked against each other.
substr(conditions$signature, 1, 24)
#> [1] "condition-space-sha256:1"
```
