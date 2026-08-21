# Identify one concrete coefficient parameterization

`map` declares `mu = map %*% beta`, where `mu` follows the semantic
condition coordinates and `beta` follows the compiled coefficient axis.
The object belongs to a design receipt, not scientific plan identity.

## Usage

``` r
coefficient_parameterization(
  map,
  conditions,
  coefficients = colnames(map),
  coding_id,
  provenance = list(),
  tolerance = sqrt(.Machine$double.eps)
)
```

## Arguments

- map:

  A finite condition-by-coefficient matrix.

- conditions:

  The bound
  [`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md).

- coefficients:

  Unique ordered coefficient names. Defaults to `map` column names.

- coding_id:

  One coding identifier recorded in the receipt.

- provenance:

  Portable compiler provenance.

- tolerance:

  Positive rank tolerance.

## Value

An `effect_coefficient_parameterization`: a list with the
condition-by-coefficient `$map`, the bound `$condition_space`, the
`$coefficients` axis, the `$coding_id`, the achieved `$semantic_rank`,
and a `$parameterization_id`.

## See also

[`lower_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/lower_effect_map.md)
to apply it to an
[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md),
and
[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md)
which carries one parameterization per partition.

Other studies and effect maps:
[`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md),
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
conditions <- condition_space(c("face", "body", "tool"))

# The compiled design adds a drift column. This declares that the three
# condition means are the first three coefficients and ignore drift, so the
# effect request survives a change of coding.
map <- cbind(diag(3), drift = 0)
dimnames(map) <- list(
  conditions$coordinates, c(conditions$coordinates, "drift")
)
coding <- coefficient_parameterization(
  map, conditions, coding_id = "cell-means-plus-drift"
)
coding
#> coefficient_parameterization<3 conditions x 4 coefficients; cell-means-plus-drift>
coding$semantic_rank
#> [1] 3

# A coding that cannot identify every condition is refused, not rounded.
collapsed <- map
collapsed["body", ] <- collapsed["face", ]
refusal <- catch_refusal(coefficient_parameterization(
  collapsed, conditions, coding_id = "face-and-body-collapsed"
))
refusal$capability
#> [1] "valid_effect_lowering"
```
