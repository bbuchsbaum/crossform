# Lower semantic effects into a concrete coefficient basis

Lower semantic effects into a concrete coefficient basis

## Usage

``` r
lower_effect_map(effects, parameterization)
```

## Arguments

- effects:

  An
  [`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md).

- parameterization:

  A
  [`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md)
  for the same semantic condition space.

## Value

An `effect_lowered_map`: a list with the effect-by-coefficient
`$target`, the `$effect_space` and `$condition_space` it came from, the
`$effect_map_id`, `$parameterization_id`, and `$coding_id` receipt
fields, `$capabilities` (all `TRUE` on this route), and a
`$lowering_id`.

## See also

[`raw_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/raw_effect_map.md)
for the degenerate route with no condition space, and
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
which lowers effects for every partition.

Other studies and effect maps:
[`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md),
[`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md),
[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md),
[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md),
[`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md),
[`observation_model()`](https://bbuchsbaum.github.io/crossform/reference/observation_model.md),
[`raw_design_model()`](https://bbuchsbaum.github.io/crossform/reference/raw_design_model.md),
[`raw_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/raw_effect_map.md),
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md),
[`study_axis()`](https://bbuchsbaum.github.io/crossform/reference/study_axis.md),
[`study_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/study_capabilities.md)

## Examples

``` r
conditions <- condition_space(c("face", "body"))
weights <- rbind(`face-body` = c(1, -1))
colnames(weights) <- conditions$coordinates
effects <- effect_map(weights, conditions)

# A design with a drift nuisance column, and the coding that names which
# coefficients carry the condition means.
map <- cbind(diag(2), drift = 0)
dimnames(map) <- list(
  conditions$coordinates, c(conditions$coordinates, "drift")
)
coding <- coefficient_parameterization(
  map, conditions, coding_id = "cell-means-plus-drift"
)

# Lowering moves the request from condition names onto the design's
# coefficient axis without changing what was asked for.
lowered <- lower_effect_map(effects, coding)
lowered$target
#>           face body drift
#> face-body    1   -1     0
```
