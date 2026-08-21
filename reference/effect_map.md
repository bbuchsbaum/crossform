# Declare effects in a semantic condition vocabulary

The rows of `weights` are named scientific effects and the columns are
bound semantic condition coordinates. The returned identity never
contains a compiled design matrix or coefficient coding.

## Usage

``` r
effect_map(
  weights,
  conditions = colnames(weights),
  effects = rownames(weights),
  units = NULL,
  scale = 1,
  component = "canonical-amplitude",
  provenance = list()
)
```

## Arguments

- weights:

  A finite effect-by-condition matrix.

- conditions:

  A
  [`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md)
  or condition names. Defaults to the column names of `weights`.

- effects:

  An
  [`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md)
  or output effect names. Defaults to the row names of `weights`.

- units:

  Output units. By default the shared condition-space unit.

- scale:

  Output scale.

- component:

  One targeted temporal-basis component.

- provenance:

  Portable provenance for the declared functional.

## Value

An `effect_condition_map`: a list with the effect-by-condition
`$weights`, the bound `$condition_space` and derived `$effect_space`,
the targeted `$component`, and an `$effect_map_id` covering all of them.

## See also

[`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md)
for the input vocabulary,
[`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md)
plus
[`lower_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/lower_effect_map.md)
to reach a coefficient axis, and
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
to request these effects from a study.

Other studies and effect maps:
[`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md),
[`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md),
[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md),
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

# Request the three condition means and one named contrast at once, so
# later views can form the contrast without refitting the run models.
weights <- rbind(
  face = c(1, 0, 0), body = c(0, 1, 0), tool = c(0, 0, 1),
  `face-body` = c(1, -1, 0)
)
colnames(weights) <- conditions$coordinates
effects <- effect_map(weights, conditions)
effects
#> effect_map<4 effects x 3 conditions; canonical-amplitude>
effects$weights["face-body", ]
#> face body tool 
#>    1   -1    0 

# The derived effect space inherits the condition units, and its identity
# records the functional rather than any design coding.
effects$effect_space$basis_id
#> [1] "condition-functional:canonical-amplitude"
```
