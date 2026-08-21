# Bind a raw coefficient-space target

This is the honest degenerate route for a supplied `T` matrix. Its
identity includes the coefficient order and target values, and it cannot
claim symbolic or coding-invariant effects.

## Usage

``` r
raw_effect_map(
  target,
  effects = rownames(target),
  coefficients = colnames(target),
  units = "arbitrary",
  scale = 1,
  provenance = list()
)
```

## Arguments

- target:

  A finite effect-by-coefficient target matrix.

- effects:

  An
  [`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md)
  or output effect names.

- coefficients:

  Unique ordered coefficient names.

- units, scale:

  Output units and scale when `effects` is a character vector.

- provenance:

  Portable target provenance.

## Value

An `effect_raw_map` (also an `effect_lowered_map`): a list with the
effect-by-coefficient `$target`, its `$effect_space`, a `$coding_id` of
`"raw-X-T"`, `$capabilities` with `symbolic_effects`,
`valid_effect_lowering`, and `coding_invariant` all `FALSE`, and a
`$lowering_id` that includes the target values themselves.

## See also

[`lower_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/lower_effect_map.md)
for the semantic route, and
[`raw_design_model()`](https://bbuchsbaum.github.io/crossform/reference/raw_design_model.md),
the design model this map must be paired with.

Other studies and effect maps:
[`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md),
[`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md),
[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md),
[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md),
[`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md),
[`lower_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/lower_effect_map.md),
[`observation_model()`](https://bbuchsbaum.github.io/crossform/reference/observation_model.md),
[`raw_design_model()`](https://bbuchsbaum.github.io/crossform/reference/raw_design_model.md),
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md),
[`study_axis()`](https://bbuchsbaum.github.io/crossform/reference/study_axis.md),
[`study_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/study_capabilities.md)

## Examples

``` r
# Use this when you have a target matrix but no condition-space meaning for
# the coefficient columns it multiplies.
target <- rbind(`face-body` = c(1, -1, 0))
colnames(target) <- c("beta_face", "beta_body", "beta_drift")
raw <- raw_effect_map(target, units = "arbitrary-BOLD")
raw
#> raw_effect_map<1 effects x 3 coefficients; raw-X-T>

# The honest cost: the numeric values are part of the identity, and no
# coding-invariance claim is made, unlike lower_effect_map().
unlist(raw$capabilities)
#>      symbolic_effects valid_effect_lowering      coding_invariant 
#>                 FALSE                 FALSE                 FALSE 
```
