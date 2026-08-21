# Construct an explicit linear effect extractor

An extractor is the declared map `E` in `B = E Y`. It contains no neural
data and may therefore be reused across feature blocks.

## Usage

``` r
effect_extractor(
  map,
  effects = rownames(map),
  estimator = "explicit",
  diagnostics = list()
)
```

## Arguments

- map:

  A finite effect-by-observation numeric matrix.

- effects:

  An
  [`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md)
  or unique names used as shorthand for an unspecified-basis effect
  space.

- estimator:

  Short estimator identity.

- diagnostics:

  Optional estimator diagnostics.

## Value

An `effect_extractor`: a list with the effect-by-observation `$map`, its
`$effect_space` and `$effects` labels, `$n_observations`, the
`$estimator` identity, and `$diagnostics`.

## See also

[`lm_extractor()`](https://bbuchsbaum.github.io/crossform/reference/lm_extractor.md)
to compile one from a design and target, and
[`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md),
which applies one extractor per partition.

Other relation planning and fitting:
[`adapter_version_certificate()`](https://bbuchsbaum.github.io/crossform/reference/adapter_version_certificate.md),
[`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md),
[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md),
[`file_matrix_source()`](https://bbuchsbaum.github.io/crossform/reference/file_matrix_source.md),
[`fmridesign_design_model()`](https://bbuchsbaum.github.io/crossform/reference/fmridesign_design_model.md),
[`fmrireg_relation()`](https://bbuchsbaum.github.io/crossform/reference/fmrireg_relation.md),
[`lm_extractor()`](https://bbuchsbaum.github.io/crossform/reference/lm_extractor.md),
[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md),
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md),
[`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md),
[`relation_block()`](https://bbuchsbaum.github.io/crossform/reference/relation_block.md),
[`relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit.md),
[`relation_fit_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit_capabilities.md),
[`relation_plan_receipts()`](https://bbuchsbaum.github.io/crossform/reference/relation_plan_receipts.md),
[`residual_block()`](https://bbuchsbaum.github.io/crossform/reference/residual_block.md),
[`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md),
[`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)

## Examples

``` r
# Average four observations, and contrast the first two with the last two.
map <- rbind(
  mean = c(0.25, 0.25, 0.25, 0.25),
  difference = c(0.5, 0.5, -0.5, -0.5)
)
extractor <- effect_extractor(map, estimator = "hand-specified")
extractor$effects
#> [1] "mean"       "difference"
extractor$n_observations
#> [1] 4

# The extractor holds no neural data, so the same map is reused across every
# feature block a relation reads.
set.seed(1)
extractor$map %*% matrix(rnorm(8), 4L, 2L)
#>                   [,1]       [,2]
#> mean        0.07921043  0.1836983
#> difference -0.60123134 -0.8583572
```
