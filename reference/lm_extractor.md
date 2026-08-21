# Compile a supplied linear model into an effect extractor

The function accepts an already constructed design matrix; it does not
own event timing, HRF construction, nuisance selection, or a formula
language. Full-rank designs use pivoted QR. Rank-deficient designs use
an SVD pseudoinverse only after every requested effect is proven
estimable.

## Usage

``` r
lm_extractor(
  design,
  effects,
  observation_whitener = NULL,
  effect_names = rownames(effects),
  tolerance = sqrt(.Machine$double.eps),
  whiten = NULL,
  solver = c("auto", "qr", "svd")
)
```

## Arguments

- design:

  Finite observation-by-coefficient design matrix.

- effects:

  Finite effect-by-coefficient target matrix.

- observation_whitener:

  Optional finite square observation whitener `L`.

- effect_names:

  Optional names or an
  [`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md)
  for target effects.

- tolerance:

  Positive rank and estimability tolerance.

- whiten:

  Deprecated alias for `observation_whitener`. Supplying both is an
  error.

- solver:

  Numerical factorization route: automatic, pivoted QR, or SVD. This
  changes execution provenance, not the requested effect.

## Value

An `effect_extractor` whose `$map` is `T (L X)^+ L`, with `$diagnostics`
recording `solver`, `solver_policy`, `observations`, `coefficients`,
`rank`, `rank_deficient`, per-effect `estimability_error`, `tolerance`,
and the `observation_whitener` descriptor.

## See also

[`effect_extractor()`](https://bbuchsbaum.github.io/crossform/reference/effect_extractor.md)
for a hand-specified map, and
[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md),
which additionally keeps the residual error channel.

Other relation planning and fitting:
[`adapter_version_certificate()`](https://bbuchsbaum.github.io/crossform/reference/adapter_version_certificate.md),
[`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md),
[`effect_extractor()`](https://bbuchsbaum.github.io/crossform/reference/effect_extractor.md),
[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md),
[`file_matrix_source()`](https://bbuchsbaum.github.io/crossform/reference/file_matrix_source.md),
[`fmridesign_design_model()`](https://bbuchsbaum.github.io/crossform/reference/fmridesign_design_model.md),
[`fmrireg_relation()`](https://bbuchsbaum.github.io/crossform/reference/fmrireg_relation.md),
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
condition <- factor(rep(c("face", "body"), each = 3L))
design <- cbind(
  stats::model.matrix(~ 0 + condition), drift = seq(-1, 1, length.out = 6L)
)
colnames(design)[1:2] <- c("face", "body")

# Ask for one contrast on the design's coefficient axis.
target <- rbind(`face-body` = c(1, -1, 0))
colnames(target) <- colnames(design)
extractor <- lm_extractor(design, target)
extractor$diagnostics$solver
#> [1] "pivoted_qr"
round(extractor$map, 3)
#>            [,1]   [,2]   [,3]  [,4]  [,5]   [,6]
#> face-body 0.417 -0.333 -1.083 1.083 0.333 -0.417

# Adding an intercept aliases the two condition columns, so a request for
# `face` alone is refused and the aliased regressors are named.
aliased <- cbind(intercept = 1, design)
request <- matrix(0, 1L, ncol(aliased),
  dimnames = list("face", colnames(aliased)))
request[, "face"] <- 1
refusal <- catch_refusal(lm_extractor(aliased, request))
refusal$capability
#> [1] "estimable_effects"
refusal$reasons
#> [1] "The design has rank 3 for 4 regressors."                
#> [2] "Aliased regressor set: intercept, face, body."          
#> [3] "Requested effect outside the estimable row space: face."
```
