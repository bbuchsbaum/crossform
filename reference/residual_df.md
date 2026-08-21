# Read residual degrees of freedom from a fitted partition

The divisor for any noise-variance estimate built from
[`residual_block()`](https://bbuchsbaum.github.io/crossform/reference/residual_block.md):
observations minus the numerical rank of the whitened design, not minus
the number of design columns.

## Usage

``` r
residual_df(x, partition)
```

## Arguments

- x:

  An `effect_relation_fit` with residual-block capability.

- partition:

  One partition name or index.

## Value

One positive integer.

## See also

[`residual_block()`](https://bbuchsbaum.github.io/crossform/reference/residual_block.md)
and
[`effect_covariance()`](https://bbuchsbaum.github.io/crossform/reference/effect_covariance.md)
for the rest of the error channel, and
[`relation_fit_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit_capabilities.md)
to check availability.

Other relation planning and fitting:
[`adapter_version_certificate()`](https://bbuchsbaum.github.io/crossform/reference/adapter_version_certificate.md),
[`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md),
[`effect_extractor()`](https://bbuchsbaum.github.io/crossform/reference/effect_extractor.md),
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
[`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)

## Examples

``` r
example <- example_fmri_effects()

# 32 trials per run minus the four estimated condition means.
residual_df(example$fit, "run1")
#> [1] 28

# Each independent run contributes its own degrees of freedom.
vapply(example$fit$relation$partitions,
  function(partition) residual_df(example$fit, partition), integer(1))
#> run1 run2 run3 run4 
#>   28   28   28   28 

# A relation built from precomputed betas has none to report.
catch_refusal(residual_df(example$fit$relation, "run1"))$capability
#> [1] "residual_df"
```
