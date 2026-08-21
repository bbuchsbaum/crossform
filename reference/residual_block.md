# Read fitted residuals for a neural feature block

Residuals are produced lazily, one feature block at a time, without ever
materializing a dense observation residualizer. This is the input a
learned neural metric or an analytic RDM covariance draws on.

## Usage

``` r
residual_block(x, partition, features)
```

## Arguments

- x:

  An `effect_relation_fit` with residual-block capability.

- partition:

  One partition name or index.

- features:

  Unique neural feature indices.

## Value

A whitened residual matrix with one row per observation in the partition
and one column per requested feature.

## See also

[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md)
which installs the residual channel,
[`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md)
for the matching degrees of freedom, and
[`noise_precision()`](https://bbuchsbaum.github.io/crossform/reference/noise_precision.md),
which learns a metric from these blocks.

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
[`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md),
[`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)

## Examples

``` r
example <- example_fmri_effects()

# Whitened residuals for the first two neural features of one run.
residuals <- residual_block(example$fit, "run1", 1:2)
dim(residuals)
#> [1] 32  2

# They are orthogonal to the fitted design, so the condition means have
# already been projected out.
round(colMeans(residuals), 8)
#> 1 2 
#> 0 0 

# Divide by residual_df(), not nrow(), to estimate the noise variance.
round(colSums(residuals^2) / residual_df(example$fit, "run1"), 3)
#>     1     2 
#> 0.507 0.317 
```
