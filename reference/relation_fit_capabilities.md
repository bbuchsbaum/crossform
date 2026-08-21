# Inspect statistical capabilities of a relation or relation fit

Call this before requesting anything that needs residuals, so a missing
error channel surfaces as an explicit `FALSE` rather than a later
refusal.

## Usage

``` r
relation_fit_capabilities(x)
```

## Arguments

- x:

  An `effect_relation` or `effect_relation_fit`.

## Value

A data frame with one row per partition: `partition` plus the logical
columns `error_model`, `residual_blocks`, `effect_covariance`,
`residual_df`, `separable_error`, `learned_metric_input`, and
`within_participant_calibration`.

## See also

[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md)
to obtain the capabilities,
[`sampling_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/sampling_capabilities.md)
for the uncertainty-side report, and
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)
to inspect the refusal raised when one is missing.

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
[`relation_plan_receipts()`](https://bbuchsbaum.github.io/crossform/reference/relation_plan_receipts.md),
[`residual_block()`](https://bbuchsbaum.github.io/crossform/reference/residual_block.md),
[`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md),
[`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)

## Examples

``` r
example <- example_fmri_effects()

# A fit built from raw responses carries the full error channel.
relation_fit_capabilities(example$fit)[
  , c("partition", "residual_blocks", "within_participant_calibration")
]
#>   partition residual_blocks within_participant_calibration
#> 1      run1            TRUE                           TRUE
#> 2      run2            TRUE                           TRUE
#> 3      run3            TRUE                           TRUE
#> 4      run4            TRUE                           TRUE

# The bare relation underneath reports every statistical capability FALSE:
# the point geometry is intact, the uncertainty channel is not there.
relation_fit_capabilities(example$fit$relation)[1L, ]
#>   partition error_model residual_blocks effect_covariance residual_df
#> 1      run1       FALSE           FALSE             FALSE       FALSE
#>   separable_error learned_metric_input within_participant_calibration
#> 1           FALSE                FALSE                          FALSE
```
