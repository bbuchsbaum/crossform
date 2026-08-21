# Inspect first-moment compiler conformance

The conformance court validates the portable receipt required by
`first-moment-relation-v1`. It reports construction evidence rather than
trusting an adapter name or package version.

## Usage

``` r
compiler_conformance(x)
```

## Arguments

- x:

  A
  [`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
  result.

## Value

A data frame with one row per partition: `partition` plus the logical
columns `semantic_identity`, `regressor_axis`, `condition_lowering`,
`row_lineage`, `rank_and_aliases`, `censor_accounting`,
`solver_diagnostics`, `whitening_provenance`, `source_revision`, and
`portable_receipt`.

## See also

[`relation_plan_receipts()`](https://bbuchsbaum.github.io/crossform/reference/relation_plan_receipts.md)
for the underlying receipts,
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
for the plan, and
[`study_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/study_capabilities.md)
for the equivalent report on the study facts.

Other relation planning and fitting:
[`adapter_version_certificate()`](https://bbuchsbaum.github.io/crossform/reference/adapter_version_certificate.md),
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
[`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md),
[`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)

## Examples

``` r
set.seed(1)
domain <- abstract_domain(3L, id = "conformance-example")
index <- observation_index(paste0("scan-", 1:4), "run-1")
facts <- study(observations(
  list(`run-1` = matrix(rnorm(12), 4L, 3L)), list(`run-1` = index), domain
))
design <- cbind(face = c(1, 0, 1, 0), body = c(0, 1, 0, 1))
rownames(design) <- paste0("scan-", 1:4)
target <- rbind(`face-body` = c(1, -1))
colnames(target) <- colnames(design)
plan <- plan_relation(
  facts, raw_design_model(list(`run-1` = design)), raw_effect_map(target),
  observation_model("ols", sampling_unit = "scan")
)

# Every field is evidence read back off the receipt, so this holds for an
# external compiler exactly as it does for a hand-built design.
conformance <- compiler_conformance(plan)
t(conformance)
#>                      [,1]   
#> partition            "run-1"
#> semantic_identity    "TRUE" 
#> regressor_axis       "TRUE" 
#> condition_lowering   "TRUE" 
#> row_lineage          "TRUE" 
#> rank_and_aliases     "TRUE" 
#> censor_accounting    "TRUE" 
#> solver_diagnostics   "TRUE" 
#> whitening_provenance "TRUE" 
#> source_revision      "TRUE" 
#> portable_receipt     "TRUE" 
all(as.matrix(conformance[-1L]))
#> [1] TRUE
```
