# Inspect portable design receipts

Returns the per-partition record of how each design was actually
compiled: the matrix used, the rows censoring retained, the achieved
rank and any aliased regressors, the solver, and the whitening
provenance.

## Usage

``` r
relation_plan_receipts(x)
```

## Arguments

- x:

  An
  [`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
  result.

## Value

A named list of `effect_design_receipt` values, one per partition. Each
carries `$design`, `$coefficient_axis`, `$lowered_target`,
`$effect_space`, `$row_lineage`, `$censoring`, `$solver`, `$rank`,
`$aliases`, `$residual_df`, `$observation_whitener`, `$capabilities`,
and a `$design_receipt_id`.

## See also

[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
for the plan,
[`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md)
for the boolean conformance summary of the same receipts, and
[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md)
to execute them.

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
[`residual_block()`](https://bbuchsbaum.github.io/crossform/reference/residual_block.md),
[`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md),
[`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)

## Examples

``` r
set.seed(1)
domain <- abstract_domain(3L, id = "receipts-example")
index <- observation_index(paste0("scan-", 1:4), "run-1")
facts <- study(observations(
  list(`run-1` = matrix(rnorm(12), 4L, 3L)), list(`run-1` = index), domain
))
design <- cbind(face = c(1, 0, 1, 0), body = c(0, 1, 0, 1))
rownames(design) <- paste0("scan-", 1:4)
target <- rbind(`face-body` = c(1, -1))
colnames(target) <- colnames(design)

# Even the raw route, which claims no semantic coding, yields a complete
# receipt: rank, aliases, censoring, solver, and whitening are all recorded.
plan <- plan_relation(
  facts, raw_design_model(list(`run-1` = design)), raw_effect_map(target),
  observation_model("ols", sampling_unit = "scan")
)
receipts <- relation_plan_receipts(plan)
names(receipts)
#> [1] "run-1"
receipts$`run-1`$rank
#> [1] 2
receipts$`run-1`$residual_df
#> [1] 2
receipts$`run-1`$aliases
#> character(0)
```
