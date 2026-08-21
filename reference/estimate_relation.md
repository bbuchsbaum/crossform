# Estimate a planned relation family

This is the authorized first-moment execution verb. It emits the
existing `effect_relation_fit` contract, so every second-moment plan and
view remains unchanged.

## Usage

``` r
estimate_relation(x, ...)
```

## Arguments

- x:

  An
  [`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
  result.

- ...:

  Reserved for future execution policies.

## Value

An `effect_relation_fit` whose `$signature` binds the relation plan,
source revisions, realized row lineage, and design receipts, and whose
`$provenance` records `relation_plan_id`, `design_receipt_ids`,
`study_id`, and `observation_model_id`.

## See also

[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
for the plan it executes,
[`fmrireg_relation()`](https://bbuchsbaum.github.io/crossform/reference/fmrireg_relation.md)
for the external point-parity adapter, and
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
for the second-moment question that follows.

Other relation planning and fitting:
[`adapter_version_certificate()`](https://bbuchsbaum.github.io/crossform/reference/adapter_version_certificate.md),
[`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md),
[`effect_extractor()`](https://bbuchsbaum.github.io/crossform/reference/effect_extractor.md),
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
domain <- abstract_domain(3L, id = "estimate-relation-example")
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

# Reading neural values happens here, and only here.
fit <- estimate_relation(plan)
round(relation_block(fit, "run-1", 1:3), 3)
#>             [,1] [,2]  [,3]
#> face-body -1.621 0.45 1.002

# A fixed observation model earns the residual channel, and the fit records
# which plan produced it.
relation_fit_capabilities(fit)$residual_blocks
#> [1] TRUE
identical(fit$provenance$relation_plan_id, plan$relation_plan_id)
#> [1] TRUE
```
