# Plan a first-moment experimental–neural relation

The plan composes factual study identity, a semantic or raw design
model, named effects, and the observation-model assumptions. Design
compilation and estimability are validated without reading neural
values.

## Usage

``` r
plan_relation(
  study,
  model,
  effects,
  observation_model,
  tolerance = sqrt(.Machine$double.eps)
)
```

## Arguments

- study:

  A
  [`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md).

- model:

  A
  [`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md)
  or
  [`raw_design_model()`](https://bbuchsbaum.github.io/crossform/reference/raw_design_model.md).

- effects:

  An
  [`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md)
  for a semantic model, or
  [`raw_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/raw_effect_map.md)
  for a raw design.

- observation_model:

  An
  [`observation_model()`](https://bbuchsbaum.github.io/crossform/reference/observation_model.md).

- tolerance:

  Positive numerical rank and estimability tolerance.

## Value

An `effect_relation_plan`: a list with the validated `$study`, `$model`,
`$effects` and `$observation_model`, the `$partitions`, the
per-partition `$lowered_effects` and portable `$design_receipts`, the
`$retained_rows` left by censoring, the resolved `$sampling_unit` and
`$whiteners`, the `$tolerance`, `$capabilities`, and a
`$relation_plan_id`.

## See also

[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md)
to execute the plan,
[`relation_plan_receipts()`](https://bbuchsbaum.github.io/crossform/reference/relation_plan_receipts.md)
and
[`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md)
to inspect it, and
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md),
[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md),
[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md),
[`observation_model()`](https://bbuchsbaum.github.io/crossform/reference/observation_model.md)
for the four inputs.

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
# Four scans, two conditions, three neural features.
set.seed(1)
domain <- abstract_domain(3L, id = "plan-relation-example")
index <- observation_index(paste0("scan-", 1:4), "run-1")
facts <- study(observations(
  list(`run-1` = matrix(rnorm(12), 4L, 3L)), list(`run-1` = index), domain
))

conditions <- condition_space(c("face", "body"))
design <- cbind(face = c(1, 0, 1, 0), body = c(0, 1, 0, 1))
rownames(design) <- paste0("scan-", 1:4)
coding <- coefficient_parameterization(
  diag(2), conditions,
  coefficients = colnames(design), coding_id = "cell-means"
)
model <- design_model(
  list(target = "condition means"), conditions,
  designs = list(`run-1` = design),
  parameterizations = list(`run-1` = coding)
)
weights <- rbind(`face-body` = c(1, -1))
colnames(weights) <- conditions$coordinates

# Compilation and estimability are checked without reading neural values.
plan <- plan_relation(
  facts, model, effect_map(weights, conditions),
  observation_model("ols", sampling_unit = "scan")
)
plan
#> relation_plan<1 partitions; semantic; uncertainty analytic>
plan$design_receipts$`run-1`$rank
#> [1] 2

# Effects and design must bind the same condition space: a different basis
# is a different scientific request, and is refused rather than coerced.
other <- condition_space(c("face", "body"), basis_id = "other-basis")
catch_refusal(plan_relation(
  facts, model, effect_map(weights, other),
  observation_model("ols", sampling_unit = "scan")
))$capability
#> [1] "valid_effect_lowering"
```
