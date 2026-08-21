# Bind already compiled raw design matrices

This degenerate route binds design values, coefficient order, row
identity, and solver into model identity. It makes no semantic
coding-invariance claim.

## Usage

``` r
raw_design_model(designs, row_ids = NULL, solver = "auto", provenance = list())
```

## Arguments

- designs:

  Named observation-by-coefficient design matrices.

- row_ids:

  Ordered observation identifiers per design.

- solver:

  Numerical route per partition.

- provenance:

  Portable provenance for the external construction.

## Value

An `effect_raw_design_model` (also an `effect_design_model`) with the
same fields as
[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md),
except that `$condition_space` and `$parameterizations` are `NULL`,
`$design_model_id` includes the design values themselves, and
`$capabilities` reports `symbolic_model` and `coding_invariant` as
`FALSE`.

## See also

[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md)
for the semantic route, and
[`raw_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/raw_effect_map.md),
the effect map this design must be paired with in
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md).

Other studies and effect maps:
[`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md),
[`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md),
[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md),
[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md),
[`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md),
[`lower_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/lower_effect_map.md),
[`observation_model()`](https://bbuchsbaum.github.io/crossform/reference/observation_model.md),
[`raw_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/raw_effect_map.md),
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md),
[`study_axis()`](https://bbuchsbaum.github.io/crossform/reference/study_axis.md),
[`study_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/study_capabilities.md)

## Examples

``` r
design <- cbind(b1 = c(1, 0, 1, 0), b2 = c(0, 1, 0, 1))
rownames(design) <- paste0("scan-", 1:4)
model <- raw_design_model(list(`run-1` = design))

# The honest cost of skipping the condition space: the numeric design is
# part of the identity, so a re-coded design is a different request.
model$capabilities$coding_invariant
#> [1] FALSE
model$capabilities$row_lineage
#> [1] TRUE

# A raw design must be paired with a raw target on the same column axis.
target <- rbind(`b1-b2` = c(1, -1))
colnames(target) <- colnames(design)
raw_effect_map(target)
#> raw_effect_map<1 effects x 2 coefficients; raw-X-T>
```
