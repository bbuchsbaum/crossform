# Read one experimental-neural relation block

This is where neural values are actually read. Only the requested
feature columns are pulled from the source and mapped through that
partition's extractor, which is what keeps out-of-memory sources
bounded.

## Usage

``` r
relation_block(x, partition, features)
```

## Arguments

- x:

  An `effect_relation` or `effect_relation_fit`.

- partition:

  One partition name or index.

- features:

  Unique neural feature indices.

## Value

A finite effect-by-feature matrix with one row per effect coordinate, in
effect-space order, and one column per requested feature.

## See also

[`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md)
and
[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md)
for the objects it reads, and
[`residual_block()`](https://bbuchsbaum.github.io/crossform/reference/residual_block.md)
for the matching residual read.

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
[`relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit.md),
[`relation_fit_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit_capabilities.md),
[`relation_plan_receipts()`](https://bbuchsbaum.github.io/crossform/reference/relation_plan_receipts.md),
[`residual_block()`](https://bbuchsbaum.github.io/crossform/reference/residual_block.md),
[`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md),
[`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)

## Examples

``` r
set.seed(1)
domain <- abstract_domain(5L, id = "relation-block-example")
betas <- matrix(rnorm(15), 3L, 5L,
  dimnames = list(c("face", "body", "tool"), NULL))
point <- relation(list(`run-1` = betas), domain = domain)

# Ask for two features only; the source is read for those columns alone.
round(relation_block(point, "run-1", c(1L, 4L)), 3)
#>        [,1]   [,2]
#> face -0.626 -0.305
#> body  0.184  1.512
#> tool -0.836  0.390

# Partitions may be named or given by index, and rows always come back in
# effect-space order regardless of how the source was stored.
rownames(relation_block(point, 1L, seq_len(5L)))
#> [1] "face" "body" "tool"
```
