# Declare source execution capabilities

`source_capabilities()` states what an out-of-memory relation source can
actually do, so the compiler can choose a bounded execution route and
record the source revision in the receipt. Declare it when supplying a
custom source such as
[`file_matrix_source()`](https://bbuchsbaum.github.io/crossform/reference/file_matrix_source.md).

## Usage

``` r
source_capabilities(
  block_read,
  reopenable = FALSE,
  thread_safe = FALSE,
  stable_revision
)
```

## Arguments

- block_read:

  Whether bounded feature-block reads are supported.

- reopenable:

  Whether a fresh read-only handle can be opened safely.

- thread_safe:

  Whether concurrent reads through one handle are supported.

- stable_revision:

  A strong immutable source revision or checksum.

## Value

An `effect_source_capabilities` value with the three logical flags
`$block_read`, `$reopenable`, `$thread_safe`, and the `$stable_revision`
identifier bound into the execution receipt.

## See also

[`file_matrix_source()`](https://bbuchsbaum.github.io/crossform/reference/file_matrix_source.md),
which carries these capabilities, and
[`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md),
which is checked against them before execution.

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
[`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md)

## Examples

``` r
# A block-readable, reopenable source identified by a content checksum.
revision <- paste0("sha256:", paste(rep("a", 64), collapse = ""))
capabilities <- source_capabilities(
  block_read = TRUE, reopenable = TRUE, stable_revision = revision
)
capabilities$block_read
#> [1] TRUE
capabilities$thread_safe
#> [1] FALSE

# The revision must be a strong identifier: a mutable label such as a file
# modification time is refused, because receipts must stay verifiable.
weak <- try(
  source_capabilities(TRUE, stable_revision = "2026-08-15"), silent = TRUE
)
conditionMessage(attr(weak, "condition"))
#> [1] "`stable_revision` must be a sha256 identifier with 64 hexadecimal digits."
```
