# Construct a lazy experimental-neural relation

Each source is read only by neural feature block. Extractors map its
observation rows into one common named experimental space. Omitting
`extract` declares that sources already contain effect-by-feature
matrices. A relation made from precomputed effects supports point
evidence but carries no residual error channel. Beta matrices alone
cannot recover within- participant analytic uncertainty. If later work
needs a learned neural metric or
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md),
start from raw responses with
[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md).
Subject-level resampling answers a population question; it does not
recreate discarded first-level residual information.

## Usage

``` r
relation(
  sources,
  extract = NULL,
  effects = NULL,
  source_dims = NULL,
  partitions = NULL,
  domain = NULL,
  domain_id = "abstract",
  capabilities = NULL,
  provenance = list()
)
```

## Arguments

- sources:

  A matrix, or a named list of matrix, function, or
  `effect_source_descriptor` response sources.

- extract:

  NULL, one `effect_extractor`, or one extractor per partition.

- effects:

  An
  [`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md)
  for already estimated sources, or unique names used as shorthand for
  an unspecified-basis effect space.

- source_dims:

  Required dimensions for function sources, as one two-element vector
  per partition.

- partitions:

  Optional partition names when not supplied by `sources`.

- domain:

  Optional `effect_domain`; its identity is recorded without resampling
  or changing neural features.

- domain_id:

  Stable neural-domain identity.

- capabilities:

  Optional source-capability values, one per partition.

- provenance:

  Optional provenance metadata.

## Value

An `effect_relation`: a list with the compiled `$sources`, one
`$extractors` entry per partition, the shared `$effect_space` and its
`$effects` labels, `$partitions`, `$n_features`, the `$domain` reference
and `$domain_id`, optional source `$capabilities`, and `$provenance`. No
neural values have been read.

## Structure

A relation declares what will be read, so its elements are identities
and shapes rather than values.

- `$effects`: the effect labels, in the order every partition estimates
  them and the order unnamed contrast weights are read in.

- `$effect_space`: the shared
  [`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md)
  those labels coordinate, carrying the basis identity, units, and
  signature every partition agrees on.

- `$partitions`: the partition names, in the order they were supplied.
  These are what
  [`cross_partitions()`](https://bbuchsbaum.github.io/crossform/reference/cross_partitions.md)
  and
  [`pairing()`](https://bbuchsbaum.github.io/crossform/reference/pairing.md)
  name.

- `$n_features`: the number of neural features every partition spans.

- `$domain`: the neural domain the source columns are bound to, carrying
  its `id`, `n_features`, and `feature_ids`; `$domain_id` repeats that
  `id`.

- `$capabilities`: one
  [`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)
  per partition when the sources declared them, otherwise `NULL`.

- `$provenance`: the metadata supplied at construction, unchanged.

`$sources` and `$extractors` are the compiled read path
[`relation_block()`](https://bbuchsbaum.github.io/crossform/reference/relation_block.md)
uses; they and any element not listed here are internal and may change.

## See also

[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md)
when raw responses are available and a residual channel is needed,
[`relation_block()`](https://bbuchsbaum.github.io/crossform/reference/relation_block.md)
to read one block, and
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
for the next step.

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
[`relation_block()`](https://bbuchsbaum.github.io/crossform/reference/relation_block.md),
[`relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit.md),
[`relation_fit_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit_capabilities.md),
[`relation_plan_receipts()`](https://bbuchsbaum.github.io/crossform/reference/relation_plan_receipts.md),
[`residual_block()`](https://bbuchsbaum.github.io/crossform/reference/residual_block.md),
[`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md),
[`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)

## Examples

``` r
# Precomputed condition betas: one effect-by-feature matrix per run.
set.seed(1)
domain <- abstract_domain(5L, id = "relation-example")
betas <- lapply(c("run-1", "run-2"), function(partition) {
  matrix(rnorm(15), 3L, 5L,
    dimnames = list(c("face", "body", "tool"), NULL))
})
names(betas) <- c("run-1", "run-2")

point <- relation(betas, domain = domain)
point$effects
#> [1] "face" "body" "tool"
point$partitions
#> [1] "run-1" "run-2"

# Betas alone carry no residual information, so anything needing the error
# channel refuses here rather than inventing a standard error.
catch_refusal(residual_df(point, "run-1"))$capability
#> [1] "residual_df"
```
