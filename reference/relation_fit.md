# Attach statistical error information to a relation

`relation_fit()` keeps the mathematical `effect_relation` intact and
adds a separate, capability-bearing error channel. A relation can remain
valid without this envelope; operations that learn a neural metric or
calibrate within-participant uncertainty must require the corresponding
capability.

## Usage

``` r
relation_fit(relation, error_models = NULL, provenance = list())
```

## Arguments

- relation:

  A validated `effect_relation`.

- error_models:

  Optional named list of internal error-model values, one per relation
  partition. Omitting it records an explicitly absent channel.

- provenance:

  Compact fit-level provenance.

## Value

An `effect_relation_fit`: a list with the untouched `$relation`, one
`$error_models` entry per partition (`NULL` where absent), the derived
per-partition `$capabilities`, `$provenance`, and a `$signature` binding
the relation and every error model.

## See also

[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md),
the route that installs a real residual channel, and
[`relation_fit_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit_capabilities.md)
to inspect the result.

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
[`relation_fit_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit_capabilities.md),
[`relation_plan_receipts()`](https://bbuchsbaum.github.io/crossform/reference/relation_plan_receipts.md),
[`residual_block()`](https://bbuchsbaum.github.io/crossform/reference/residual_block.md),
[`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md),
[`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)

## Examples

``` r
set.seed(1)
domain <- abstract_domain(4L, id = "relation-fit-example")
betas <- matrix(rnorm(8), 2L, 4L,
  dimnames = list(c("face", "body"), NULL))
point <- relation(list(`run-1` = betas), domain = domain)

# Wrapping a point relation records an explicitly absent error channel,
# which is different from an unstated one.
fit <- relation_fit(point)
relation_fit_capabilities(fit)
#>   partition error_model residual_blocks effect_covariance residual_df
#> 1     run-1       FALSE           FALSE             FALSE       FALSE
#>   separable_error learned_metric_input within_participant_calibration
#> 1           FALSE                FALSE                          FALSE

# The relation itself is unchanged, so point geometry still works.
identical(fit$relation, point)
#> [1] TRUE
```
