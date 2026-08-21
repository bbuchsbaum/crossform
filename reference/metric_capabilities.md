# Inspect exact neural-metric capabilities

`metric_capabilities()` answers what a metric or recipe admits before a
plan provokes a refusal: whether it is diagonal, feature-additive,
positive definite, already materialized, and whether inverse quadratic
forms are available.

## Usage

``` r
metric_capabilities(x)
```

## Arguments

- x:

  A
  [`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md)
  or an on-demand metric recipe. A cross-space measurement bridge is
  refused because cross-space bridges are not same-space metrics.

## Value

An `effect_metric_capabilities` record of logical flags, including
`$identity`, `$native_diagonal`, `$feature_additive`, `$support_dense`,
`$positive_definite`, `$fixed`, `$learned_recipe`, `$materialized`, and
`$inverse_quadratic`.

## See also

[`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md)
and
[`noise_precision()`](https://bbuchsbaum.github.io/crossform/reference/noise_precision.md)
for fixed metrics;
[`identity_metric()`](https://bbuchsbaum.github.io/crossform/reference/identity_metric.md),
[`diagonal_precision()`](https://bbuchsbaum.github.io/crossform/reference/diagonal_precision.md),
and
[`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)
for on-demand recipes.

Other neural metrics:
[`diagonal_precision()`](https://bbuchsbaum.github.io/crossform/reference/diagonal_precision.md),
[`gaussian_covariance_model()`](https://bbuchsbaum.github.io/crossform/reference/gaussian_covariance_model.md),
[`identity_metric()`](https://bbuchsbaum.github.io/crossform/reference/identity_metric.md),
[`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md),
[`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md),
[`noise_precision()`](https://bbuchsbaum.github.io/crossform/reference/noise_precision.md),
[`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)

## Examples

``` r
# A materialized diagonal metric is feature-additive, so it restricts
# cleanly to any searchlight support.
domain <- abstract_domain(3, id = "capability-example")
fixed <- metric_capabilities(neural_metric(diag(c(1, 2, 3)), domain))
fixed[c("native_diagonal", "feature_additive", "materialized")]
#> $native_diagonal
#> [1] TRUE
#> 
#> $feature_additive
#> [1] TRUE
#> 
#> $materialized
#> [1] TRUE
#> 

# A recipe is not yet a matrix: it promises a local metric per support.
recipe <- metric_capabilities(shrinkage_precision(0.2))
recipe[c("learned_recipe", "materialized", "support_dense")]
#> $learned_recipe
#> [1] TRUE
#> 
#> $materialized
#> [1] FALSE
#> 
#> $support_dense
#> [1] TRUE
#> 

# Anything that is not a metric or a recipe is refused rather than
# coerced, so a bare matrix or a cross-space object cannot pass as one.
refused <- try(metric_capabilities(diag(c(1, 2, 3))), silent = TRUE)
conditionMessage(attr(refused, "condition"))
#> [1] "`x` must be a neural metric or metric recipe."
```
