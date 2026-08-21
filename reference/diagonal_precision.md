# Specify on-demand diagonal residual-variance precision

The local metric is the inverse of the residual variance of each
feature, floored so that a near-silent feature cannot dominate. Use it
for univariate noise normalization when a full local covariance is not
wanted or not estimable.

## Usage

``` r
diagonal_precision(
  relative_variance_floor = 1e-08,
  absolute_variance_floor = 0,
  domain = NULL
)
```

## Arguments

- relative_variance_floor:

  Positive floor relative to the mean positive local residual variance.

- absolute_variance_floor:

  Nonnegative floor in squared response units.

- domain:

  Optional exact neural domain. Omitting it defers domain binding until
  schedule compilation.

## Value

An `effect_metric_recipe` whose `$hyperparameters` record the
`residual_diagonal_inverse` estimator and both variance floors, with
`$capabilities$native_diagonal` true and `$capabilities$materialized`
false until a schedule binds it.

## See also

[`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)
for a full local covariance shrunk to its diagonal, and
[`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md),
which consumes the recipe.

Other neural metrics:
[`gaussian_covariance_model()`](https://bbuchsbaum.github.io/crossform/reference/gaussian_covariance_model.md),
[`identity_metric()`](https://bbuchsbaum.github.io/crossform/reference/identity_metric.md),
[`metric_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/metric_capabilities.md),
[`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md),
[`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md),
[`noise_precision()`](https://bbuchsbaum.github.io/crossform/reference/noise_precision.md),
[`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)

## Examples

``` r
# The estimator and both floors are fixed by the recipe, so they are part
# of the plan identity rather than a runtime tuning choice.
recipe <- diagonal_precision(relative_variance_floor = 1e-6)
recipe$kind
#> [1] "diagonal_variance_precision"
recipe$hyperparameters[c("estimator", "relative_variance_floor")]
#> $estimator
#> [1] "residual_diagonal_inverse"
#> 
#> $relative_variance_floor
#> [1] 1e-06
#> 

# The floor is a guard against dividing by a near-zero residual variance.
diagonal_precision(absolute_variance_floor = 0.01)$
  hyperparameters$absolute_variance_floor
#> [1] 0.01

# A nonpositive relative floor would remove that guard, so it is refused.
refused <- try(diagonal_precision(relative_variance_floor = 0), silent = TRUE)
conditionMessage(attr(refused, "condition"))
#> [1] "Metric variance floors must be finite, with a positive relative floor and a nonnegative absolute floor."
```
