# Specify on-demand shrinkage-to-diagonal precision

The local covariance estimator is
`(1 - shrinkage) * S + shrinkage * diag(diag(S))`, followed by the
declared variance and spectral floors. `shrinkage` is fixed by the
recipe; it is not tuned on evaluation effects.

## Usage

``` r
shrinkage_precision(
  shrinkage = 0.1,
  relative_variance_floor = 1e-08,
  absolute_variance_floor = 0,
  relative_spectral_floor = 1e-10,
  domain = NULL
)
```

## Arguments

- shrinkage:

  Fixed number in `(0, 1]`.

- relative_variance_floor:

  Positive floor relative to the mean positive local residual variance.

- absolute_variance_floor:

  Nonnegative floor in squared response units.

- relative_spectral_floor:

  Positive minimum eigenvalue relative to the local covariance scale.

- domain:

  Optional exact neural domain. Omitting it defers domain binding until
  schedule compilation.

## Value

An `effect_metric_recipe` whose `$hyperparameters` record the
`fixed_shrinkage_to_residual_diagonal` estimator, the fixed `shrinkage`,
and the variance and spectral floors. It is not diagonal, so it needs a
dense local support.

## See also

[`diagonal_precision()`](https://bbuchsbaum.github.io/crossform/reference/diagonal_precision.md)
for the diagonal-only recipe,
[`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md)
for which partitions may train it, and
[`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md)
to compile it.

Other neural metrics:
[`diagonal_precision()`](https://bbuchsbaum.github.io/crossform/reference/diagonal_precision.md),
[`gaussian_covariance_model()`](https://bbuchsbaum.github.io/crossform/reference/gaussian_covariance_model.md),
[`identity_metric()`](https://bbuchsbaum.github.io/crossform/reference/identity_metric.md),
[`metric_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/metric_capabilities.md),
[`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md),
[`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md),
[`noise_precision()`](https://bbuchsbaum.github.io/crossform/reference/noise_precision.md)

## Examples

``` r
# The default shrinks the local residual covariance 10% toward its own
# diagonal, which keeps small searchlights invertible.
recipe <- shrinkage_precision()
recipe$hyperparameters$shrinkage
#> [1] 0.1

# Shrinkage is fixed by the recipe, never tuned on evaluation effects, so
# changing it names a different estimand.
shrinkage_precision(0.3)$hyperparameters$shrinkage
#> [1] 0.3

# Unlike a diagonal recipe, this one needs the full local support.
metric_capabilities(recipe)[c("native_diagonal", "support_dense")]
#> $native_diagonal
#> [1] FALSE
#> 
#> $support_dense
#> [1] TRUE
#> 

# Shrinkage must lie in (0, 1]; zero would be an unregularized covariance.
refused <- try(shrinkage_precision(0), silent = TRUE)
conditionMessage(attr(refused, "condition"))
#> [1] "`shrinkage` must be one finite number in (0, 1]."
```
