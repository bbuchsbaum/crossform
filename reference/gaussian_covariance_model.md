# Declare a joint Gaussian covariance interpretation

This explicit declaration is required before canonical correlations are
transformed into Gaussian mutual information.

## Usage

``` r
gaussian_covariance_model(provenance = list())
```

## Arguments

- provenance:

  Named metadata describing the fixed model assumption.

## Value

An `effect_gaussian_covariance_model` declaration recording its
`$family`, that it is `$fixed`, the stated `$provenance`, and a
`$signature` carried into the result identity. It performs no fitting
and no goodness-of-fit test.

## See also

[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md)
with `view = "gaussian_information"`, the only place this declaration is
used, and
[`canonical_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
for the correlation spectrum the information is computed from.

Other neural metrics:
[`diagonal_precision()`](https://bbuchsbaum.github.io/crossform/reference/diagonal_precision.md),
[`identity_metric()`](https://bbuchsbaum.github.io/crossform/reference/identity_metric.md),
[`metric_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/metric_capabilities.md),
[`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md),
[`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md),
[`noise_precision()`](https://bbuchsbaum.github.io/crossform/reference/noise_precision.md),
[`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)

## Examples

``` r
# The declaration is required so that a Gaussian information number
# cannot be produced without someone stating the model it rests on.
model <- gaussian_covariance_model(
  list(assumption = "joint Gaussian time observations")
)
model$family
#> [1] "joint_gaussian_covariance"
model$provenance$assumption
#> [1] "joint Gaussian time observations"

# It records the assumption; it does not test it, so the signature is a
# provenance trail rather than evidence of fit.
model$fixed
#> [1] TRUE
```
