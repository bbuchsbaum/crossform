# Construct a fixed neural noise-precision metric

`noise_precision()` is the semantically specific constructor for a known
or externally estimated precision that is fixed before effect
evaluation. It uses the same
[`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md)
representation while recording that the operator has
inverse-noise-covariance meaning.

## Usage

``` r
noise_precision(
  value,
  domain,
  support = NULL,
  covariance = NULL,
  tolerance = 1e-10,
  provenance = list()
)
```

## Arguments

- value:

  A finite symmetric positive-semidefinite precision matrix.

- domain:

  Exact neural feature domain.

- support:

  Ordered support identities for a local matrix.

- covariance:

  Optional retained inverse covariance.

- tolerance:

  Positive numerical validation tolerance.

- provenance:

  Additional compact provenance. Reserved semantic fields cannot be
  replaced.

## Value

An `effect_neural_metric` whose `$provenance$metric_role` is
`"noise_precision"`, carrying the canonicalized `$value`, its `$domain`
and `$support`, `$capabilities`, and a `$signature`. The recorded role
is what lets
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md)
name the result a Mahalanobis distance.

## See also

[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md)
and
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
(its `metric` argument);
[`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)
when the precision must instead be learned from residuals.

Other neural metrics:
[`diagonal_precision()`](https://bbuchsbaum.github.io/crossform/reference/diagonal_precision.md),
[`gaussian_covariance_model()`](https://bbuchsbaum.github.io/crossform/reference/gaussian_covariance_model.md),
[`identity_metric()`](https://bbuchsbaum.github.io/crossform/reference/identity_metric.md),
[`metric_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/metric_capabilities.md),
[`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md),
[`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md),
[`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)

## Examples

``` r
# A precision estimated outside crossform and fixed before evaluation.
domain <- abstract_domain(3, id = "noise-precision-example")
covariance <- matrix(c(1, 0.4, 0.1, 0.4, 1, 0.3, 0.1, 0.3, 1), 3, 3)
metric <- noise_precision(
  solve(covariance), domain, covariance = covariance,
  provenance = list(source = "resting-run residual covariance")
)
metric$provenance$metric_role
#> [1] "noise_precision"
metric$capabilities$positive_definite
#> [1] TRUE

# The role is what a crossnobis plan checks for, so a plain
# `neural_metric()` with the same numbers is deliberately not equivalent.
identical(
  neural_metric(solve(covariance), domain)$provenance$metric_role,
  metric$provenance$metric_role
)
#> [1] FALSE

# The semantic provenance fields are reserved and cannot be overwritten.
refused <- try(
  noise_precision(diag(3), domain, provenance = list(metric_role = "other")),
  silent = TRUE
)
conditionMessage(attr(refused, "condition"))
#> [1] "Noise-precision semantic provenance fields are reserved."
```
