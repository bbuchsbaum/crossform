# Construct a support-local same-space neural metric

A neural metric is the PSD same-space role of the evidence-pairing
operator `K`. It is distinct from a cross-space measurement bridge. The
stored matrix has width equal to its local support, never the full
neural domain unless the metric is genuinely global. Supplying `inverse`
retains a mathematically equivalent inverse action, such as the
covariance from which a precision metric was derived, so inverse
quadratic forms need no second factorization.

## Usage

``` r
neural_metric(
  value,
  domain,
  support = NULL,
  inverse = NULL,
  estimation = c("fixed", "learned_frozen"),
  tolerance = 1e-10,
  provenance = list()
)
```

## Arguments

- value:

  A finite symmetric positive-semidefinite local metric matrix.

- domain:

  Exact neural feature domain.

- support:

  Ordered feature identities for the local matrix coordinates.

- inverse:

  Optional finite symmetric inverse metric.

- estimation:

  Whether the materialized metric was fixed a priori or learned and
  frozen before evaluation.

- tolerance:

  Positive numerical validation tolerance.

- provenance:

  Compact metric provenance. Learned-frozen metrics require
  `frozen = TRUE` and a strong `training_signature`.

## Value

An `effect_neural_metric` carrying the canonicalized `$value`, its
`$domain` and `$support`, the `$estimation` status, any retained
`$inverse_representation`, the derived `$capabilities` record, and a
content-addressed `$signature` bound into every plan that uses it.

## See also

[`noise_precision()`](https://bbuchsbaum.github.io/crossform/reference/noise_precision.md)
for the semantically specific inverse-noise-covariance constructor,
[`metric_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/metric_capabilities.md)
to inspect what a metric admits, and
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
to compile it into a plan.

Other neural metrics:
[`diagonal_precision()`](https://bbuchsbaum.github.io/crossform/reference/diagonal_precision.md),
[`gaussian_covariance_model()`](https://bbuchsbaum.github.io/crossform/reference/gaussian_covariance_model.md),
[`identity_metric()`](https://bbuchsbaum.github.io/crossform/reference/identity_metric.md),
[`metric_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/metric_capabilities.md),
[`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md),
[`noise_precision()`](https://bbuchsbaum.github.io/crossform/reference/noise_precision.md),
[`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)

## Examples

``` r
# A fixed diagonal metric: each feature is reweighted independently, so
# the metric stays feature-additive and survives searchlight restriction.
domain <- abstract_domain(3, id = "metric-example")
metric <- neural_metric(diag(c(1, 2, 3)), domain)
metric_capabilities(metric)$feature_additive
#> [1] TRUE

# Retaining the inverse lets inverse quadratic forms be read without a
# second factorization.
full <- neural_metric(
  matrix(c(2, 0.5, 0, 0.5, 2, 0.5, 0, 0.5, 2), 3, 3), domain,
  inverse = solve(matrix(c(2, 0.5, 0, 0.5, 2, 0.5, 0, 0.5, 2), 3, 3))
)
full$inverse_representation$kind
#> [1] "retained_inverse_metric"

# A same-space metric must be positive semidefinite.
indefinite <- try(
  neural_metric(diag(c(1, -1, 1)), domain), silent = TRUE
)
conditionMessage(attr(indefinite, "condition"))
#> [1] "A same-space neural metric must be positive semidefinite."
```
