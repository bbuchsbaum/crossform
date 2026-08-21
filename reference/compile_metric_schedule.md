# Compile a frozen-provenance, on-demand metric schedule

Compilation binds a metric recipe to canonical residual sufficient
statistics, spatial supports, and evaluation edges. What is frozen is
the estimator identity, training assignment, source revisions, and
regularization policy. Local covariance and precision matrices are
derived when requested and are never retained as a node-by-fold table.

## Usage

``` r
compile_metric_schedule(
  recipe,
  statistics,
  at,
  over,
  training = metric_training_policy("exclude_evaluation")
)
```

## Arguments

- recipe:

  A recipe from
  [`identity_metric()`](https://bbuchsbaum.github.io/crossform/reference/identity_metric.md),
  [`diagonal_precision()`](https://bbuchsbaum.github.io/crossform/reference/diagonal_precision.md),
  or
  [`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md).

- statistics:

  Canonical
  [`residual_pair_statistics()`](https://bbuchsbaum.github.io/crossform/reference/residual_pair_statistics.md).

- at:

  The same support-index-backed frame used to accumulate `statistics`.

- over:

  An evaluation
  [`pairing()`](https://bbuchsbaum.github.io/crossform/reference/pairing.md).

- training:

  A
  [`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md).

## Value

An `effect_frozen_metric_schedule`.
