# Materialize one metric from an on-demand schedule

This is an inspection operation. Execution kernels should request one
edge-scoped provider and reuse its combined sufficient statistics while
deriving local solve handles support by support.

## Usage

``` r
materialize_metric(schedule, node, edge = 1L)
```

## Arguments

- schedule:

  A compiled metric schedule.

- node:

  One support position or node identifier.

- edge:

  One evaluation-edge position or name.

## Value

A support-local
[`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md).
