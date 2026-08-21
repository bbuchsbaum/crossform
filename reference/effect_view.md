# Construct an explicit query-only effect view

An `effect_view` contains derived values but not the geometry from which
they arose. It is never returned where a complete `effect_geometry` is
promised.

## Usage

``` r
effect_view(
  values,
  query,
  component,
  receipt,
  index = NULL,
  metadata = list(),
  effects = NULL,
  left_space = NULL,
  right_space = NULL
)
```

## Arguments

- values:

  A finite numeric measurement-by-view matrix.

- query:

  The compiled query matrix or a descriptive query value.

- component:

  Geometry component used to create the view.

- receipt:

  The
  [`execution_receipt()`](https://bbuchsbaum.github.io/crossform/reference/execution_receipt.md)
  proving how the view was made.

- index:

  Optional measurement index.

- metadata:

  Optional compact semantic metadata.

- effects:

  Compatibility self-space binding. Prefer `left_space` and
  `right_space` for rectangular views.

- left_space, right_space:

  Optional ordered effect-space bindings.

## Value

A query-only `effect_view`.
