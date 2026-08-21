# Construct a semantically complete effect geometry

An `effect_geometry` is the symmetric packed self-form specialization of
a complete `effect_form`. Its storage may be in memory or block-backed,
while configuration remains exactly total minus coherent.

## Usage

``` r
effect_geometry(
  total,
  coherent,
  marginals,
  effects,
  receipt,
  index = NULL,
  metadata = list()
)
```

## Arguments

- total, coherent:

  Packed geometry matrices, with measurements in rows, or internal
  geometry stores having the same dimensions.

- marginals:

  Pairing-appropriate signed marginals. Undirected pairings contain
  `endpoint`; directed pairings contain `left` and `right`.

- effects:

  An
  [`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md)
  whose dimension must match the triangular packed-geometry width and
  marginal columns. Unique names are accepted as shorthand for an
  unspecified-basis space.

- receipt:

  The
  [`execution_receipt()`](https://bbuchsbaum.github.io/crossform/reference/execution_receipt.md)
  proving how the result was made.

- index:

  Optional measurement index with one entry per geometry row.

- metadata:

  Optional compact semantic metadata.

## Value

A complete `effect_geometry` and `effect_form`.
