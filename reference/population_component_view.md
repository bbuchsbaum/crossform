# Build an equal-axis population component view

Converts a
[`population_decomposition()`](https://bbuchsbaum.github.io/crossform/reference/population_decomposition.md)
object into directly inspectable plotted data for total, coherent, and
configuration coefficients. One symmetric axis limit is computed across
every selected panel, so equal signed magnitudes always receive equal
visual magnitude.

## Usage

``` r
population_component_view(x, term, query = NULL)
```

## Arguments

- x:

  An `effect_population_decomposition`.

- term:

  One population model coefficient name.

- query:

  Optional query names. The default keeps every query.

## Value

An `effect_population_component_view` with `data`, `coverage`, one
`axis` contract, and a plotting receipt.
