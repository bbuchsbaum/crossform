# Build a population scale profile with explicit pointwise uncertainty

Build a population scale profile with explicit pointwise uncertainty

## Usage

``` r
population_scale_profile(
  x,
  term,
  query = NULL,
  interval = c("HC3", "classical", "wild_bootstrap"),
  bootstrap = NULL
)
```

## Arguments

- x:

  An `effect_population_decomposition`.

- term:

  One population coefficient name.

- query:

  Optional query names.

- interval:

  Explicit interval method: `"HC3"`, `"classical"`, or
  `"wild_bootstrap"`.

- bootstrap:

  For `interval = "wild_bootstrap"`, a named list containing total,
  coherent, and configuration bootstrap objects made from the three
  component results and the same coefficient contrast.

## Value

An `effect_population_scale_profile`. Bands are pointwise; no
simultaneous or maxT interpretation is supplied.
