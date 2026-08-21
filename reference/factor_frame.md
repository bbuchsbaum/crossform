# Describe a fixed or locally estimated factor frame

Factor frames use a separate contraction lowering. A locally estimated
factor retains location-dependent fitting and is never treated as an
additive-frame collapse.

## Usage

``` r
factor_frame(
  factors,
  locally_estimated = FALSE,
  domain_id = "abstract",
  domain = NULL
)
```

## Arguments

- factors:

  A nonempty list of numeric factor matrices.

- locally_estimated:

  Whether the factors are estimated independently at each location.

- domain_id:

  Stable identity of the neural feature domain.

- domain:

  Optional exact `effect_domain` or internal domain reference.

## Value

A declarative frame value.
