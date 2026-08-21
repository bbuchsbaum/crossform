# Identify an oriented coherent neural functional

The vector `a` defines the scientifically oriented amplitude `B a`. It
is separate from the metric that normalizes the energy of that
amplitude.

## Usage

``` r
coherent_functional(
  value,
  domain,
  support = NULL,
  label = "raw_weighted_mean",
  provenance = list()
)
```

## Arguments

- value:

  One finite nonzero local covector.

- domain:

  Exact neural feature domain.

- support:

  Ordered feature identities for `value`.

- label:

  Nonempty scientific identity for the functional.

- provenance:

  Compact provenance.

## Value

An immutable `effect_coherent_functional`.
