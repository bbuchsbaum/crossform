# Compute pairing-appropriate signed relation marginals

Compute pairing-appropriate signed relation marginals

## Usage

``` r
pairing_marginals(local, over, mass = 1)
```

## Arguments

- local:

  A measurement-by-effect-by-partition numeric array containing
  spatially weighted relation sums.

- over:

  An `effect_pairing`.

- mass:

  Positive frame mass for each measurement.

## Value

For undirected pairings, a list containing only `endpoint`; for directed
pairings, a list containing `left` and `right`.
