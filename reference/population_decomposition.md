# Enforce the population coefficient decomposition law

Checks `beta_total = beta_coherent + beta_configuration` for every
estimable population coefficient, node, and query. The three results
must come from the same plan, model, query bank, and cellwise subject
sets.

## Usage

``` r
population_decomposition(
  total,
  coherent,
  configuration,
  estimator = c("HC3", "classical"),
  tolerance = 1e-10
)
```

## Arguments

- total, coherent, configuration:

  Comparable `effect_population_result` objects for the named
  components.

- estimator:

  `"HC3"` or `"classical"` for the joint component coefficient
  covariance.

- tolerance:

  Absolute numerical tolerance for the conservation law.

## Value

An `effect_population_decomposition` carrying the coefficient gap, joint
component covariance, uncertainty objects, and comparability receipt.
