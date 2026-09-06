# Diagnose leave-one-subject population influence

Computes descriptive leave-one-subject coefficient changes from the
exact retained transported values in a
[`population_decomposition()`](https://bbuchsbaum.github.io/crossform/reference/population_decomposition.md).
It also binds optional cross-fitted heterogeneity output and
transport/coverage provenance.

## Usage

``` r
population_influence(
  x,
  heterogeneity = NULL,
  mode = c("bounded", "deep"),
  max_subjects = 50L,
  max_cells = 250000L
)
```

## Arguments

- x:

  An `effect_population_decomposition`.

- heterogeneity:

  Optional `effect_population_heterogeneity` from the same planned
  subjects, normally with `estimator = "cross_fit"`.

- mode:

  `"bounded"` refuses work beyond `max_subjects` or `max_cells`;
  `"deep"` is the explicit opt-in for a larger run.

- max_subjects, max_cells:

  Positive bounds for the default mode.

## Value

An `effect_population_influence` with a long influence table, subject
provenance, optional heterogeneity summary, and descriptive scope.
