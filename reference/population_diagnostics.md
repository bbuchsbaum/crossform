# Diagnose population coverage, sink exposure, and transport sensitivity

Produces descriptive diagnostics for a fitted query-bank population
result. It never refits or mutates the primary estimand. Thresholded
rows are explicitly labelled sensitivity summaries, and associations are
descriptive warnings rather than causal corrections.

## Usage

``` r
population_diagnostics(
  x,
  minimum_coverage = 0.8,
  minimum_transport_quality = 0.7,
  material_change = 0.2
)
```

## Arguments

- x:

  An `effect_population_result` returned by
  [`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md).

- minimum_coverage:

  Minimum planned-subject fraction for a warning.

- minimum_transport_quality:

  Minimum retained-territory fraction for a subject to enter the
  sensitivity summary.

- material_change:

  Fraction of contributors removed at which a composition-change warning
  is raised.

## Value

An `effect_population_diagnostics` object with subject provenance, cell
summaries, descriptive associations, sensitivity summaries, and
warnings.
