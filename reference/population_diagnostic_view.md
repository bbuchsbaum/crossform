# Bind population effects to coverage and transport diagnostics

Bind population effects to coverage and transport diagnostics

## Usage

``` r
population_diagnostic_view(effect_view, diagnostics, node = NULL, query = NULL)
```

## Arguments

- effect_view:

  An `effect_population_component_view` or
  `effect_population_scale_profile`.

- diagnostics:

  An `effect_population_diagnostics` for the total result underlying the
  same decomposition.

- node, query:

  Optional synchronized selections.

## Value

An `effect_population_diagnostic_view` with exact-key effect and support
rows, warnings, filter provenance, and separate diagnostic axes.
