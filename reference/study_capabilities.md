# Inspect factual study capabilities

Reports what the binding actually established, so downstream calls can
require a guarantee instead of assuming it.

## Usage

``` r
study_capabilities(x)
```

## Arguments

- x:

  An
  [`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md).

## Value

A one-row data frame with logical columns `aligned_observations`,
`timing_resolved`, `partition_hierarchy`, and `stable_source_revision`.

## See also

[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md),
and
[`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md)
for the equivalent report on a compiled relation plan.

Other studies and effect maps:
[`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md),
[`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md),
[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md),
[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md),
[`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md),
[`lower_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/lower_effect_map.md),
[`observation_model()`](https://bbuchsbaum.github.io/crossform/reference/observation_model.md),
[`raw_design_model()`](https://bbuchsbaum.github.io/crossform/reference/raw_design_model.md),
[`raw_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/raw_effect_map.md),
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md),
[`study_axis()`](https://bbuchsbaum.github.io/crossform/reference/study_axis.md)

## Examples

``` r
set.seed(1)
domain <- abstract_domain(4L, id = "study-capabilities-example")
index <- observation_index(
  1:6, "run-1", time = seq(0, 10, by = 2), units = "seconds"
)
record <- observations(
  list(`run-1` = matrix(rnorm(24), 6L, 4L)), list(`run-1` = index), domain
)

# With no event record there is no timing to resolve, and the report says so
# rather than defaulting to TRUE.
study_capabilities(study(record))
#>   aligned_observations timing_resolved partition_hierarchy
#> 1                 TRUE           FALSE                TRUE
#>   stable_source_revision
#> 1                   TRUE

# Binding timed events on the same clock earns `timing_resolved`.
events <- observation_events(data.frame(
  partition = "run-1", event_id = c("e1", "e2"), onset = c(0, 6),
  duration = 0.5
))
study_capabilities(study(record, events))$timing_resolved
#> [1] TRUE
```
