# Bind observations, events, confounds, clocks, and partition axes

`study()` records factual correspondence only. It does not assign model
roles, choose a sampling unit, infer independence, or select a
`generalizes_over` axis.

## Usage

``` r
study(
  observations,
  events = NULL,
  confounds = NULL,
  hierarchy = NULL,
  clock_tolerance = 0,
  provenance = list()
)
```

## Arguments

- observations:

  An
  [`observations()`](https://bbuchsbaum.github.io/crossform/reference/observations.md)
  object.

- events:

  Optional
  [`observation_events()`](https://bbuchsbaum.github.io/crossform/reference/observation_events.md)
  record.

- confounds:

  Optional
  [`observation_confounds()`](https://bbuchsbaum.github.io/crossform/reference/observation_confounds.md)
  record.

- hierarchy:

  Optional
  [`partition_hierarchy()`](https://bbuchsbaum.github.io/crossform/reference/partition_hierarchy.md).
  A leaf-only hierarchy is constructed from observation partitions when
  omitted.

- clock_tolerance:

  Nonnegative finite tolerance for clock coverage.

- provenance:

  Portable binding provenance.

## Value

An `effect_study`: a list with the validated `$observations`, `$events`,
`$confounds` and `$hierarchy`, the ordered `$partitions`, a
per-partition `$lineage` data frame recording which observations
censoring retained, `$clock_coverage`, `$clock_tolerance`,
`$capabilities`, `$provenance`, and a `$study_id`.

## See also

[`observations()`](https://bbuchsbaum.github.io/crossform/reference/observations.md),
[`observation_events()`](https://bbuchsbaum.github.io/crossform/reference/observation_events.md),
[`observation_confounds()`](https://bbuchsbaum.github.io/crossform/reference/observation_confounds.md),
[`partition_hierarchy()`](https://bbuchsbaum.github.io/crossform/reference/partition_hierarchy.md)
for the inputs;
[`study_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/study_capabilities.md)
to inspect what was established; and
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
for the next step.

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
[`study_axis()`](https://bbuchsbaum.github.io/crossform/reference/study_axis.md),
[`study_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/study_capabilities.md)

## Examples

``` r
set.seed(1)
domain <- abstract_domain(4L, id = "study-example")
indexes <- list(
  `run-1` = observation_index(
    1:6, "run-1", time = seq(0, 10, by = 2), units = "seconds"
  ),
  `run-2` = observation_index(
    1:6, "run-2", time = seq(0, 10, by = 2), units = "seconds"
  )
)
record <- observations(
  lapply(indexes, function(index) matrix(rnorm(24), 6L, 4L)), indexes, domain
)
events <- observation_events(data.frame(
  partition = rep(names(indexes), each = 2L),
  event_id = paste0("e", 1:4), onset = c(0, 6, 0, 6), duration = 0.5,
  condition = c("face", "body", "face", "body")
))

# Binding checks correspondence only; it assigns no model roles.
facts <- study(record, events)
facts
#> study<2 partitions; timing resolved>
study_capabilities(facts)
#>   aligned_observations timing_resolved partition_hierarchy
#> 1                 TRUE            TRUE                TRUE
#>   stable_source_revision
#> 1                   TRUE

# An event outside the acquired window is refused, not silently clipped.
late <- observation_events(data.frame(
  partition = "run-1", event_id = "e5", onset = 40, duration = 1
))
catch_refusal(study(record, late))$capability
#> [1] "timing_resolved"
```
