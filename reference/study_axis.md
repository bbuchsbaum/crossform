# Select a typed study axis

This accessor supplies vocabulary only. The returned object does not
select a generalization target or assert independence.

## Usage

``` r
study_axis(x, name)
```

## Arguments

- x:

  An
  [`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md).

- name:

  One hierarchy axis name.

## Value

An `effect_study_axis`: a list with the axis `$name`, its distinct
`$levels`, the enclosing `$parent` axis (`NULL` at the outermost axis),
and a `$signature` tied to the study hierarchy.

## See also

[`partition_hierarchy()`](https://bbuchsbaum.github.io/crossform/reference/partition_hierarchy.md)
which declares the available axes, and
[`cross_partitions()`](https://bbuchsbaum.github.io/crossform/reference/cross_partitions.md),
where a generalization axis is actually chosen.

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
[`study_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/study_capabilities.md)

## Examples

``` r
set.seed(1)
domain <- abstract_domain(4L, id = "study-axis-example")
partitions <- c("run-1", "run-2")
indexes <- lapply(partitions, function(partition) {
  observation_index(1:6, partition)
})
names(indexes) <- partitions
record <- observations(
  lapply(indexes, function(index) matrix(rnorm(24), 6L, 4L)), indexes, domain
)
facts <- study(record, hierarchy = partition_hierarchy(data.frame(
  partition = partitions, session = c("ses-1", "ses-2")
)))

# Vocabulary only: naming an axis does not assert independence across it.
axis <- study_axis(facts, "session")
axis
#> study_axis<session; 2 levels>
axis$levels
#> [1] "ses-1" "ses-2"

# An axis that was never declared is refused rather than invented.
catch_refusal(study_axis(facts, "subject"))$capability
#> [1] "declared_generalization_axis"
```
