# Declare observation-level confounds and censor facts

Records one row per observation. Censoring is never inferred from motion
or outlier columns: to exclude observations you must name an explicit
logical retain column.

## Usage

``` r
observation_confounds(
  data,
  partition = "partition",
  observation_id = "observation_id",
  censor = NULL,
  provenance = list()
)
```

## Arguments

- data:

  A nonempty data frame with one row per observation.

- partition, observation_id:

  Columns binding rows to observation indexes.

- censor:

  Optional logical censor column. `TRUE` means retained.

- provenance:

  Portable confound provenance.

## Value

An `effect_observation_confounds` fact object: a list with the canonical
`$data` and its `$schema`, the `$partition_column`,
`$observation_id_column` and `$censor_column` role names, `$provenance`,
and a `$confounds_id`.

## See also

[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md),
which joins these rows to the observation axis and records the resulting
row lineage, and
[`bids_confounds()`](https://bbuchsbaum.github.io/crossform/reference/bids_confounds.md)
for fMRIPrep-style TSV input.

Other typed observation facts:
[`bids_study()`](https://bbuchsbaum.github.io/crossform/reference/bids_study.md),
[`observation_events()`](https://bbuchsbaum.github.io/crossform/reference/observation_events.md),
[`observation_index()`](https://bbuchsbaum.github.io/crossform/reference/observation_index.md),
[`observations()`](https://bbuchsbaum.github.io/crossform/reference/observations.md),
[`partition_hierarchy()`](https://bbuchsbaum.github.io/crossform/reference/partition_hierarchy.md)

## Examples

``` r
confounds <- data.frame(
  partition = "run-1",
  observation_id = 1:6,
  framewise_displacement = c(0.1, 0.2, 0.9, 0.1, 0.1, 0.3),
  retained = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)
)

# Naming the censor column is what makes the exclusion a declared fact.
record <- observation_confounds(confounds, censor = "retained")
record$censor_column
#> [1] "retained"
sum(record$data$retained)
#> [1] 5

# Without `censor`, the motion column is kept but nothing is excluded: no
# censor policy is inferred from it.
is.null(observation_confounds(confounds)$censor_column)
#> [1] TRUE
```
