# Define one partition's observation axis

Records what was acquired, in what order, and on which clock, for a
single partition such as a run. Build one index per response matrix
before calling
[`observations()`](https://bbuchsbaum.github.io/crossform/reference/observations.md).

## Usage

``` r
observation_index(
  observation_id,
  partition,
  time = NULL,
  units = NULL,
  provenance = list()
)
```

## Arguments

- observation_id:

  Unique ordered observation identifiers.

- partition:

  One partition identifier.

- time:

  Optional strictly increasing finite observation times.

- units:

  One physical time unit when `time` is supplied.

- provenance:

  Portable acquisition provenance.

## Value

An `effect_observation_index`: a list with `$observation_id`,
`$partition`, `$time` and `$units` (both `NULL` when untimed), the
`$timing` flag, `$provenance`, and a `$signature`.

## See also

[`observations()`](https://bbuchsbaum.github.io/crossform/reference/observations.md)
to bind indexes to response sources, and
[`observation_events()`](https://bbuchsbaum.github.io/crossform/reference/observation_events.md)
whose clock must use the same units.

Other typed observation facts:
[`bids_study()`](https://bbuchsbaum.github.io/crossform/reference/bids_study.md),
[`observation_confounds()`](https://bbuchsbaum.github.io/crossform/reference/observation_confounds.md),
[`observation_events()`](https://bbuchsbaum.github.io/crossform/reference/observation_events.md),
[`observations()`](https://bbuchsbaum.github.io/crossform/reference/observations.md),
[`partition_hierarchy()`](https://bbuchsbaum.github.io/crossform/reference/partition_hierarchy.md)

## Examples

``` r
# One acquisition axis: six scans on a 2 s clock.
index <- observation_index(
  seq_len(6L), partition = "run-1",
  time = seq(0, by = 2, length.out = 6L), units = "seconds"
)
index$partition
#> [1] "run-1"
index$timing
#> [1] TRUE

# The clock must be strictly increasing, so a repeated acquisition time is
# rejected here rather than silently reordered later.
try(observation_index(1:3, "run-1", time = c(0, 2, 2), units = "seconds"))
#> Error : `time` must be finite, strictly increasing, and match observations.
```
