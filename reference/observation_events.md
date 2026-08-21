# Declare a typed event record

Records the experimental events as facts, keeping every supplied column.
It assigns no model roles: whether a column is a target, a nuisance
term, or unused is decided later by the design model.

## Usage

``` r
observation_events(
  data,
  partition = "partition",
  event_id = "event_id",
  onset = "onset",
  duration = "duration",
  units = "seconds",
  provenance = list()
)
```

## Arguments

- data:

  A nonempty event data frame.

- partition:

  Column naming the observation partition.

- event_id:

  Column containing event identifiers.

- onset, duration:

  Optional timing columns. Supply both or neither.

- units:

  Physical time unit when timing columns are supplied.

- provenance:

  Portable event provenance.

## Value

An `effect_events` fact object: a list with the canonical `$data` and
its `$schema`, the `$partition_column`, `$event_id_column`,
`$onset_column` and `$duration_column` role names, `$units`, the
`$timing` flag, `$provenance`, and an `$events_id`. Column roles remain
model-relative.

## See also

[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md)
to bind events to observations,
[`bids_events()`](https://bbuchsbaum.github.io/crossform/reference/bids_events.md)
to read them from BIDS TSV files, and
[`observation_confounds()`](https://bbuchsbaum.github.io/crossform/reference/observation_confounds.md)
for the observation-level table.

Other typed observation facts:
[`bids_study()`](https://bbuchsbaum.github.io/crossform/reference/bids_study.md),
[`observation_confounds()`](https://bbuchsbaum.github.io/crossform/reference/observation_confounds.md),
[`observation_index()`](https://bbuchsbaum.github.io/crossform/reference/observation_index.md),
[`observations()`](https://bbuchsbaum.github.io/crossform/reference/observations.md),
[`partition_hierarchy()`](https://bbuchsbaum.github.io/crossform/reference/partition_hierarchy.md)

## Examples

``` r
events <- data.frame(
  partition = rep(c("run-1", "run-2"), each = 2L),
  event_id = paste0("e", 1:4),
  onset = c(0, 6, 0, 6),
  duration = 0.5,
  condition = c("face", "body", "face", "body")
)
record <- observation_events(events)
record$timing
#> [1] TRUE
record$units
#> [1] "seconds"

# `condition` is preserved but carries no model role yet; the design model,
# not this record, decides what is a target or a nuisance term.
names(record$data)
#> [1] "partition" "event_id"  "onset"     "duration"  "condition"

# Event ids must be unique within a partition.
duplicated_ids <- events
duplicated_ids$event_id <- c("e1", "e1", "e3", "e4")
try(observation_events(duplicated_ids))
#> Error : Event identifiers must be unique within partition.
```
