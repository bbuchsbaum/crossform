# Bind BIDS files into a generic study

BIDS is an adapter boundary, not the crossform object model: the result
is the same
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md)
that can be built from any event and confound source.

## Usage

``` r
bids_study(
  observations,
  event_files,
  confound_files = NULL,
  partitions = observations$partitions,
  observation_ids = NULL,
  censor = NULL,
  hierarchy = NULL,
  units = "seconds"
)
```

## Arguments

- observations:

  A declared
  [`observations()`](https://bbuchsbaum.github.io/crossform/reference/observations.md)
  record.

- event_files:

  Event-TSV paths in observation partition order.

- confound_files:

  Optional confound-TSV paths in the same order.

- partitions:

  Explicit ordered partition identifiers.

- observation_ids:

  Optional confound-row identifiers. Defaults to the identifiers already
  declared by `observations`.

- censor:

  Optional explicit logical retain-column name.

- hierarchy:

  Optional
  [`partition_hierarchy()`](https://bbuchsbaum.github.io/crossform/reference/partition_hierarchy.md).

- units:

  Physical event-time unit.

## Value

A generic
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md)
object, identical in structure to one built from any other event and
confound source, with `$provenance$adapter` recording the BIDS route.

## See also

[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md)
for the object returned,
[`bids_events()`](https://bbuchsbaum.github.io/crossform/reference/bids_events.md)
and
[`bids_confounds()`](https://bbuchsbaum.github.io/crossform/reference/bids_confounds.md)
for the individual file adapters, and
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
for the next step.

Other typed observation facts:
[`observation_confounds()`](https://bbuchsbaum.github.io/crossform/reference/observation_confounds.md),
[`observation_events()`](https://bbuchsbaum.github.io/crossform/reference/observation_events.md),
[`observation_index()`](https://bbuchsbaum.github.io/crossform/reference/observation_index.md),
[`observations()`](https://bbuchsbaum.github.io/crossform/reference/observations.md),
[`partition_hierarchy()`](https://bbuchsbaum.github.io/crossform/reference/partition_hierarchy.md)

## Examples

``` r
set.seed(1)
domain <- abstract_domain(3L, id = "bids-study-example")
index <- observation_index(
  1:4, "run-1", time = seq(0, by = 2, length.out = 4L), units = "seconds"
)
record <- observations(
  list(`run-1` = matrix(rnorm(12), 4L, 3L)), list(`run-1` = index), domain
)

path <- tempfile(fileext = ".tsv")
utils::write.table(
  data.frame(
    onset = c(0, 4), duration = 0.5, trial_type = c("face", "body")
  ),
  path, sep = "\t", row.names = FALSE, quote = FALSE
)

# BIDS stays at the file boundary: what comes back is the ordinary study
# object, with the event clock checked against the acquisition clock.
facts <- bids_study(record, event_files = c(`run-1` = path))
class(facts)
#> [1] "effect_study"
study_capabilities(facts)$timing_resolved
#> [1] TRUE
unlink(path)
```
