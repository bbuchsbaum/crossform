# Bind raw observation sources to indexes and a neural domain

The constructor validates only source metadata and axes. Neural values
are read later by
[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md).

## Usage

``` r
observations(
  sources,
  index,
  domain,
  source_dims = NULL,
  partitions = NULL,
  capabilities = NULL,
  provenance = list()
)
```

## Arguments

- sources:

  A matrix or named list of matrix, function, or
  `effect_source_descriptor` sources.

- index:

  One
  [`observation_index()`](https://bbuchsbaum.github.io/crossform/reference/observation_index.md)
  or a named list, one per partition.

- domain:

  The exact neural domain.

- source_dims:

  Dimensions required for function sources.

- partitions:

  Optional partition order.

- capabilities:

  Optional
  [`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)
  per partition. Function sources require explicit capabilities.

- provenance:

  Portable observation provenance.

## Value

An `effect_observations` fact object: a list with the compiled
`$sources`, `$indexes`, `$partitions`, the `$domain` reference,
`$n_features`, per-partition `$capabilities`, `$provenance`, and an
`$observations_id` covering all of them.

## See also

[`observation_index()`](https://bbuchsbaum.github.io/crossform/reference/observation_index.md)
for the axes,
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md)
to bind events and confounds to these observations, and
[`file_matrix_source()`](https://bbuchsbaum.github.io/crossform/reference/file_matrix_source.md)
for out-of-memory sources.

Other typed observation facts:
[`bids_study()`](https://bbuchsbaum.github.io/crossform/reference/bids_study.md),
[`observation_confounds()`](https://bbuchsbaum.github.io/crossform/reference/observation_confounds.md),
[`observation_events()`](https://bbuchsbaum.github.io/crossform/reference/observation_events.md),
[`observation_index()`](https://bbuchsbaum.github.io/crossform/reference/observation_index.md),
[`partition_hierarchy()`](https://bbuchsbaum.github.io/crossform/reference/partition_hierarchy.md)

## Examples

``` r
set.seed(1)
domain <- abstract_domain(4L, id = "observations-example")
indexes <- list(
  `run-1` = observation_index(
    1:6, "run-1", time = seq(0, 10, by = 2), units = "seconds"
  ),
  `run-2` = observation_index(
    1:6, "run-2", time = seq(0, 10, by = 2), units = "seconds"
  )
)
responses <- lapply(indexes, function(index) matrix(rnorm(24), 6L, 4L))

record <- observations(responses, indexes, domain)
record$partitions
#> [1] "run-1" "run-2"
record$n_features
#> [1] 4

# No neural values were read: only shapes, axes, and source revisions were
# checked, so a mis-shaped run is caught before any fit is attempted.
short <- responses
short$`run-2` <- short$`run-2`[1:5, , drop = FALSE]
try(observations(short, indexes, domain))
#> Error : Partition `run-2` source dimensions are 5 x 4; its observation index and neural domain require 6 x 4.
```
