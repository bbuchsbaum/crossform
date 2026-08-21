# Import BIDS task events as typed event facts

The adapter preserves arbitrary BIDS columns and adds private partition
and event-key columns required by the generic
[`observation_events()`](https://bbuchsbaum.github.io/crossform/reference/observation_events.md)
contract. Partition identity is explicit rather than inferred from
filenames.

## Usage

``` r
bids_events(files, partitions = names(files), units = "seconds")
```

## Arguments

- files:

  Character event-TSV paths, one per partition.

- partitions:

  Explicit ordered partition identifiers. Named `files` may supply these
  names.

- units:

  Physical onset and duration units.

## Value

An
[`observation_events()`](https://bbuchsbaum.github.io/crossform/reference/observation_events.md)
fact object whose `$data` holds every original BIDS column plus the
private `.bids_partition` and `.bids_event_id` key columns, and whose
`$provenance` records each file's basename and content revision.

## See also

[`observation_events()`](https://bbuchsbaum.github.io/crossform/reference/observation_events.md)
for the generic contract,
[`bids_confounds()`](https://bbuchsbaum.github.io/crossform/reference/bids_confounds.md)
for the confound side, and
[`bids_study()`](https://bbuchsbaum.github.io/crossform/reference/bids_study.md)
to bind both into a study.

## Examples

``` r
# Stand in for one run's events.tsv.
path <- tempfile(fileext = ".tsv")
utils::write.table(
  data.frame(
    onset = c(0, 6), duration = 0.5, trial_type = c("face", "body")
  ),
  path, sep = "\t", row.names = FALSE, quote = FALSE
)

# Partition identity is declared by the caller, never parsed from the
# filename, and every original column survives the import.
record <- crossform:::bids_events(c(`run-1` = path))
record$timing
#> [1] TRUE
names(record$data)
#> [1] "onset"           "duration"        "trial_type"      ".bids_partition"
#> [5] ".bids_event_id" 
record$data$.bids_event_id
#> [1] "event-000001" "event-000002"
unlink(path)
```
