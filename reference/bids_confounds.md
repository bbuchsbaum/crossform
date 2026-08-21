# Import fMRIPrep-style confounds as observation facts

All columns are preserved without assigning them target or nuisance
roles. If `censor` is supplied it must name a complete logical retain
column; no censor policy is inferred from motion or outlier columns.

## Usage

``` r
bids_confounds(
  files,
  partitions = names(files),
  observation_ids = NULL,
  censor = NULL
)
```

## Arguments

- files:

  Character confound-TSV paths, one per partition.

- partitions:

  Explicit ordered partition identifiers.

- observation_ids:

  Optional named list of observation identifiers. Row numbers are used
  when omitted.

- censor:

  Optional logical retain-column name already present in every input
  table.

## Value

An
[`observation_confounds()`](https://bbuchsbaum.github.io/crossform/reference/observation_confounds.md)
fact object whose `$data` holds every original confound column plus the
private `.bids_partition` and `.bids_observation_id` key columns, with
`$censor_column` set only when `censor` was supplied.

## See also

[`observation_confounds()`](https://bbuchsbaum.github.io/crossform/reference/observation_confounds.md)
for the generic contract,
[`bids_events()`](https://bbuchsbaum.github.io/crossform/reference/bids_events.md)
for the event side, and
[`bids_study()`](https://bbuchsbaum.github.io/crossform/reference/bids_study.md)
to bind both.

## Examples

``` r
# Stand in for one run's fMRIPrep confounds table.
path <- tempfile(fileext = ".tsv")
utils::write.table(
  data.frame(
    framewise_displacement = c(0, 0.1, 0.9, 0.2),
    retained = c(TRUE, TRUE, FALSE, TRUE)
  ),
  path, sep = "\t", row.names = FALSE, quote = FALSE
)

# Naming the retain column is required: no censor policy is inferred from
# the motion columns, however large they are.
record <- crossform:::bids_confounds(c(`run-1` = path), censor = "retained")
record$censor_column
#> [1] "retained"
sum(record$data$retained)
#> [1] 3

# Naming a column that is absent or not logical refuses explicitly.
catch_refusal(
  crossform:::bids_confounds(
    c(`run-1` = path), censor = "framewise_displacement"
  )
)$capability
#> [1] "censoring_declared"
unlink(path)
```
