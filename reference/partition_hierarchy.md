# Declare nested partition axes

Columns are ordered from the leaf partition outward. Each child level
must map to exactly one parent level.

## Usage

``` r
partition_hierarchy(data, leaf = names(data)[[1L]], provenance = list())
```

## Arguments

- data:

  A data frame with one row per leaf partition and one column per nested
  axis.

- leaf:

  Name of the leaf partition column; defaults to the first column.

- provenance:

  Portable hierarchy provenance.

## Value

An `effect_partition_hierarchy`: a list with the canonical `$data`, the
`$leaf` column name, the ordered `$axes`, the distinct `$levels` per
axis, the child-to-parent `$parent_maps`, `$provenance`, and a
`$signature`.

## See also

[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md),
which binds the hierarchy to observations, and
[`study_axis()`](https://bbuchsbaum.github.io/crossform/reference/study_axis.md)
to select one axis from it.

Other typed observation facts:
[`bids_study()`](https://bbuchsbaum.github.io/crossform/reference/bids_study.md),
[`observation_confounds()`](https://bbuchsbaum.github.io/crossform/reference/observation_confounds.md),
[`observation_events()`](https://bbuchsbaum.github.io/crossform/reference/observation_events.md),
[`observation_index()`](https://bbuchsbaum.github.io/crossform/reference/observation_index.md),
[`observations()`](https://bbuchsbaum.github.io/crossform/reference/observations.md)

## Examples

``` r
# Runs nested in sessions nested in one subject, leaf column first.
hierarchy <- partition_hierarchy(data.frame(
  partition = c("run-1", "run-2", "run-3", "run-4"),
  session = c("ses-1", "ses-1", "ses-2", "ses-2"),
  subject = "sub-01"
))
hierarchy$axes
#> [1] "partition" "session"   "subject"  
hierarchy$parent_maps$partition
#>   run-1   run-2   run-3   run-4 
#> "ses-1" "ses-1" "ses-2" "ses-2" 

# Nesting must be exact: a session that belongs to two subjects is refused
# rather than quietly flattened.
crossed <- data.frame(
  partition = c("run-1", "run-2"),
  session = c("ses-1", "ses-1"),
  subject = c("sub-01", "sub-02")
)
try(partition_hierarchy(crossed))
#> Error : Hierarchy axis `session` does not nest uniquely within `subject`.
```
