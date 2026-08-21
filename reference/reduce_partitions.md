# Reduce normalized and transformed partition edges

This is the default estimand: normalize each declared edge, transform
it, and only then apply the partition weights. Choose it when the
quantity you mean is the average of per-partition normalized edges.

## Usage

``` r
reduce_partitions()
```

## Value

An `effect_partition_reducer` recording `$kind`, `$weight_convention`,
and the stage `$order` (`"edge_first"`).

## See also

[`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md),
the other stage order, and
[`pairing()`](https://bbuchsbaum.github.io/crossform/reference/pairing.md),
which supplies the partition weights this reducer applies.

Other generalization pairings:
[`cross_partitions()`](https://bbuchsbaum.github.io/crossform/reference/cross_partitions.md),
[`pairing()`](https://bbuchsbaum.github.io/crossform/reference/pairing.md)

## Examples

``` r
# Edge-first: each partition pair is normalized and transformed before the
# pairing weights are applied.
reducer <- reduce_partitions()
reducer$order
#> [1] "edge_first"

# The two reducers name different estimands, so they are not
# interchangeable defaults.
c(edge_first = reduce_partitions()$order,
  aggregate_first = aggregate_first()$order)
#>        edge_first   aggregate_first 
#>      "edge_first" "aggregate_first" 
```
