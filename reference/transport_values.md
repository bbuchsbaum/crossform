# Carry native node values onto the group nodes

`transport_values()` applies a
[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md)
to one value per native node — a queried ledger, or a matrix of packed
forms with one row per native node — and returns the group-node values
with the sink last.

## Usage

``` r
transport_values(transport, values)
```

## Arguments

- transport:

  A
  [`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md).

- values:

  One finite value per native node: a numeric vector, or a matrix with
  one row per native node.

## Value

A named numeric vector of length `m + 1` when `values` is a vector,
otherwise an `(m + 1)` by `ncol(values)` matrix. The last entry or row
is the sink, labelled `<sink>`.

## Details

Under `"budget"` semantics the result is `t(P) %*% values`, and the
transported total including the sink equals the native total exactly.
The argument is `P 1 = 1` and never `values >= 0`, so this is a
signed-sum law: it holds for the signed, crossvalidated ledgers a
conservative geometry actually produces. Under `"density"` semantics
each group entry is divided by the transported row mass, group nodes
reached by no native mass come back `NA` rather than `0`, and the sink
row stays in budget units.

## See also

[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md).

Other population transports:
[`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md),
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md),
[`external_transport()`](https://bbuchsbaum.github.io/crossform/reference/external_transport.md),
[`heterogeneity()`](https://bbuchsbaum.github.io/crossform/reference/heterogeneity.md),
[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md),
[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md),
[`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md),
[`population_prevalence()`](https://bbuchsbaum.github.io/crossform/reference/population_prevalence.md),
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md),
[`population_views`](https://bbuchsbaum.github.io/crossform/reference/population_views.md)

## Examples

``` r
transport <- anatomical_transport(
  native_coords = cbind(c(0, 1, 4, 5, 9)),
  group_coords = cbind(c(0, 5)),
  semantics = "budget", radius = 2
)

# Native node 5 sits further than the radius from either group centre, so
# its whole budget lands in the sink and the total still closes.
ledger <- c(2, -1, 0.5, 1.5, 3)
transport_values(transport, ledger)
#> group1 group2 <sink> 
#>      1      2      3 
sum(transport_values(transport, ledger)) - sum(ledger)
#> [1] 0

# A matrix of packed forms transports column by column.
forms <- cbind(a = ledger, b = rev(ledger))
transport_values(transport, forms)
#>        a    b
#> group1 1  4.5
#> group2 2 -0.5
#> <sink> 3  2.0
```
