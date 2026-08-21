# Declare uncentered inner-product edge geometry

`inner_product()` names the default edge normalization: raw uncentered
products, with no centering, scaling, or correlation normalization
applied to an edge before it is reduced. Use it wherever an operation
stage must record that the geometry is the plain bilinear form.

## Usage

``` r
inner_product()
```

## Value

An `effect_edge_normalizer` recording `$kind` and its `$zero_policy`. It
is a declaration, not a computation.

## See also

[`reduce_partitions()`](https://bbuchsbaum.github.io/crossform/reference/reduce_partitions.md)
and
[`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md)
for the partition stage that follows normalization;
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md)
for the normalized views that are gated behind explicit capabilities.

## Examples

``` r
# The normalizer only names the stage; nothing is computed here.
normalizer <- crossform:::inner_product()
normalizer$kind
#> [1] "inner_product"

# Because no denominator is involved, there is no zero-norm case to
# resolve at execution time.
normalizer$zero_policy
#> [1] "error"
```
