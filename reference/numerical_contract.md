# Define crossform's numerical reproducibility contract

The contract separates three claims. Reordering completion of fixed
tasks is bitwise invariant because reduction order is canonical.
Changing feature blocks or numerical platform is guaranteed only within
the declared combined absolute and relative tolerance. Bitwise equality
across block partitions or platforms is deliberately not promised.

## Usage

``` r
numerical_contract(atol = 1e-12, rtol = 1e-10)
```

## Arguments

- atol, rtol:

  Nonnegative finite absolute and relative tolerances.

## Value

An `effect_numerical_contract`: the declared `$atol` and `$rtol` plus
one guarantee record per claim (`$scheduling`, `$block_partition`,
`$cross_platform`) and the two explicit non-promises
`$bitwise_across_blocking` and `$bitwise_across_platforms`.

## See also

[`numerical_agreement()`](https://bbuchsbaum.github.io/crossform/reference/numerical_agreement.md),
which tests two results against one of these named guarantees.

Other numerical contracts and receipts:
[`numerical_agreement()`](https://bbuchsbaum.github.io/crossform/reference/numerical_agreement.md)

## Examples

``` r
# The default contract. Read the guarantee you intend to rely on rather
# than assuming bitwise equality everywhere.
contract <- numerical_contract()
contract$scheduling
#> $guarantee
#> [1] "bitwise"
#> 
#> $condition
#> [1] "same tasks and canonical reduction order"
#> 
contract$block_partition
#> $guarantee
#> [1] "tolerance"
#> 
#> $condition
#> [1] "same estimand and precision"
#> 

# Bitwise reproducibility across block partitions is deliberately not
# promised, so a tolerance is the only portable claim.
contract$bitwise_across_blocking
#> [1] FALSE

# Tighten the tolerance when a downstream check needs it.
numerical_contract(atol = 1e-14, rtol = 1e-12)$rtol
#> [1] 1e-12
```
