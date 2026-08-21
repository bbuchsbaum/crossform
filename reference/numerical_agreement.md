# Assess results under a named numerical guarantee

Tolerance agreement uses
`abs(x - y) <= atol + rtol * max(abs(x), abs(y))` elementwise.

## Usage

``` r
numerical_agreement(
  x,
  y,
  guarantee = c("scheduling", "block_partition", "cross_platform"),
  contract = numerical_contract()
)
```

## Arguments

- x, y:

  Numeric objects with identical dimensions.

- guarantee:

  One of `scheduling`, `block_partition`, or `cross_platform`.
  `scheduling` is checked bitwise; the other two are checked against the
  contract tolerance.

- contract:

  An `effect_numerical_contract`.

## Value

An `effect_numeric_agreement` with `$passed`, the `$guarantee` and
`$comparison` used (`"bitwise"` or `"tolerance"`),
`$max_absolute_error`, `$max_allowed_error`, and the `$atol`/`$rtol` in
force.

## See also

[`numerical_contract()`](https://bbuchsbaum.github.io/crossform/reference/numerical_contract.md)
for the guarantees themselves.

Other numerical contracts and receipts:
[`numerical_contract()`](https://bbuchsbaum.github.io/crossform/reference/numerical_contract.md)

## Examples

``` r
# Re-running the same fixed plan under a different feature-block size is a
# `block_partition` change: crossform promises tolerance, not bitwise
# equality.
baseline <- c(1, 2, 3) / 7
reblocked <- baseline + c(0, 1e-15, -2e-15)
agreement <- numerical_agreement(baseline, reblocked, "block_partition")
agreement$passed
#> [1] TRUE
agreement$max_absolute_error
#> [1] 1.998401e-15

# The same pair fails the strict bitwise scheduling guarantee, which is
# reserved for reordering completion of identical tasks.
numerical_agreement(baseline, reblocked, "scheduling")$passed
#> [1] FALSE
numerical_agreement(baseline, baseline, "scheduling")$passed
#> [1] TRUE
```
