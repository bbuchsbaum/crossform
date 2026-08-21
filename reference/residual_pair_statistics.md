# Accumulate canonical residual pair sufficient statistics

Computes residual cross-products only for neural feature pairs that
coexist in at least one requested support. The GEMM tile shape is
derived from the support topology and is independent of
`workspace_bytes`; the budget changes only how many already-computed
residual tiles may be cached. Consequently, changing the workspace plan
cannot change the accumulated floating-point values.

## Usage

``` r
residual_pair_statistics(
  x,
  at,
  partitions = NULL,
  workspace_bytes = 512 * 1024^2
)
```

## Arguments

- x:

  An `effect_relation_fit` with residual-block capability.

- at:

  A compiled frame carrying an explicit support index, such as a
  compiled searchlight frame. An internal `effect_support_index` is also
  accepted.

- partitions:

  Optional relation partitions. They are canonicalized to relation
  order.

- workspace_bytes:

  Positive crossform-owned workspace budget.

## Value

An `effect_residual_pair_statistics` object. `$pair_i`/`$pair_j` list
the coexisting feature pairs, `$partitions` names the partitions,
`$atomic` holds one `$cross_products` vector and `$residual_df` per
partition, and `$numerical_contract` records the fixed tile shape.
`$execution` diagnostics are excluded from the scientific `$signature`.

## See also

[`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md),
which compiles these statistics into a learned metric schedule, and
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md),
which can reuse them through
`residual_strategy = "shared_pair_statistics"`.

## Examples

``` r
# Residual cross-products are accumulated only for feature pairs that
# coexist in some searchlight, so the cost tracks support topology rather
# than the square of the feature count.
example <- example_fmri_effects()
statistics <- crossform:::residual_pair_statistics(
  example$fit, example$frame
)
length(statistics$pair_i)
#> [1] 2763
statistics$partitions
#> [1] "run1" "run2" "run3" "run4"

# One atomic record per partition, each carrying its own residual df.
statistics$atomic[["run1"]]$residual_df
#> [1] 28

# The workspace budget is a cache size, not part of the numerical shape,
# so shrinking it cannot change the accumulated values.
frugal <- crossform:::residual_pair_statistics(
  example$fit, example$frame, workspace_bytes = 4 * 1024^2
)
identical(frugal$signature, statistics$signature)
#> [1] TRUE
```
