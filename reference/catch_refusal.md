# Inspect a capability refusal

Every contract-level refusal in crossform signals a condition of class
`effect_capability_refusal`. The condition carries the missing
`$capability`, its `$namespace`, all unmet `$reasons`, and suggested
`$remedies`, so callers can branch on the *cause* rather than matching
message prose.

## Usage

``` r
catch_refusal(expr)
```

## Arguments

- expr:

  An expression expected to refuse.

## Value

The captured `effect_capability_refusal` condition, or `NULL` if `expr`
succeeded.

## See also

[crossform_conditions](https://bbuchsbaum.github.io/crossform/reference/crossform_conditions.md)
for the other condition classes crossform raises, and how to branch on
them.

Other conditions:
[`crossform_conditions`](https://bbuchsbaum.github.io/crossform/reference/crossform_conditions.md)

## Examples

``` r
domain <- abstract_domain(2, id = "refusal-example")
relation <- relation(
  list(a = matrix(1:4, 2), b = matrix(2:5, 2)),
  effects = c("x", "y"), domain = domain
)
refusal <- catch_refusal(
  rdm_sampling_covariance(
    plan_geometry(relation, compile_frame(whole_brain(), domain),
      cross_partitions(relation, independence = "independent")),
    relation, target = "null", at = 1L
  )
)
refusal$capability
#> [1] "sampling_covariance"
refusal$reasons
#> [1] "missing_error_channel"                
#> [2] "sampling_axis_missing_or_inconsistent"
```
