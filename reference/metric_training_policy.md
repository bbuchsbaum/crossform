# Declare which residual partitions may train a metric

`metric_training_policy()` fixes the partition discipline used when a
metric recipe is estimated from residuals, so metric training cannot
silently borrow the partitions whose products are being evaluated. Pass
it to
[`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md).

## Usage

``` r
metric_training_policy(
  kind = c("exclude_evaluation", "all_partitions_residual_orthogonality"),
  justification = NULL
)
```

## Arguments

- kind:

  `"exclude_evaluation"` estimates each edge metric without its two
  evaluation partitions. `"all_partitions_residual_orthogonality"`
  admits evaluation-partition GLM residuals under an explicit
  orthogonality justification.

- justification:

  Required for the all-partitions policy. It records why residual reuse
  is scientifically admitted; it is not treated as proof.

## Value

An `effect_metric_training_policy` recording `$kind`,
`$includes_evaluation_residuals`, the named `$assumption` it rests on,
the stored `$justification`, and a `$signature` bound into plan
identity.

## Refusal

Requesting `"all_partitions_residual_orthogonality"` without a
`justification` signals an `effect_capability_refusal` carrying
capability `"evaluation_residual_reuse"` in namespace
`"metric_learning"`, with reason
`"residual_reuse_justification_absent"`. Inspect it with
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md).

## See also

[`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md),
which enforces this policy, and
[`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)
for the recipe it governs.

Other neural metrics:
[`diagonal_precision()`](https://bbuchsbaum.github.io/crossform/reference/diagonal_precision.md),
[`gaussian_covariance_model()`](https://bbuchsbaum.github.io/crossform/reference/gaussian_covariance_model.md),
[`identity_metric()`](https://bbuchsbaum.github.io/crossform/reference/identity_metric.md),
[`metric_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/metric_capabilities.md),
[`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md),
[`noise_precision()`](https://bbuchsbaum.github.io/crossform/reference/noise_precision.md),
[`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)

## Examples

``` r
# The default trains each edge's metric without its two evaluation
# partitions, so the metric and the products stay disjoint.
policy <- metric_training_policy()
policy$kind
#> [1] "exclude_evaluation"
c(reuses_evaluation = policy$includes_evaluation_residuals,
  assumption = policy$assumption)
#>                    reuses_evaluation                           assumption 
#>                              "FALSE" "partition_disjoint_metric_training" 

# Reusing evaluation-partition residuals is admitted only with an explicit
# written justification, which is recorded, not verified.
permissive <- metric_training_policy(
  "all_partitions_residual_orthogonality",
  justification = "GLM residuals are orthogonal to the fitted effects."
)
permissive$includes_evaluation_residuals
#> [1] TRUE

# Omitting that justification is refused, and the refusal is classed, so a
# caller can branch on the capability rather than on the prose.
refusal <- catch_refusal(
  metric_training_policy("all_partitions_residual_orthogonality")
)
refusal$capability
#> [1] "evaluation_residual_reuse"
refusal$reasons
#> [1] "residual_reuse_justification_absent"
```
