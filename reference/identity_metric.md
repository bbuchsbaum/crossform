# Specify an on-demand identity metric

These constructors describe how a local neural metric is to be produced.
They do not allocate a domain-wide matrix. Use
[`compile_metric_schedule()`](https://bbuchsbaum.github.io/crossform/reference/compile_metric_schedule.md)
to bind a recipe to residual sufficient statistics, spatial supports,
and evaluation edges.

## Usage

``` r
identity_metric(domain = NULL)
```

## Arguments

- domain:

  Optional exact neural domain. Omitting it defers domain binding until
  schedule compilation.

## Value

An `effect_metric_recipe` with `$kind`, the optional bound `$domain`, a
`$capabilities` record, the fixed `$hyperparameters` (estimator,
randomness, seed), and a `$signature`. It allocates no matrix.

## See also

[`diagonal_precision()`](https://bbuchsbaum.github.io/crossform/reference/diagonal_precision.md)
and
[`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)
for residual-derived recipes,
[`metric_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/metric_capabilities.md)
to inspect one, and
[`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md),
which binds a recipe to residual statistics.

Other neural metrics:
[`diagonal_precision()`](https://bbuchsbaum.github.io/crossform/reference/diagonal_precision.md),
[`gaussian_covariance_model()`](https://bbuchsbaum.github.io/crossform/reference/gaussian_covariance_model.md),
[`metric_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/metric_capabilities.md),
[`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md),
[`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md),
[`noise_precision()`](https://bbuchsbaum.github.io/crossform/reference/noise_precision.md),
[`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)

## Examples

``` r
# A recipe is a declaration, not a matrix: no domain-wide operator exists
# until a schedule binds it to supports.
recipe <- identity_metric()
recipe$kind
#> [1] "identity"
recipe$hyperparameters$estimator
#> [1] "fixed_identity"
metric_capabilities(recipe)$materialized
#> [1] FALSE

# Bind the domain up front when you want a domain mismatch caught before
# schedule compilation rather than during it.
domain <- abstract_domain(3, id = "identity-metric-example")
identity_metric(domain)$domain$id
#> [1] "identity-metric-example"

# Unlike the residual-derived recipes, this one estimates nothing, so it
# stays diagonal and needs no residual channel.
c(identity = metric_capabilities(recipe)$native_diagonal,
  shrinkage = metric_capabilities(shrinkage_precision())$native_diagonal)
#>  identity shrinkage 
#>      TRUE     FALSE 
```
