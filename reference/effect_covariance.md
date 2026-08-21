# Read the unscaled effect-coordinate covariance of a fitted partition

This is the design-side factor of the separable error model. By
convention the neural residual covariance is excluded, so the value
depends on the design and whitener but not on the data.

## Usage

``` r
effect_covariance(x, partition)
```

## Arguments

- x:

  An `effect_relation_fit` with residual-block capability.

- partition:

  One partition name or index.

## Value

A symmetric effect-by-effect matrix, with rows and columns named by the
relation's effect coordinates, excluding the neural residual covariance
factor.

## See also

[`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md)
and
[`residual_block()`](https://bbuchsbaum.github.io/crossform/reference/residual_block.md)
for the other two pieces of the error channel, and
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md),
which combines them.

## Examples

``` r
example <- example_fmri_effects()

# Four condition means estimated from eight trials each: the design factor
# is diagonal with entries 1/8, because the conditions are orthogonal.
covariance <- crossform:::effect_covariance(example$fit, "run1")
round(covariance, 4)
#>        face  body house  tool
#> face  0.125 0.000 0.000 0.000
#> body  0.000 0.125 0.000 0.000
#> house 0.000 0.000 0.125 0.000
#> tool  0.000 0.000 0.000 0.125

# Scaling it by a residual variance gives an effect standard error.
residuals <- residual_block(example$fit, "run1", 1L)
variance <- sum(residuals^2) / residual_df(example$fit, "run1")
round(sqrt(diag(covariance) * variance), 3)
#>  face  body house  tool 
#> 0.252 0.252 0.252 0.252 
```
